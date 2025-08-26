# Usage: python gnps_download_cluster.py -m new
import os
import pandas as pd
import shutil
import shlex
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, TimeoutError
from queue import Queue
from threading import Thread
from collections import defaultdict
import argparse
from get_ms2_from_raw_data import convert_raw_data_to_mgf


MAX_DOWNLOADED_DATASETS = 25
PROCESSED_DATASETS_FILE = "processed_datasets.txt"
ERROR_DATASETS_FILE = "error_datasets.txt"

download_queue = Queue()
processing_queue = Queue()


def get_file_polarity_centroid(file_name, timeout=10):
    def process_file():
        ion_mode = 'unknown'
        centroided = 'unknown'
        ext = os.path.splitext(file_name)[1].lower()

        try:
            with open(file_name, 'r', encoding='iso-8859-1') as f:
                if ext == '.mzml':
                    for line_num, line in enumerate(f, 1):
                        if line_num > 300:
                            break

                        lower_line = line.lower()
                        if ion_mode == 'unknown':
                            if 'positive' in lower_line:
                                ion_mode = 'positive'
                            elif 'negative' in lower_line:
                                ion_mode = 'negative'

                        if centroided == 'unknown':
                            if 'centroid spectrum' in lower_line:
                                centroided = 'centroided'
                            elif 'profile spectrum' in lower_line:
                                centroided = 'profile'

                        if ion_mode != 'unknown' and centroided != 'unknown':
                            break

                elif ext == '.mzxml':
                    for line_num, line in enumerate(f, 1):
                        if line_num > 300:
                            break

                        lower_line = line.lower()
                        if ion_mode == 'unknown':
                            if 'polarity="+"' in lower_line:
                                ion_mode = 'positive'
                            elif 'polarity="-"' in lower_line:
                                ion_mode = 'negative'

                        if centroided == 'unknown':
                            if 'centroided="1"' in lower_line:
                                centroided = 'centroided'
                            elif 'centroided="0"' in lower_line:
                                centroided = 'profile'

                        if ion_mode != 'unknown' and centroided != 'unknown':
                            break

        except Exception as e:
            print(f"Error processing {file_name}: {str(e)}")

        return file_name, ion_mode, centroided

    with ThreadPoolExecutor() as executor:
        future = executor.submit(process_file)
        try:
            return future.result(timeout=timeout)
        except TimeoutError:
            print(f"Processing of {file_name} timed out")
            return file_name, 'unknown', 'unknown'


def prepare_datasets_to_download():

    df = pd.read_csv('data/gnps_new_datasets_filtered.tsv', sep='\t')
    unique_msv_ids = df['dataset'].unique()

    downloaded_msvs = [d for d in os.listdir(DOWNLOAD_DIR) if d.startswith('MSV')]

    to_download = [msv for msv in unique_msv_ids if msv not in downloaded_msvs]

    return to_download


def run_lftp_command(command):
    full_command = f'lftp -c "open massive-ftp.ucsd.edu; {command}"'
    result = subprocess.run(full_command, shell=True, capture_output=True, text=True)
    return result.stdout.strip().split('\n')


def find_and_download_files(msv_id):
    base_dirs = ['z01', 'v09', 'v08', 'v07', 'v06', 'v05', 'v04', 'v03', 'v02', 'v01', 'x01']
    base_dir = None
    for dir in base_dirs:
        all_dirs = run_lftp_command(f'ls {dir}')
        for entry in all_dirs:
            if msv_id in entry:
                base_dir = dir
                break
        if base_dir:
            break

    if not base_dir:
        print(f"Could not find base directory for {msv_id}")
        return

    msv_path = f"{base_dir}/{msv_id}"
    all_files = collect_files(msv_path)
    deduplicated_files = deduplicate_files(all_files)
    download_files(msv_id, deduplicated_files)


def collect_files(remote_dir):
    all_files = []

    def process_dir(current_dir):
        files = run_lftp_command(f'ls -l {shlex.quote(current_dir)}')
        for file in files:
            file_parts = file.split(maxsplit=8)
            if len(file_parts) < 9:
                continue
            file_type = file_parts[0][0]
            file_name = file_parts[8]

            if file_type == 'd':
                process_dir(f"{current_dir}/{shlex.quote(file_name)}")
            elif file_name.lower().endswith(('.mzml', '.mzxml')):
                all_files.append(f"{current_dir}/{file_name}")

    process_dir(remote_dir)
    return all_files


def deduplicate_files(file_list):
    file_dict = defaultdict(list)
    for file_path in file_list:
        base_name = os.path.basename(file_path)
        file_dict[base_name].append(file_path)

    return [paths[0] for paths in file_dict.values()]


def download_files(msv_id, file_list):
    local_msv_path = os.path.join(DOWNLOAD_DIR, msv_id)
    os.makedirs(local_msv_path, exist_ok=True)

    for remote_file_path in file_list:
        file_name = os.path.basename(remote_file_path)
        local_file_path = os.path.join(local_msv_path, file_name)  # DOWNLOAD_DIR/MSVXXXXX/file_name

        cmd = f'lftp massive-ftp.ucsd.edu -e "get {shlex.quote(remote_file_path)} -o {shlex.quote(local_file_path)}; exit"'
        try:
            subprocess.run(cmd, shell=True, check=True)
            print(f"Downloaded: {file_name}")
        except subprocess.CalledProcessError:
            print(f"Failed to download: {remote_file_path}")


def get_msv_count():
    return len([d for d in os.listdir(DOWNLOAD_DIR) if d.startswith('MSV')])


def downloader_thread():
    while True:
        try:
            if download_queue.qsize() > 0:
                msv_count = get_msv_count()
                if msv_count < MAX_DOWNLOADED_DATASETS:
                    msv_id = download_queue.get()
                    find_and_download_files(msv_id)
                    processing_queue.put(msv_id)
                    print(f"Downloaded {msv_id}. Current MSV count: {get_msv_count()}")
            time.sleep(60)  # Check every 1 minute
        except Exception as e:
            print(f"Error in downloader thread: {str(e)}")


def get_processed_datasets():
    if os.path.exists(PROCESSED_DATASETS_FILE):
        with open(PROCESSED_DATASETS_FILE, 'r') as f:
            return set(line.strip() for line in f)
    return set()


def add_processed_dataset(msv_id):
    with open(PROCESSED_DATASETS_FILE, 'a') as f:
        f.write(f"{msv_id}\n")


def add_error_dataset(msv_id, error_message):
    with open(ERROR_DATASETS_FILE, 'a') as f:
        f.write(f"{msv_id}: {error_message}\n")


def initialize_queues():
    downloaded = set(d for d in os.listdir(DOWNLOAD_DIR) if d.startswith('MSV'))
    processed = get_processed_datasets()

    # Add downloaded but unprocessed datasets to processing queue
    for msv_id in downloaded:
        if msv_id not in processed:
            processing_queue.put(msv_id)

    # Prepare datasets to download
    to_download = prepare_datasets_to_download()

    # Add datasets that are neither downloaded nor processed to download queue
    for msv_id in to_download:
        if msv_id not in downloaded and msv_id not in processed:

            ###############
            # skip some huge datasets which were examined and found to be not useful/proteomics
            ###############
            if msv_id in ['MSV000096103', 'MSV000095534', 'MSV000095357', 'MSV000095360', 'MSV000095324']:
                continue
            ###############

            download_queue.put(msv_id)


def processor_thread():
    while True:
        try:
            if not processing_queue.empty():
                msv_id = processing_queue.get()
                msv_path = os.path.join(DOWNLOAD_DIR, msv_id)

                if not os.path.exists(msv_path):
                    print(f"MSV directory {msv_id} not found. Marking as error and skipping.")
                    add_error_dataset(msv_id, "Directory not found")
                    continue

                # All files in the MSV directory
                all_files = [f for f in os.listdir(msv_path) if f.lower().endswith(('.mzml', '.mzxml'))]
                all_files = [os.path.join(msv_path, f) for f in all_files]

                print(f'Converting raw data to MGF for {msv_id}...')
                # convert raw data to MGF
                _mgf_file = os.path.join(RESULT_DIR, f'{msv_id}.mgf')

                try:
                    convert_raw_data_to_mgf(all_files, msv_id, _mgf_file, chunk_size=10)

                    # Mark as processed
                    add_processed_dataset(msv_id)
                    print(f"Processed {msv_id}. MGF file: {_mgf_file}")

                except Exception as e:
                    error_message = f"Error processing {msv_id}: {str(e)}"
                    print(error_message)
                    add_error_dataset(msv_id, error_message)

                # Remove the raw data files and the directory to save space
                shutil.rmtree(msv_path)
                print(f"Removed raw files and directory for {msv_id}.")

            else:
                time.sleep(10)  # Wait for 10 seconds if the queue is empty
        except Exception as e:
            error_message = f"Error in processor thread: {str(e)}"
            print(error_message)
            add_error_dataset(msv_id, error_message)


def main():

    os.makedirs(DOWNLOAD_DIR, exist_ok=True)
    os.makedirs(RESULT_DIR, exist_ok=True)

    initialize_queues()

    # Start downloader thread
    Thread(target=downloader_thread, daemon=True).start()

    # Start processor thread
    Thread(target=processor_thread, daemon=True).start()

    # Main loop to keep the script running
    while True:
        try:
            time.sleep(10000)
            print(f"Download queue size: {download_queue.qsize()}")
            print(f"Processing queue size: {processing_queue.qsize()}")
            print(f"Downloaded datasets: {get_msv_count()}")
        except Exception as e:
            print(f"Error in main loop: {str(e)}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process GNPS datasets.')
    parser.add_argument('--local', action='store_true', help='Run in local mode')
    args = parser.parse_args()

    if args.local:
        DOWNLOAD_DIR = "/clustering/data/gnps"
        RESULT_DIR = "/clustering/data/gnps_result"
    else:
        DOWNLOAD_DIR = "/home/shipei/projects/cluster/msv_datasets"
        RESULT_DIR = "/home/shipei/projects/cluster/new_unclustered_msv"

    main()
