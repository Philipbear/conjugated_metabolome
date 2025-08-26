import argparse
import multiprocessing as mp
import os
import shlex
import time
from typing import List

import numpy as np
import pandas as pd
from pyteomics import mzml, mzxml
from tqdm import tqdm


def process_file(args):
    try:
        file_path, title_name = args
        file_format = os.path.splitext(file_path)[1].lower()
        file_name = os.path.splitext(os.path.basename(file_path))[0]

        output = []
        if file_format == '.mzml':
            reader = mzml.MzML(file_path)
            output = process_mzml(reader, title_name, file_name)
        elif file_format == '.mzxml':
            reader = mzxml.MzXML(file_path)
            output = process_mzxml(reader, title_name, file_name)
        return output
    except:
        return []


def process_mzml(reader, title_name, file_name):
    output = []
    for spectrum in reader:
        if spectrum['ms level'] == 2:  # MS2 scans only
            precursor_mz = round(float(
                spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z']), 4)

            if precursor_mz > 2000:
                continue

            mz_array = np.array(spectrum['m/z array'])
            intensity_array = np.array(spectrum['intensity array'])

            if len(mz_array) < 5:
                continue

            # >1% intensity
            mask = (intensity_array > np.max(intensity_array) * 0.01) & (mz_array < precursor_mz + 0.5)
            mz_array = mz_array[mask]
            intensity_array = intensity_array[mask]

            if len(mz_array) < 5:
                continue

            peaks = np.column_stack((mz_array, intensity_array))
            # Limit the number of peaks
            max_peaks = 50
            if len(peaks) > max_peaks:
                peaks = peaks[peaks[:, 1].argsort()[-max_peaks:][::-1]]
            # sort by mz
            peaks = peaks[peaks[:, 0].argsort()]
            # normalize intensity
            peaks[:, 1] = peaks[:, 1] / np.max(peaks[:, 1]) * 999

            scan_number = spectrum['index'] + 1
            try:
                charge = int(
                    spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['charge state'])
            except:
                charge = 1

            if charge != 1:  # Only consider singly charged spectra
                continue

            polarity = '+' if 'positive scan' in spectrum else '-'
            charge_str = f"{charge}{polarity}"

            output.append({
                'title': f"{title_name}:{file_name}:scan:{scan_number}",
                'pepmass': precursor_mz,
                'charge': charge_str,
                'peaks': peaks
            })

    return output


def process_mzxml(reader, title_name, file_name):
    output = []
    for spectrum in reader:
        if spectrum['msLevel'] == 2:  # MS2 scans only
            precursor_mz = round(float(spectrum['precursorMz'][0]['precursorMz']), 4)

            if precursor_mz > 2000:
                continue

            mz_array = np.array(spectrum['m/z array'])
            intensity_array = np.array(spectrum['intensity array'])

            if len(mz_array) < 5:
                continue

            # >1% intensity
            mask = (intensity_array > np.max(intensity_array) * 0.01) & (mz_array < precursor_mz + 0.5)
            mz_array = mz_array[mask]
            intensity_array = intensity_array[mask]

            if len(mz_array) < 5:
                continue

            peaks = np.column_stack((mz_array, intensity_array))
            # Limit the number of peaks
            max_peaks = 50
            if len(peaks) > max_peaks:
                peaks = peaks[peaks[:, 1].argsort()[-max_peaks:][::-1]]
            # sort by mz
            peaks = peaks[peaks[:, 0].argsort()]
            # normalize intensity
            peaks[:, 1] = peaks[:, 1] / np.max(peaks[:, 1]) * 999

            scan_number = spectrum['num']

            polarity = spectrum['polarity']
            if polarity not in ['+', '-']:
                polarity = '+'

            try:
                charge = int(spectrum['precursorMz'][0]['precursorCharge'])
            except:
                charge = 1

            if charge != 1:  # Only consider singly charged spectra
                continue

            charge_str = f"{charge}{polarity}"

            output.append({
                'title': f"{title_name}:{file_name}:scan:{scan_number}",
                'pepmass': precursor_mz,
                'charge': charge_str,
                'peaks': peaks
            })

    return output


def write_to_mgf(spectra, output_file, mode='a'):
    """Write spectra to an MGF file, either appending or writing a new file."""
    with open(output_file, mode) as mgf_out:
        for spectrum in spectra:
            mgf_out.write("BEGIN IONS\n")
            mgf_out.write(f"TITLE={spectrum['title']}\n")
            mgf_out.write(f"PEPMASS={spectrum['pepmass']}\n")
            mgf_out.write(f"CHARGE={spectrum['charge']}\n")

            np.savetxt(mgf_out, spectrum['peaks'], fmt="%.4f %.0f")

            mgf_out.write("END IONS\n\n")


def chunk_files(input_files: List[str], chunk_size: int) -> List[List[str]]:
    """Split the input files into chunks of specified size."""
    return [input_files[i:i + chunk_size] for i in range(0, len(input_files), chunk_size)]


def send_mgf_to_server(file_path, server_path):
    """Send a file to the server using rsync and remove it locally."""
    cmd = f'rsync -avz -e "ssh -i /home/user/.ssh/id_rsa" --progress {file_path} {server_path} && rm {file_path}'
    os.system(cmd)
    return


def convert_raw_data_to_mgf(input_files: List[str], title_name: str, output_file: str, chunk_size: int = 10,
                            single_file: bool = False, send_to_server: bool = False, server_path: str = None):
    """
    Convert a list of mzML or mzXML files to MGF file(s) using parallel processing and chunking.

    Args:
    input_files (list): List of input mzML or mzXML file paths.
    title_name (str): Name to be used in the title of each spectrum.
    output_file (str): Path to the output MGF file (or base name for multiple files).
    chunk_size (int): Number of files to process in each chunk.
    single_file (bool): If True, write all results to a single file. If False, write a separate file for each chunk.
    send_to_server (bool): If True, send each chunk file to the server immediately after generation.
    server_path (str): The path on the server where files should be sent.
    """
    output_dir = os.path.dirname(output_file)
    os.makedirs(output_dir, exist_ok=True)

    if single_file:
        open(output_file, 'w').close()

    chunked_files = chunk_files(input_files, chunk_size)

    total_chunks = len(chunked_files)
    for i, chunk in enumerate(chunked_files, 1):
        print(f"Processing chunk {i}/{total_chunks} (contains {len(chunk)} files)...")
        if single_file:
            process_chunk(chunk, title_name, output_file)
        else:
            chunk_output_file = f"{os.path.splitext(output_file)[0]}_chunk{i}.mgf"
            process_chunk(chunk, title_name, chunk_output_file, mode='w')
            if send_to_server:
                try:
                    send_mgf_to_server(chunk_output_file, server_path)
                except:
                    print(f"Error sending {chunk_output_file} to server.")

    if single_file:
        print(f"All chunks processed. Output is stored in {output_file}")
    else:
        print(
            f"All chunks processed. Output files are stored in {output_dir} with prefix {os.path.basename(output_file)}")
        if send_to_server:
            print("All chunk files have been sent to the server and removed locally.")


def process_chunk(chunk: List[str], title_name: str, output_file: str, mode='a'):
    """Process a chunk of files and write the results to an MGF file."""
    pool = mp.Pool(processes=mp.cpu_count())

    args = [(file_path, title_name) for file_path in chunk]

    results = []
    for result in pool.imap_unordered(process_file, args):
        results.extend(result)

    pool.close()
    pool.join()

    write_to_mgf(results, output_file, mode=mode)


def main(repo='ST', send_to_server=True):
    """
    :param repo: ST, MTBLS, MassIVE
    :param send_to_server: If True, send each chunk file to the server immediately after generation.
    """

    start_time = time.time()

    output_dir = f'/home/user/shipei/clustering/{repo}'
    os.makedirs(output_dir, exist_ok=True)

    # Create or load the processed file list
    processed_list_file = os.path.join(output_dir, f"{repo}_processed_files.txt")
    if os.path.exists(processed_list_file):
        with open(processed_list_file, 'r') as f:
            processed_files = set(f.read().splitlines())
    else:
        processed_files = set()

    # Create error file list
    error_list_file = os.path.join(output_dir, f"{repo}_error_files.txt")

    df = pd.read_csv('all_filtered.tsv', sep='\t')
    df = df[df['Repo'] == repo].reset_index(drop=True)

    # get unique datasets
    datasets = df['Dataset'].unique()

    server_path = "shipei@137.110.133.146:/home/shipei/projects/cluster/unclustered_all"

    for dataset in tqdm(datasets, desc='Datasets', total=len(datasets)):
        if dataset in processed_files:
            print(f"Skipping {dataset} as it has already been processed.")
            continue

        dataset_df = df[df['Dataset'] == dataset].reset_index(drop=True)
        all_files = [shlex.quote(filename) for filename in dataset_df['Filename']]

        print(f'Converting raw data to MGF for {dataset}...')
        _mgf_file = f'{output_dir}/{dataset}.mgf'
        convert_raw_data_to_mgf(all_files, dataset, _mgf_file, chunk_size=10,
                                send_to_server=send_to_server, server_path=server_path)

        # Add to processed files list if success
        processed_files.add(dataset)
        with open(processed_list_file, 'a') as f:
            f.write(f"{dataset}\n")
        print(f"Processed {dataset}")

    print(f"Total time: {(time.time() - start_time) / 60} minutes")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Get MS2 spectra from raw data files')
    parser.add_argument('--repo', '-r', type=str, default='ST', help='Repository name, ST, MTBLS, or MassIVE')
    args = parser.parse_args()

    main(args.repo)

    # convert_raw_data_to_mgf(['/Users/shipei/Documents/projects/playground/falcon/22_BB4_01_27191.d.mzML'], 'ST000001', '/Users/shipei/Documents/projects/playground/falcon/22.mgf')
