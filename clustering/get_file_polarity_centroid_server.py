import os
import pandas as pd
import argparse
from multiprocessing import Pool, cpu_count as get_cpu_count
from tqdm import tqdm
import time
import signal


def timeout_handler(signum, frame):
    raise TimeoutError("File processing timed out")


def get_file_polarity_centroid(file_name, timeout=10):
    """
    Find the polarity (ion mode) and centroiding status of the file.
    Stops processing after 400 lines. Uses ISO-8859-1 encoding.
    Implements a timeout mechanism.
    """
    ion_mode = 'unknown'
    centroided = 'unknown'
    ext = os.path.splitext(file_name)[1].lower()

    signal.signal(signal.SIGALRM, timeout_handler)
    signal.alarm(timeout)

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
    except TimeoutError:
        print(f"Processing of {file_name} timed out")
    except Exception as e:
        print(f"Error processing {file_name}: {str(e)}")
    finally:
        signal.alarm(0)

    return file_name, ion_mode, centroided


def process_chunk(chunk_data):
    chunk, timeout = chunk_data
    return [get_file_polarity_centroid(file_name, timeout) for file_name in chunk]


def main(repo='msv', cpu_count=None, chunk_size=1000, timeout=10):
    """
    Process the repository data with multiprocessing and a progress bar.
    """
    df = pd.read_csv(f'{repo}_results.tsv', sep='\t')
    df['Filename'] = df['Filename'].apply(lambda x: '/data/datasets/server/' + str(x))
    df['Polarity'] = 'unknown'
    df['Centroided'] = 'unknown'

    # Determine the number of CPUs to use
    max_cpu_count = get_cpu_count()
    if cpu_count is None or cpu_count <= 0:
        cpu_count = max_cpu_count
    else:
        cpu_count = min(cpu_count, max_cpu_count)

    print(f"Using {cpu_count} CPU cores")

    # Prepare data for multiprocessing
    file_names = df['Filename'].tolist()
    chunks = [file_names[i:i + chunk_size] for i in range(0, len(file_names), chunk_size)]
    chunk_data = [(chunk, timeout) for chunk in chunks]  # Include timeout in chunk data

    # Set up multiprocessing pool
    with Pool(processes=cpu_count) as pool:
        results = []
        with tqdm(total=len(file_names), desc="Processing files") as pbar:
            for chunk_result in pool.imap_unordered(process_chunk, chunk_data):
                results.extend(chunk_result)
                pbar.update(len(chunk_result))

    # Update DataFrame with results
    result_dict = {file_name: (polarity, centroided) for file_name, polarity, centroided in results}
    df['Polarity'] = df['Filename'].map(lambda x: result_dict.get(x, ('unknown', 'unknown'))[0])
    df['Centroided'] = df['Filename'].map(lambda x: result_dict.get(x, ('unknown', 'unknown'))[1])

    df.to_csv(f'{repo}_results_with_polarity_centroided.tsv', sep='\t', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Get file polarity and centroid information with multiprocessing')
    parser.add_argument('--repo', '-r', type=str, default='msv', help='Repository name, wmb, mtbls, or msv')
    parser.add_argument('--cpu', '-c', type=int, default=None, help='Number of CPU cores to use')
    parser.add_argument('--chunk-size', '-s', type=int, default=1000, help='Chunk size for processing')
    parser.add_argument('--timeout', '-t', type=int, default=30, help='Timeout for processing each file (in seconds)')
    args = parser.parse_args()

    main(args.repo, args.cpu, args.chunk_size, args.timeout)