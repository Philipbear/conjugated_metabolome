
import pandas as pd
import os
import multiprocessing as mp
from tqdm import tqdm
from functools import partial
import pickle


def process_file(file):
    """Process a single TSV file and return mz list."""
    try:
        df = pd.read_csv(file, sep='\t', low_memory=False)

        # dereplicate by qry_id
        df = df.drop_duplicates(subset=['qry_id']).reset_index(drop=True)

        mz_ls = df['qry_mz'].tolist()
        return mz_ls

    except Exception as e:
        return []


def main(input_dir, out_name, n_processes=8):
    """
    Process all TSV files in the input directory using parallel processing.

    Args:
        input_dir: Directory containing TSV files
        n_processes: Number of processes to use (defaults to CPU count)
    """

    # Get list of TSV files
    files = [os.path.join(input_dir, x) for x in os.listdir(input_dir)
             if x.endswith('.tsv') and not x.startswith('.')]

    print(f"Found {len(files)} TSV files to process")

    # Create process pool
    pool = mp.Pool(processes=n_processes)

    # Process files in parallel with progress bar
    all_mzs = []

    try:
        # Use partial function for processing
        process_func = partial(process_file)

        # Process files in parallel with progress bar
        for result in tqdm(pool.imap_unordered(process_func, files),
                           total=len(files),
                           desc="Processing files"):
            all_mzs.extend(result)

    finally:
        # Clean up
        pool.close()
        pool.join()

    with open(out_name, 'wb') as f:
        pickle.dump(all_mzs, f)



if __name__ == '__main__':

    main('analysis/data/pos', 'overall_analysis/mz/data/pos_annotated_mz.pkl')

    main('analysis/data/neg', 'overall_analysis/mz/data/neg_annotated_mz.pkl')
