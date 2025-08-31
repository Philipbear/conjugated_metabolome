
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

        # sort by ref_score
        df = df.sort_values(by=['ref_1_score', 'ref_2_score'], ascending=False).reset_index(drop=True)
        # Drop duplicate qry_id
        df = df.drop_duplicates(subset=['qry_id'], keep='first').reset_index(drop=True)
    
        # delta
        df = df[df['ref_2_inchi'].isna()].reset_index(drop=True)
        
        mz_ls = df['delta_mass'].tolist()
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

    main('analysis/data/pos', 'overall_analysis/delta_mass/data/pos_delta_mass.pkl')

    main('analysis/data/neg', 'overall_analysis/delta_mass/data/neg_delta_mass.pkl')
