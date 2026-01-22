"""
Count annotated spectra and unique annotations using parallel processing
"""

import pandas as pd
import os
import multiprocessing as mp
from tqdm import tqdm
from functools import partial


def process_file(file):
    """Process a single TSV file and return unique qry_ids and annotations."""

    df = pd.read_csv(file, sep='\t', low_memory=False)

    # Get unique query IDs
    qry_ids = set(df['qry_id'].tolist())

    # Create unique annotations
    df['delta_mass_round2'] = df['delta_mass'].apply(lambda x: str(round(x, 2)))
    df['annotations'] = df.apply(lambda x: f"{x['ref_1_inchi']}_{x['delta_mass_round2']}", axis=1)
    annotations = set(df['annotations'].tolist())

    # sort by ref_2_score
    df = df.sort_values(by='ref_2_score', ascending=False).drop_duplicates(subset=['qry_id'], keep='first').reset_index(drop=True)
    
    # count how many ref_2_score > 0.5
    high_score_ids = df[df['ref_2_score'] > 0.5]['qry_id'].tolist()

    return qry_ids, annotations, high_score_ids


def main(input_dir, n_processes=10):
    """
    Process all TSV files in the input directory using parallel processing.

    Args:
        input_dir: Directory containing TSV files
        n_processes: Number of processes to use (defaults to CPU count)
    """
    # Set default number of processes
    if n_processes is None:
        n_processes = max(1, mp.cpu_count() - 1)  # Leave one CPU free

    print(f"Using {n_processes} processes")

    # Get list of TSV files
    files = [os.path.join(input_dir, x) for x in os.listdir(input_dir)
             if x.endswith('.tsv') and not x.startswith('.')]

    print(f"Found {len(files)} TSV files to process")

    # Create process pool
    pool = mp.Pool(processes=n_processes)

    # Process files in parallel with progress bar
    all_qry_ids = set()
    all_annotations = set()
    ref_2_anno_ids = set()

    # Use partial function for processing
    process_func = partial(process_file)

    # Process files in parallel with progress bar
    results = []
    for result in tqdm(pool.imap_unordered(process_func, files),
                        total=len(files),
                        desc="Processing files"):
        qry_ids, annotations, anno_cnt = result
        all_qry_ids.update(qry_ids)
        all_annotations.update(annotations)
        ref_2_anno_ids.update(anno_cnt)
        # Update progress stats
        results.append((len(qry_ids), len(annotations)))

    # Print results
    print(f"Total number of spectra: {len(all_qry_ids)}")
    print(f"Total number of unique annotations: {len(all_annotations)}")
    print(f"Total ref_2 annotated: {len(ref_2_anno_ids)}")

    # Clean up
    pool.close()
    pool.join()

    return

    
if __name__ == '__main__':

    main('analysis/data/pos')
    """
    Total number of spectra: 22431333
    Total number of unique annotations: 9375256
    Total ref_2 annotated: 3448465
    """

    main('analysis/data/neg')
    """
    Total number of spectra: 1796106
    Total number of unique annotations: 766298
    Total ref_2 annotated: 97878
    """
    
    '''
    22431333 + 1796106 = 24227439
    '''
    