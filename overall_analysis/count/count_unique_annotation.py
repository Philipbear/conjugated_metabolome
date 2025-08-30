"""
Count annotated spectra and unique annotations using parallel processing
"""

import pandas as pd
import os
import multiprocessing as mp
from tqdm import tqdm
from functools import partial
import pickle


def process_file(file):
    """Process a single TSV file and return unique qry_ids and annotations."""

    df = pd.read_csv(file, sep='\t', low_memory=False)
    
    # sort by ref_2_score then ref_1_score
    df = df.sort_values(by=['ref_2_score', 'ref_1_score'], ascending=False).reset_index(drop=True)
    # Drop duplicate qry_id
    df = df.drop_duplicates(subset=['qry_id'], keep='first').reset_index(drop=True)

    # Create unique annotations
    df['delta_mass_round2'] = df['delta_mass'].apply(lambda x: str(round(x, 2)))
    
    df_1 = df[df['ref_2_score'] > 0.1].reset_index(drop=True)
    df_2 = df[df['ref_2_score'] <= 0.1].reset_index(drop=True)
    
    df_1['annotations'] = df_1.apply(lambda x: f"{x['ref_1_inchi']}_{x['ref_2_inchi']}", axis=1)
    df_2['annotations'] = df_2.apply(lambda x: f"{x['ref_1_inchi']}_{x['delta_mass_round2']}", axis=1)
    
    spec_spec_annotations = df_1['annotations'].value_counts().to_dict()
    spec_delta_annotations = df_2['annotations'].value_counts().to_dict()

    return spec_spec_annotations, spec_delta_annotations

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
    all_spec_spec_annotations = {}
    all_spec_delta_annotations = {}

    # Use partial function for processing
    process_func = partial(process_file)

    # Process files in parallel with progress bar
    for result in tqdm(pool.imap_unordered(process_func, files),
                        total=len(files),
                        desc="Processing files"):
        spec_spec_annotations, spec_delta_annotations = result
        
        # Merge spec-spec annotations
        for annotation, count in spec_spec_annotations.items():
            if annotation in all_spec_spec_annotations:
                all_spec_spec_annotations[annotation] += count
            else:
                all_spec_spec_annotations[annotation] = count
        
        # Merge spec-delta annotations
        for annotation, count in spec_delta_annotations.items():
            if annotation in all_spec_delta_annotations:
                all_spec_delta_annotations[annotation] += count
            else:
                all_spec_delta_annotations[annotation] = count

    # Clean up
    pool.close()
    pool.join()
    
    # Save results    
    results = {'spec_spec': all_spec_spec_annotations, 'spec_delta': all_spec_delta_annotations}
    out_name = 'pos' if 'pos' in input_dir else 'neg'
    out_name = f'{out_name}_unique_annotation.pkl'
    out_path = f'overall_analysis/count/data/{out_name}'
    with open(out_path, 'wb') as f:
        pickle.dump(results, f)
        
    return results


def merge(pos, neg):
    # merge pos and neg
    all_spec_spec_annotations = pos['spec_spec']
    all_spec_delta_annotations = pos['spec_delta']
    
    for annotation, count in neg['spec_spec'].items():
        if annotation in all_spec_spec_annotations:
            all_spec_spec_annotations[annotation] += count
        else:
            all_spec_spec_annotations[annotation] = count
    
    for annotation, count in neg['spec_delta'].items():
        if annotation in all_spec_delta_annotations:
            all_spec_delta_annotations[annotation] += count
        else:
            all_spec_delta_annotations[annotation] = count
    
    # Save merged results
    results = {'spec_spec': all_spec_spec_annotations, 'spec_delta': all_spec_delta_annotations}
    out_path = 'overall_analysis/count/data/all_unique_annotation.pkl'
    with open(out_path, 'wb') as f:
        pickle.dump(results, f)


if __name__ == '__main__':

    # Process positive mode data
    print("Processing positive mode data...")
    pos = main('analysis/data/pos')
    
    print("\nProcessing negative mode data...")
    neg = main('analysis/data/neg')

    merged_results = merge(pos, neg)
    print("\nProcessing complete. Merged results saved.")
