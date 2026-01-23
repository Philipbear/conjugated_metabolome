"""
Count annotated spectra and unique annotations using parallel processing
"""

from httpx import get
import pandas as pd
import os
import multiprocessing as mp
from tqdm import tqdm
from functools import partial


def process_file(file):
    """Process a single TSV file and return unique qry_ids and annotations."""

    df = pd.read_csv(file, sep='\t', low_memory=False)
    
    # sort by ref_score
    df = df.sort_values(by=['ref_1_score', 'ref_2_score'], ascending=False).reset_index(drop=True)
    # Drop duplicate qry_id
    df = df.drop_duplicates(subset=['qry_id'], keep='first').reset_index(drop=True)

    # Create unique annotations
    df['delta_mass_round2'] = df['delta_mass'].apply(lambda x: str(round(x, 2)))

    def add_annotation(row):
        if pd.notnull(row['ref_2_inchi']):
            # sort 2 inchis alphabetically
            inchi_1 = row['ref_1_inchi']
            inchi_2 = row['ref_2_inchi']
            if inchi_1 < inchi_2:
                return f"{inchi_1}_{inchi_2}"
            else:
                return f"{inchi_2}_{inchi_1}"
        else:
            return f"{row['ref_1_inchi']}_{row['delta_mass_round2']}"
    df['annotation'] = df.apply(add_annotation, axis=1)

    # add dataset
    df['dataset'] = df['qry_id'].apply(lambda x: x.split(':')[0])

    # df = df.drop_duplicates(subset=['annotation'], keep='first').reset_index(drop=True)

    return df[['ref_2_inchi', 'annotation', 'dataset']]


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
    all_dfs = pd.DataFrame(columns=['ref_2_inchi', 'annotation', 'dataset'])

    # Use partial function for processing
    process_func = partial(process_file)

    # Process files in parallel with progress bar
    for result in tqdm(pool.imap_unordered(process_func, files),
                        total=len(files),
                        desc="Processing files"):
        all_dfs = pd.concat([all_dfs, result], ignore_index=True)

    # Clean up
    pool.close()
    pool.join()
    
    # save
    ion_mode = 'pos' if 'pos' in input_dir else 'neg'
    out_path = f'overall_analysis/count/data/all_unique_annotation_{ion_mode}.pkl'
    all_dfs.to_pickle(out_path)


def merge_results():
    pos = pd.read_pickle('overall_analysis/count/data/all_unique_annotation_pos.pkl')
    neg = pd.read_pickle('overall_analysis/count/data/all_unique_annotation_neg.pkl')
    merged = pd.concat([pos, neg], ignore_index=True)

    all_spec_spec_annotations = merged[merged['ref_2_inchi'].notnull()].copy()
    all_spec_delta_annotations = merged[merged['ref_2_inchi'].isnull()].copy()

    # dereplicate by annotation and dataset
    all_spec_spec_annotations = all_spec_spec_annotations.drop_duplicates(subset=['annotation', 'dataset'], keep='first').reset_index(drop=True)
    all_spec_delta_annotations = all_spec_delta_annotations.drop_duplicates(subset=['annotation', 'dataset'], keep='first').reset_index(drop=True)
    
    # save
    all_spec_spec_annotations.to_pickle('overall_analysis/count/data/all_unique_spec_annotation_merged.pkl')
    all_spec_delta_annotations.to_pickle('overall_analysis/count/data/all_unique_delta_annotation_merged.pkl')
    
    # further dereplicate by annotation
    all_spec_spec_annotations = all_spec_spec_annotations.drop_duplicates(subset=['annotation'], keep='first').reset_index(drop=True)
    all_spec_delta_annotations = all_spec_delta_annotations.drop_duplicates(subset=['annotation'], keep='first').reset_index(drop=True)
    
    print(f"Total unique spec-spec annotations: {all_spec_spec_annotations.shape[0]}")
    print(f"Total unique spec-delta annotations: {all_spec_delta_annotations.shape[0]}")
    

def get_dataset_count():
    """Get dataset count for each annotation"""
    
    all_spec_spec_annotations = pd.read_pickle('overall_analysis/count/data/all_unique_spec_annotation_merged.pkl')
    # group by annotation and count unique datasets
    spec_spec_counts = all_spec_spec_annotations.groupby('annotation')['dataset'].nunique().reset_index()
    spec_spec_counts = spec_spec_counts.rename(columns={'dataset': 'dataset_count'})
    # sort
    spec_spec_counts = spec_spec_counts.sort_values(by='dataset_count', ascending=False).reset_index(drop=True)
    spec_spec_counts.to_csv('overall_analysis/count/data/spec_spec_annotation_dataset_count.tsv', sep='\t', index=False)
    # summarize the value counts of dataset_count
    value_counts = spec_spec_counts['dataset_count'].value_counts().sort_index()
    # save as tsv
    value_counts.to_csv('overall_analysis/count/data/spec_spec_annotation_dataset_count_summary.tsv', sep='\t', header=['count'])
    
    all_spec_delta_annotations = pd.read_pickle('overall_analysis/count/data/all_unique_delta_annotation_merged.pkl')
    # group by annotation and count unique datasets
    spec_delta_counts = all_spec_delta_annotations.groupby('annotation')['dataset'].nunique().reset_index()
    spec_delta_counts = spec_delta_counts.rename(columns={'dataset': 'dataset_count'})
    # sort
    spec_delta_counts = spec_delta_counts.sort_values(by='dataset_count', ascending=False).reset_index(drop=True)
    spec_delta_counts.to_csv('overall_analysis/count/data/spec_delta_annotation_dataset_count.tsv', sep='\t', index=False)
    # summarize the value counts of dataset_count
    value_counts = spec_delta_counts['dataset_count'].value_counts().sort_index()
    # save as tsv
    value_counts.to_csv('overall_analysis/count/data/spec_delta_annotation_dataset_count_summary.tsv', sep='\t', header=['count'])


if __name__ == '__main__':

    # # Process positive mode data
    # print("Processing positive mode data...")
    # main('analysis/data/pos')
    
    # print("\nProcessing negative mode data...")
    # main('analysis/data/neg')

    # merge_results()
    # print("\nProcessing complete. Merged results saved.")
    
    # """
    # Total unique spec-spec annotations: 217291
    # Total unique spec-delta annotations: 3412720
    # """
    
    get_dataset_count()
