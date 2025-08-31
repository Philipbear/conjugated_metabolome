"""
Count annotated spectra and unique annotations using parallel processing
"""

from matplotlib.pyplot import ion
import pandas as pd
import os
import multiprocessing as mp
from tqdm import tqdm
from functools import partial
import pickle


def process_file(file):
    """Process a single TSV file and return unique qry_ids and annotations."""

    df = pd.read_csv(file, sep='\t', low_memory=False)
    
    # sort by ref_score
    df = df.sort_values(by=['ref_1_score', 'ref_2_score'], ascending=False).reset_index(drop=True)
    # Drop duplicate qry_id
    df = df.drop_duplicates(subset=['qry_id'], keep='first').reset_index(drop=True)

    df['dataset'] = df['qry_id'].apply(lambda x: x.split(':')[0])

    # Create unique annotations
    df['delta_mass_round2'] = df['delta_mass'].apply(lambda x: str(round(x, 2)))

    df_1 = df[df['ref_2_inchi'].notna()].reset_index(drop=True)
    df_2 = df[df['ref_2_inchi'].isna()].reset_index(drop=True)

    def get_annotation(row):
        # sort ref_1_inchi and ref_2_inchi alphabetically
        inchi_1 = row['ref_1_inchi']
        inchi_2 = row['ref_2_inchi']
        try:
            if inchi_1 > inchi_2:
                inchi_1, inchi_2 = inchi_2, inchi_1
        except Exception as e:
            print(f"Error processing row {row.name}: {e}")
            print(f"Values: ref_1_inchi={inchi_1}, ref_2_inchi={inchi_2}")
        return f"{inchi_1}_{inchi_2}"

    df_1['annotations'] = df_1.apply(get_annotation, axis=1)
    df_2['annotations'] = df_2.apply(lambda x: f"{x['ref_1_inchi']}_{x['delta_mass_round2']}", axis=1)

    df_1 = df_1[['dataset','annotations']]
    df_2 = df_2[['dataset', 'annotations']]
    
    df_1 = df_1.drop_duplicates(subset=['dataset', 'annotations'], keep='first').reset_index(drop=True)
    df_2 = df_2.drop_duplicates(subset=['dataset', 'annotations'], keep='first').reset_index(drop=True)

    return df_1, df_2


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
    all_spec_spec_annotations = pd.DataFrame(columns=['dataset', 'annotations'])
    all_spec_delta_annotations = pd.DataFrame(columns=['dataset', 'annotations'])

    # Use partial function for processing
    process_func = partial(process_file)

    # Process files in parallel with progress bar
    for result in tqdm(pool.imap_unordered(process_func, files),
                        total=len(files),
                        desc="Processing files"):
        spec_spec_annotations, spec_delta_annotations = result
        
        # Merge
        all_spec_spec_annotations = pd.concat([all_spec_spec_annotations, spec_spec_annotations], ignore_index=True)
        all_spec_delta_annotations = pd.concat([all_spec_delta_annotations, spec_delta_annotations], ignore_index=True)

    # Clean up
    pool.close()
    pool.join()
    
    # dereplicate
    all_spec_spec_annotations = all_spec_spec_annotations.drop_duplicates(subset=['dataset', 'annotations'], keep='first').reset_index(drop=True)
    all_spec_delta_annotations = all_spec_delta_annotations.drop_duplicates(subset=['dataset', 'annotations'], keep='first').reset_index(drop=True)

    # save
    ion_mode = 'pos' if 'pos' in input_dir else 'neg'
    out_path_1 = f'overall_analysis/count/data/all_unique_spec_annotation_{ion_mode}.pkl'
    all_spec_spec_annotations.to_pickle(out_path_1)
    out_path_2 = f'overall_analysis/count/data/all_unique_delta_annotation_{ion_mode}.pkl'
    all_spec_delta_annotations.to_pickle(out_path_2)


def merge_results():
    pos = pd.read_pickle('overall_analysis/count/data/all_unique_spec_annotation_pos.pkl')
    neg = pd.read_pickle('overall_analysis/count/data/all_unique_spec_annotation_neg.pkl')
    merged = pd.concat([pos, neg], ignore_index=True)

    # save
    out_path = 'overall_analysis/count/data/all_unique_spec_annotation_merged.pkl'
    merged.to_pickle(out_path)

    pos = pd.read_pickle('overall_analysis/count/data/all_unique_delta_annotation_pos.pkl')
    neg = pd.read_pickle('overall_analysis/count/data/all_unique_delta_annotation_neg.pkl')
    merged = pd.concat([pos, neg], ignore_index=True)

    # save
    out_path = 'overall_analysis/count/data/all_unique_delta_annotation_merged.pkl'
    merged.to_pickle(out_path)


if __name__ == '__main__':

    # Process positive mode data
    print("Processing positive mode data...")
    main('analysis/data/pos')
    
    print("\nProcessing negative mode data...")
    main('analysis/data/neg')

    merge_results()
    print("\nProcessing complete. Merged results saved.")
