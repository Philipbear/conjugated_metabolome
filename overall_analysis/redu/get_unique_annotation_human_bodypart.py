import os
import pandas as pd
import multiprocessing as mp
from tqdm import tqdm
import pickle
from functools import partial


def process_file(file, redu):
    """Process a single TSV file and return unique qry_ids and annotations."""

    df = pd.read_csv(file, sep='\t', low_memory=False)

    # sort by ref_1_score then ref_2_score
    df = df.sort_values(by=['ref_1_score', 'ref_2_score'], ascending=False).reset_index(drop=True)
    # Drop duplicate qry_id
    df = df.drop_duplicates(subset=['qry_id'], keep='first').reset_index(drop=True)
    
    df['mri'] = df['qry_id'].apply(lambda x: x.split(':scan')[0])
    # merge with redu
    df = df.merge(redu, on='mri', how='left')
    # Filter for human samples only
    df = df[df['NCBITaxonomy'].str.contains('9606|Homo sapiens', na=False)]
    
    # Filter out missing values from UBERONBodyPartName
    df = df[df['UBERONBodyPartName'] != 'missing value']
    df = df.dropna(subset=['UBERONBodyPartName'])

    # Create unique annotations - combine all without splitting by ref_2_score
    df['delta_mass_round2'] = df['delta_mass'].apply(lambda x: str(round(x, 2)))
    
    # Create annotation
    def create_annotation(row):
        if pd.notnull(row['ref_2_inchi']):
            return f"{row['ref_1_inchi']}_{row['ref_2_inchi']}"
        else:
            return f"{row['ref_1_inchi']}_{row['delta_mass_round2']}"
    
    df['annotation'] = df.apply(create_annotation, axis=1)
    
    # group by annotations and UBERONBodyPartName, count unique annotations
    df_combined = df.groupby(['annotation', 'UBERONBodyPartName']).size().reset_index(name='count')

    return df_combined


def process_ion_mode(input_dir, n_processes=10):
    
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
    
    redu = load_redu()

    # Process files in parallel with progress bar
    all_annotations_df = pd.DataFrame()

    # Use partial function for processing
    process_func = partial(process_file, redu=redu)

    # Process files in parallel with progress bar
    for result in tqdm(pool.imap_unordered(process_func, files),
                        total=len(files),
                        desc="Processing files"):
        annotations_df = result

        all_annotations_df = pd.concat([all_annotations_df, annotations_df], ignore_index=True)
    
    # group by annotations and UBERONBodyPartName, sum count
    all_annotations_df = all_annotations_df.groupby(['annotation', 'UBERONBodyPartName']).sum().reset_index()

    # Clean up
    pool.close()
    pool.join()

    return all_annotations_df


def merge(pos, neg):
    # merge pos and neg
    all_annotations_df = pd.concat([pos, neg], ignore_index=True)
    
    # Group by annotation and body part again to combine pos/neg counts
    all_annotations_df = all_annotations_df.groupby(['annotation', 'UBERONBodyPartName']).sum().reset_index()

    return all_annotations_df
    

def get_raw_data_main():

    # Process positive mode data
    print("Processing positive mode data...")
    pos = process_ion_mode('analysis/data/pos')
    
    print("\nProcessing negative mode data...")
    neg = process_ion_mode('analysis/data/neg')

    merged_results = merge(pos, neg)
    
    # Save results - now just one combined dataset
    out_name = 'all_unique_annotation_human_bodypart.pkl'
    out_path = f'overall_analysis/redu/data/{out_name}'
    with open(out_path, 'wb') as f:
        pickle.dump(merged_results, f)


def load_redu():
    print('Loading REDU data...')
    redu = pd.read_csv('overall_analysis/redu/data/redu.tsv', sep='\t', low_memory=False)
    return redu

    
if __name__ == '__main__':
    
    get_raw_data_main()