"""
'qry_id', 'qry_mz',
'ref_1_id', 'ref_1_score', 'ref_1_peak', 'ref_1_usage', 'ref_1_prec_int',
'ref_2_id', 'ref_2_score', 'ref_2_peak', 'ref_2_usage'
'ref_1_inchi', 'ref_1_name', 'ref_1_formula', 'ref_1_mass', 'ref_1_cls_superclass', 'ref_1_cls_class', 'ref_1_cls_subclass', 'ref_1_np_superclass',
'ref_1_np_class', 'ref_1_np_pathway',
'ref_2_inchi', 'ref_2_name', 'ref_2_formula', 'ref_2_mass', 'ref_2_cls_superclass', 'ref_2_cls_class', 'ref_2_cls_subclass', 'ref_2_np_superclass',
'ref_2_np_class', 'ref_2_np_pathway',
'delta_mass'
"""

import pandas as pd
import numpy as np
from tqdm import tqdm
import glob
from multiprocessing import Pool
import os

# Create the delta mass reference array (10 to 300, step 0.01)
DELTA_MASS_ARRAY = np.arange(10, 300.01, 0.01)
DELTA_MASS_MIN = 10.0
DELTA_MASS_MAX = 300.0
DELTA_MASS_STEP = 0.01

def delta_mass_to_index(delta_mass):
    """Convert delta mass to array index"""
    if delta_mass < DELTA_MASS_MIN or delta_mass > DELTA_MASS_MAX:
        return None
    return int(round((delta_mass - DELTA_MASS_MIN) / DELTA_MASS_STEP))

def index_to_delta_mass(index):
    """Convert array index back to delta mass"""
    return DELTA_MASS_MIN + index * DELTA_MASS_STEP


def get_all_results():
    """
    Process all chemicals by processing each file once and writing results for each compound.
    """

    # Create output directory if it doesn't exist
    output_folder = 'data'
    os.makedirs(output_folder, exist_ok=True)

    # Process files
    print("Processing files...")
    pos_results = process_mode_files(mode='pos')
    neg_results = process_mode_files(mode='neg')

    # Combine results
    all_results = []
    if pos_results is not None:
        all_results.append(pos_results)
    if neg_results is not None:
        all_results.append(neg_results)
    
    if all_results:
        # Concatenate all results and group again
        final_df = pd.concat(all_results, ignore_index=True)
        final_grouped = final_df.groupby(['ref_1_inchi', 'delta_mass_idx'])['count'].sum().reset_index()
        
        # Save final results
        final_grouped.to_csv(os.path.join(output_folder, 'combined_results.tsv'), sep='\t', index=False)
        print(f"Final results saved with {len(final_grouped)} unique combinations")


def process_mode_files(mode):
    """
    Process all files for a given mode in parallel.
    """
    input_folder = f'/Users/shipei/Documents/projects/conjugated_metabolome/analysis/data/{mode}'

    # Get all TSV files in the input folder
    all_files = glob.glob(os.path.join(input_folder, '*.tsv'))

    if not all_files:
        print(f"No TSV files found in {input_folder}")
        return None

    print(f"Found {len(all_files)} files to process in {input_folder}")

    # Process files in parallel
    num_workers = min(11, len(all_files))
    print(f"Using {num_workers} workers for parallel processing")

    with Pool(num_workers) as pool:
        # Create a list of arguments for each file
        args = [(file_path,) for file_path in all_files]

        # Process files in parallel with progress bar
        results = list(tqdm(pool.imap(apply_process_single_file, args), total=len(args)))

    # Filter out None results and concatenate
    valid_results = [result for result in results if result is not None and not result.empty]
    
    if not valid_results:
        print(f"No valid results found for mode {mode}")
        return None
    
    # Concatenate all results from this mode
    mode_df = pd.concat(valid_results, ignore_index=True)
    
    # Group by ref_1_inchi and delta_mass_idx again to combine counts from all files
    mode_grouped = mode_df.groupby(['ref_1_inchi', 'delta_mass_idx'])['count'].sum().reset_index()
    
    print(f"Mode {mode}: {len(mode_grouped)} unique combinations")
    return mode_grouped


def apply_process_single_file(args):
    return process_single_file(*args)


def process_single_file(file_path):
    """
    Process a single TSV file for all target compounds and write individual results.
    """
    try:
        # Read the file
        df = pd.read_csv(file_path, sep='\t', low_memory=False)

        # valid ref_1_inchi
        df = df[df['ref_1_inchi'].notna() & (df['ref_1_inchi'] != '')].reset_index(drop=True)
        
        if df.empty:
            return pd.DataFrame()
        
        # Convert delta_mass to index
        df['delta_mass_idx'] = df['delta_mass'].apply(delta_mass_to_index)
        
        # Remove rows where delta_mass is out of range
        df = df[df['delta_mass_idx'].notna()].reset_index(drop=True)
        df['delta_mass_idx'] = df['delta_mass_idx'].astype(int)

        if df.empty:
            return pd.DataFrame()

        # group by ref_1_inchi and delta_mass_idx, add count column
        grouped_df = df.groupby(['ref_1_inchi', 'delta_mass_idx']).size().reset_index(name='count')

        return grouped_df
    
    except Exception as e:
        print(f"Error processing file {file_path}: {e}")
        return pd.DataFrame()


def summarize_results():
    """
    Prepare data for UMAP plots by creating vectors for each unique ref_1_inchi.
    Each vector indicates the count for each unique delta_mass index.
    """
    output_folder = 'data'
    combined_file = os.path.join(output_folder, 'combined_results.tsv')
    
    if not os.path.exists(combined_file):
        print(f"Combined results file not found: {combined_file}")
        return
    
    # Read the combined results
    df = pd.read_csv(combined_file, sep='\t')
    print(f"Loaded {len(df)} combinations from combined results")
    
    # Get unique compounds
    unique_inchi = df['ref_1_inchi'].unique()
    print(f"Found {len(unique_inchi)} unique compounds")
    
    # Create feature matrix: rows = compounds, columns = delta_mass indices
    num_compounds = len(unique_inchi)
    num_features = len(DELTA_MASS_ARRAY)
    
    feature_matrix = np.zeros((num_compounds, num_features), dtype=int)
    
    # Create mapping from inchi to row index
    inchi_to_row = {inchi: i for i, inchi in enumerate(unique_inchi)}
    
    # Fill the feature matrix
    for _, row in df.iterrows():
        inchi = row['ref_1_inchi']
        delta_idx = int(row['delta_mass_idx'])
        count = row['count']
        
        row_idx = inchi_to_row[inchi]
        if 0 <= delta_idx < num_features:
            feature_matrix[row_idx, delta_idx] += count
    
    print(f"Created feature matrix with shape {feature_matrix.shape}")
    
    # Save the feature matrix
    feature_matrix_file = os.path.join(output_folder, 'umap_feature_matrix.npy')
    np.save(feature_matrix_file, feature_matrix)
    print(f"Feature matrix saved to {feature_matrix_file}")
    
    # Save the compound identifiers
    compound_ids_file = os.path.join(output_folder, 'compound_ids.txt')
    with open(compound_ids_file, 'w') as f:
        for inchi in unique_inchi:
            f.write(f"{inchi}\n")
    print(f"Compound IDs saved to {compound_ids_file}")
    
    # Save delta mass array for reference
    delta_mass_file = os.path.join(output_folder, 'delta_mass_array.npy')
    np.save(delta_mass_file, DELTA_MASS_ARRAY)
    print(f"Delta mass array saved to {delta_mass_file}")
    
    # Save metadata
    metadata = {
        'delta_mass_min': DELTA_MASS_MIN,
        'delta_mass_max': DELTA_MASS_MAX,
        'delta_mass_step': DELTA_MASS_STEP,
        'num_compounds': num_compounds,
        'num_features': num_features
    }
    
    metadata_file = os.path.join(output_folder, 'metadata.txt')
    with open(metadata_file, 'w') as f:
        for key, value in metadata.items():
            f.write(f"{key}: {value}\n")
    print(f"Metadata saved to {metadata_file}")
    
    # Print summary statistics
    print(f"\nSummary:")
    print(f"- Number of compounds: {num_compounds}")
    print(f"- Number of delta_mass features: {num_features}")
    print(f"- Delta mass range: {DELTA_MASS_MIN} to {DELTA_MASS_MAX} (step {DELTA_MASS_STEP})")
    print(f"- Matrix shape: {feature_matrix.shape}")
    print(f"- Total non-zero entries: {np.count_nonzero(feature_matrix)}")
    print(f"- Sparsity: {1 - np.count_nonzero(feature_matrix) / feature_matrix.size:.2%}")
    print(f"- Max count in matrix: {feature_matrix.max()}")
    print(f"- Total counts: {feature_matrix.sum()}")
    
    return feature_matrix, unique_inchi, DELTA_MASS_ARRAY


if __name__ == "__main__":
    
    import os
    os.chdir(os.path.dirname(__file__))
    
    get_all_results()
    feature_matrix, compound_ids, delta_mass_array = summarize_results()