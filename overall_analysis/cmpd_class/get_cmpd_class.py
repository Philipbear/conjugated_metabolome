import pandas as pd
import os
import multiprocessing as mp
from tqdm import tqdm
from functools import partial
import pickle
import numpy as np

# Define the class order for natural product pathways
CLASS_ORDER = [
    'Fatty acids',
    'Shikimates and Phenylpropanoids',
    'Terpenoids',
    'Alkaloids',
    'Amino acids and Peptides',
    'Carbohydrates',
    'Polyketides',
    'Others'
]


def get_primary_pathway(pathway_str):
    """Extract the primary pathway from a comma-separated list."""
    if pd.isna(pathway_str):
        return 'Others'

    # Split by comma and get the first pathway
    pathways = [p.strip() for p in pathway_str.split(',')]
    primary = pathways[0]

    # Check if primary pathway is in our defined classes
    if primary in CLASS_ORDER:
        return primary
    else:
        return 'Others'


def process_file(file):
    """Process a single TSV file and return categorized data."""
    df = pd.read_csv(file, sep='\t', low_memory=False)

    # ref_2_id not null
    df = df[df['ref_2_id'].notnull()].reset_index(drop=True)

    df = df[df['ref_1_np_pathway'].notnull()].reset_index(drop=True)
    df = df[df['ref_2_np_pathway'].notnull()].reset_index(drop=True)

    # dereplicate by qry_id
    df = df.drop_duplicates(subset=['qry_id']).reset_index(drop=True)

    # Extract primary pathways
    df['class_1'] = df['ref_1_np_pathway'].apply(get_primary_pathway)
    df['class_2'] = df['ref_2_np_pathway'].apply(get_primary_pathway)

    # Create combined class name
    df['class_class'] = df['class_1'] + '_' + df['class_2']

    # Count occurrences of each combination
    result = df.groupby(['class_1', 'class_2']).size().reset_index(name='count')

    return result


def main(input_dir, out_name, n_processes=8):
    """
    Process all TSV files in the input directory using parallel processing.

    Args:
        input_dir: Directory containing TSV files
        out_name: Output filename for the pickle file
        n_processes: Number of processes to use (defaults to 8)
    """

    # Get list of TSV files
    files = [os.path.join(input_dir, x) for x in os.listdir(input_dir)
             if x.endswith('.tsv') and not x.startswith('.')]

    print(f"Found {len(files)} TSV files to process")

    # Create process pool
    pool = mp.Pool(processes=n_processes)

    # Process files in parallel with progress bar
    all_results = []

    try:
        # Use partial function for processing
        process_func = partial(process_file)

        # Process files in parallel with progress bar
        for result in tqdm(pool.imap_unordered(process_func, files),
                           total=len(files),
                           desc="Processing files"):
            all_results.append(result)

    finally:
        # Clean up
        pool.close()
        pool.join()

    # Combine all results
    combined_df = pd.concat(all_results)

    # Sum up counts for the same combinations
    summary_df = combined_df.groupby(['class_1', 'class_2'])['count'].sum().reset_index()

    # Create a cross-tabulation of the data
    cross_tab = pd.pivot_table(
        summary_df,
        values='count',
        index='class_1',
        columns='class_2',
        fill_value=0
    )

    # Ensure all classes are present in both rows and columns
    for cls in CLASS_ORDER:
        if cls not in cross_tab.index:
            cross_tab.loc[cls] = 0
        if cls not in cross_tab.columns:
            cross_tab[cls] = 0

    # Reorder rows and columns according to CLASS_ORDER
    cross_tab = cross_tab.loc[CLASS_ORDER][CLASS_ORDER]

    # Convert to the format shown in the image
    result_df = []

    for i, class_1 in enumerate(CLASS_ORDER):
        for j, class_2 in enumerate(CLASS_ORDER):
            # Skip combinations where j < i to avoid duplicates
            if j < i:
                continue

            count = cross_tab.loc[class_1, class_2]
            if count > 0:
                class_class = f"{class_1}_{class_2}"
                result_df.append({
                    'smaller_class_idx': i,
                    'larger_class_idx': j,
                    'count': int(count),
                    'class_1': class_1,
                    'class_2': class_2,
                    'class_class': class_class
                })

    result_df = pd.DataFrame(result_df)

    # Save both the raw data and the formatted result
    with open(out_name, 'wb') as f:
        pickle.dump({
            'cross_tab': cross_tab,
            'summary_table': result_df
        }, f)

    # Also save as CSV for easy viewing
    result_df.to_csv(f"{out_name.replace('.pkl', '.csv')}", index=False)

    print(f"Saved results to {out_name} and {out_name.replace('.pkl', '.csv')}")
    return result_df


if __name__ == '__main__':
    import os
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    
    # Process positive mode data
    pos_result = main('/Users/shipei/Documents/projects/conjugated_metabolome/analysis/data/pos',
                      'data/pos.pkl')

    # Process negative mode data
    neg_result = main('/Users/shipei/Documents/projects/conjugated_metabolome/analysis/data/neg',
                      'data/neg.pkl')

    # combine the results
    combined_result = pd.concat([pos_result, neg_result])
    combined_result = combined_result.groupby(['smaller_class_idx', 'larger_class_idx', 'class_1', 'class_2'])[
        'count'].sum().reset_index()
    combined_result['class_class'] = combined_result['class_1'] + '_' + combined_result['class_2']
    combined_result.to_csv('data/combined_summary.csv', index=False)

    print("Processing complete. Results saved to data/pos.pkl, data/neg.pkl, and data/combined_summary.csv")