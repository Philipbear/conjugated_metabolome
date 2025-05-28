import pandas as pd
from msbuddy import Msbuddy
from tqdm import tqdm
import multiprocessing as mp
from functools import partial


def process_chunk(chunk, engine, dm_df, mode='pos'):
    """
    Process a chunk of rows with MSBuddy and delta mass annotations
    """
    for idx in chunk.index:
        row = chunk.loc[idx]

        # Skip if ref_2_id already exists
        if pd.notnull(row['ref_2_id']):
            continue

        # Calculate target mass based on mode
        if mode == 'pos':
            delta_mass = row['qry_mz'] - row['ref_1_mass'] - 1.007276466812
        else:
            delta_mass = row['qry_mz'] - row['ref_1_mass'] + 1.007276466812

        target_mass = delta_mass + 18.010565
        mz_tol = max(0.004, 20 * target_mass / 1e6)

        # MSBuddy processing
        formula_list = engine.mass_to_formula(target_mass, mz_tol, False, False, 0.)

        if not formula_list:
            continue

        # Get the first formula from the list
        chunk.loc[idx, 'ref_2_formula'] = formula_list[0].formula
        chunk.loc[idx, 'ref_2_mass'] = target_mass + formula_list[0].mass_error

        # Search in delta mass df
        sub_db = dm_df[dm_df['formula'] == formula_list[0].formula]

        if not sub_db.empty:
            matched_row = sub_db.iloc[0]
            chunk.loc[idx, 'ref_2_name'] = '; '.join(matched_row['name'])

    return chunk


def main(tsv_path, mode='pos', chunk_size=50000, n_processes=None):
    """
    Add metadata, delta mass annotations from ms2db & lipidmaps, and msbuddy annotations
    using parallel processing
    """
    if n_processes is None:
        n_processes = max(1, mp.cpu_count() - 1)  # Leave one CPU free

    print(f"Using {n_processes} processes")

    # Read data
    df = pd.read_csv(tsv_path, sep='\t', low_memory=False)

    ################
    # Load metadata
    (db_id_to_name, db_id_to_formula, db_id_to_mass, db_id_to_inchi14,
     db_id_to_cls_superclass, db_id_to_cls_class, db_id_to_cls_subclass,
     db_id_to_np_superclass, db_id_to_np_class, db_id_to_np_pathway,
     db_id_to_cooh, db_id_to_free_oh, db_id_to_nh2, db_id_to_nh) = load_ms2db_metadata()

    # Fill in metadata
    print("Adding metadata...")
    metadata_mappings = {
        'ref_1_name': db_id_to_name,
        'ref_1_formula': db_id_to_formula,
        'ref_1_mass': db_id_to_mass,
        'ref_1_cls_superclass': db_id_to_cls_superclass,
        'ref_1_cls_class': db_id_to_cls_class,
        'ref_1_cls_subclass': db_id_to_cls_subclass,
        'ref_1_np_superclass': db_id_to_np_superclass,
        'ref_1_np_class': db_id_to_np_class,
        'ref_1_np_pathway': db_id_to_np_pathway,
        'ref_2_inchi': db_id_to_inchi14,
        'ref_2_name': db_id_to_name,
        'ref_2_formula': db_id_to_formula,
        'ref_2_mass': db_id_to_mass,
        'ref_2_cls_superclass': db_id_to_cls_superclass,
        'ref_2_cls_class': db_id_to_cls_class,
        'ref_2_cls_subclass': db_id_to_cls_subclass,
        'ref_2_np_superclass': db_id_to_np_superclass,
        'ref_2_np_class': db_id_to_np_class,
        'ref_2_np_pathway': db_id_to_np_pathway
    }

    for col, mapping in metadata_mappings.items():
        ref_col = 'ref_1_id' if col.startswith('ref_1_') else 'ref_2_id'
        df[col] = df[ref_col].map(mapping)

    # Calculate delta mass
    if mode == 'pos':
        df['delta_mass'] = df['qry_mz'] - df['ref_1_mass'] - 1.007276466812
    else:
        df['delta_mass'] = df['qry_mz'] - df['ref_1_mass'] + 1.007276466812

    # Filter based on delta mass
    df = df[df['delta_mass'] > 1.5].reset_index(drop=True)

    # Load delta mass data
    dm_df = load_delta_mass()

    # Initialize MSBuddy engine for each process
    engine = Msbuddy()

    # Split dataframe into chunks
    chunks = [df[i:i + chunk_size] for i in range(0, len(df), chunk_size)]
    print(f"Processing {len(chunks)} chunks of size {chunk_size}...")

    # Create partial function with fixed arguments
    process_chunk_partial = partial(process_chunk, engine=engine, dm_df=dm_df, mode=mode)

    # Process chunks in parallel with a single progress bar
    processed_chunks = []
    with mp.Pool(processes=n_processes) as pool:
        for chunk_result in tqdm(
                pool.imap_unordered(process_chunk_partial, chunks),
                total=len(chunks),
                desc="Processing chunks",
                unit="chunk"
        ):
            processed_chunks.append(chunk_result)

    # Combine processed chunks
    result_df = pd.concat(processed_chunks, axis=0)

    # Save results
    output_tsv = tsv_path.replace('raw.tsv', 'annotated.tsv')
    result_df.to_csv(output_tsv, sep='\t', index=False)
    print(f"Results saved to {output_tsv}")


def load_ms2db_metadata(
        # path='/Users/shipei/Documents/projects/chemical_conjugate_discovery/db/ms2db/all/all_ms2db_metadata.pkl'
        path='/home/shipei/projects/revcos/search/all_ms2db_metadata.pkl'
):
    print("Loading MS2DB metadata...")

    metadata = pd.read_pickle(path)

    db_id_to_name = metadata.set_index('db_id')['name'].to_dict()
    db_id_to_formula = metadata.set_index('db_id')['formula'].to_dict()
    db_id_to_mass = metadata.set_index('db_id')['monoisotopic_mass'].to_dict()
    db_id_to_inchi14 = metadata.set_index('db_id')['inchikey_14'].to_dict()
    db_id_to_cls_superclass = metadata.set_index('db_id')['classyfire_superclass'].to_dict()
    db_id_to_cls_class = metadata.set_index('db_id')['classyfire_class'].to_dict()
    db_id_to_cls_subclass = metadata.set_index('db_id')['classyfire_subclass'].to_dict()
    db_id_to_np_superclass = metadata.set_index('db_id')['npclassifier_superclass_results'].to_dict()
    db_id_to_np_class = metadata.set_index('db_id')['npclassifier_class_results'].to_dict()
    db_id_to_np_pathway = metadata.set_index('db_id')['npclassifier_pathway_results'].to_dict()
    db_id_to_cooh = metadata.set_index('db_id')['cooh'].to_dict()
    db_id_to_free_oh = metadata.set_index('db_id')['free_oh'].to_dict()
    db_id_to_nh2 = metadata.set_index('db_id')['primary_nh2'].to_dict()
    db_id_to_nh = metadata.set_index('db_id')['secondary_nh'].to_dict()

    return (db_id_to_name, db_id_to_formula, db_id_to_mass, db_id_to_inchi14,
            db_id_to_cls_superclass, db_id_to_cls_class, db_id_to_cls_subclass,
            db_id_to_np_superclass, db_id_to_np_class, db_id_to_np_pathway,
            db_id_to_cooh, db_id_to_free_oh, db_id_to_nh2, db_id_to_nh)


def load_delta_mass(
        # path='/Users/shipei/Documents/projects/chemical_conjugate_discovery/analysis/data/delta_mass.pkl'
        path='/home/shipei/projects/revcos/search/delta_mass.pkl'
):
    """
    'formula', 'exact_mass', 'unique_inchikeys', 'inchikey_14', 'name'
    """
    print("Loading delta mass data...")
    delta = pd.read_pickle(path)
    return delta


def filter(tsv_path):
    """
    Filter the annotated results to keep only the relevant columns
    """
    df = pd.read_csv(tsv_path, sep='\t', low_memory=False)

    df = df[(df['ref_1_prec_frag_int'] >= 0.05) | (df['ref_1_prec_frag_water_loss_int'] >= 0.05)].reset_index(drop=True)

    df.to_csv(tsv_path.replace('.tsv', '_final.tsv'), sep='\t', index=False)
    # df.to_pickle(tsv_path.replace('.tsv', '_final.pkl'))
    print(f"Filtered results saved")


if __name__ == "__main__":

    # main('/home/shipei/projects/revcos/search/results/neg_stitched/neg_best_raw.tsv',
    #      mode='neg')
    filter('/home/shipei/projects/revcos/search/results/neg_stitched/neg_best_annotated.tsv')
    
    # main('/home/shipei/projects/revcos/search/results/pos_stitched/pos_best_raw.tsv',
    #      mode='pos')
    filter('/home/shipei/projects/revcos/search/results/pos_stitched/pos_best_annotated.tsv')