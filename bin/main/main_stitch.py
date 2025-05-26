"""
main stitching code
"""

import pickle
import pandas as pd
import numpy as np
from numba import njit
import os
from tqdm import tqdm
import time
from multiprocessing import Pool, cpu_count


def main(pkl_folder='/home/shipei/projects/revcos/search/results/pos_recheck',
         out_folder='/home/shipei/projects/revcos/search/results/pos_stitched',
         mode='pos', n_jobs=None):

    os.makedirs(out_folder, exist_ok=True)

    # make separate out folders
    # all_out = os.path.join(out_folder, 'all')
    best_out = os.path.join(out_folder, 'best')
    # os.makedirs(all_out, exist_ok=True)
    os.makedirs(best_out, exist_ok=True)

    # Write parameters to file
    write_params_to_file(pkl_folder, out_folder, mode, n_jobs)

    all_pkls = os.listdir(pkl_folder)
    all_pkls = [x for x in all_pkls if x.endswith('.pkl')]

    all_pkls = [os.path.join(pkl_folder, x) for x in all_pkls]

    # Check which files need processing
    files_to_process = []
    for pkl in all_pkls:
        basename = os.path.basename(pkl).split('_results')[0]
        out_path_best = os.path.join(os.path.join(out_folder, 'best'), f"{basename}_stitched_best.tsv")
        # out_path_all = os.path.join(os.path.join(out_folder, 'all'), f"{basename}_stitched_all.pkl")

        # # If either output file doesn't exist, add to processing list
        # if not (os.path.exists(out_path_best) and os.path.exists(out_path_all)):
        #     files_to_process.append(pkl)

        if not (os.path.exists(out_path_best)):
            files_to_process.append(pkl)

    print(f"Found {len(all_pkls)} files, {len(files_to_process)} need processing.")

    if not files_to_process:
        print("All files have already been processed. Skipping to merge step.")
    else:
        # load ms2db metadata
        db_id_to_mass, db_id_to_inchi14 = load_ms2db_metadata()

        # Determine number of processes
        if n_jobs is None:
            n_jobs = min(max(1, cpu_count() - 1), len(files_to_process))

        # Create process pool and process files in parallel
        print(f"Processing {len(files_to_process)} files using {n_jobs} processes...")
        with Pool(n_jobs) as pool:
            # Create argument list for each file
            args = [(pkl, out_folder, mode, db_id_to_mass, db_id_to_inchi14) for pkl in files_to_process]

            # Use imap to process files with progress bar
            for _ in tqdm(pool.imap_unordered(process_file_wrapper, args),
                          total=len(files_to_process),
                          desc="Processing pkl files"):
                pass
    ######
    # merge
    individual_best_files = os.listdir(best_out)
    individual_best_files = [x for x in individual_best_files if x.endswith('.tsv')]
    individual_best_files = [os.path.join(best_out, x) for x in individual_best_files]

    all_best_results = []
    for file in tqdm(individual_best_files, desc='Merging best results'):
        df = pd.read_csv(file, sep='\t')
        all_best_results.append(df)

    all_best_results = pd.concat(all_best_results)
    all_best_results.to_csv(os.path.join(out_folder, f'{mode}_best_raw.tsv'), sep='\t', index=False)


def write_params_to_file(pkl_folder, out_folder, mode, n_jobs):
    """
    Write all parameters to a text file in the output folder.
    """
    params_file = os.path.join(out_folder, 'search_parameters.txt')
    with open(params_file, 'w') as f:
        f.write("Search Parameters\n")
        f.write("================\n\n")
        f.write(f"Input Folder: {pkl_folder}\n")
        f.write(f"Output Folder: {out_folder}\n")
        f.write(f"Ion Mode: {mode}\n")
        f.write(f"Number of Parallel Jobs: {n_jobs if n_jobs else 'Auto'}\n")
        f.write(f"\nRun Date: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")


def load_ms2db_metadata(path=None):
    if path is None:
        path = 'all_ms2db_metadata.pkl'

    with open(path, 'rb') as f:
        metadata = pickle.load(f)

    db_id_to_mass = metadata.set_index('db_id')['monoisotopic_mass'].to_dict()
    db_id_to_inchi14 = metadata.set_index('db_id')['inchikey_14'].to_dict()

    return db_id_to_mass, db_id_to_inchi14


def process_file_wrapper(args):
    """Wrapper function to unpack arguments for pool.imap"""
    return process_file(*args)


def process_file(pkl_path, out_folder, mode, db_id_to_mass, db_id_to_inchi14):
    # basic search parameters
    ms1_tol_ppm = 20  # as ref mzs are recalculated

    basename = os.path.basename(pkl_path).split('_results')[0]

    out_path_best = os.path.join(os.path.join(out_folder, 'best'), f"{basename}_stitched_best.tsv")
    # out_path_all = os.path.join(os.path.join(out_folder, 'all'), f"{basename}_stitched_all.pkl")

    # load search results
    with open(pkl_path, 'rb') as f:
        search_results = pickle.load(f)

    results = []
    for spec in search_results:
        q_mass = spec['qry_mz']
        if mode == 'pos':
            q_mass -= 1.007276
        else:
            q_mass += 1.007276

        if spec['delta_result'] is None:
            this_spec_results = gen_from_empty_delta_results(spec)
            results.extend(this_spec_results)
            continue

        ##################
        # Filter ref_2 results before processing
        if 'score_arr' in spec['delta_result']:
            # Create mask for filtering based on score, peak, and usage thresholds
            mask = (spec['delta_result']['score_arr'] >= 0.6) & \
                   (spec['delta_result']['peak_arr'] >= 3) & \
                   (spec['delta_result']['usage_arr'] >= 0.10)

            filtered_indices = np.where(mask)[0]

            # If no results pass the filter, treat as empty delta result
            if len(filtered_indices) == 0:
                this_spec_results = gen_from_empty_delta_results(spec)
                results.extend(this_spec_results)
                continue

            # Create a filtered delta_result
            filtered_delta_result = {
                'ref_id_arr': [spec['delta_result']['ref_id_arr'][i] for i in filtered_indices],
                'score_arr': [spec['delta_result']['score_arr'][i] for i in filtered_indices],
                'peak_arr': [spec['delta_result']['peak_arr'][i] for i in filtered_indices],
                'usage_arr': [spec['delta_result']['usage_arr'][i] for i in filtered_indices]
            }

            # Create a copy of spec with filtered delta_result
            filtered_spec = spec.copy()
            filtered_spec['delta_result'] = filtered_delta_result

            this_spec_results = gen_from_delta_results(filtered_spec, q_mass, db_id_to_mass, ms1_tol_ppm)
        else:
            # Handle case where delta_result doesn't have expected arrays
            this_spec_results = gen_from_empty_delta_results(spec)

        results.extend(this_spec_results)

    if len(results) == 0:
        print(f"No results for {basename}")
        return

    # # save all results
    # with open(out_path_all, 'wb') as f:
    #     pickle.dump(results, f)

    #####
    # filter best results
    # for each unique ref_1 inchi, keep the best ref_2
    results_df = pd.DataFrame(results)

    results_df['ref_1_id'] = results_df['ref_1_id'].astype('str')
    results_df['ref_1_inchi'] = results_df['ref_1_id'].map(db_id_to_inchi14)

    float_columns = ['ref_1_score', 'ref_2_score', 'ref_1_usage', 'ref_2_usage', 'total_usage', 'ref_1_prec_frag_int', 
                     'ref_1_prec_frag_water_loss_int', 'ref_1_prec_nl_int', 'ref_1_prec_nl_water_loss_int']
    results_df[float_columns] = results_df[float_columns].astype('float32')

    # Sort by qry_id, ref_1_inchikey, total_usage, ref_1_score, ref_2_score
    results_df = results_df.sort_values(['qry_id', 'ref_1_inchi', 'total_usage', 'ref_1_score', 'ref_2_score'],
                                        ascending=[True, True, False, False, False])

    # Keep the highest scoring ref_2 for each qry_id + ref_1_inchikey
    results_df = results_df.drop_duplicates(subset=['qry_id', 'ref_1_inchi'], keep='first')

    # save best results
    results_df.to_csv(out_path_best, sep='\t', index=False)


def gen_from_empty_delta_results(spec):
    """
    Generate delta results from empty results
    """
    this_results = []
    for idx, ref_id in enumerate(spec['ref_id_arr']):
        this_results.append({
            'qry_id': spec['qry_id'],
            'qry_mz': spec['qry_mz'],
            'ref_1_id': ref_id,
            'ref_1_score': spec['score_arr'][idx],
            'ref_1_peak': spec['peak_arr'][idx],
            'ref_1_usage': spec['usage_arr'][idx],
            'ref_1_prec_frag_int': spec['ref_prec_int_frag_arr'][idx],
            'ref_1_prec_frag_water_loss_int': spec['ref_prec_int_frag_water_loss_arr'][idx],
            'ref_1_prec_nl_int': spec['ref_prec_int_nl_arr'][idx],
            'ref_1_prec_nl_water_loss_int': spec['ref_prec_int_nl_water_loss_arr'][idx],
            'ref_2_id': None,
            'ref_2_score': 0,
            'ref_2_peak': 0,
            'ref_2_usage': 0,
            'total_usage': spec['usage_arr'][idx]
        })
    return this_results


def gen_from_delta_results(spec, q_mass, db_id_to_mass, ms1_tol_ppm):

    # ref_2_mass_arr (delta mass)
    ref_2_mass_arr = np.empty(len(spec['delta_result']['ref_id_arr']), dtype=np.float64)
    for idx, ref_id in enumerate(spec['delta_result']['ref_id_arr']):
        ref_2_mass_arr[idx] = db_id_to_mass.get(ref_id, 0.0)

    ref_2_usage_arr = np.array(spec['delta_result']['usage_arr'], dtype=np.float32)

    this_results = []
    for idx, ref_id in enumerate(spec['ref_id_arr']):
        ref_1_mass = db_id_to_mass.get(ref_id, 0.0)
        ref_1_usage = spec['usage_arr'][idx]

        # find which ref_2_mass is within mz_tol and total usage > min_total_usage
        mask1 = np.abs(ref_2_mass_arr + ref_1_mass - 18.010565 - q_mass) < q_mass * ms1_tol_ppm / 1e6
        mask2 = ref_1_usage + ref_2_usage_arr >= min_total_usage

        mask = np.logical_and(mask1, mask2)
        mask_idx = np.where(mask)[0]

        if len(mask_idx) == 0:
            this_results.append({
                'qry_id': spec['qry_id'],
                'qry_mz': spec['qry_mz'],
                'ref_1_id': ref_id,
                'ref_1_score': spec['score_arr'][idx],
                'ref_1_peak': spec['peak_arr'][idx],
                'ref_1_usage': spec['usage_arr'][idx],
                'ref_1_prec_frag_int': spec['ref_prec_int_frag_arr'][idx],
                'ref_1_prec_frag_water_loss_int': spec['ref_prec_int_frag_water_loss_arr'][idx],
                'ref_1_prec_nl_int': spec['ref_prec_int_nl_arr'][idx],
                'ref_1_prec_nl_water_loss_int': spec['ref_prec_int_nl_water_loss_arr'][idx],
                'ref_2_id': None,
                'ref_2_score': 0,
                'ref_2_peak': 0,
                'ref_2_usage': 0,
                'total_usage': spec['usage_arr'][idx]
            })
        else:
            for idx2 in mask_idx:
                this_results.append({
                    'qry_id': spec['qry_id'],
                    'qry_mz': spec['qry_mz'],
                    'ref_1_id': ref_id,
                    'ref_1_score': spec['score_arr'][idx],
                    'ref_1_peak': spec['peak_arr'][idx],
                    'ref_1_usage': spec['usage_arr'][idx],
                    'ref_1_prec_frag_int': spec['ref_prec_int_frag_arr'][idx],
                    'ref_1_prec_frag_water_loss_int': spec['ref_prec_int_frag_water_loss_arr'][idx],
                    'ref_1_prec_nl_int': spec['ref_prec_int_nl_arr'][idx],
                    'ref_1_prec_nl_water_loss_int': spec['ref_prec_int_nl_water_loss_arr'][idx],
                    'ref_2_id': spec['delta_result']['ref_id_arr'][idx2],
                    'ref_2_score': spec['delta_result']['score_arr'][idx2],
                    'ref_2_peak': spec['delta_result']['peak_arr'][idx2],
                    'ref_2_usage': spec['delta_result']['usage_arr'][idx2],
                    'total_usage': spec['usage_arr'][idx] + spec['delta_result']['usage_arr'][idx2]
                })
    return this_results


if __name__ == "__main__":

    min_total_usage = 0.80

    start_time = time.time()
    main('/home/shipei/projects/revcos/search/results/pos_recheck',
         '/home/shipei/projects/revcos/search/results/pos_stitched', 'pos', None)
    print(f"Pos total time: {(time.time() - start_time) / 60} minutes")

    start_time = time.time()
    main('/home/shipei/projects/revcos/search/results/neg_recheck',
         '/home/shipei/projects/revcos/search/results/neg_stitched', 'neg', None)
    print(f"Neg total time: {(time.time() - start_time) / 60} minutes")


