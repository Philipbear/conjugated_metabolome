import pickle
import pandas as pd
import numpy as np
from numba import njit
import os
from _centroid_data import centroid_spectrum_for_search
from tqdm import tqdm
import time
from ms_entropy import read_one_spectrum
from multiprocessing import Pool, cpu_count

ms1_tol_ppm = 20
mz_tol = 0.025


def process_chunk(args):
    chunk_spectra, search_eng, mode, score_cutoff, min_matched_peak, min_spec_usage, min_prec_int, chunk_id = args

    chunk_results = []

    for spec in chunk_spectra:
        try:
            centroided_peaks = centroid_spectrum_for_search(spec['peaks'], width_da=mz_tol * 2.015)

            qry_mz = float(spec.get('pepmass') or spec.get('precursor_mz') or spec['precursormz'])

            matching_result = search_eng.search(
                precursor_mz=qry_mz,
                peaks=centroided_peaks,
                ms1_tolerance_in_da=ms1_tol_ppm * qry_mz / 1e6,
                ms2_tolerance_in_da=mz_tol,
                method='hybrid',
                precursor_ions_removal_da=-0.5,
                noise_threshold=0.0,
                min_ms2_difference_in_da=mz_tol * 2.015,
                reverse=True
            )

            _score_arr, _matched_peak_arr, _spec_usage_arr = matching_result['hybrid_search']
            v = np.where((_score_arr >= score_cutoff) &
                         (_matched_peak_arr >= min_matched_peak) &
                         (_spec_usage_arr >= min_spec_usage))[0]

            if len(v) == 0:
                continue

            for idx in v:
                ref_mz = search_eng[idx]['precursor_mz']
                ref_prec_int_frag1 = _get_fragment_intensity(centroided_peaks, ref_mz, mz_tol / 2)
                ref_prec_int_frag2 = _get_fragment_intensity(centroided_peaks, ref_mz - 18.010565, mz_tol / 2)

                new_target_mz = qry_mz - ref_mz + 18.010565
                if mode == 'pos':
                    new_target_mz += 1.007276
                else:
                    new_target_mz -= 1.007276
                ref_prec_int_nl1 = _get_fragment_intensity(centroided_peaks, new_target_mz, mz_tol / 2)
                ref_prec_int_nl2 = _get_fragment_intensity(centroided_peaks, new_target_mz - 18.010565, mz_tol / 2)

                ref_prec_int = max(ref_prec_int_frag1, ref_prec_int_frag2, ref_prec_int_nl1, ref_prec_int_nl2)

                if ref_prec_int < min_prec_int - 1e-4:
                    continue

                # taurine peak
                ref_prec_int_taurine = _get_fragment_intensity(centroided_peaks, 126.02194, mz_tol / 2)
                # glycine peak
                ref_prec_int_glycine = _get_fragment_intensity(centroided_peaks, 76.03930, mz_tol / 2)

                chunk_results.append({
                    'qry_name': spec.get('name'),
                    'qry_id': spec.get('spectrumid'),
                    'qry_mz': round(qry_mz, 4),
                    'ref_id': search_eng[idx]['id'],
                    'score': _score_arr[idx],
                    'peak': _matched_peak_arr[idx],
                    'usage': _spec_usage_arr[idx],
                    'ref_peak': ref_prec_int_frag1,
                    'ref_water_loss_peak': ref_prec_int_frag2,
                    'ref_nl_peak': ref_prec_int_nl1,
                    'ref_water_loss_nl_peak': ref_prec_int_nl2,
                    'tau_peak': ref_prec_int_taurine,
                    'gly_peak': ref_prec_int_glycine,
                })

        except Exception as e:
            print(f"Error processing spectrum: {str(e)}")

    # Save chunk results to a temporary file
    df = pd.DataFrame(chunk_results)
    out_folder = os.path.dirname(os.environ.get('MGF_PATH', 'input.mgf'))
    out_basename = os.path.basename(os.environ.get('MGF_PATH', 'input.mgf')).replace('.mgf', '')
    chunk_file = os.path.join(out_folder, f"{out_basename}_chunk_{chunk_id}_results.tsv")
    df.to_csv(chunk_file, sep='\t', index=False)

    return chunk_file


def add_metadata_to_chunk(args):
    chunk_file, ms2db_metadata_path = args
    try:
        # Load the chunk file
        df = pd.read_csv(chunk_file, sep='\t')

        if len(df) == 0:
            print(f"Chunk file {chunk_file} is empty. Skipping metadata addition.")
            return chunk_file

        # Load metadata
        with open(ms2db_metadata_path, 'rb') as f:
            metadata = pickle.load(f)

        metadata = metadata.rename(columns={'db_id': 'ref_id', 'recalc_prec_mz': 'ref_mz'})
        metadata = metadata[['ref_id', 'ref_mz', 'name', 'prec_type', 'inchikey', 'ion_mode',
                             'instrument_type', 'db', 'smiles', 'inchi', 'formula', 'monoisotopic_mass',
                             'npclassifier_superclass_results', 'npclassifier_class_results',
                             'npclassifier_pathway_results',
                             'npclassifier_isglycoside', 'classyfire_superclass', 'classyfire_class',
                             'classyfire_subclass']]

        # Merge with metadata
        df = df.merge(metadata, on='ref_id', how='left')

        # sort by score
        df = df.sort_values('score', ascending=False).reset_index(drop=True)

        # dereplicate by qry_id, inchikey
        df['inchikey_14'] = df['inchikey'].apply(lambda x: x[:14] if pd.notnull(x) else x)
        df = df.drop_duplicates(subset=['qry_id', 'inchikey_14', 'db'], keep='first').reset_index(drop=True)

        # Add OH count
        df['oh_count'] = df['qry_name'].apply(lambda x: _map_OH_count(x))

        df['ba_core_mass'] = 360.302825 + df['oh_count'] * 15.994915

        # possible adducts: M+H, M+Na, M+K, M+NH4, M+H-H2O, M+H-2H2O
        df['stitch_type'], df['adduct'] = zip(*df.apply(lambda x: _stitch(x['qry_mz'], x['ba_core_mass'],
                                                                        x['monoisotopic_mass'], x['tau_peak'],
                                                                        x['gly_peak']), axis=1))

        # filter out rows with None
        df = df[~df['stitch_type'].isnull()].reset_index(drop=True)

        # Save the updated chunk file
        df.to_csv(chunk_file, sep='\t', index=False)

        return chunk_file
    except Exception as e:
        print(f"Error adding metadata to chunk {chunk_file}: {str(e)}")
        return None


adduct_ls = ['M+H', 'M+Na', 'M+K', 'M+NH4', 'M+H-H2O', 'M+H-2H2O']
adduct_mass_ls = [1.007276, 22.989218, 38.963158, 18.0338254, 1.007276 - 18.010565, 1.007276 - 2*18.010565]

def _stitch(qry_mz, ba_core_mass, ref_exact_mass, tau_peak, gly_peak):
    # if observe tau or gly, try tau/gly + adducts, try ref + tau/gly + adducts
    # if not observe tau or gly, try ref + adducts

    # tau mass: 125.014664; gly mass: 75.032028

    if tau_peak > 0:
        # try tau + adducts
        for i, adduct_mass in enumerate(adduct_mass_ls):
            if abs(ba_core_mass + 125.014664 + adduct_mass - 18.010565 - qry_mz) < ms1_tol_ppm * qry_mz / 1e6:
                return 'Tau', adduct_ls[i]
        # try ref + tau + adducts
        for i, adduct_mass in enumerate(adduct_mass_ls):
            if abs(ba_core_mass + ref_exact_mass + 125.014664 + adduct_mass - 2*18.010565 - qry_mz) < ms1_tol_ppm * qry_mz / 1e6:
                return 'Tau_ref', adduct_ls[i]

    if gly_peak > 0:
        # try gly + adducts
        for i, adduct_mass in enumerate(adduct_mass_ls):
            if abs(ba_core_mass + 75.032028 + adduct_mass - 18.010565 - qry_mz) < ms1_tol_ppm * qry_mz / 1e6:
                return 'Gly', adduct_ls[i]
        # try ref + gly + adducts
        for i, adduct_mass in enumerate(adduct_mass_ls):
            if abs(ba_core_mass + ref_exact_mass + 75.032028 + adduct_mass - 2*18.010565 - qry_mz) < ms1_tol_ppm * qry_mz / 1e6:
                return 'Gly_ref', adduct_ls[i]

    # try ref + adducts
    for i, adduct_mass in enumerate(adduct_mass_ls):
        if abs(ba_core_mass + ref_exact_mass + adduct_mass - 18.010565 - qry_mz) < ms1_tol_ppm * qry_mz / 1e6:
            return 'Ref', adduct_ls[i]

    return None, None


def _map_OH_count(name):
    if 'Nonhydroxylated' in name:
        return 0
    elif 'Monohydroxylated' in name:
        return 1
    elif 'Dihydroxylated' in name:
        return 2
    elif 'Trihydroxylated' in name:
        return 3
    elif 'Tetrahydroxylated' in name:
        return 4
    elif 'Pentahydroxylated' in name:
        return 5
    else:
        return None


def chunk_list(lst, chunk_size):
    """Split list into chunks of specified size."""
    return [lst[i:i + chunk_size] for i in range(0, len(lst), chunk_size)]


def main(mgf_path='input.mgf',
         ms2db_metadata_path=None,
         mode='pos', score_cutoff=0.6, min_matched_peak=3,
         min_spec_usage=0.10, min_prec_int=0.0, n_jobs=None):
    # Store mgf_path in environment for child processes
    os.environ['MGF_PATH'] = mgf_path

    if n_jobs is None:
        n_jobs = max(1, cpu_count() - 1)

    write_params_to_file(mgf_path, mode, score_cutoff, min_matched_peak,
                         min_spec_usage, min_prec_int, n_jobs)

    with open(f'indexed_lib_{mode}.pkl', 'rb') as f:
        search_engine = pickle.load(f)


    print("Reading all spectra...")
    spectra = list(read_one_spectrum(mgf_path))

    # chunk size
    chunk_size = 1000
    chunks = chunk_list(spectra, chunk_size)

    # Prepare arguments for parallel processing
    process_args = [(chunk, search_engine, mode, score_cutoff, min_matched_peak,
                     min_spec_usage, min_prec_int, i) for i, chunk in enumerate(chunks)]

    print(f"Processing {len(chunks)} chunks in parallel...")
    with Pool(n_jobs) as pool:
        chunk_files = list(tqdm(pool.imap(process_chunk, process_args),
                                total=len(chunks), desc="Processing chunks"))

    # Filter out any None values
    chunk_files = [f for f in chunk_files if f is not None]

    if ms2db_metadata_path:
        print("Adding metadata to chunk files in parallel...")
        metadata_args = [(chunk_file, ms2db_metadata_path) for chunk_file in chunk_files]

        with Pool(n_jobs) as pool:
            updated_files = list(tqdm(pool.imap(add_metadata_to_chunk, metadata_args),
                                      total=len(chunk_files), desc="Adding metadata"))

        # Create a manifest file of all chunk files
        out_folder = os.path.dirname(mgf_path)
        out_basename = os.path.basename(mgf_path).replace('.mgf', '')
        manifest_file = os.path.join(out_folder, f"{out_basename}_chunks_manifest.txt")

        with open(manifest_file, 'w') as f:
            f.write("# Chunk files manifest\n")
            f.write(f"# Generated on: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"# Total chunks: {len(updated_files)}\n\n")
            for chunk_file in updated_files:
                if chunk_file:
                    f.write(f"{chunk_file}\n")

        print(f"Manifest file created: {manifest_file}")
        print(f"Total processed chunks: {len(updated_files)}")


@njit
def _get_fragment_intensity(peaks, mz, tol):
    max_intensity = np.max(peaks[:, 1])
    idx = np.abs(peaks[:, 0] - mz) < tol
    if not np.any(idx):
        return 0.0
    return np.max(peaks[idx, 1]) / max_intensity


def write_params_to_file(mgf, mode, score_cutoff, min_matched_peak,
                         min_spec_usage, min_prec_int, n_jobs):
    folder = os.path.dirname(mgf)
    params_file = os.path.join(folder, 'search_parameters.txt')
    with open(params_file, 'w') as f:
        f.write("Search Parameters\n")
        f.write("================\n\n")
        f.write(f"Input MGF: {mgf}\n")
        f.write(f"Ion Mode: {mode}\n")
        f.write(f"Score Cutoff: {score_cutoff}\n")
        f.write(f"Minimum Matched Peaks: {min_matched_peak}\n")
        f.write(f"Minimum Spectral Usage: {min_spec_usage}\n")
        f.write(f"Minimum Precursor Intensity: {min_prec_int}\n")
        f.write(f"Number of Parallel Jobs: {n_jobs if n_jobs else 'Auto'}\n")
        f.write(f"\nRun Date: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")


def combine_results():
    tsv_files = os.listdir('ba')
    tsv_files = [x for x in tsv_files if x.endswith('.tsv') and not x.startswith('.')]
    tsv_files = [os.path.join('ba', x) for x in tsv_files]

    all_results = []
    for file in tqdm(tsv_files, desc='Combining results'):
        df = pd.read_csv(file, sep='\t', low_memory=False)
        if len(df) == 0:
            continue
        all_results.append(df)

    all_results = pd.concat(all_results)
    all_results.to_csv('ba/all_ba_results.tsv', sep='\t', index=False)


if __name__ == '__main__':

    start_time = time.time()
    main('ba/GNPS-BILE-ACID-MODIFICATIONS.mgf',
         'all_ms2db_metadata.pkl',
         'pos', 0.6, 2, 0.05, 0.0, 18)

    combine_results()

    print(f"Total time: {(time.time() - start_time) / 60} minutes")
