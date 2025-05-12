import pickle
import pandas as pd
import numpy as np
from numba import njit
import os
from _centroid_data import centroid_spectrum_for_search
from tqdm import tqdm
import argparse
import time
from ms_entropy import read_one_spectrum
from multiprocessing import Pool, cpu_count


def process_chunk(args):
    chunk_spectra, search_eng, mode, score_cutoff, min_matched_peak, min_spec_usage, min_prec_int = args
    ms1_tol_ppm = 20
    mz_tol = 0.025
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
                ref_prec_int_frag = _get_fragment_intensity(centroided_peaks, ref_mz, mz_tol / 2)

                new_target_mz = qry_mz - ref_mz + 18.010565
                if mode == 'pos':
                    new_target_mz += 1.007276
                else:
                    new_target_mz -= 1.007276
                ref_prec_int_nl = _get_fragment_intensity(centroided_peaks, new_target_mz, mz_tol / 2)

                ref_prec_int = max(ref_prec_int_frag, ref_prec_int_nl)
                if ref_prec_int < min_prec_int - 1e-4:
                    continue

                qry_id = (spec.get('title') or spec.get('spectrum_id') or
                          spec.get('spectrumid') or spec.get('feature_id') or
                          spec.get('scans') or len(chunk_results))

                chunk_results.append({
                    'qry_id': qry_id,
                    'qry_mz': round(qry_mz, 4),
                    'ref_id': search_eng[idx]['id'],
                    'score': _score_arr[idx],
                    'peak': _matched_peak_arr[idx],
                    'usage': _spec_usage_arr[idx],
                    'ref_prec_int_in_qry_spec': ref_prec_int
                })

        except Exception as e:
            print(f"Error processing spectrum: {str(e)}")

    return chunk_results


def chunk_list(lst, chunk_size):
    """Split list into chunks of specified size."""
    return [lst[i:i + chunk_size] for i in range(0, len(lst), chunk_size)]


def main(mgf_path='input.mgf',
         ms2db_metadata_path=None,
         mode='pos', score_cutoff=0.6, min_matched_peak=3,
         min_spec_usage=0.10, min_prec_int=0.0, n_jobs=None):
    write_params_to_file(mgf_path, mode, score_cutoff, min_matched_peak,
                         min_spec_usage, min_prec_int, n_jobs)

    with open(f'data/indexed_lib_gnps_{mode}.pkl', 'rb') as f:
        search_engine = pickle.load(f)

    # with open(f'data/indexed_lib_{mode}.pkl', 'rb') as f:
    #     search_engine = pickle.load(f)

    if n_jobs is None:
        n_jobs = max(1, cpu_count() - 1)

    print("Reading all spectra...")
    spectra = list(read_one_spectrum(mgf_path))

    # chunk size
    chunk_size = 200
    chunks = chunk_list(spectra, chunk_size)

    # Prepare arguments for parallel processing
    process_args = [(chunk, search_engine, mode, score_cutoff, min_matched_peak,
                     min_spec_usage, min_prec_int) for chunk in chunks]

    print(f"Processing {len(chunks)} chunks in parallel...")
    with Pool(n_jobs) as pool:
        all_results = list(tqdm(pool.imap(process_chunk, process_args),
                                total=len(chunks), desc="Processing chunks"))

    # Flatten results
    flat_results = [item for chunk_result in all_results for item in chunk_result]

    # Save results
    out_folder = os.path.dirname(mgf_path)
    out_basename = os.path.basename(mgf_path).replace('.mgf', '')
    out_basename = f"{out_basename}_search_results.tsv"
    out_path = os.path.join(out_folder, out_basename)

    df = pd.DataFrame(flat_results)

    if ms2db_metadata_path:
        print("Adding metadata...")
        with open(ms2db_metadata_path, 'rb') as f:
            metadata = pickle.load(f)

        metadata = metadata.rename(columns={'db_id': 'ref_id', 'recalc_prec_mz': 'ref_mz'})
        metadata = metadata[['ref_id', 'ref_mz', 'name', 'prec_type', 'inchikey', 'ion_mode',
                             'instrument_type', 'db', 'smiles', 'inchi', 'formula', 'monoisotopic_mass',
                             'npclassifier_superclass_results', 'npclassifier_class_results',
                             'npclassifier_pathway_results',
                             'npclassifier_isglycoside', 'classyfire_superclass', 'classyfire_class',
                             'classyfire_subclass']]
        df = df.merge(metadata, on='ref_id', how='left')

    df.to_csv(out_path, sep='\t', index=False)


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


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Chunked parallel reverse search')
    parser.add_argument('--input_path', '-i', help='MGF file path', default=None)
    parser.add_argument('--mode', '-m', help='Ion mode: pos or neg', default='pos')
    parser.add_argument('--score', '-sc', type=float, help='Score cutoff', default=0.70)
    parser.add_argument('--peak', '-p', type=int, help='Minimum number of matched peaks', default=3)
    parser.add_argument('--usage', '-u', type=float, help='Minimum spectral usage', default=0.10)
    parser.add_argument('--prec_int', '-pi', type=float,
                        help='Minimum precursor intensity (fragment or neutral loss)', default=0.05)
    parser.add_argument('--n_jobs', '-j', type=int, help='Number of parallel processes', default=None)

    args = parser.parse_args()
    start_time = time.time()
    main(args.input_path,
         '//db/ms2db/all/all_ms2db_metadata.pkl',
         args.mode, args.score, args.peak, args.usage, args.prec_int, args.n_jobs)
    print(f"Total time: {(time.time() - start_time) / 60} minutes")