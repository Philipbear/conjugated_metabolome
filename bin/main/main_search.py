"""
Unified reverse cosine search script combining both rounds into a single step.
Some improvement: water loss precursor (eg, cholic acid water loss precursor instead of cholic acid, as -OH is already removed after conjugation)
"""

import argparse
import os
import pickle
import time
from multiprocessing import Pool, cpu_count

import numpy as np
from numba import njit
from tqdm import tqdm

from _centroid_data import centroid_spectrum_for_search


def main(folder='/home/shipei/projects/cluster/unclustered_all',
         out_folder='/home/shipei/projects/revcos/search/search_results',
         mode='pos',
         top_score_cutoff=0.8, top_min_matched_peak=3, top_min_spec_usage=0.30, top_min_prec_int=0.05,
         delta_score_cutoff=0.5, delta_min_matched_peak=2, delta_min_spec_usage=0.10,
         n_jobs=None):
    """
    Main function to run the unified search process
    """
    os.makedirs(out_folder, exist_ok=True)

    # Write parameters to file
    write_params_to_file(out_folder, folder, mode,
                         top_score_cutoff, top_min_matched_peak, top_min_spec_usage, top_min_prec_int,
                         delta_score_cutoff, delta_min_matched_peak, delta_min_spec_usage,
                         n_jobs)

    all_mgfs = os.listdir(folder)
    all_mgfs = [x for x in all_mgfs if x.endswith('.mgf') and not x.startswith('.')]
    all_mgfs = [os.path.join(folder, x) for x in all_mgfs]

    # load search engines and metadata
    with open(f'indexed_lib_{mode}.pkl', 'rb') as f:
        search_engine = pickle.load(f)

    with open('db_id_to_mass.pkl', 'rb') as f:
        db_id_to_mass = pickle.load(f)

    # Determine number of processes
    if n_jobs is None:
        n_jobs = max(1, cpu_count() - 1)

    # skip files that have been processed
    all_results = os.listdir(out_folder)
    all_results = [x for x in all_results if x.endswith('.pkl')]

    mgfs_to_process = []
    for _mgf in all_mgfs:
        out_basename = os.path.basename(_mgf).replace('.mgf', '') + '_results.pkl'
        if out_basename not in all_results:
            mgfs_to_process.append(_mgf)
        else:
            print(f"File already processed: {_mgf}")

    # Create process pool and process files in parallel
    print(f"Processing {len(mgfs_to_process)} files using {n_jobs} processes...")
    with Pool(n_jobs) as pool:
        args = [(mgf, search_engine, db_id_to_mass, out_folder, mode,
                 top_score_cutoff, top_min_matched_peak, top_min_spec_usage, top_min_prec_int,
                 delta_score_cutoff, delta_min_matched_peak, delta_min_spec_usage)
                for mgf in mgfs_to_process]

        for _ in tqdm(pool.imap_unordered(process_mgf_wrapper, args),
                      total=len(mgfs_to_process),
                      desc="Processing MGF files"):
            pass


def write_params_to_file(out_folder, folder, mode,
                         top_score_cutoff, top_min_matched_peak, top_min_spec_usage, top_min_prec_int,
                         delta_score_cutoff, delta_min_matched_peak, delta_min_spec_usage,
                         n_jobs):
    """Write search parameters to a file"""
    params_file = os.path.join(out_folder, 'search_parameters.txt')
    with open(params_file, 'w') as f:
        f.write("Search Parameters\n")
        f.write("================\n\n")
        f.write(f"Input Folder: {folder}\n")
        f.write(f"Output Folder: {out_folder}\n")
        f.write(f"Ion Mode: {mode}\n")
        f.write("\nFirst Pass Parameters:\n")
        f.write("\nTop Hit Parameters:\n")
        f.write(f"Top Score Cutoff: {top_score_cutoff}\n")
        f.write(f"Top Minimum Matched Peaks: {top_min_matched_peak}\n")
        f.write(f"Top Minimum Spectral Usage: {top_min_spec_usage}\n")
        f.write(f"Top Minimum Precursor Intensity: {top_min_prec_int}\n")
        f.write("\nDelta Search Parameters:\n")
        f.write(f"Delta Score Cutoff: {delta_score_cutoff}\n")
        f.write(f"Delta Minimum Matched Peaks: {delta_min_matched_peak}\n")
        f.write(f"Delta Minimum Spectral Usage: {delta_min_spec_usage}\n")
        f.write(f"Number of Parallel Jobs: {n_jobs if n_jobs else 'Auto'}\n")
        f.write(f"\nRun Date: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")


def process_mgf_wrapper(args):
    """Wrapper function to unpack arguments for pool.imap"""
    return process_mgf(*args)


def process_mgf(mgf, search_eng, db_id_to_mass, out_folder, mode='pos',
                top_score_cutoff=0.8, top_min_matched_peak=4, top_min_spec_usage=0.40, top_min_prec_int=0.05,
                delta_score_cutoff=0.5, delta_min_matched_peak=2, delta_min_spec_usage=0.10):
    """Process a single MGF file with unified search"""

    out_basename = os.path.basename(mgf).replace('.mgf', '') + '_results.pkl'
    out_path = os.path.join(out_folder, out_basename)
    
    ms1_tol_ppm = 25
    mz_tol = 0.025
    all_results = []

    try:
        for spec in iterate_mgf(mgf):
            centroided_peaks = centroid_spectrum_for_search(spec['peaks'], width_da=mz_tol * 2.015)
            qry_mz = float(spec['pepmass'])

            # search
            matching_result = search_eng.search(
                precursor_mz=qry_mz,
                peaks=centroided_peaks,
                ms1_tolerance_in_da=ms1_tol_ppm * qry_mz / 1e6,
                ms2_tolerance_in_da=mz_tol,
                method='open',
                precursor_ions_removal_da=-0.5,
                noise_threshold=0.0,
                min_ms2_difference_in_da=mz_tol * 2.015,
                reverse=True
            )

            _score_arr, _matched_peak_arr, _spec_usage_arr = matching_result['open_search']

            # Filter by top matching cutoffs
            v = np.where((_score_arr >= top_score_cutoff) &
                         (_matched_peak_arr >= top_min_matched_peak) &
                         (_spec_usage_arr >= top_min_spec_usage))[0]

            if len(v) == 0:
                continue

            # Process top hits
            ref_id_arr = []
            ref_prec_int_frag_arr = []
            ref_prec_int_frag_water_loss_arr = []
            ref_prec_int_nl_arr = []
            ref_prec_int_nl_water_loss_arr = []
            score_arr = []
            matched_peak_arr = []
            spec_usage_arr = []
            delta_mass_arr = []

            q_mass = qry_mz
            if mode == 'pos':
                q_mass -= 1.007276
            else:
                q_mass += 1.007276

            for idx in v:
                ref_mz = search_eng[idx]['precursor_mz']
                if ref_mz > qry_mz - 10:
                    continue

                # Find precursor intensity
                ref_prec_int_frag1 = _get_fragment_intensity(centroided_peaks, ref_mz, 0.01)
                ref_prec_int_frag2 = _get_fragment_intensity(centroided_peaks, ref_mz - 18.010565, 0.01)

                new_target_mz = qry_mz - ref_mz + 18.010565
                if mode == 'pos':
                    new_target_mz += 1.007276
                else:
                    new_target_mz -= 1.007276
                ref_prec_int_nl1 = _get_fragment_intensity(centroided_peaks, new_target_mz, 0.01)
                ref_prec_int_nl2 = _get_fragment_intensity(centroided_peaks, new_target_mz - 18.010565, 0.01)

                ref_prec_int = max(ref_prec_int_frag1, ref_prec_int_frag2, ref_prec_int_nl1, ref_prec_int_nl2)

                if ref_prec_int < top_min_prec_int:
                    continue

                ref_id = search_eng[idx]['id']
                ref_id_arr.append(ref_id)
                ref_prec_int_frag_arr.append(ref_prec_int_frag1)
                ref_prec_int_frag_water_loss_arr.append(ref_prec_int_frag2)
                ref_prec_int_nl_arr.append(ref_prec_int_nl1)
                ref_prec_int_nl_water_loss_arr.append(ref_prec_int_nl2)
                score_arr.append(_score_arr[idx])
                matched_peak_arr.append(_matched_peak_arr[idx])
                spec_usage_arr.append(_spec_usage_arr[idx])
                delta_mass_arr.append(q_mass - db_id_to_mass.get(ref_id, 0.0))

            if not ref_id_arr:
                continue

            # Convert lists to numpy arrays
            ref_id_arr = np.array(ref_id_arr)
            ref_prec_int_frag_arr = np.array(ref_prec_int_frag_arr, dtype=np.float16)
            ref_prec_int_frag_water_loss_arr = np.array(ref_prec_int_frag_water_loss_arr, dtype=np.float16)
            ref_prec_int_nl_arr = np.array(ref_prec_int_nl_arr, dtype=np.float16)
            ref_prec_int_nl_water_loss_arr = np.array(ref_prec_int_nl_water_loss_arr, dtype=np.float16)
            score_arr = np.array(score_arr, dtype=np.float16)
            matched_peak_arr = np.array(matched_peak_arr, dtype=np.int8)
            spec_usage_arr = np.array(spec_usage_arr, dtype=np.float16)
            delta_mass_arr = np.array(delta_mass_arr, dtype=np.float32)

            # for delta masses

            # Filter delta search results
            v = np.where((_score_arr >= delta_score_cutoff) &
                         (_matched_peak_arr >= delta_min_matched_peak) &
                         (_spec_usage_arr >= delta_min_spec_usage))[0]

            if len(v) > 0:
                delta_ref_id_arr = []
                delta_score_arr = []
                delta_matched_peak_arr = []
                delta_spec_usage_arr = []

                for idx in v:
                    ref_id = search_eng[idx]['id']
                    ref_mass = db_id_to_mass.get(ref_id, 0.0)

                    # Check if ref mass matches any delta mass
                    if np.any(np.abs(delta_mass_arr + 18.010565 - ref_mass) < ms1_tol_ppm * delta_mass_arr / 1e6):
                        delta_ref_id_arr.append(ref_id)
                        delta_score_arr.append(_score_arr[idx])
                        delta_matched_peak_arr.append(_matched_peak_arr[idx])
                        delta_spec_usage_arr.append(_spec_usage_arr[idx])

                if delta_ref_id_arr:
                    delta_result = {
                        'ref_id_arr': np.array(delta_ref_id_arr),
                        'score_arr': np.array(delta_score_arr, dtype=np.float16),
                        'peak_arr': np.array(delta_matched_peak_arr, dtype=np.int8),
                        'usage_arr': np.array(delta_spec_usage_arr, dtype=np.float16)
                    }
                else:
                    delta_result = None
            else:
                delta_result = None

            # Store results
            result = {
                'qry_id': spec['title'],
                'qry_mz': round(float(spec['pepmass']), 4),
                'ref_id_arr': ref_id_arr,
                'score_arr': score_arr,
                'peak_arr': matched_peak_arr,
                'usage_arr': spec_usage_arr,
                'ref_prec_int_frag_arr': ref_prec_int_frag_arr,
                'ref_prec_int_frag_water_loss_arr': ref_prec_int_frag_water_loss_arr,
                'ref_prec_int_nl_arr': ref_prec_int_nl_arr,
                'ref_prec_int_nl_water_loss_arr': ref_prec_int_nl_water_loss_arr,
                'delta_mass_arr': delta_mass_arr,
                'delta_result': delta_result
            }

            all_results.append(result)

        # Save results
        with open(out_path, 'wb') as f:
            pickle.dump(all_results, f)

    except Exception as e:
        print(f"Error processing {mgf}: {str(e)}")
        return


def read_mgf_spectrum(file_obj):
    """Read a single spectrum block from an open MGF file"""
    spectrum = {
        'title': '',
        'pepmass': 0.0,
        'charge': '',
        'peaks': []
    }

    for line in file_obj:
        if line.strip() == 'BEGIN IONS':
            break
    else:
        return None

    for line in file_obj:
        line = line.strip()

        if not line:  # Skip empty lines
            continue

        if line == 'END IONS':
            if spectrum['peaks']:  # Only return if we have peaks
                return spectrum
            break

        if line.startswith('TITLE='):
            spectrum['title'] = line[6:]
        elif line.startswith('PEPMASS='):
            spectrum['pepmass'] = float(line[8:].split()[0])  # Handle additional intensity value
        elif line.startswith('CHARGE='):
            spectrum['charge'] = line[7:]
        elif line and not line.startswith(('BEGIN', 'END')):  # Should be a peak line
            mz, intensity = line.split()
            spectrum['peaks'].append((float(mz), float(intensity)))

    return None


def iterate_mgf(mgf_path, buffer_size=8192):
    """Iterate through spectra in an MGF file efficiently using buffered reading"""
    with open(mgf_path, 'r', buffering=buffer_size) as f:
        while True:
            spectrum = read_mgf_spectrum(f)
            if spectrum is None:
                break
            yield spectrum


@njit
def _get_fragment_intensity(peaks, mz, tol):
    """Get the maximum intensity of fragments within tolerance of target m/z"""
    max_intensity = np.max(peaks[:, 1])

    idx = np.abs(peaks[:, 0] - mz) < tol
    if not np.any(idx):
        return 0.0
    return np.max(peaks[idx, 1]) / max_intensity


if __name__ == '__main__':

    start_time = time.time()
    
    main(folder='/home/shipei/projects/revcos/cluster/dataset_clustered/pos',
         out_folder='/home/shipei/projects/revcos/search/results/pos',
         mode='pos',
         top_score_cutoff=0.7, top_min_matched_peak=3, top_min_spec_usage=0.20, top_min_prec_int=0.05,
         delta_score_cutoff=0.7, delta_min_matched_peak=2, delta_min_spec_usage=0.05,
         n_jobs=40)
    
    # main(folder='/home/shipei/projects/revcos/cluster/dataset_clustered/neg',
    #      out_folder='/home/shipei/projects/revcos/search/results/neg',
    #      mode='neg',
    #      top_score_cutoff=0.7, top_min_matched_peak=3, top_min_spec_usage=0.20, top_min_prec_int=0.05,
    #      delta_score_cutoff=0.7, delta_min_matched_peak=2, delta_min_spec_usage=0.05,
    #      n_jobs=40)

    print(f"Total time: {(time.time() - start_time) / 60:.2f} minutes")
    