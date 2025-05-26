import os
import pickle
import time
from multiprocessing import Pool, cpu_count

import numpy as np
from numba import njit
from tqdm import tqdm

from _centroid_data import centroid_spectrum_for_search


def main(folder='/home/shipei/projects/cluster/unclustered_all',
         results_folder='/home/shipei/projects/revcos/search/search_results',
         out_folder=None,
         mode='pos',
         n_jobs=None):
    """
    Recalculate precursor intensities in existing search results
    """
    # If no output folder is specified, use the results folder
    if out_folder is None:
        out_folder = results_folder
    
    os.makedirs(out_folder, exist_ok=True)

    # Map MGF files to result files
    all_mgfs = os.listdir(folder)
    all_mgfs = [x for x in all_mgfs if x.endswith('.mgf') and not x.startswith('.')]
    
    all_results = os.listdir(results_folder)
    all_results = [x for x in all_results if x.endswith('.pkl') and not x.startswith('search_parameters')]
    
    # Create mapping from MGF filenames to results filenames
    mgf_to_result = {}
    for mgf in all_mgfs:
        result_basename = mgf.replace('.mgf', '') + '_results.pkl'
        if result_basename in all_results:
            mgf_to_result[os.path.join(folder, mgf)] = os.path.join(results_folder, result_basename)
    
    print(f"Found {len(mgf_to_result)} matching MGF and result files")
    
    # load search engines and metadata
    with open(f'indexed_lib_{mode}.pkl', 'rb') as f:
        search_engine = pickle.load(f)
        
    # Create mapping from ref_id to index in search_eng
    ref_id_to_idx = {}
    for idx in range(len(search_engine.precursor_mz_array)):
        ref_id_to_idx[search_engine[idx]['id']] = idx

    # Determine number of processes
    if n_jobs is None:
        n_jobs = max(1, cpu_count() - 1)

    # Create process pool and process files in parallel
    print(f"Processing {len(mgf_to_result)} files using {n_jobs} processes...")
    with Pool(n_jobs) as pool:
        args = [(mgf_path, result_path, search_engine, ref_id_to_idx, out_folder, mode)
                for mgf_path, result_path in mgf_to_result.items()]

        for _ in tqdm(pool.imap_unordered(process_file_wrapper, args),
                      total=len(mgf_to_result),
                      desc="Recalculating precursor intensities"):
            pass


def process_file_wrapper(args):
    """Wrapper function to unpack arguments for pool.imap"""
    return process_file(*args)


def process_file(mgf_path, result_path, search_eng, ref_id_to_idx, out_folder, mode='pos'):
    """Process a single MGF-result file pair to recalculate precursor intensities"""
    
    out_path = os.path.join(out_folder, os.path.basename(result_path))
    
    try:
        # Load the results
        with open(result_path, 'rb') as f:
            results = pickle.load(f)
        
        # Load all spectra from MGF into a dictionary for fast access
        spectra = {}
        for spec in iterate_mgf(mgf_path):
            spectra[spec['title']] = {
                'mz': float(spec['pepmass']),
                'peaks': centroid_spectrum_for_search(spec['peaks'], width_da=0.025 * 2.015)
            }
        
        # Process each result
        for idx, result in enumerate(results):
            qry_id = result['qry_id']
            if qry_id not in spectra:
                print(f"Warning: Could not find spectrum for {qry_id}")
                continue
            
            spec_data = spectra[qry_id]
            qry_mz = spec_data['mz']
            centroided_peaks = spec_data['peaks']
            
            # Get the number of ref IDs
            n_refs = len(result['ref_id_arr'])
            
            # Initialize the new arrays
            ref_prec_int_frag1_arr = np.zeros(n_refs, dtype=np.float16)
            ref_prec_int_frag2_arr = np.zeros(n_refs, dtype=np.float16)
            ref_prec_int_nl1_arr = np.zeros(n_refs, dtype=np.float16)
            ref_prec_int_nl2_arr = np.zeros(n_refs, dtype=np.float16)
            
            for i in range(n_refs):
                ref_id = result['ref_id_arr'][i]
                
                # Get the reference precursor m/z
                ref_spec_idx = ref_id_to_idx.get(ref_id)
                if ref_spec_idx is None:
                    continue
                
                ref_mz = search_eng[ref_spec_idx]['precursor_mz']
                
                # Calculate the four precursor intensities
                ref_prec_int_frag1 = _get_fragment_intensity(centroided_peaks, ref_mz, 0.01)
                ref_prec_int_frag2 = _get_fragment_intensity(centroided_peaks, ref_mz - 18.010565, 0.01)
                
                new_target_mz = qry_mz - ref_mz + 18.010565
                if mode == 'pos':
                    new_target_mz += 1.007276
                else:
                    new_target_mz -= 1.007276
                    
                ref_prec_int_nl1 = _get_fragment_intensity(centroided_peaks, new_target_mz, 0.01)
                ref_prec_int_nl2 = _get_fragment_intensity(centroided_peaks, new_target_mz - 18.010565, 0.01)
                
                # Store the values
                ref_prec_int_frag1_arr[i] = ref_prec_int_frag1
                ref_prec_int_frag2_arr[i] = ref_prec_int_frag2
                ref_prec_int_nl1_arr[i] = ref_prec_int_nl1
                ref_prec_int_nl2_arr[i] = ref_prec_int_nl2
            
            # Update the result with new arrays and remove the old one
            result['ref_prec_int_frag_arr'] = ref_prec_int_frag1_arr
            result['ref_prec_int_frag_water_loss_arr'] = ref_prec_int_frag2_arr
            result['ref_prec_int_nl_arr'] = ref_prec_int_nl1_arr
            result['ref_prec_int_nl_water_loss_arr'] = ref_prec_int_nl2_arr
            
            if 'ref_prec_int_arr' in result:
                del result['ref_prec_int_arr']
        
        # Save the updated results
        with open(out_path, 'wb') as f:
            pickle.dump(results, f)
            
        return True
        
    except Exception as e:
        print(f"Error processing {mgf_path} / {result_path}: {str(e)}")
        return False


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
    
    # main('/home/shipei/projects/revcos/cluster/dataset_clustered/pos',
    #      '/home/shipei/projects/revcos/search/results/pos',
    #      '/home/shipei/projects/revcos/search/results/pos_recheck',
    #      'pos')
    
    main('/home/shipei/projects/revcos/cluster/dataset_clustered/neg',
         '/home/shipei/projects/revcos/search/results/neg',
         '/home/shipei/projects/revcos/search/results/neg_recheck',
         'neg')

    print(f"Total time: {(time.time() - start_time) / 60:.2f} minutes")