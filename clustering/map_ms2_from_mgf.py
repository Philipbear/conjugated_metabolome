"""
this is to give MS2 given a USI, for analysis purpose (filter out annotations based on MassQL)
"""
from collections import defaultdict
import os
import pickle
from multiprocessing import Pool, cpu_count

import pandas as pd
from tqdm import tqdm


def read_mgf_spectrum(file_obj):
    """Read a single spectrum block from an open MGF file.

    Args:
        file_obj: An opened file object positioned at the start of a spectrum

    Returns:
        dict: Spectrum information containing title, pepmass, charge and peaks
        or None if end of file is reached
    """
    spectrum = {
        'title': '',
        'pepmass': 0.0,
        'charge': '',
        'peaks': []
    }

    # Skip any empty lines before BEGIN IONS
    for line in file_obj:
        if line.strip() == 'BEGIN IONS':
            break
    else:  # EOF reached
        return None

    # Read spectrum metadata and peaks
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
    """Iterate through spectra in an MGF file efficiently using buffered reading.

    Args:
        mgf_path: Path to the MGF file
        buffer_size: Read buffer size in bytes

    Yields:
        dict: Spectrum information containing title, pepmass, charge and peaks
    """
    with open(mgf_path, 'r', buffering=buffer_size) as f:
        while True:
            spectrum = read_mgf_spectrum(f)
            if spectrum is None:
                break
            yield spectrum


def process_mgf(mgf, usi_list):
    """Process a single MGF file and return a dictionary of spectra.

    Args:
        mgf: Path to the MGF file

    Returns:
        dict: Dictionary mapping spectrum titles to peak data
    """
    out = {}
    try:
        for spec in iterate_mgf(mgf):
            if spec['title'] in usi_list:
                new_spec = {
                    'precmz': round(spec['pepmass'], 5),
                    'peaks': spec['peaks']
                }
                out[spec['title']] = new_spec
    except Exception as e:
        print(f"Error processing {mgf}: {str(e)}")
        return {}
    return out


def main(folder, out_name):
    """Process all MGF files in a folder and save results to a pickle file.

    Args:
        folder: Path to folder containing MGF files
        out_name: Path where to save the output pickle file
    """
    # Get list of MGF files
    all_mgfs = [os.path.join(folder, x) for x in os.listdir(folder)
                if x.endswith('.mgf') and not x.startswith('.')]

    if not all_mgfs:
        print(f"No MGF files found in {folder}")
        return

    # load usi list with annotations
    if 'pos' in folder:
        with open('data/pos_usis.pkl', 'rb') as f:
            usi_list = pickle.load(f)
    else:
        with open('data/neg_usis.pkl', 'rb') as f:
            usi_list = pickle.load(f)

    # Process files in parallel and collect results
    combined_dict = {}
    for mgf_file in tqdm(all_mgfs, total=len(all_mgfs), desc='Processing MGFs'):
        result = process_mgf(mgf_file, usi_list)

        # Combine results
        combined_dict.update(result)

    # Save combined dictionary
    print(f"Saving {len(combined_dict)} spectra to {out_name}")
    with open(out_name, 'wb') as f:
        pickle.dump(combined_dict, f)


def gen_usi_list():

    pos_files = os.listdir('/cmpd_analysis/data/pos')
    pos_files = [os.path.join('/cmpd_analysis/data/pos', x) for x in pos_files if x.endswith('.tsv') and not x.startswith('.')]

    all_pos_usi_ls = []
    for file in tqdm(pos_files):
        df = pd.read_csv(file, sep='\t', low_memory=False)
        list_to_add = df['qry_id'].tolist()
        list_to_add = list(set(list_to_add))
        all_pos_usi_ls.extend(list_to_add)
    all_pos_usi_ls = set(all_pos_usi_ls)
    # save as pickle
    with open('data/pos_usis.pkl', 'wb') as f:
        pickle.dump(all_pos_usi_ls, f)

    df = pd.read_csv('/cmpd_analysis/data/neg/neg_best_annotated.tsv', sep='\t', low_memory=False)
    usi_list = df['qry_id'].tolist()
    usi_list = set(usi_list)
    # save as pickle
    with open('data/neg_usis.pkl', 'wb') as f:
        pickle.dump(usi_list, f)


if __name__ == '__main__':

    ########## 1. generate USIs with annotations ##########
    # gen_usi_list()

    ########## 2. process_mgf ##########

    main('/home/shipei/projects/revcos/cluster/dataset_clustered/pos',
         '/home/shipei/projects/revcos/cluster/dataset_clustered/pos_all.pkl')

    main('/home/shipei/projects/revcos/cluster/dataset_clustered/neg',
         '/home/shipei/projects/revcos/cluster/dataset_clustered/neg_all.pkl')
