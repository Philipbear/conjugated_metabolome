"""
get all mz from clustered mgfs
"""
import argparse
import os
import pickle
import time
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


def process_mgf(mgf):
    """Process a single MGF file and return a dictionary of spectra.

    Args:
        mgf: Path to the MGF file

    Returns:
        dict: Dictionary mapping spectrum titles to peak data
    """
    mz_ls = []
    try:
        for spec in iterate_mgf(mgf):
            mz_ls.append(spec['pepmass'])
    except Exception as e:
        print(f"Error processing {mgf}: {str(e)}")
        return []
    return mz_ls


def main(folder, out_name):
    """Process all MGF files in a folder and save results to a pickle file.

    Args:
        folder: Path to folder containing MGF files
        out_name: Path where to save the output pickle file
    """
    # Get list of MGF files
    all_mgfs = [os.path.join(folder, x) for x in os.listdir(folder)
                if x.endswith('.mgf') and not x.startswith('.')]

    # Set up multiprocessing pool
    num_processes = min(cpu_count(), len(all_mgfs))
    print(f"Using {num_processes} processes to process {len(all_mgfs)} MGF files")

    # Process files in parallel
    all_mz_lists = []
    with Pool(processes=num_processes) as pool:
        results = list(tqdm(
            pool.imap(process_mgf, all_mgfs),
            total=len(all_mgfs),
            desc="Processing MGF files"
        ))

        # Combine all mz lists
        all_mz_lists = [mz for sublist in results for mz in sublist]

    # Save combined list
    print(f"Saving {len(all_mz_lists)} precursor m/z values to {out_name}")
    with open(out_name, 'wb') as f:
        pickle.dump(all_mz_lists, f)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process MGF files into a combined dictionary')
    parser.add_argument('--input_folder', '-i',
                        help='Folder containing the MGF files',
                        default='/home/shipei/projects/revcos/cluster/dataset_clustered/pos')
    parser.add_argument('--out_name', '-o',
                        help='Path for output pickle file',
                        default='/home/shipei/projects/revcos/cluster/datad/pos_all_prec_mz.pkl')

    args = parser.parse_args()

    start_time = time.time()
    main(args.input_folder, args.out_name)
    print(f"Total time: {(time.time() - start_time) / 60:.2f} minutes")
