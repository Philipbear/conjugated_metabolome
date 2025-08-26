from collections import defaultdict
from multiprocessing import Pool

import os
from pathlib import Path
import numpy as np
from tqdm import tqdm


def read_mgf(input_file):
    """Read precursor m/z values and charges from an MGF file."""
    pos_mzs = []
    neg_mzs = []

    with open(input_file, 'r') as mgf_in:
        current_mz = None
        current_charge = None

        for line in mgf_in:
            line = line.strip()

            if line == "BEGIN IONS":
                current_mz = None
                current_charge = None
            elif line == "END IONS":
                if current_mz is not None and current_charge is not None:
                    if current_charge == '1+':
                        pos_mzs.append(current_mz)
                    elif current_charge == '1-':
                        neg_mzs.append(current_mz)
            elif line.startswith("PEPMASS="):
                current_mz = float(line[8:])
            elif line.startswith("CHARGE="):
                current_charge = line[7:]

    return pos_mzs, neg_mzs


def get_all_prec_mzs():
    # Get list of MGF files
    mgfs = [f for f in os.listdir('unclustered_all') if f.endswith('.mgf')]
    mgfs = [os.path.join('unclustered_all', f) for f in mgfs]

    # Use parallel processing with progress bar
    with Pool() as pool:
        results = list(tqdm(pool.imap(read_mgf, mgfs), total=len(mgfs),
                            desc="Processing MGF files"))

    # Combine results from all files
    all_pos_mzs = []
    all_neg_mzs = []
    for pos_mzs, neg_mzs in results:
        all_pos_mzs.extend(pos_mzs)
        all_neg_mzs.extend(neg_mzs)

    # Convert to numpy arrays
    pos_array = np.array(all_pos_mzs)
    neg_array = np.array(all_neg_mzs)

    # Save arrays
    np.save('data/pos_mzs.npy', pos_array)
    np.save('data/neg_mzs.npy', neg_array)

    print(f"Total positive mode spectra: {len(pos_array)}")
    print(f"Total negative mode spectra: {len(neg_array)}")


def get_all_prec_mzs_for_a_dataset(input_dir, all_mgfs, dataset_name):
    # Get list of MGF files
    mgfs = [f for f in all_mgfs if f.split('_')[0] == dataset_name]
    mgfs = [os.path.join(input_dir, f) for f in mgfs]

    # Combine results from all files
    all_pos_mzs = []
    all_neg_mzs = []
    for mgf in mgfs:
        pos_mzs, neg_mzs = read_mgf(mgf)
        all_pos_mzs.extend(pos_mzs)
        all_neg_mzs.extend(neg_mzs)

    # Convert to numpy arrays
    pos_array = np.array(all_pos_mzs)
    neg_array = np.array(all_neg_mzs)

    return pos_array, neg_array


def find_group_id(mz, bin_size=50.0, start_mz=50.0):
    """
    Find the group ID for a given m/z value, handling both small and large bin sizes.
    Uses integer arithmetic for small bins (<1 Da) and direct division for large bins (≥1 Da).

    Parameters:
    -----------
    mz : float
        The m/z value to bin
    bin_size : float, optional
        Size of each bin in Da (default: 0.5)
    start_mz : float, optional
        Starting m/z value (default: 50.0)

    Returns:
    --------
    int
        Bin index
    """
    if bin_size < 1:
        # For small bins, use scaling to avoid floating point arithmetic errors
        scale = int(1 / bin_size)
        bin_index = int((mz * scale) - (start_mz * scale))
    else:
        # For large bins, use direct division
        bin_index = int((mz - start_mz) / bin_size)

    return bin_index


def process_mgf_file(input_file, output_dir):
    """Process a single MGF file and split its spectra into groups."""

    # Initialize dictionaries to store spectra for each group
    pos_spectra = defaultdict(list)
    neg_spectra = defaultdict(list)

    # Read and group spectra as before
    current_spectrum = []
    current_mz = None
    current_charge = None

    with open(input_file, 'r') as mgf_in:
        for line in mgf_in:
            line = line.strip()
            current_spectrum.append(line)

            if line.startswith("PEPMASS="):
                current_mz = float(line[8:])
            elif line.startswith("CHARGE="):
                current_charge = line[7:]
            elif line == "END IONS":
                if current_mz < 50. or np.abs(current_mz - 1000.0000) < 1e-4:
                    current_spectrum = []
                    current_mz = None
                    current_charge = None
                    continue

                if current_charge == '1+':
                    group_id = find_group_id(current_mz)
                    pos_spectra[group_id].extend(current_spectrum)
                    pos_spectra[group_id].append('')
                elif current_charge == '1-':
                    group_id = find_group_id(current_mz)
                    neg_spectra[group_id].extend(current_spectrum)
                    neg_spectra[group_id].append('')

                current_spectrum = []
                current_mz = None
                current_charge = None

    # Write spectra to files
    base_name = os.path.splitext(os.path.basename(input_file))[0].split('_')[0]

    # Create output directories
    os.makedirs(os.path.join(output_dir, 'pos'), exist_ok=True)
    os.makedirs(os.path.join(output_dir, 'neg'), exist_ok=True)

    # Write positive mode spectra
    for group_id, spectra in pos_spectra.items():
        output_file = os.path.join(output_dir, 'pos', f'{base_name}_{group_id}.mgf')
        with open(output_file, 'a') as f:
            f.write('\n'.join(spectra) + '\n')

    # Write negative mode spectra
    for group_id, spectra in neg_spectra.items():
        output_file = os.path.join(output_dir, 'neg', f'{base_name}_{group_id}.mgf')
        with open(output_file, 'a') as f:
            f.write('\n'.join(spectra) + '\n')

    return


def find_split_points(mz_array, chunk_size=100000):
    """
    Find split points in the mz array to create chunks of roughly equal size.
    Returns an array of mz values to use as boundaries between groups.
    """
    if len(mz_array) <= chunk_size:
        return np.array([])

    # Sort the array
    sorted_mz = np.sort(mz_array)

    # Calculate number of chunks needed
    n_chunks = max(1, len(sorted_mz) // chunk_size)

    # Get split points
    split_indices = np.linspace(0, len(sorted_mz), n_chunks + 1)[1:-1]
    split_points = sorted_mz[split_indices.astype(int)]

    # Round to integers
    split_points = np.round(split_points).astype(int)

    # Remove duplicates
    split_points = np.unique(split_points)

    return split_points


def get_mz_group(mz, split_points):
    """
    Determine which group a given m/z value belongs to based on split points.
    """
    if len(split_points) == 0:
        return 0
    return np.searchsorted(split_points, mz)


def process_mgf_file_with_groups(input_file, output_dir, pos_split_points, neg_split_points, dataset_name):
    """Process a single MGF file and split its spectra into groups based on pre-calculated split points."""
    # Initialize dictionaries to store spectra for each group
    pos_spectra = defaultdict(list)
    neg_spectra = defaultdict(list)

    current_spectrum = []
    current_mz = None
    current_charge = None

    with open(input_file, 'r') as mgf_in:
        for line in mgf_in:
            line = line.strip()
            current_spectrum.append(line)

            if line.startswith("PEPMASS="):
                current_mz = float(line[8:])
            elif line.startswith("CHARGE="):
                current_charge = line[7:]
            elif line == "END IONS":
                if current_mz is None or current_mz < 50. or np.abs(current_mz - 1000.0000) < 1e-4:
                    current_spectrum = []
                    current_mz = None
                    current_charge = None
                    continue

                if current_charge == '1+':
                    group_id = get_mz_group(current_mz, pos_split_points)
                    pos_spectra[group_id].extend(current_spectrum)
                    pos_spectra[group_id].append('')
                elif current_charge == '1-':
                    group_id = get_mz_group(current_mz, neg_split_points)
                    neg_spectra[group_id].extend(current_spectrum)
                    neg_spectra[group_id].append('')

                current_spectrum = []
                current_mz = None
                current_charge = None

    # Write positive mode spectra
    os.makedirs(os.path.join(output_dir, 'pos'), exist_ok=True)
    for group_id, spectra in pos_spectra.items():
        output_file = os.path.join(output_dir, 'pos', f'{dataset_name}_{group_id}.mgf')
        with open(output_file, 'a') as f:
            f.write('\n'.join(spectra) + '\n')

    # Write negative mode spectra
    os.makedirs(os.path.join(output_dir, 'neg'), exist_ok=True)
    for group_id, spectra in neg_spectra.items():
        output_file = os.path.join(output_dir, 'neg', f'{dataset_name}_{group_id}.mgf')
        with open(output_file, 'a') as f:
            f.write('\n'.join(spectra) + '\n')


def process_dataset(args):
    """Process a single dataset with all its MGF files."""
    input_dir, output_dir, dataset_name, mgf_files, chunk_size = args

    # Get precursor m/z values for this dataset
    pos_array, neg_array = get_all_prec_mzs_for_a_dataset(input_dir, mgf_files, dataset_name)

    # Find split points for positive and negative mode
    pos_split_points = find_split_points(pos_array, chunk_size)
    neg_split_points = find_split_points(neg_array, chunk_size)

    # Get files for this dataset
    dataset_files = [f for f in mgf_files if f.split('_')[0] == dataset_name]

    # Process each file
    for mgf_file in dataset_files:
        input_file = os.path.join(input_dir, mgf_file)
        process_mgf_file_with_groups(input_file, output_dir, pos_split_points,
                                     neg_split_points, dataset_name)

    return f"Completed processing dataset: {dataset_name}"


def split_mgf_files(input_dir, output_dir, chunk_size=35000):
    """Split all MGF files in input directory using parallel processing at dataset level."""
    print("Initializing parallel processing...")

    os.makedirs(output_dir, exist_ok=True)

    # Get all MGF files and unique datasets
    mgf_files = [f for f in os.listdir(input_dir) if f.endswith('.mgf')]
    all_datasets = list(set(f.split('_')[0] for f in mgf_files))

    # Prepare arguments for parallel processing
    process_args = [(input_dir, output_dir, dataset_name, mgf_files, chunk_size)
                    for dataset_name in all_datasets]

    # Use parallel processing with progress bar
    with Pool() as pool:
        list(tqdm(
            pool.imap(process_dataset, process_args),
            total=len(all_datasets),
            desc="Processing datasets in parallel"
        ))

    print(f"\nCompleted processing {len(all_datasets)} datasets")


import math


def split_large_mgf_files(input_folder: str, output_folder: str, spectra_per_file: int = 120000) -> None:
    """
    Split large MGF files into multiple files, each containing a specified number of spectra.

    Args:
        input_folder (str): Path to folder containing MGF files
        output_folder (str): Path to folder where split files will be saved
        spectra_per_file (int): Number of spectra per output file (default: 120000)
    """
    # Create output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)

    # Process all MGF files in input folder
    for mgf_file in Path(input_folder).glob('*.mgf'):
        # Initialize variables
        spectrum_count = 0
        current_spectrum = []
        current_file_index = 1
        current_output_file = None

        # First count total spectra to know how many files we'll create
        total_spectra = 0
        with open(mgf_file, 'r') as input_file:
            for line in input_file:
                if line.strip() == "END IONS":
                    total_spectra += 1

        total_files = math.ceil(total_spectra / spectra_per_file)
        print(f"\nProcessing {mgf_file.name}:")
        print(f"Total spectra: {total_spectra}")
        print(f"Will create {total_files} files with {spectra_per_file} spectra each\n")

        # Now split the file
        current_output_file = open(os.path.join(output_folder, f"{mgf_file.stem}_part{current_file_index}.mgf"), 'w')

        with open(mgf_file, 'r') as input_file:
            for line in input_file:
                current_spectrum.append(line)

                # When we hit "END IONS", we've completed a spectrum
                if line.strip() == "END IONS":
                    spectrum_count += 1

                    # Write the spectrum to current file
                    current_output_file.writelines(current_spectrum)

                    # If we've reached the spectra limit for this file, start a new file
                    if spectrum_count % spectra_per_file == 0 and spectrum_count < total_spectra:
                        current_output_file.close()
                        current_file_index += 1
                        current_output_file = open(os.path.join(output_folder,
                                                                f"{mgf_file.stem}_part{current_file_index}.mgf"), 'w')
                        print(f"Created part{current_file_index - 1}.mgf with {spectra_per_file} spectra")

                    # Reset current spectrum buffer
                    current_spectrum = []

            # Close the last file
            if current_output_file:
                current_output_file.close()

        # Print final stats
        final_file_spectra = total_spectra % spectra_per_file
        if final_file_spectra > 0:
            print(f"Created part{total_files}.mgf with {final_file_spectra} spectra")
        print(f"\nCompleted processing {mgf_file.name}")
        print(f"Total spectra processed: {spectrum_count}")



if __name__ == "__main__":
    # get_all_prec_mzs()  # server
    '''
    Total positive mode spectra: 429964017
    Total negative mode spectra: 109240316 
    '''

    # split_mgf_files('unclustered_all', 'unclustered_mzgrouped')  # server

    split_large_mgf_files('/Users/shipei/Documents/projects/chemical_conjugate_discovery/clustering/data/tmp_data/neg_unclustered',
                          '/Users/shipei/Documents/projects/chemical_conjugate_discovery/clustering/data/tmp_data/neg_unclustered_small')
