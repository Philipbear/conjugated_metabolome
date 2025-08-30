"""
after clustering, how many spectra?
"""

import os
import time
from tqdm import tqdm


def read_spec(file_obj):
    """Read a single spectrum block from an open MGF file and return 1 if valid.

    Args:
        file_obj: An opened file object positioned at the start of a spectrum

    Returns:
        int: 1 if a valid spectrum was found, 0 if end of file is reached
    """
    # Skip any empty lines before BEGIN IONS
    for line in file_obj:
        if line.strip() == 'BEGIN IONS':
            break
    else:  # EOF reached
        return 0

    has_peaks = False

    # Read until END IONS
    for line in file_obj:
        line = line.strip()

        if line == 'END IONS':
            return 1 if has_peaks else 0  # Only count if we have peaks

        # Check if this is a peak line (contains two float values)
        if line and not line.startswith(('BEGIN', 'END', 'TITLE=', 'PEPMASS=', 'CHARGE=')):
            try:
                # Just check if it's a valid peak line (two floats)
                mz, intensity = line.split()
                float(mz)
                float(intensity)
                has_peaks = True
            except (ValueError, IndexError):
                # Not a valid peak line, ignore
                pass

    return 0  # EOF reached without END IONS


def process_mgf(mgf_path, buffer_size=8192):
    """Count spectra in an MGF file efficiently.

    Args:
        mgf_path: Path to the MGF file
        buffer_size: Read buffer size in bytes

    Returns:
        int: Number of valid spectra in the file
    """
    count = 0
    with open(mgf_path, 'r', buffering=buffer_size) as f:
        while True:
            result = read_spec(f)
            if result == 0:
                break
            count += 1
    return count


def main(folder):
    """
    Process all MGF files in a folder and print results
    """
    # Get list of MGF files
    all_mgfs = [os.path.join(folder, x) for x in os.listdir(folder)
                if x.endswith('.mgf') and not x.startswith('.')]

    # Process files
    total = 0
    for mgf_file in tqdm(all_mgfs, total=len(all_mgfs), desc='Processing MGFs'):
        result = process_mgf(mgf_file)
        total += result

    print(f"Total spectra in {folder}: {total}")


if __name__ == '__main__':

    start_time = time.time()
    main('/home/shipei/projects/revcos/cluster/dataset_clustered/pos')
    print(f"Total time: {(time.time() - start_time) / 60:.2f} minutes")

    start_time = time.time()
    main('/home/shipei/projects/revcos/cluster/dataset_clustered/neg')
    print(f"Total time: {(time.time() - start_time) / 60:.2f} minutes")

    """
    Processing MGFs: 100%|██████████████████████████████████████████| 4336/4336 [25:17<00:00,  2.86it/s]
    Total spectra in /home/shipei/projects/revcos/cluster/dataset_clustered/pos: 118764042
    Total time: 25.29 minutes
    Processing MGFs: 100%|██████████████████████████████████████████| 2675/2675 [06:02<00:00,  7.37it/s]
    Total spectra in /home/shipei/projects/revcos/cluster/dataset_clustered/neg: 31181284
    Total time: 6.05 minutes
    """
