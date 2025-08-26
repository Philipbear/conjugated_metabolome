import os
import time
import argparse
import tqdm
import tempfile
import shutil
import pickle
import pandas as pd


def process_file(file_path, mgf_name, work_dir, remove_unclustered=False, keep_csv=False):

    # clean work space
    os.system(f'rm -rf {work_dir}/nn/ 2>/dev/null; rm -rf {work_dir}/spectra/ 2>/dev/null')

    # main command
    cmd = f'falcon {file_path} {mgf_name}_clustered --work_dir {work_dir}\
    --export_representatives --export_include_singletons \
    --usi_pxd {mgf_name} --precursor_tol 20 ppm --fragment_tol 0.05 --eps 0.05 \
    --min_peaks 4 --min_mz_range 0 --min_mz 0 --max_mz 2000 --min_intensity 0.01 --max_peaks_used 50'

    exit_code = os.system(cmd)
    if exit_code != 0:
        print(f"Error: {mgf_name}.mgf")

        # shutil.copy(file_path, work_dir)
    else:
        print(f"Success: {mgf_name}.mgf")

        # clean work space
        cmd_2 = f'rm -rf {work_dir}/nn/ 2>/dev/null; rm -rf {work_dir}/spectra/ 2>/dev/null'
        os.system(cmd_2)

        if keep_csv:
            cmd_3 = f'mv {mgf_name}_clustered.csv {work_dir}/ 2>/dev/null; mv {mgf_name}_clustered.mgf {work_dir}/ 2>/dev/null'
            os.system(cmd_3)

            new_csv_path = os.path.join(work_dir, f'{mgf_name}_clustered.csv')
            if not os.path.exists(new_csv_path):
                print(f"Files not found: {new_csv_path}")
                return
        else:
            cmd_3 = f'rm {mgf_name}_clustered.csv 2>/dev/null; mv {mgf_name}_clustered.mgf {work_dir}/ 2>/dev/null'
            os.system(cmd_3)

        print('Refining files...')
        new_mgf_path = os.path.join(work_dir, f'{mgf_name}_clustered.mgf')

        if not os.path.exists(new_mgf_path):
            print(f"Files not found: {new_mgf_path}")
            return

        metadata_dict = gen_metadata_from_mgf(file_path, mgf_name)

        # modify the newly generated mgf file
        all_title_list = modify_and_write_mgf(new_mgf_path, metadata_dict)

        if keep_csv:
            modify_and_write_pkl(new_csv_path, all_title_list, metadata_dict)

        if remove_unclustered:
            os.remove(file_path)


def gen_metadata_from_mgf(mgf_path, mgf_name):
    """
    Generate a dataframe from mgf files
    """
    all_out_dict = {}

    with open(mgf_path, 'r') as file:
        out_dict = {}
        index = 0
        for line in file:
            # empty line
            _line = line.strip()
            if not _line:
                continue
            elif line.startswith('BEGIN IONS'):
                spectrum = {}
            elif line.startswith('END IONS'):
                key_name = mgf_name + f':index:{index}'
                out_dict[key_name] = spectrum['TITLE']
                index += 1
                continue
            else:
                # if line contains '=', it is a key-value pair
                if '=' in _line:
                    # split by first '='
                    key, value = _line.split('=', 1)
                    spectrum[key] = value
                else:
                    # if no '=', it is a spectrum pair
                    continue
    all_out_dict.update(out_dict)

    return all_out_dict


def format_peak_line(line):
    parts = line.strip().split()
    if len(parts) == 2:
        mz, intensity = parts
        try:
            formatted_mz = f"{float(mz):.4f}"
            formatted_intensity = f"{int(float(intensity))}"
            return f"{formatted_mz} {formatted_intensity}\n"
        except ValueError:
            return line  # Return original line if conversion fails
    return line  # Return original line if it doesn't have exactly two parts


def modify_and_write_mgf(mgf, metadata_dict):
    """
    return all_title_ls: list of all titles in the mgf file
    """

    all_title_ls = []
    temp_file = tempfile.NamedTemporaryFile(mode='w+', delete=False)
    try:
        with open(mgf, 'r') as infile:
            in_peaks = False
            for line in infile:
                if line.strip() == "BEGIN IONS":
                    in_peaks = False
                elif line.strip() == "END IONS":
                    in_peaks = False

                if line.startswith('TITLE'):
                    try:
                        key_name = ":".join(line.strip().rsplit(":", 3)[-3:])
                        title = metadata_dict[key_name]
                        temp_file.write(f'TITLE={title}\n')
                        all_title_ls.append(title)
                    except (AttributeError, ValueError, KeyError) as e:
                        print(f"Error processing TITLE line: {e}")
                        temp_file.write(line)  # Keep original line if there's an error
                elif line.startswith('RTINSECONDS') or line.startswith('CLUSTER'):
                    continue  # Skip these lines
                elif in_peaks and not line.startswith('END IONS'):
                    temp_file.write(format_peak_line(line))
                else:
                    temp_file.write(line)

                if line.startswith('CHARGE'):
                    in_peaks = True

        temp_file.close()

        # Replace the original file with the temporary file
        shutil.move(temp_file.name, mgf)
    except Exception as e:
        print(f"An error occurred: {e}")
        os.unlink(temp_file.name)  # Delete the temporary file if an error occurs
        raise

    return all_title_ls


def modify_and_write_pkl(csv, all_title_list, metadata_dict):
    df = pd.read_csv(csv, comment='#')

    df = df[df['cluster'] != -1].reset_index(drop=True)

    # drop cols of precursor_charge, retention_time, precursor_mz
    df = df.drop(['precursor_charge', 'retention_time', 'precursor_mz'], axis=1)

    df['new_title'] = df['identifier'].apply(lambda x: ":".join(x.rsplit(":", 3)[-3:]))
    df['new_title'] = df['new_title'].apply(lambda x: metadata_dict[x])

    df['representative'] = df['new_title'].apply(lambda x: x in all_title_list)

    # dict to store the representative: [list of cluster members]
    out_dict = {}

    # Group by cluster
    for cluster_id, group in df.groupby('cluster'):
        # Find the representative for this cluster
        rep_row = group[group['representative']].iloc[0]
        rep_title = rep_row['new_title']

        # Get all member titles except the representative
        member_titles = group['new_title'].tolist()

        # Store in dictionary
        out_dict[rep_title] = member_titles

    # remove the old csv file
    os.remove(csv)

    # save out_dict as pickle
    pkl_name = csv.replace('.csv', '.pkl')
    with open(pkl_name, 'wb') as f:
        pickle.dump(out_dict, f)


def main(file_folder, output_folder, remove_unclustered=False):

    os.makedirs(output_folder, exist_ok=True)

    start_time = time.time()

    all_mgf_files = os.listdir(file_folder)
    all_mgf_files = [x for x in all_mgf_files if x.endswith('.mgf') and not x.startswith('.')]
    all_mgf_files = [os.path.join(file_folder, x) for x in all_mgf_files]
    all_mgf_files = [x for x in all_mgf_files if os.path.getsize(x) > 1]

    for mgf_file in tqdm.tqdm(all_mgf_files, desc='Processing files', total=len(all_mgf_files)):

        base_name = os.path.splitext(os.path.basename(mgf_file))[0]

        # if file is in the output folder, skip
        if f'{base_name}_clustered.mgf' in os.listdir(output_folder):
            print(f"File already processed: {base_name}.mgf")
            continue

        print('='*50)
        print(f"Processing file: {base_name}.mgf")
        process_file(mgf_file, base_name, output_folder, remove_unclustered)

    # print time in minutes
    print(f"Total time: {(time.time() - start_time) / 60} minutes")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run falcon for clustering')
    parser.add_argument('--input_folder', '-i', help='Folder containing the MGF files',
                        default='unclustered_all')
    parser.add_argument('--output_folder', '-o', help='Folder to store the results',
                        default='clustered_all')
    parser.add_argument('--remove_unclustered', '-rm',
                        action='store_true', help='Remove unclustered MGF files')
    args = parser.parse_args()

    main(args.input_folder, args.output_folder, args.remove_unclustered)

    # main('/Users/shipei/Documents/projects/chemical_conjugate_discovery/clustering/data/MTBLS/test',
    #      '/Users/shipei/Documents/projects/chemical_conjugate_discovery/clustering/data/MTBLS/test_out')
