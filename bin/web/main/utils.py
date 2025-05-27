import pandas as pd
import urllib.parse

'''
# ref_1_id, ref_2_id, delta_mass, count, dataset, file_path, file_scan, full_file_path
ms2db_df: name, 'db', 'db_id', 'inchikey_14', monoisotopic_mass'
'''

def filter_search_results(df, ms2db_df, inchikey, mono_mass, min_count=3, match_filter_ls=['spec', 'delta']):
    """
    Filter the DataFrame by 2D InChIKey and monoisotopic mass.
    """
    
    dbid_to_inchi_dict = ms2db_df.set_index('db_id')['inchikey_14'].to_dict()
    
    result_dfs = []
    if 'spec' in match_filter_ls:
        # from ms2db_df, get rows with given inchikey
        ms2db_filtered = ms2db_df[(ms2db_df['inchikey_14'] == inchikey)].reset_index(drop=True)
        if not ms2db_filtered.empty:
            # Get db_id from the filtered ms2db_df
            db_ids = ms2db_filtered['db_id'].tolist()
            del ms2db_filtered
        
            df_filtered1 = df[(df['ref_1_id'].isin(db_ids)) & (df['count'] >= min_count)].reset_index(drop=True)
            df_filtered1['Match type'] = 'Spectral match (Ref 1)'
            df_filtered2 = df[(df['ref_2_id'].isin(db_ids)) & (df['count'] >= min_count)].reset_index(drop=True)
            df_filtered2['Match type'] = 'Spectral match (Ref 2)'
            
            
            if not df_filtered1.empty:
                df_filtered1['conjugate_inchikey'] = df_filtered1['ref_2_id'].apply(lambda x: dbid_to_inchi_dict.get(x, None))
                
                result_dfs.append(df_filtered1)
            
            if not df_filtered2.empty:
                df_filtered2['conjugate_inchikey'] = df_filtered2['ref_1_id'].apply(lambda x: dbid_to_inchi_dict.get(x, None))
                
                db_to_mass_dict = ms2db_df.set_index('db_id')['monoisotopic_mass'].to_dict()
                df_filtered2['delta_mass'] = df_filtered2['ref_1_id'].apply(lambda x: round(db_to_mass_dict.get(x, 0) - 18.0106, 2))
                del db_to_mass_dict
                
                result_dfs.append(df_filtered2)
    
    if 'delta' in match_filter_ls:
        target_mass = round(mono_mass - 18.0106, 2)
        # Filter by delta mass
        df_filtered3 = df[(df['delta_mass'] == target_mass) & (df['count'] >= min_count)].reset_index(drop=True)
        if not df_filtered3.empty:
            df_filtered3['Match type'] = 'Delta mass'
            # Add conjugate inchikey
            df_filtered3['conjugate_inchikey'] = df_filtered3['ref_1_id'].apply(lambda x: dbid_to_inchi_dict.get(x, None))
            del dbid_to_inchi_dict
            # Add delta mass results
            result_dfs.append(df_filtered3)
    
    if not result_dfs:
        return pd.DataFrame()
    
    # Concatenate the filtered DataFrames
    df_filtered = pd.concat(result_dfs, ignore_index=True)
    del result_dfs
    
    # fill in conjugate names
    inchikey_to_name = ms2db_df.set_index('inchikey_14')['name'].to_dict()
    df_filtered['Conjugate name'] = df_filtered['conjugate_inchikey'].apply(
        lambda x: inchikey_to_name.get(x, None) if pd.notna(x) else None
    )
    # drop col
    df_filtered = df_filtered.drop(columns=['conjugate_inchikey'])
    del inchikey_to_name
    
    def _get_qry_usi(row):
        return 'mzspec:' + row['dataset'] + ':' + row['full_file_path'] + ':scan:' + row['file_scan']
    
    def _get_ref_usi(row, ref_col='ref_1'):
        if pd.isna(row[f'{ref_col}_id']):
            return None
        if row[f'{ref_col}_db'] == 'gnps':
            return 'mzspec:GNPS:GNPS-LIBRARY:accession:CCMSLIB000' + row[f'{ref_col}_id']
        elif row[f'{ref_col}_db'] == 'mona':
            return 'mzspec:MASSBANK::accession:' + row[f'{ref_col}_id']
        else:
            return 'No valid USI'
    
    df_filtered['qry_usi'] = df_filtered.apply(lambda x: _get_qry_usi(x), axis=1)
    df_filtered = df_filtered.drop(columns=['dataset', 'file_scan', 'full_file_path'])
    
    # Add names for ref_1 and ref_2
    db_dict = ms2db_df.set_index('db_id')['db'].to_dict()
    df_filtered['ref_1_db'] = df_filtered['ref_1_id'].apply(lambda x: db_dict.get(x, ''))
    df_filtered['ref_2_db'] = df_filtered['ref_2_id'].apply(lambda x: db_dict.get(x, ''))
    del db_dict
    
    df_filtered['usi1'] = df_filtered.apply(lambda x: _get_ref_usi(x, 'ref_1'), axis=1)
    df_filtered['usi2'] = df_filtered.apply(lambda x: _get_ref_usi(x, 'ref_2'), axis=1)
    df_filtered = df_filtered.drop(columns=['ref_1_id', 'ref_2_id', 'ref_1_db', 'ref_2_db'])
    
    # sort by count
    df_filtered = df_filtered.sort_values(by='count', ascending=False).reset_index(drop=True)
    
    return df_filtered


def add_mirror_plot_urls(df):
    """
    Add mirror plot URLs to the DataFrame.
    """
    # Initialize columns for URLs and display text
    df['Mirror plot (Ref 1)'] = df.apply(lambda x: gen_mirror_plot_url(x['qry_usi'], x['usi1']), axis=1)
    df['Mirror plot (Ref 2)'] = df.apply(lambda x: gen_mirror_plot_url(x['qry_usi'], x['usi2']), axis=1)
    
    return df


def gen_mirror_plot_url(usi1, usi2):
    """
    Generate a URL for the metabolomics-usi.gnps2.org mirror plot.
    
    Parameters:
    -----------
    usi1 : str
        First USI for the mirror plot (top spectrum)
    usi2 : str, optional
        Second USI for the mirror plot (bottom spectrum). If None, only usi1 will be displayed.
    width : float
        Width of the plot in inches
    height : float
        Height of the plot in inches
    mz_min : float or None
        Minimum m/z to display
    mz_max : float or None
        Maximum m/z to display
    max_intensity : int
        Maximum intensity as percentage
        
    Returns:
    --------
    str : URL for the mirror plot
    """
    
    if usi2 is None:
        return None

    # Create the base URL
    base_url = "https://metabolomics-usi.gnps2.org/dashinterface/"
    
    # Define parameters
    params = {
        "usi1": usi1,
        "usi2": usi2,
        # "width": 10.0,
        # "height": 6.0,
        # "mz_min": None,
        # "mz_max": None,
        "max_intensity": 150,
        # "annotate_precision": 4,
        # "annotation_rotation": 90,
        "cosine": "shifted",
        "fragment_mz_tolerance": 0.05,
        # "grid": True,

        # "annotate_peaks": "[[],[]]"
    }
    if usi2 == 'No valid USI':
        # If usi2 is not valid, remove it from the parameters
        params.pop("usi2")
    
    # URL encode the parameters
    query_string = urllib.parse.urlencode(params)
    
    # Create the full URL
    url = f"{base_url}?{query_string}"
    
    return url


if __name__ == "__main__":
    
    # df = pd.read_parquet("/Users/shipei/Documents/projects/conjugated_metabolome/bin/web/main/neg_web_full_path.parquet")
    # print(df.shape)
    
    # Example usage
    usi1 = "mzspec:GNPS:GNPS-LIBRARY:accession:CCMSLIB00005436077"
    usi2 = "mzspec:GNPS:GNPS-LIBRARY:accession:CCMSLIB00005436078"
    url = gen_mirror_plot_url(usi1, usi2)
    print(url)  # This will print the generated URL for the mirror plot
