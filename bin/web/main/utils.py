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
            df_filtered1['Match type'] = 'spec (ref 1)'
            df_filtered2 = df[(df['ref_2_id'].isin(db_ids)) & (df['count'] >= min_count)].reset_index(drop=True)
            df_filtered2['Match type'] = 'spec (ref 2)'
            
            
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
        df_filtered3 = df[(abs(df['delta_mass'] - target_mass) <= 0.01)& (df['count'] >= min_count)].reset_index(drop=True)
        if not df_filtered3.empty:
            df_filtered3['Match type'] = 'delta'
            
            # Add conjugate inchikey
            df_filtered3['conjugate_inchikey'] = df_filtered3['ref_1_id'].apply(lambda x: dbid_to_inchi_dict.get(x, None))
            del dbid_to_inchi_dict
            
            db_to_mass_dict = ms2db_df.set_index('db_id')['monoisotopic_mass'].to_dict()
            df_filtered3['delta_mass'] = df_filtered3['ref_1_id'].apply(lambda x: round(db_to_mass_dict.get(x, 0) - 18.0106, 2))
            del db_to_mass_dict
            
            # Add results
            result_dfs.append(df_filtered3)
    
    if not result_dfs:
        return pd.DataFrame()
    
    # Concatenate the filtered DataFrames
    df_filtered = pd.concat(result_dfs, ignore_index=True)
    del result_dfs
    
    df_filtered['Annotation type'] = df_filtered['ref_2_id'].apply(
        lambda x: 'spec_spec' if pd.notna(x) else 'spec_delta'
    )
    
    # fill in conjugate names
    inchikey_to_name = ms2db_df.set_index('inchikey_14')['name'].to_dict()
    df_filtered['Conjugate name'] = df_filtered['conjugate_inchikey'].apply(
        lambda x: inchikey_to_name.get(x, None) if pd.notna(x) else None
    )
    # drop col
    df_filtered = df_filtered.drop(columns=['conjugate_inchikey'])
    del inchikey_to_name
    
    def _get_qry_usi(row):
        try:
            return 'mzspec:' + row['dataset'] + ':' + row['full_file_path'] + ':scan:' + row['file_scan']
        except:
            return None
    
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


def prepare_delta_mass_plot(df_filtered):
    """
    Efficiently prepare data for delta mass distribution plot.
    Groups data by delta mass and ion polarity, combining entries that appear in both modes.
    """
    # Group by delta mass and ion polarity, summing the counts
    delta_mass_counts = df_filtered.groupby(['Conjugate delta mass', 'Ion polarity'], as_index=False)['Count'].sum()
    
    # Split by ion polarity
    pos_data = delta_mass_counts[delta_mass_counts['Ion polarity'] == '+']
    neg_data = delta_mass_counts[delta_mass_counts['Ion polarity'] == '-']
    
    # Get unique delta masses for each polarity
    pos_deltas = set(pos_data['Conjugate delta mass'])
    neg_deltas = set(neg_data['Conjugate delta mass'])
    
    # Find delta masses that appear in both polarities
    common_deltas = pos_deltas.intersection(neg_deltas)
    
    # For delta masses in both polarities, create a single 'both' entry
    if common_deltas:
        # Create a list to hold our processed data
        result_rows = []
        
        # Process entries unique to positive mode
        unique_pos = pos_data[~pos_data['Conjugate delta mass'].isin(common_deltas)]
        result_rows.extend(unique_pos.to_dict('records'))
        
        # Process entries unique to negative mode
        unique_neg = neg_data[~neg_data['Conjugate delta mass'].isin(common_deltas)]
        result_rows.extend(unique_neg.to_dict('records'))
        
        # Process common entries - combine into 'both' entries using groupby
        # First, filter the original dataframe to only include rows with common delta masses
        common_data = delta_mass_counts[delta_mass_counts['Conjugate delta mass'].isin(common_deltas)]
        
        # Group by delta mass and sum the counts
        combined_counts = common_data.groupby('Conjugate delta mass', as_index=False)['Count'].sum()
        
        # Add the 'both' polarity to each row
        combined_counts['Ion polarity'] = 'both'
        
        # Add the combined rows to our result
        result_rows.extend(combined_counts.to_dict('records'))
        
        # Create a new DataFrame from the processed rows
        delta_mass_counts = pd.DataFrame(result_rows)
    
    # Sort by delta mass for consistent display
    delta_mass_counts = delta_mass_counts.sort_values('Conjugate delta mass')
    
    return delta_mass_counts


def add_urls(df):
    """
    Add mirror plot URLs to the DataFrame.
    """
    # Initialize columns for URLs and display text
    df['Mirror plot (Ref 1)'] = df.apply(lambda x: gen_mirror_plot_url(x['qry_usi'], x['usi1']), axis=1)
    df['Mirror plot (Ref 2)'] = df.apply(lambda x: gen_mirror_plot_url(x['qry_usi'], x['usi2']), axis=1)
    df['MASST'] = df.apply(lambda x: gen_fasst_url(x['qry_usi']), axis=1)
    
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
    
    if usi1 is None or usi2 is None:
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


def gen_fasst_url(usi):
    """
    Generate a URL for fasst MASST
    """
    

    # Create the base URL
    base_url = "https://fasst.gnps2.org/fastsearch/"
    
    # Define parameters
    params = {
        "library_select": "metabolomicspanrepo_index_latest",
        "usi1": usi
    }
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
