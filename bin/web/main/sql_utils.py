import sqlite3
import pandas as pd
import urllib.parse


def get_git_short_rev():
    try:
        with open('.git/logs/HEAD', 'r') as f:
            last_line = f.readlines()[-1]
            hash_val = last_line.split()[1]
        return hash_val[:7]
    except Exception:
        return ".git/ not found"


def filter_search_results(db_path, inchikey_14, mono_mass, min_count=1):
    """
    Main function to filter search results by InChIKey and monoisotopic mass
    """
    conn = create_production_connection(db_path)

    # Search by InChIKey
    df_spec = search_by_inchikey(conn, inchikey_14, min_count)

    # Search by delta mass (monoisotopic mass - water loss)
    target_mass = round(mono_mass - 18.0106, 2)
    df_mass = search_by_delta_mass(conn, target_mass, min_count)
    conn.close()
    
    # Combine results and drop duplicates
    df = pd.concat([df_spec, df_mass], ignore_index=True)
    
    # If still empty, return empty DataFrame
    if df.empty:
        return df
    
    # add annotation type column
    df['annotation_type'] = df['ref_2_id'].apply(
        lambda x: 'spec_spec' if pd.notna(x) else 'spec_delta'
    )
    
    # add USI columns
    df['qry_id'] = 'mzspec:' + df['qry_id']
    df['ref_1_id'] = df.apply(_get_ref_usi, axis=1, ref_col='ref_1')
    df['ref_2_id'] = df.apply(_get_ref_usi, axis=1, ref_col='ref_2')
    
    df['mirror_plot_ref_1'] = df.apply(
        lambda row: gen_mirror_plot_url(row['qry_id'], row['ref_1_id']), axis=1
    )
    df['mirror_plot_ref_2'] = df.apply(
        lambda row: gen_mirror_plot_url(row['qry_id'], row['ref_2_id']), axis=1
    )
    df['masst'] = df['qry_id'].apply(gen_fasst_url)

    return df


def prepare_delta_mass_plot(df_filtered):
    """
    Efficiently prepare data for delta mass distribution plot.
    Groups data by delta mass and ion polarity, combining entries that appear in both modes.
    """
    # Group by delta mass and ion polarity, summing the counts
    delta_mass_counts = df_filtered.groupby(['conjugate_delta_mass', 'ion_polarity'], as_index=False)['count'].sum()
    
    # Split by ion polarity
    pos_data = delta_mass_counts[delta_mass_counts['ion_polarity'] == '+']
    neg_data = delta_mass_counts[delta_mass_counts['ion_polarity'] == '-']

    # Get unique delta masses for each polarity
    pos_deltas = set(pos_data['conjugate_delta_mass'])
    neg_deltas = set(neg_data['conjugate_delta_mass'])
    
    # Find delta masses that appear in both polarities
    common_deltas = pos_deltas.intersection(neg_deltas)
    
    # For delta masses in both polarities, create a single 'both' entry
    if common_deltas:
        # Create a list to hold our processed data
        result_rows = []
        
        # Process entries unique to positive mode
        unique_pos = pos_data[~pos_data['conjugate_delta_mass'].isin(common_deltas)]
        result_rows.extend(unique_pos.to_dict('records'))
        
        # Process entries unique to negative mode
        unique_neg = neg_data[~neg_data['conjugate_delta_mass'].isin(common_deltas)]
        result_rows.extend(unique_neg.to_dict('records'))
        
        # Process common entries - combine into 'both' entries using groupby
        # First, filter the original dataframe to only include rows with common delta masses
        common_data = delta_mass_counts[delta_mass_counts['conjugate_delta_mass'].isin(common_deltas)]
        
        # Group by delta mass and sum the counts
        combined_counts = common_data.groupby('conjugate_delta_mass', as_index=False)['count'].sum()

        # Add the 'both' polarity to each row
        combined_counts['ion_polarity'] = 'both'
        
        # Add the combined rows to our result
        result_rows.extend(combined_counts.to_dict('records'))
        
        # Create a new DataFrame from the processed rows
        delta_mass_counts = pd.DataFrame(result_rows)
    
    # Sort by delta mass for consistent display
    delta_mass_counts = delta_mass_counts.sort_values('conjugate_delta_mass')

    return delta_mass_counts


def _get_ref_usi(row, ref_col='ref_1'):
    if pd.isna(row[f'{ref_col}_id']):
        return None
    if row[f'{ref_col}_db'] == 'gnps':
        return 'mzspec:GNPS:GNPS-LIBRARY:accession:CCMSLIB000' + row[f'{ref_col}_id']
    elif row[f'{ref_col}_db'] == 'mona':
        return 'mzspec:MASSBANK::accession:' + row[f'{ref_col}_id']
    else:
        return 'No valid USI'


def create_production_connection(db_path):
    """
    Create optimized connection for production queries
    """
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    # Optimal settings for <1GB RAM
    cursor.execute("PRAGMA cache_size=10000")
    cursor.execute("PRAGMA mmap_size=134217728")     # 128MB (your optimal setting)
    cursor.execute("PRAGMA temp_store=memory")
    cursor.execute("PRAGMA query_only=1")            # Read-only safety
    
    return conn


def search_by_inchikey(conn, target_inchikey, min_count):
    """
    Search for metabolites by InChIKey with optimized performance
    Returns results with all relevant database IDs and metadata
    """
    
    # Optimized query using UNION for better index usage
    query = '''
        SELECT DISTINCT
            dl.name as qry_dataset, 
            fl.name as qry_file_name, 
            sr.qry_scan_id, 
            sr.qry_mz,
            sr.ref_1_id, 
            sr.ref_2_id, 
            sr.delta_mass / 100.0 as delta_mass,
            sr.count, 
            CASE sr.ion_polarity WHEN 1 THEN '+' ELSE '-' END as ion_polarity,
            -- Reference 1 info
            ms1.inchikey_14 as ref_1_inchikey, 
            c1.name as ref_1_name,
            c1.monoisotopic_mass as ref_1_mono_mass,
            ms1.db as ref_1_db,
            -- Reference 2 info  
            ms2.inchikey_14 as ref_2_inchikey, 
            c2.name as ref_2_name,
            c2.monoisotopic_mass as ref_2_mono_mass,
            ms2.db as ref_2_db
        FROM (
            -- Query for ref_1 matches
            SELECT sr.*
            FROM search_results sr
            JOIN ms2db ms1 ON sr.ref_1_id = ms1.db_id
            JOIN chem_db c1 ON ms1.inchikey_14 = c1.inchikey_14
            WHERE c1.inchikey_14 = ?
            
            UNION
            
            -- Query for ref_2 matches  
            SELECT sr.*
            FROM search_results sr
            JOIN ms2db ms2 ON sr.ref_2_id = ms2.db_id
            JOIN chem_db c2 ON ms2.inchikey_14 = c2.inchikey_14
            WHERE c2.inchikey_14 = ?
        ) sr
        LEFT JOIN dataset_lookup dl ON sr.qry_dataset = dl.id
        LEFT JOIN filename_lookup fl ON sr.qry_file_name = fl.id
        LEFT JOIN ms2db ms1 ON sr.ref_1_id = ms1.db_id
        LEFT JOIN chem_db c1 ON ms1.inchikey_14 = c1.inchikey_14
        LEFT JOIN db_lookup dbl1 ON ms1.db = dbl1.id
        LEFT JOIN ms2db ms2 ON sr.ref_2_id = ms2.db_id
        LEFT JOIN chem_db c2 ON ms2.inchikey_14 = c2.inchikey_14
        LEFT JOIN db_lookup dbl2 ON ms2.db = dbl2.id
        ORDER BY sr.count DESC
    '''
    
    df = pd.read_sql_query(query, conn, params=[target_inchikey, target_inchikey])
    
    if df.empty:
        return df
    
    # Filter by minimum count
    df = df[df['count'] >= min_count].reset_index(drop=True)
    if df.empty:
        return df

    # Reconstruct qry_id
    df['qry_id'] = df.apply(lambda row: reconstruct_qry_id(
        row['qry_dataset'], row['qry_file_name'], row['qry_scan_id']
    ), axis=1)
    
    # Add match type indicator
    def get_match_type_and_conjugate_info(row):
        ref1_match = row['ref_1_inchikey'] == target_inchikey
        ref2_match = row['ref_2_inchikey'] == target_inchikey
        
        if ref1_match:
            return 'spec (ref 1)', row['ref_2_name'], row['delta_mass']
        elif ref2_match:
            return 'spec (ref 2)', row['ref_1_name'], round(row['ref_1_mono_mass'] - 18.0106, 2)
        else:
            raise ValueError("InChIKey not found in either reference")

    df['match_type'], df['conjugate_name'], df['conjugate_delta_mass'] = zip(*df.apply(get_match_type_and_conjugate_info, axis=1))

    # Clean column order with all database info including mono masses
    column_order = [
        'qry_id', 'qry_mz', 'ion_polarity',
        'ref_1_id', 'ref_1_name', 'ref_1_mono_mass', 'ref_1_db',
        'ref_2_id', 'ref_2_name', 'ref_2_mono_mass', 'ref_2_db',
        'match_type', 'delta_mass', 'count', 'conjugate_name', 'conjugate_delta_mass'
    ]
    
    return df[column_order]


def search_by_delta_mass(conn, target_delta_mass, min_count, mass_tolerance=0.01):
    """
    Search for metabolites by delta mass with specified tolerance
    Returns results with all relevant database IDs and metadata
    
    Args:
        db_path: Path to SQLite database
        target_delta_mass: Target delta mass (in Da)
        mass_tolerance: Mass tolerance (default 0.01 Da)
    """
    
    # Convert target delta mass to stored format (multiplied by 100)
    target_delta_stored = round(target_delta_mass * 100)
    tolerance_stored = round(mass_tolerance * 100)
    
    # Single query with tolerance range - more efficient than multiple queries
    query = '''
        SELECT DISTINCT
            dl.name as qry_dataset, 
            fl.name as qry_file_name, 
            sr.qry_scan_id, 
            sr.qry_mz,
            sr.ref_1_id, 
            sr.ref_2_id, 
            sr.delta_mass / 100.0 as delta_mass,
            sr.count, 
            CASE sr.ion_polarity WHEN 1 THEN '+' ELSE '-' END as ion_polarity,
            -- Reference 1 info
            ms1.inchikey_14 as ref_1_inchikey, 
            c1.name as ref_1_name,
            c1.monoisotopic_mass as ref_1_mono_mass,
            ms1.db as ref_1_db,
            -- Reference 2 info  
            ms2.inchikey_14 as ref_2_inchikey, 
            c2.name as ref_2_name,
            c2.monoisotopic_mass as ref_2_mono_mass,
            dbl2.name as ref_2_db
        FROM search_results sr
        LEFT JOIN dataset_lookup dl ON sr.qry_dataset = dl.id
        LEFT JOIN filename_lookup fl ON sr.qry_file_name = fl.id
        LEFT JOIN ms2db ms1 ON sr.ref_1_id = ms1.db_id
        LEFT JOIN chem_db c1 ON ms1.inchikey_14 = c1.inchikey_14
        LEFT JOIN db_lookup dbl1 ON ms1.db = dbl1.id
        LEFT JOIN ms2db ms2 ON sr.ref_2_id = ms2.db_id
        LEFT JOIN chem_db c2 ON ms2.inchikey_14 = c2.inchikey_14
        LEFT JOIN db_lookup dbl2 ON ms2.db = dbl2.id
        WHERE sr.delta_mass BETWEEN ? AND ?
        ORDER BY ABS(sr.delta_mass - ?) ASC, sr.count DESC
    '''
    
    # Calculate range parameters
    min_delta = target_delta_stored - tolerance_stored
    max_delta = target_delta_stored + tolerance_stored
    
    df = pd.read_sql_query(query, conn, params=[min_delta, max_delta, target_delta_stored])
    
    if df.empty:
        return df
    
    # Filter by minimum count
    df = df[df['count'] >= min_count].reset_index(drop=True)
    if df.empty:
        return df
    
    # Reconstruct qry_id
    df['qry_id'] = df.apply(lambda row: reconstruct_qry_id(
        row['qry_dataset'], row['qry_file_name'], row['qry_scan_id']
    ), axis=1)
    
    # Add match type indicator
    df['match_type'] = 'delta mass'
    df['conjugate_name'] = df['ref_1_name']  # For delta mass matches, use ref_1_name as conjugate name
    # conjugate_delta_mass is ref_1_mono_mass - water loss, rounded to 2 decimal places
    df['conjugate_delta_mass'] = df['ref_1_mono_mass'].apply(lambda x: round(x - 18.0106, 2))

    # Clean column order with all database info including mono masses
    column_order = [
        'qry_id', 'qry_mz', 'ion_polarity',
        'ref_1_id', 'ref_1_name', 'ref_1_mono_mass', 'ref_1_db',
        'ref_2_id', 'ref_2_name', 'ref_2_mono_mass', 'ref_2_db',
        'match_type', 'delta_mass', 'count', 'conjugate_name', 'conjugate_delta_mass'
    ]
    
    return df[column_order]


def reconstruct_qry_id(dataset, file_name, scan_id):
    """
    Reconstruct qry_id from components
    """
    return f"{dataset}:{file_name}:scan:{scan_id}"



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
        # "max_intensity": 150,
        # "annotate_precision": 4,
        # "annotation_rotation": 90,
        "cosine": "shifted",
        # "fragment_mz_tolerance": 0.05,
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
        "library_select": "metabolomicspanrepo_index_nightly",
        "usi1": usi
    }
    # URL encode the parameters
    query_string = urllib.parse.urlencode(params)
    
    # Create the full URL
    url = f"{base_url}?{query_string}"
    
    return url

if __name__ == "__main__":
    
    db_path = '/Users/shipei/Documents/projects/conjugated_metabolome/bin/web/main/conjugated_metabolome.db'
    target_inchikey_14 = 'OUYCCCASQSFEME'
    target_mono_mass = 181.073893
    min_count = 1
    
    results_df = filter_search_results(db_path, target_inchikey_14, target_mono_mass, min_count)
    print(results_df)
