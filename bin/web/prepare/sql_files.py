import sqlite3
import pandas as pd
import json
import os
from pathlib import Path


def optimize_qry_id(qry_id):
    """
    Convert qry_id to a more compact format by extracting components
    Format: '{dataset}:{file_name}:scan:{scan_id}'
    Returns: (dataset, file_name, scan_id)
    """
    if pd.isna(qry_id) or not isinstance(qry_id, str):
        raise ValueError(f"Invalid qry_id (not a string): {qry_id}")
    
    parts = qry_id.split(':')
    if len(parts) == 4 and parts[2] == 'scan':
        dataset = parts[0]
        file_name = parts[1]
        scan_id = int(parts[3])  # Convert scan_id to integer
        return dataset, file_name, scan_id
    else:
        raise ValueError(f"Invalid qry_id format (expected 'dataset:file:scan:id'): {qry_id}")


def decompress_qry_id(dataset, file_name, scan_id):
    """
    Reconstruct qry_id from components
    """
    if pd.isna(dataset) or pd.isna(file_name) or pd.isna(scan_id):
        return None
    return f"{dataset}:{file_name}:scan:{scan_id}"


def create_database_schema(db_path):
    """
    Create schema optimized for <1GB size and <1GB RAM with fast InChIKey queries
    """
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    # OPTIMAL PRAGMA settings based on your benchmarks
    cursor.execute("PRAGMA page_size=65536")
    cursor.execute("PRAGMA journal_mode=DELETE")
    cursor.execute("PRAGMA synchronous=OFF")
    cursor.execute("PRAGMA cache_size=10000")
    cursor.execute("PRAGMA temp_store=memory")
    cursor.execute("PRAGMA auto_vacuum=FULL")
    cursor.execute("PRAGMA mmap_size=134217728")     # 128MB for I/O
    
    # Create tables with maximum size optimization
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS chem_db (
            inchikey_14 TEXT(14) PRIMARY KEY,
            name TEXT NOT NULL,
            monoisotopic_mass REAL NOT NULL
        ) WITHOUT ROWID
    ''')
    
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS ms2db (
            db_id TEXT PRIMARY KEY,
            db INTEGER NOT NULL,
            inchikey_14 TEXT(14) NOT NULL,
            FOREIGN KEY (inchikey_14) REFERENCES chem_db (inchikey_14)
        ) WITHOUT ROWID
    ''')
    
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS search_results (
            id INTEGER PRIMARY KEY,
            qry_dataset INTEGER NOT NULL,
            qry_file_name INTEGER NOT NULL,
            qry_scan_id INTEGER NOT NULL,
            qry_mz REAL NOT NULL,
            ref_1_id TEXT,
            ref_2_id TEXT,
            delta_mass INTEGER NOT NULL,
            count INTEGER NOT NULL,
            ion_polarity INTEGER NOT NULL
        )
    ''')
    
    # Create lookup tables
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS db_lookup (
            id INTEGER PRIMARY KEY,
            name TEXT(10) UNIQUE
        ) WITHOUT ROWID
    ''')
    
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS dataset_lookup (
            id INTEGER PRIMARY KEY,
            name TEXT(15) UNIQUE
        ) WITHOUT ROWID
    ''')
    
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS filename_lookup (
            id INTEGER PRIMARY KEY,
            name TEXT(150) UNIQUE
        ) WITHOUT ROWID
    ''')
    
    # Insert lookup values
    cursor.execute("INSERT OR IGNORE INTO db_lookup (id, name) VALUES (1, 'gnps'), (2, 'massbank'), (3, 'mona'), (4, 'nist20')")
    
    # CRITICAL: Add the missing indexes that will fix your performance
    print("Creating essential indexes for InChIKey queries...")
    
    # These are ESSENTIAL for your query pattern:
    cursor.execute('CREATE INDEX IF NOT EXISTS idx_search_ref1 ON search_results(ref_1_id)')
    cursor.execute('CREATE INDEX IF NOT EXISTS idx_search_ref2 ON search_results(ref_2_id)')
    cursor.execute('CREATE INDEX IF NOT EXISTS idx_ms2db_inchikey ON ms2db(inchikey_14)')
    
    # Composite index for better join performance
    cursor.execute('CREATE INDEX IF NOT EXISTS idx_ms2db_composite ON ms2db(db_id, inchikey_14)')
    
    print("Essential indexes created")
    
    conn.commit()
    return conn

def process_ms2db(parquet_path, conn):
    """
    Process and insert ms2db data with optimized schema
    Creates both chem_db (deduplicated) and ms2db tables
    """
    print(f"Processing ms2db from {parquet_path}")
    df = pd.read_parquet(parquet_path)
    
    # Validate and clean data
    df = df.dropna(subset=['db_id', 'inchikey_14'])  # Remove rows without required fields
    df['monoisotopic_mass'] = df['monoisotopic_mass'].round(4)
    
    # Ensure db values are valid
    valid_dbs = ['gnps', 'massbank', 'mona', 'nist20']
    df = df[df['db'].isin(valid_dbs)]
    
    print(f"Processing {len(df)} records from ms2db")
    
    # Create chemical database (deduplicated by inchikey_14)
    print("Creating deduplicated chemical database...")
    
    # For each inchikey_14, keep the row with the shortest name
    def get_shortest_name(group):
        # Sort by name length, then alphabetically for consistency
        sorted_group = group.sort_values(['name'], key=lambda x: x.str.len())
        return sorted_group.iloc[0]
    
    chem_df = df.groupby('inchikey_14').apply(get_shortest_name).reset_index(drop=True)
    chem_df = chem_df[['inchikey_14', 'name', 'monoisotopic_mass']].copy()
    
    print(f"Deduplicated to {len(chem_df)} unique chemicals")
    
    # Insert chemical database in chunks
    chunk_size = 50000
    total_inserted = 0
    
    for i in range(0, len(chem_df), chunk_size):
        chunk = chem_df.iloc[i:i + chunk_size]
        chunk.to_sql('chem_db', conn, if_exists='append', index=False)
        total_inserted += len(chunk)
        print(f"Inserted {total_inserted}/{len(chem_df)} chemical records...")
    
    # Create ms2db table (only db_id, db, inchikey_14)
    print("Creating ms2db table...")
    ms2_df = df[['db_id', 'db', 'inchikey_14']].copy()
    
    total_inserted = 0
    for i in range(0, len(ms2_df), chunk_size):
        chunk = ms2_df.iloc[i:i + chunk_size]
        chunk.to_sql('ms2db', conn, if_exists='append', index=False)
        total_inserted += len(chunk)
        print(f"Inserted {total_inserted}/{len(ms2_df)} ms2db records...")
    
    print(f"Created chem_db with {len(chem_df)} unique chemicals")
    print(f"Created ms2db with {len(ms2_df)} spectrum records")
    
    return len(chem_df), len(ms2_df)

def process_search_results_with_maps(parquet_path, conn, ion_polarity, dataset_map, file_map):
    """
    Process search results using pre-created lookup maps
    """
    print(f"Processing search results from {parquet_path} (polarity: {ion_polarity})")
    df = pd.read_parquet(parquet_path)
    
    df = df[df['delta_mass'] > 1.5].reset_index(drop=True)
    
    # Process qry_id into components
    qry_components = df['qry_id'].apply(optimize_qry_id)
    df[['qry_dataset', 'qry_file_name', 'qry_scan_id']] = pd.DataFrame(qry_components.tolist(), index=df.index)
    
    # Apply the pre-created mappings
    df['qry_dataset'] = df['qry_dataset'].map(dataset_map)
    df['qry_file_name'] = df['qry_file_name'].map(file_map)
    
    # Check for unmapped values
    if df['qry_dataset'].isna().any():
        unmapped = df[df['qry_dataset'].isna()]['qry_dataset'].unique()
        print(f"WARNING: Found unmapped datasets: {unmapped}")
    
    if df['qry_file_name'].isna().any():
        unmapped = df[df['qry_file_name'].isna()]['qry_file_name'].unique()
        print(f"WARNING: Found unmapped filenames: {unmapped}")
    
    # Convert numeric data
    df['qry_mz'] = df['qry_mz'].round(4)
    df['delta_mass'] = (df['delta_mass'].round(2) * 100).astype(int)
    df['ion_polarity'] = 1 if ion_polarity == '+' else 0
    
    # Insert data
    columns_to_keep = [
        'qry_dataset', 'qry_file_name', 'qry_scan_id', 'qry_mz',
        'ref_1_id', 'ref_2_id', 'delta_mass', 'count', 'ion_polarity'
    ]
    df = df[columns_to_keep]
    
    chunk_size = 150000
    for i in range(0, len(df), chunk_size):
        chunk = df.iloc[i:i + chunk_size]
        chunk.to_sql('search_results', conn, if_exists='append', index=False)
    
    conn.commit()
    return len(df)


def create_lookup_tables(conn, pos_df, neg_df):
    """
    Create lookup tables with consistent IDs across both positive and negative data
    """
    cursor = conn.cursor()
    
    # Get all unique datasets and filenames from both dataframes
    all_datasets = pd.concat([pos_df['qry_dataset'], neg_df['qry_dataset']]).unique()
    all_filenames = pd.concat([pos_df['qry_file_name'], neg_df['qry_file_name']]).unique()
    
    # Create dataset lookup
    dataset_map = {}
    for i, dataset in enumerate(sorted(all_datasets), 1):  # Sort for consistency
        cursor.execute("INSERT OR IGNORE INTO dataset_lookup (id, name) VALUES (?, ?)", (i, dataset))
        dataset_map[dataset] = i
    
    # Create filename lookup
    file_map = {}
    for i, filename in enumerate(sorted(all_filenames), 1):  # Sort for consistency
        cursor.execute("INSERT OR IGNORE INTO filename_lookup (id, name) VALUES (?, ?)", (i, filename))
        file_map[filename] = i
    
    conn.commit()
    return dataset_map, file_map

def optimize_database(conn):
    """
    Optimize the database after all inserts
    """
    print("Final database compression...")
    cursor = conn.cursor()
    
    # Change to DELETE mode for final compression
    cursor.execute("PRAGMA journal_mode=DELETE")
    
    # Analyze tables
    cursor.execute("ANALYZE")
    
    # Full vacuum for maximum compression
    cursor.execute("VACUUM")
    
    # Set read-only optimizations for production use
    cursor.execute("PRAGMA optimize")
    
    conn.commit()
    print("Database optimization complete")


def convert_parquet_to_sqlite(input_dir, output_db_path):
    """
    Updated main function with consistent lookup table handling
    """
    input_path = Path(input_dir)
    
    # Define file paths
    ms2db_path = input_path / "ms2db.parquet"
    pos_path = input_path / "pos_web.parquet"
    neg_path = input_path / "neg_web.parquet"
    
    # Check if files exist
    missing_files = []
    for file_path, name in [(ms2db_path, "ms2db.parquet"), 
                           (pos_path, "pos_web.parquet"), 
                           (neg_path, "neg_web.parquet")]:
        if not file_path.exists():
            missing_files.append(name)
    
    if missing_files:
        print(f"Missing files: {missing_files}")
        return False
        
    # Remove existing database if it exists
    if os.path.exists(output_db_path):
        os.remove(output_db_path)
    
    # Create database and schema
    print(f"Creating SQLite database: {output_db_path}")
    conn = create_database_schema(output_db_path)
    
    try:
        # Process ms2db first
        chem_count, ms2db_count = process_ms2db(ms2db_path, conn)
        
        # Load both search result files to create consistent lookup tables
        print("Loading search result files to create lookup tables...")
        pos_df = pd.read_parquet(pos_path)
        neg_df = pd.read_parquet(neg_path)
        
        # Extract qry_id components from both files
        pos_components = pos_df['qry_id'].apply(optimize_qry_id)
        pos_df[['qry_dataset', 'qry_file_name', 'qry_scan_id']] = pd.DataFrame(pos_components.tolist(), index=pos_df.index)
        
        neg_components = neg_df['qry_id'].apply(optimize_qry_id)
        neg_df[['qry_dataset', 'qry_file_name', 'qry_scan_id']] = pd.DataFrame(neg_components.tolist(), index=neg_df.index)
        
        # Create consistent lookup tables
        print("Creating consistent lookup tables...")
        dataset_map, file_map = create_lookup_tables(conn, pos_df, neg_df)
        
        print(f"Created mappings for {len(dataset_map)} datasets and {len(file_map)} filenames")
        
        # Process search results with consistent mappings
        neg_count = process_search_results_with_maps(neg_path, conn, '-', dataset_map, file_map)
        pos_count = process_search_results_with_maps(pos_path, conn, '+', dataset_map, file_map)
        
        # Optimize database
        optimize_database(conn)
        
        # Print summary
        print(f"\n=== Conversion Summary ===")
        print(f"Unique chemicals (chem_db): {chem_count:,}")
        print(f"MS2 spectrum records (ms2db): {ms2db_count:,}")
        print(f"Positive search results: {pos_count:,}")
        print(f"Negative search results: {neg_count:,}")
        print(f"Total search results: {pos_count + neg_count:,}")
        
        # Check file size
        db_size = os.path.getsize(output_db_path)
        print(f"Database size: {db_size / (1024*1024):.2f} MB")
        
        return True
        
    except Exception as e:
        print(f"Error during conversion: {e}")
        import traceback
        traceback.print_exc()
        return False
        
    finally:
        conn.close()
    
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

def search_by_inchikey(db_path, target_inchikey="OUYCCCASQSFEME"):
    """
    Search for metabolites by InChIKey with optimized performance
    Returns results with all relevant database IDs and metadata
    """
    conn = create_production_connection(db_path)
    
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
    conn.close()
    
    if df.empty:
        print(f"No results found for InChIKey: {target_inchikey}")
        return df
    
    # Reconstruct qry_id
    df['qry_id'] = df.apply(lambda row: decompress_qry_id(
        row['qry_dataset'], row['qry_file_name'], row['qry_scan_id']
    ), axis=1)
    
    # Add match type indicator
    def get_match_type(row):
        ref1_match = row['ref_1_inchikey'] == target_inchikey
        ref2_match = row['ref_2_inchikey'] == target_inchikey
        
        if ref1_match and ref2_match:
            return 'both'
        elif ref1_match:
            return 'ref_1'
        elif ref2_match:
            return 'ref_2'
        else:
            return 'none'
    
    df['match_type'] = df.apply(get_match_type, axis=1)
    
    # Clean column order with all database info including mono masses
    column_order = [
        'qry_id', 'qry_mz', 'ion_polarity',
        'ref_1_id', 'ref_1_name', 'ref_1_mono_mass', 'ref_1_db', 'ref_1_inchikey',
        'ref_2_id', 'ref_2_name', 'ref_2_mono_mass', 'ref_2_db', 'ref_2_inchikey',
        'match_type', 'delta_mass', 'count'
    ]
    
    return df[column_order]


def search_by_delta_mass(db_path, target_delta_mass, mass_tolerance=0.01):
    """
    Search for metabolites by delta mass with specified tolerance
    Returns results with all relevant database IDs and metadata
    
    Args:
        db_path: Path to SQLite database
        target_delta_mass: Target delta mass (in Da)
        mass_tolerance: Mass tolerance (default 0.01 Da)
    """
    conn = create_production_connection(db_path)
    
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
    conn.close()
    
    if df.empty:
        print(f"No results found for delta mass: {target_delta_mass} ± {mass_tolerance} Da")
        return df
    
    # Reconstruct qry_id
    df['qry_id'] = df.apply(lambda row: decompress_qry_id(
        row['qry_dataset'], row['qry_file_name'], row['qry_scan_id']
    ), axis=1)

    # Add delta mass error, round to 2 decimal places
    df['delta_mass_error'] = abs(df['delta_mass'] - target_delta_mass).round(2)

    # Clean column order with all database info including mono masses
    column_order = [
        'qry_id', 'qry_mz', 'ion_polarity',
        'ref_1_id', 'ref_1_name', 'ref_1_mono_mass', 'ref_1_db', 'ref_1_inchikey',
        'ref_2_id', 'ref_2_name', 'ref_2_mono_mass', 'ref_2_db', 'ref_2_inchikey',
        'delta_mass', 'delta_mass_error', 'count'
    ]
    
    print(f"Found {len(df)} results for delta mass {target_delta_mass} ± {mass_tolerance} Da")
    return df[column_order]

if __name__ == "__main__":
    import time
    input_directory = "/Users/shipei/Documents/projects/conjugated_metabolome/bin/web/prepare"
    output_database = "/Users/shipei/Documents/projects/conjugated_metabolome/bin/web/main/conjugated_metabolome.db"
    
    # Convert parquet files to SQLite (comment out if already done)
    success = convert_parquet_to_sqlite(input_directory, output_database)

    ################
    test_inchikeys = ["XYNPYHXGMWJBLV", "BSYNRYMUTXBXSQ", "KZSNJWFQEVHDMF", 'PSIFNNKUMBGKDQ', 'JIVPVXMEBJLZRO']
    start_time = time.time()    
    for inchikey in test_inchikeys:
        this_start = time.time()
        inchikey_results = search_by_inchikey(output_database, inchikey)
        this_end = time.time()
        print(f"InChIKey {inchikey}: {len(inchikey_results)} results in {(this_end - this_start)*1000:.1f} ms")    
        
        # # save
        # out_path = f"/Users/shipei/Documents/projects/conjugated_metabolome/bin/web/prepare/inchikey_{inchikey}_results.csv"
        # inchikey_results.to_csv(out_path, index=False)
    end_time = time.time()
    print(f"Total test time: {(end_time - start_time)*1000:.1f} ms")

    ################
    # test_delta_masses = [176.0321, 162.0528, 80.9642, 79.9568, 97.9769]
    # start_time = time.time()    
    # for delta_mass in test_delta_masses:
    #     this_start = time.time()
    #     delta_mass_results = search_by_delta_mass(output_database, delta_mass, mass_tolerance=0.01)
    #     this_end = time.time()
    #     print(f"Delta mass {delta_mass}: {len(delta_mass_results)} results in {(this_end - this_start)*1000:.1f} ms")
        
    #     # save
    #     out_path = f"/Users/shipei/Documents/projects/conjugated_metabolome/bin/web/prepare/delta_mass_{delta_mass:.4f}_results.csv"
    #     delta_mass_results.to_csv(out_path, index=False)
        
    # end_time = time.time()
    # print(f"Total test time: {(end_time - start_time)*1000:.1f} ms")

