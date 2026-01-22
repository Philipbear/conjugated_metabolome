import streamlit as st
import pandas as pd
import sqlite3
import os
import urllib.parse


def render_delta_mass_browser():
    """Render the delta mass browser interface"""
    
    st.title("Browse by Delta Mass")
    st.markdown("Explore metabolite conjugates by their conjugation mass differences")
    
    st.markdown("---")
    
    # Common conjugate presets
    st.subheader("Common Conjugation Patterns")
    
    presets = {
        "Sulfate": 79.96,
        "Glucuronide": 176.03,
        "Hexose": 162.05,
        # "Deoxyhexose": 146.06,
        "Pentose": 132.04,
        "Phosphate": 79.97,
        # "Deoxypentose": 116.05,
        "Glycine": 57.02,
        "Taurine": 107.00,
        # "Glutamine": 128.06,
        "Proline": 97.05,
        # "Phenylalanine": 147.07,
        # "Aspartate": 115.03,
        # "Tryptophan": 186.08,
        "C16:0": 256.24,
        "C16:1": 254.24,
        "C18:0": 282.26,
        "C18:1": 280.26
    }
    
    # Initialize session state for delta mass if not exists
    if 'selected_delta_mass' not in st.session_state:
        st.session_state.selected_delta_mass = 162.05
    
    # Create preset buttons
    cols = st.columns(6)
    
    for idx, (name, mass) in enumerate(presets.items()):
        col = cols[idx % 6]
        with col:
            if st.button(f"{name}\n({mass:.2f} Da)", use_container_width=True, key=f"preset_{name}"):
                st.session_state.selected_delta_mass = mass
                st.rerun()
    
    st.markdown("---")
    
    # Input section
    st.subheader("Search Parameters")
    
    col1, col2, col3 = st.columns([2, 1, 1])
    
    with col1:
        delta_mass = st.number_input(
            "Delta Mass (Da):",
            min_value=0.0,
            max_value=1000.0,
            value=float(st.session_state.selected_delta_mass),
            step=0.01,
            format="%.2f"
        )
        # Update session state when user manually changes the input
        st.session_state.selected_delta_mass = delta_mass
    
    with col2:
        tolerance = st.number_input(
            "Tolerance (Da):",
            min_value=0.0,
            max_value=0.10,
            value=0.01,
            step=0.005,
            format="%.3f"
        )
    
    with col3:
        min_count = st.number_input(
            "Min dataset count:",
            min_value=1,
            max_value=100,
            value=1,
            step=1
        )
    
    # Search button
    btn_col, _ = st.columns([1, 3])
    with btn_col:
        if st.button("🔍 Search", type="primary", use_container_width=True):
            with st.spinner("Searching database..."):
                results = search_by_delta_mass_range(delta_mass, tolerance, min_count)
                st.session_state['delta_mass_results'] = results
    
    # Display results
    if 'delta_mass_results' in st.session_state and st.session_state['delta_mass_results'] is not None:
        display_delta_mass_results(st.session_state['delta_mass_results'], delta_mass, tolerance)


def search_by_delta_mass_range(target_mass, tolerance, min_count):
    """Search for conjugates within a delta mass range"""
    
    db_path = get_db_path()
    if not db_path:
        return None
    
    conn = sqlite3.connect(db_path)
    
    # Convert to stored format (multiplied by 100)
    target_stored = round(target_mass * 100)
    tolerance_stored = tolerance * 100  # Don't round tolerance to preserve precision
    min_delta = target_stored - tolerance_stored
    max_delta = target_stored + tolerance_stored
    
    query = '''
        SELECT DISTINCT
            dl.name as qry_dataset, 
            fl.name as qry_file_name, 
            sr.qry_scan_id, 
            sr.qry_mz,
            sr.ref_1_id, 
            sr.delta_mass / 100.0 as delta_mass,
            sr.count, 
            CASE sr.ion_polarity WHEN 1 THEN '+' ELSE '-' END as ion_polarity,
            ms1.inchikey_14 as ref_1_inchikey, 
            c1.name as ref_1_name,
            c1.monoisotopic_mass as ref_1_mono_mass,
            ms1.db as ref_1_db
        FROM search_results sr
        LEFT JOIN dataset_lookup dl ON sr.qry_dataset = dl.id
        LEFT JOIN filename_lookup fl ON sr.qry_file_name = fl.id
        LEFT JOIN ms2db ms1 ON sr.ref_1_id = ms1.db_id
        LEFT JOIN chem_db c1 ON ms1.inchikey_14 = c1.inchikey_14
        WHERE sr.delta_mass BETWEEN ? AND ?
          AND sr.count >= ?
        ORDER BY sr.count DESC
    '''
    
    df = pd.read_sql_query(query, conn, params=[min_delta, max_delta, min_count])
    conn.close()
    
    if df.empty:
        return df
    
    # Reconstruct qry_id
    df['qry_usi'] = df.apply(lambda row: f"mzspec:{row['qry_dataset']}:{row['qry_file_name']}:scan:{row['qry_scan_id']}", axis=1)
    
    # Add reference USI column
    df['ref_1_usi'] = df.apply(_get_ref_usi, axis=1)
    
    # Add mirror plot and MASST URLs
    df['mirror_plot'] = df.apply(
        lambda row: gen_mirror_plot_url(row['qry_usi'], row['ref_1_usi']), axis=1
    )
    df['masst'] = df['qry_usi'].apply(gen_fasst_url)
    
    return df


def _get_ref_usi(row):
    """Generate USI for reference spectrum"""
    if pd.isna(row['ref_1_id']):
        return None
    ref_id = str(row['ref_1_id'])
    if row['ref_1_db'] == 'gnps':
        return 'mzspec:GNPS:GNPS-LIBRARY:accession:CCMSLIB000' + ref_id
    elif row['ref_1_db'] == 'mona':
        return 'mzspec:MASSBANK::accession:' + ref_id
    else:
        return None


def gen_mirror_plot_url(usi1, usi2):
    """Generate a URL for the metabolomics-usi.gnps2.org mirror plot"""
    
    if usi1 is None or usi2 is None:
        return None

    base_url = "https://metabolomics-usi.gnps2.org/dashinterface/"
    
    params = {
        "usi1": usi1,
        "usi2": usi2,
        "cosine": "shifted",
    }
    
    query_string = urllib.parse.urlencode(params)
    url = f"{base_url}?{query_string}"
    
    return url


def gen_fasst_url(usi):
    """Generate a URL for fasst MASST"""
    
    base_url = "https://fasst.gnps2.org/fastsearch/"
    
    params = {
        "library_select": "metabolomicspanrepo_index_nightly",
        "usi1": usi
    }
    
    query_string = urllib.parse.urlencode(params)
    url = f"{base_url}?{query_string}"
    
    return url


def display_delta_mass_results(df, target_mass, tolerance):
    """Display delta mass search results"""
    
    if df.empty:
        st.warning(f"❌ No conjugates found with delta mass {target_mass} ± {tolerance} Da")
        return
    
    # Summary statistics
    st.success(f"✅ Found {len(df)} matches")
    
    # Center the metrics
    st.markdown("""
        <style>
        [data-testid="stMetricValue"], [data-testid="stMetricLabel"] {
            text-align: center;
            justify-content: center;
        }
        </style>
    """, unsafe_allow_html=True)
    
    # Calculate unique compounds by polarity
    total_unique = df['ref_1_inchikey'].nunique()
    pos_unique = df[df['ion_polarity'] == '+']['ref_1_inchikey'].nunique()
    neg_unique = df[df['ion_polarity'] == '-']['ref_1_inchikey'].nunique()
    
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("Unique Compounds (Total)", total_unique)
    with col2:
        st.metric("Unique Compounds (positive polarity data)", pos_unique)
    with col3:
        st.metric("Unique Compounds (negative polarity data)", neg_unique)
    with col4:
        st.metric("Total Observations", df['count'].sum())
    
    st.markdown("---")
    
    # Results table with filters
    st.subheader("Search Results")
    
    # Prepare data for display
    df_display = df[['ion_polarity', 'count', 'qry_mz', 'ref_1_name', 
                     'ref_1_inchikey', 'ref_1_mono_mass', 'mirror_plot', 'masst']].copy()
    
    # Filter controls
    _, filter_col1, filter_col2, _ = st.columns([1, 2, 2, 1])
    
    with filter_col1:
        polarity_filter = st.selectbox('Ion polarity:', ['All', '+', '-'])
    
    with filter_col2:
        sort_by = st.selectbox('Sort by:', ['Dataset count', 'Precursor m/z', 'Compound Name'])
    
    # Apply filters
    filtered_df = df_display.copy()
    
    if polarity_filter != 'All':
        filtered_df = filtered_df[filtered_df['ion_polarity'] == polarity_filter]
    
    # Apply sorting
    if sort_by == 'Dataset count':
        filtered_df = filtered_df.sort_values('count', ascending=False)
    elif sort_by == 'Precursor m/z':
        filtered_df = filtered_df.sort_values('qry_mz', ascending=True)
    elif sort_by == 'Compound Name':
        filtered_df = filtered_df.sort_values('ref_1_name', ascending=True)
    
    # Display info
    _, info_col, _ = st.columns([1, 8, 1])
    with info_col:
        st.info(f"Showing {len(filtered_df)} of {len(df_display)} results")
    
    # Display table
    if not filtered_df.empty:
        st.dataframe(
            filtered_df,
            column_config={
                "ion_polarity": st.column_config.TextColumn(
                    "Ion polarity", 
                    width="small",
                    help="Ionization polarity of the query MS/MS"
                ),
                "count": st.column_config.NumberColumn(
                    "Dataset count", 
                    width="small",
                    help="Dataset count in public domains", 
                    format="%d"
                ),
                "qry_mz": st.column_config.NumberColumn(
                    "Precursor m/z", 
                    width="small",
                    help="m/z of the precursor ion", 
                    format="%.4f"
                ),
                "ref_1_name": st.column_config.TextColumn(
                    "Compound Name", 
                    width="medium",
                    help="Name of the matched compound"
                ),
                "ref_1_inchikey": st.column_config.TextColumn(
                    "Compound InChIKey (planar)", 
                    width="small",
                    help="InChIKey-14 of the matched compound"
                ),
                "ref_1_mono_mass": st.column_config.NumberColumn(
                    "Compound exact mass", 
                    width="small",
                    help="Monoisotopic mass of the compound", 
                    format="%.4f"
                ),
                "mirror_plot": st.column_config.LinkColumn(
                    "Mirror plot", 
                    width="small",
                    help="Click to view mirror plot between query MS/MS and reference",
                    display_text="View", 
                    required=False
                ),
                "masst": st.column_config.LinkColumn(
                    "MASST", 
                    width="small",
                    help="Link to the fast MASST search results",
                    display_text="🔍 MASST", 
                    required=False
                )
            },
            hide_index=True,
            use_container_width=True
        )


def get_db_path():
    """Get database path"""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    db_path = os.path.join(script_dir, "conjugated_metabolome.db")
    if os.path.exists(db_path):
        return db_path
    else:
        st.error(f"❌ Database SQLite file not found: {db_path}")
        return None
