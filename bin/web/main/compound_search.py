import streamlit as st
import os
from sql_utils import filter_search_results, prepare_delta_mass_plot
from chem_utils import smiles_to_formula_inchikey, calc_monoisotopic_mass, inchikey_to_common_name, get_structure_image_gnps2, get_compound_description_pubchem
from pubchem_utils import pubchem_autocomplete, name_to_cid, cid_to_canonical_smiles


def render_compound_search():
    """Render the compound search interface"""
    
    # Initialize session state for compound search
    initialize_search_session_state()
    
    db_path = get_db_path()
    if not db_path:
        return
    
    st.title("Compound Search")
    st.markdown("Search for conjugated metabolites by compound name or SMILES string.")
    
    # st.header("Search")
    
    # Handle search inputs and get effective SMILES
    effective_smiles, compound_name = handle_search_interface()
    
    # Process search if we have valid input
    if effective_smiles:
        process_search(effective_smiles, compound_name, db_path)
    
    # Display existing results (independent of new searches)
    display_stored_results()


def initialize_search_session_state():
    """Initialize all session state variables for compound search"""
    session_vars = {
        'search_history': [],
        'name_suggestions': [],
        'show_suggestions': False,
        'selected_smiles': "",
        'selected_compound_name': "",
        'current_search_results': None,
        'current_compound_info': None
    }
    
    for var, default_value in session_vars.items():
        if var not in st.session_state:
            st.session_state[var] = default_value


def handle_search_interface():
    """Handle all search interface logic and return effective SMILES and compound name"""
    # Add tabs for different search methods
    search_tab1, search_tab2 = st.tabs(["📚 Search by Name", "📝 Search by SMILES"])
    
    effective_smiles = ""
    compound_name = ""
    
    with search_tab1:
        effective_smiles, compound_name = handle_name_search()
    
    with search_tab2:
        smiles_result = handle_smiles_search()
        if smiles_result:
            effective_smiles = smiles_result
            compound_name = "Chemical compound"  # Default name for SMILES search
    
    return effective_smiles, compound_name


def handle_name_search():
    """Handle name search interface"""
    st.markdown("Search for compounds by name using PubChem:")
    
    # Demo compound buttons
    st.markdown("**Try demo compounds:**")
    demo_col1, demo_col2, _ = st.columns([2, 2, 8])
    
    with demo_col1:
        if st.button("Ferulic Acid", use_container_width=True, key="demo_ferulic_acid"):
            st.session_state.selected_smiles = "COC1=C(C=CC(=C1)/C=C/C(=O)O)O"
            st.session_state.selected_compound_name = "Ferulic acid"
            st.rerun()
    
    with demo_col2:
        if st.button("Phenylacetic Acid", use_container_width=True, key="demo_phenylacetic_acid"):
            st.session_state.selected_smiles = "C1=CC=C(C=C1)CC(=O)O"
            st.session_state.selected_compound_name = "Phenylacetic acid"
            st.rerun()
    
    st.markdown("---")
    
    # Name input
    compound_name = st.text_input(
        "Enter compound name:",
        placeholder=f"e.g., tyrosine, phenylacetic acid, ferulic acid",
        value=st.session_state.get('compound_name_input', ''),
        key="compound_name_input"
    )
    compound_name = compound_name.strip()
    
    # Buttons
    col1, _ = st.columns([3, 10])
    
    with col1:
        search_pressed = st.button(
            "**Retrieve Candidate Compounds**",
            type="primary",
            use_container_width=True,
            icon=":material/search:",
            key="search_name_button"
        )
    
    # Handle search button
    if search_pressed and compound_name:
        with st.spinner("Searching PubChem..."):
            suggestions = pubchem_autocomplete(compound_name)
        
        if suggestions:
            st.session_state.name_suggestions = suggestions
            st.session_state.show_suggestions = True
        else:
            st.error(f"No compounds found for '{compound_name}'.")
            st.session_state.show_suggestions = False
    
    # Handle compound selection
    if st.session_state.get('show_suggestions') and st.session_state.get('name_suggestions'):
        return handle_compound_selection()
    
    # Handle selected SMILES from name search
    if st.session_state.get('selected_smiles'):
        smiles = st.session_state.selected_smiles
        name = st.session_state.get('selected_compound_name', 'Selected compound')
        # Clear after use
        st.session_state.selected_smiles = ""
        st.session_state.selected_compound_name = ""
        return smiles, name
    
    return "", ""


def handle_compound_selection():
    """Handle compound selection dropdown and return SMILES if selected"""
    st.subheader("Select a compound:")
    
    compound_options = ['Select a compound...'] + st.session_state.name_suggestions
    selected_compound = st.selectbox(
        "Available compounds:",
        options=compound_options,
        key="compound_selector"
    )
    
    if selected_compound and selected_compound != 'Select a compound...':
        with st.spinner(f"Getting SMILES for {selected_compound}..."):
            cid = name_to_cid(selected_compound)
            if cid:
                smiles = cid_to_canonical_smiles(cid)
                if smiles:
                    st.success(f"✅ **{selected_compound}**")
                    st.code(smiles, language=None)
                    
                    if st.button(
                        f"Search for conjugates of {selected_compound}",
                        type="primary",
                        use_container_width=True,
                        icon=":material/search:",
                        key="proceed_with_compound"
                    ):
                        # Clear selection state
                        st.session_state.show_suggestions = False
                        if 'compound_selector' in st.session_state:
                            del st.session_state.compound_selector
                        return smiles, selected_compound
                else:
                    st.error(f"Could not get SMILES for {selected_compound}")
            else:
                st.error(f"Could not find compound ID for {selected_compound}")
    
    return "", ""


def handle_smiles_search():
    """Handle SMILES search interface"""
    # Determine default value
    default_value = ""
    if st.session_state.selected_smiles:
        default_value = st.session_state.selected_smiles
    
    # SMILES input
    smiles_input = st.text_input(
        "Enter SMILES:",
        placeholder="Enter a valid SMILES string",
        value=default_value,
        key="smiles_input_field"
    )
    smiles_input = smiles_input.strip()
    
    # Buttons
    col1, _ = st.columns([3, 10])
    
    with col1:
        search_pressed = st.button(
            "**Search for conjugates**",
            type="primary",
            use_container_width=True,
            icon=":material/search:",
            key="search_smiles_button"
        )
    
    # Handle search button
    if search_pressed and smiles_input:
            return smiles_input
    
    return ""


def process_search(effective_smiles, compound_name, db_path):
    """Process a new search and store results"""
    # Convert SMILES to formula and InChIKey
    formula, inchikey = smiles_to_formula_inchikey(effective_smiles)
    
    if inchikey is None:
        st.error("❌ Invalid SMILES string. Please double-check your input.")
        return
    
    # Get compound information
    common_names = inchikey_to_common_name(inchikey)
    if not compound_name or compound_name == "Chemical compound":
        compound_name = common_names[0] if common_names else "Chemical structure"
    
    mono_mass = calc_monoisotopic_mass(formula)
    
    # Store compound info
    st.session_state.current_compound_info = {
        'smiles': effective_smiles,
        'compound_name': compound_name,
        'common_names': common_names,
        'formula': formula,
        'inchikey': inchikey,
        'mono_mass': mono_mass
    }
    
    # Add to search history
    add_to_search_history(compound_name, effective_smiles)
    
    # Perform search
    inchikey_14 = inchikey[:14]
    df_results = filter_search_results(db_path, inchikey_14, mono_mass)
    
    # Store results
    st.session_state.current_search_results = df_results


def display_stored_results():
    """Display stored search results and compound info"""
    if st.session_state.current_compound_info is None:
        return
    
    # Display compound information
    display_compound_info()
    
    # Display search results
    if st.session_state.current_search_results is not None:
        display_search_results()


def display_compound_info():
    """Display compound information section"""
    compound_info = st.session_state.current_compound_info
    
    st.subheader("Compound Information")
    
    _, structure_col, _, info_col, description_col = st.columns([1, 3, 1, 7, 7])
    
    with structure_col:
        image_url = get_structure_image_gnps2(compound_info['smiles'])
        st.image(image_url, caption=compound_info['compound_name'], use_container_width=True)
    
    with info_col:
        if compound_info['common_names']:
            names_str = '| '.join(compound_info['common_names'][:3])
            st.markdown(f"**Common names:** {names_str}")
        st.markdown(f"**SMILES:**\n```\n{compound_info['smiles']}\n```")
        st.markdown(f"**Formula:** {compound_info['formula']}")
        st.markdown(f"**InChIKey:** {compound_info['inchikey']}")
        st.markdown(f"**Monoisotopic mass:** {compound_info['mono_mass']:.4f}")
    
    with description_col:
        st.markdown("**Compound description (from PubChem):**")
        description = get_compound_description_pubchem(compound_info['smiles'])
        if description:
            st.write(description)


def display_search_results():
    """Display search results with filtering"""
    df_results = st.session_state.current_search_results
    
    if df_results.empty:
        st.warning("❌ No matches found for the searched compound.")
        return
    
    # Results summary
    total_matches = len(df_results)
    st.info(f"✅ Found {total_matches} matches for the target compound.")
    
    # Chart for multiple results
    if len(df_results) > 1:
        display_delta_mass_chart(df_results)
    
    # Results table with filtering
    display_results_table(df_results)


def display_delta_mass_chart(df_results):
    """Display delta mass distribution chart"""
    st.subheader("Distribution of Conjugate Delta Masses")
    
    delta_mass_counts = prepare_delta_mass_plot(df_results)
    
    _, chart_col, _ = st.columns([1, 12, 1])
    with chart_col:
        st.scatter_chart(
            data=delta_mass_counts,
            x='conjugate_delta_mass',
            y='count',
            x_label='Conjugate delta mass (Da)',
            y_label='Dataset frequency',
            size='count',
            color='ion_polarity',
            height=450,
            use_container_width=True
        )


def display_results_table(df_results):
    """Display results table with filtering controls"""
    st.subheader("Result Table")
    
    # Prepare data for display
    df_display = df_results[['ion_polarity', 'annotation_type', 'count', 'qry_mz',
                            'conjugate_delta_mass', 'conjugate_name', 'mirror_plot_ref_1',
                            'mirror_plot_ref_2', 'match_type', 'masst']].copy()
    
    # Filter controls
    _, filter_col1, filter_col2, filter_col3, filter_col4, _ = st.columns([1, 2, 2, 2, 2, 1])
    
    with filter_col1:
        polarity_filter = st.selectbox('Ion polarity:', ['+/-', '+', '-'])
    
    with filter_col2:
        annotation_filter = st.selectbox('Annotation type:', ['All matches', 'spec_spec', 'spec_delta'], index=0)
    
    with filter_col3:
        name_filter = st.selectbox('Conjugate name:', ['All matches', 'With name (annotated)', 'Without name (unannotated)'], index=0)
    
    with filter_col4:
        match_filter = st.selectbox('Match type:', 
                                   ['All matches', 'spec (ref 1)', 'spec (ref 2)', 'spec (ref 1) or spec (ref 2)', 'delta mass'], index=3)
    
    # Apply filters
    filtered_results = apply_filters(df_display, polarity_filter, annotation_filter, name_filter, match_filter)
    
    if filtered_results.empty:
        # display all results if no matches after filtering
        filtered_results = apply_filters(df_display, '+/-', 'All matches', 'All matches', 'All matches')
        # show warning
        st.warning("No results match the current filter criteria. Displaying all results instead.")
    
    # Always sort by dataset count in descending order
    filtered_results = filtered_results.sort_values('count', ascending=False)
    
    # Display filter results summary
    _, info_col, _ = st.columns([1, 8, 1])
    with info_col:
        st.info(f"Showing {len(filtered_results)} of {len(df_display)} results")
    
    # Display table
    if not filtered_results.empty:
        st.dataframe(
            filtered_results,
            column_config=get_column_config(),
            hide_index=True,
            use_container_width=True
        )
    else:
        st.warning("No results match the current filter criteria.")


def apply_filters(df, polarity_filter, annotation_filter, name_filter, match_filter):
    """Apply all filters to the dataframe"""
    filtered_df = df.copy()
    
    if polarity_filter != '+/-':
        filtered_df = filtered_df[filtered_df['ion_polarity'] == polarity_filter]
    
    if annotation_filter != 'All matches':
        filtered_df = filtered_df[filtered_df['annotation_type'] == annotation_filter]
    
    if match_filter != 'All matches':
        if match_filter == 'spec (ref 1) or spec (ref 2)':
            filtered_df = filtered_df[filtered_df['match_type'] != 'delta mass']
        else:
            filtered_df = filtered_df[filtered_df['match_type'] == match_filter]
    
    if name_filter == 'With name (annotated)':
        filtered_df = filtered_df[filtered_df['conjugate_name'].notna() & 
                                 (filtered_df['conjugate_name'] != '')]
    elif name_filter == 'Without name (unannotated)':
        filtered_df = filtered_df[filtered_df['conjugate_name'].isna() | 
                                 (filtered_df['conjugate_name'] == '')]
    
    return filtered_df


def get_column_config():
    """Get column configuration for the results table"""
    return {
        "ion_polarity": st.column_config.TextColumn(
            "Ion polarity", width="small",
            help="Ionization polarity of the query MS/MS"
        ),
        "annotation_type": st.column_config.TextColumn(
            "Annotation type", width="small",
            help="How query MS/MS are annotated (spec_spec: both components have spectral matches; spec_delta: one component has spectral match)"
        ),
        "count": st.column_config.NumberColumn(
            "Dataset count", width="small",
            help="Dataset count in public domains", format="%d"
        ),
        "qry_mz": st.column_config.NumberColumn(
            "Precursor m/z", width="small",
            help="m/z of the precursor ion", format="%.4f"
        ),
        "conjugate_delta_mass": st.column_config.NumberColumn(
            "Conjugate delta mass", width="small",
            help="Mass of the conjugate component", format="%.2f"
        ),
        "conjugate_name": st.column_config.TextColumn(
            "Conjugate name", width="medium",
            help="Name of the conjugate component"
        ),
        "mirror_plot_ref_1": st.column_config.LinkColumn(
            "Mirror plot (Ref 1)", width="small",
            help="Click to view mirror plot between query MS/MS and reference 1",
            display_text="View", required=False
        ),
        "mirror_plot_ref_2": st.column_config.LinkColumn(
            "Mirror plot (Ref 2)", width="small",
            help="Click to view mirror plot between query MS/MS and reference 2",
            display_text="View", required=False
        ),
        "match_type": st.column_config.TextColumn(
            "Match type", width="small",
            help="How the target compound is found (spectral match or delta mass search)"
        ),
        "masst": st.column_config.LinkColumn(
            "MASST", width="small",
            help="Link to the fast MASST search results",
            display_text="🔍 MASST", required=False
        )
    }


# Helper functions
def get_db_path():
    """Get database path"""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    db_path = os.path.join(script_dir, "conjugated_metabolome.db")
    if os.path.exists(db_path):
        return db_path
    else:
        st.error(f"❌ Database SQLite file not found: {db_path}")
        return None


def add_to_search_history(compound_name, smiles):
    """Add compound name and SMILES to search history"""
    if smiles and compound_name:
        search_entry = (compound_name, smiles)
        existing_smiles = [entry[1] for entry in st.session_state.search_history]
        if smiles not in existing_smiles:
            st.session_state.search_history.insert(0, search_entry)
            if len(st.session_state.search_history) > 10:
                st.session_state.search_history = st.session_state.search_history[:10]


def render_search_history():
    """Render search history sidebar component"""
    st.header("🕒 Search History")
    display_search_history()
    
    if st.session_state.search_history:
        if st.button("🗑️ Clear All History", type="secondary"):
            st.session_state.search_history = []
            st.rerun()


def display_search_history():
    """Display search history with clickable compound names"""
    if st.session_state.search_history:
        for i, (compound_name, hist_smiles) in enumerate(st.session_state.search_history):
            col1, col2 = st.columns([5, 1])
            
            with col1:
                if st.button(
                    f"🔍 {compound_name}",
                    key=f"history_{i}",
                    help=f"Click to search: {compound_name}"
                ):
                    st.session_state.selected_smiles = hist_smiles
                    st.session_state.selected_compound_name = compound_name
                    st.rerun()
            
            with col2:
                if st.button("❌", key=f"delete_{i}", help="Remove from history", type="secondary"):
                    st.session_state.search_history.pop(i)
                    st.rerun()
    else:
        st.info("No search history yet. Start searching to see your history here!")
