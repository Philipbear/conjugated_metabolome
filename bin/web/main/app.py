import streamlit as st
import pandas as pd
import os
from sql_utils import filter_search_results, prepare_delta_mass_plot, get_git_short_rev
from chem_utils import smiles_to_formula_inchikey, calc_monoisotopic_mass, inchikey_to_common_name, get_structure_image_gnps2, get_compound_description_pubchem
from pubchem_utils import pubchem_autocomplete, name_to_cid, cid_to_canonical_smiles


DEMO_SMILES = 'C1=CC(=CC=C1C[C@@H](C(=O)O)N)O'


def main():
    # Set the page configuration
    app_version = "2025-10-01"
    try:
        git_hash = get_git_short_rev()
    except:
        git_hash = "unknown"
    repo_link = "https://github.com/Philipbear/conjugated_metabolome"

    st.set_page_config(page_title="Conjugated Metabolome Explorer (under development)", layout="wide",
                       menu_items={"About": (f"**App Version**: {app_version} | "
                                             f"[**Git Hash**: {git_hash}]({repo_link}/commit/{git_hash})")}
                       )
    
    initialize_search_history()

    # Left panel with app info
    with st.sidebar:
        st.title("Conjugated Metabolome Explorer (under development)")
        st.image("https://ccms-ucsd.github.io/GNPSDocumentation/img/logo/GNPS_logo_original_transparent.png", width=150)
        
        st.markdown("""
        ### 📖 About
        This app allows you to explore potential metabolite conjugations using compound names or SMILES strings.
        
        This webpage does not include all conjugation results. For more comprehensive results, please refer to [our paper](https://github.com/Philipbear/conjugated_metabolome) and [Zenodo repository](https://github.com/Philipbear/conjugated_metabolome).
            
        ### 📝 Citation
        Please use it responsibly and cite [our work](https://github.com/Philipbear/conjugated_metabolome) if you find it useful:
        - S. Xing, V. Charron-Lamoureux, A. Patan, ..... [Navigating the underexplored conjugated metabolome](https://github.com/Philipbear/conjugated_metabolome). 2025
        
        ### 📧 Contact
        For questions or feedback, please contact Shipei Xing at
        [philipxsp@hotmail.com](mailto:philipxsp@hotmail.com)
        """)
        
        st.header("🕒 Search History")
        display_search_history()
        
        # Clear all history button
        if st.session_state.search_history:
            if st.button("🗑️ Clear All History", type="secondary"):
                st.session_state.search_history = []
                st.rerun()

    # Create a layout
    _, main_col, _ = st.columns([1, 8, 1])
    # Main panel with search functionality
    with main_col:       
        # get database path
        db_path = get_db_path()
        
        initialize_demo_smiles()
        initialize_name_search()
        
        st.title("Conjugated Metabolome Explorer")
        # Create a search section
        st.header("Search")
        
        # Add tabs for different search methods (Name search prioritized)
        search_tab1, search_tab2 = st.tabs(["🔍 Search by Name", "🧪 Search by SMILES"])

        with search_tab1:
            st.markdown("Search for compounds by name using PubChem:")
            
            # Create columns for name search
            name_col1, name_col2 = st.columns([3, 1])
            
            with name_col1:
                # Check if we should load demo name
                if 'demo_name' in st.session_state and st.session_state.demo_name:
                    default_name_value = st.session_state.demo_name
                    # Clear the demo flag after using it
                    st.session_state.demo_name = ""
                else:
                    default_name_value = st.session_state.get('compound_name_input', '')
                
                compound_name = st.text_input(
                    "Enter compound name:",
                    placeholder="e.g., aspirin, glucose, caffeine",
                    value=default_name_value,
                    key="compound_name_input"
                )
                compound_name = compound_name.strip()
            
            with name_col2:
                # Use the same min_count as SMILES search
                min_count_name = st.number_input("Min frequency:", min_value=1, max_value=50, value=3, step=1,
                                                help="Minimum frequency of a conjugation observed in public LC-MS/MS datasets to be included.",
                                                format="%d", key="min_count_name")
            
            # Create a row for the buttons
            name_button_col1, _, name_button_col2, _ = st.columns([2, 1, 2, 7])
            
            # Search by name button
            with name_button_col1:
                name_search_button = st.button(
                    "**Search Name**",
                    type="primary",
                    use_container_width=True,
                    icon=":material/search:",
                    help="Search PubChem for compound names",
                    key="search_name_button"
                )
            
            # Add Demo button for name search
            with name_button_col2:
                name_demo_button = st.button(
                    "**Load Demo**",
                    type="secondary",
                    use_container_width=True,
                    icon=":material/login:",
                    help="Load example: Tyrosine",
                    key="name_demo_button"
                )
            
            # Process the name demo button
            if name_demo_button:
                # Set demo name in session state and trigger rerun
                st.session_state.demo_name = "Tyrosine"
                st.rerun()
            
            # Handle name search
            if name_search_button and compound_name:
                with st.spinner("Searching PubChem for compound names..."):
                    suggestions = pubchem_autocomplete(compound_name)
                
                if suggestions:
                    st.session_state.name_suggestions = suggestions
                    st.session_state.show_suggestions = True
                else:
                    st.error(f"No compounds found for '{compound_name}'. Please try a different name.")
                    st.session_state.show_suggestions = False
            
            # Display suggestions if available
            if st.session_state.get('show_suggestions', False) and st.session_state.get('name_suggestions', []):
                st.subheader("Select a compound:")
                
                # Create a dropdown for compound selection
                compound_options = ['Select a compound...'] + st.session_state.name_suggestions
                selected_compound = st.selectbox(
                    "Available compounds:",
                    options=compound_options,
                    key="compound_selector",
                    help="Choose a compound from the search results"
                )
                
                # Handle compound selection
                if selected_compound and selected_compound != 'Select a compound...':
                    with st.spinner(f"Getting SMILES for {selected_compound}..."):
                        cid = name_to_cid(selected_compound)
                        if cid:
                            smiles = cid_to_canonical_smiles(cid)
                            if smiles:
                                # Display the SMILES for the selected compound
                                st.success(f"✅ **{selected_compound}**")
                                st.code(smiles, language=None)
                                
                                # Add a button to proceed with this compound
                                if st.button(
                                    f"🔍 Search for conjugates of {selected_compound}",
                                    type="primary",
                                    use_container_width=True,
                                    key="proceed_with_compound"
                                ):
                                    st.session_state.selected_smiles = smiles
                                    st.session_state.selected_compound_name = selected_compound
                                    st.session_state.show_suggestions = False
                                    # Clear the dropdown selection for next time
                                    if 'compound_selector' in st.session_state:
                                        del st.session_state.compound_selector
                                    st.info("Proceeding with search...")
                                    st.rerun()
                            else:
                                st.error(f"Could not get SMILES for {selected_compound}")
                        else:
                            st.error(f"Could not find compound ID for {selected_compound}")

        with search_tab2:
            # Create two columns for the input fields
            input_col1, input_col2 = st.columns([3, 1])
            
            # Add SMILES input in the first column
            with input_col1:
                # Check if we should load demo SMILES
                if 'demo_smiles' in st.session_state and st.session_state.demo_smiles:
                    default_value = st.session_state.demo_smiles
                    # Clear the demo flag after using it
                    st.session_state.demo_smiles = ""
                elif 'selected_smiles' in st.session_state and st.session_state.selected_smiles:
                    default_value = st.session_state.selected_smiles
                    st.session_state.selected_smiles = ""  # Clear after using
                else:
                    default_value = st.session_state.get('smiles_input', '')
                
                smiles_input = st.text_input(
                    "Enter SMILES:",
                    placeholder="Enter a valid SMILES string",
                    value=default_value,
                    key="smiles_input_field"
                )
                smiles_input = smiles_input.strip()

            # Add min_count input in the second column
            with input_col2:
                min_count = st.number_input("Min frequency:", min_value=3, max_value=100, value=3, step=1,
                                            help="Minimum frequency of a conjugation observed in public LC-MS/MS datasets to be included.",
                                            format="%d")
                
            # Create a row for the buttons
            button_col1, _, button_col2, _ = st.columns([2, 1, 2, 7])
            
            # Add the Search button in the first column
            with button_col1:
                search_button = st.button(
                    "**Search**", 
                    type="primary", 
                    use_container_width=True,
                    icon=":material/search:",
                    help="Search for conjugated metabolites",
                    key="search_smiles_button"
                )
            
            # Add the Demo button in the second column (to the right of Search)
            with button_col2:
                demo_button = st.button(
                    "**Load Demo**", 
                    type="secondary", 
                    use_container_width=True,
                    icon=":material/login:",
                    help="Load example: Tyrosine",
                    key="demo_button"
                )

            # Process the demo button
            if demo_button:
                # Set demo SMILES in session state and trigger rerun
                st.session_state.demo_smiles = DEMO_SMILES
                st.rerun()

        # ===== MAIN SEARCH LOGIC (OUTSIDE THE TABS) =====
        # Determine effective input for processing (common for both tabs)
        effective_smiles = ""
        effective_min_count = 3
        
        # Handle name search selection
        if st.session_state.get('selected_smiles'):
            effective_smiles = st.session_state.selected_smiles
            effective_min_count = min_count_name
            # Clear the selected SMILES after using it
            st.session_state.selected_smiles = ""
        # Handle SMILES search
        elif search_button and smiles_input:
            effective_smiles = smiles_input
            effective_min_count = min_count

        # Process the search (common for both tabs)
        if effective_smiles:
            # Convert SMILES to formula and InChIKey
            formula, inchikey = smiles_to_formula_inchikey(effective_smiles)
            
            if inchikey is None:
                st.error("❌ Invalid SMILES string. Please double-check your input.")
            else:
                ###################
                # Debug information
                st.write(f"DEBUG: effective_smiles = {effective_smiles}")
                st.write(f"DEBUG: effective_min_count = {effective_min_count}")
                st.write(f"DEBUG: inchikey = {inchikey}")
                st.write(f"DEBUG: inchikey_14 = {inchikey[:14]}")
                st.write(f"DEBUG: mono_mass = {calc_monoisotopic_mass(formula)}")
                
                # Add to search history when a search is performed
                add_to_search_history(effective_smiles)

                # Create container for results
                results_container = st.container()
                
                with results_container:
                    st.subheader("Compound Information")
                    
                    # Results section with columns for info and structure
                    _, structure_col, _, info_col, description_col = st.columns([1, 3, 1, 7, 7])
                    
                    with structure_col:
                        common_names = inchikey_to_common_name(inchikey)
                        # Use selected compound name if available
                        if st.session_state.get('selected_compound_name'):
                            image_caption = st.session_state.selected_compound_name
                            st.session_state.selected_compound_name = ""  # Clear after using
                        else:
                            image_caption = common_names[0] if common_names else "Chemical structure"
                        
                        # Display the chemical structure image
                        image_url = get_structure_image_gnps2(effective_smiles)
                        st.image(image_url, caption=image_caption, use_container_width=True)
                    
                    with info_col:  
                        if common_names:
                            names_str = ', '.join(common_names[:3])
                            st.markdown(f"**Common names:** {names_str}")
                        # Use st.code to display SMILES as plain text without Markdown interpretation
                        st.markdown(f"**SMILES:**\n```\n{effective_smiles}\n```")
                        st.markdown(f"**Formula:** {formula}")
                        st.markdown(f"**InChIKey:** {inchikey}")
                        mono_mass = calc_monoisotopic_mass(formula)                
                        st.markdown(f"**Monoisotopic mass:** {mono_mass:.4f}")
                        
                    with description_col:
                        st.markdown("**Compound description (from PubChem):**")
                        # Fetch and display the description from PubChem
                        description = get_compound_description_pubchem(effective_smiles)
                        if description:
                            st.write(description)                    
                
                # Process search
                inchikey_14 = inchikey[:14]  # 2D InChIKey (first 14 characters)
                
                # Search both positive and negative modes
                df_filtered = filter_search_results(db_path, inchikey_14, mono_mass, effective_min_count)
                st.write(f"DEBUG: df_filtered shape = {df_filtered.shape}")

                # Display results
                if df_filtered.empty:
                    st.warning(f"❌ No matches found for SMILES: {effective_smiles}")
                else:
                    # total matches
                    total_matches = len(df_filtered)
                    st.info(f"✅ Found {total_matches} matches for the target SMILES.")
                    
                    # Add a bar chart showing the frequency of delta masses
                    if len(df_filtered) > 1:  # Only show chart if there are multiple results
                        st.subheader("Distribution of Conjugate Delta Masses")                       
                        
                        delta_mass_counts = prepare_delta_mass_plot(df_filtered)
                                                
                        # Create a column with specific width to control the chart size
                        _, chart_col, _ = st.columns([1, 12, 1])
                        
                        with chart_col:
                            # Create the bar chart in the middle column
                            chart = st.scatter_chart(
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
                    
                    # Select columns to display
                    df_filtered = df_filtered[['ion_polarity', 'annotation_type', 'count', 'conjugate_delta_mass', 'conjugate_name', 
                                               'mirror_plot_ref_1', 'mirror_plot_ref_2', 'match_type', 'masst']]
                                                                    
                    # After preparing df_filtered but before displaying it
                    if not df_filtered.empty:
                        st.subheader("Result Table")
                        
                        _, filter_col1, filter_col2, filter_col3, filter_col4, _ = st.columns([1, 2, 2, 2, 2, 1])
                        
                        # Filter by Ion polarity
                        with filter_col1:
                            available_polarities = ['All', '+', '-']
                            polarity_filter = st.selectbox('Ion polarity:', available_polarities)
                            
                        # Filter by Annotation type
                        with filter_col2:
                            available_annotations = ['All', 'spec_spec', 'spec_delta']
                            annotation_filter = st.selectbox('Annotation type:', available_annotations)
                        
                        # Filter by Conjugate name (has name or not)
                        with filter_col3:
                            name_filter = st.selectbox('Conjugate name:', ['All', 'With name (annotated)', 'Without name (unannotated)'])
                        
                        # Filter by Match type
                        with filter_col4:
                            available_matches = ['All', 'spec (ref 1)', 'spec (ref 2)', 'spec (ref 1) or spec (ref 2)', 'delta mass']
                            match_filter = st.selectbox('Match type:', available_matches, index=3)  # Default to 'spec (ref 1) or spec (ref 2)'
                        
                        # Apply the filters
                        filtered_results = df_filtered.copy()
                        
                        if polarity_filter != 'All':
                            filtered_results = filtered_results[filtered_results['ion_polarity'] == polarity_filter]
                            
                        if annotation_filter != 'All':
                            filtered_results = filtered_results[filtered_results['annotation_type'] == annotation_filter]
                            
                        if match_filter != 'All':
                            if match_filter == 'spec (ref 1) or spec (ref 2)':
                                filtered_results = filtered_results[filtered_results['match_type'] != 'delta mass']
                            else:
                                filtered_results = filtered_results[filtered_results['match_type'] == match_filter]
                            
                        if name_filter == 'With name (annotated)':
                            filtered_results = filtered_results[filtered_results['conjugate_name'].notna() & 
                                                                (filtered_results['conjugate_name'] != '')]
                        elif name_filter == 'Without name (unannotated)':
                            filtered_results = filtered_results[filtered_results['conjugate_name'].isna() | 
                                                                (filtered_results['conjugate_name'] == '')]

                        _, info_col, _ = st.columns([1, 8, 1])
                        with info_col:
                            # Show how many results are displayed after filtering
                            st.info(f"Showing {len(filtered_results)} of {len(df_filtered)} results")
                        
                        # Display the filtered dataframe with clickable links
                        st.dataframe(
                            filtered_results,
                            column_config={
                                "ion_polarity": st.column_config.TextColumn(
                                    "Ion polarity",
                                    width="small",
                                    help="Ionization polarity of the query MS/MS"
                                ),
                                "annotation_type": st.column_config.TextColumn(
                                    "Annotation type",
                                    width="small",
                                    help="How query MS/MS are annotated"
                                ),
                                "count": st.column_config.NumberColumn(
                                    "Count",
                                    width="small",
                                    help="Frequency in public LC-MS/MS datasets",
                                    format="%d"
                                ),
                                "conjugate_delta_mass": st.column_config.NumberColumn(
                                    "Conjugate delta mass",
                                    width="small",
                                    help="Mass of the conjugate component",
                                    format="%.2f"
                                ),
                                "conjugate_name": st.column_config.TextColumn(
                                    "Conjugate name",
                                    width="large",
                                    help="Name of the conjugate component"
                                ),
                                "mirror_plot_ref_1": st.column_config.LinkColumn(
                                    "Mirror plot (Ref 1)",
                                    width="small",
                                    help="Click to view mirror plot between query MS/MS and reference 1",
                                    display_text="View",
                                    required=False
                                ),
                                "mirror_plot_ref_2": st.column_config.LinkColumn(
                                    "Mirror plot (Ref 2)",
                                    width="small",
                                    help="Click to view mirror plot between query MS/MS and reference 2",
                                    display_text="View",
                                    required=False
                                ),
                                "match_type": st.column_config.TextColumn(
                                    "Match type",
                                    width="small",
                                    help="How the target compound is found"
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

        # Handle cases where no search is triggered
        elif not effective_smiles and (search_button or name_search_button):
            st.warning("⚠️ Please enter a SMILES string or search for a compound name.")
                            
        # Add final notes and copyright at the bottom of the page
        st.markdown("---")
        st.markdown("""
        ### Notes
        1. This web app does not include all conjugation results. For more comprehensive results, please refer to [our paper](https://doi.org/10.1101/2025.01.01.123456) and [Zenodo repository](https://zenodo.org/record/1234567). Due to memory limit, for each conjugation, we only reserve one representative query MS/MS and its corresponding reference MS/MS spectra in the result table.
        2. Reference spectra from [GNPS](https://external.gnps2.org/gnpslibrary) have USI links. Some spectra from [MassBank](https://github.com/MassBank/MassBank-data/releases) and all NIST20 spectra (commercially available) do not have USIs available. In spectral matches where MassBank or NIST20 spectra are involved, only the query MS/MS will be shown in the mirror plot viewer.
        3. All search results are based on 2D chemical structure.
        4. Column descriptions:
        - **Ion polarity**: The ion polarity of the query MS/MS.
        - **Annotation type**: how query MS/MS are annotated in the search results.
            - spec_spec: Query MS/MS is explained as a conjugate of two component molecules, and both components are explained by reference MS/MS via spectral matching.
            - spec_delta: Query MS/MS is explained as a conjugate of a reference MS/MS via spectral matching and a delta mass.
        - **Count**: The frequency of the conjugation in public LC-MS/MS datasets.
        - **Conjugate delta mass**: The mass of the conjugate component.
        - **Conjugate name**: The name of the conjugate component, if available.
        - **Mirror plot (Ref 1)**: Link to the mirror plot between the query MS/MS and reference 1.
        - **Mirror plot (Ref 2)**: Link to the mirror plot between the query MS/MS and reference 2.
        - **Match type**: how the target compound is found in the search results.
            - spec: Spectral match with a reference MS/MS.
            - delta: Delta mass match.
            
        © All rights reserved, Shipei Xing 2025
        """)


def get_db_path():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    db_path = os.path.join(script_dir, "conjugated_metabolome.db")
    if os.path.exists(db_path):
        return db_path
    else:
        st.error(f"❌ Database SQLite file not found: {db_path}. ")
        return None


def initialize_demo_smiles():
    """Initialize demo_smiles in session state"""
    if 'demo_smiles' not in st.session_state:
        st.session_state.demo_smiles = ""


def initialize_name_search():
    """Initialize name search related session state"""
    if 'name_suggestions' not in st.session_state:
        st.session_state.name_suggestions = []
    if 'show_suggestions' not in st.session_state:
        st.session_state.show_suggestions = False
    if 'selected_smiles' not in st.session_state:
        st.session_state.selected_smiles = ""
    if 'selected_compound_name' not in st.session_state:
        st.session_state.selected_compound_name = ""
    if 'demo_name' not in st.session_state:
        st.session_state.demo_name = ""


def initialize_search_history():
    """Initialize search history in session state"""
    if 'search_history' not in st.session_state:
        st.session_state.search_history = []


def add_to_search_history(smiles):
    """Add SMILES to search history if not already present"""
    if smiles and smiles not in st.session_state.search_history:
        st.session_state.search_history.insert(0, smiles)  # Add to beginning of list
        # Keep only the last 10 searches to avoid clutter
        if len(st.session_state.search_history) > 10:
            st.session_state.search_history = st.session_state[:10]


def display_search_history():
    """Display search history with clickable SMILES"""
    if st.session_state.search_history:
        st.markdown("### 📚 Recent Searches")
        
        # Create a container for the history
        history_container = st.container()
        
        with history_container:
            # Display each SMILES as a clickable button
            for i, hist_smiles in enumerate(st.session_state.search_history):
                col1, col2 = st.columns([5, 1])
                
                with col1:
                    # Create a button for each historical SMILES
                    if st.button(
                        f"🔍 {hist_smiles[:50]}{'...' if len(hist_smiles) > 50 else ''}", 
                        key=f"history_{i}",
                        help=f"Click to search: {hist_smiles}"
                    ):
                        # Set the SMILES input and trigger search
                        st.session_state.selected_smiles = hist_smiles
                        st.rerun()
                
                with col2:
                    # Add a small delete button for each entry
                    if st.button("❌", key=f"delete_{i}", help="Remove from history", type="secondary"):
                        st.session_state.search_history.pop(i)
                        st.rerun()
    else:
        st.info("No search history yet. Start searching to see your history here!")

    
if __name__ == "__main__":
    main()