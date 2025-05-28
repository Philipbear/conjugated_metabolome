import streamlit as st
import pandas as pd
import os
from utils import filter_search_results, prepare_delta_mass_plot, add_mirror_plot_urls
from chem_utils import smiles_to_formula_inchikey, calc_monoisotopic_mass, inchikey_to_common_name, get_structure_image_pubchem, get_structure_image_gnps2, get_compound_description_pubchem


st.set_page_config(page_title="Conjugated Metabolome Explorer", layout="wide")

# Left panel with app info
with st.sidebar:
    st.title("Conjugated Metabolome Explorer")
    st.image("https://ccms-ucsd.github.io/GNPSDocumentation/img/logo/GNPS_logo_original_transparent.png", width=150)
    
    st.markdown("""
    ### About
    This app allows you to explore potential metabolite conjugations using SMILES strings.
    
    This webpage does not include all conjugation results. For more comprehensive results, please refer to [our paper](https://doi.org/10.1101/2025.01.01.123456) and [Zenodo repository](https://zenodo.org/record/1234567).
        
    ### Citation
    Please use it responsibly and cite [our work](https://doi.org/10.1101/2025.01.01.123456) if you find it useful:
    - S. Xing, V. Charron-Lamoureux, A. Patan, ..... [Navigating the underexplored conjugated metabolome](https://doi.org/10.1101/2025.01.01.123456). 2025
    
    ### Contact
    For questions or feedback, please contact Shipei Xing at
    [philipxsp@hotmail.com](mailto:philipxsp@hotmail.com)
    """)

# Create a layout
col1, col2, col3 = st.columns([1, 8, 1])
# Main panel with search functionality
with col2:
    # Load the data
    @st.cache_data
    def load_data():
        script_dir = os.path.dirname(os.path.abspath(__file__))
        pos_path = os.path.join(script_dir, "pos_refined.parquet")
        if os.path.exists(pos_path):
            # st.info(f"✅ Loading files from script directory: {script_dir}")
            pos_df = pd.read_parquet(pos_path)
            neg_df = pd.read_parquet(os.path.join(script_dir, "neg_refined.parquet"))
            ms2db_df = pd.read_parquet(os.path.join(script_dir, "ms2db.parquet"))
            return pos_df, neg_df, ms2db_df
        else:
            st.error(f"❌ Data files not found in script directory: {script_dir}. ")
            return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()            

    # Load the data
    pos_df, neg_df, ms2db_df = load_data()
    
    st.title("Conjugated Metabolome Explorer")
    # Create a search section
    st.header("Search")
    
    # Create two columns for the input fields
    input_col1, input_col2, input_col3 = st.columns([3, 1, 1])
    
    # Initialize demo SMILES in session state if demo button is clicked
    if 'demo_smiles' not in st.session_state:
        st.session_state.demo_smiles = ""
    
    # Add SMILES input in the first column
    with input_col1:
        smiles_input = st.text_input("Enter a SMILES string:", value=st.session_state.demo_smiles)
        # strip
        smiles_input = smiles_input.strip()

    # Add filter dropdown in the third column
    with input_col2:
        match_filter = st.selectbox(
            "Target compound is found by:",
            ["Spectral match", "Delta mass", "Spectral match or delta mass"],
            help="Select the type of match to filter results. 'Spectral match' requires a high MS/MS similarity score, 'Delta mass' filters based on mass difference.",
            index=2
        )
        # map match_filter to integer
        if match_filter == "Spectral match":
            match_filter_ls = ['spec']
        elif match_filter == "Delta mass":
            match_filter_ls = ['delta']
        else:
            match_filter_ls = ['spec', 'delta']
    
    # Add min_count input in the second column
    with input_col3:
        min_count = st.number_input("Min frequency:", min_value=2, max_value=100, value=2, step=1,
                                    help="Minimum frequency of a conjugation in public LC-MS/MS datasets to be included in the results.",
                                    format="%d")
        
    # Create a row for the buttons
    button_col1, button_col2, button_col3 = st.columns([1, 1, 8])
    
    # Add the Search button in the first column
    with button_col1:
        search_button = st.button("Search")
    
    # Add the Demo button in the second column (to the right of Search)
    with button_col2:
        demo_button = st.button("Load Demo")

    # Process the demo button
    if demo_button:
        # Set demo SMILES in session state and trigger rerun
        st.session_state.demo_smiles = 'C1=CC=C(C=C1)C[C@@H](C(=O)O)N'
        st.rerun()
    
    # Clear demo SMILES from session state if user manually changes the input
    if smiles_input != st.session_state.demo_smiles:
        st.session_state.demo_smiles = ""

    # Process the input when the user submits
    if search_button or smiles_input:
        if not smiles_input:
            st.error("❌ Please enter a SMILES string.")
        else:
            # Convert SMILES to formula and InChIKey
            formula, inchikey = smiles_to_formula_inchikey(smiles_input)
            
            if inchikey is None:
                st.error("❌ Invalid SMILES string. Please double-check your input.")
            else:
                # Create container for results
                results_container = st.container()
                
                with results_container:
                    st.subheader("Compound Information")
                    
                    # Results section with columns for info and structure
                    _, structure_col, _, info_col, description_col = st.columns([1, 3, 1, 7, 7])
                    
                    with structure_col:
                        common_names = inchikey_to_common_name(inchikey)
                        image_caption = common_names[0] if common_names else "Chemical structure"
                        
                        # Display the chemical structure image
                        image_url = get_structure_image_gnps2(smiles_input)
                        st.image(image_url, caption=image_caption, use_container_width=True)
                    
                    with info_col:  
                        if common_names:
                            names_str = ', '.join(common_names[:3])
                            st.markdown(f"**Common names:** {names_str}")
                        # Use st.code to display SMILES as plain text without Markdown interpretation
                        st.markdown(f"**SMILES:**\n```\n{smiles_input}\n```")
                        st.markdown(f"**Formula:** {formula}")
                        st.markdown(f"**InChIKey:** {inchikey}")
                        mono_mass = calc_monoisotopic_mass(formula)                
                        st.markdown(f"**Monoisotopic mass:** {mono_mass:.4f}")
                        
                    with description_col:
                        # Fetch and display the description from PubChem
                        description = get_compound_description_pubchem(smiles_input)
                        if description:
                            st.write(description)                    
                
                # Process search
                inchikey_14 = inchikey[:14]  # Use the first 14 characters of the InChIKey
                
                # Search both positive and negative modes
                pos_filtered = filter_search_results(pos_df, ms2db_df, inchikey_14, mono_mass, min_count, match_filter_ls)
                if not pos_filtered.empty:
                    pos_filtered['Ion polarity'] = '+'
                
                neg_filtered = filter_search_results(neg_df, ms2db_df, inchikey_14, mono_mass, min_count, match_filter_ls)
                if not neg_filtered.empty:
                    neg_filtered['Ion polarity'] = '-'
                
                # Display results
                if pos_filtered.empty and neg_filtered.empty:
                    st.warning(f"❌ No matches found for SMILES: {smiles_input}")
                else:
                    # total matches
                    total_matches = len(pos_filtered) + len(neg_filtered)                    
                    st.info(f"✅ Found {total_matches} matches for the target SMILES.")                 
                    
                    # Concatenate the filtered DataFrames
                    df_filtered = pd.concat([pos_filtered, neg_filtered], ignore_index=True)
                    # Rename columns for clarity
                    df_filtered = df_filtered.rename(columns={
                        'count': 'Count',
                        'delta_mass': 'Conjugate delta mass'
                    })
                    
                    # Add a bar chart showing the frequency of delta masses
                    if len(df_filtered) > 1:  # Only show chart if there are multiple results
                        st.subheader("Distribution of Conjugate Delta Masses")                       
                        
                        delta_mass_counts = prepare_delta_mass_plot(df_filtered)
                                                
                        # Create a column with specific width to control the chart size
                        chart_col1, chart_col2, chart_col3 = st.columns([1, 12, 1])
                        
                        with chart_col2:
                            # Create the bar chart in the middle column
                            chart = st.scatter_chart(
                                data=delta_mass_counts,
                                x='Conjugate delta mass',
                                y='Count',
                                x_label='Conjugate delta mass (Da)',
                                y_label='Observed frequency',
                                size='Count',
                                color='Ion polarity',
                                height=450,
                                use_container_width=True
                            )
                    
                    # Add mirror plot URLs
                    df_filtered = add_mirror_plot_urls(df_filtered)
                    
                    # Select columns to display
                    df_filtered = df_filtered[['Ion polarity', 'Annotation type', 'Count', 'Conjugate delta mass', 'Conjugate name', 
                                               'Mirror plot (Ref 1)', 'Mirror plot (Ref 2)', 'Match type']]
                    
                    
                    
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
                            name_filter = st.selectbox('Conjugate name:', ['All', 'Has name', 'No name'])
                        
                        # Filter by Match type
                        with filter_col4:
                            available_matches = ['All', 'spec (ref 1)', 'spec (ref 2)', 'spec (ref 1) or spec (ref 2)', 'delta']
                            match_filter = st.selectbox('Match type:', available_matches)
                        
                        
                        # Apply the filters
                        filtered_results = df_filtered.copy()
                        
                        if polarity_filter != 'All':
                            filtered_results = filtered_results[filtered_results['Ion polarity'] == polarity_filter]
                            
                        if annotation_filter != 'All':
                            filtered_results = filtered_results[filtered_results['Annotation type'] == annotation_filter]
                            
                        if match_filter != 'All':
                            if match_filter == 'spec (ref 1) or spec (ref 2)':
                                filtered_results = filtered_results[(filtered_results['Match type'] == 'spec (ref 1)') | 
                                                                   (filtered_results['Match type'] == 'spec (ref 2)')]
                            else:
                                filtered_results = filtered_results[filtered_results['Match type'] == match_filter]
                            
                        if name_filter == 'Has name':
                            filtered_results = filtered_results[filtered_results['Conjugate name'].notna() & 
                                                                (filtered_results['Conjugate name'] != '')]
                        elif name_filter == 'No name':
                            filtered_results = filtered_results[filtered_results['Conjugate name'].isna() | 
                                                                (filtered_results['Conjugate name'] == '')]
                        
                        _, info_col, _ = st.columns([1, 8, 1])
                        with info_col:
                            # Show how many results are displayed after filtering
                            st.info(f"Showing {len(filtered_results)} of {len(df_filtered)} results")
                        
                        # # Display the filtered dataframe
                        # st.subheader(f"Result Table")
                        
                        # Display the filtered dataframe with clickable links
                        st.dataframe(
                            filtered_results,
                            column_config={
                                "Ion polarity": st.column_config.TextColumn(
                                    "Ion polarity",
                                    width="small",
                                    help="Ionization polarity of the query MS/MS"
                                ),
                                "Annotation type": st.column_config.TextColumn(
                                    "Annotation type",
                                    width="small",
                                    help="How query MS/MS are annotated"
                                ),
                                "Count": st.column_config.NumberColumn(
                                    "Count",
                                    width="small",
                                    help="Frequency in public LC-MS/MS datasets",
                                    format="%d"
                                ),
                                "Conjugate delta mass": st.column_config.NumberColumn(
                                    "Conjugate delta mass",
                                    width="small",
                                    help="Mass of the conjugate component",
                                    format="%.2f"
                                ),
                                "Conjugate name": st.column_config.TextColumn(
                                    "Conjugate name",
                                    width="large",
                                    help="Name of the conjugate component"
                                ),
                                "Mirror plot (Ref 1)": st.column_config.LinkColumn(
                                    "Mirror plot (Ref 1)",
                                    width="small",
                                    help="Click to view mirror plot between query MS/MS and reference 1",
                                    display_text="View",
                                    required=False
                                ),
                                "Mirror plot (Ref 2)": st.column_config.LinkColumn(
                                    "Mirror plot (Ref 2)",
                                    width="small",
                                    help="Click to view mirror plot between query MS/MS and reference 2",
                                    display_text="View",
                                    required=False
                                ),
                                "Match type": st.column_config.TextColumn(
                                    "Match type",
                                    width="small",
                                    help="How the target compound is found"
                                ),
                            },
                            hide_index=True,
                            use_container_width=True
                        )
                    
    # Add final notes and copyright at the bottom of the page
    st.markdown("---")
    st.markdown("""
    ### Notes
    1. This app does not include all conjugation results. For more comprehensive results, please refer to [our paper](https://doi.org/10.1101/2025.01.01.123456) and [Zenodo repository](https://zenodo.org/record/1234567). Due to memory limit, for each conjugation, we only reserve one representative query MS/MS and its corresponding reference MS/MS spectra in the result table.
    2. Reference spectra from [MassBank](https://github.com/MassBank/MassBank-data/releases) and NIST20 (commercially available) do not have USIs. In spectral matches where MassBank or NIST20 spectra are involved, only the query MS/MS will be shown in the mirror plot viewer.
    3. All search results are based on 2D chemical structure.
    4. Column descriptions:
    - **Ion polarity**: The ion polarity of the query MS/MS.
    - **Annotation type**: how query MS/MS are annotated in the search results.
        - spec_spec: Query MS/MS is explained as a conjugate of two component molecules, and both components are explained by reference MS/MS via spectral matching.
        - spec_delta: Query MS/MS is explained as a conjugate of a reference MS/MS via spectral matching and a delta mass (one component unannotated).
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