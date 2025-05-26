import streamlit as st
import pandas as pd
import os
from utils import filter_by_inchikey, add_mirror_plot_urls
from chem_utils import smiles_to_formula_inchikey, calc_monoisotopic_mass, inchikey_to_common_name, get_structure_image_pubchem, get_compound_description_pubchem


st.set_page_config(page_title="Conjugated Metabolome Explorer", layout="wide")

# Left panel with app info
with st.sidebar:
    st.title("Conjugated Metabolome Explorer")
    st.image("https://ccms-ucsd.github.io/GNPSDocumentation/img/logo/GNPS_logo_original_transparent.png", width=150)
    
    st.markdown("""
    ### About
    This app allows you to explore potential metabolite conjugations using SMILES strings.
    
        
    ### Note
    - This webpage does not include all conjugation results. For more comprehensive results, please refer to [our paper](https://doi.org/10.1101/2025.01.01.123456) and [Zenodo repository](https://zenodo.org/record/1234567).
    - If the reference spectra are not from GNPS or MassBank, only the query MS/MS will be shown in the mirror plot viewer.
    - Due to memory usage, for each conjugation, we only reserve one representative query MS/MS and its corresponding reference MS/MS spectra in the result table.
        
    ### Citation
    Please use it responsibly and cite [our work](https://doi.org/10.1101/2025.01.01.123456) if you find it useful:
    - Xing S, et al. "Conjugated Metabolome..."  
      <https://doi.org/10.1101/2025.01.01.123456>  
      (Accessed: 2025-01-01)
    
    ### Contact
    For questions or feedback, please contact Shipei Xing at
    [philipxsp@hotmail.com](mailto:philipxsp@hotmail.com)
    """)

# Create a layout
col1, col2, col3 = st.columns([1, 7, 1])
# Main panel with search functionality
with col2:
    # Load the data
    @st.cache_data
    def load_data():
        try:
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
            
        except Exception as e:
            st.info(f"❌ Error loading data: {str(e)}")
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

    # Add filter dropdown in the third column
    with input_col2:
        match_filter = st.selectbox(
            "Target compound is found by:",
            ["Spectral match", "Delta mass", "Spectral match or delta mass"],
            help="Select the type of match to filter results. 'Spectral match' requires a high MS/MS similarity score, 'Delta mass' filters based on mass difference.",
            index=2
        )
    
    # Add min_count input in the second column
    with input_col3:
        min_count = st.number_input("Min frequency:", min_value=1, max_value=100, value=3, step=1,
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
                    info_col, structure_col, description_col = st.columns([3, 1, 3])
                    
                    with structure_col:
                        common_names = inchikey_to_common_name(inchikey)
                        image_caption = common_names[0] if common_names else "Chemical structure"
                        
                        # Display the chemical structure image
                        image_url = get_structure_image_pubchem(smiles_input)
                        st.image(image_url, caption=image_caption, use_container_width=True)
                    
                    with info_col:                        
                        # Use st.code to display SMILES as plain text without Markdown interpretation
                        st.markdown(f"**SMILES:**\n```\n{smiles_input}\n```")
                        if common_names:
                            names_str = ', '.join(common_names[:3])
                            st.markdown(f"**Common names:** {names_str}")
                        st.markdown(f"**Formula:** {formula}")
                        st.markdown(f"**InChIKey:** {inchikey}")                        
                        st.markdown(f"**Monoisotopic mass:** {calc_monoisotopic_mass(formula):.4f}")
                        
                    with description_col:
                        # Fetch and display the description from PubChem
                        description = get_compound_description_pubchem(smiles_input)
                        if description:
                            st.write(description)                    
                
                # Process search
                inchikey_14 = inchikey[:14]  # Use the first 14 characters of the InChIKey
                
                # Search both positive and negative modes
                pos_filtered = filter_by_inchikey(pos_df, ms2db_df, inchikey_14, min_count)
                if not pos_filtered.empty:
                    pos_filtered['Ion polarity'] = '+'
                
                neg_filtered = filter_by_inchikey(neg_df, ms2db_df, inchikey_14, min_count)
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
                    
                    # Add a bar chart showing the frequency of delta masses
                    if len(df_filtered) > 1:  # Only show chart if there are multiple results
                        st.subheader("Distribution of Delta Masses")
                        
                        # Group by delta_mass and sum the counts
                        delta_mass_counts = df_filtered.groupby('delta_mass')['count'].sum().reset_index()
                        delta_mass_counts.columns = ['Delta mass', 'Total count']
                        delta_mass_counts = delta_mass_counts.sort_values('Delta mass')
                        
                        # Create a column with specific width to control the chart size
                        chart_col1, chart_col2, chart_col3 = st.columns([1, 10, 1])
                        
                        with chart_col2:
                            # Create the bar chart in the middle column
                            chart = st.bar_chart(
                                data=delta_mass_counts,
                                x='Delta mass',
                                y='Total count',
                                use_container_width=True
                            )
                    
                    # Add mirror plot URLs
                    df_filtered = add_mirror_plot_urls(df_filtered)
                    
                    # Select columns to display
                    display_cols = ['Ion polarity', 'count', 'delta_mass', 'Conjugate name',
                                'Mirror plot (Ref 1)', 'Mirror plot (Ref 2)', 'Match type']
                    df_filtered = df_filtered[display_cols]
                    # Rename columns for clarity
                    df_filtered = df_filtered.rename(columns={
                        'count': 'Count',
                        'delta_mass': 'Delta mass'
                    })
                    
                    # Display results
                    st.subheader(f"Result Table")
                    # Display the dataframe with clickable links
                    st.dataframe(
                        df_filtered,
                        column_config={
                            "Mirror plot (Ref 1)": st.column_config.LinkColumn(
                                "Mirror plot (Ref 1)",
                                width="medium",
                                help="Click to view mirror plot between query MS/MS and reference 1",
                                display_text="View",
                                required=False
                            ),
                            "Mirror plot (Ref 2)": st.column_config.LinkColumn(
                                "Mirror plot (Ref 2)",
                                width="medium",
                                help="Click to view mirror plot between query MS/MS and reference 2",
                                display_text="View",
                                required=False
                            ),
                        },
                        hide_index=True,
                    )