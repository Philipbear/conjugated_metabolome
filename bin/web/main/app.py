import streamlit as st
import pandas as pd
import os
from utils import filter_by_inchikey, add_mirror_plot_urls
from chem_utils import smiles_to_formula_inchikey, calc_monoisotopic_mass, inchikey_to_common_name, get_structure_image_pubchem


st.set_page_config(page_title="Conjugated Metabolome Explorer", layout="wide")

# Left panel with app info
with st.sidebar:
    st.title("Conjugated Metabolome Explorer")
    st.image("https://ccms-ucsd.github.io/GNPSDocumentation/img/logo/GNPS_logo_original_transparent.png", width=150)
    
    st.markdown("""
    ### About
    This webpage allows you to search for potential metabolite conjugations using SMILES strings.
    
        
    ### Note
    - This webpage does not cover all conjugation results. For more comprehensive results, please refer to the [original publication](https://doi.org/10.1101/2025.01.01.123456) and [Zenodo repository](https://zenodo.org/record/1234567).
    - If the reference spectra are not from GNPS or MassBank, the mirror plots will only show the query MS/MS.
    - For each conjugation, we only show one representative query MS/MS and its corresponding reference MS/MS spectra for mirror plots.
        
    ### References
    Please use it responsibly and cite the [original work](https://doi.org/10.1101/2025.01.01.123456) if you find it useful:
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
            # Try different path strategies and report which one works
            
            # Strategy 1: Current directory
            if os.path.exists("pos_refined.parquet"):
                st.success("✅ Loading files from current directory")
                pos_df = pd.read_parquet("pos_refined.parquet")
                neg_df = pd.read_parquet("neg_refined.parquet")
                ms2db_df = pd.read_parquet("ms2db.parquet")
                return pos_df, neg_df, ms2db_df
            
            # Strategy 2: Script directory
            script_dir = os.path.dirname(os.path.abspath(__file__))
            pos_path = os.path.join(script_dir, "pos_refined.parquet")
            if os.path.exists(pos_path):
                st.success(f"✅ Loading files from script directory: {script_dir}")
                pos_df = pd.read_parquet(pos_path)
                neg_df = pd.read_parquet(os.path.join(script_dir, "neg_refined.parquet"))
                ms2db_df = pd.read_parquet(os.path.join(script_dir, "ms2db.parquet"))
                return pos_df, neg_df, ms2db_df
            
            # Strategy 3: Look in parent directories
            for i in range(3):  # Check up to 3 levels up
                test_dir = script_dir
                for _ in range(i):
                    test_dir = os.path.dirname(test_dir)
                
                pos_path = os.path.join(test_dir, "pos_refined.parquet")
                if os.path.exists(pos_path):
                    st.success(f"✅ Loading files from parent directory (level {i}): {test_dir}")
                    pos_df = pd.read_parquet(pos_path)
                    neg_df = pd.read_parquet(os.path.join(test_dir, "neg_refined.parquet"))
                    ms2db_df = pd.read_parquet(os.path.join(test_dir, "ms2db.parquet"))
                    return pos_df, neg_df, ms2db_df
            
            # If we get here, we couldn't find the files
            st.error("❌ Could not locate parquet files. Please check file locations.")
            st.write(f"Script directory: {script_dir}")
            st.write(f"Current working directory: {os.getcwd()}")
            st.write(f"Files in script directory: {os.listdir(script_dir) if os.path.exists(script_dir) else 'Directory not found'}")
            
            # Return empty dataframes as fallback
            return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
            
        except Exception as e:
            st.error(f"❌ Error loading data: {str(e)}")
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
    
    # Add min_count input in the second column
    with input_col2:
        min_count = st.number_input("Min frequency:", min_value=1, max_value=100, value=3, step=1)

    # Create a row for the buttons
    button_col1, button_col2, button_col3 = st.columns([1, 1, 6])
    
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
            st.error("Please enter a SMILES string.")
        else:
            # Convert SMILES to formula and InChIKey
            formula, inchikey = smiles_to_formula_inchikey(smiles_input)
            
            if inchikey is None:
                st.error("Invalid SMILES string. Could not generate InChIKey.")
            else:
                # Create container for results
                results_container = st.container()
                
                with results_container:
                    st.subheader("Results")
                    
                    # Results section with columns for info and structure
                    info_col, structure_col = st.columns([2, 1])
                    
                    with structure_col:
                        common_names = inchikey_to_common_name(inchikey)
                        image_caption = common_names[0] if common_names else "Chemical structure"
                        
                        # Display the chemical structure image
                        image_url = get_structure_image_pubchem(smiles_input)
                        st.image(image_url, caption=image_caption, width=250)
                    
                    with info_col:
                        # Display the chemical info
                        if common_names:
                            names_str = ', '.join(common_names[:3])
                            st.markdown(f"**Common names:** {names_str}")
                        
                        # Use st.code to display SMILES as plain text without Markdown interpretation
                        st.markdown("**SMILES:**")
                        st.code(smiles_input, language=None)
                        
                        st.markdown(f"**InChIKey:** {inchikey}; **Formula:** {formula}; **Monoisotopic mass:** {calc_monoisotopic_mass(formula):.4f}")
                
                # Process search
                results = []
                inchikey_14 = inchikey[:14]  # Use the first 14 characters of the InChIKey
                
                # Search both positive and negative modes
                pos_filtered = filter_by_inchikey(pos_df, ms2db_df, inchikey_14, min_count)
                if pos_filtered is not None:
                    pos_filtered['Ion polarity'] = '+'
                    results.append(pos_filtered)
                
                neg_filtered = filter_by_inchikey(neg_df, ms2db_df, inchikey_14, min_count)
                if neg_filtered is not None:
                    neg_filtered['Ion polarity'] = '-'
                    results.append(neg_filtered)
                
                # Display results
                if not results:
                    st.warning(f"No matches found for SMILES: {smiles_input}")
                else:
                    # Combine results
                    df_filtered = pd.concat(results, ignore_index=True)
                    
                    st.info(f"Found {len(df_filtered)} matches for SMILES: {smiles_input}")
                    
                    # Add mirror plot URLs
                    df_filtered = add_mirror_plot_urls(df_filtered)                    
                    
                    # Add a bar chart showing the frequency of delta masses
                    if len(df_filtered) > 1:  # Only show chart if there are multiple results
                        st.subheader("Distribution of Delta Masses")
                        
                        # Group by delta_mass and sum the counts
                        delta_mass_counts = df_filtered.groupby('delta_mass')['count'].sum().reset_index()
                        delta_mass_counts.columns = ['Delta mass', 'Total count']
                        delta_mass_counts = delta_mass_counts.sort_values('Delta mass')
                        
                        # Create a column with specific width to control the chart size
                        chart_col1, chart_col2, chart_col3 = st.columns([1, 8, 1])
                        
                        with chart_col2:
                            # Create the bar chart in the middle column
                            chart = st.bar_chart(
                                data=delta_mass_counts,
                                x='Delta mass',
                                y='Total count',
                                use_container_width=True
                            )
                    
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
                    st.subheader(f"Table of matches:")
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