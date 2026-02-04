import streamlit as st


def render_homepage():
    """Render the homepage with project information"""
    
    # Title
    st.title("Pan-repository conjugated metabolome navigation")
    
    st.markdown("---")
    
    # Introduction
    st.markdown("""
    <div style="font-size: 1.1em; line-height: 1.6;">
    Life's chemical diversity far exceeds current biochemical maps. While metabolomics has catalogued tens of thousands of small molecules, <strong>conjugated metabolites</strong>, formed when two or more molecular entities fuse through amidation, esterification, or related linkages, remain largely unexplored. These molecules can act as microbial signals, detoxification intermediates, or endogenous regulators, yet their global diversity is poorly characterized.
    </div>
    """, unsafe_allow_html=True)
    
    st.write("")  # Blank line
    
    st.markdown("""
    <div style="font-size: 1.1em; line-height: 1.6;">
    Here, we mined 1.32 billion MS/MS spectra across public metabolomics repositories using reverse spectral searching and delta-mass analysis to infer conjugation events. We generated structural hypotheses for 24,227,439 spectra, encompassing 217,291 substructure pairs with dual spectral support and 3,412,720 additional candidates with single-match support. Predictions include host–microbe co-metabolites, diet-derived conjugates, and drug-derived species, including drug-ethanolamine and creatinine conjugates that may alter biological activity, and reveal steroid-phosphoethanolamine conjugates. We synthesized and confirmed 55 conjugates by MS/MS, 27 of which were validated by retention time. These results provide a pan-repository map of conjugation chemistry, establish a resource for structural discovery, and offer a framework to further explore the potential scale and diversity of the conjugated metabolome.
    </div>
    """, unsafe_allow_html=True)
    
    st.markdown("---")
    
    # for images
    _, image_col, _ = st.columns([8, 84, 8])
    
    with image_col:
        # Workflow image
        st.image("https://raw.githubusercontent.com/Philipbear/conjugated_metabolome/main/gitfigs/workflow.svg", 
                use_container_width=True)
        st.markdown("<p style='text-align: center;'>Workflow for pan-repository conjugated metabolome discovery</p>", unsafe_allow_html=True)
    
    st.write("")  # Blank line
    
    # multi-column layout for additional images (also centered)
    _, col1, col2, col3, _ = st.columns([8, 38.5, 18.5, 27, 8])
    
    with col1:
        st.image("https://raw.githubusercontent.com/Philipbear/conjugated_metabolome/main/gitfigs/revcos.svg",
                use_container_width=True)
    
    with col2:
        st.image("https://raw.githubusercontent.com/Philipbear/conjugated_metabolome/main/gitfigs/class.svg",
                use_container_width=True)
    
    with col3:
        st.image("https://raw.githubusercontent.com/Philipbear/conjugated_metabolome/main/gitfigs/distributions.png",
                use_container_width=True)
    
    _, col1, col2, col3, _ = st.columns([8, 38.5, 18.5, 27, 8])
    with col1:
        st.markdown("<p style='text-align: center;'>Reverse cosine similarity analysis</p>", unsafe_allow_html=True)
    
    with col2:
        st.markdown("<p style='text-align: center;'>Chemical class linkages</p>", unsafe_allow_html=True)
    
    with col3:
        st.markdown("<p style='text-align: center;'>Broad distributions across biological systems</p>", unsafe_allow_html=True)
    
    st.markdown("---")
    
    # Citation
    st.header("📄 Citation")
    st.info("""
    **S. Xing et al.** Navigation of the conjugated metabolome. 
    https://github.com/Philipbear/conjugated_metabolome
    """)
    
    # License
    st.header("⚖️ License")
    st.markdown("This work is licensed under the **Apache License 2.0**.")
    
    # Contributors
    st.header("👥 Contributors")
    st.markdown("""
    - **Shipei Xing, PhD** - University of California San Diego
    - **Wilhan Nunes, PhD** - University of California San Diego
    - **Mingxun Wang, PhD** - University of California Riverside
    """)
    
    # Contact
    st.markdown("""
    ### 📧 Contact
    For questions or feedback, please contact Shipei Xing at [philipxsp@hotmail.com](mailto:philipxsp@hotmail.com).
    """)
    
    st.markdown("© All rights reserved, Shipei Xing 2026")

