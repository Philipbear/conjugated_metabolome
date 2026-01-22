import streamlit as st
from sql_utils import get_git_short_rev
from homepage import render_homepage
from compound_search import render_compound_search, render_search_history, initialize_search_session_state
from delta_mass_browser import render_delta_mass_browser


def initialize_session_state():
    """Initialize all session state variables at app start"""
    initialize_search_session_state()


def main():
    # Set the page configuration
    app_version = "2026-01-22"
    try:
        git_hash = get_git_short_rev()
    except:
        git_hash = "unknown"
    repo_link = "https://github.com/Philipbear/conjugated_metabolome"

    st.set_page_config(
        page_title="Conjugated Metabolome Explorer", 
        layout="wide",
        menu_items={
            "About": (f"**App Version**: {app_version} | "
                     f"[**Git Hash**: {git_hash}]({repo_link}/commit/{git_hash})")
        }
    )
    
    # Initialize navigation state
    if 'current_page' not in st.session_state:
        st.session_state.current_page = "Homepage"
    
    # Initialize session state variables for compound search
    initialize_session_state()
    
    # Render sidebar navigation
    render_sidebar(app_version)
    
    # Main content
    _, main_col, _ = st.columns([1, 8, 1])
    with main_col:
        # Route to appropriate page
        if st.session_state.current_page == "Homepage":
            render_homepage()
        elif st.session_state.current_page == "Compound Search":
            render_compound_search()
            # Footer for Compound Search page only
            render_footer()
        elif st.session_state.current_page == "Browse by Delta Mass":
            render_delta_mass_browser()


def render_sidebar(app_version):
    """Render the sidebar with navigation and info"""
    with st.sidebar:
        st.image("https://ccms-ucsd.github.io/GNPSDocumentation/img/logo/GNPS_logo_original_transparent.png", width=150)
        st.title("Conjugated Metabolome Explorer")
                
        st.markdown("---")
        
        # Navigation
        st.header("📑 Navigation")
        
        # Homepage button
        if st.button(
            "🏠 Homepage",
            use_container_width=True,
            type="primary" if st.session_state.current_page == "Homepage" else "secondary"
        ):
            st.session_state.current_page = "Homepage"
            st.rerun()
        
        # Compound Search button
        if st.button(
            "🔍 Compound Search",
            use_container_width=True,
            type="primary" if st.session_state.current_page == "Compound Search" else "secondary"
        ):
            st.session_state.current_page = "Compound Search"
            st.rerun()
        
        # Delta Mass Browser button
        if st.button(
            "⚖️ Browse by Delta Mass",
            use_container_width=True,
            type="primary" if st.session_state.current_page == "Browse by Delta Mass" else "secondary"
        ):
            st.session_state.current_page = "Browse by Delta Mass"
            st.rerun()
        
        st.markdown("---")
        
        # Show search history only on Compound Search page
        if st.session_state.current_page == "Compound Search":
            if 'search_history' in st.session_state:  # Check before accessing
                render_search_history()
                st.markdown("---")
        
        # About section
        st.markdown("""
        ### 📖 About
        This app allows you to explore the pan-repository conjugated metabolome and search for potential metabolite conjugations.
        
        ### 📧 Contact
        For questions or feedback, please contact Shipei Xing at [philipxsp@hotmail.com](mailto:philipxsp@hotmail.com).
        """)
        
        # App version info
        st.markdown(f"""
        ---
        **Version:** {app_version}
        """)


def render_footer():
    """Render the footer with notes"""
    st.markdown("---")
    st.markdown("""
    ### Notes
    1. This web application does not include all conjugation results across repositories. For more comprehensive results, please refer to our paper.
    2. Mirror plots will be correctly shown when reference MS/MS spectra with valid universal spectrum identifiers (USIs).
    3. All search results are based on 2D chemical structure.
    """)


if __name__ == "__main__":
    main()
