import streamlit as st
import session
import ui

st.set_page_config(
    page_title="COBRApy GUI",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

def main():
    """Main application entry point."""
    session.init_session_state()
    
    ui.render_intro_content()

    st.sidebar.header("Choose Analysis Mode")
    st.session_state.analysis_mode = st.sidebar.radio(
        "Select Mode",
        ["Single Model Analysis", "Multi Model Analysis"],
        key="analysis_mode_selector"
    )

    if st.session_state.analysis_mode == "Single Model Analysis":
        ui.single_model_analysis_tabs()
    elif st.session_state.analysis_mode == "Multi Model Analysis":
        ui.multi_model_analysis_tabs()

if __name__ == "__main__":
    main()
