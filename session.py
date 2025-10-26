import streamlit as st

def init_session_state():
    """Initialize all session state variables."""
    if 'models' not in st.session_state:
        st.session_state.models = {}
    if 'current_model_name' not in st.session_state:
        st.session_state.current_model_name = None
    if 'fba_solutions' not in st.session_state:
        st.session_state.fba_solutions = {}
    if 'analysis_mode' not in st.session_state:
        st.session_state.analysis_mode = "Single Model Analysis"
    if 'summary_display_counter' not in st.session_state:
        st.session_state.summary_display_counter = 0
