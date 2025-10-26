import streamlit as st
import cobra
import pandas as pd
import tempfile
import os
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from cobra.flux_analysis import (
    pfba, flux_variability_analysis,
    single_gene_deletion, single_reaction_deletion,
    double_gene_deletion, double_reaction_deletion,
    find_blocked_reactions
)
from cobra.io import load_model, save_json_model
import utils

COMPARTMENT_NAMES = {
    'c': 'Cytosol',
    'e': 'Extracellular',
    'p': 'Periplasm',
    'm': 'Mitochondria',
    'l': 'Lysosome',
    'n': 'Nucleus',
    'r': 'Endoplasmic Reticulum',
    'g': 'Golgi Apparatus',
    'v': 'Vacuole',
    'x': 'Peroxisome',
    'i': 'Inner Membrane',
    'o': 'Outer Membrane',
    'h': 'Hydrogenosome',
    'f': 'Flagellum',
    'w': 'Cell Wall',
    'u': 'Unknown/Unassigned',
}

def safe_set_bounds(model, rxn_id: str, lb: float = None, ub: float = None):
    """Safely set reaction bounds with validation."""
    if rxn_id not in model.reactions:
        raise ValueError(f"Reaction ID not in model: {rxn_id}")
    r = model.reactions.get_by_id(rxn_id)
    if lb is not None:
        r.lower_bound = lb
    if ub is not None:
        r.upper_bound = ub
    return r

def load_model_from_file(uploaded_file):
    """Load model from uploaded file with format detection."""
    try:
        with tempfile.NamedTemporaryFile(delete=False, suffix=uploaded_file.name) as tmp_file:
            tmp_file.write(uploaded_file.getvalue())
            tmp_path = tmp_file.name

        if uploaded_file.name.endswith('.xml'):
            model = cobra.io.read_sbml_model(tmp_path)
        elif uploaded_file.name.endswith('.json'):
            model = cobra.io.load_json_model(tmp_path)
        elif uploaded_file.name.endswith('.mat'):
            model = cobra.io.load_matlab_model(tmp_path)
        else:
            st.error("Unsupported file format. Please use SBML (.xml), JSON (.json), or MATLAB (.mat)")
            return None

        os.unlink(tmp_path)
        st.success("Model loaded successfully!")
        return model

    except Exception as e:
        st.error(f"Error loading model: {str(e)}")
        return None

def display_model_summary(model, unique_key_suffix=""):
    """Display comprehensive model summary."""
    col1, col2, col3 = st.columns(3)
    col1.metric("Reactions", len(model.reactions))
    col2.metric("Metabolites", len(model.metabolites))
    col3.metric("Genes", len(model.genes))

    objective = model.objective
    if isinstance(objective, str):
        st.write(f"**Objective Function:** {objective}")
    elif hasattr(objective, 'expression'):
        st.write(f"**Objective Function:** {objective.expression}")
    else:
        st.write(f"**Objective Function:** Not explicitly set or complex.")

    with st.expander("Reactions Overview"):
        reaction_data = []
        for rxn in model.reactions[:10]:
            reaction_data.append({
                "ID": rxn.id,
                "Name": rxn.name if rxn.name else "N/A",
                "Equation": rxn.build_reaction_string(),
                "Bounds": f"{rxn.lower_bound} to {rxn.upper_bound}"
            })
        st.dataframe(pd.DataFrame(reaction_data), hide_index=True)

    with st.expander("Metabolites Overview"):
        all_compartments = sorted(list(set(m.compartment for m in model.metabolites if m.compartment)))
        display_compartment_options = ["All"] + [COMPARTMENT_NAMES.get(c, c) for c in all_compartments]
        
        selected_compartment_display = st.selectbox(
            "View metabolites by compartment",
            display_compartment_options,
            key=f"met_compartment_select_{model.id}_{unique_key_suffix}"
        )
        
        selected_compartment_filter = selected_compartment_display
        if selected_compartment_display != "All":
            found_original_id = False
            for comp_id, comp_name in COMPARTMENT_NAMES.items():
                if comp_name == selected_compartment_display:
                    selected_compartment_filter = comp_id
                    found_original_id = True
                    break
            if not found_original_id:
                selected_compartment_filter = selected_compartment_display

        met_data = []
        for met in model.metabolites:
            if selected_compartment_filter == "All" or met.compartment == selected_compartment_filter:
                met_data.append({
                    "ID": met.id,
                    "Name": met.name if met.name else "N/A",
                    "Compartment": COMPARTMENT_NAMES.get(met.compartment, met.compartment) if met.compartment else "N/A"
                })
        st.dataframe(pd.DataFrame(met_data).head(100), hide_index=True)

    with st.expander("Genes Overview"):
        gene_data = []
        for gene in model.genes[:10]:
            gene_data.append({
                "ID": gene.id,
                "Name": gene.name if gene.name else "N/A",
                "Reactions": len(gene.reactions)
            })
        st.dataframe(pd.DataFrame(gene_data), hide_index=True)

def get_reaction_options(model):
    """Get formatted reaction options for dropdowns."""
    return [f"{rxn.name} ({rxn.id})" if rxn.name else rxn.id for rxn in model.reactions]

def get_metabolite_options(model):
    """Get formatted metabolite options for dropdowns."""
    return [f"{met.name} ({met.id})" if met.name else met.id for met in model.metabolites]

def get_gene_options(model):
    """Get formatted gene options for dropdowns."""
    return [f"{gene.name} ({gene.id})" if gene.name else gene.id for gene in model.genes]

def validate_genes(model, gene_ids):
    """Validate gene IDs against model."""
    valid_ids = []
    for gid in gene_ids:
        try:
            gene_obj = model.genes.get_by_id(gid)
            valid_ids.append(gene_obj.id)
        except KeyError:
            st.warning(f"Gene {gid} not found in model. Skipping...")
    return valid_ids

def perform_fba(model, objective_reaction_id=None):
    """Perform Flux Balance Analysis."""
    try:
        temp_model = model.copy()

        if objective_reaction_id:
            try:
                resolved_obj_id = utils.resolve_reaction_id(objective_reaction_id, temp_model)
                temp_model.objective = resolved_obj_id
            except ValueError as e:
                st.warning(f"Objective reaction '{objective_reaction_id}' not found or not specified. Using model's current objective: {temp_model.objective.expression}")
        else:
            st.warning(f"No objective reaction specified. Using model's current objective: {temp_model.objective.expression}")
            
        solution = temp_model.optimize()
        
        # FIX: Store using current_model_name instead of model.id
        if st.session_state.current_model_name:
            st.session_state.fba_solutions[st.session_state.current_model_name] = solution
        else:
            st.session_state.fba_solutions[model.id] = solution

        if solution.status == "optimal":
            st.success(f"Optimization successful! Objective value: {solution.objective_value:.4f}")

            fluxes = solution.fluxes
            nonzero_fluxes = fluxes[abs(fluxes) > 1e-6].sort_values(key=abs, ascending=False)

            st.subheader("Top Active Fluxes")
            flux_df_display = nonzero_fluxes.reset_index()
            flux_df_display.columns = ['Reaction ID', 'Flux']
            flux_df_display['Reaction ID'] = flux_df_display['Reaction ID'].astype(str)
            flux_df_display['Reaction Name'] = flux_df_display['Reaction ID'].apply(
                lambda x: temp_model.reactions.get_by_id(x).name if x in temp_model.reactions else "N/A"
            )
            st.dataframe(flux_df_display[['Reaction ID', 'Reaction Name', 'Flux']].head(10), hide_index=True)

            fig = px.histogram(
                fluxes,
                nbins=50,
                title="Flux Distribution",
                labels={'value': 'Flux Value'}
            )
            fig.update_layout(
                xaxis_title="Flux Value",
                yaxis_title="Count",
                template='plotly_white'
            )
            st.plotly_chart(fig, use_container_width=True)

            csv = flux_df_display.to_csv(index=False)
            st.download_button(
                label="Download Full Flux Distribution",
                data=csv,
                file_name=f'{model.id}_fba_fluxes.csv',
                mime='text/csv'
            )
        else:
            st.error(f"Optimization failed: {solution.status}")

        return solution

    except Exception as e:
        st.error(f"FBA failed: {str(e)}")
        return None

def plot_growth_impact(results_df, title, id_col_name="Gene ID"):
    """Plot growth impact for deletion analyses."""
    if results_df.empty:
        return

    results_df = results_df.copy()
    results_df.sort_values("Growth Rate", ascending=False, inplace=True)

    colors = ['red' if rate < 1e-6 else 'green' for rate in results_df["Growth Rate"]]

    plot_x_axis_col = None
    if "Gene Display" in results_df.columns:
        plot_x_axis_col = "Gene Display"
    elif "Reaction Display" in results_df.columns:
        plot_x_axis_col = "Reaction Display"
    elif id_col_name in results_df.columns:
        plot_x_axis_col = id_col_name
    elif f"{id_col_name.replace(' ID', '')} Display" in results_df.columns:
        plot_x_axis_col = f"{id_col_name.replace(' ID', '')} Display"
    else:
        st.warning(f"Could not find suitable display column for plotting. Using '{id_col_name}'.")
        plot_x_axis_col = id_col_name

    plot_x_axis = results_df[plot_x_axis_col]

    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            x=plot_x_axis,
            y=results_df["Growth Rate"],
            marker_color=colors,
            text=results_df["Growth Rate"].round(6),
            textposition='auto'
        )
    )

    fig.update_layout(
        title=title,
        xaxis_title=plot_x_axis_col.replace(' Display', '').replace(' ID', ''),
        yaxis_title="Growth Rate",
        template='plotly_white',
        height=500,
        hovermode='x'
    )
    fig.update_xaxes(tickangle=45)
    fig.add_hline(y=1e-6, line_dash="dash", line_color="gray")

    st.plotly_chart(fig, use_container_width=True)
