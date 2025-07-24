import streamlit as st
import cobra
import pandas as pd
import tempfile
import os
import re
import numpy as np
from cobra.flux_analysis import (
    pfba, flux_variability_analysis,
    single_gene_deletion, single_reaction_deletion,
    double_gene_deletion, double_reaction_deletion,
    find_blocked_reactions
)
from cobra.io import load_model, save_json_model
import plotly.express as px
import plotly.graph_objects as go
import json
import itertools
import traceback

st.set_page_config(
    page_title="COBRApy GUI",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

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

def load_model_from_file(uploaded_file):
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
    return [f"{rxn.name} ({rxn.id})" if rxn.name else rxn.id for rxn in model.reactions]

def get_metabolite_options(model):
    return [f"{met.name} ({met.id})" if met.name else met.id for met in model.metabolites]

def get_gene_options(model):
    return [f"{gene.name} ({gene.id})" if gene.name else gene.id for gene in model.genes]

def extract_id(option_string):
    match = re.search(r'\((.*?)\)', option_string)
    if match:
        return match.group(1)
    return option_string

def validate_genes(model, gene_ids):
    valid_ids = []
    for gid in gene_ids:
        try:
            gene_obj = model.genes.get_by_id(gid)
            valid_ids.append(gene_obj.id)
        except KeyError:
            st.warning(f"Gene {gid} not found in model. Skipping...")
    return valid_ids

def perform_fba(model, objective_reaction_id=None):
    try:
        temp_model = model.copy()

        if objective_reaction_id and objective_reaction_id in temp_model.reactions:
            temp_model.objective = objective_reaction_id
        else:
            st.warning(f"Objective reaction '{objective_reaction_id}' not found or not specified. Using model's current objective: {temp_model.objective.expression}")
            
        solution = temp_model.optimize()
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


def single_model_analysis_tabs():
    model_names = list(st.session_state.models.keys())
    
    if model_names:
        st.session_state.current_model_name = st.sidebar.selectbox(
            "Select Active Model",
            model_names,
            index=model_names.index(st.session_state.current_model_name) if st.session_state.current_model_name in model_names else 0,
            key="single_model_selector"
        )
    else:
        st.session_state.current_model_name = None
        st.warning("Please load a model first in the 'Model Loading' tab.")

    if st.session_state.current_model_name:
        model = st.session_state.models[st.session_state.current_model_name]
    else:
        tabs = ["Model Loading"]
        selected_tab = st.sidebar.selectbox("Navigation", tabs, key="single_model_nav_no_model")
        if selected_tab == "Model Loading":
            render_model_loading_section_single_mode()
        return

    tabs = [
        "Model Loading", "Reaction Details", "Flux Analysis", "Genetic Manipulations",
        "Model Diagnostics", "Export & Visualization"
    ]
    selected_tab = st.sidebar.selectbox("Navigation", tabs, key="single_model_nav")
    
    if selected_tab == "Model Loading":
        render_model_loading_section_single_mode()
            
    elif selected_tab == "Reaction Details":
        st.header("🔬 Reaction Details")
        st.info("Select a reaction to view its equation, involved metabolites, and associated genes.")

        reaction_options = get_reaction_options(model)
        selected_reaction_display = st.selectbox(
            "Select a reaction",
            reaction_options,
            key="reaction_details_select"
        )

        if selected_reaction_display:
            rxn_id = extract_id(selected_reaction_display)
            try:
                reaction = model.reactions.get_by_id(rxn_id)

                st.subheader(f"Details for: {reaction.name} ({reaction.id})")
                st.write(f"**Equation:** `{reaction.build_reaction_string()}`")
                st.write(f"**Bounds:** `{reaction.lower_bound} <= Flux <= {reaction.upper_bound}`")
                st.write(f"**Subsystem:** `{reaction.subsystem if reaction.subsystem else 'N/A'}`")
                st.write(f"**Gene-Protein-Reaction (GPR) Rule:** `{reaction.gene_reaction_rule if reaction.gene_reaction_rule else 'N/A'}`")

                st.markdown("#### Involved Metabolites")
                if reaction.metabolites:
                    metabolite_data = []
                    for met, coeff in reaction.metabolites.items():
                        metabolite_data.append({
                            "Metabolite ID": met.id,
                            "Metabolite Name": met.name if met.name else "N/A",
                            "Coefficient": coeff
                        })
                    st.dataframe(pd.DataFrame(metabolite_data), hide_index=True)
                else:
                    st.info("No metabolites found for this reaction.")

                st.markdown("#### Associated Genes")
                if reaction.genes:
                    gene_data = []
                    for gene in reaction.genes:
                        associated_reactions = [f"{r.name} ({r.id})" if r.name else r.id for r in gene.reactions]
                        gene_data.append({
                            "Gene ID": str(gene.id),
                            "Gene Name": gene.name if gene.name else "N/A",
                            "Gene Display": f"{gene.name} ({gene.id})" if gene.name else gene.id,
                            "Associated Reactions": ", ".join(associated_reactions)
                        })
                    st.dataframe(pd.DataFrame(gene_data), hide_index=True)
                else:
                    st.info("No genes associated with this reaction.")

            except KeyError:
                st.error(f"Reaction '{rxn_id}' not found in the model.")
            except Exception as e:
                st.error(f"Error retrieving reaction details: {str(e)}")
                st.exception(e)
                
    elif selected_tab == "Flux Analysis":
        st.header("📊 Flux Analysis")
        
        reaction_options = get_reaction_options(model)

        with st.expander("Flux Balance Analysis (FBA)"):
            st.info("FBA calculates the optimal flux distribution to maximize/minimize the objective function.")

            current_objective_id = None
            for rxn in model.reactions:
                if rxn.objective_coefficient != 0:
                    current_objective_id = rxn.id
                    break

            current_objective_display = None
            if current_objective_id:
                current_objective_rxn = model.reactions.get_by_id(current_objective_id)
                current_objective_display = f"{current_objective_rxn.name} ({current_objective_rxn.id})" if current_objective_rxn.name else current_objective_rxn.id
            else:
                if reaction_options:
                    current_objective_display = reaction_options[0]

            default_objective_index = 0
            if current_objective_display and current_objective_display in reaction_options:
                default_objective_index = reaction_options.index(current_objective_display)

            objective_display = st.selectbox(
                "Objective Reaction",
                reaction_options,
                index=default_objective_index,
                key="fba_objective_select",
                help="Select the reaction to optimize"
            )
            objective_rxn_id = extract_id(objective_display)

            st.subheader("Adjust Reaction Bounds")
            rxn_display = st.selectbox("Select reaction to modify", reaction_options, key="set_bounds_rxn_select_fba")
            rxn_id_for_bounds = extract_id(rxn_display)
            selected_rxn = model.reactions.get_by_id(rxn_id_for_bounds)

            col_lb, col_ub = st.columns(2)
            with col_lb:
                new_lb = st.number_input(
                    "Lower Bound",
                    value=float(selected_rxn.lower_bound),
                    format="%f",
                    key=f"lb_input_{rxn_id_for_bounds}_fba"
                )
            with col_ub:
                new_ub = st.number_input(
                    "Upper Bound",
                    value=float(selected_rxn.upper_bound),
                    format="%f",
                    key=f"ub_input_{rxn_id_for_bounds}_fba"
                )

            if st.button("Apply Bounds", key=f"apply_bounds_btn_{rxn_id_for_bounds}_fba"):
                selected_rxn.bounds = (new_lb, new_ub)
                st.success(f"Updated bounds for {rxn_id_for_bounds}: [{new_lb}, {new_ub}]")

            if st.button("Run FBA", type="primary", key="run_fba_btn_main"):
                with st.spinner("Performing Flux Balance Analysis..."):
                    perform_fba(model, objective_rxn_id)

        with st.expander("Parsimonious FBA (pFBA)"):
            st.info("pFBA finds the minimal total flux solution that achieves optimal growth.")

            if st.button("Run pFBA", key="run_pfba_btn"):
                with st.spinner("Performing Parsimonious FBA..."):
                    try:
                        p_solution = pfba(model)
                        st.session_state.fba_solutions[model.id] = p_solution

                        st.success(f"pFBA completed! Total flux: {p_solution.fluxes.abs().sum():.4f}")

                        if model.id in st.session_state.fba_solutions:
                            fba_solution = st.session_state.fba_solutions[model.id]
                            fba_fluxes = fba_solution.fluxes
                            p_fluxes = p_solution.fluxes

                            nonzero_fba = sum(abs(fba_fluxes) > 1e-6)
                            nonzero_pfba = sum(abs(p_fluxes) > 1e-6)

                            st.write(f"**FBA Active Reactions:** {nonzero_fba}")
                            st.write(f"**pFBA Active Reactions:** {nonzero_pfba}")

                            st.subheader("Top pFBA Fluxes")
                            nonzero_pfluxes = p_fluxes[abs(p_fluxes) > 1e-6].sort_values(key=abs, ascending=False)
                            
                            flux_df_display = nonzero_pfluxes.reset_index()
                            flux_df_display.columns = ['Reaction ID', 'Flux']
                            flux_df_display['Reaction ID'] = flux_df_display['Reaction ID'].astype(str)
                            flux_df_display['Reaction Name'] = flux_df_display['Reaction ID'].apply(
                                lambda x: model.reactions.get_by_id(x).name if x in model.reactions else "N/A"
                            )
                            st.dataframe(flux_df_display[['Reaction ID', 'Reaction Name', 'Flux']].head(10), hide_index=True)

                            fig = go.Figure()
                            fig.add_trace(
                                go.Scatter(
                                    x=fba_fluxes,
                                    y=p_fluxes,
                                    mode='markers',
                                    marker=dict(
                                        size=8,
                                        color='rgba(55, 126, 184, 0.7)',
                                        line=dict(width=1, color='DarkSlateGrey')
                                )
                            ))

                            max_val = max(fba_fluxes.max(), p_fluxes.max())
                            min_val = min(fba_fluxes.min(), p_fluxes.min())
                            fig.add_trace(
                                go.Scatter(
                                    x=[min_val, max_val],
                                    y=[min_val, max_val],
                                    mode='lines',
                                    line=dict(color='red', dash='dash'),
                                    name='y=x'
                                )
                            )

                            fig.update_layout(
                                title="FBA vs pFBA Flux Comparison",
                                xaxis_title="FBA Flux",
                                yaxis_title="pFBA Flux",
                                template='plotly_white',
                                showlegend=False
                            )
                            st.plotly_chart(fig, use_container_width=True)
                    except Exception as e:
                        st.error(f"pFBA failed: {str(e)}")

        with st.expander("Flux Variability Analysis (FVA)"):
            st.info("FVA calculates the minimum and maximum possible flux for each reaction.")

            reactions = st.multiselect(
                "Select reactions for FVA (leave empty for all)",
                reaction_options,
                help="Analyze specific reactions or all reactions in model",
                key="fva_rxn_select"
            )
            reaction_ids = [extract_id(r) for r in reactions] if reactions else []

            fraction = st.slider(
                "Fraction of Optimum",
                min_value=0.5, max_value=1.0, value=0.95, step=0.01,
                help="Percentage of optimal growth rate to consider",
                key="fva_fraction_slider"
            )

            if st.button("Run FVA", key="run_fva_btn"):
                if not reaction_ids:
                    reaction_list = model.reactions
                else:
                    reaction_list = [model.reactions.get_by_id(rxn_id) for rxn_id in reaction_ids]

                with st.spinner("Performing Flux Variability Analysis (this may take a while)..."):
                    try:
                        fva_results = flux_variability_analysis(
                            model,
                            reaction_list=reaction_list,
                            fraction_of_optimum=fraction
                        )

                        st.success("FVA completed successfully!")
                        
                        fva_results_display = fva_results.reset_index()
                        fva_results_display.columns = ['Reaction ID', 'minimum', 'maximum']
                        fva_results_display['Reaction ID'] = fva_results_display['Reaction ID'].astype(str)
                        fva_results_display['Reaction Name'] = fva_results_display['Reaction ID'].apply(
                            lambda x: model.reactions.get_by_id(x).name if x in model.reactions else "N/A"
                        )
                        st.dataframe(fva_results_display[['Reaction ID', 'Reaction Name', 'minimum', 'maximum']], hide_index=True)

                        fig = go.Figure()

                        fig.add_trace(
                            go.Bar(
                                x=fva_results_display['Reaction ID'],
                                y=fva_results_display['minimum'],
                                name='Minimum Flux',
                                marker_color='#636EFA'
                            )
                        )

                        fig.add_trace(
                            go.Bar(
                                x=fva_results_display['Reaction ID'],
                                y=fva_results_display['maximum'] - fva_results_display['minimum'],
                                name='Flux Range',
                                base=fva_results_display['minimum'],
                                marker_color='#EF553B'
                            )
                        )

                        fig.update_layout(
                            title="Flux Variability Analysis",
                            xaxis_title="Reaction",
                            yaxis_title="Flux Range",
                            barmode='stack',
                            template='plotly_white',
                            height=600,
                            legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1)
                        )
                        fig.update_xaxes(tickangle=45)
                        st.plotly_chart(fig, use_container_width=True)

                        csv = fva_results_display.to_csv(index=False)
                        st.download_button(
                            label="Download FVA Results",
                            data=csv,
                            file_name=f'{model.id}_fva_results.csv',
                            mime='text/csv'
                        )
                    except Exception as e:
                        st.error(f"FVA failed: {str(e)}")

    elif selected_tab == "Genetic Manipulations":
        st.header("🧬 Genetic Manipulations")
        
        gene_options = get_gene_options(model)
        reaction_options = get_reaction_options(model)

        tab_gene_del, tab_rxn_del, tab_essential_genes, tab_set_sim = st.tabs([
            "Gene Deletion Analysis", "Reaction Deletion Analysis", "Find Essential Genes", "Gene/Reaction Set Simulation"
        ])

        with tab_gene_del:
            st.subheader("Gene Deletion Analysis")
            deletion_type = st.radio("Deletion Type", ["Single", "Double"], index=0, key="gene_del_type_radio")

            selected_genes = st.multiselect(
                "Select genes to delete",
                gene_options,
                key="single_double_gene_select"
            )
            gene_ids_raw = [extract_id(g) for g in selected_genes]

            if st.button("Simulate Gene Deletions", key="run_gene_deletion_btn"):
                if not gene_ids_raw:
                    st.warning("Please select at least one gene.")
                else:
                    gene_ids = validate_genes(model, gene_ids_raw)
                    if not gene_ids:
                        st.error("No valid genes selected for deletion that exist in the model.")
                    else:
                        with st.spinner("Analyzing gene deletions..."):
                            try:
                                if deletion_type == "Single":
                                    deletion_results = single_gene_deletion(model, gene_ids)

                                    results_df = pd.DataFrame({
                                        "Gene ID": [str(g_id) for g_id in deletion_results.index.get_level_values(0)],
                                        "Growth Rate": deletion_results['growth']
                                    })
                                    results_df["Gene Name"] = results_df["Gene ID"].apply(
                                        lambda g_id: model.genes.get_by_id(g_id).name if g_id in model.genes else "N/A"
                                    )
                                    results_df["Gene Display"] = results_df.apply(
                                        lambda row: f"{row['Gene Name']} ({row['Gene ID']})" if row['Gene Name'] != "N/A" else row['Gene ID'], axis=1
                                    )
                                    results_df["Lethal"] = results_df["Growth Rate"] < 1e-6
                                    st.dataframe(results_df[["Gene Display", "Growth Rate", "Lethal", "Gene ID", "Gene Name"]], hide_index=True)
                                    plot_growth_impact(results_df, "Single Gene Deletion Impact", id_col_name="Gene Display")

                                else:
                                    if len(gene_ids) < 2:
                                        st.error("Please select at least 2 genes for double deletion")
                                        return
                                    deletion_results = double_gene_deletion(model, gene_ids)

                                    results_df = pd.DataFrame({
                                        "Gene 1 ID": [str(g_id) for g_id in deletion_results.index.get_level_values(0)],
                                        "Gene 2 ID": [str(g_id) for g_id in deletion_results.index.get_level_values(1)],
                                        "Growth Rate": deletion_results['growth']
                                    })
                                    results_df["Gene 1 Name"] = results_df["Gene 1 ID"].apply(
                                        lambda g_id: model.genes.get_by_id(g_id).name if g_id in model.genes else "N/A"
                                    )
                                    results_df["Gene 1 Display"] = results_df.apply(
                                        lambda row: f"{row['Gene 1 Name']} ({row['Gene 1 ID']})" if row['Gene 1 Name'] != "N/A" else row['Gene 1 ID'], axis=1
                                    )
                                    results_df["Gene 2 Name"] = results_df["Gene 2 ID"].apply(
                                        lambda g_id: model.genes.get_by_id(g_id).name if g_id in model.genes else "N/A"
                                    )
                                    results_df["Gene 2 Display"] = results_df.apply(
                                        lambda row: f"{row['Gene 2 Name']} ({row['Gene 2 ID']})" if row['Gene 2 Name'] != "N/A" else row['Gene 2 ID'], axis=1
                                    )
                                    results_df["Lethal"] = results_df["Growth Rate"] < 1e-6
                                    st.dataframe(results_df[["Gene 1 Display", "Gene 2 Display", "Growth Rate", "Lethal", "Gene 1 ID", "Gene 1 Name", "Gene 2 ID", "Gene 2 Name"]], hide_index=True)
                                
                                lethal_count = sum(results_df["Lethal"])
                                st.info(f"**{lethal_count} lethal deletion(s) found**")

                                csv = results_df.to_csv(index=False)
                                st.download_button(
                                    label="Download Gene Deletion Results",
                                    data=csv,
                                    file_name=f'{st.session_state.current_model_name}_gene_deletion_results.csv',
                                    mime='text/csv'
                                )
                            except Exception as e:
                                st.error(f"Gene deletion analysis failed: {str(e)}")
                                st.exception(e)

        with tab_rxn_del:
            st.subheader("Reaction Deletion Analysis")
            deletion_type_rxn = st.radio("Deletion Type", ["Single", "Double"], index=0, key="rxn_del_type_radio")

            selected_reactions = st.multiselect(
                "Select reactions to delete",
                reaction_options,
                key="single_double_rxn_select"
            )
            reaction_ids = [extract_id(r) for r in selected_reactions]

            if st.button("Simulate Reaction Deletions", key="run_rxn_deletion_btn"):
                if not reaction_ids:
                    st.warning("Please select at least one reaction")
                else:
                    with st.spinner("Analyzing reaction deletions..."):
                        try:
                            if deletion_type_rxn == "Single":
                                deletion_results = single_reaction_deletion(model, reaction_ids)
                                results_df = pd.DataFrame({
                                    "Reaction ID": deletion_results.index.get_level_values(0).astype(str),
                                    "Growth Rate": deletion_results['growth']
                                })
                                results_df["Reaction Name"] = results_df["Reaction ID"].apply(
                                    lambda r_id: model.reactions.get_by_id(r_id).name if r_id in model.reactions else "N/A"
                                )
                                results_df["Reaction Display"] = results_df.apply(
                                    lambda row: f"{row['Reaction Name']} ({row['Reaction ID']})" if row['Reaction Name'] != "N/A" else row['Reaction ID'], axis=1
                                )
                                results_df["Lethal"] = results_df["Growth Rate"] < 1e-6
                                st.dataframe(results_df[["Reaction Display", "Growth Rate", "Lethal", "Reaction ID", "Reaction Name"]], hide_index=True)
                                plot_growth_impact(results_df, "Single Reaction Deletion Impact", id_col_name="Reaction Display")
                            else:
                                if len(reaction_ids) < 2:
                                    st.error("Please select at least 2 reactions for double deletion")
                                    return
                                deletion_results = double_reaction_deletion(model, reaction_ids)

                                results_df = pd.DataFrame({
                                    "Reaction 1 ID": deletion_results.index.get_level_values(0).astype(str),
                                    "Reaction 2 ID": deletion_results.index.get_level_values(1).astype(str),
                                    "Growth Rate": deletion_results['growth']
                                })
                                results_df["Reaction 1 Name"] = results_df["Reaction 1 ID"].apply(
                                    lambda r_id: model.reactions.get_by_id(r_id).name if r_id in model.reactions else "N/A"
                                )
                                results_df["Reaction 1 Display"] = results_df.apply(
                                    lambda row: f"{row['Reaction 1 Name']} ({row['Reaction 1 ID']})" if row['Reaction 1 Name'] != "N/A" else row['Reaction 1 ID'], axis=1
                                )
                                results_df["Reaction 2 Name"] = results_df["Reaction 2 ID"].apply(
                                    lambda r_id: model.reactions.get_by_id(r_id).name if r_id in model.reactions else "N/A"
                                )
                                results_df["Reaction 2 Display"] = results_df.apply(
                                    lambda row: f"{row['Reaction 2 Name']} ({row['Reaction 2 ID']})" if row['Reaction 2 Name'] != "N/A" else row['Reaction 2 ID'], axis=1
                                )
                                results_df["Lethal"] = results_df["Growth Rate"] < 1e-6
                                st.dataframe(results_df[["Reaction 1 Display", "Reaction 2 Display", "Growth Rate", "Lethal", "Reaction 1 ID", "Reaction 1 Name", "Reaction 2 ID", "Reaction 2 Name"]], hide_index=True)

                            lethal_count_rxn = sum(results_df["Lethal"])
                            st.info(f"**{lethal_count_rxn} lethal deletion(s) found**")

                            csv = results_df.to_csv(index=False)
                            st.download_button(
                                label="Download Reaction Deletion Results",
                                data=csv,
                                file_name=f'{st.session_state.current_model_name}_reaction_deletion_results.csv',
                                mime='text/csv'
                            )
                        except Exception as e:
                            st.error(f"Reaction deletion analysis failed: {str(e)}")
                            st.exception(e)

        with tab_essential_genes:
            st.subheader("Find Essential Genes")
            st.info("Identifies genes whose deletion leads to a significant drop in growth rate (lethal).")

            if st.button("Run Essential Gene Analysis", key="run_essential_gene_btn"):
                with st.spinner("Finding essential genes..."):
                    try:
                        gene_objects = [model.genes.get_by_id(gid) for gid in model.genes.list_attr("id")]
                        essential_results = single_gene_deletion(model, gene_list=gene_objects)
                        gene_ids = [gene.id for gene in gene_objects]
                        essential_results.index = gene_ids
                        essential_results.index.name = "Gene ID"
                        gene_ids_str = [str(g) for g in gene_ids]
                        essential_df = pd.DataFrame({
                            "Gene ID": gene_ids_str,
                            "Growth Rate": essential_results['growth'].values
                        })
                        essential_df["Gene Name"] = essential_df["Gene ID"].apply(
                            lambda g_id: model.genes.get_by_id(g_id).name if g_id in model.genes else "N/A"
                        )
                        essential_df["Gene Display"] = essential_df.apply(
                            lambda row: f"{row['Gene Name']} ({row['Gene ID']})" if row['Gene Name'] != "N/A" else row["Gene ID"],
                            axis=1
                        )
                        essential_df["Status"] = essential_df["Growth Rate"].apply(
                            lambda x: "Lethal" if x < 1e-6 else "Non-Lethal"
                        )
                        lethal_genes = essential_df[essential_df["Status"] == "Lethal"]
                        non_lethal_genes = essential_df[essential_df["Status"] == "Non-Lethal"]

                        st.success("Essential gene analysis completed!")
                        st.info(f"Found **{len(lethal_genes)} lethal genes** and **{len(non_lethal_genes)} non-lethal genes**.")

                        st.subheader("Lethal Genes (Growth Rate < 1e-6)")
                        if not lethal_genes.empty:
                            st.dataframe(lethal_genes[["Gene Display", "Growth Rate", "Status", "Gene ID", "Gene Name"]], hide_index=True)
                        else:
                            st.info("No lethal genes found under current conditions.")

                        st.subheader("Non-Lethal Genes")
                        if not non_lethal_genes.empty:
                            st.dataframe(non_lethal_genes[["Gene Display", "Growth Rate", "Status", "Gene ID", "Gene Name"]], hide_index=True)
                        else:
                            st.info("All genes are lethal or no genes found.")

                        csv = essential_df.to_csv(index=False)
                        st.download_button(
                            label="Download Essential Gene Results",
                            data=csv,
                            file_name=f'{st.session_state.current_model_name}_essential_gene_results.csv',
                            mime='text/csv'
                        )

                    except Exception as e:
                        st.error(f"Essential gene analysis failed: {str(e)}")
                        st.exception(e)

        with tab_set_sim:
            st.subheader("Gene Set Deletion Simulation")
            st.info("Simulate the deletion of multiple gene sets.")

            num_gene_sets = st.number_input("Number of Gene Sets to Define", min_value=1, value=1, key="num_gene_sets")
            gene_sets = []
            for i in range(num_gene_sets):
                set_name = st.text_input(f"Gene Set {i+1} Name", value=f"GeneSet_{i+1}", key=f"gene_set_name_{i}")
                selected_genes_in_set = st.multiselect(
                    f"Select genes for {set_name}",
                    gene_options,
                    key=f"gene_set_select_{i}"
                )
                gene_sets.append((set_name, [extract_id(g) for g in selected_genes_in_set]))

            if st.button("Run Gene Set Deletions", key="run_gene_set_deletion_btn"):
                if not gene_sets:
                    st.warning("Please define at least one gene set.")
                else:
                    all_results = []
                    for set_name, gene_ids in gene_sets:
                        if not gene_ids:
                            st.warning(f"Gene Set '{set_name}' is empty. Skipping.")
                            continue
                        gene_ids_valid = validate_genes(model, gene_ids)
                        if not gene_ids_valid:
                            st.error(f"No valid genes in set '{set_name}'. Skipping.")
                            continue

                        st.write(f"**Simulating deletion for Gene Set: {set_name}**")
                        with st.spinner(f"Analyzing gene set '{set_name}'..."):
                            try:
                                with model:
                                    for gene_id in gene_ids_valid:
                                        model.genes.get_by_id(gene_id).knock_out()
                                    
                                    solution = model.optimize()
                                    growth_rate = solution.objective_value if solution.status == "optimal" else 0.0
                                    display_genes_names = []
                                    display_genes_ids = []
                                    display_genes_combined = []
                                    for g_id in gene_ids_valid:
                                        gene_obj = model.genes.get_by_id(g_id)
                                        gene_name = gene_obj.name if gene_obj.name else "N/A"
                                        display_genes_names.append(gene_name)
                                        display_genes_ids.append(str(g_id))
                                        display_genes_combined.append(f"{gene_name} ({g_id})" if gene_name != "N/A" else g_id)

                                    all_results.append({
                                        "Set Name": set_name,
                                        "Genes Deleted (Display)": ", ".join(display_genes_combined),
                                        "Growth Rate": growth_rate,
                                        "Status": "Lethal" if growth_rate < 1e-6 else "Non-Lethal",
                                        "Gene IDs Deleted": ", ".join(display_genes_ids),
                                        "Gene Names Deleted": ", ".join(display_genes_names)
                                    })
                            except Exception as e:
                                st.error(f"Error simulating gene set '{set_name}': {str(e)}")
                                all_results.append({
                                    "Set Name": set_name,
                                    "Genes Deleted (Display)": ", ".join(gene_ids_valid),
                                    "Growth Rate": "Error",
                                    "Status": "Error",
                                    "Gene IDs Deleted": ", ".join(gene_ids_valid),
                                    "Gene Names Deleted": "Error"
                                })
                    
                    if all_results:
                        results_df = pd.DataFrame(all_results)
                        st.subheader("Gene Set Deletion Results")
                        st.dataframe(results_df, hide_index=True)

                        csv = results_df.to_csv(index=False)
                        st.download_button(
                            label="Download Gene Set Deletion Results",
                            data=csv,
                            file_name=f'{st.session_state.current_model_name}_gene_set_deletion_results.csv',
                            mime='text/csv'
                        )

            st.subheader("Reaction Set Deletion Simulation")
            st.info("Simulate the deletion of multiple reaction sets (pathway blocks).")

            num_rxn_sets = st.number_input("Number of Reaction Sets to Define", min_value=1, value=1, key="num_rxn_sets")
            rxn_sets = []
            for i in range(num_rxn_sets):
                set_name = st.text_input(f"Reaction Set {i+1} Name", value=f"ReactionSet_{i+1}", key=f"rxn_set_name_{i}")
                selected_rxns_in_set = st.multiselect(
                    f"Select reactions for {set_name}",
                    reaction_options,
                    key=f"rxn_set_select_{i}"
                )
                rxn_sets.append((set_name, [extract_id(r) for r in selected_rxns_in_set]))

            if st.button("Run Reaction Set Deletions", key="run_rxn_set_deletion_btn"):
                if not rxn_sets:
                    st.warning("Please define at least one reaction set.")
                else:
                    all_results = []
                    for set_name, rxn_ids in rxn_sets:
                        if not rxn_ids:
                            st.warning(f"Reaction Set '{set_name}' is empty. Skipping.")
                            continue
                        
                        valid_rxn_ids = []
                        for rid in rxn_ids:
                            if rid in model.reactions:
                                valid_rxn_ids.append(rid)
                            else:
                                st.warning(f"Reaction {rid} not found in model for set '{set_name}'. Skipping.")
                        
                        if not valid_rxn_ids:
                            st.error(f"No valid reactions in set '{set_name}'. Skipping.")
                            continue

                        st.write(f"**Simulating deletion for Reaction Set: {set_name}**")
                        with st.spinner(f"Analyzing reaction set '{set_name}'..."):
                            try:
                                with model:
                                    for rxn_id in valid_rxn_ids:
                                        model.reactions.get_by_id(rxn_id).knock_out()
                                    
                                    solution = model.optimize()
                                    growth_rate = solution.objective_value if solution.status == "optimal" else 0.0

                                    display_rxns_names = []
                                    display_rxns_ids = []
                                    display_rxns_combined = []
                                    for r_id in valid_rxn_ids:
                                        rxn_obj = model.reactions.get_by_id(r_id)
                                        rxn_name = rxn_obj.name if rxn_obj.name else "N/A"
                                        display_rxns_names.append(rxn_name)
                                        display_rxns_ids.append(str(r_id))
                                        display_rxns_combined.append(f"{rxn_name} ({r_id})" if rxn_name != "N/A" else r_id)

                                    all_results.append({
                                        "Set Name": set_name,
                                        "Reactions Deleted (Display)": ", ".join(display_rxns_combined),
                                        "Growth Rate": growth_rate,
                                        "Status": "Lethal" if growth_rate < 1e-6 else "Non-Lethal",
                                        "Reaction IDs Deleted": ", ".join(display_rxns_ids),
                                        "Reaction Names Deleted": ", ".join(display_rxns_names)
                                    })
                            except Exception as e:
                                st.error(f"Error simulating reaction set '{set_name}': {str(e)}")
                                all_results.append({
                                    "Set Name": set_name,
                                    "Reactions Deleted (Display)": ", ".join(valid_rxn_ids),
                                    "Growth Rate": "Error",
                                    "Status": "Error",
                                    "Reaction IDs Deleted": ", ".join(valid_rxn_ids),
                                    "Reaction Names Deleted": "Error"
                                })
                    
                    if all_results:
                        results_df = pd.DataFrame(all_results)
                        st.subheader("Reaction Set Deletion Results")
                        st.dataframe(results_df, hide_index=True)

                        csv = results_df.to_csv(index=False)
                        st.download_button(
                            label="Download Reaction Set Deletion Results",
                            data=csv,
                            file_name=f'{st.session_state.current_model_name}_reaction_set_deletion_results.csv',
                            mime='text/csv'
                        )

    elif selected_tab == "Model Diagnostics":
        st.header("🔍 Model Diagnostics")
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.subheader("Blocked Reactions")
            st.info("Reactions that cannot carry any flux under any condition")
            
            if st.button("Find Blocked Reactions", key="find_blocked_rxns_btn"):
                with st.spinner("Identifying blocked reactions..."):
                    try:
                        blocked = find_blocked_reactions(model)
                        if len(blocked) > 0:
                            st.warning(f"Found {len(blocked)} blocked reactions")
                            
                            blocked_data = []
                            for rxn_id in blocked:
                                rxn = model.reactions.get_by_id(rxn_id)
                                blocked_data.append({
                                    "ID": str(rxn.id),
                                    "Name": rxn.name if rxn.name else "N/A",
                                    "Equation": rxn.build_reaction_string()
                                })
                            st.dataframe(pd.DataFrame(blocked_data), hide_index=True)
                            
                            if st.button("Remove Blocked Reactions", key="remove_blocked_rxns_btn"):
                                model.remove_reactions(blocked)
                                st.success(f"Removed {len(blocked)} blocked reactions")
                        else:
                            st.success("No blocked reactions found")
                    except Exception as e:
                        st.error(f"Error finding blocked reactions: {str(e)}")
        
        with col2:
            st.subheader("Model Consistency Check")
            st.info("Verify model structural and numerical consistency")
            
            if st.button("Check Model Consistency", key="check_model_consistency_btn"):
                with st.spinner("Checking model consistency..."):
                    try:
                        solution = model.optimize()
                        if solution.status == "optimal":
                            st.success("Model is consistent and solvable")
                        else:
                            st.error(f"Model inconsistency detected: {solution.status}")
                    except Exception as e:
                        st.error(f"Model consistency check failed: {str(e)}")

    elif selected_tab == "Export & Visualization":
        st.header("💾 Export & Visualization")
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.subheader("Export Model")
            export_format = st.selectbox(
                "Select export format",
                ["SBML", "JSON", "MATLAB"],
                key="export_format_select"
            )
            
            export_filename = st.text_input(
                "Filename", 
                value=f"{model.id}_export",
                help="Filename without extension",
                key="export_filename_input"
            )
            
            if st.button("Export Model", key="export_model_btn"):
                try:
                    with tempfile.TemporaryDirectory() as tmpdir:
                        filepath = os.path.join(tmpdir, export_filename)
                        
                        if export_format == "SBML":
                            filepath += ".xml"
                            cobra.io.write_sbml_model(model, filepath)
                        elif export_format == "JSON":
                            filepath += ".json"
                            save_json_model(model, filepath)
                        elif export_format == "MATLAB":
                            filepath += ".mat"
                            cobra.io.save_matlab_model(model, filepath)
                        
                        with open(filepath, "rb") as f:
                            st.download_button(
                                label=f"Download {export_format} File",
                                data=f,
                                file_name=os.path.basename(filepath),
                                key="download_exported_model_btn"
                            )
                except Exception as e:
                    st.error(f"Export failed: {str(e)}")
        
        with col2:
            st.subheader("Flux Visualization")
            
            if st.session_state.current_model_name in st.session_state.fba_solutions:
                st.success("Flux solution available for visualization")
                fba_solution = st.session_state.fba_solutions[st.session_state.current_model_name]
                flux_df = fba_solution.fluxes.reset_index()
                flux_df.columns = ['Reaction ID', 'Flux']
                
                flux_df['Reaction ID'] = flux_df['Reaction ID'].astype(str)
                flux_df['Reaction Name'] = flux_df['Reaction ID'].apply(
                    lambda x: model.reactions.get_by_id(x).name if x in model.reactions else "N/A"
                )
                
                csv = flux_df.to_csv(index=False)
                st.download_button(
                    label="Download Flux Data",
                    data=csv,
                    file_name=f'{model.id}_flux_data.csv',
                    mime='text/csv',
                    key="download_flux_data_btn"
                )
                
                top_fluxes = flux_df[abs(flux_df['Flux']) > 1e-6].sort_values(
                    by='Flux', key=abs, ascending=False
                ).head(10)
                
                if not top_fluxes.empty:
                    fig = px.bar(
                        top_fluxes,
                        x='Reaction Name',
                        y='Flux',
                        color='Flux',
                        color_continuous_scale='RdBu',
                        title="Top 10 Fluxes",
                        text='Flux'
                    )
                    fig.update_traces(texttemplate='%{text:.4f}', textposition='outside')
                    fig.update_layout(
                        xaxis_title="Reaction",
                        yaxis_title="Flux Value",
                        template='plotly_white',
                        height=500
                    )
                    fig.update_xaxes(tickangle=45)
                    st.plotly_chart(fig, use_container_width=True)
            else:
                st.warning("Run FBA first to generate flux solution")

def render_model_loading_section_single_mode():
    st.header("📤 Upload a metabolic model or create a new one.")

    col1, col2 = st.columns(2)
    with col1:
        model_source = st.radio("Model source:", ["Built-in", "Upload File", "URL", "Create Blank Model"], key="single_load_source_radio")

        if model_source == "Built-in":
            builtin_models = {
                "E. coli core": "e_coli_core",
                "Textbook E. coli": "textbook",
                "Salmonella core": "salmonella",
                "iJO1366 (E. coli)": "iJO1366",
                "Recon3D (Human)": "Recon3D"
            }
            selected_model_key = st.selectbox("Select model", list(builtin_models.keys()), key="builtin_model_select_single")
            model_label = st.text_input("Label for this model", value=selected_model_key, key="builtin_model_label_single")

            if st.button("Load Built-in Model", key="load_builtin_btn_single"):
                if model_label in st.session_state.models:
                    st.warning(f"Model with label '{model_label}' already exists. Please choose a different label.")
                else:
                    try:
                        new_model = load_model(builtin_models[selected_model_key])
                        st.session_state.models[model_label] = new_model
                        st.session_state.current_model_name = model_label
                        st.rerun()
                    except Exception as e:
                        st.error(f"Error loading model: {str(e)}")

        elif model_source == "Upload File":
            uploaded_file = st.file_uploader(
                "Upload model (SBML/JSON/MAT)",
                type=["xml", "json", "mat"],
                help="Supported formats: SBML (.xml), JSON (.json), MATLAB (.mat)",
                key="uploaded_file_uploader_single"
            )
            model_label = st.text_input("Label for this model", value=uploaded_file.name.split('.')[0] if uploaded_file else "uploaded_model", key="uploaded_model_label_single")

            if st.button("Load Uploaded Model", key="load_uploaded_btn_single") and uploaded_file:
                if model_label in st.session_state.models:
                    st.warning(f"Model with label '{model_label}' already exists. Please choose a different label.")
                else:
                    new_model = load_model_from_file(uploaded_file)
                    if new_model:
                        st.session_state.models[model_label] = new_model
                        st.session_state.current_model_name = model_label
                        st.rerun()

        elif model_source == "URL":
            url = st.text_input(
                "Model URL",
                placeholder="e.g., http://bigg.ucsd.edu/models/e_coli_core.xml",
                help="Enter direct download URL for model file",
                key="url_input_single"
            )
            model_label = st.text_input("Label for this model", value="url_model", key="url_model_label_single")

            if st.button("Load from URL", key="load_url_btn_single"):
                if model_label in st.session_state.models:
                    st.warning(f"Model with label '{model_label}' already exists. Please choose a different label.")
                else:
                    st.warning("URL loading requires additional implementation")

        elif model_source == "Create Blank Model":
            blank_model_label = st.text_input("Label for new blank model", value="new_model", key="blank_model_label_input_single")
            if st.button("Create Blank Model", key="create_blank_model_btn_single"):
                if blank_model_label in st.session_state.models:
                    st.warning(f"Model with label '{blank_model_label}' already exists. Please choose a different label.")
                else:
                    new_blank_model = cobra.Model(id=blank_model_label)
                    st.session_state.models[blank_model_label] = new_blank_model
                    st.session_state.current_model_name = blank_model_label
                    st.success(f"Blank model '{blank_model_label}' created. Go to 'Model Editing' to build it out.")
                    st.rerun()

    if st.session_state.current_model_name and st.session_state.current_model_name in st.session_state.models:
        st.subheader(f"Summary for Active Model: {st.session_state.current_model_name}")
        display_model_summary(st.session_state.models[st.session_state.current_model_name], unique_key_suffix="active")
    
def render_model_loading_section_multi_mode():
    st.header("📤 Upload a metabolic model or create a new one.")

    col1, col2 = st.columns(2)
    with col1:
        model_source = st.radio("Model source:", ["Built-in", "Upload File", "URL", "Create Blank Model"], key="multi_model_source_radio")

        if model_source == "Built-in":
            builtin_models = {
                "E. coli core": "e_coli_core",
                "Textbook E. coli": "textbook",
                "Salmonella core": "salmonella",
                "iJO1366 (E. coli)": "iJO1366",
                "Recon3D (Human)": "Recon3D"
            }
            selected_model_key = st.selectbox("Select model", list(builtin_models.keys()), key="multi_builtin_model_select")
            model_label = st.text_input("Label for this model", value=selected_model_key, key="multi_builtin_model_label")

            if st.button("Load Built-in Model", key="multi_load_builtin_btn"):
                if model_label in st.session_state.models:
                    st.warning(f"Model with label '{model_label}' already exists. Please choose a different label.")
                else:
                    try:
                        new_model = load_model(builtin_models[selected_model_key])
                        st.session_state.models[model_label] = new_model
                        st.success(f"Loaded built-in model: {selected_model_key} as '{model_label}'")
                        st.rerun()
                    except Exception as e:
                        st.error(f"Error loading model: {str(e)}")

        elif model_source == "Upload File":
            uploaded_file = st.file_uploader(
                "Upload model (SBML/JSON/MAT)",
                type=["xml", "json", "mat"],
                help="Supported formats: SBML (.xml), JSON (.json), MATLAB (.mat)",
                key="multi_uploaded_file_uploader"
            )
            model_label = st.text_input("Label for this model", value=uploaded_file.name.split('.')[0] if uploaded_file else "uploaded_model", key="multi_uploaded_model_label")

            if st.button("Load Uploaded Model", key="multi_load_uploaded_btn") and uploaded_file:
                if model_label in st.session_state.models:
                    st.warning(f"Model with label '{model_label}' already exists. Please choose a different label.")
                else:
                    new_model = load_model_from_file(uploaded_file)
                    if new_model:
                        st.session_state.models[model_label] = new_model
                        st.rerun()

        elif model_source == "URL":
            url = st.text_input(
                "Model URL",
                placeholder="e.g., http://bigg.ucsd.edu/models/e_coli_core.xml",
                help="Enter direct download URL for model file",
                key="multi_url_input"
            )
            model_label = st.text_input("Label for this model", value="url_model", key="multi_url_model_label")

            if st.button("Load from URL", key="multi_load_url_btn"):
                if model_label in st.session_state.models:
                    st.warning(f"Model with label '{model_label}' already exists. Please choose a different label.")
                else:
                    st.warning("URL loading requires additional implementation")

        elif model_source == "Create Blank Model":
            blank_model_label = st.text_input("Label for new blank model", value="new_model", key="multi_blank_model_label_input")
            if st.button("Create Blank Model", key="multi_create_blank_model_btn"):
                if blank_model_label in st.session_state.models:
                    st.warning(f"Model with label '{blank_model_label}' already exists. Please choose a different label.")
                else:
                    new_blank_model = cobra.Model(id=blank_model_label)
                    st.session_state.models[blank_model_label] = new_blank_model
                    st.success(f"Blank model '{blank_model_label}' created.")
                    st.rerun()

    st.subheader("Currently Loaded Models")
    if st.session_state.models:
        for name, model_obj in st.session_state.models.items():
            st.write(f"- **{name}**: {len(model_obj.reactions)} reactions, {len(model_obj.metabolites)} metabolites, {len(model_obj.genes)} genes")
            with st.expander(f"View Summary for {name}", expanded=False):
                display_model_summary(model_obj, unique_key_suffix=f"loaded_{name}")
            if st.button(f"Remove {name}", key=f"remove_model_{name}"):
                del st.session_state.models[name]
                if name in st.session_state.fba_solutions:
                    del st.session_state.fba_solutions[name]
                if st.session_state.current_model_name == name:
                    st.session_state.current_model_name = None
                st.success(f"Model '{name}' removed.")
                st.rerun()
    else:
        st.info("No models loaded yet.")

def multi_model_analysis_tabs():
    model_names = list(st.session_state.models.keys())
    
    tabs = [
        "Model Loading", "Compare Models"
    ]
    selected_tab = st.sidebar.selectbox("Navigation", tabs, key="multi_model_nav")

    if selected_tab == "Model Loading":
        render_model_loading_section_multi_mode()

    elif selected_tab == "Compare Models":
        st.header("📊 Compare Models")
        
        if len(model_names) < 2:
            st.warning("Please load at least two models to compare.")
            return

        col_m1, col_m2 = st.columns(2)
        with col_m1:
            model1_name = st.selectbox("Select Model 1", model_names, key="compare_model1")
        with col_m2:
            model2_options = [name for name in model_names if name != model1_name]
            if not model2_options:
                st.warning("Please load another model to compare.")
                return
            model2_name = st.selectbox("Select Model 2", model2_options, key="compare_model2")

        if model1_name and model2_name:
            model1 = st.session_state.models[model1_name]
            model2 = st.session_state.models[model2_name]

            st.subheader("Objective Values Comparison")
            if st.button("Run FBA and Compare Objectives", key="compare_objectives_btn"):
                with st.spinner("Running FBA for both models..."):
                    sol1 = perform_fba(model1.copy(), model1_name)
                    sol2 = perform_fba(model2.copy(), model2_name)

                if sol1 and sol2 and sol1.status == "optimal" and sol2.status == "optimal":
                    st.write(f"**{model1_name} Objective Value:** {sol1.objective_value:.4f}")
                    st.write(f"**{model2_name} Objective Value:** {sol2.objective_value:.4f}")

                    obj_data = pd.DataFrame({
                        "Model": [model1_name, model2_name],
                        "Objective Value": [sol1.objective_value, sol2.objective_value]
                    })
                    fig = px.bar(
                        obj_data,
                        x="Model",
                        y="Objective Value",
                        title="Objective Value Comparison",
                        color="Model",
                        text="Objective Value"
                    )
                    fig.update_traces(texttemplate='%{text:.4f}', textposition='outside')
                    st.plotly_chart(fig, use_container_width=True)
                else:
                    st.warning("Could not obtain optimal solutions for both models.")

            st.subheader("Top 10 Fluxes Comparison")
            if st.button("Compare Top 10 Fluxes", key="compare_top_fluxes_btn"):
                if model1_name not in st.session_state.fba_solutions or model2_name not in st.session_state.fba_solutions:
                    st.warning("Please run FBA for both models first (in Objective Values Comparison section).")
                else:
                    sol1 = st.session_state.fba_solutions[model1_name]
                    sol2 = st.session_state.fba_solutions[model2_name]

                    if sol1.status == "optimal" and sol2.status == "optimal":
                        flux1_df = sol1.fluxes.reset_index()
                        flux1_df.columns = ['Reaction ID', 'Flux']
                        flux1_df['Reaction ID'] = flux1_df['Reaction ID'].astype(str)
                        flux1_df['Reaction Name'] = flux1_df['Reaction ID'].apply(
                            lambda x: model1.reactions.get_by_id(x).name if x in model1.reactions else "N/A"
                        )
                        top_flux1 = flux1_df[abs(flux1_df['Flux']) > 1e-6].sort_values(
                            by='Flux', key=abs, ascending=False
                        ).head(10)
                        top_flux1['Model'] = model1_name

                        flux2_df = sol2.fluxes.reset_index()
                        flux2_df.columns = ['Reaction ID', 'Flux']
                        flux2_df['Reaction ID'] = flux2_df['Reaction ID'].astype(str)
                        flux2_df['Reaction Name'] = flux2_df['Reaction ID'].apply(
                            lambda x: model2.reactions.get_by_id(x).name if x in model2.reactions else "N/A"
                        )
                        top_flux2 = flux2_df[abs(flux2_df['Flux']) > 1e-6].sort_values(
                            by='Flux', key=abs, ascending=False
                        ).head(10)
                        top_flux2['Model'] = model2_name

                        combined_fluxes = pd.concat([top_flux1, top_flux2])

                        st.write(f"**Top 10 Fluxes for {model1_name}**")
                        st.dataframe(top_flux1[['Reaction ID', 'Reaction Name', 'Flux']], hide_index=True)
                        st.write(f"**Top 10 Fluxes for {model2_name}**")
                        st.dataframe(top_flux2[['Reaction ID', 'Reaction Name', 'Flux']], hide_index=True)

                        fig = px.bar(
                            combined_fluxes,
                            x='Reaction Name',
                            y='Flux',
                            color='Model',
                            barmode='group',
                            title="Top 10 Fluxes Comparison",
                            color_discrete_map={model1_name: 'blue', model2_name: 'orange'},
                            text='Flux'
                        )
                        fig.update_traces(texttemplate='%{text:.4f}', textposition='outside')
                        fig.update_layout(
                            xaxis_title="Reaction",
                            yaxis_title="Flux Value",
                            template='plotly_white',
                            height=600
                        )
                        fig.update_xaxes(tickangle=45)
                        st.plotly_chart(fig, use_container_width=True)
                    else:
                        st.warning("Optimal FBA solutions not available for both models.")

            st.subheader("Essential Gene Overlap")
            if st.button("Compare Essential Genes", key="compare_essential_genes_btn"):
                with st.spinner("Finding essential genes for both models..."):
                    
                    gene_objs1 = [model1.genes.get_by_id(g.id) for g in model1.genes]
                    essential_results1 = single_gene_deletion(model1.copy(), gene_list=gene_objs1)
                    essential_results1.index = [gene.id for gene in gene_objs1]
                    essential_results1.index.name = "Gene ID"

                    gene_names1 = [gene.name if gene.name else "N/A" for gene in gene_objs1]
                    essential_df1 = pd.DataFrame({
                        "Gene ID": essential_results1.index,
                        "Growth Rate": essential_results1['growth'],
                        "Gene Name": gene_names1
                    })
                    essential_df1["Status"] = essential_df1["Growth Rate"].apply(lambda x: "Lethal" if abs(x) < 1e-6 else "Non-Lethal")
                    essential_df1["Gene Display"] = essential_df1.apply(
                        lambda row: f"{row['Gene Name']} ({row['Gene ID']})" if row['Gene Name'] != "N/A" else row['Gene ID'], axis=1
                    )
                    lethal_genes1 = set(essential_df1[essential_df1["Status"] == "Lethal"]["Gene ID"])

                    gene_objs2 = [model2.genes.get_by_id(g.id) for g in model2.genes]
                    essential_results2 = single_gene_deletion(model2.copy(), gene_list=gene_objs2)
                    essential_results2.index = [gene.id for gene in gene_objs2]
                    essential_results2.index.name = "Gene ID"

                    gene_names2 = [gene.name if gene.name else "N/A" for gene in gene_objs2]
                    essential_df2 = pd.DataFrame({
                        "Gene ID": essential_results2.index,
                        "Growth Rate": essential_results2['growth'],
                        "Gene Name": gene_names2
                    })
                    essential_df2["Status"] = essential_df2["Growth Rate"].apply(lambda x: "Lethal" if abs(x) < 1e-6 else "Non-Lethal")
                    essential_df2["Gene Display"] = essential_df2.apply(
                        lambda row: f"{row['Gene Name']} ({row['Gene ID']})" if row['Gene Name'] != "N/A" else row['Gene ID'], axis=1
                    )
                    lethal_genes2 = set(essential_df2[essential_df2["Status"] == "Lethal"]["Gene ID"])

                common_lethal = lethal_genes1.intersection(lethal_genes2)
                unique_lethal1 = lethal_genes1.difference(lethal_genes2)
                unique_lethal2 = lethal_genes2.difference(lethal_genes1)

                st.write(f"**Common Lethal Genes ({model1_name} & {model2_name}):** {len(common_lethal)} genes")
                if common_lethal:
                    common_lethal_data = []
                    for g_id in common_lethal:
                        name1 = model1.genes.get_by_id(g_id).name if g_id in model1.genes else "N/A"
                        name2 = model2.genes.get_by_id(g_id).name if g_id in model2.genes else "N/A"
                        common_lethal_data.append({"Gene ID": str(g_id), "Gene Name (Model 1)": name1, "Gene Name (Model 2)": name2, "Gene Display (Model 1)": f"{name1} ({g_id})" if name1 != "N/A" else g_id, "Gene Display (Model 2)": f"{name2} ({g_id})" if name2 != "N/A" else g_id})
                    st.dataframe(pd.DataFrame(common_lethal_data), hide_index=True)
                else:
                    st.info("No common lethal genes found.")

                st.write(f"**Lethal Only in {model1_name}:** {len(unique_lethal1)} genes")
                if unique_lethal1:
                    unique_lethal1_data = []
                    for g_id in unique_lethal1:
                        name1 = model1.genes.get_by_id(g_id).name if g_id in model1.genes else "N/A"
                        unique_lethal1_data.append({"Gene ID": str(g_id), "Gene Name": name1, "Gene Display": f"{name1} ({g_id})" if name1 != "N/A" else g_id})
                    st.dataframe(pd.DataFrame(unique_lethal1_data), hide_index=True)
                else:
                    st.info(f"No genes lethal only in {model1_name}.")

                st.write(f"**Lethal Only in {model2_name}:** {len(unique_lethal2)} genes")
                if unique_lethal2:
                    unique_lethal2_data = []
                    for g_id in unique_lethal2:
                        name2 = model2.genes.get_by_id(g_id).name if g_id in model2.genes else "N/A"
                        unique_lethal2_data.append({"Gene ID": str(g_id), "Gene Name": name2, "Gene Display": f"{name2} ({g_id})" if name2 != "N/A" else g_id})
                    st.dataframe(pd.DataFrame(unique_lethal2_data), hide_index=True)
                else:
                    st.info(f"No genes lethal only in {model2_name}.")

            st.subheader("Model Content Differences")
            if st.button("Show Differences in Reactions/Metabolites/Genes", key="show_diff_btn"):
                st.write(f"**Reactions in {model1_name} but not {model2_name}:**")
                rxns1_ids = set(r.id for r in model1.reactions)
                rxns2_ids = set(r.id for r in model2.reactions)
                diff_rxns1 = rxns1_ids.difference(rxns2_ids)
                if diff_rxns1:
                    diff_rxn_data = []
                    for rxn_id in diff_rxns1:
                        rxn = model1.reactions.get_by_id(rxn_id)
                        diff_rxn_data.append({"ID": str(rxn.id), "Name": rxn.name if rxn.name else "N/A"})
                    st.dataframe(pd.DataFrame(diff_rxn_data), hide_index=True)
                else:
                    st.info("No unique reactions in Model 1.")

                st.write(f"**Reactions in {model2_name} but not {model1_name}:**")
                diff_rxns2 = rxns2_ids.difference(rxns1_ids)
                if diff_rxns2:
                    diff_rxn_data = []
                    for rxn_id in diff_rxns2:
                        rxn = model2.reactions.get_by_id(rxn_id)
                        diff_rxn_data.append({"ID": str(rxn.id), "Name": rxn.name if rxn.name else "N/A"})
                    st.dataframe(pd.DataFrame(diff_rxn_data), hide_index=True)
                else:
                    st.info("No unique reactions in Model 2.")
                
                st.write(f"**Metabolites in {model1_name} but not {model2_name}:**")
                mets1_ids = set(m.id for m in model1.metabolites)
                mets2_ids = set(m.id for m in model2.metabolites)
                diff_mets1 = mets1_ids.difference(mets2_ids)
                if diff_mets1:
                    diff_met_data = []
                    for met_id in diff_mets1:
                        met = model1.metabolites.get_by_id(met_id)
                        diff_met_data.append({"ID": str(met.id), "Name": met.name if met.name else "N/A"})
                    st.dataframe(pd.DataFrame(diff_met_data), hide_index=True)
                else:
                    st.info("No unique metabolites in Model 1.")

                st.write(f"**Metabolites in {model2_name} but not {model1_name}:**")
                diff_mets2 = mets2_ids.difference(mets1_ids)
                if diff_mets2:
                    diff_met_data = []
                    for met_id in diff_mets2:
                        met = model2.metabolites.get_by_id(met_id)
                        diff_met_data.append({"ID": str(met.id), "Name": met.name if met.name else "N/A"})
                    st.dataframe(pd.DataFrame(diff_met_data), hide_index=True)
                else:
                    st.info("No unique metabolites in Model 2.")

                st.write(f"**Genes in {model1_name} but not {model2_name}:**")
                genes1_ids = set(g.id for g in model1.genes)
                genes2_ids = set(g.id for g in model2.genes)
                diff_genes1 = genes1_ids.difference(genes2_ids)
                if diff_genes1:
                    diff_gene_data = []
                    for gene_id in diff_genes1:
                        gene = model1.genes.get_by_id(gene_id)
                        diff_gene_data.append({"ID": str(gene.id), "Name": gene.name if gene.name else "N/A"})
                    st.dataframe(pd.DataFrame(diff_gene_data), hide_index=True)
                else:
                    st.info("No unique genes in Model 1.")

                st.write(f"**Genes in {model2_name} but not {model1_name}:**")
                diff_genes2 = genes2_ids.difference(genes1_ids)
                if diff_genes2:
                    diff_gene_data = []
                    for gene_id in diff_genes2:
                        gene = model2.genes.get_by_id(gene_id)
                        diff_gene_data.append({"ID": str(gene.id), "Name": gene.name if gene.name else "N/A"})
                    st.dataframe(pd.DataFrame(diff_gene_data), hide_index=True)
                else:
                    st.info("No unique genes in Model 2.")

def render_intro_content():
    st.title("🧬 Welcome to COBRApy GUI")
    st.markdown("This Streamlit-based graphical user interface (GUI) is designed to simplify and enhance the interaction with `COBRApy`, a powerful Python package for constraint-based reconstruction and analysis of biological networks.")

    with st.expander("About"):
        st.markdown("""
        ### **Project Overview**
        `COBRApy_GUI` is a Streamlit-based graphical user interface (GUI) designed to simplify and enhance the interaction with `COBRApy`, a widely used Python package for constraint-based reconstruction and analysis of biological networks. This tool provides an visual environment for performing core metabolic modeling tasks, making complex analyses more accessible to researchers who may prefer a GUI over direct command-line interaction.

        ### **Features**
        This application offers a suite of functionalities, with key `COBRApy` capabilities:

        * **Model Management:**
            * **Loading:** Supports loading metabolic models from SBML (.xml), JSON (.json), and MATLAB (.mat) formats. Includes built-in access to common `cobra` models (e.g., E. coli core, iJO1366, Recon3D).
            * **Summary & Inspection:** Provides a detailed overview of loaded models, including counts of reactions, metabolites, and genes, along with objective function details and snippets of model components.
            * **Export:** Allows users to export modified models back into SBML, JSON, or MATLAB formats.
        * **Flux Analysis:**
            * **Flux Balance Analysis (FBA):** Perform FBA to determine optimal flux distributions, with options to set objective functions and adjust reaction bounds. Visualizes flux distributions and active pathways.
            * **Parsimonious FBA (pFBA):** Execute pFBA to identify the most parsimonious (minimal total flux) solution achieving optimal growth, offering insights into metabolic efficiency. Includes FBA vs. pFBA flux comparison.
            * **Flux Variability Analysis (FVA):** Calculate the minimum and maximum possible flux for individual or all reactions, providing insights into reaction flexibility and robustness. Visualizes flux ranges.
            * **Set Reaction Bounds:** Directly adjust the lower and upper bounds of individual reactions to simulate specific environmental conditions or genetic modifications.
        * **Genetic & Reaction Manipulations:**
            * **Gene Deletion Analysis:** Simulate the effects of single or double gene knockouts on the model's objective (e.g., growth rate), identifying essential genes and synthetic lethals.
            * **Reaction Deletion Analysis:** Perform single or double reaction deletions to assess their impact on model performance.
            * **Find Essential Genes:** A dedicated tool to identify genes whose deletion leads to a significant drop in growth rate (lethal).
            * **Gene/Reaction Set Simulation:** Allows defining and simulating the deletion of multiple custom gene sets or reaction sets (pathway blocks) in batch.
            * Visualizes growth impact of deletions and identifies lethal knockouts.
        * **Model Diagnostics:**
            * **Blocked Reactions:** Identify and optionally remove reactions that cannot carry any flux under defined conditions, aiding in model curation.
            * **Consistency Check:** Verify the solvability and consistency of the metabolic model.
        * **Reaction Details:**
            * **Explore individual reactions:** View detailed information about any reaction, including its equation, involved metabolites with their stoichiometry, and associated genes with their GPR rules.
        * **Data Export & Visualization:**
            * Download FBA flux distributions, FVA results, and gene/reaction deletion analysis reports as CSV files for external analysis.
            * Integrated Plotly visualizations for flux distributions, FVA ranges, and growth rate impacts.
        * **Multi-Model Support:**
            * Load and manage multiple models simultaneously.
            * **Compare Models:** Select two loaded models and compare their objective values, top 10 fluxes, essential gene overlap, and content differences (reactions, metabolites, genes).
        """)

def main():
    render_intro_content()

    st.sidebar.header("Choose Analysis Mode")
    st.session_state.analysis_mode = st.sidebar.radio(
        "Select Mode",
        ["Single Model Analysis", "Multi Model Analysis"],
        key="analysis_mode_selector"
    )

    if st.session_state.analysis_mode == "Single Model Analysis":
        single_model_analysis_tabs()
    elif st.session_state.analysis_mode == "Multi Model Analysis":
        multi_model_analysis_tabs()

if __name__ == "__main__":
    main()

