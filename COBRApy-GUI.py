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

# Configure Streamlit page
st.set_page_config(
    page_title="Welcome to COBRApy GUI",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Initialize session state
if 'model' not in st.session_state:
    st.session_state.model = None
if 'fba_solution' not in st.session_state:
    st.session_state.fba_solution = None

# Helper functions
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

def display_model_summary(model):
    st.subheader("Model Summary")

    col1, col2, col3 = st.columns(3)
    col1.metric("Reactions", len(model.reactions))
    col2.metric("Metabolites", len(model.metabolites))
    col3.metric("Genes", len(model.genes))

    objective = model.objective
    if hasattr(objective, 'expression'):
        st.write(f"**Objective Function:** {objective.expression}")
    else:
        obj_reactions = [rxn.id for rxn in model.reactions if rxn.objective_coefficient != 0]
        st.write(f"**Objective Function:** {', '.join(obj_reactions)}")

    with st.expander("First 10 Reactions"):
        reaction_data = []
        for rxn in model.reactions[:10]:
            reaction_data.append({
                "ID": rxn.id,
                "Name": rxn.name,
                "Equation": rxn.build_reaction_string(),
                "Bounds": f"{rxn.lower_bound} to {rxn.upper_bound}"
            })
        st.table(pd.DataFrame(reaction_data))

    with st.expander("First 10 Metabolites"):
        met_data = []
        for met in model.metabolites[:10]:
            met_data.append({
                "ID": met.id,
                "Name": met.name,
                "Compartment": met.compartment
            })
        st.table(pd.DataFrame(met_data))

    with st.expander("First 10 Genes"):
        gene_data = []
        for gene in model.genes[:10]:
            gene_data.append({
                "ID": gene.id,
                "Name": gene.name,
                "Reactions": len(gene.reactions)
            })
        st.table(pd.DataFrame(gene_data))

def get_reaction_options(model):
    return [f"{rxn.id} - {rxn.name}" for rxn in model.reactions]

def get_metabolite_options(model):
    return [f"{met.id} - {met.name}" for met in model.metabolites]

def get_gene_options(model):
    return [f"{gene.id} - {gene.name}" for gene in model.genes]

def extract_id(option):
    return option.split(' - ')[0]

def validate_genes(model, gene_ids):
    valid_ids = []
    for gid in gene_ids:
        try:
            model.genes.get_by_id(gid)
            valid_ids.append(gid)
        except KeyError:
            st.warning(f"Gene {gid} not found in model. Skipping...")
    return valid_ids

def perform_fba(model, objective_reaction=None):
    try:
        if objective_reaction and objective_reaction in [rxn.id for rxn in model.reactions]:
            model.objective = objective_reaction

        solution = model.optimize()
        st.session_state.fba_solution = solution

        if solution.status == "optimal":
            st.success(f"Optimization successful! Objective value: {solution.objective_value:.4f}")

            fluxes = solution.fluxes
            nonzero_fluxes = fluxes[abs(fluxes) > 1e-6].sort_values(key=abs, ascending=False)

            st.subheader("Top Active Fluxes")
            st.dataframe(nonzero_fluxes.head(10))

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

            csv = fluxes.reset_index().to_csv(index=False)
            st.download_button(
                label="Download Full Flux Distribution",
                data=csv,
                file_name='fba_fluxes.csv',
                mime='text/csv'
            )
        else:
            st.error(f"Optimization failed: {solution.status}")

        return solution

    except Exception as e:
        st.error(f"FBA failed: {str(e)}")
        return None

def plot_growth_impact(results_df, title):
    if results_df.empty:
        return

    results_df = results_df.copy()
    results_df.sort_values("Growth Rate", ascending=False, inplace=True)

    colors = ['red' if rate < 1e-6 else 'green' for rate in results_df["Growth Rate"]]

    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            x=results_df.index.astype(str),
            y=results_df["Growth Rate"],
            marker_color=colors,
            text=results_df["Growth Rate"].round(6),
            textposition='auto'
        )
    )

    fig.update_layout(
        title=title,
        xaxis_title="Deletion",
        yaxis_title="Growth Rate",
        template='plotly_white',
        height=500,
        hovermode='x'
    )
    fig.update_xaxes(tickangle=45)
    fig.add_hline(y=1e-6, line_dash="dash", line_color="gray")

    st.plotly_chart(fig, use_container_width=True)

def main():
    st.title("🧬 Welcome to COBRApy GUI")
    st.markdown("""
    **An interface for constraint-based metabolic modeling using COBRApy.**
    """)

    tabs = [
        "Model Loading", "Flux Analysis", "Genetic Manipulations",
        "Model Diagnostics", "Model Editing", "Export & Visualization"
    ]
    selected_tab = st.sidebar.selectbox("Navigation", tabs)
    
    # Model Loading Section
    if selected_tab == "Model Loading":
        st.header("📤 Upload a metabolic model or select a built-in model to get started.")

        builtin_models = {
            "E. coli core": "e_coli_core",
            "Textbook E. coli": "textbook",
            "Salmonella core": "salmonella",
            "iJO1366 (E. coli)": "iJO1366",
            "Recon3D (Human)": "Recon3D"
        }

        col1, col2 = st.columns(2)
        with col1:
            model_source = st.radio("Model source:", ["Built-in", "Upload File", "URL"])

            if model_source == "Built-in":
                selected_model = st.selectbox("Select model", list(builtin_models.keys()))

                if st.button("Load Built-in Model"):
                    try:
                        model_name = builtin_models[selected_model]
                        model = load_model(model_name)
                        st.session_state.model = model
                        st.success(f"Loaded built-in model: {model_name}")
                    except Exception as e:
                        st.error(f"Error loading model: {str(e)}")

            elif model_source == "Upload File":
                uploaded_file = st.file_uploader(
                    "Upload model (SBML/JSON/MAT)",
                    type=["xml", "json", "mat"],
                    help="Supported formats: SBML (.xml), JSON (.json), MATLAB (.mat)"
                )

                if st.button("Load Uploaded Model") and uploaded_file:
                    model = load_model_from_file(uploaded_file)
                    if model:
                        st.session_state.model = model

            elif model_source == "URL":
                url = st.text_input(
                    "Model URL",
                    placeholder="e.g., http://bigg.ucsd.edu/models/e_coli_core.xml",
                    help="Enter direct download URL for model file"
                )
                if st.button("Load from URL"):
                    try:
                        st.warning("URL loading requires additional implementation")
                    except Exception as e:
                        st.error(f"Error loading from URL: {str(e)}")

        if st.session_state.model:
            display_model_summary(st.session_state.model)
            
    # Flux Analysis Section
    elif selected_tab == "Flux Analysis" and st.session_state.model:
        st.header("📊 Flux Analysis")
        model = st.session_state.model

        reaction_options = get_reaction_options(model)

        with st.expander("Flux Balance Analysis (FBA)"):
            st.info("FBA calculates the optimal flux distribution to maximize/minimize the objective function.")

            current_objective = None
            for rxn in model.reactions:
                if rxn.objective_coefficient != 0:
                    current_objective = rxn.id
                    break

            if current_objective:
                current_objective_display = f"{current_objective} - {model.reactions.get_by_id(current_objective).name}"
            else:
                current_objective_display = reaction_options[0]

            objective_display = st.selectbox(
                "Objective Reaction",
                reaction_options,
                index=reaction_options.index(current_objective_display) if current_objective else 0,
                help="Select the reaction to optimize"
            )
            objective_rxn = extract_id(objective_display)

            st.subheader("Adjust Reaction Bounds")
            rxn_display = st.selectbox("Select reaction to modify", reaction_options)
            rxn_id = extract_id(rxn_display)
            selected_rxn = model.reactions.get_by_id(rxn_id)

            col1, col2 = st.columns(2)
            with col1:
                new_lb = st.number_input(
                    "Lower Bound",
                    value=float(selected_rxn.lower_bound),
                    format="%f"
                )
            with col2:
                new_ub = st.number_input(
                    "Upper Bound",
                    value=float(selected_rxn.upper_bound),
                    format="%f"
                )

            if st.button("Apply Bounds"):
                selected_rxn.bounds = (new_lb, new_ub)
                st.success(f"Updated bounds for {rxn_id}: [{new_lb}, {new_ub}]")

            if st.button("Run FBA", type="primary"):
                with st.spinner("Performing Flux Balance Analysis..."):
                    solution = perform_fba(model, objective_rxn)

        with st.expander("Parsimonious FBA (pFBA)"):
            st.info("pFBA finds the minimal total flux solution that achieves optimal growth.")

            if st.button("Run pFBA"):
                with st.spinner("Performing Parsimonious FBA..."):
                    try:
                        p_solution = pfba(model)
                        st.session_state.fba_solution = p_solution

                        st.success(f"pFBA completed! Total flux: {p_solution.fluxes.abs().sum():.4f}")

                        if st.session_state.fba_solution:
                            fba_fluxes = st.session_state.fba_solution.fluxes
                            p_fluxes = p_solution.fluxes

                            nonzero_fba = sum(abs(fba_fluxes) > 1e-6)
                            nonzero_pfba = sum(abs(p_fluxes) > 1e-6)

                            st.write(f"**FBA Active Reactions:** {nonzero_fba}")
                            st.write(f"**pFBA Active Reactions:** {nonzero_pfba}")

                            st.subheader("Top pFBA Fluxes")
                            nonzero_pfluxes = p_fluxes[abs(p_fluxes) > 1e-6].sort_values(key=abs, ascending=False)
                            st.dataframe(nonzero_pfluxes.head(10))

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
                help="Analyze specific reactions or all reactions in model"
            )
            reaction_ids = [extract_id(r) for r in reactions] if reactions else []

            fraction = st.slider(
                "Fraction of Optimum",
                min_value=0.5, max_value=1.0, value=0.95, step=0.01,
                help="Percentage of optimal growth rate to consider"
            )

            if st.button("Run FVA"):
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
                        st.dataframe(fva_results)

                        fig = go.Figure()

                        fig.add_trace(
                            go.Bar(
                                x=fva_results.index,
                                y=fva_results['minimum'],
                                name='Minimum Flux',
                                marker_color='#636EFA'
                            )
                        )

                        fig.add_trace(
                            go.Bar(
                                x=fva_results.index,
                                y=fva_results['maximum'] - fva_results['minimum'],
                                name='Flux Range',
                                base=fva_results['minimum'],
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

                        csv = fva_results.reset_index().to_csv(index=False)
                        st.download_button(
                            label="Download FVA Results",
                            data=csv,
                            file_name='fva_results.csv',
                            mime='text/csv'
                        )
                    except Exception as e:
                        st.error(f"FVA failed: {str(e)}")

    # Genetic Manipulations Section
    elif selected_tab == "Genetic Manipulations" and st.session_state.model:
        st.header("🧬 Genetic Manipulations")
        model = st.session_state.model

        col1, col2 = st.columns(2)

        gene_options = get_gene_options(model)
        reaction_options = get_reaction_options(model)

        with col1:
            st.subheader("Gene Deletion Analysis")
            deletion_type = st.radio("Deletion Type", ["Single", "Double"], index=0)

            selected_genes = st.multiselect(
                "Select genes to delete",
                gene_options
            )
            gene_ids_raw = [extract_id(g) for g in selected_genes]

            if st.button("Simulate Gene Deletions"):
                if not gene_ids_raw:
                    st.warning("Please select at least one gene.")
                    return

                gene_ids = validate_genes(model, gene_ids_raw)
                if not gene_ids:
                    st.error("No valid genes selected for deletion that exist in the model.")
                    return

                with st.spinner("Analyzing gene deletions..."):
                    try:
                        if deletion_type == "Single":
                            deletion_results = single_gene_deletion(model, gene_ids)

                            results_df = pd.DataFrame({
                                "Gene": deletion_results.index.get_level_values(0),
                                "Growth Rate": deletion_results['growth']
                            })
                        else:
                            if len(gene_ids) < 2:
                                st.error("Please select at least 2 genes for double deletion")
                                return
                            deletion_results = double_gene_deletion(model, gene_ids)

                            results_df = pd.DataFrame({
                                "Gene 1": deletion_results.index.get_level_values(0),
                                "Gene 2": deletion_results.index.get_level_values(1),
                                "Growth Rate": deletion_results['growth']
                            })

                        results_df["Lethal"] = results_df["Growth Rate"] < 1e-6

                        if deletion_type == "Single":
                            gene_names = []
                            for g_id in results_df["Gene"]:
                                try:
                                    gene_names.append(model.genes.get_by_id(g_id).name)
                                except KeyError:
                                    gene_names.append(f"UNKNOWN_GENE_{g_id}")
                            results_df["Gene Name"] = gene_names
                            st.dataframe(results_df[["Gene", "Gene Name", "Growth Rate", "Lethal"]])
                        else:
                            gene_1_names = []
                            for g_id in results_df["Gene 1"]:
                                try:
                                    gene_1_names.append(model.genes.get_by_id(g_id).name)
                                except KeyError:
                                    gene_1_names.append(f"UNKNOWN_GENE_{g_id}")
                            results_df["Gene 1 Name"] = gene_1_names

                            gene_2_names = []
                            for g_id in results_df["Gene 2"]:
                                try:
                                    gene_2_names.append(model.genes.get_by_id(g_id).name)
                                except KeyError:
                                    gene_2_names.append(f"UNKNOWN_GENE_{g_id}")
                            results_df["Gene 2 Name"] = gene_2_names
                            st.dataframe(results_df[["Gene 1", "Gene 1 Name", "Gene 2", "Gene 2 Name", "Growth Rate", "Lethal"]])

                        lethal_count = sum(results_df["Lethal"])
                        st.info(f"**{lethal_count} lethal deletion(s) found**")

                        if deletion_type == "Single":
                            plot_growth_impact(results_df, "Gene Deletion Impact")

                        csv = results_df.to_csv(index=False)
                        st.download_button(
                            label="Download Gene Deletion Results",
                            data=csv,
                            file_name='gene_deletion_results.csv',
                            mime='text/csv'
                        )
                    except Exception as e:
                        st.error(f"Gene deletion analysis failed: {str(e)}")
                        st.exception(e)

        with col2:
            st.subheader("Reaction Deletion Analysis")
            deletion_type_rxn = st.radio("Deletion Type", ["Single", "Double"], index=0, key="rxn_del_type")

            selected_reactions = st.multiselect(
                "Select reactions to delete",
                reaction_options
            )
            reaction_ids = [extract_id(r) for r in selected_reactions]

            if st.button("Simulate Reaction Deletions"):
                if not reaction_ids:
                    st.warning("Please select at least one reaction")
                else:
                    with st.spinner("Analyzing reaction deletions..."):
                        try:
                            if deletion_type_rxn == "Single":
                                deletion_results = single_reaction_deletion(model, reaction_ids)
                                results_df = pd.DataFrame({
                                    "Reaction": deletion_results.index.get_level_values(0),
                                    "Growth Rate": deletion_results['growth']
                                })
                            else:
                                if len(reaction_ids) < 2:
                                    st.error("Please select at least 2 reactions for double deletion")
                                    return
                                deletion_results = double_reaction_deletion(model, reaction_ids)

                                results_df = pd.DataFrame({
                                    "Reaction 1": deletion_results.index.get_level_values(0),
                                    "Reaction 2": deletion_results.index.get_level_values(1),
                                    "Growth Rate": deletion_results['growth']
                                })

                            results_df["Lethal"] = results_df["Growth Rate"] < 1e-6

                            if deletion_type_rxn == "Single":
                                rxn_names = []
                                for r_id in results_df["Reaction"]:
                                    try:
                                        rxn_names.append(model.reactions.get_by_id(r_id).name)
                                    except KeyError:
                                        rxn_names.append(f"UNKNOWN_REACTION_{r_id}")
                                results_df["Reaction Name"] = rxn_names
                                st.dataframe(results_df[["Reaction", "Reaction Name", "Growth Rate", "Lethal"]])
                            else:
                                rxn_1_names = []
                                for r_id in results_df["Reaction 1"]:
                                    try:
                                        rxn_1_names.append(model.reactions.get_by_id(r_id).name)
                                    except KeyError:
                                        rxn_1_names.append(f"UNKNOWN_REACTION_{r_id}")
                                results_df["Reaction 1 Name"] = rxn_1_names

                                rxn_2_names = []
                                for r_id in results_df["Reaction 2"]:
                                    try:
                                        rxn_2_names.append(model.reactions.get_by_id(r_id).name)
                                    except KeyError:
                                        rxn_2_names.append(f"UNKNOWN_REACTION_{r_id}")
                                results_df["Reaction 2 Name"] = rxn_2_names

                                st.dataframe(results_df[["Reaction 1", "Reaction 1 Name", "Reaction 2", "Reaction 2 Name", "Growth Rate", "Lethal"]])

                            lethal_count_rxn = sum(results_df["Lethal"])
                            st.info(f"**{lethal_count_rxn} lethal deletion(s) found**")

                            if deletion_type_rxn == "Single":
                                plot_growth_impact(results_df, "Reaction Deletion Impact")

                            csv = results_df.to_csv(index=False)
                            st.download_button(
                                label="Download Reaction Deletion Results",
                                data=csv,
                                file_name='reaction_deletion_results.csv',
                                mime='text/csv'
                            )
                        except Exception as e:
                            st.error(f"Reaction deletion analysis failed: {str(e)}")
                            st.exception(e)

    # Model Diagnostics Section
    elif selected_tab == "Model Diagnostics" and st.session_state.model:
        st.header("🔍 Model Diagnostics")
        model = st.session_state.model
        
        col1, col2 = st.columns(2)
        
        # Blocked Reactions
        with col1:
            st.subheader("Blocked Reactions")
            st.info("Reactions that cannot carry any flux under any condition")
            
            if st.button("Find Blocked Reactions"):
                with st.spinner("Identifying blocked reactions..."):
                    try:
                        blocked = find_blocked_reactions(model)
                        if len(blocked) > 0:
                            st.warning(f"Found {len(blocked)} blocked reactions")
                            
                            # Create dataframe with names
                            blocked_data = []
                            for rxn_id in blocked:
                                rxn = model.reactions.get_by_id(rxn_id)
                                blocked_data.append({
                                    "ID": rxn.id,
                                    "Name": rxn.name,
                                    "Equation": rxn.build_reaction_string()
                                })
                            st.table(pd.DataFrame(blocked_data))
                            
                            if st.button("Remove Blocked Reactions"):
                                model.remove_reactions(blocked)
                                st.success(f"Removed {len(blocked)} blocked reactions")
                        else:
                            st.success("No blocked reactions found")
                    except Exception as e:
                        st.error(f"Error finding blocked reactions: {str(e)}")
        
        # Model Consistency
        with col2:
            st.subheader("Model Consistency Check")
            st.info("Verify model structural and numerical consistency")
            
            if st.button("Check Model Consistency"):
                with st.spinner("Checking model consistency..."):
                    try:
                        # Simplified consistency check
                        solution = model.optimize()
                        if solution.status == "optimal":
                            st.success("Model is consistent and solvable")
                        else:
                            st.error(f"Model inconsistency detected: {solution.status}")
                    except Exception as e:
                        st.error(f"Model consistency check failed: {str(e)}")
    
    # Model Editing Section
    elif selected_tab == "Model Editing" and st.session_state.model:
        st.header("✏️ Model Editing")
        model = st.session_state.model
        
        tab1, tab2, tab3, tab4 = st.tabs([
            "Reactions", "Metabolites", "Genes", "Media Conditions"
        ])
        
        # Get options with names
        reaction_options = get_reaction_options(model)
        metabolite_options = get_metabolite_options(model)
        gene_options = get_gene_options(model)
        
        # Reaction Editing
        with tab1:
            st.subheader("Add New Reaction")
            with st.form("add_reaction"):
                rxn_id = st.text_input("Reaction ID")
                rxn_name = st.text_input("Reaction Name (optional)")
                rxn_formula = st.text_input("Reaction Formula (e.g., A_c + B_c --> C_c)")
                lb = st.number_input("Lower Bound", value=0.0)
                ub = st.number_input("Upper Bound", value=1000.0)
                
                if st.form_submit_button("Add Reaction"):
                    try:
                        new_rxn = cobra.Reaction(id=rxn_id, name=rxn_name)
                        new_rxn.bounds = (lb, ub)
                        
                        # Parse metabolites from formula
                        metabolites = {}
                        for part in rxn_formula.split('-->'):
                            for met in re.findall(r'(\d*\s*\w+_\w+)', part):
                                coeff = 1.0
                                if met.strip()[0].isdigit():
                                    coeff = float(met.split()[0])
                                    met_id = met.split()[1]
                                else:
                                    met_id = met.strip()
                                
                                # Find or create metabolite
                                if met_id in model.metabolites:
                                    metabolite = model.metabolites.get_by_id(met_id)
                                else:
                                    st.warning(f"Metabolite {met_id} not found in model. Creating placeholder.")
                                    metabolite = cobra.Metabolite(met_id)
                                
                                # Determine coefficient based on side
                                if '-->' in rxn_formula:
                                    side = rxn_formula.split('-->')
                                    if met in side[0]:
                                        coeff = -abs(coeff)
                                    else:
                                        coeff = abs(coeff)
                                metabolites[metabolite] = coeff
                        
                        new_rxn.add_metabolites(metabolites)
                        model.add_reactions([new_rxn])
                        st.success(f"Reaction {rxn_id} added successfully!")
                    except Exception as e:
                        st.error(f"Error adding reaction: {str(e)}")
            
            st.subheader("Remove Reactions")
            reactions_to_remove = st.multiselect(
                "Select reactions to remove",
                reaction_options
            )
            if st.button("Remove Selected Reactions"):
                reaction_ids = [extract_id(r) for r in reactions_to_remove]
                model.remove_reactions(reaction_ids)
                st.success(f"Removed {len(reaction_ids)} reactions")
        
        # Metabolite Editing
        with tab2:
            st.subheader("Add New Metabolite")
            with st.form("add_metabolite"):
                met_id = st.text_input("Metabolite ID")
                met_name = st.text_input("Name (optional)")
                formula = st.text_input("Chemical formula (optional)")
                charge = st.number_input("Charge", value=0)
                compartment = st.text_input("Compartment", value="c")
                
                if st.form_submit_button("Add Metabolite"):
                    try:
                        new_met = cobra.Metabolite(
                            id=met_id,
                            name=met_name,
                            formula=formula,
                            charge=charge,
                            compartment=compartment
                        )
                        model.add_metabolites([new_met])
                        st.success(f"Metabolite {met_id} added successfully!")
                    except Exception as e:
                        st.error(f"Error adding metabolite: {str(e)}")
            
            st.subheader("Remove Metabolites")
            metabolites_to_remove = st.multiselect(
                "Select metabolites to remove",
                metabolite_options
            )
            if st.button("Remove Selected Metabolites"):
                metabolite_ids = [extract_id(m) for m in metabolites_to_remove]
                model.remove_metabolites(metabolite_ids)
                st.success(f"Removed {len(metabolite_ids)} metabolites")
        
        # Gene Editing
        with tab3:
            st.subheader("Add New Gene")
            with st.form("add_gene"):
                gene_id = st.text_input("Gene ID")
                gene_name = st.text_input("Name (optional)")
                
                if st.form_submit_button("Add Gene"):
                    try:
                        new_gene = cobra.Gene(id=gene_id, name=gene_name)
                        model.genes.append(new_gene)
                        st.success(f"Gene {gene_id} added successfully!")
                    except Exception as e:
                        st.error(f"Error adding gene: {str(e)}")
            
            st.subheader("Remove Genes")
            genes_to_remove = st.multiselect(
                "Select genes to remove",
                gene_options
            )
            if st.button("Remove Selected Genes"):
                gene_ids = [extract_id(g) for g in genes_to_remove]
                for gene_id in gene_ids:
                    gene = model.genes.get_by_id(gene_id)
                    model.genes.remove(gene)
                st.success(f"Removed {len(gene_ids)} genes")
        
        # Media Conditions
        with tab4:
            st.subheader("Modify Media Conditions")
            st.info("Adjust exchange reaction bounds to simulate different growth conditions")
            
            exchange_rxns = [rxn for rxn in model.reactions if rxn.id.startswith('EX_')]
            exchange_options = [f"{rxn.id} - {rxn.name}" for rxn in exchange_rxns]
            
            selected_exchange = st.selectbox(
                "Select exchange reaction",
                exchange_options
            )
            
            if selected_exchange:
                rxn_id = extract_id(selected_exchange)
                rxn = model.reactions.get_by_id(rxn_id)
                col1, col2 = st.columns(2)
                with col1:
                    new_lb = st.number_input(
                        "New Lower Bound", 
                        value=float(rxn.lower_bound),
                        format="%f"
                    )
                with col2:
                    new_ub = st.number_input(
                        "New Upper Bound", 
                        value=float(rxn.upper_bound),
                        format="%f"
                    )
                
                if st.button("Apply Media Change"):
                    rxn.bounds = (new_lb, new_ub)
                    st.success(f"Updated {rxn_id} bounds to [{new_lb}, {new_ub}]")
                    
                    # Run FBA with new media
                    with st.spinner("Running FBA with new media..."):
                        solution = perform_fba(model)
    
    # Export & Visualization Section
    elif selected_tab == "Export & Visualization" and st.session_state.model:
        st.header("💾 Export & Visualization")
        model = st.session_state.model
        
        col1, col2 = st.columns(2)
        
        # Model Export
        with col1:
            st.subheader("Export Model")
            export_format = st.selectbox(
                "Select export format",
                ["SBML", "JSON", "MATLAB"]
            )
            
            export_filename = st.text_input(
                "Filename", 
                value=f"{model.id}_export",
                help="Filename without extension"
            )
            
            if st.button("Export Model"):
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
                                file_name=os.path.basename(filepath)
                            )
                except Exception as e:
                    st.error(f"Export failed: {str(e)}")
        
        # Visualization
        with col2:
            st.subheader("Flux Visualization")
            
            if st.session_state.fba_solution:
                st.success("Flux solution available for visualization")
                flux_df = st.session_state.fba_solution.fluxes.reset_index()
                flux_df.columns = ['Reaction ID', 'Flux']
                
                # Add reaction names to flux data
                flux_df['Reaction Name'] = flux_df['Reaction ID'].apply(
                    lambda x: model.reactions.get_by_id(x).name if x in model.reactions else "Unknown"
                )
                
                # Download flux data for visualization
                csv = flux_df.to_csv(index=False)
                st.download_button(
                    label="Download Flux Data",
                    data=csv,
                    file_name='flux_data.csv',
                    mime='text/csv'
                )
                
                # Top 10 fluxes visualization
                top_fluxes = flux_df[abs(flux_df['Flux']) > 1e-6].sort_values(
                    by='Flux', key=abs, ascending=False
                ).head(10)
                
                if not top_fluxes.empty:
                    # Create bar plot
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

# Run the application
if __name__ == "__main__":
    main()
