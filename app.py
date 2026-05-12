import streamlit as st
import streamlit.components.v1 as components
import scanpy as sc
import squidpy as sq
import anndata as ad
import pandas as pd
import numpy as np
import plotly.express as px
from io import BytesIO
import json
from html import escape
from scipy import sparse as sp

st.set_page_config(page_title="Spatial Transcriptomics Analyzer", layout="wide")

st.title("🧬 Spatial Transcriptomics Analyzer")

# Sidebar: Data ingestion
st.sidebar.header("📊 Data Loading")

# Data type selector
data_type = st.sidebar.radio("Data Type", ["Single Cell", "Spatial"], index=0)

@st.cache_data
def load_sample_data():
    """Load sample PBMC dataset from Scanpy."""
    adata = sc.datasets.pbmc68k_reduced()
    return adata

@st.cache_data
def load_sample_spatial_data():
    """Load sample spatial dataset from Squidpy."""
    adata = sq.datasets.visium_hne_adata()
    return adata

@st.cache_data
def load_uploaded_file(uploaded_file):
    """Load AnnData file from upload."""
    try:
        adata = sc.read_h5ad(uploaded_file)
        return adata
    except Exception as e:
        st.error(f"Error loading file: {e}")
        return None

# File upload or sample data
uploaded_file = st.sidebar.file_uploader("Upload .h5ad file", type=["h5ad"])

if uploaded_file is not None:
    adata = load_uploaded_file(uploaded_file)
else:
    if data_type == "Single Cell":
        if st.sidebar.button("Load Sample Dataset (PBMC 68K)"):
            adata = load_sample_data()
            st.sidebar.success("✅ Sample single-cell dataset loaded!")
    else:
        if st.sidebar.button("Load Sample Dataset (Visium HnE)"):
            adata = load_sample_spatial_data()
            st.sidebar.success("✅ Sample spatial dataset loaded!")

# Initialize adata in session state
if "adata" not in st.session_state:
    st.session_state.adata = None

if "adata_processed" not in st.session_state:
    st.session_state.adata_processed = None

if "adata_spatial" not in st.session_state:
    st.session_state.adata_spatial = None

# Update session state if adata was loaded
if 'adata' in locals():
    st.session_state.adata = adata

# Scanpy Preprocessing Function
@st.cache_data
def run_scanpy_pipeline(_adata):
    """Run standard Scanpy preprocessing pipeline."""
    adata_copy = _adata.copy()

    def _sanitize_x(_adata_obj):
        if sp.issparse(_adata_obj.X):
            _adata_obj.X = _adata_obj.X.tocsr(copy=True)
            _adata_obj.X.data = np.nan_to_num(_adata_obj.X.data, nan=0.0, posinf=0.0, neginf=0.0)
        else:
            _adata_obj.X = np.nan_to_num(np.asarray(_adata_obj.X), nan=0.0, posinf=0.0, neginf=0.0)

    def _set_manual_hvg(_adata_for_hvg):
        """Fallback HVG stats when scanpy binning fails."""
        if sp.issparse(_adata_for_hvg.X):
            x = _adata_for_hvg.X.toarray()
        else:
            x = np.asarray(_adata_for_hvg.X)

        x = np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)
        gene_means = x.mean(axis=0)
        gene_vars = x.var(axis=0)

        n_vars = _adata_for_hvg.n_vars
        top_k = min(2000, n_vars)
        ranked_idx = np.argsort(gene_vars)[::-1]
        hv_idx = ranked_idx[:top_k]

        highly_variable = np.zeros(n_vars, dtype=bool)
        highly_variable[hv_idx] = True

        hvg_rank = np.full(n_vars, np.nan)
        hvg_rank[hv_idx] = np.arange(top_k, dtype=float)

        _adata_for_hvg.var["highly_variable"] = highly_variable
        _adata_for_hvg.var["highly_variable_rank"] = hvg_rank
        _adata_for_hvg.var["means"] = gene_means
        _adata_for_hvg.var["dispersions_norm"] = gene_vars

    # Sanitize expression matrix and remove all-zero observations/features.
    _sanitize_x(adata_copy)
    sc.pp.filter_cells(adata_copy, min_counts=1)
    sc.pp.filter_genes(adata_copy, min_counts=1)

    if adata_copy.n_obs == 0 or adata_copy.n_vars == 0:
        raise ValueError("Dataset has no non-zero cells/genes after filtering; cannot run PCA/UMAP.")

    # Preprocessing steps
    sc.pp.normalize_total(adata_copy, inplace=True)
    _sanitize_x(adata_copy)
    sc.pp.log1p(adata_copy)
    _sanitize_x(adata_copy)
    try:
        sc.pp.highly_variable_genes(adata_copy)
    except Exception:
        # Retry with fewer bins for low-variance or degenerate datasets.
        try:
            sc.pp.highly_variable_genes(adata_copy, n_bins=3)
        except Exception:
            _set_manual_hvg(adata_copy)

    _sanitize_x(adata_copy)

    # PCA, neighbors, UMAP
    sc.tl.pca(adata_copy)
    sc.pp.neighbors(adata_copy, n_neighbors=15)
    sc.tl.umap(adata_copy)

    # Leiden clustering for annotations
    sc.tl.leiden(adata_copy, flavor="igraph")

    return adata_copy

# Vitessce Visualization Function
def create_vitessce_html(_adata):
    """Create Vitessce HTML visualization from AnnData."""
    try:
        from vitessce import (
            VitessceConfig,
            Component as VitessceComponent,
            AnnDataWrapper,
        )
        try:
            from vitessce import CoordinationType as CT
        except Exception:
            CT = None
        import tempfile
        import os
        import inspect

        # Create a temporary directory for file storage
        temp_dir = tempfile.mkdtemp()
        h5ad_path = os.path.join(temp_dir, "adata.h5ad")

        # Prepare AnnData for robust H5AD export (nullable strings can fail on older anndata)
        adata_export = _adata.copy()
        for df in (adata_export.obs, adata_export.var):
            for col in df.columns:
                if isinstance(df[col].dtype, pd.StringDtype):
                    df[col] = df[col].astype(object)

        # Save AnnData to h5ad file
        prev_nullable_setting = ad.settings.allow_write_nullable_strings
        ad.settings.allow_write_nullable_strings = True
        try:
            adata_export.write_h5ad(h5ad_path)
        finally:
            ad.settings.allow_write_nullable_strings = prev_nullable_setting

        # Create Vitessce config (compatible with versions where schema_version is required)
        config_kwargs = {
            "name": "Spatial Transcriptomics",
            "description": "Interactive UMAP and spatial visualization",
        }
        vc_init_sig = inspect.signature(VitessceConfig.__init__)
        if "schema_version" in vc_init_sig.parameters:
            schema_version = None
            try:
                from vitessce.constants import DEFAULT_SCHEMA_VERSION
                schema_version = DEFAULT_SCHEMA_VERSION
            except Exception:
                schema_version = "1.0.15"
            config_kwargs["schema_version"] = schema_version

        vc = VitessceConfig(**config_kwargs)

        # Create dataset with AnnDataWrapper pointing to h5ad file
        dataset = vc.add_dataset(name="Spatial Data", uid="spatial_data")

        ann_wrapper_base_kwargs = {
            "obs_set_paths": ["leiden"],
            "obsm_names": ["X_umap"],
        }
        ann_wrapper_sig = inspect.signature(AnnDataWrapper.__init__)
        has_var_kwargs = any(
            p.kind == inspect.Parameter.VAR_KEYWORD
            for p in ann_wrapper_sig.parameters.values()
        )
        param_names = set(ann_wrapper_sig.parameters.keys())

        source_candidates = [
            ("_path", h5ad_path),
            ("path", h5ad_path),
            ("adata_path", h5ad_path),
            ("_url", f"file://{h5ad_path}"),
            ("url", f"file://{h5ad_path}"),
            ("adata_url", f"file://{h5ad_path}"),
        ]

        wrapper_errors = []
        ann_wrapper = None
        for source_key, source_value in source_candidates:
            if not has_var_kwargs and source_key not in param_names:
                continue
            try:
                ann_wrapper = AnnDataWrapper(
                    **ann_wrapper_base_kwargs,
                    **{source_key: source_value},
                )
                dataset.add_object(ann_wrapper)
                break
            except Exception as wrapper_error:
                wrapper_errors.append(f"{source_key}: {wrapper_error}")

        if ann_wrapper is None:
            raise RuntimeError(
                "Failed to initialize AnnDataWrapper with a valid data input parameter. "
                f"Tried keys: {[k for k, _ in source_candidates]}. "
                f"Errors: {' | '.join(wrapper_errors) if wrapper_errors else 'No matching constructor keys.'}"
            )

        def add_view_compat(component, view_name, dataset_uid):
            add_view_attempts = [
                {"name": view_name, "dataset_uid": dataset_uid},
                {"dataset_uid": dataset_uid},
                {"name": view_name, "dataset": dataset},
                {"dataset": dataset},
                {"name": view_name},
                {},
            ]
            last_error = None
            for kwargs in add_view_attempts:
                try:
                    return vc.add_view(component, **kwargs)
                except TypeError as e:
                    last_error = e
                    continue

            raise RuntimeError(f"Unable to create Vitessce view for {component}: {last_error}")

        def resolve_component_compat(*component_candidates):
            for candidate in component_candidates:
                if hasattr(VitessceComponent, candidate):
                    return getattr(VitessceComponent, candidate)
            return None

        # Add scatterplot view for UMAP visualization
        scatterplot = add_view_compat(
            VitessceComponent.SCATTERPLOT,
            view_name="UMAP",
            dataset_uid="spatial_data",
        )

        # Configure embedding if supported by the installed Vitessce API.
        if CT is not None and hasattr(scatterplot, "set_coordination_value"):
            scatterplot.set_coordination_value(CT.EMBEDDING_TYPE, "X_umap")
        elif hasattr(scatterplot, "set_props"):
            try:
                scatterplot.set_props(embeddingType="X_umap")
            except Exception:
                pass

        # Add obs/cell sets view for cluster annotations (schema-version compatible)
        obs_sets_component = resolve_component_compat("OBS_SETS", "CELL_SETS")
        if obs_sets_component is None:
            obs_sets_component = "obsSets"

        cell_sets = add_view_compat(
            obs_sets_component,
            view_name="Clusters",
            dataset_uid="spatial_data",
        )

        # Create side-by-side layout
        vc.layout((scatterplot | cell_sets))

        def as_inline_html(render_output):
            if isinstance(render_output, str) and render_output.strip().lower().startswith(("http://", "https://")):
                safe_url = escape(render_output.strip(), quote=True)
                return (
                    f"<iframe src='{safe_url}' style='width:100%;height:780px;border:none;' "
                    "allow='clipboard-read; clipboard-write'></iframe>"
                )
            return render_output

        # Prefer HTML-producing APIs and always convert URL outputs to inline iframe HTML.
        if hasattr(vc, "to_html"):
            try:
                return as_inline_html(vc.to_html())
            except Exception:
                pass

        if hasattr(vc, "widget"):
            try:
                return as_inline_html(vc.widget(theme="light").to_html())
            except Exception:
                pass

        # Fallback for versions exposing web_app only.
        if hasattr(vc, "web_app"):
            try:
                return as_inline_html(vc.web_app(theme="light"))
            except TypeError:
                return as_inline_html(vc.web_app())

        raise RuntimeError("Unable to render Vitessce inline: no supported HTML/widget/web_app render method found.")

    except ModuleNotFoundError as e:
        if "anywidget" in str(e):
            return (
                "<div style='padding: 16px; font-family: sans-serif;'>"
                "Vitessce is available, but the optional widget dependency <b>anywidget</b> is missing "
                "for this render path. The app is using compatibility fallbacks; if you still see this message, "
                "install <code>anywidget</code> in the same Python environment used by Streamlit."
                "</div>"
            )
        import traceback
        error_msg = f"{str(e)}\n\n{traceback.format_exc()}"
        return f"<div style='color: red; padding: 20px; font-family: monospace; font-size: 12px; white-space: pre-wrap; max-height: 400px; overflow-y: auto;'>Vitessce Error: {error_msg}</div>"
    except Exception as e:
        import traceback
        error_msg = f"{str(e)}\n\n{traceback.format_exc()}"
        return f"<div style='color: red; padding: 20px; font-family: monospace; font-size: 12px; white-space: pre-wrap; max-height: 400px; overflow-y: auto;'>Vitessce Error: {error_msg}</div>"

# Spatial Analysis Function (Squidpy)
@st.cache_data
def run_spatial_analysis(_adata):
    """Run Squidpy spatial transcriptomics analysis."""
    try:
        adata_spatial = _adata.copy()

        # Check if spatial coordinates exist
        if "spatial" not in adata_spatial.obsm:
            adata_spatial.uns["spatial_analysis_status"] = "failed"
            adata_spatial.uns["spatial_analysis_message"] = "No spatial coordinates found in adata.obsm['spatial']."
            return adata_spatial

        # Compute spatial neighbors
        coord_type = "grid" if "spatial" in adata_spatial.uns else "generic"
        sq.gr.spatial_neighbors(adata_spatial, coord_type=coord_type)

        # Compute spatial autocorrelation (Moran's I)
        sq.gr.spatial_autocorr(
            adata_spatial,
            mode="moran",
            genes=adata_spatial.var_names[:100] if adata_spatial.n_vars > 100 else adata_spatial.var_names
        )

        adata_spatial.uns["spatial_analysis_status"] = "success"
        adata_spatial.uns["spatial_analysis_message"] = "Spatial neighbors and Moran's I computed successfully."

        return adata_spatial

    except Exception as e:
        adata_fallback = _adata.copy()
        adata_fallback.uns["spatial_analysis_status"] = "failed"
        adata_fallback.uns["spatial_analysis_message"] = f"Spatial analysis failed: {e}"
        return adata_fallback

# Main content area
if st.session_state.adata is not None:
    adata = st.session_state.adata

    st.subheader("Dataset Overview")

    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Total Cells", adata.n_obs)
    with col2:
        st.metric("Total Genes", adata.n_vars)
    with col3:
        st.metric("Observations", len(adata.obs.columns))

    # Display cell metadata
    st.subheader("Cell Metadata (adata.obs)")
    st.dataframe(adata.obs.head(10), width="stretch")

    # Display gene metadata
    st.subheader("Gene Metadata (adata.var)")
    st.dataframe(adata.var.head(10), width="stretch")

    # Sidebar: Run Scanpy Pipeline
    st.sidebar.header("🔬 Analysis")
    if st.sidebar.button("Run Scanpy Pipeline"):
        with st.spinner("Running preprocessing and UMAP..."):
            try:
                st.session_state.adata_processed = run_scanpy_pipeline(adata)
                st.sidebar.success("✅ Pipeline complete!")
            except Exception as e:
                st.session_state.adata_processed = None
                st.sidebar.error(f"Scanpy pipeline failed: {e}")

    # Spatial Analysis (if applicable)
    if data_type == "Spatial":
        if st.sidebar.button("Run Spatial Analysis (Squidpy)"):
            if st.session_state.adata_processed is not None:
                with st.spinner("Computing spatial neighbors and Moran's I..."):
                    st.session_state.adata_spatial = run_spatial_analysis(st.session_state.adata_processed)
                    if st.session_state.adata_spatial.uns.get("spatial_analysis_status") == "success":
                        st.sidebar.success("✅ Spatial analysis complete!")
                    else:
                        st.sidebar.error(st.session_state.adata_spatial.uns.get("spatial_analysis_message", "Spatial analysis failed."))
            else:
                st.sidebar.warning("⚠️ Please run the Scanpy pipeline first!")

    # Visualization section
    if st.session_state.adata_processed is not None:
        adata_processed = st.session_state.adata_processed

        # View selection tabs
        tab1, tab2 = st.tabs(["📊 Quick UMAP (Plotly)", "🔬 Advanced Viewer (Vitessce)"])

        with tab1:
            st.subheader("UMAP Visualization")

            # Extract UMAP coordinates and cell type info
            umap_df = pd.DataFrame(
                adata_processed.obsm["X_umap"],
                columns=["UMAP-1", "UMAP-2"],
                index=adata_processed.obs_names
            )
            umap_df["Cluster"] = adata_processed.obs["leiden"].values

            # Create interactive Plotly scatter plot
            fig = px.scatter(
                umap_df.reset_index(),
                x="UMAP-1",
                y="UMAP-2",
                color="Cluster",
                hover_name="index",
                title="UMAP: Cell Clustering",
                labels={"index": "Cell ID"},
                color_discrete_sequence=px.colors.qualitative.Light24
            )
            fig.update_layout(height=600, width=800)
            st.plotly_chart(fig, width = "stretch")

            # Expression heatmap for top genes
            st.subheader("Top Highly Variable Genes")
            required_hvg_cols = {"highly_variable", "highly_variable_rank", "means", "dispersions_norm"}
            if required_hvg_cols.issubset(set(adata_processed.var.columns)):
                top_hvg = adata_processed.var.dropna(subset=["highly_variable_rank"]).nlargest(20, "highly_variable_rank")
                st.write(f"Top 20 highly variable genes (out of {adata_processed.n_vars} total):")
                st.dataframe(top_hvg[["highly_variable", "highly_variable_rank", "means", "dispersions_norm"]], width = "stretch")
            else:
                st.info("Highly variable gene statistics are unavailable for this dataset.")

        with tab2:
            st.subheader("Vitessce Advanced Viewer")
            st.info("High-performance interactive visualization with linked views")

            with st.spinner("Loading Vitessce viewer..."):
                vitessce_html = create_vitessce_html(adata_processed)
                if isinstance(vitessce_html, str) and vitessce_html.strip().lower().startswith(("http://", "https://")):
                    components.iframe(vitessce_html.strip(), height=800)
                else:
                    components.html(vitessce_html, height=800)

        # Spatial Analysis Results
        if data_type == "Spatial" and "adata_spatial" in st.session_state and st.session_state.adata_spatial is not None:
            st.divider()
            st.subheader("🗺️ Spatial Transcriptomics Analysis")

            adata_spatial = st.session_state.adata_spatial

            status_msg = adata_spatial.uns.get("spatial_analysis_message") if hasattr(adata_spatial, "uns") else None
            if status_msg:
                if adata_spatial.uns.get("spatial_analysis_status") == "success":
                    st.success(status_msg)
                elif adata_spatial.uns.get("spatial_analysis_status") == "failed":
                    st.error(status_msg)

            # Display Moran's I results if available
            morans_df = None

            if "moranI" in adata_spatial.uns and isinstance(adata_spatial.uns["moranI"], pd.DataFrame):
                morans_df = adata_spatial.uns["moranI"].copy()
                if "I" in morans_df.columns:
                    morans_df = morans_df.rename(columns={"I": "morans_i"})
            elif "morans_i" in adata_spatial.var.columns:
                morans_df = adata_spatial.var[["morans_i"]].copy()

            if morans_df is not None and "morans_i" in morans_df.columns:
                st.subheader("Spatial Autocorrelation (Moran's I)")

                morans_df = morans_df.sort_values("morans_i", ascending=False)
                st.dataframe(morans_df.head(20), width = "stretch")

                # Plot Moran's I scores
                fig = px.bar(
                    x=morans_df.head(10)["morans_i"],
                    y=morans_df.head(10).index,
                    orientation="h",
                    title="Top 10 Spatially Variable Genes (Moran's I)",
                    labels={"x": "Moran's I Score", "y": "Gene"}
                )
                st.plotly_chart(fig, width = "stretch")

                st.info("Moran's I measures spatial autocorrelation. Higher values indicate stronger spatial clustering of gene expression.")
            else:
                # Fallback: show spatial neighbors graph if available
                if "spatial_distances" in adata_spatial.obsp:
                    st.success("✓ Spatial neighbor graph computed")
                    st.write(f"Spatial graph shape: {adata_spatial.obsp['spatial_distances'].shape}")
                else:
                    st.info("Spatial analysis not yet performed. Run 'Run Spatial Analysis (Squidpy)' button in the sidebar to compute spatial autocorrelation.")


else:
    st.info("👈 Please upload a .h5ad file or load a sample dataset from the sidebar to get started.")
