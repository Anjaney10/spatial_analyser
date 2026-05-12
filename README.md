# Spatial-SC Explorer 🧬

**Spatial-SC Explorer** is a high-performance, interactive web application built for the rapid analysis and visualization of single-cell and spatial transcriptomics data.

By bridging the gap between the **scverse** (Scanpy/Squidpy) backend and the **Vitessce** high-fidelity visualization framework, this tool allows researchers to move seamlessly from raw `.h5ad` files to interactive "brushing and linking" discovery—all within a local, offline-capable Streamlit environment.

---

## 🚀 Key Features

* **Interactive Brushing & Linking:** Seamlessly select clusters in UMAP/Latent space and instantly see their physical distribution in Tissue space.
* **Integrated scverse Pipeline:** Run standard preprocessing (Normalization, Log1p, PCA, Neighbors, UMAP, and Leiden clustering) directly through the UI.
* **Spatial Autocorrelation:** Calculate and visualize spatial variable genes using Squidpy (Moran's I).
* **Advanced WebGL Rendering:** Support for large datasets (100k+ cells) using Vitessce's GPU-accelerated viewer.
* **Local & Offline:** Designed to run on your local machine using your system's resources, ensuring data privacy and offline accessibility.

---

## 🛠️ Tech Stack

* **Language:** Python 3.11 (Recommended)
* **Analysis:** `scanpy`, `squidpy`, `anndata`
* **Frontend:** `streamlit`
* **Visualization:** `vitessce`, `plotly`
* **Environment:** `mamba` / `miniforge`

---

## 📦 Installation

Due to the complex C-dependencies in bioinformatics libraries and recent breaking changes in `pandas 3.0` and `zarr 3.0`, please follow these steps precisely to ensure a stable build.

### 1. Prerequisites

We recommend using **Miniforge** (which includes `mamba`) for faster and more reliable dependency resolution.

### 2. Create the Environment

```bash
# Create a dedicated Python 3.11 environment
mamba create -n spatial_sc python=3.11 -y
mamba activate spatial_sc

```

### 3. Install Core Dependencies

We install the heavy-lifting libraries via `mamba` to ensure pre-compiled binary compatibility.

```bash
mamba install "zarr<3.0.0" "pandas<3.0.0" "squidpy<1.8.0" "python-igraph" "leidenalg" scanpy anndata plotly -c conda-forge -y

```

### 4. Install UI Components

```bash
pip install streamlit vitessce

```

---

## 🖥️ How to Use

1. **Clone the Repository:**
```bash
git clone https://github.com/yourusername/spatial-sc-explorer.git
cd spatial-sc-explorer

```


2. **Launch the App:**

```bash
    streamlit run app.py
    ```

3.  **Analyze Your Data:**
    *   **Upload:** Drag and drop your `.h5ad` file into the sidebar.
    *   **Process:** Click **"Run Pipeline"** to compute embeddings and clusters.
    *   **Explore:** Toggle between Plotly and Vitessce viewers to interrogate your tissue samples.

---

## ⚠️ Troubleshooting & Known Issues

### Python Version
**Do not use Python 3.12.** Several core dependencies (like `numba` and `vitessce`) currently rely on modules (e.g., `distutils`, `pkgutil.ImpImporter`) that were removed in 3.12. Stick to **3.11**.

### Zarr Conflict
The app uses a "manual export" bridge to handle Vitessce visualizations. If you see Zarr-related import errors, ensure you have forced the Zarr 2.x downgrade:
`mamba install "zarr<3.0.0"`

### Memory Management
Single-cell data is RAM-intensive. The app uses `@st.cache_resource` to keep data in memory. If the app crashes on a laptop with <16GB RAM when using large datasets, try subsetting your data before uploading.

---

## 🗺️ Roadmap

- [ ] **WebAssembly Deployment:** Porting the logic to `stlite` for zero-install browser execution.
- [ ] **Multi-modal Support:** Adding support for H&E image overlays in the spatial view.
- [ ] **Differential Expression:** Integrated UI for volcano plots and marker gene identification.

---

## 📄 License
Distributed under the MIT License. See `LICENSE` for more information.

> **Note:** This tool is a Proof of Concept (PoC) intended for research purposes. It is designed to demonstrate the power of interactive spatial transcriptomics and is not a substitute for a full production bioinformatics platform.

```
