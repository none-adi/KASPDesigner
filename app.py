# app.py – KASP Primer Designer (visible larger text, matching footer background)
import streamlit as st
import pandas as pd
import concurrent.futures
from io import BytesIO
from kasp_core import design_one_locus

# ---------------------------------------------------------------------------
# Page configuration
# ---------------------------------------------------------------------------
st.set_page_config(
    page_title="KASP Primer Designer",
    page_icon="🧬",
    layout="centered",
    initial_sidebar_state="collapsed",
)

# ---------------------------------------------------------------------------
# Custom CSS – high specificity to override Streamlit defaults
# ---------------------------------------------------------------------------
st.markdown("""
<style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap');

    /* Set the entire app background */
    html, body, .stApp, .main, .block-container {
        background-color: #F8FAFC !important;
        font-family: 'Inter', sans-serif !important;
    }

    /* Increase base font size for all text */
    .stApp {
        font-size: 1.05rem !important;
    }

    /* Main title – using h1 for maximum size, override Streamlit's h1 style */
    h1.main-title {
        font-size: 3.8rem !important;
        font-weight: 700 !important;
        color: #1E293B !important;
        margin-bottom: 0.2rem !important;
        line-height: 1.2 !important;
    }

    /* Subtitle */
    .subtitle {
        font-size: 1.2rem !important;
        color: #64748B !important;
        margin-bottom: 2rem !important;
        line-height: 1.5 !important;
    }

    /* Section labels */
    .section-label {
        font-size: 1.3rem !important;
        font-weight: 600 !important;
        color: #334155 !important;
        margin-top: 2rem !important;
        margin-bottom: 0.5rem !important;
    }

    /* Descriptive text */
    .desc-text {
        font-size: 1.05rem !important;
        color: #475569 !important;
        line-height: 1.6 !important;
    }

    /* Soft blue divider (gradient) */
    .soft-divider {
        height: 1px;
        background: linear-gradient(to right, #93B4D7, transparent) !important;
        margin: 1.8rem 0 !important;
    }

    /* Fixed footer – background matches page, no white */
    .fixed-footer {
        position: fixed;
        bottom: 0;
        left: 0;
        right: 0;
        background-color: #F8FAFC !important;
        border-top: 1px solid #E2E8F0 !important;
        padding: 0.8rem 1rem;
        text-align: center;
        font-size: 0.95rem !important;
        color: #94A3B8 !important;
        z-index: 1000;
    }

    /* Enlarge buttons */
    .stButton>button, .stDownloadButton>button {
        font-size: 1rem !important;
        padding: 0.5rem 1.2rem !important;
    }

    /* Make the upload zone border consistent */
    .upload-zone {
        border: 2px dashed #CBD5E1 !important;
        border-radius: 8px !important;
        padding: 1rem 1.5rem !important;
        background-color: #FFFFFF !important;
    }
</style>
""", unsafe_allow_html=True)

# ---------------------------------------------------------------------------
# Helper – sample CSV
# ---------------------------------------------------------------------------
def get_sample_csv_bytes() -> BytesIO:
    sample = (
        "rs123,A,G,ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT[A/G]ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n"
        "rs456,T,C,TAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC[T/C]GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA\n"
    )
    return BytesIO(sample.encode("utf-8"))

# ===========================================================================
# PAGE CONTENT
# ===========================================================================

# --- Header (h1 with class) ---
st.markdown('<h1 class="main-title">🧬 KASP Primer Designer</h1>', unsafe_allow_html=True)
st.markdown(
    '<p class="subtitle">Design allele‑specific & common primers for high‑throughput genotyping — upload, run, download.</p>',
    unsafe_allow_html=True,
)

# --- Brief about & sample download ---
col1, col2 = st.columns([3, 1])
with col1:
    st.markdown(
        """
        <p class="desc-text">
        Implements the <strong>EasyKASP</strong> algorithm (BMC Bioinformatics 2025) with nearest‑neighbour Tm,
        primer‑dimer/hairpin checks, and multi‑threaded processing.
        </p>
        """,
        unsafe_allow_html=True,
    )
with col2:
    st.download_button(
        label="📄 Get sample CSV",
        data=get_sample_csv_bytes(),
        file_name="sample_loci.csv",
        mime="text/csv",
    )

st.markdown('<div class="soft-divider"></div>', unsafe_allow_html=True)

# --- Threads ---
st.markdown('<p class="section-label">⚙️ Processing Threads</p>', unsafe_allow_html=True)
workers = st.slider(
    "Number of threads to use",
    min_value=1,
    max_value=8,
    value=4,
    help="More threads process large files faster, but consume more CPU.",
)

st.markdown('<div class="soft-divider"></div>', unsafe_allow_html=True)

# --- Upload ---
st.markdown('<p class="section-label">📁 Upload Your CSV</p>', unsafe_allow_html=True)
st.markdown(
    """
    <p class="desc-text">
    Format: <code>Locus, Allele1, Allele2, Sequence</code> (no header). Variant inside <code>[ ]</code> with <code>/</code> (e.g., ...A/G...). Flanks ≥40 bp.
    </p>
    """,
    unsafe_allow_html=True,
)

uploaded_file = st.file_uploader("", type=["csv"], label_visibility="collapsed")

df_input = None
if uploaded_file is not None:
    try:
        df_input = pd.read_csv(uploaded_file, header=None)
        if df_input.shape[1] < 4:
            st.error("❌ CSV must have exactly 4 columns.")
            df_input = None
        else:
            st.success(f"✅ {len(df_input)} loci loaded.")
    except Exception as e:
        st.error(f"Error reading file: {e}")
        df_input = None

st.markdown('<div class="soft-divider"></div>', unsafe_allow_html=True)

# --- Design button & progress ---
if df_input is not None:
    st.markdown('<p class="section-label">🚀 Design Primers</p>', unsafe_allow_html=True)
    if st.button("Start Design", type="primary"):
        rows = [(i, row.tolist()) for i, row in df_input.iterrows()]
        progress_bar = st.progress(0)
        status_text = st.empty()

        results = []
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
            futures = {executor.submit(design_one_locus, row): row for row in rows}
            total = len(rows)
            completed = 0
            for future in concurrent.futures.as_completed(futures):
                try:
                    res = future.result()
                    results.append(res)
                except Exception as e:
                    row = futures[future]
                    idx = row[0]
                    results.append({
                        "locus": f"Row_{idx}",
                        "fam": "",
                        "hex": "",
                        "common": "",
                        "status": "Failed",
                        "reason": f"Exception: {str(e)}",
                        "_index": idx,
                    })
                completed += 1
                progress_bar.progress(completed / total)
                status_text.text(f"Processing {completed} of {total}…")

        # Sort by original index
        results.sort(key=lambda r: r["_index"])
        for r in results:
            del r["_index"]

        df_out = pd.DataFrame(results)
        success_count = (df_out["status"] == "Success").sum()
        fail_count = len(df_out) - success_count

        c1, c2, c3 = st.columns(3)
        c1.metric("🧬 Total loci", len(df_out))
        c2.metric("✅ Success", success_count)
        c3.metric("❌ Failed", fail_count)

        st.session_state["results_df"] = df_out
        st.session_state["results_csv"] = df_out.to_csv(index=False)

        if fail_count > 0:
            with st.expander(f"⚠️ {fail_count} failures – click for details"):
                st.dataframe(df_out[df_out["status"] == "Failed"][["locus", "reason"]], use_container_width=True)

        progress_bar.empty()
        status_text.empty()
        st.success("Design complete.")

    # --- Results ---
    if "results_df" in st.session_state:
        st.markdown('<p class="section-label">📊 Results</p>', unsafe_allow_html=True)
        st.dataframe(st.session_state["results_df"], use_container_width=True)
        st.download_button(
            label="📥 Download CSV",
            data=st.session_state["results_csv"],
            file_name="kasp_primers.csv",
            mime="text/csv",
        )

# ---------------------------------------------------------------------------
# Fixed footer – background matches page, no white strip
# ---------------------------------------------------------------------------
st.markdown(
    """
    <div class="fixed-footer">
        Built with ❤️ using Streamlit · EasyKASP methodology (BMC Bioinformatics 2025)
    </div>
    """,
    unsafe_allow_html=True,
)