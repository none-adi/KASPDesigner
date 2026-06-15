import streamlit as st
import pandas as pd
from kasp_core import process_rows

st.set_page_config(page_title="KASP Primer Designer", layout="wide")
st.title("🧬 KASP Primer Design Tool")
st.markdown("Upload a CSV file with your loci (four columns: Locus, Allele1, Allele2, Sequence).")

uploaded_file = st.file_uploader("Choose a CSV file", type=["csv"])

if uploaded_file is not None:
    # Read the uploaded CSV (no header)
    df = pd.read_csv(uploaded_file, header=None)
    if df.shape[1] < 4:
        st.error("The CSV must have at least 4 columns: LocusName, Allele1, Allele2, Sequence.")
    else:
        st.success(f"Loaded {len(df)} loci.")
        # Convert DataFrame to list of (index, row) tuples
        rows = [(i, row.tolist()) for i, row in df.iterrows()]

        if st.button("🚀 Design Primers"):
            with st.spinner("Designing primers... This may take a moment."):
                progress_bar = st.progress(0)
                # process_rows handles threading and returns sorted results
                # We'll manually update progress because process_rows is synchronous.
                # A better approach: run process_rows in a thread and poll? 
                # Simpler: we trust the user waits. For better UX, we could
                # iterate with ThreadPoolExecutor ourselves and update progress.
                # Let's do that here for a smoother experience.

                import concurrent.futures
                from kasp_core import design_one_locus

                results = []
                max_workers = 4
                with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
                    futures = {executor.submit(design_one_locus, row): row for row in rows}
                    for i, future in enumerate(concurrent.futures.as_completed(futures)):
                        res = future.result()
                        results.append(res)
                        progress_bar.progress((i+1)/len(rows))

                # Sort by original index
                results.sort(key=lambda r: r["_index"])
                # Remove the internal index before display
                for r in results:
                    del r["_index"]

                df_out = pd.DataFrame(results)
                st.subheader("Results")
                st.dataframe(df_out)

                # Download button
                csv = df_out.to_csv(index=False)
                st.download_button("📥 Download results as CSV", csv, "primers.csv", "text/csv")