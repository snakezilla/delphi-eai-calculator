"""
Delphi EAI Calculator - Streamlit Web Application

A simple, user-friendly interface for calculating SET2,3 scores
from XPonent and Accession files.

Usage:
    streamlit run app.py
"""

import streamlit as st
import pandas as pd
import io
from datetime import datetime

from parsers import parse_xponent_csv, parse_accession_xlsx, XPonentData, AccessionData
from algorithm import process_batch, results_to_dataframe, ALGO_VERSION

# Constants
MAX_FILE_SIZE_MB = 50
MAX_FILE_SIZE_BYTES = MAX_FILE_SIZE_MB * 1024 * 1024

# Page configuration
st.set_page_config(
    page_title="Delphi EAI Calculator",
    page_icon="🧬",
    layout="wide"
)


def main():
    # Header
    st.title("🧬 Delphi EAI Calculator")
    st.markdown(f"**SET2,3 Score Calculator for HR+/HER2- Breast Cancer** | Algorithm v{ALGO_VERSION}")
    st.markdown("---")

    # Instructions
    with st.expander("📖 Instructions (Click to expand)", expanded=False):
        st.markdown("""
        ### How to Use This Tool

        1. **Upload your XPonent file** (.csv)
           - This is the export from your Luminex instrument
           - Must contain the "Net MFI" data section
           - Should include negative control samples (labeled "nc")

        2. **Upload your Accession file** (.xlsx)
           - Contains patient/specimen information
           - Must have columns for: Sample ID, Sample Type, Tumor Stage, Tumor Node
           - For pathological samples: also needs Tumor Size and Positive Lymph Nodes

        3. **Click "Calculate SET2,3 Scores"**

        4. **Review results and download**

        ### Required File Formats

        **XPonent CSV:**
        - Standard Luminex xPONENT export format
        - Must contain "Net MFI" section with gene expression data

        **Accession XLSX:**
        - Must have "Sample Id" column matching XPonent sample names
        - Must have "Clinical or Pathology" column
        - Must have "Tumor Stage" column (cT0-cT4 or pT0-pT4)
        - Must have "Tumor Node" column (cN0-cN3 or pN0-pN3)
        - For pathological: "Tumor Size" and "Number of Positive Lymph Nodes"
        """)

    # File upload section
    st.header("📁 Upload Files")

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("XPonent File")
        xponent_file = st.file_uploader(
            "Upload XPonent CSV file",
            type=['csv'],
            help="The CSV export from your Luminex instrument containing Net MFI data"
        )

    with col2:
        st.subheader("Accession File")
        accession_file = st.file_uploader(
            "Upload Accession XLSX file",
            type=['xlsx', 'xls'],
            help="The specimen information file with sample type and staging data"
        )

    st.markdown("---")

    # Process button
    if xponent_file is not None and accession_file is not None:
        if st.button("🔬 Calculate SET2,3 Scores", type="primary", use_container_width=True):
            process_files(xponent_file, accession_file)
    elif xponent_file is not None or accession_file is not None:
        st.warning("⚠️ Please upload both files to continue.")
    else:
        st.info("👆 Upload your XPonent CSV and Accession XLSX files above to get started.")


def process_files(xponent_file, accession_file):
    """Process uploaded files and display results."""

    # File size validation
    if xponent_file.size > MAX_FILE_SIZE_BYTES:
        st.error(f"❌ XPonent file too large ({xponent_file.size / 1024 / 1024:.1f}MB). "
                f"Maximum allowed size is {MAX_FILE_SIZE_MB}MB.")
        return

    if accession_file.size > MAX_FILE_SIZE_BYTES:
        st.error(f"❌ Accession file too large ({accession_file.size / 1024 / 1024:.1f}MB). "
                f"Maximum allowed size is {MAX_FILE_SIZE_MB}MB.")
        return

    # Progress indicator
    progress_bar = st.progress(0)
    status_text = st.empty()

    try:
        # Step 1: Parse XPonent file
        status_text.text("📊 Parsing XPonent file...")
        progress_bar.progress(10)

        xponent_content = xponent_file.read()
        xponent_data = parse_xponent_csv(xponent_content)

        st.success(f"✅ XPonent file loaded: {len(xponent_data.samples)} samples, "
                   f"{len(xponent_data.negative_controls)} negative control(s)")

        progress_bar.progress(30)

        # Step 2: Parse Accession file
        status_text.text("📋 Parsing Accession file...")

        accession_content = accession_file.read()
        accession_data = parse_accession_xlsx(accession_content)

        st.success(f"✅ Accession file loaded: {len(accession_data.samples)} records")

        progress_bar.progress(50)

        # Step 3: Match samples
        status_text.text("🔗 Matching samples...")

        xponent_samples = set(xponent_data.samples) - set(xponent_data.negative_controls)
        accession_samples = set(accession_data.samples['sample_id'].tolist())

        matched = xponent_samples & accession_samples
        unmatched_xponent = xponent_samples - accession_samples
        unmatched_accession = accession_samples - xponent_samples

        if len(matched) == 0:
            st.error("❌ No matching samples found between files! Please check that Sample IDs match.")
            st.write("**XPonent samples:**", list(xponent_samples)[:10])
            st.write("**Accession samples:**", list(accession_samples)[:10])
            return

        if len(unmatched_xponent) > 0:
            st.warning(f"⚠️ {len(unmatched_xponent)} samples in XPonent not found in Accession: "
                      f"{list(unmatched_xponent)[:5]}{'...' if len(unmatched_xponent) > 5 else ''}")

        progress_bar.progress(60)

        # Step 4: Run algorithm
        status_text.text("🧮 Calculating SET2,3 scores...")

        results = process_batch(
            net_mfi_df=xponent_data.net_mfi,
            counts_df=xponent_data.counts,
            accession_df=accession_data.samples,
            nc_sample_ids=xponent_data.negative_controls
        )

        progress_bar.progress(90)

        # Step 5: Format results
        status_text.text("📊 Formatting results...")

        results_df = results_to_dataframe(results)

        progress_bar.progress(100)
        status_text.text("✅ Complete!")

        # Display results
        st.markdown("---")
        st.header("📊 Results")

        # Summary statistics
        col1, col2, col3, col4 = st.columns(4)

        with col1:
            st.metric("Total Samples", len(results))

        with col2:
            qc_passed = sum(1 for r in results if r.qc_passed)
            st.metric("QC Passed", f"{qc_passed}/{len(results)}")

        with col3:
            high_risk = sum(1 for r in results if r.set23_category == "High")
            st.metric("High Risk", high_risk)

        with col4:
            low_risk = sum(1 for r in results if r.set23_category == "Low")
            st.metric("Low Risk", low_risk)

        st.markdown("---")

        # Results table
        st.subheader("📋 Full Results Table")

        # Display key columns first
        display_cols = [
            'Sample', 'PASSED (Y/N)', 'Sample Type', 'SET2,3', 'SET2,3 Category',
            'SET ER/PR', 'ESR1', 'ERBB2', 'PGR', 'AURKA',
            'RNA4', 'Stage Risk', 'Node Risk', 'BPI'
        ]

        # Filter to columns that exist
        display_cols = [c for c in display_cols if c in results_df.columns]

        st.dataframe(
            results_df[display_cols],
            use_container_width=True,
            hide_index=True
        )

        # Show QC failures
        qc_failures = results_df[results_df['PASSED (Y/N)'] == 'FAIL']
        if len(qc_failures) > 0:
            st.subheader("⚠️ QC Failures")
            st.dataframe(
                qc_failures[['Sample', 'PASSED (Y/N)', 'Mean Ref Genes', 'Bead Count', 'QC Fail Reason']],
                use_container_width=True,
                hide_index=True
            )

        st.markdown("---")

        # Download section
        st.subheader("📥 Download Results")

        # Prepare CSV download
        csv_buffer = io.StringIO()
        results_df.to_csv(csv_buffer, index=False)
        csv_data = csv_buffer.getvalue()

        # Prepare Excel download
        excel_buffer = io.BytesIO()
        with pd.ExcelWriter(excel_buffer, engine='openpyxl') as writer:
            results_df.to_excel(writer, sheet_name='Results', index=False)

            # Add summary sheet
            summary_data = {
                'Metric': ['Total Samples', 'QC Passed', 'QC Failed', 'High Risk', 'Low Risk', 'Algorithm Version'],
                'Value': [
                    len(results),
                    qc_passed,
                    len(results) - qc_passed,
                    high_risk,
                    low_risk,
                    ALGO_VERSION
                ]
            }
            pd.DataFrame(summary_data).to_excel(writer, sheet_name='Summary', index=False)

        excel_data = excel_buffer.getvalue()

        col1, col2 = st.columns(2)

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

        with col1:
            st.download_button(
                label="📄 Download CSV",
                data=csv_data,
                file_name=f"EAI_Results_{timestamp}.csv",
                mime="text/csv",
                use_container_width=True
            )

        with col2:
            st.download_button(
                label="📊 Download Excel",
                data=excel_data,
                file_name=f"EAI_Results_{timestamp}.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                use_container_width=True
            )

        # Show full details expander
        with st.expander("🔍 View All Columns"):
            st.dataframe(results_df, use_container_width=True, hide_index=True)

    except ValueError as e:
        st.error(f"❌ File format error: {str(e)}")
        st.info("Please check that you've uploaded the correct file types and formats.")

    except Exception as e:
        st.error(f"❌ Error processing files: {str(e)}")
        with st.expander("🔍 Technical details (for support)"):
            st.exception(e)

    finally:
        progress_bar.empty()
        status_text.empty()


# Footer
def show_footer():
    st.markdown("---")
    st.markdown(
        f"""
        <div style="text-align: center; color: #888; font-size: 0.8em;">
            Delphi EAI Calculator v{ALGO_VERSION} |
            For research and clinical use |
            © 2026 Delphi Diagnostics
        </div>
        """,
        unsafe_allow_html=True
    )


if __name__ == "__main__":
    main()
    show_footer()
