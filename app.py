# ============ Import Libraries ==============
import streamlit as st
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from memvax_core import *
import tempfile
import os
import textwrap


# ========== Display Formatting Helpers ===========
def format_regions_for_display(regions):
    """
    Convert internal TM-region list into a readable string for
    table/CSV display. Does not affect scoring logic.
    """
    if regions is None:
        return "Not analyzed"
    if len(regions) == 0:
        return "None detected"
    return "; ".join(f"{start}-{end}" for start, end in regions)


def format_lengths_for_display(lengths):
    """
    Convert internal TM-length list into a readable string for
    table/CSV display. Does not affect scoring logic.
    """
    if lengths is None:
        return "Not analyzed"
    if len(lengths) == 0:
        return "None detected"
    return ", ".join(f"{length} aa" for length in lengths)


# ============ Page Setup ==============
st.set_page_config(page_title="MemVax-Py", layout="wide")

st.markdown(
    textwrap.dedent("""
    <style>
    html, body, [data-testid="stAppViewContainer"] {
        height: 100%;
    }

    [data-testid="stAppViewContainer"] > .main {
        display: flex;
        flex-direction: column;
        min-height: 100vh;
    }

    .memvax-footer {
        margin-top: auto;
    }

    .stMarkdown, .stCaption, p {
        font-size: 1.05rem !important;
    }
    </style>
    """),
    unsafe_allow_html=True
)

st.title("MemVax-Py")

st.markdown(
    "**An Integrative Computational Pipeline for Vaccine Antigen Candidate Prioritization Across Diverse Pathogens**"
)


# ============ List of organisms that tool supports =============
SUPPORTED_ORGANISMS = {
    "SARS-CoV-2": "sars",
    "Neisseria meningitidis": "neisseria",
    "Bordetella pertussis": "bordetella",
    "Streptococcus pneumoniae": "streptococcus",
    "Mycobacterium tuberculosis": "mycobacterium",
    "Hepatitis B virus": "hepatitis",
    "Respiratory syncytial virus": "rsv",
    "Human papillomavirus": "hpv",
}


# ============= SIDEBAR =============
st.sidebar.header("Upload & Settings")
organism = st.sidebar.selectbox("Select organism",list(SUPPORTED_ORGANISMS.keys()))
uploaded_files = st.sidebar.file_uploader("Upload FASTA file(s)", type=["fasta", "fa"], accept_multiple_files=True)
run_btn = st.sidebar.button("Run Analysis", type="primary")


# ============= MAIN AREA =============
st.divider()
col1, col2, col3, col4, col5 = st.columns(5)

with col1:
    st.markdown("**1. Antigenicity**")
    st.caption(
        "Estimates intrinsic antigenicity using the IAPred model."
    )

with col2:
    st.markdown("**2. Topology/Accessibility**")
    st.caption(
        "Uses hydropathy and signal-peptide features to estimate topology-related accessibility."
    )

with col3:
    st.markdown("**3. Host Similarity**")
    st.caption(
        "Screens for detectable sequence similarity to human proteins using BLASTp."
    )

with col4:
    st.markdown("**4. Conservation**")
    st.caption(
        "Assesses sequence conservation across the selected pathogen strains."
    )

with col5:
    st.markdown("**5. Usable Protein Region**")
    st.caption(
        "Estimates the largest contiguous region remaining outside predicted signal and TM-like regions."
    )

st.divider()

st.markdown(
    "Upload one or more FASTA protein sequences from the sidebar, "
    "select the source organism, and click **Run Analysis** to calculate "
    "the composite Antigen Candidacy Score (ACS) for each protein."
)


# =============================================================================
# RUN ANALYSIS
# =============================================================================

if run_btn and uploaded_files:
    all_results = []
    with st.spinner("Running analysis..."):
        for uploaded_file in uploaded_files:

            # Save uploaded FASTA file to a temporary location
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp:
                tmp.write(uploaded_file.getvalue().decode('utf-8'))
                tmp_path = Path(tmp.name)

            try:
                # Read the first FASTA sequence
                record = next(SeqIO.parse(tmp_path, "fasta"))
                protein_name = Path(uploaded_file.name).stem.replace("_", " ")
                sequence = record.seq

                # =================================================================
                # MODULE 1: ANTIGENICITY
                # =================================================================

                score, category = antigenicity_score(tmp_path)
                if score is None:
                    st.warning(f"Skipping {protein_name}:"
                        "IAPred failed. Check terminal for details."
                    )
                    continue

                antigenic_norm = antigenicity_normalization(score)


                # =================================================================
                # MODULE 2: TOPOLOGY / ACCESSIBILITY
                # =================================================================

                # Run topology analysis ONCE.
                topology = analyze_topology(sequence)

                signalp_bool = topology["signal_peptide"]
                cleavage_index = topology["cleavage_index"]

                tm_positions = topology["tm_positions"]
                tm_lengths = topology["tm_lengths"]
                tmhelix_result = topology["tm_count"]

                surface_score = topology["surface_score"]
                surface_norm = topology["surface_normalized"]


                # =================================================================
                # MODULE 3: HOST SIMILARITY
                # =================================================================

                homology_score = score_homology(str(sequence))

                if homology_score is None:
                    st.warning(
                        f"Skipping {protein_name}: "
                        "human BLASTp analysis failed."
                    )
                    continue


                # =================================================================
                # MODULE 4: PATHOGEN CONSERVATION
                # =================================================================

                keyword = SUPPORTED_ORGANISMS[organism]
                strain_list = STRAIN_MAP.get( keyword, [])
                conservation = conservation_score(sequence,strain_list)

                if conservation is None:
                    st.warning(
                        f"Skipping {protein_name}: "
                        "conservation analysis could not be completed."
                    )
                    continue


                # =================================================================
                # MODULE 5: USABLE PROTEIN REGION
                # =================================================================

                # Module 5 reuses Module 2 output.
                usable_region = calculate_module5_from_module2(sequence, topology)
                usable_start = usable_region["usable_start"]
                usable_end = usable_region["usable_end"]
                usable_length = usable_region["usable_length"]
                usable_fraction = usable_region["usable_fraction"]
                usable_score = usable_region["usable_region_score"]


                # =================================================================
                # COMPOSITE ANTIGEN CANDIDACY SCORE
                # =================================================================

                acs = calculate_acs(
                    antigenic_norm,
                    surface_norm,
                    homology_score,
                    conservation,
                    usable_score
                )


                # =================================================================
                # SAVE RESULT
                # =================================================================

                all_results.append({
                    "Protein": protein_name,
                    "Antigenicity": round(float(score), 3),
                    "Antigenicity Category": category,
                    "Antigenicity Score (0-25)": antigenic_norm,
                    "Signal Peptide": signalp_bool,
                    "Cleavage Index": cleavage_index,
                    "TM-like Count": tmhelix_result,
                    "TM-like Regions": format_regions_for_display(tm_positions),
                    "TM-like Lengths": format_lengths_for_display(tm_lengths),
                    "Topology / Accessibility Score": surface_norm,
                    "Host Similarity Score": homology_score,
                    "Conservation Score": conservation,
                    "Usable Start": usable_start,
                    "Usable End": usable_end,
                    "Usable Length": usable_length,
                    "Usable Fraction": usable_fraction,
                    "Usable Region Score": usable_score,
                    "ACS": acs,
                })


            except StopIteration:
                st.warning(
                    f"Skipping {protein_name}: "
                    "no valid FASTA sequence was found."
                )

            except Exception as e:
                st.warning(
                    f"Error while analyzing {protein_name}: {e}"
                )

            finally:

                # Remove temporary FASTA file
                if os.path.exists(tmp_path):
                    os.remove(tmp_path)


    # =========================================================================
    # SHOW RESULTS
    # =========================================================================

    if all_results:
        df = pd.DataFrame(all_results)
        # Rank by final ACS
        df = df.sort_values("ACS", ascending=False).reset_index(drop=True)
        df["Rank"] = df.index + 1
        st.success(f"Done! {len(df)} proteins analyzed.")
        st.markdown("### Antigen Candidacy Ranking")
        st.dataframe(df, use_container_width=True, hide_index=True)


        # =====================================================================
        # SCORE INTERPRETATION
        # =====================================================================

        st.markdown("### Score interpretation")
        st.caption(
            "The ACS integrates five computational criteria on a common "
            "0-25 scale using the initial model weights: Antigenicity 30%, "
            "Topology/Accessibility 25%, Conservation 20%, Host Similarity 15%, "
            "and Usable Protein Region 10%. The weights represent the design "
            "logic of the initial prioritization model and should be evaluated "
            "through sensitivity and ablation analysis."
        )


        # =====================================================================
        # DOWNLOAD
        # =====================================================================

        csv = df.to_csv(index=False)
        st.download_button("Download CSV", csv, "memvax_results.csv", "text/csv")

    else:
        st.error("No valid results.")


elif run_btn and not uploaded_files:
    st.error("Please upload at least one FASTA file.")


# =============================================================================
# FOOTER
# =============================================================================

st.markdown(
    '<div class="memvax-footer" style="background-color:#161B26;border-top:2px solid #444444;'
    'margin-top:3rem;padding:2rem 1rem;text-align:center;color:#B0B0B0;">'
    '<p style="font-size:1.6rem;font-weight:700;color:#FFFFFF;margin-bottom:0.4rem;">MemVax-Py</p>'
    '<p style="font-size:1.3rem;color:#D0D0D0;margin-bottom:0.8rem;">'
    'An Integrative Computational Pipeline for Vaccine Antigen Candidate Prioritization<br>'
    'Across Diverse Pathogens</p>'
    '<p style="font-size:1.1rem;color:#999999;margin-bottom:1rem;">Developed by Sana Aziz Sial • Version 1.0</p>'
    '<p style="font-size:1.1rem;color:#4A9D7F;font-weight:600;margin-bottom:1rem;">'
    'Antigenicity • Topology / Accessibility • Host Similarity • Conservation • Usable Protein Region</p>'
    '<p style="font-size:1rem;font-style:italic;color:#999999;margin-bottom:0.5rem;">'
    'Computational predictions are intended for research and require experimental validation.</p>'
    '<p style="font-size:1rem;color:#777777;">© 2026 MemVax-Py</p>'
    '</div>',
    unsafe_allow_html=True
)