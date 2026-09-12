"""
MemVax-Py: Integrative Reverse Vaccinology Candidate Prioritization Pipeline

Research question
-----------------
Can a transparent, multi-criteria computational framework integrating
antigenicity, pathogen conservation, reduced host similarity, and
topology-derived accessibility features prioritize pathogen proteins with
characteristics associated with experimentally characterized vaccine antigens?

The pipeline evaluates protein sequences across five modules:

    1. Antigenicity
       IAPred-based intrinsic antigenicity scoring.

    2. Topology and Surface-Localization Features
       Kyte-Doolittle hydropathy-based prediction of TM-like hydrophobic
       regions and a transparent N-terminal signal-peptide heuristic.

    3. Host Similarity
       BLASTp against the human RefSeq protein database using sequence
       identity and query coverage.

    4. Pathogen Conservation
       BLASTp across predefined pathogen strain proteomes using sequence
       identity and query coverage.

    5. Usable Protein Region
       Estimates the largest contiguous non-signal, non-TM region using
       topology information already calculated in Module 2.

The final output is a ranked candidate list based on the Antigen Candidacy
Score (ACS). The pipeline is intended for candidate prioritization and
requires experimental validation.
"""

from pathlib import Path
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils.ProtParamData import kd
from Bio.Blast import NCBIXML
import sys
import io
import os
import tempfile
import pandas as pd
import subprocess


# -----------------------------------------------------------------------------
# CONFIGURATION & PATH SETUP
# -----------------------------------------------------------------------------

current_dir = Path(__file__).parent

# Path to the BLAST human protein database
HUMAN_DATABASE = current_dir / "data" / "blast_db" / "human_refseq"

# Base directory containing strain-specific proteome databases
STRAINS_DIR = current_dir / "data" / "proteomes"

# Mapping of pathogen keywords to their corresponding strain database folders
STRAIN_MAP = {
    "sars": [
        STRAINS_DIR / "SARS_CoV_2_Wuhan",
        STRAINS_DIR / "SARS_CoV_2_Delta",
    ],
    "neisseria": [
        STRAINS_DIR / "Neisseria_MC58",
        STRAINS_DIR / "Neisseria_8013",
        STRAINS_DIR / "Neisseria_alpha14",
    ],
    "bordetella": [
        STRAINS_DIR / "Bordetella_Tohama",
        STRAINS_DIR / "Bordetella_18323",
    ],
    "streptococcus": [
        STRAINS_DIR / "Streptococcus_pneumoniae_TIGR4",
        STRAINS_DIR / "Streptococcus_pneumoniae_R6",
        STRAINS_DIR / "Streptococcus_pneumoniae_D39",
    ],
    "mycobacterium": [
        STRAINS_DIR / "Mycobacterium_tuberculosis_H37Rv",
        STRAINS_DIR / "Mycobacterium_tuberculosis_CDC1551",
        STRAINS_DIR / "Mycobacterium_tuberculosis_H37Ra",
    ],
    "hepatitis": [
        STRAINS_DIR / "Hepatitis_B_genotype_D",
        STRAINS_DIR / "Hepatitis_B_genotype_C",
        STRAINS_DIR / "Hepatitis_B_genotype_B",
    ],
    "rsv": [
        STRAINS_DIR / "RSV_A2",
        STRAINS_DIR / "RSV_B1",
    ],
    "hpv": [
        STRAINS_DIR / "HPV16",
        STRAINS_DIR / "HPV18",
        STRAINS_DIR / "HPV31",
    ],
}


# =============================================================================
# MODULE 1: ANTIGENICITY
# =============================================================================

def antigenicity_score(fasta_file):
    """
    Run the external IAPred tool to compute intrinsic antigenicity.

    IAPred is executed as a subprocess. It reads the input FASTA and writes
    a single-row CSV containing the antigenicity score and category.

    Args:
        fasta_file (Path): Path to the query protein FASTA file.

    Returns:
        tuple: (score, category) where score is a float and category is a
               string (e.g., "High", "Moderate", "Low"). Returns (None, None)
               if IAPred fails or the output file is not created.
    """

    iapred_dir = current_dir / "IAPred"
    iapred_script = iapred_dir / "IApred.py"
    output_csv = current_dir / f"{fasta_file.stem}.csv"

    # Execute IAPred in its own working directory
    result = subprocess.run(
        [sys.executable, str(iapred_script), str(fasta_file), str(output_csv)],
        cwd=str(iapred_dir),
        capture_output=True,
        text=True
    )

    # Handle subprocess failure
    if result.returncode != 0:
        print("IAPred ERROR:")
        print(result.stderr)
        return None, None

    # Handle missing output file
    if not output_csv.exists():
        print("Output file was not created")
        return None, None

    # Parse the CSV output produced by IAPred
    df = pd.read_csv(output_csv)
    score = df["Intrinsic_Antigenicity_Score"].iloc[0]
    category = df["Antigenicity_Category"].iloc[0]

    return score, category


def antigenicity_normalization(score):
    """
    Linearly rescale an IAPred raw score from [-3, 3] to [0, 25].

    This maps the antigenicity metric into the pipeline's common
    0-25 scoring scale for weighted aggregation.

    Args:
        score (float): Raw IAPred intrinsic antigenicity score.

    Returns:
        float: Normalized score on the 0-25 scale.
    """

    normalized = ((score + 3) / 6) * 25

    # Keep the score within the intended 0-25 range
    normalized = max(0, min(25, normalized))

    return round(normalized, 2)


# =============================================================================
# MODULE 2: TOPOLOGY AND SURFACE-LOCALIZATION FEATURES
# =============================================================================

def get_tm_helix_positions(sequence, window_size=19, threshold=1.6, min_gap=5):
    """
    Predict TM-like hydrophobic regions using a Kyte-Doolittle
    sliding-window scan.

    This is a transparent hydropathy heuristic rather than TMHMM.

    Positive hydrophobic windows are converted into actual residue
    coordinates and overlapping or closely spaced windows are merged.

    Args:
        sequence (Seq or str): Amino acid sequence.
        window_size (int): Length of the sliding window.
        threshold (float): Kyte-Doolittle hydropathy threshold.
        min_gap (int): Maximum gap allowed when merging regions.

    Returns:
        list: Merged (start, end) residue coordinates.

              Coordinates are 0-based and inclusive.
    """

    seq_str = str(sequence)

    # Too short to contain a complete hydrophobic window
    if len(seq_str) < window_size:
        return None

    analysis = ProteinAnalysis(seq_str)
    scores = analysis.protein_scale(param_dict=kd, window=window_size)

    # Identify hydrophobic windows
    stretches = []

    for idx, score in enumerate(scores):

        if score > threshold:

            # The hydropathy score at idx represents:
            # idx ... idx + window_size - 1
            start = idx
            end = idx + window_size - 1

            stretches.append((start, end))

    if not stretches:
        return []

    # Merge overlapping or closely spaced windows
    merged = [stretches[0]]

    for start, end in stretches[1:]:

        previous_start, previous_end = merged[-1]
        gap = start - previous_end - 1

        if gap <= min_gap:
            merged[-1] = (previous_start, max(previous_end, end))
        else:
            merged.append((start, end))

    return merged


def tmhelix_count(sequence, window_size=19, threshold=1.6, min_gap=5):
    """
    Count predicted TM-like hydrophobic regions.

    This function is retained from the original pipeline.

    Args:
        sequence (Seq or str): Amino acid sequence.
        window_size (int): Length of sliding window.
        threshold (float): Hydropathy threshold.
        min_gap (int): Maximum gap allowed when merging regions.

    Returns:
        int: Number of predicted TM-like regions.
    """

    tm_positions = get_tm_helix_positions(sequence, window_size, threshold, min_gap)

    if tm_positions is None:
        return None

    return len(tm_positions)


def predict_signal_peptide(sequence, window_size=7, threshold=1.6):
    """
    Predict a possible N-terminal signal peptide using a hydropathy heuristic.

    This is NOT SignalP.

    The heuristic scans the first 40 residues for a hydrophobic core and
    then checks the downstream region for a simple cleavage-site pattern.

    Args:
        sequence (Seq or str): Amino acid sequence.
        window_size (int): Sliding window for hydrophobic core detection.
        threshold (float): Hydropathy threshold.

    Returns:
        tuple: (signal_peptide_bool, cleavage_index)

               cleavage_index is the number of residues before the
               predicted mature protein region.
    """

    seq_str = str(sequence)
    n_term = seq_str[:40]

    # Sequence too short to contain a meaningful signal peptide
    if len(n_term) < 10:
        return False, 0

    analysis = ProteinAnalysis(n_term)
    scores = analysis.protein_scale(param_dict=kd, window=window_size)

    if not scores:
        return False, 0

    max_score = max(scores)

    # No sufficiently hydrophobic core found
    if max_score <= threshold:
        return False, 0

    # Locate the strongest hydrophobic core
    core_start_index = scores.index(max_score)
    core_end_index = core_start_index + window_size

    # Inspect downstream residues for a possible cleavage site
    cleavage_zone = n_term[core_end_index:core_end_index + 10]

    for offset, aa in enumerate(cleavage_zone):

        if aa in "AGSC":
            cleavage_index = core_end_index + offset + 1
            return True, cleavage_index

    # Hydrophobic core present but no recognizable cleavage pattern
    return False, 0


def calculate_surface_score(tmhelix_result, signalp_bool):
    """
    Compute a heuristic topology/accessibility score.

    This score is NOT a direct measurement of experimentally confirmed
    surface exposure.

    The scoring favors proteins with a predicted signal peptide and
    fewer TM-like regions.

    Args:
        tmhelix_result (int): Number of predicted TM-like regions.
        signalp_bool (bool): True if a signal peptide is predicted.

    Returns:
        int: Raw topology/accessibility score.
    """

    if tmhelix_result == 0 and signalp_bool:
        return 100
    elif tmhelix_result == 1 and signalp_bool:
        return 85
    elif tmhelix_result == 2 and signalp_bool:
        return 70
    elif tmhelix_result == 3 and signalp_bool:
        return 50
    elif tmhelix_result > 3 and signalp_bool:
        return 35
    elif tmhelix_result == 0 and not signalp_bool:
        return 20
    elif tmhelix_result >= 1 and not signalp_bool:
        return 5


def surface_normalization(surface_score):
    """
    Rescale the raw topology/accessibility score to the 0-25 pipeline scale.

    Args:
        surface_score (int): Raw topology/accessibility score.

    Returns:
        float: Normalized score from 0 to 25.
    """

    normalized = surface_score / 4

    return round(max(0, min(25, normalized)), 2)


def analyze_topology(sequence, window_size=19, threshold=1.6, min_gap=5):
    """
    Perform the complete topology analysis once.

    Module 2 calculates all topology-related features that are needed
    elsewhere in the pipeline. Module 5 reuses these results instead
    of running another TM/signal-peptide prediction.

    Args:
        sequence (Seq or str): Amino acid sequence.
        window_size (int): TM hydropathy window.
        threshold (float): Hydropathy threshold.
        min_gap (int): Merge gap.

    Returns:
        dict: Topology and surface-localization features.
    """

    seq_str = str(sequence)

    # Predict TM-like regions once
    tm_positions = get_tm_helix_positions(seq_str, window_size, threshold, min_gap)

    # Calculate length of each merged region
    tm_lengths = []

    for start, end in tm_positions:
        length = end - start + 1
        tm_lengths.append(length)

    # Predict signal peptide once
    signalp_bool, cleavage_index = predict_signal_peptide(seq_str)

    # Count predicted regions
    tm_count = len(tm_positions)

    # Total hydrophobic-region length
    total_tm_length = sum(tm_lengths)

    # Fraction of sequence occupied by predicted hydrophobic regions
    if len(seq_str) > 0:
        tm_fraction = total_tm_length / len(seq_str)
    else:
        tm_fraction = 0

    # Calculate topology/accessibility score
    surface_score = calculate_surface_score(tm_count, signalp_bool)
    surface_normalized = surface_normalization(surface_score)

    return {
        "tm_positions": tm_positions,
        "tm_lengths": tm_lengths,
        "tm_count": tm_count,
        "total_tm_length": total_tm_length,
        "tm_fraction": round(tm_fraction, 4),
        "signal_peptide": signalp_bool,
        "cleavage_index": cleavage_index,
        "surface_score": surface_score,
        "surface_normalized": surface_normalized
    }


# =============================================================================
# MODULE 3: HOST SIMILARITY
# =============================================================================

def score_homology(sequence):
    """
    Calculate a host-similarity score using BLASTp against human RefSeq.

    The score considers both sequence identity and query coverage so that
    a short partial match does not receive the same penalty as a
    full-length human homolog.

    Higher scores indicate lower detectable sequence similarity to human
    proteins.

    IMPORTANT:
        This is a sequence-based host-similarity filter. It does NOT prove
        absence of autoimmune risk.

    Args:
        sequence (str): Amino acid sequence.

    Returns:
        float: Host-similarity score from 0 to 25.
               Returns None if BLAST fails.
    """

    temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False)
    temp_file.write(f'>query\n{sequence}\n')
    temp_file.close()

    temp_path = temp_file.name

    try:

        blast_output = subprocess.run(
            args=[
                'blastp',
                '-query', temp_path,
                '-db', 'human_refseq',
                '-outfmt', '5'
            ],
            cwd=str(HUMAN_DATABASE.parent),
            text=True,
            check=True,
            capture_output=True
        )

        data = blast_output.stdout
        xml_data = io.StringIO(data)
        blast_record = NCBIXML.read(xml_data)
        alignments = blast_record.alignments

        # No human hits
        if not alignments:
            return 25

        # Find significant hits
        significant_hits = []

        for alignment in alignments:

            for hsp in alignment.hsps:

                if hsp.expect <= 0.001:
                    significant_hits.append((alignment, hsp))

        # No significant human similarity detected
        if not significant_hits:
            return 25

        # Select the strongest significant HSP
        alignment, hsp = max(significant_hits, key=lambda x: getattr(x[1], "bits", 0))

        # Calculate percent identity
        percent_identity = hsp.identities / hsp.align_length

        # Calculate query coverage
        query_length = len(str(sequence))
        query_coverage = hsp.align_length / query_length
        query_coverage = min(1, query_coverage)

        # Combine identity and coverage
        similarity_burden = percent_identity * query_coverage

        # Higher similarity = lower score
        score = 25 * (1 - similarity_burden)

        # Keep score within 0-25
        score = max(0, min(25, score))

        return round(score, 2)

    except subprocess.CalledProcessError as e:

        print("BLAST ERROR:", e.stderr)
        return None

    finally:
        if os.path.exists(temp_path):
            os.remove(temp_path)


# =============================================================================
# MODULE 4: PATHOGEN CONSERVATION
# =============================================================================

def conservation_score(sequence, strain_list):
    """
    Calculate conservation across pathogen strains using BLASTp.

    For each strain:

        Significant hit
            -> identity and query coverage contribute to conservation.

        No significant hit
            -> zero contribution.

    A missing biological hit is therefore different from a failed BLAST
    calculation.

    Args:
        sequence (Seq or str): Amino acid sequence.
        strain_list (list): List of strain database folders.

    Returns:
        float: Conservation score from 0 to 25.

        None: If conservation cannot be calculated reliably or is not
              applicable.
    """

    # Conservation is not applicable when no strain databases are supplied.
    if not strain_list:
        return None

    temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False)
    temp_file.write(f'>query\n{sequence}\n')
    temp_file.close()

    temp_path = temp_file.name
    conservation_values = []

    try:

        for strain_db in strain_list:

            try:

                blast_output = subprocess.run(
                    args=[
                        'blastp',
                        '-query', temp_path,
                        '-db', strain_db.name,
                        '-outfmt', '5'
                    ],
                    cwd=str(STRAINS_DIR),
                    text=True,
                    check=True,
                    capture_output=True
                )

                data = blast_output.stdout
                xml_data = io.StringIO(data)
                blast_record = NCBIXML.read(xml_data)
                alignments = blast_record.alignments

                # No BLAST hit
                if not alignments:
                    conservation_values.append(0)
                    continue

                # Find significant hits
                significant_hits = []

                for alignment in alignments:

                    for hsp in alignment.hsps:

                        if hsp.expect <= 0.001:
                            significant_hits.append((alignment, hsp))

                # No significant hit in this strain
                if not significant_hits:
                    conservation_values.append(0)
                    continue

                # Select strongest significant HSP
                alignment, hsp = max(significant_hits, key=lambda x: getattr(x[1], "bits", 0))

                # Percent identity
                percent_identity = hsp.identities / hsp.align_length

                # Query coverage
                query_length = len(str(sequence))
                query_coverage = hsp.align_length / query_length
                query_coverage = min(1, query_coverage)

                # Conservation contribution for this strain
                strain_conservation = percent_identity * query_coverage

                conservation_values.append(strain_conservation)

            except subprocess.CalledProcessError as e:

                print(f"BLAST ERROR for {strain_db.name}:", e.stderr)

                # Technical failure is not treated as biological non-conservation
                conservation_values.append(None)

        # Keep only strains for which BLAST completed successfully
        usable_values = [value for value in conservation_values if value is not None]

        # All BLAST calculations failed
        if not usable_values:
            return None

        # Average across successfully analyzed strains
        average_conservation = sum(usable_values) / len(usable_values)

        # Convert 0-1 score to 0-25 scale
        score = average_conservation * 25
        score = max(0, min(25, score))

        return round(score, 2)

    finally:

        if os.path.exists(temp_path):
            os.remove(temp_path)


# =============================================================================
# MODULE 5: USABLE PROTEIN REGION
# =============================================================================

def find_usable_region(sequence, cleavage_index, tm_positions):
    """
    Determine the largest contiguous non-signal, non-TM region.

    Module 5 does not perform another topology prediction. It uses the
    TM-like regions already calculated by Module 2.

    The function evaluates all regions remaining after excluding the
    predicted signal peptide and TM-like regions and selects the longest
    contiguous region.

    Args:
        sequence (Seq or str): Full protein sequence.
        cleavage_index (int): Predicted signal peptide cleavage position.
        tm_positions (list): Merged TM-like regions from Module 2.

    Returns:
        tuple: (usable_start, usable_end)

               Coordinates are 0-based and inclusive.

               Returns (-1, -1) if no usable region exists.
    """

    seq_str = str(sequence)
    total_length = len(seq_str)

    if total_length == 0:
        return -1, -1

    # Start after the predicted signal peptide
    usable_start_boundary = max(0, min(cleavage_index, total_length))

    if usable_start_boundary >= total_length:
        return -1, -1

    # Keep only valid TM coordinates
    valid_tm_positions = []

    for start, end in tm_positions:

        start = max(0, min(start, total_length - 1))
        end = max(0, min(end, total_length - 1))

        if end >= start:
            valid_tm_positions.append((start, end))

    # Sort regions by sequence position
    valid_tm_positions.sort()

    # Merge overlapping TM-like regions
    merged_tm = []

    for start, end in valid_tm_positions:

        if not merged_tm:
            merged_tm.append((start, end))

        else:

            previous_start, previous_end = merged_tm[-1]

            if start <= previous_end + 1:
                merged_tm[-1] = (previous_start, max(previous_end, end))
            else:
                merged_tm.append((start, end))

    # Identify soluble/non-TM segments
    soluble_segments = []
    current_start = usable_start_boundary

    for tm_start, tm_end in merged_tm:

        # Ignore TM region completely before mature protein
        if tm_end < usable_start_boundary:
            continue

        tm_start = max(tm_start, usable_start_boundary)

        # Region before this TM-like region is usable
        if tm_start > current_start:
            soluble_segments.append((current_start, tm_start - 1))

        # Continue after this TM-like region
        current_start = max(current_start, tm_end + 1)

    # Add final C-terminal region
    if current_start <= total_length - 1:
        soluble_segments.append((current_start, total_length - 1))

    # No usable region remains
    if not soluble_segments:
        return -1, -1

    # Select the largest contiguous region
    best_region = max(soluble_segments, key=lambda x: x[1] - x[0] + 1)

    return best_region


def calculate_construct_score(usable_start, usable_end, total_length):
    """
    Calculate the fraction represented by the largest contiguous usable region.

    Args:
        usable_start (int): Start coordinate.
        usable_end (int): End coordinate.
        total_length (int): Total protein length.

    Returns:
        float: Usable fraction between 0 and 1.
    """

    if total_length <= 0:
        return 0

    # No usable region
    if usable_start < 0 or usable_end < usable_start:
        return 0

    usable_length = usable_end - usable_start + 1
    usable_fraction = usable_length / total_length

    # Keep fraction within 0-1
    usable_fraction = max(0, min(1, usable_fraction))

    return usable_fraction


def construct_normalization(fraction):
    """
    Convert the usable protein fraction into the 0-25 score scale.

    Args:
        fraction (float): Usable fraction between 0 and 1.

    Returns:
        float: Normalized score from 0 to 25.
    """

    score = fraction * 25
    score = max(0, min(25, score))

    return round(score, 2)


def calculate_module5_from_module2(sequence, topology_result):
    """
    Calculate Module 5 directly from Module 2 output.

    Module 2 predicts the topology.

    Module 5 uses those predictions to estimate the largest contiguous
    usable protein region.

    No second TM or signal-peptide prediction is performed.

    Args:
        sequence (Seq or str): Protein sequence.
        topology_result (dict): Output from analyze_topology().

    Returns:
        dict: Usable-region coordinates, length, fraction, and score.
    """

    seq_str = str(sequence)

    cleavage_index = topology_result["cleavage_index"]
    tm_positions = topology_result["tm_positions"]

    usable_start, usable_end = find_usable_region(seq_str, cleavage_index, tm_positions)

    # Calculate usable length
    if usable_start >= 0 and usable_end >= usable_start:
        usable_length = usable_end - usable_start + 1
    else:
        usable_length = 0

    # Calculate usable fraction
    fraction = calculate_construct_score(usable_start, usable_end, len(seq_str))

    # Normalize to 0-25
    normalized_score = construct_normalization(fraction)

    return {
        "usable_start": usable_start,
        "usable_end": usable_end,
        "usable_length": usable_length,
        "usable_fraction": round(fraction, 4),
        "usable_region_score": normalized_score
    }


# =============================================================================
# COMPOSITE ANTIGEN CANDIDACY SCORE
# =============================================================================

def calculate_acs(
    antigenicity_normalized,
    surface_normalized,
    host_similarity_score,
    conservation_score_value,
    usable_region_score,
    weights=None
):
    """
    Calculate the Antigen Candidacy Score (ACS).

    All five modules are normalized to the same 0-25 scale before
    weighted aggregation.

    The initial weights reflect the intended biological importance of
    each criterion within the candidate-prioritization framework:

        Antigenicity           = 30%
        Topology/Accessibility = 25%
        Conservation           = 20%
        Host Similarity        = 15%
        Usable Region          = 10%

    Antigenicity, topology/accessibility, and conservation receive greater
    weight because they are central to the biological objective of
    prioritizing promising vaccine-antigen candidates.

    Host similarity acts as a supporting host-related filter, while usable
    region is treated as a lower-weight downstream feature.

    These are model-design weights and should be evaluated through
    sensitivity and ablation analyses during validation.

    Args:
        antigenicity_normalized (float): Module 1 score, 0-25.
        surface_normalized (float): Module 2 score, 0-25.
        host_similarity_score (float): Module 3 score, 0-25.
        conservation_score_value (float): Module 4 score, 0-25.
        usable_region_score (float): Module 5 score, 0-25.
        weights (dict, optional): Custom module weights.

    Returns:
        float: ACS on a 0-100 scale.
    """

    # Initial model-design weights
    if weights is None:
        weights = {
            "antigenicity": 0.30,
            "surface": 0.25,
            "conservation": 0.20,
            "host_similarity": 0.15,
            "usable_region": 0.10
        }

    # Store module scores together
    values = {
        "antigenicity": antigenicity_normalized,
        "surface": surface_normalized,
        "conservation": conservation_score_value,
        "host_similarity": host_similarity_score,
        "usable_region": usable_region_score
    }

    # A complete ACS requires all module scores
    if any(value is None for value in values.values()):
        return None

    # Check weight total
    weight_sum = sum(weights.values())

    if weight_sum <= 0:
        raise ValueError("ACS weights must have a positive sum.")

    # Weighted score on the 0-25 scale
    weighted_score = 0

    for module_name, value in values.items():
        normalized_weight = weights[module_name] / weight_sum
        weighted_score += float(value) * normalized_weight

    # Convert 0-25 scale to 0-100
    acs = weighted_score * 4
    acs = max(0, min(100, acs))

    return round(acs, 2)