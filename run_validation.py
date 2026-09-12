import os
import sys
import tempfile
from pathlib import Path
import pandas as pd
from Bio import SeqIO

# Import your core pipeline functions from memvax_core.py
from memvax_core import (
    antigenicity_score,
    antigenicity_normalization,
    analyze_topology,
    score_homology,
    conservation_score,
    calculate_module5_from_module2,
    calculate_acs,
    current_dir,
    STRAIN_MAP
)

# Map the exact pathogen names from your CSV to the keys in your STRAIN_MAP
PATHOGEN_TO_STRAIN_KEY = {
    "Bordetella pertussis": "bordetella",
    "Hepatitis B virus": "hepatitis",
    "Human papillomavirus": "hpv",
    "Mycobacterium tuberculosis": "mycobacterium",
    "Neisseria meningitidis": "neisseria",
    "Respiratory syncytial virus": "rsv",
    "SARS-CoV-2": "sars",
    "Streptococcus pneumoniae": "streptococcus"
}

def run_validation_pipeline(fasta_path, metadata_path, output_csv_path):
    print("Loading metadata and sequences")
    meta_df = pd.read_csv(metadata_path)
    
    # Load all sequences into a dictionary using the Validation_ID
    records = {rec.id: str(rec.seq) for rec in SeqIO.parse(fasta_path, "fasta")}
    results = []

    print(f"Beginning pipeline execution for {len(meta_df)} proteins\n")
    
    for index, row in meta_df.iterrows():
        val_id = row['Validation_ID']
        pathogen = row['Pathogen']
        sequence = records[val_id]

        print(f"[{index+1}/{len(meta_df)}] Processing {val_id} ({pathogen})")

        # 1. Prepare temp file for Module 1 (IAPred requires a file path)
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp_fasta:
            tmp_fasta.write(f">{val_id}\n{sequence}\n")
            tmp_fasta_path = Path(tmp_fasta.name)

        try:
            # MODULE 1: Antigenicity
            raw_ag, cat = antigenicity_score(tmp_fasta_path)
            mod1 = antigenicity_normalization(raw_ag) if raw_ag is not None else None

            # MODULE 2: Topology
            top_res = analyze_topology(sequence)
            mod2 = top_res["surface_normalized"]

            # MODULE 3: Host Similarity
            mod3 = score_homology(sequence)

            # MODULE 4: Conservation
            strain_key = PATHOGEN_TO_STRAIN_KEY.get(pathogen)
            strain_list = STRAIN_MAP.get(strain_key, [])
            mod4 = conservation_score(sequence, strain_list)

            # MODULE 5: Usable Region
            mod5_res = calculate_module5_from_module2(sequence, top_res)
            mod5 = mod5_res["usable_region_score"]

            # FINAL ACS
            acs = calculate_acs(mod1, mod2, mod3, mod4, mod5)

            results.append({
                "Validation_ID": val_id,
                "Mod1_Antigenicity": mod1,
                "Mod2_Topology": mod2,
                "Mod3_HostSim": mod3,
                "Mod4_Conservation": mod4,
                "Mod5_UsableRegion": mod5,
                "ACS": acs
            })

        finally:
            # Cleanup temp files
            if tmp_fasta_path.exists():
                os.remove(tmp_fasta_path)
            iapred_output = current_dir / f"{tmp_fasta_path.stem}.csv"
            if iapred_output.exists():
                os.remove(iapred_output)

    # Merge results with metadata and save
    results_df = pd.DataFrame(results)
    final_df = pd.merge(meta_df, results_df, on="Validation_ID")
    final_df.to_csv(output_csv_path, index=False)
    
    print("\n========================================")
    print(f"VALIDATION COMPLETE! Output saved to:\n{output_csv_path}")
    print("========================================\n")

if __name__ == "__main__":
    val_fasta = current_dir / "validation" / "output" / "validation_proteins.fasta"
    val_meta = current_dir / "validation" / "output" / "validation_metadata.csv"
    output_csv = current_dir / "validation" / "output" / "validation_results_scored.csv"
    
    run_validation_pipeline(val_fasta, val_meta, output_csv)