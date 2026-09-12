from pathlib import Path
from Bio import SeqIO
import csv
import random
import re

BASE_DIR = Path(__file__).resolve().parent
ANTIGEN_DIR = BASE_DIR / "known_antigens"
PROTEOME_DIR = BASE_DIR / "proteomes"
OUTPUT_DIR = BASE_DIR / "output"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

CONTROL_RATIO = 5
RANDOM_SEED = 42


def extract_accession(record):
    parts = record.id.split("|")

    if len(parts) >= 2 and parts[0] == "sp":
        return parts[1]

    match = re.search(r"\b([A-Z][A-Z0-9]{5,9})\b", record.description)

    return match.group(1) if match else record.id


def identify_pathogen(filename):
    name = filename.lower()

    if "bordetella" in name:
        return "Bordetella pertussis"

    if "hepatitis_b" in name:
        return "Hepatitis B virus"

    if "hpv" in name:
        return "Human papillomavirus"

    if "mycobacterium_tuberculosis" in name:
        return "Mycobacterium tuberculosis"

    if "neisseria" in name:
        return "Neisseria meningitidis"

    if "rsv" in name:
        return "Respiratory syncytial virus"

    if "sars_cov_2" in name:
        return "SARS-CoV-2"

    if "streptococcus_pneumoniae" in name:
        return "Streptococcus pneumoniae"

    return "Unknown"


def read_single_fasta(path):
    records = list(SeqIO.parse(path, "fasta"))

    if not records:
        raise ValueError(f"No FASTA sequence found in: {path.name}")

    if len(records) > 1:
        raise ValueError(f"Expected one protein in antigen file, but found multiple: {path.name}")

    return records[0]


def read_proteome(path):
    records = list(SeqIO.parse(path, "fasta"))

    if not records:
        raise ValueError(f"No FASTA sequences found in proteome: {path.name}")

    return records


def main():
    antigen_files = sorted(ANTIGEN_DIR.glob("*.fasta"))
    proteome_files = sorted(PROTEOME_DIR.glob("*.fasta"))

    if not antigen_files:
        raise FileNotFoundError(f"No antigen FASTA files found in {ANTIGEN_DIR}")

    if not proteome_files:
        raise FileNotFoundError(f"No proteome FASTA files found in {PROTEOME_DIR}")

    positives = []
    positive_ids = set()

    print("\nReading known antigen FASTA files...")

    for path in antigen_files:
        record = read_single_fasta(path)
        accession = extract_accession(record)
        pathogen = identify_pathogen(path.name)

        if pathogen == "Unknown":
            print(f"Warning: pathogen could not be identified from {path.name}")
            continue

        if accession in positive_ids:
            print(f"Skipping duplicate positive: {accession}")
            continue

        positive_ids.add(accession)

        positives.append({
            "record": record,
            "protein_id": accession,
            "protein_name": record.description,
            "pathogen": pathogen,
            "source": path.name,
            "label": 1,
            "evidence_status": "Experimentally characterized antigen set"
        })

    print(f"Known antigen proteins: {len(positives)}")

    proteome_by_pathogen = {}

    print("\nReading pathogen proteomes...")

    for path in proteome_files:
        pathogen = identify_pathogen(path.name)

        if pathogen == "Unknown":
            print(f"Warning: pathogen could not be identified from {path.name}")
            continue

        records = read_proteome(path)

        if pathogen not in proteome_by_pathogen:
            proteome_by_pathogen[pathogen] = {}

        for record in records:
            accession = extract_accession(record)

            if accession in positive_ids:
                continue

            if accession not in proteome_by_pathogen[pathogen]:
                proteome_by_pathogen[pathogen][accession] = {
                    "record": record,
                    "protein_id": accession,
                    "protein_name": record.description,
                    "pathogen": pathogen,
                    "source": path.name,
                    "label": 0,
                    "evidence_status": "Background protein without qualifying positive-set evidence"
                }

    positive_counts = {}

    for item in positives:
        positive_counts[item["pathogen"]] = positive_counts.get(item["pathogen"], 0) + 1

    controls = []
    random.seed(RANDOM_SEED)

    print("\nSelecting background proteins...")

    for pathogen, count in positive_counts.items():
        candidates = list(proteome_by_pathogen.get(pathogen, {}).values())

        if not candidates:
            print(f"Warning: no background proteins found for {pathogen}")
            continue

        target_controls = min(count * CONTROL_RATIO, len(candidates))
        selected = random.sample(candidates, target_controls)

        controls.extend(selected)

        print(f"{pathogen}: {count} positives, {target_controls} controls")

    dataset = positives + controls

    if not dataset:
        raise RuntimeError("Validation dataset is empty.")

    random.shuffle(dataset)

    fasta_path = OUTPUT_DIR / "validation_proteins.fasta"
    metadata_path = OUTPUT_DIR / "validation_metadata.csv"

    with open(fasta_path, "w", encoding="utf-8") as fasta_file:
        for index, item in enumerate(dataset, start=1):
            record = item["record"]
            new_id = f"VAL_{index:04d}_{item['protein_id']}"

            sequence = str(record.seq).replace(" ", "").replace("\n", "")

            fasta_file.write(f">{new_id} {item['protein_name']}\n")
            fasta_file.write(f"{sequence}\n")

    with open(metadata_path, "w", newline="", encoding="utf-8") as csv_file:
        writer = csv.writer(csv_file)

        writer.writerow([
            "Validation_ID",
            "Protein_ID",
            "Protein_Name",
            "Pathogen",
            "Source",
            "Label",
            "Evidence_Status"
        ])

        for index, item in enumerate(dataset, start=1):
            writer.writerow([
                f"VAL_{index:04d}_{item['protein_id']}",
                item["protein_id"],
                item["protein_name"],
                item["pathogen"],
                item["source"],
                item["label"],
                item["evidence_status"]
            ])

    print("\n========================================")
    print("VALIDATION DATASET CREATED")
    print("========================================")
    print(f"Total proteins: {len(dataset)}")
    print(f"Positive proteins: {len(positives)}")
    print(f"Control/background proteins: {len(controls)}")

    print("\nPositive proteins by pathogen:")

    for pathogen, count in sorted(positive_counts.items()):
        print(f"  {pathogen}: {count}")

    print(f"\nFASTA: {fasta_path}")
    print(f"Metadata: {metadata_path}")
    print("========================================\n")


if __name__ == "__main__":
    main()