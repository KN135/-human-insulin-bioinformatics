
from pathlib import Path
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt

print("HUMAN INSULIN PROTEIN ANALYSIS")
print("==============================")

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data" / "raw"
RESULTS_DIR = BASE_DIR / "results"
FIGURES_DIR = RESULTS_DIR / "figures"
REPORTS_DIR = RESULTS_DIR / "reports"
PROCESSED_DIR = BASE_DIR / "data" / "processed"

for folder in [FIGURES_DIR, REPORTS_DIR, PROCESSED_DIR]:
    folder.mkdir(parents=True, exist_ok=True)

# Read annotated GenBank record
input_file = DATA_DIR / "insulin_mrna.gb"

if not input_file.exists():
    raise FileNotFoundError("Missing file: " + str(input_file))

record = SeqIO.read(str(input_file), "genbank")

# Find the annotated coding sequence
cds_features = [
    feature for feature in record.features
    if feature.type == "CDS"
]

if len(cds_features) != 1:
    raise ValueError(
        "Expected one CDS feature, found: " + str(len(cds_features))
    )

cds_feature = cds_features[0]

# Extract and translate annotated CDS
coding_sequence = cds_feature.extract(record.seq)
protein = coding_sequence.translate()

# Remove terminal stop symbol, if present
if str(protein).endswith("*"):
    protein = protein[:-1]

protein_string = str(protein)

print("Sequence ID:", record.id)
print("CDS length:", len(coding_sequence))
print("Protein length:", len(protein_string))
print("Protein sequence:", protein_string)

# Check expected reference protein length
if len(protein_string) != 110:
    print("WARNING: Protein length differs from expected 110 amino acids.")
else:
    print("Protein length check: PASSED")

# Save protein FASTA
from Bio.SeqRecord import SeqRecord

protein_record = SeqRecord(
    protein,
    id="human_insulin_precursor",
    description="Human insulin precursor translated from annotated CDS"
)

protein_file = PROCESSED_DIR / "insulin_protein.fasta"
SeqIO.write(protein_record, str(protein_file), "fasta")

# Amino acid composition
amino_acids = sorted(set(protein_string))

composition = {
    aa: protein_string.count(aa)
    for aa in amino_acids
}

df = pd.DataFrame(
    list(composition.items()),
    columns=["Amino_Acid", "Count"]
)

csv_file = RESULTS_DIR / "amino_acid_composition.csv"
df.to_csv(str(csv_file), index=False)

# Generate graph
plt.figure(figsize=(12, 6))
plt.bar(composition.keys(), composition.values())
plt.title("Human Insulin Precursor Amino Acid Composition")
plt.xlabel("Amino Acid")
plt.ylabel("Count")
plt.tight_layout()

graph_file = FIGURES_DIR / "amino_acid_composition.png"
plt.savefig(str(graph_file), dpi=300)
plt.show()
plt.close()

# Save report
report_file = REPORTS_DIR / "protein_analysis_report.txt"

with open(str(report_file), "w") as report:
    report.write("HUMAN INSULIN PROTEIN ANALYSIS\n")
    report.write("Source: " + record.id + "\n")
    report.write("CDS length: " + str(len(coding_sequence)) + "\n")
    report.write("Protein length: " + str(len(protein_string)) + "\n")
    report.write("Protein sequence:\n" + protein_string + "\n\n")
    report.write("Amino acid composition:\n")

    for aa, count in composition.items():
        report.write(aa + ": " + str(count) + "\n")

print("Protein FASTA saved:", protein_file)
print("CSV saved:", csv_file)
print("Graph saved:", graph_file)
print("Report saved:", report_file)
print("ANALYSIS COMPLETED")


