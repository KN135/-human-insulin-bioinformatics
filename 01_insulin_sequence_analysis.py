from Bio import SeqIO
from pathlib import Path

# Project folder
BASE_DIR = Path(
    r"C:\Users\USER\insulin-bioinformatics-project"
)

# FASTA file location
input_file = Path(
    r"C:\Users\USER\insulin-bioinformatics-project\data\raw\insulin_mrna.fasta"
) 

# Read sequence
record = SeqIO.read(str(input_file), "fasta")

sequence = str(record.seq).upper()

# Calculate sequence length
length = len(sequence)

# Nucleotide composition
composition = {
    base: sequence.count(base)
    for base in "ATGC"
}

# GC content
gc_content = (
    (composition["G"] + composition["C"]) / length
) * 100

# Display results
print("HUMAN INSULIN GENE ANALYSIS")
print("===========================")

print("Sequence ID:", record.id)
print("Sequence Length:", length)
print("Nucleotide Composition:", composition)
print("GC Content:", round(gc_content, 2), "%")


from pathlib import Path
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt

# ==========================================
# HUMAN INSULIN mRNA SEQUENCE ANALYSIS
# ==========================================

print("HUMAN INSULIN mRNA ANALYSIS")
print("===========================")

# Project directories
BASE_DIR = Path(__file__).resolve().parent.parent

DATA_DIR = BASE_DIR / "data" / "raw"
RESULTS_DIR = BASE_DIR / "results"
FIGURES_DIR = RESULTS_DIR / "figures"
REPORTS_DIR = RESULTS_DIR / "reports"

# Create output directories if they do not exist
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
REPORTS_DIR.mkdir(parents=True, exist_ok=True)

# Input FASTA file
input_file = DATA_DIR / "insulin_mrna.fasta"

# Check whether input file exists
if not input_file.exists():
    raise FileNotFoundError(
        "FASTA file not found: " + str(input_file)
    )

# Read FASTA sequence
record = SeqIO.read(str(input_file), "fasta")

sequence = str(record.seq).upper()

# Count nucleotides
nucleotides = ["A", "T", "G", "C"]

composition = {
    nucleotide: sequence.count(nucleotide)
    for nucleotide in nucleotides
}

# Calculate sequence length
sequence_length = len(sequence)

# Calculate GC content using unambiguous A, T, G and C
valid_bases = sum(composition.values())

if valid_bases == 0:
    raise ValueError("No valid A, T, G or C nucleotides found.")

gc_content = (
    (composition["G"] + composition["C"])
    / valid_bases
) * 100

# Display results
print("Sequence ID:", record.id)
print("Sequence Length:", sequence_length)
print("Nucleotide Composition:", composition)
print("GC Content:", round(gc_content, 2), "%")

# ==========================================
# SAVE NUCLEOTIDE COMPOSITION CSV
# ==========================================

df = pd.DataFrame(
    list(composition.items()),
    columns=["Nucleotide", "Count"]
)

csv_file = RESULTS_DIR / "nucleotide_composition.csv"
df.to_csv(str(csv_file), index=False)

print("\nCSV saved:", csv_file)

# ==========================================
# GENERATE NUCLEOTIDE COMPOSITION GRAPH
# ==========================================

plt.figure(figsize=(8, 5))

plt.bar(
    list(composition.keys()),
    list(composition.values())
)

plt.title("Human Insulin mRNA Nucleotide Composition")
plt.xlabel("Nucleotide")
plt.ylabel("Count")

plt.tight_layout()

graph_file = FIGURES_DIR / "nucleotide_composition.png"
plt.savefig(str(graph_file), dpi=300)

plt.show()
plt.close()

print("Graph saved:", graph_file)

# ==========================================
# SAVE ANALYSIS REPORT
# ==========================================

report_file = REPORTS_DIR / "sequence_analysis_report.txt"

with open(str(report_file), "w") as report:
    report.write("HUMAN INSULIN mRNA ANALYSIS REPORT\n")
    report.write("==================================\n\n")

    report.write("Sequence ID: " + record.id + "\n")
    report.write("Sequence Length: " + str(sequence_length) + "\n\n")

    report.write("NUCLEOTIDE COMPOSITION\n")
    report.write("----------------------\n")

    for nucleotide, count in composition.items():
        report.write(
            nucleotide + ": " + str(count) + "\n"
        )

    report.write(
        "\nGC Content: " + str(round(gc_content, 2)) + "%\n"
    )

    report.write(
        "\nGC content is calculated using unambiguous A, T, G and C bases.\n"
    )

print("Report saved:", report_file)

print("\n===================================")
print("ANALYSIS COMPLETED SUCCESSFULLY!")
print("===================================")



