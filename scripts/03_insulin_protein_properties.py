
from pathlib import Path
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd

print("HUMAN INSULIN PROTEIN PROPERTIES")
print("================================")

# Project folders
BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data" / "processed"
RESULTS_DIR = BASE_DIR / "results"
REPORTS_DIR = RESULTS_DIR / "reports"

REPORTS_DIR.mkdir(parents=True, exist_ok=True)

# Read translated protein FASTA
protein_file = DATA_DIR / "insulin_protein.fasta"

if not protein_file.exists():
    raise FileNotFoundError("Protein FASTA file not found.")

record = SeqIO.read(str(protein_file), "fasta")
protein_sequence = str(record.seq)

# Analyze protein
analysis = ProteinAnalysis(protein_sequence)

molecular_weight = analysis.molecular_weight()
isoelectric_point = analysis.isoelectric_point()
instability_index = analysis.instability_index()
aromaticity = analysis.aromaticity()
gravy = analysis.gravy()

# Display results
print("Protein ID:", record.id)
print("Protein Length:", len(protein_sequence))
print("Molecular Weight:", round(molecular_weight, 2), "Da")
print("Theoretical pI:", round(isoelectric_point, 2))
print("Instability Index:", round(instability_index, 2))
print("Aromaticity:", round(aromaticity, 4))
print("GRAVY:", round(gravy, 4))

# Save results
results = {
    "Property": [
        "Protein Length",
        "Molecular Weight (Da)",
        "Theoretical pI",
        "Instability Index",
        "Aromaticity",
        "GRAVY"
    ],
    "Value": [
        len(protein_sequence),
        round(molecular_weight, 2),
        round(isoelectric_point, 2),
        round(instability_index, 2),
        round(aromaticity, 4),
        round(gravy, 4)
    ]
}

df = pd.DataFrame(results)

csv_file = RESULTS_DIR / "insulin_protein_properties.csv"
df.to_csv(str(csv_file), index=False)

# Save text report
report_file = REPORTS_DIR / "insulin_protein_properties_report.txt"

with open(str(report_file), "w") as report:
    report.write("HUMAN INSULIN PROTEIN PROPERTIES\n")
    report.write("================================\n\n")

    for index, row in df.iterrows():
        report.write(
            str(row["Property"]) + ": "
            + str(row["Value"]) + "\n"
        )

print("\nCSV saved:", csv_file)
print("Report saved:", report_file)
print("\nPROTEIN PROPERTY ANALYSIS COMPLETED!")
