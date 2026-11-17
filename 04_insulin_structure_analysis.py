print("INSULIN PROJECT TEST")
print("Python is working!")
input("Press Enter to finish...")

from pathlib import Path
from Bio.PDB import PDBParser
import pandas as pd

print("HUMAN INSULIN 3D STRUCTURE ANALYSIS")
print("===================================")

# Project directories
BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data" / "raw"
RESULTS_DIR = BASE_DIR / "results"
REPORTS_DIR = RESULTS_DIR / "reports"

REPORTS_DIR.mkdir(parents=True, exist_ok=True)

# PDB input file
pdb_file = DATA_DIR / "insulin_4EXX.pdb"

# Read structure
parser = PDBParser(QUIET=True)
structure = parser.get_structure("insulin", str(pdb_file))

# Analyze first model
model = structure[0]

chain_data = []

for chain in model:
    residues = [residue for residue in chain if residue.id[0] == " "]

    atom_count = sum(
        len(list(residue.get_atoms()))
        for residue in residues
    )

    chain_data.append({
        "Chain": chain.id,
        "Residue_Count": len(residues),
        "Atom_Count": atom_count
    })
    
df = pd.DataFrame(chain_data)

print("\nCHAIN INFORMATION")
print(df.to_string(index=False))

print("\nNumber of chains:", len(chain_data))

total_residues = sum(
    item["Residue_Count"] for item in chain_data
)

print("Total protein residues:", total_residues)

# Save CSV
csv_file = RESULTS_DIR / "insulin_structure_chains.csv"
df.to_csv(str(csv_file), index=False)

# Save report
report_file = REPORTS_DIR / "insulin_structure_report.txt"

with open(str(report_file), "w") as report:
    report.write("HUMAN INSULIN STRUCTURE ANALYSIS\n")
    report.write("PDB ID: 4EXX\n\n")
    report.write("Number of chains: " + str(len(chain_data)) + "\n")
    report.write("Total protein residues: " + str(total_residues) + "\n\n")
    report.write(df.to_string(index=False))

print("\nCSV saved:", csv_file)
print("Report saved:", report_file)

print("\nSTRUCTURE ANALYSIS COMPLETED!")


