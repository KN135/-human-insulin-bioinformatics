# 🧬 Human Insulin Bioinformatics Analysis

### Computational Analysis of Human Insulin at the Nucleotide, Protein and Structural Levels

![Python](https://img.shields.io/badge/Python-3.7%2B-blue)
![Biopython](https://img.shields.io/badge/Biopython-1.81-green)
![Pandas](https://img.shields.io/badge/Pandas-1.3.5-blue)
![Matplotlib](https://img.shields.io/badge/Matplotlib-3.2.2-orange)
![Domain](https://img.shields.io/badge/Domain-Bioinformatics-purple)

---

## 🔬 Project Overview

This project presents a computational bioinformatics analysis of human insulin using Python-based biological data analysis tools.

The study explores human insulin at multiple biological levels, beginning with nucleotide sequence analysis, followed by protein sequence translation, physicochemical property prediction and three-dimensional protein structure analysis.

The project demonstrates the practical application of bioinformatics programming, biological sequence analysis, data processing and scientific visualization.

## 🎯 Project Objectives

- Analyze human insulin mRNA nucleotide composition.
- Calculate nucleotide frequency and GC content.
- Extract and translate the annotated coding sequence.
- Determine amino acid composition.
- Calculate important physicochemical properties of the insulin precursor.
- Examine the chain composition of a human insulin 3D structural model.
- Generate scientific graphs, CSV datasets and analysis reports.

## 🛠️ Technologies and Tools

| Technology | Application |
|---|---|
| Python | Bioinformatics programming |
| Biopython | Sequence processing and PDB analysis |
| Pandas | Data processing and CSV generation |
| Matplotlib | Scientific visualization |
| NCBI RefSeq | Reference nucleotide sequence |
| RCSB PDB | Protein structural data |

## 🧪 Analysis Workflow

### 01. Nucleotide Sequence Analysis

**Reference:** NCBI RefSeq — `NM_000207.3`

The human insulin mRNA sequence was analyzed to determine nucleotide composition and GC content.

| Parameter | Result |
|---|---:|
| Sequence length | 465 nt |
| Adenine (A) | 91 |
| Thymine (T) | 77 |
| Guanine (G) | 141 |
| Cytosine (C) | 156 |
| GC content | 63.87% |

**Output files:**
- `nucleotide_composition.csv`
- `nucleotide_composition.png`
- `sequence_analysis_report.txt`

### 02. Protein Sequence Analysis

The annotated coding sequence was extracted and translated into its corresponding amino acid sequence.

| Parameter | Result |
|---|---:|
| CDS length | 333 nt |
| Protein precursor length | 110 amino acids |

The resulting protein sequence represents the human insulin precursor.

**Output files:**
- `insulin_protein.fasta`
- `amino_acid_composition.csv`
- `amino_acid_composition.png`
- `protein_analysis_report.txt`

### 03. Protein Physicochemical Properties

Biopython's ProteinAnalysis module was used to calculate important theoretical physicochemical properties.

| Property | Result |
|---|---:|
| Molecular weight | 11,980.79 Da |
| Theoretical isoelectric point | 5.22 |
| Instability index | 40.33 |
| Aromaticity | 0.0818 |
| GRAVY | 0.1927 |

These values describe the sequence-based properties of the 110-amino-acid insulin precursor.

**Output files:**
- `insulin_protein_properties.csv`
- `insulin_protein_properties_report.txt`

### 04. Three-Dimensional Structure Analysis

**PDB ID:** `4EXX`

The downloaded structural model was parsed using Biopython's PDBParser.

The analysis examined protein chain distribution, residue counts and atom counts.

| Chain | Residues | Atoms |
|---|---:|---:|
| A | 21 | 163 |
| B | 30 | 241 |
| C | 21 | 163 |
| D | 30 | 241 |
| **Total** | **102** | **808** |

**Output files:**
- `insulin_structure_chains.csv`
- `insulin_structure_report.txt`
- `insulin_3d_structure.png` (structure screenshot, if included)

## 📂 Project Directory Structure

```text
human-insulin-bioinformatics/
│
├── data/
│   ├── raw/
│   │   ├── insulin_mrna.fasta
│   │   └── insulin_4EXX.pdb
│   │
│   └── processed/
│       └── insulin_protein.fasta
│
├── scripts/
│   ├── 01_insulin_sequence_analysis.py
│   ├── 02_insulin_protein_analysis.py
│   ├── 03_insulin_protein_properties.py
│   └── 04_insulin_structure_analysis.py
│
├── results/
│   ├── nucleotide_composition.csv
│   ├── amino_acid_composition.csv
│   ├── insulin_protein_properties.csv
│   ├── insulin_structure_chains.csv
│   │
│   ├── figures/
│   │   ├── nucleotide_composition.png
│   │   ├── amino_acid_composition.png
│   │   └── insulin_3d_structure.png
│   │
│   └── reports/
│       ├── sequence_analysis_report.txt
│       ├── protein_analysis_report.txt
│       ├── insulin_protein_properties_report.txt
│       └── insulin_structure_report.txt
│
├── README.md
└── .gitignore
```

*The directory tree describes the intended organization. Include only files that are actually present in your repository.*

## ⚙️ Installation and Requirements

The project uses the following Python libraries:

biopython==1.81
pandas==1.3.5
matplotlib==3.2.2

Install the dependencies using:

pip install biopython==1.81 pandas==1.3.5 matplotlib==3.2.2

▶️ Running the Project

Run the scripts in the following sequence from the project directory:

python scripts/01_insulin_sequence_analysis.py

python scripts/02_insulin_protein_analysis.py

python scripts/03_insulin_protein_properties.py

python scripts/04_insulin_structure_analysis.py

Ensure the required input files are placed in data/raw/ before running the relevant scripts.

📊 Key Results and Deliverables

This project produces:

Nucleotide composition analysis.

GC content calculation.

Translated insulin precursor FASTA sequence.

Amino acid composition dataset.

Protein physicochemical property calculations.

Protein chain and atom-count analysis.

Scientific visualization images.

Automated text reports and CSV result files.

🧠 Scientific Interpretation and Limitations

The 110-amino-acid sequence represents the insulin precursor rather than mature insulin alone.

The structural analysis describes the chains and atoms present in the downloaded PDB model.

Sequence and structural residue counts are not directly equivalent because the biological forms and processing stages differ.

Physicochemical properties are computational predictions and should not be interpreted as experimental measurements.

Structural analysis using Biopython parses atomic coordinates but does not itself perform experimental structure determination.

📚 Data Sources

NCBI RefSeq: Human insulin reference transcript, NM_000207.3

RCSB Protein Data Bank: Human insulin structural entry, 4EXX

🚀 Skills Demonstrated

Biological sequence analysis

DNA/RNA nucleotide composition

Coding sequence extraction and translation

Protein sequence analysis

Physicochemical property prediction

Protein structure data parsing

Python programming

Data visualization

Scientific reporting

Bioinformatics workflow development

🔭 Future Improvements

Protein sequence alignment and comparison.

Identification and analysis of conserved regions.

Disulfide bond analysis.

Additional protein structural measurements.

Interactive molecular visualization.

Automated integration of sequence and structural analysis.

👨‍🔬 Author
K.Nissi
Bioinformatics & Biotechnology Enthusiast

This project was developed as a practical learning and portfolio project to demonstrate computational approaches to biological sequence and protein structure analysis.
