# *ExoproteoML* Ensemble learning framework for eukaryotic exoproteome prediction
 
ExoproteoML is an interactive R/Shiny platform designed for the high-throughput analysis of the eukaryotic exoproteome. 
It leverages a robust Ensemble architecture (XGBoost, Random Forest, and SVM) to classify proteins into classical and non-conventional secretion pathways with state-of-the-art accuracy.

![Graphical Abstract](img/ExoproteoML_GraphicalAbstract.jpg)

## **Installation & Setup**

# *All the files and manuals are in the zip file that cam be download from the main branch of the repository.*

Clone the Repository: 

```bash
git clone https://github.com/usuario/ExoproteoML.git
```

**Download the Model:** Due to its size (700MB), the trained ensemble bundle is hosted on Zenodo.
	
Download: https://zenodo.org/uploads/record/18754701 
	
**Placement:** Place the Ensemble_Bundle.rds file inside the models/ folder.

**Install dependencies**

```R
# 1. CRAN 
dependencies_cran <- c("shiny", "shinythemes", "tidymodels", "DT", "tidyverse", 
                   "Peptides", "plotly", "htmlwidgets", "cluster", "factoextra")

install.packages(dependencies_cran, dependencies = TRUE)

# 2. Bioconductor 
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("protr", "Biostrings"))
```

**Run the App:** Open app.R in RStudio and click "Run App", 
or run:

```R
shiny::runApp()
```

Ensure all dependencies (shiny, tidymodels, protr, Biostrings, plotly) are installed.

## 📖 **User Manual & Modules**

### **Data Input**
• FASTA Upload: Supports large proteome files (up to 50MB).

• Manual Paste: Quickly test single or multiple sequences in FASTA format.

• Validation: The app automatically cleans sequences, removing non-standard amino acids and gaps to ensure descriptor integrity.
    
### **Full Results Table**

This module provides a detailed breakdown of every protein:

• **Veredict:** Final classification—Conventional (SP+), Non-Conventional (UPS), Non-Secreted, or SP+ (ER-Retained) (detected via C-terminal KDEL/HDEL motifs).

• **Confidence:** The consensus probability of the Ensemble (0-100%).

• **Is_SSP:** Flags candidates for Small Secreted Proteins/Effectors.

• **Physicochemical Data:** Automated calculation of Molecular Weight (kDa), Charge (pH 7), and Hydrophobicity.

    
### **Global Analytics**

• **CAZy Potential:** Visualizes the proportion of the secretome with high enzymatic affinity.

• **SSP/Effector Distribution:** Highlighting the "Small Secreted Protein" space, crucial for plant-pathogen interaction studies.
    
### **3D Secretome Space (Bio-Discovery)**
	
Using Principal Component Analysis (PCA) and K-means Clustering, the app maps your proteins into a 3D interactive universe based on their physicochemical properties.
This allows for the identification of functional clusters and "dark proteome" outliers.

## **Interpretation of Results**

To help you navigate the biological significance of the output, use the following thresholds as a guide:
 
### **Secretion Pathways**

• **Conventional (SP+):** Proteins containing a classical N-terminal Signal Peptide recognized by the Sec translocon.

• **Non-Conventional (UPS):** Proteins secreted via Unconventional Protein Secretion. These lack a signal peptide but are identified by the Ensemble's global feature analysis.

• **ER-Retained:** Proteins with a Signal Peptide that also carry a retention motif, indicating they likely remain in the Endoplasmic Reticulum. Functional Indicators • CAZy Index (> 25): A high aromatic residue density (F, W, Y) often correlates with carbohydrate-binding modules. A score above 25 suggests a "High Potential" for being a Carbohydrate-Active Enzyme.

• **Candidate SSP/Effector:** Proteins flagged with a "★" meet the "Small Secreted Protein" criteria: Length < 250 aa and Cysteine content > 3%. These are prime candidates for virulence factors or signaling molecules.

• **Confidence Scores:**
o > 85%: High reliability; likely candidates for experimental validation.
o 50% - 70%: Putative candidates; may represent proteins with "moonlighting" functions or complex membrane topologies.

## **Exporting Data**

• **TSV Report:** A "Safe-Format" tab-separated file that prevents protein IDs with commas from breaking your columns in Excel/LibreOffice.

• **Plots:** All charts can be exported as Interative HTML (to keep the 3D rotation and tooltips) or SVG (for publication-ready figures). 

## **Citation** If you use **ExoproteoML** in your research, please cite:

**Lavin, JL and Oguiza, JA (2026).** ExoproteoML: An Integrated Ensemble Learning Framework for Eukaryotic Secretome and Effector Discovery. [Journal Name / DOI]
