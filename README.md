# Understanding Rearrangements in CNS Tumors  

**Author:** Areeba Salman  
**Degree:** MSc Bioinformatics – Individual Project (2025)  
**Institution:** University of Nottingham  

---

## 📌 Overview
This repository accompanies my MSc thesis *“Understanding Rearrangements in CNS Tumors.”*  

It provides **scripts and workflows** for reproducing the code and figures (Circos plots, CNVs, IGV screenshots) described in the thesis.  

The project compares two Nanopore sequencing workflows:  
- **ROBIN** → baseline workflow (methylation classification, SV/CNV calling)  
- **BATMAN** → extended workflow (adaptive sampling, breakpoint enrichment)  

---

## 🧰 Tools & Environments

<details>
<summary><b>Click to expand</b></summary>

### SV Calling
- Sniffles2  
- cuteSV  
- SVIM  
- Epi2me `wf-human-variation`

### SV Merging & Filtering
- SURVIVOR  

### Visualization
- IGV  
- Samplot  
- Ribbon  
- Circos-style Python scripts (in this repo)  

### CNV Analysis
- CNVkit  

### Annotation
- Bedtools  
- COSMIC CGC v96  
- Ensembl biomaRt (R)  

### Supporting
- Samtools  
- Bcftools  
- R (`VariantAnnotation`, `ggplot2`, `dplyr`, `eulerr`)  

</details>

---

## 📂 Repository Structure


