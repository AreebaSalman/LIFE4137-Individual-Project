# Understanding Rearrangements in CNS Tumors  

**Author:** Areeba Salman  
**Degree:** MSc Bioinformatics – Individual Project (2025)  
**Institution:** University of Nottingham  

---

## Overview
This repository accompanies my MSc thesis *“Understanding Rearrangements in CNS Tumors.”*  

It provides **scripts and workflows** for reproducing the code and figures (Circos plots, CNVs, IGV screenshots) described in the thesis.  

The project compares two Nanopore sequencing workflows:  
- **ROBIN** → baseline workflow (methylation classification, SV/CNV calling)  
- **BATMAN** → extended workflow (adaptive sampling, breakpoint enrichment)  

---

## Tools & Environments

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

.
├── BATMAN/
│ ├── make_circos_batman.py # Circos plot with SV links
│ └── example_batman_sv.csv # Example input for Fig. 8a
│
├── ROBIN/
│ ├── make_circos_robin_blank.py # Chromosome circle only
│ └── example_robin.csv # Empty placeholder
│
├── data/ # Example SV/CNV input tables
├── results/ # Output directory for PNG/PDF figures
├── envs/ # Conda environment YAMLs
├── README.md # This file

## Setup

<details>
<summary><b>Click to expand installation steps</b></summary>

### Clone
```bash
git clone https://github.com/AreebaSalman/LIFE4137-Individual-Project.git
cd LIFE4137-Individual-Project
Conda environment (minimal)

conda create -n SV
conda activate SV

conda install bioconda::sniffles= 2.6.3
conda install bioconda::cutesv= 2.1.2
conda install bioconda::svim= 2.0.0
conda install bioconda::survivor= 1.0.3
conda install bioconda::samtools= 1.21
conda install bioconda::samplot= 1.3.0
conda install bioconda::circus= v1.10.0
conda install bioconda::cnvkit= 0.9.10
conda install bioconda::bedtools= v2.31.1

</details>
Workflow Steps
<details> <summary><b>1. Reference Genome Preparation</b></summary>

wget ftp://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
</details> <details> <summary><b>2. Structural Variant Calling</b></summary>

# Sniffles2
sniffles --input sample.bam --vcf sniffles.vcf --minsupport 5 --output-rnames

# cuteSV
cuteSV sample.bam ref.fa cutesv.vcf ./cutesv_output

# SVIM
svim alignment svim_out/ sample.bam ref.fa
</details> <details> <summary><b>3. Merging & Filtering with SURVIVOR</b></summary>

SURVIVOR merge input_files.txt 500 2 1 1 0 30 merged.vcf
bcftools index merged.vcf
</details> <details> <summary><b>4. Visualization & Validation</b></summary>
IGV: inspect breakpoints (coverage dips, split/supplementary reads).

Samplot:
samplot plot -n SAMPLE -b sample.bam -c chr3 -s 75376331 -e 75591822 -t DEL -o sv.png
Ribbon: load BAM at genomeribbon.com.

Circos plots:
python BATMAN/make_circos_batman.py --csv BATMAN/example_batman_sv.csv --out_prefix results/chromosomal_map_BATMAN
python ROBIN/make_circos_robin_blank.py --out_prefix results/chromosomal_map_ROBIN
</details> <details> <summary><b>5. CNV Analysis (CNVkit)</b></summary>

cnvkit.py batch sample.bam --fasta ref.fa --output-dir cnvkit_out/ --diagram --scatter
</details> <details> <summary><b>6. Annotation</b></summary>

bedtools intersect -a merged.vcf -b COSMIC_CGCv96.bed -wa -wb > SVs_with_genes.tsv
Annotation also performed in R with VariantAnnotation and biomaRt.

</details>
 Input Data
Full BAM/VCF files are not included (size/privacy).

Example CSVs with SV coordinates are provided for reproducing figures.

All data preparation steps are detailed in the thesis.

Figures Reproduced
BATMAN Circos plot (SV links)

ROBIN Circos plot (ideogram only)

CNV plots – CNVkit 

IGV screenshots – deletions, duplications, inversions, translocations

Samplot – standardized SV visualizations

Ribbon – read-level translocations

📜 Data Availability
Scripts and summary tables are hosted here.
Full BAM/VCF files are stored on Ada HPC (University of Nottingham) and available upon request.

🙏 Acknowledgements
Supervisor: Prof. Matthew Loose

Mentors: Thomas Murray, Simon Deacon

Ada HPC team (University of Nottingham)

📌 Citation
Salman, A. (2025). Understanding Rearrangements in CNS Tumors. MSc Thesis, University of Nottingham.
GitHub: https://github.com/AreebaSalman/LIFE4137-Individual-Project




