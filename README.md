# Understanding Rearrangements in CNS Tumors  

**Author:** Areeba Salman  
**Degree:** MSc Bioinformatics – Individual Project (2025)  
**Institution:** University of Nottingham  

---

## Background  

This project investigates **structural variants (SVs)** and **copy number variations (CNVs)** in CNS tumor genomes using Oxford Nanopore long-read sequencing.  

Two workflows were compared:  
- **ROBIN** → baseline pipeline (methylation classification, SV & CNV calling)  
- **BATMAN** → extended pipeline (adaptive sampling with real-time breakpoint enrichment)  

This study aims to investigate structural variants (SVs) (Hurles et al. 2008) in central nervous system (CNS) tumors using long read nanopore sequencing, and to evaluate the performance and effectiveness of two analytical pipelines, ROBIN(Deacon et al. 2025) and BATMAN, in detecting these variants.   

---

## Input Data – BAM Files  

- Basecalled Nanopore reads were processed using **ROBIN** and **BATMAN** workflows.  
- This produced sorted **BAM files**, which served as input for SV detection.  
- BAM files are **not uploaded** here due to size/privacy.  

---

## Reference Genome Preparation  

The **GRCh38.p14** reference genome was downloaded from Ensembl and indexed for downstream analysis.  

```bash
# Download reference genome
wget ftp://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Unzip
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Index with samtools
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
```

## Structural Variant Detection

SVs were identified using Sniffles2, cuteSV, and SVIM.

### Installation
```bash
conda create -n sv-callers -y
conda activate sv-callers

conda install -c bioconda sniffles= 2.6.3
conda install -c bioconda cutesv= 2.1.2
conda install -c bioconda svim= 2.0.0
```

### Running
```bash
# Sniffles2
sniffles --input sample.bam --vcf sniffles.vcf --minsupport 5 --output-rnames

# cuteSV
cuteSV sample.bam Homo_sapiens.GRCh38.dna.primary_assembly.fa cutesv.vcf ./cutesv_output

# SVIM
svim alignment svim_out/ sample.bam Homo_sapiens.GRCh38.dna.primary_assembly.fa
```

## Merging and Filtering of SVs

SV calls from all three tools were merged using SURVIVOR.
We required ≥2 callers to support an event within a 500 bp window.

### Installation
```bash
conda create -n survivor -y
conda activate survivor
conda install -c bioconda survivor= 1.0.3
```

Running
```bash
# input_files.txt contains sniffles.vcf, cutesv.vcf, svim.vcf
SURVIVOR merge input_files.txt 500 2 1 1 0 30 merged.vcf
```

## Variant Filtering and Prioritization

To reduce false positives and highlight biologically relevant events:

Variants with QUAL > 40 were retained.

For each SV type (DEL, DUP, INV, TRA), the top 10 events by size were selected.

This filtering was performed using downstream R/Python scripts

## Visualization and Validation

SVs were validated using multiple visualization tools.

### IGV

Used for manual inspection of breakpoints.

Deletions → coverage drop + split reads

Duplications → coverage spike

Inversions → orientation changes

Translocations → reads mapping to multiple chromosomes

### Prepare files for IGV:
```bash
conda create -n samtools -y
conda activate samtools
conda install bioconda::samtools= 1.21

samtools index sample.bam
bcftools index merged.vcf
```

### Samplot
```bash
conda create -n samplot -y
conda activate samplot
conda install -c bioconda samplot
samplot plot -n SAMPLE -b sample.bam -c chr3 -s 75376331 -e 75591822 -t DEL -o sv.png
```

### Ribbon
Upload BAM files to Ribbon to visualize long-read evidence of translocations.

### Circos-style Plots
Python scripts included in this repo were used for genome-wide SV mapping:
```bash
# BATMAN Circos plot
python BATMAN/make_circos_batman.py --csv BATMAN/example_batman_sv.csv --out_prefix results/chromosomal_map_BATMAN

# ROBIN Circos plot
python ROBIN/make_circos_robin_blank.py --out_prefix results/chromosomal_map_ROBIN.
```

## Comparison of Variant Callers

**Sniffles2**, **cuteSV**, and **SVIM** results were compared.

**SURVIVOR** was used for merged calls.

**Venn diagrams** were generated in R using the eulerr package.

ROBIN vs BATMAN outputs were compared to evaluate the effect of adaptive sampling.

## Annotation

SVs were annotated to identify overlaps with cancer-related genes.

### Installation
```bash
conda create -n bedtools -y
conda activate bedtools
conda install -c bioconda bedtools
```
### Running
```bash
bedtools intersect -a merged.vcf -b COSMIC_CGCv96.bed -wa -wb > SVs_with_genes.tsv
```

Annotation was further refined in R using:

biomaRt for gene annotation

VariantAnnotation for VCF parsing

ggplot2 and dplyr for visualization

## Data Availability

Scripts and codes to reproduce are provided here.

Full BAM/VCF files are large and stored on Ada HPC (University of Nottingham). Available upon request.

## Acknowledgements

Supervisor: Prof. Matthew Loose

Mentors: Thomas Murray, Simon Deacon

Ada HPC team (University of Nottingham)

## Citation
Salman, A. (2025). Understanding Rearrangements in CNS Tumors. MSc Thesis, University of Nottingham.
GitHub: https://github.com/AreebaSalman/LIFE4137-Individual-Project
