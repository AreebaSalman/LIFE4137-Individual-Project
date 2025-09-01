# Understanding Rearrangements in CNS Tumors  

**Author:** Areeba Salman  
**Degree:** MSc Bioinformatics – Individual Project (2025)  
**Institution:** University of Nottingham  

---

## 📌 Background  

This project investigates **structural variants (SVs)** and **copy number variations (CNVs)** in CNS tumor genomes using Oxford Nanopore long-read sequencing.  

Two workflows were compared:  
- **ROBIN** → baseline pipeline (methylation classification, SV & CNV calling)  
- **BATMAN** → extended pipeline (adaptive sampling with real-time breakpoint enrichment)  

The aim was to identify genomic rearrangements and validate them through multiple variant callers and visualization tools.  

---

## 📂 Input Data – BAM Files  

- Basecalled Nanopore reads were processed using **ROBIN** and **BATMAN** workflows.  
- This produced sorted **BAM files**, which served as input for SV detection.  
- BAM files are **not uploaded** here due to size/privacy.  

---

## 🧬 Reference Genome Preparation  

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

conda install -c bioconda sniffles=2.0
conda install -c bioconda cutesv
conda install -c bioconda svim
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
