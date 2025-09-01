Understanding Rearrangements in CNS Tumors

Author: Areeba Salman
Degree: MSc Bioinformatics – Individual Project (2025)
Institution: University of Nottingham

📌 Background

This project investigates structural variants (SVs) and copy number variations (CNVs) in CNS tumor genomes using Oxford Nanopore long-read sequencing.

Two workflows were compared:

ROBIN → baseline pipeline (methylation classification, SV & CNV calling)

BATMAN → extended pipeline (adaptive sampling with real-time breakpoint enrichment)

The aim was to identify rearrangements and validate them through multiple variant callers and visualization tools.

📂 Input Data – Obtaining BAM files

Basecalled Nanopore reads were processed through the ROBIN and BATMAN workflows.

Final BAM files for each sample were used as input for downstream SV calling.

These BAMs are not hosted here due to size restrictions, but all subsequent scripts assume you start from a sorted BAM.

🧬 Reference Genome Preparation

We used GRCh38.p14 from Ensembl as the reference genome.

# Download reference genome from Ensembl
wget ftp://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Unzip
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Index reference for samtools
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa

🔎 Structural Variant Detection

Three SV callers were used: Sniffles2, cuteSV, and SVIM.

Each was installed and run in its own conda environment.

Install environments
# Create environment
conda create -n sv-callers -y
conda activate sv-callers

# Install tools
conda install -c bioconda sniffles=2.0
conda install -c bioconda cutesv
conda install -c bioconda svim

Run commands
# Sniffles2
sniffles --input sample.bam --vcf sniffles.vcf --minsupport 5 --output-rnames

# cuteSV
cuteSV sample.bam Homo_sapiens.GRCh38.dna.primary_assembly.fa cutesv.vcf ./cutesv_output

# SVIM
svim alignment svim_out/ sample.bam Homo_sapiens.GRCh38.dna.primary_assembly.fa

🔗 Merging and Filtering of SVs

We merged VCF outputs from the three callers using SURVIVOR, requiring at least 2 callers to support an event and a maximum distance of 500 bp.

Install SURVIVOR
conda install -c bioconda survivor

Run command
# input_files.txt contains the VCFs (sniffles.vcf, cutesv.vcf, svim.vcf)
SURVIVOR merge input_files.txt 500 2 1 1 0 30 merged.vcf

🎯 Variant Filtering and Prioritization

To reduce false positives and focus on biologically meaningful events:

Variants with QUAL > 40 were kept.

For each SV type (DEL, DUP, INV, TRA), the top 10 events by size were retained.

This filtering was performed downstream in R/Python scripts (provided separately).

👁 Visualization and Validation

Variants were validated using multiple visualization tools.

IGV

Manual inspection of BAM alignments for:

Deletions → drop in coverage + split reads

Duplications → spike in coverage

Inversions → inward/outward facing reads

Translocations → split reads mapping to other chromosomes

Prepare files for IGV:

samtools index sample.bam
bcftools index merged.vcf

Samplot
conda install -c bioconda samplot
samplot plot -n SAMPLE -b sample.bam -c chr3 -s 75376331 -e 75591822 -t DEL -o sv.png

Ribbon

Upload BAM files to Ribbon
 for read-level visualization of complex events.

Circos-style Plots

Python scripts included here were used for genome-wide SV mapping:

# BATMAN Circos plot
python BATMAN/make_circos_batman.py --csv BATMAN/example_batman_sv.csv --out_prefix results/chromosomal_map_BATMAN

# ROBIN Circos plot
python ROBIN/make_circos_robin_blank.py --out_prefix results/chromosomal_map_ROBIN

📊 Comparison of Variant Callers

We compared Sniffles2, cuteSV, and SVIM:

Venn diagrams generated in R (eulerr package) to show overlaps.

Concordance evaluated at merged events (SURVIVOR).

ROBIN vs BATMAN outputs were compared to highlight the effect of adaptive sampling.

🧬 Annotation

To interpret the biological relevance, SVs were annotated against cancer gene databases.

Install bedtools
conda install -c bioconda bedtools

Run annotation
bedtools intersect -a merged.vcf -b COSMIC_CGCv96.bed -wa -wb > SVs_with_genes.tsv

R-based annotation

biomaRt for gene annotation

VariantAnnotation for VCF parsing

ggplot2/dplyr for visualization

📜 Data Availability

Scripts and summary tables are available in this repository.

BAM/VCF files are too large and not uploaded, but are stored on Ada HPC and available upon request.

🙏 Acknowledgements

Supervisor: Prof. Matthew Loose

Mentors: Thomas Murray, Simon Deacon

Ada HPC team (University of Nottingham)

📌 Citation

Salman, A. (2025). Understanding Rearrangements in CNS Tumors. MSc Thesis, University of Nottingham.
GitHub: https://github.com/AreebaSalman/LIFE4137-Individual-Project
