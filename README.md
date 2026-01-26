# Introduction

TORBiT (TCR Reconstruction-Based Immunophenotyping Tool) is a computational workflow specifically designed for reconstructing T cell receptor sequences from single-cell transcriptome sequencing data and analyzing TCR repertoire dynamics. This workflow integrates multiple bioinformatics tools and methods, enabling the complete reconstruction of TCR α/β/γ/δ chain sequences at the single-cell level based on scRNA-seq raw sequencing data.

# Clone the repository

git clone https://github.com/yourusername/fastq-processor.git
cd TORBiT

# Install dependencies

pip install -r requirements.txt(未编写)

# Installation of bioinformatics tool dependencies

mamba install -c bioconda bwa samtools Trinity TRUST4 # 如果使用conda

#### **Running example**

# Single-End Sequencing

python main.py -rf ./ref/hg_tcr.fa -r1 ./test_data/soloseq_f1.fastq -o ./results -i ./ref/ -t 8

# Paired-End Sequencing

python main.py -rf ./ref/hg_tcr.fa -r1 ./test_data/soloseq_f1.fastq -r2 ./test_data/soloseq_r2.fastq -pe -o ./results -i ./ref/ -t 8

# Single-Cell Mode (Single-End)

python main.py -rf ./ref/hg_tcr.fa -r1 ./test_data/soloseq_f1.fastq -sc -bc 16 -o ./results -i ./ref/ -t 8

