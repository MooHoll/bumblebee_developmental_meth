#!/bin/bash

#SBATCH --job-name=TE_annot
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=110:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=hollie_marshall@hotmail.co.uk
#SBATCH --account=evo-epi

# Run script in the working directory it was submitted in
cd $SLURM_SUBMIT_DIR

### Before running this script, need to install EDTA
#git clone https://github.com/oushujun/EDTA.git ./edta/
#cd edta
#conda env create -f EDTA.yml 
#conda update --all

### Also need to make the correct files for EDTA, cds fasta: (NOT WORKING)
# Install from: https://github.com/gpertea/gffread
#gffread -g GCF_000214255.1_Bter_1.0_genomic.fa \
#-x bumblebee_cds.fa ref_Bter_1.0_top_level.gff3

### Alternative to get cds fasta
#grep "CDS" ref_Bter_1.0_top_level.gff3 > cds.txt
#cut -f 1,4,5 cds.txt > cds.bed
#wget https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz
#tar -zxvf bedtools-2.29.1.tar.gz
#cd bedtools2
#make
#/scratch/evo-epi/hjm32/bin/bedtools2/bin/fastaFromBed -fi GCF_000214255.1_Bter_1.0_genomic.fa -bed cds.bed -fo bumblebee_cds.fa


### Make the genes.bed
#awk '{if ($3 == "gene") print $0;}' ref_Bter_1.0_top_level.gff3 > genes
#cut -f1,4,5,9 genes > gene_cut
#sed -i 's/ID=.*gene=//g' gene_cut
#sed 's/;.*$//g' gene_cut > bumblebee_genes.bed

# Activate environment for EDTA
source ~/miniconda3/bin/activate EDTA

export BLASTDB_LMDB_MAP_SIZE=100000000

# Run EDTA
EDTA.pl \
--genome GCF_000214255.1_Bter_1.0_genomic.fa \
--cds bumblebee_cds.fa \
--anno 1 \
--sensitive 1 \
--exclude bumblebee_genes.bed \
--threads 20