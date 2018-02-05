#!/bin/bash

# Bash script to download
#  (1) GRCm38 (mm10) reference genome and annotation from Ensembl
#  (2) BK000964.3 rDNA repeat sequence and annotation from Dropbox
# and construct combined reference sequence and annotation.
#
# Author: Maurits Evers (maurits.evers@anu.edu.au)
# Last changed: 27-07-2017


# Download reference genome
if [ ! -f GRCm38+rDNA_repeat.fa ]; then
    # Download individual chromosome files
    wget ftp://ftp.ensembl.org/pub/release-89/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_rm.chromosome.1.fa.gz
    wget ftp://ftp.ensembl.org/pub/release-89/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_rm.chromosome.2.fa.gz
    wget ftp://ftp.ensembl.org/pub/release-89/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_rm.chromosome.3.fa.gz
    wget ftp://ftp.ensembl.org/pub/release-89/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_rm.chromosome.4.fa.gz
    wget ftp://ftp.ensembl.org/pub/release-89/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_rm.chromosome.5.fa.gz
    wget ftp://ftp.ensembl.org/pub/release-89/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_rm.chromosome.6.fa.gz
    wget ftp://ftp.ensembl.org/pub/release-89/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_rm.chromosome.7.fa.gz
    wget ftp://ftp.ensembl.org/pub/release-89/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_rm.chromosome.8.fa.gz
    wget ftp://ftp.ensembl.org/pub/release-89/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_rm.chromosome.9.fa.gz
    wget ftp://ftp.ensembl.org/pub/release-89/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_rm.chromosome.10.fa.gz
    wget ftp://ftp.ensembl.org/pub/release-89/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_rm.chromosome.11.fa.gz
    wget ftp://ftp.ensembl.org/pub/release-89/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_rm.chromosome.12.fa.gz
    wget ftp://ftp.ensembl.org/pub/release-89/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_rm.chromosome.13.fa.gz
    wget ftp://ftp.ensembl.org/pub/release-89/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_rm.chromosome.14.fa.gz
    wget ftp://ftp.ensembl.org/pub/release-89/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_rm.chromosome.15.fa.gz
    wget ftp://ftp.ensembl.org/pub/release-89/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_rm.chromosome.16.fa.gz
    wget ftp://ftp.ensembl.org/pub/release-89/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_rm.chromosome.17.fa.gz
    wget ftp://ftp.ensembl.org/pub/release-89/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_rm.chromosome.18.fa.gz
    wget ftp://ftp.ensembl.org/pub/release-89/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_rm.chromosome.19.fa.gz
    wget ftp://ftp.ensembl.org/pub/release-89/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_rm.chromosome.X.fa.gz
    wget ftp://ftp.ensembl.org/pub/release-89/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_rm.chromosome.Y.fa.gz
    wget ftp://ftp.ensembl.org/pub/release-89/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_rm.chromosome.MT.fa.gz

    # Download BK000964.3 rDNA repeat sequence from Dropbox
    wget https://www.dropbox.com/s/8lan5aqt1m47y52/rDNA_repeat_BK000964.3.fa

    # Unzip
    gunzip Mus_musculus.GRCm38.dna_rm.chromosome.*

    # Concatenate
    cat Mus_musculus.GRCm38.dna_rm.chromosome.1.fa \
    Mus_musculus.GRCm38.dna_rm.chromosome.2.fa \
    Mus_musculus.GRCm38.dna_rm.chromosome.3.fa \
    Mus_musculus.GRCm38.dna_rm.chromosome.4.fa \
    Mus_musculus.GRCm38.dna_rm.chromosome.5.fa \
    Mus_musculus.GRCm38.dna_rm.chromosome.6.fa \
    Mus_musculus.GRCm38.dna_rm.chromosome.7.fa \
    Mus_musculus.GRCm38.dna_rm.chromosome.8.fa \
    Mus_musculus.GRCm38.dna_rm.chromosome.9.fa \
    Mus_musculus.GRCm38.dna_rm.chromosome.10.fa \
    Mus_musculus.GRCm38.dna_rm.chromosome.11.fa \
    Mus_musculus.GRCm38.dna_rm.chromosome.12.fa \
    Mus_musculus.GRCm38.dna_rm.chromosome.13.fa \
    Mus_musculus.GRCm38.dna_rm.chromosome.14.fa \
    Mus_musculus.GRCm38.dna_rm.chromosome.15.fa \
    Mus_musculus.GRCm38.dna_rm.chromosome.16.fa \
    Mus_musculus.GRCm38.dna_rm.chromosome.17.fa \
    Mus_musculus.GRCm38.dna_rm.chromosome.18.fa \
    Mus_musculus.GRCm38.dna_rm.chromosome.19.fa \
    Mus_musculus.GRCm38.dna_rm.chromosome.X.fa \
    Mus_musculus.GRCm38.dna_rm.chromosome.Y.fa \
    Mus_musculus.GRCm38.dna_rm.chromosome.MT.fa \
    rDNA_repeat_BK000964.3.fa > GRCm38+rDNA_repeat.fa

    # Remove individual files
    rm -f Mus_musculus.GRCm38.dna_rm.chromosome.*.fa
    rm -f rDNA_repeat_BK000964.3.fa
fi

# Download gene annotation
if [ ! -f GRCm38+rDNA_repeat.gtf ]; then
    # Download gtf file
    wget ftp://ftp.ensembl.org/pub/release-89/gtf/mus_musculus/Mus_musculus.GRCm38.89.gtf.gz

    # Unzip
    gunzip Mus_musculus.GRCm38.89.gtf.gz

    # Remove comment lines
    sed -i '/^#/d' Mus_musculus.GRCm38.89.gtf

    # Download BK000964.3 rDNA repeat annotation from Dropbox
    wget https://www.dropbox.com/s/1gnk4lzo0rp49c9/rDNA_repeat_BK000964.3.gtf

    # Concatenate
    cat Mus_musculus.GRCm38.89.gtf \
    rDNA_repeat_BK000964.3.gtf > GRCm38+rDNA_repeat.gtf

    # Remove individual files
    rm -f Mus_musculus.GRCm38.89.gtf
    rm -f rDNA_repeat_BK000964.3.gtf
fi
