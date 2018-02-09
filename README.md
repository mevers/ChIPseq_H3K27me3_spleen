# ChIPseq_H3K27me3_spleen

This repository contains snakemake workflow rules and configuration files for the analysis of H3K27me3 ChIP-seq data from the [ENCODE Consortium](https://www.encodeproject.org/reference-epigenomes/ENCSR902FHL).

## Raw sequencing data

Raw sequencing data can be downloaded from the ENCODE Experiment summary website.

1. [H3K27me3 ChIP-seq data for replicate 1](https://www.encodeproject.org/files/ENCFF001KVG/)
2. [H3K27me3 ChIP-seq data for replicate 2](https://www.encodeproject.org/files/ENCFF001KVH/)
3. [Input data for replicate 1](https://www.encodeproject.org/files/ENCFF001KVR)
4. [Input data for replicate 2](https://www.encodeproject.org/files/ENCFF001KWE)

The gzipped FASTQ files must be placed in rolder `./rawData`.

## Requirements

1. Raw sequencing data must be downloaded manually and place in folder `./rawData`.

2. The full analysis was tested using the following program/package versions:

    * bedtools 2.25.0-24-g3d31735-dirty
    * bowtie2 and bowtie2-build 2.2.6
    * deeptools 2.4.1
    * FastQC 0.10.1
    * picard-tools 2.7.1-SNAPSHOT
    * samtools 1.3.1
    * snakemake 3.9.0


## Author

In case of questions please contact [Maurits Evers](mailto:maurits.evers@anu.edu.au).
