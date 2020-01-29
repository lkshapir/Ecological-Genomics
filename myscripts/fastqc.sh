#!/bin/bash
cd ~/Ecological-Genomics/myresults/

#I am creating a new directory to store my results
mkdir fastqc

#Performing fastqc on each file labeled with "XDS"

for file in /data/project_data/RS_ExomeSeq/fastq/edge_fastq/XDS*fastq.gz

do

fastqc ${file} -o fastqc/

done

#Escape insert mode and escape file: hit escape


