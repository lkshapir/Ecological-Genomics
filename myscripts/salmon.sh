#!/bin/bash
cd /data/project_data/RS_RNASeq/fastq/cleanreads/

for file in ESC_01_C*R1.cl.fq
do
  salmon quant -i /data/project_data/RS_RNASeq/ReferenceTranscriptome/Pabies_cds_index \
   -l A \
   -r ${file} \
   --validateMappings \
   --seqBias \
   -o /data/project_data/RS_RNASeq/salmon/allmapping/${file}
done
