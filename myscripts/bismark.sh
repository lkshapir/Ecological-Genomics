
bismark -bowtie -multicore 1 \
--genome /data/project_data/epigenetics/reference_genome \
--output_dir /data/project_data/epigenetics/aligned_output \
-1 /data/project_data/epigenetics/trimmed_fastq/HA_F25_3_1.fq.gz \
-2 /data/project_data/epigenetics/trimmed_fastq/HA_F25_3_2.fq.gz \
-rg_tag --rg_id HA_F25_3 --rg_sample HA_F25_3 --gzip --local --maxins 1000

#Last line adds read group to .bam file (adds sample ID), gzip compresses temporary
#files, local tells bowtie 2(which does the alignment) to use soft clippings
#to increase mapping rates, maxins is max insert size
