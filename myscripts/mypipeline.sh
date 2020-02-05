#!/bin/bash

# We'll use this as a wrapper to run our different scripts

#Path to my repo
myrepo="/users/l/k/lkshapir/Ecological-Genomics"

# My population
mypop="XDS"

#Directory to cleaned and paired reads

input="/data/project_data/RS_ExomeSeq/fastq/edge_fastq/pairedcleanreads/${mypop}"

#Directory to store outputs of our mapping

output="/data/project_data/RS_ExomeSeq/mapping"

# Run mapping.sh

source ./mapping.sh

#Run the post-processing myscripts

source ./process_bam.sh
