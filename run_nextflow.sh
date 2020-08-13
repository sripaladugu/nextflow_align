#!/usr/bin/env bash

#State the custom file working on
echo "Working on $1"

#Load Singularity
module load singularity

#Run nextflow
nextflow run -with-dag dag.png ./main.nf -c $1 -resume -bg > .information
