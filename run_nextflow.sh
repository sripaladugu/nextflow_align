#!/usr/bin/env bash

#State the custom file working on
echo "Working on $1"

#Run nextflow
nextflow run -with-dag dag.png ./main.nf \
         -c $1 \
         -w /vol1/work/ \
         -resume -bg > .information
