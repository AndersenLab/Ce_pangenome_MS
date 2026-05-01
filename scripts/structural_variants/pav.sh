#!/bin/bash

# Container execution of PAV - make sure to run with "--nt" Snakemake flag otherwise PAV might fail because of run competition issues due to intermediate files being deleted
## Follow PAV GitHub directions for setting up PAV directory structure and pull PAV container from Docker Hub with tag "latest"
singularity run --bind $analysis_dir:/analysis_dir --bind $pav_dir:/pav_dir $container_image -c 48 --latency-wait 900 --nt

    
