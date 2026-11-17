#!/bin/bash

source activate trees

### Run with: sbatch --array=1-5000%50 protein_msa.sh
set -euo pipefail

OG=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ../../processed_data/orthology/orthofinder/orthofinder_output/Orthogroups_SingleCopyOrthologues.txt)

IN="/vast/eande106/projects/Lance/THESIS_WORK/assemblies/orthology/elegans/orthofinder/64_core/OrthoFinder/Results_Dec07/Single_Copy_Orthologue_Sequences/${OG}.fa"
OUT="../../processed_data/divergence_estimates/single_copy_msa_orthogroups/${OG}.MSA.fa"

famsa -t 2 $IN $OUT
