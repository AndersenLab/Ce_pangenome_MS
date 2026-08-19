#!/bin/bash

#SBATCH -J famsa_MSA
#SBATCH -A eande106          # Allocation name
#SBATCH -p parallel          # Partition/Queue name
#SBATCH -t 10:00:00          # Job walltime/duration (hh:mm:ss)
#SBATCH -N 1                 # Number of nodes
#SBATCH -c 2                 # Number of cores
#SBATCH --output=slurm_logs/famsa_%A_%a.out
#SBATCH --error=slurm_logs/famsa_%A_%a.err

source activate trees

### Run with: sbatch --array=1-5000%50 protein_msa.sh
set -euo pipefail

OG=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ../../processed_data/orthology/orthofinder/orthofinder_output/Orthogroups_SingleCopyOrthologues.txt)

IN="/vast/eande106/projects/Lance/THESIS_WORK/assemblies/orthology/elegans/orthofinder/64_core/OrthoFinder/Results_Dec07/Single_Copy_Orthologue_Sequences/${OG}.fa"
OUT="../../processed_data/divergence_estimates/single_copy_msa_orthogroups/${OG}.MSA.fa"

famsa -t 2 $IN $OUT
