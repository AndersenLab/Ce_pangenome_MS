#!/bin/bash

source activate trees 

# Run OrthoFinder on the core proteome set selected for genetic diversity
orthofinder -f ../../processed_data/orthology/orthofinder/proteomes_HC_CE/ -t 48

