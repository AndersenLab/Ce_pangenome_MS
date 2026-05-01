#!/bin/bash

# Run OrthoFinder on the core proteome set selected for genetic diversity
orthofinder -f ../../processed_data/orthology/orthofinder/proteomes/core_64 -t 48

# Finish running OrthoFinder on the rest of the 78 C. elegans wild strains
orthofinder --assign ../../processed_data/orthology/orthofinder/proteomes/accessory_78 --core ../../processed_data/orthology/orthofinder/proteomes/core_64/OrthoFinder/Results_Dec06
