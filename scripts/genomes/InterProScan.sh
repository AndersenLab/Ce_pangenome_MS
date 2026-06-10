#!/bin/bash

base=$(basename $file)
strain=${base%%.*}

interproscan.sh \
	--formats TSV \
	--input $file \
	--goterms \
	--cpu 12 \
	--iprlookup \
	--disable-precalc \
	--output-file-base output/${strain}_allApps \
	--tempdir $TMP
