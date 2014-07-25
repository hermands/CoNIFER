#!/bin/bash

#Pull down some RPKMs
dir="$HOME/repos/CoNIFER/targeted_data"
if [ ! -d $dir ]; then
	mkdir $dir
fi
cd $dir
rsync -v narwhal:/home/data/Conifer-Coseq21/no_dup_bed/RPKM/*.txt .
rsync -v narwhal:/home/genetics/Genomes/BED_Files/BROCA_v6_0465501.bed .
