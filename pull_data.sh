#!/bin/bash

#Pull down some RPKMs
cd ~/repos/CoNIFER/targeted_data
rsync -v hermands@narwhal:/home/data/Conifer-Coseq21/no_dup_bed/RPKM/*.txt .
rsync -v hermands@narwhal:/home/genetics/Genomes/BED_Files/BROCA_v6_0465501.bed .
