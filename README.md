# Codes for DGRP plasticity eQTL paper

This repository contains various codes and data (summary level data for figures only, raw data in public repository). It is organized in three directories.

### analysisCodes

This directory contains scripts (bash, R, perl) used to perform the analysese. The analyses were performed in the order based on six bash scripts whose file names start with numbers (e.g. 1rnaSeq.bash, ... etc.). These bash scripts call other scripts that are also stored in this directory.

### figureCodes

This directory contains R scripts to make figures and tables in the paper.

### figureData

This directory contains data summarized from analyses as documented in the analysisCodes directory. These data were used to generate figures and tables. Two large files could not be hosted on github.com. However, the script (runNetwork.R) to generate them and the necessary input files are provided.

The figures and tables can be reproduced by running the Makefile contained in the figureCodes directory.
