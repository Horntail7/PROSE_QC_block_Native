#!/bin/bash

#Change NCORES (default = 40) depending on workload and resources
snakemake -s snakefile --cores 40 --printshellcmds --verbose --keep-going --use-conda

