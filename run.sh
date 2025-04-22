#!/bin/bash

#Change NCORES depending on workload and resources
snakemake -s snakefile --cores 80 --printshellcmds --verbose --keep-going --use-conda

