import os
import random

configfile: "config.yaml"

datadir = config["location"]

if not os.path.isdir(datadir):
    sys.exit(f"\n❌ Data directory '{datadir}' not found! Please create it or check the path.\n")

# TBD
CHIP_DIRS = config["samples"]

# List all Native directories for each Chip
NATIVE_DIRS = {}
for chip in CHIP_DIRS:
    chip_dir=f"{datadir}/{chip}"
    subfolder = next(os.path.join(chip_dir, d) for d in os.listdir(chip_dir) if os.path.isdir(os.path.join(chip_dir, d)))
    fastq_dir = os.path.join(subfolder, "fastq_pass")
    NATIVE_DIRS[chip] = [name for name in os.listdir(fastq_dir) if os.path.isdir(os.path.join(fastq_dir, name)) and 'barcode' in name]


chip_native_pairs = [
    (chip, native)
    for chip in CHIP_DIRS
    for native in NATIVE_DIRS[chip]
]

# TBD
BLOCKS = {}
for c,n in chip_native_pairs:
    BLOCKS[(c,n)] = [block for block in range(config["prep"]["nblocks"])]

chip_native_block_trips = [
    (chip, native, block)
    for chip in CHIP_DIRS
    for native in NATIVE_DIRS[chip]
    for block in BLOCKS[(chip,native)]
]

chips, natives, blocks = zip(*chip_native_block_trips)
# Define rule all to specify final outputs for each chip
rule all:
    input:
        expand(
            "results/{chip}/{native}/{block}/fastq_files.lst",
            zip,
            chip=chips,
            native=natives,
            block=blocks
        ),
        expand(
            "results/{chip}/{native}/{block}/merged.fastq.gz",
            zip,
            chip=chips,
            native=natives,
            block=blocks
        ),
	expand(
            "results/{chip}/{native}/{block}/sequential_summary.csv",
            zip,
            chip=chips,
            native=natives,
            block=blocks
        ),
	expand(
            "results/{chip}/{native}/{block}/non_sequential_summary.csv",
            zip,
            chip=chips,
            native=natives,
            block=blocks
        )

# Rule to create the blocks with fastq_files.lst
rule create_blocks:
    output:
        "results/{chip}/{native}/{block}/fastq_files.lst"
    params:
        nsample=config["prep"]["nsample"],
        nblocks=config["prep"]["nblocks"],
        datadir=config["location"]
    shell:
        """
        realpath {params.datadir}/{wildcards.chip}/*/fastq_pass/{wildcards.native}/*.fastq.gz \
            | head -n {params.nsample} > {output}.temp

        nfiles=$(wc -l < {output}.temp)
        nfiles_per_block=$((nfiles / {params.nblocks}))
        nhead=$((nfiles_per_block * ({wildcards.block}+1)))

        head -n $nhead {output}.temp | tail -n $nfiles_per_block > {output}
        rm -f {output}.temp
        """


# Rule to merge the FASTQ files
rule merge_fastq:
    input:
        "results/{chip}/{native}/{block}/fastq_files.lst"
    output:
        "results/{chip}/{native}/{block}/merged.fastq.gz"
    shell:
        """
        for f in `cat {input}`; do
            cat "$f" >> {output}.temp;
        done
        mv {output}.temp {output}
        """


rule align:
    input:
        Templates=config["alignment"]["reference"],
	script=config["alignment"]["align"],
	fastq="results/{chip}/{native}/{block}/merged.fastq.gz"		
    output:
        Non_sequential_summary="results/{chip}/{native}/{block}/non_sequential_summary.csv",
        Sequential_summary="results/{chip}/{native}/{block}/sequential_summary.csv",
        Aligned_reads="results/{chip}/{native}/{block}/aligned.fastq.gz",
	Errors_out="results/{chip}/{native}/{block}/inserts.out"
    params:
        output_dir="results/{chip}/{native}/{block}",
	cut_off=config["alignment"]["cutoff"]
    shell:
        """
        python3 {input.script} --iq {input.fastq} --ia {input.Templates} --od {params.output_dir} --score_cut_off {params.cut_off} 
        """
