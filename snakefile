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
    NATIVE_DIRS[chip] = [
    name for name in os.listdir(fastq_dir)
    if os.path.isdir(os.path.join(fastq_dir, name))
    and 'barcode' in name
    and len(os.listdir(os.path.join(fastq_dir, name))) >= config["prep"]["nblocks"]  # ensures the directory is not empty
    ]


chip_native_pairs = [
    (chip, native)
    for chip in CHIP_DIRS
    for native in NATIVE_DIRS[chip]
]

chps, ntvs = zip(*chip_native_pairs)

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
        ),
	expand(
            "AGGRE/{chip}/{native}/non_sequential_summary.csv",
	    zip,
	    chip=chps,
	    native=ntvs
        ),
	expand(
            "AGGRE/{chip}/{native}/sequential_summary.csv",
	    zip,
	    chip=chps,
	    native=ntvs
        ),
	expand(
            "AGGRE/{chip}/{native}/aligned.fastq.gz",
	    zip,
	    chip=chps,
	    native=ntvs
        ),
	expand(
            "AGGRE/{chip}/{native}/inserts.out",
	    zip,
	    chip=chps,
	    native=ntvs
        )
	

# Rule to create the blocks with fastq_files.lst
rule create_blocks:
    output:
        "results/{chip}/{native}/{block}/fastq_files.lst"
    params:
        nsample=config["prep"]["nsample"],
        nblocks=config["prep"]["nblocks"],
        datadir=config["location"]
    resources:
        mem_mb=4000
    run:
        import os
        import glob

        pattern = os.path.join(params.datadir, wildcards.chip, "*", "fastq_pass", wildcards.native, "*.fastq.gz")
        files = sorted(glob.glob(pattern))

        files = files[:params.nsample]

        nfiles = len(files)
        nfiles_per_block = nfiles // params.nblocks

        block = int(wildcards.block)
        start = block * nfiles_per_block
        end = start + nfiles_per_block

        selected_files = files[start:end]

        with open(output[0], "w") as out_f:
            out_f.write("\n".join(os.path.realpath(f) for f in selected_files))


# Rule to merge the FASTQ files
rule merge_fastq:
    input:
        "results/{chip}/{native}/{block}/fastq_files.lst"
    output:
        "results/{chip}/{native}/{block}/merged.fastq.gz"
    resources:
        mem_mb=4000
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
	Inserts_out="results/{chip}/{native}/{block}/inserts.out"
    params:
        output_dir="results/{chip}/{native}/{block}",
	cut_off=config["alignment"]["cutoff"]
    resources:
        mem_mb=4000
    shell:
        """
        python3 {input.script} --iq {input.fastq} --ia {input.Templates} --od {params.output_dir} --score_cut_off {params.cut_off} 
        """

def get_nonseq_summary_inputs(wildcards):
    try:
        chip = wildcards.chip
        native = wildcards.native
        nblocks = int(config["prep"]["nblocks"])
        return [f"results/{chip}/{native}/{b}/non_sequential_summary.csv" for b in range(nblocks)]
    except Exception as e:
        print(f"⚠️ Error in get_nonseq_summary_inputs: {e}")
        return []  # <-- this would cause Snakemake failure if not handled

def get_seq_summary_inputs(wildcards):
    try:
        chip = wildcards.chip
        native = wildcards.native
        nblocks = int(config["prep"]["nblocks"])
        return [f"results/{chip}/{native}/{b}/sequential_summary.csv" for b in range(nblocks)]
    except Exception as e:
        print(f"⚠️ Error in get_seq_summary_inputs: {e}")
        return []  # <-- this would cause Snakemake failure if not handled

def get_aligned_inputs(wildcards):
    try:
        chip = wildcards.chip
        native = wildcards.native
        nblocks = int(config["prep"]["nblocks"])
        return [f"results/{chip}/{native}/{b}/aligned.fastq.gz" for b in range(nblocks)]
    except Exception as e:
        print(f"⚠️ Error in get_aligned_inputs: {e}")
        return []  # <-- this would cause Snakemake failure if not handled

def get_inserts_inputs(wildcards):
    try:
        chip = wildcards.chip
        native = wildcards.native
        nblocks = int(config["prep"]["nblocks"])
        return [f"results/{chip}/{native}/{b}/inserts.out" for b in range(nblocks)]
    except Exception as e:
        print(f"⚠️ Error in get_inserts_inputs: {e}")
        return []  # <-- this would cause Snakemake failure if not handled



# Rule to aggregate the block data
rule aggregate_data_summary:
    input:
        non_sequential_summary=get_nonseq_summary_inputs,
        sequential_summary=get_seq_summary_inputs,
        aligned=get_aligned_inputs,
	inserts=get_inserts_inputs
    output:
        Aggregate_non_sequential_summary="AGGRE/{chip}/{native}/non_sequential_summary.csv",
        Aggregate_sequential_summary="AGGRE/{chip}/{native}/sequential_summary.csv",
	Aggregate_aligned="AGGRE/{chip}/{native}/aligned.fastq.gz",
	Aggregate_inserts="AGGRE/{chip}/{native}/inserts.out"
    params:
    	templates=config["alignment"]["reference"],
    	nblocks=config["prep"]["nblocks"],
    	datadir=config["location"]
    resources:
        mem_mb=4000
    shell:
        """
    	echo "Seq,Read_count,Percentage" > {output.Aggregate_non_sequential_summary}
    	echo "Seq,Read_count,Percentage" > {output.Aggregate_sequential_summary}

    	total=0
    	for f in {input.sequential_summary}; do
            c=$(grep Total "$f" | sed "s|,| |g" | awk '{{print $2}}')
            total=$(echo $total $c | awk '{{print $1+$2}}')
    	done

    	nlen=0
    	for f in {input.non_sequential_summary}; do
            c=$(grep Reads_not_meeting_length_criteria "$f" | sed "s|,| |g" | awk '{{print $2}}')
            nlen=$(echo $nlen $c | awk '{{print $1+$2}}')
    	done

    	nunmatched=0
    	for f in {input.non_sequential_summary}; do
            c=$(grep Reads_unmatched "$f" | sed "s|,| |g" | awk '{{print $2}}')
            nunmatched=$(echo $nunmatched $c | awk '{{print $1+$2}}')
    	done

    	for temp in $(grep '^>' {params.templates} | sed 's/>//' ); do
            count=0
            for f in {input.non_sequential_summary}; do
            	c=$(grep "$temp" "$f" | sed "s|,| |g" | awk '{{print $2}}')
            	count=$(echo $count $c | awk '{{print $1+$2}}')
            done
            printf "%s,%d,%6.2f\\n" "$temp" "$count" "$(awk -v a=$count -v b=$total 'BEGIN{{print (a/b)*100}}')" >> {output.Aggregate_non_sequential_summary}

            count=0
            for f in {input.sequential_summary}; do
            	c=$(grep "$temp" "$f" | sed "s|,| |g" | awk '{{print $2}}')
            	count=$(echo $count $c | awk '{{print $1+$2}}')
            done
            printf "%s,%d,%6.2f\\n" "$temp" "$count" "$(awk -v a=$count -v b=$total 'BEGIN{{print (a/b)*100}}')" >> {output.Aggregate_sequential_summary}
        done

    	printf "Reads_not_meeting_length_criteria,%d,%6.2f\\n" "$nlen" "$(awk -v a=$nlen -v b=$total 'BEGIN{{print (a/b)*100}}')" >> {output.Aggregate_non_sequential_summary}
    	printf "Reads_unmatched,%d,%6.2f\\n" "$nunmatched" "$(awk -v a=$nunmatched -v b=$total 'BEGIN{{print (a/b)*100}}')" >> {output.Aggregate_non_sequential_summary}
    	printf "Total,%d,100.00\\n" "$total" >> {output.Aggregate_sequential_summary}

	for f in {input.aligned}; do
            cat "$f" >> {output.Aggregate_aligned}.temp;
        done
        mv {output.Aggregate_aligned}.temp {output.Aggregate_aligned}

	for f in {input.inserts}; do
            cat "$f" >> {output.Aggregate_inserts}.temp;
        done
        mv {output.Aggregate_inserts}.temp {output.Aggregate_inserts}
       

   	"""
	
