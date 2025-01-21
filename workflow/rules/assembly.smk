ASSE_RUN = os.path.join(WORKDIR, config["path"]["assembly_path"])

def megahit_input(input):
    Nfiles = len(input)

    if Nfiles == 1:
        out = f"--read {input[0]}"
    else:
        out = f"-1 {input[0]} -2 {input[1]} "

        if Nfiles == 3:
            out += f"--read {input[2]}"
    return out

rule run_megahit:
    input:
        expand(CLEAN_RUN + "/{{sra_run}}/{{sra_run}}_kneaddata{frac}.fastq.gz", frac = Sra_frac),
    output:
        temp(ASSE_RUN + "/megahit/{sra_run}/{sra_run}_prefilter.contigs.fa"),
    benchmark:
        "logs/benchmarks/assembly/megahit/{sra_run}.tsv",
    log:
        "logs/assembly/megahit/{sra_run}.log",
    params:
        min_count=config["megahit_min_count"],
        k_min=config["megahit_k_min"],
        k_max=config["megahit_k_max"],
        k_step=config["megahit_k_step"],
        merge_level=config["megahit_merge_level"],
        prune_level=config["megahit_prune_level"],
        low_local_ratio=config["megahit_low_local_ratio"],
        outdir=lambda wc, output: os.path.dirname(output[0]),
        inputs=lambda wc, input: megahit_input(input),
    conda:
        "../envs/megahit.yaml"
    threads: config["assembly_threads"]
    resources:
        mem_mb=config["assembly_memory"] * 1000,
        time_min=60 * config["assembly_runtime"],
    shell:
        """
        rm -rf {params.outdir}
        
        megahit \
        {params.inputs} \
        --num-cpu-threads {threads} \
        --k-min {params.k_min} \
        --k-max {params.k_max} \
        --k-step {params.k_step} \
        --out-dir {params.outdir} \
        --out-prefix {wildcards.sra_run}_prefilter \
        --min-count {params.min_count} \
        --merge-level {params.merge_level} \
        --prune-level {params.prune_level} \
        --low-local-ratio {params.low_local_ratio} \
        --memory {resources.mem_mb}000000 >> {log} 2>&1
        """

localrules: rename_megahit_output,

rule rename_megahit_output:
    input:
        ASSE_RUN + "/megahit/{sra_run}/{sra_run}_prefilter.contigs.fa",
    output:
        ASSE_RUN + "/megahit/{sra_run}/{sra_run}_raw_contigs.fasta",
    shell:
        "cp {input} {output}"

rule rename_contigs:
    input:
        ASSE_RUN + "/megahit/{sra_run}/{sra_run}_raw_contigs.fasta",
    output:
        fasta=temp(ASSE_RUN + "/megahit/{sra_run}/{sra_run}_prefilter_contigs.fasta"),
        mapping_table=ASSE_RUN + "/megahit/{sra_run}/old2new_contig_names.tsv",
    threads: config.get("simplejob_threads", 1)
    resources:
        mem_mb=config["simplejob_mem"] * 1000,
        time_min=60 * config["runtime"]["default"],
    log:
        "logs/assembly/megahit/{sra_run}_rename.log",
    params:
        minlength=config["minimum_contig_length"],
    conda:
        "../envs/fasta.yaml"
    script:
        "../scripts/rename_assembly.py"

localrules:
    finalize_contigs,

rule finalize_contigs:
    input:
        ASSE_RUN + "/megahit/{sra_run}/{sra_run}_prefilter_contigs.fasta",
    output:
        ASSE_RUN + "/megahit/{sra_run}/{sra_run}_final_contigs.fasta",
    threads: 1
    shell:
        "cp {input} {output}"


rule calculate_contigs_stats:
    input:
        ASSE_RUN + "/megahit/{sra_run}/{sra_run}_final_contigs.fasta",
    output:
        ASSE_RUN + "/megahit/{sra_run}/final_contig_stats.txt",
    conda:
        "../envs/assembly_rep.yaml"
    log:
        "logs/assembly/megahit/{sra_run}_calculate.log",
    threads: 1
    resources:
        mem_mb=1000,
    shell:
        "stats.sh in={input} format=3 out={output} &> {log}"

def get_align_input(input):
    Nfiles = len(input.query)

    if Nfiles == 1:
        out = f"-U {input.query[0]}"
    else:
        out = f"-1 {input.query[0]} -2 {input.query[1]} "
    return out

rule align_reads_to_contigs:
    input:
        query=expand(CLEAN_RUN + "/{{sra_run}}/{{sra_run}}_kneaddata{frac}.fastq.gz", frac = Sra_frac),
        target=ASSE_RUN + "/megahit/{sra_run}/{sra_run}_final_contigs.fasta",
    output:
        ASSE_RUN + "/megahit/{sra_run}/sequence_alignment/{sra_run}_sort.bam",
    params:
        ref_out=lambda wc, output: os.path.dirname(output[0]),
        input=lambda wc, input: get_align_input(input),
    benchmark:
        "logs/benchmarks/assembly/align/{sra_run}.tsv"
    log:
        "logs/assembly/megahit/{sra_run}.log",
    conda:
        "../envs/align.yaml",
    threads: 16
    resources:
        mem_mb=config["mem"] * 1000,
    shell:
        """
        echo "Starting index building..." 2> {log} 
        bowtie2-build \
        --threads {threads} \
        -f {input.target} \
        {params.ref_out}/index \
        2>> {log}
        
        echo "Starting alignment and sorting..." 2>> {log} 
        bowtie2 \
        --threads {threads} \
        -x {params.ref_out}/index \
        {params.input} \
        | samtools sort \
        -O bam \
        --threads {threads} \
        -o - > {output} \
        2>> {log}
        """

rule pileup_contigs_sample:
    input:
        fasta=ASSE_RUN + "/megahit/{sra_run}/{sra_run}_final_contigs.fasta",
        bam=ASSE_RUN + "/megahit/{sra_run}/sequence_alignment/{sra_run}_sort.bam",
        bai=ASSE_RUN + "/megahit/{sra_run}/sequence_alignment/{sra_run}_sort.bam.bai",
    output:
        covhist=ASSE_RUN + "/megahit/{sra_run}/contig_stats/postfilter_coverage_histogram.txt",
        covstats=ASSE_RUN + "/megahit/{sra_run}/contig_stats/postfilter_coverage_stats.txt",
        bincov=ASSE_RUN + "/megahit/{sra_run}/contig_stats/postfilter_coverage_binned.txt",
    params:
        pileup_secondary=(
            "t"
            if config.get("count_multi_mapped_reads", CONTIG_COUNT_MULTI_MAPPED_READS)
            else "f"
        ),
        minmapq=config["minimum_map_quality"],
    benchmark:
        "logs/benchmarks/assembly/calculate_coverage/pileup/{sra_run}.txt"
    log:
        "logs/assembly/calculate_coverage/{sra_run}_pilup_final_contigs.log",  # This log file is uesd for report
    conda:
        "../envs/assembly_rep.yaml"
    threads: 16
    resources:
        mem_mb=config["mem"] * 1000,
        java_mem=int(config["mem"] * JAVA_MEM_FRACTION),
    shell:
        "pileup.sh "
        " ref={input.fasta} "
        " in={input.bam} "
        " threads={threads} "
        " -Xmx{resources.java_mem}G "
        " out={output.covstats} "
        " hist={output.covhist} "
        " concise=t "
        " minmapq={params.minmapq} "
        " secondary={params.pileup_secondary} "
        " bincov={output.bincov} "
        " 2> {log} "

rule create_bam_index:
    input:
        "{file}.bam",
    output:
        "{file}.bam.bai",
    conda:
        "../envs/assembly_rep.yaml"
    log: 
        "logs/{file}.log"
    threads: 4
    resources:
        mem_mb=8000,
    shell:
        "samtools index {input} 2> {log}"

rule predict_genes:
    input:
        ASSE_RUN + "/megahit/{sra_run}/{sra_run}_final_contigs.fasta",
    output:
        fna=ASSE_RUN + "/predicted_genes/{sra_run}/{sra_run}.fna",
        faa=ASSE_RUN + "/predicted_genes/{sra_run}/{sra_run}.faa",
        gff=ASSE_RUN + "/predicted_genes/{sra_run}/{sra_run}.gff",
    conda:
        "../envs/prodigal.yaml"
    log:
        "logs/assembly/gene_annotation/{sra_run}_prodigal.log",
    benchmark:
        "logs/benchmarks/assembly/gene_annotation/{sra_run}_prodigal.txt",
    threads: 4
    resources:
        mem_mb=config["simplejob_mem"] * 1000,
        time_min=60 * config["runtime"]["simplejob"],
    shell:
        """
        prodigal -i {input} -o {output.gff} -d {output.fna} \
            -a {output.faa} -p meta -f gff 2> {log}
        """

localrules:
    get_contigs_from_gene_names,

rule get_contigs_from_gene_names:
    input:
        faa=ASSE_RUN + "/predicted_genes/{sra_run}/{sra_run}.faa",
    output:
        tsv=ASSE_RUN + "/predicted_genes/{sra_run}/{sra_run}.tsv",
    run:
        header = [
            "gene_id",
            "Contig",
            "Gene_nr",
            "Start",
            "Stop",
            "Strand",
            "Annotation",
        ]
        with open(output.tsv, "w") as tsv:
            tsv.write("\t".join(header) + "\n")
            with open(input.faa) as fin:
                gene_idx = 0
                for line in fin:
                    if line[0] == ">":
                        text = line[1:].strip().split(" # ")
                        old_gene_name = text[0]
                        text.remove(old_gene_name)
                        old_gene_name_split = old_gene_name.split("_")
                        gene_nr = old_gene_name_split[-1]
                        contig_nr = old_gene_name_split[-2]
                        sample = "_".join(
                            old_gene_name_split[: len(old_gene_name_split) - 2]
                        )
                        tsv.write(
                            "{gene_id}\t{sample}_{contig_nr}\t{gene_nr}\t{text}\n".format(
                                text="\t".join(text),
                                gene_id=old_gene_name,
                                i=gene_idx,
                                sample=sample,
                                gene_nr=gene_nr,
                                contig_nr=contig_nr,
                            )
                        )
                        gene_idx += 1

localrules:
    build_assembly_report,
    combine_contig_stats,

# log待调整
rule combine_contig_stats:
    input:
        contig_stats=expand(
            ASSE_RUN + "/megahit/{sra_run}/final_contig_stats.txt", sra_run=IDS
        ),
        gene_tables=expand(
            ASSE_RUN + "/predicted_genes/{sra_run}/{sra_run}.tsv", sra_run=IDS
        ),
        mapping_logs=expand(
            "logs/assembly/calculate_coverage/{sra_run}_pilup_final_contigs.log",
            sra_run=IDS,
        ),
        # mapping logs will be incomplete unless we wait on alignment to finish
        bams=expand(ASSE_RUN + "/megahit/{sra_run}/sequence_alignment/{sra_run}_sort.bam", sra_run=IDS),
    output:
        combined_contig_stats=WORKDIR + "stats/combined_contig_stats.tsv",
    params:
        samples=IDS,
    log:
        "logs/assembly/combine_contig_stats.log",
    script:
        "../scripts/combine_contig_stats.py"


# log待调整
rule build_assembly_report:
    input:
        combined_contig_stats=WORKDIR + "stats/combined_contig_stats.tsv",
    output:
        report=WORKDIR + "reports/assembly_report.html",
    conda:
        "../envs/report.yaml"
    log:
        "logs/assembly/report.log",
    script:
        "../scripts/assembly_report.py"


