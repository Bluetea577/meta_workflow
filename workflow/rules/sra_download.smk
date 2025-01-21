localrules:
    prefetch,

wildcard_constraints:
    sra_run="[S,E,D]RR[0-9]+",

SRA_RUN = os.path.join(WORKDIR, config["path"]["sra_path"])

rule prefetch:
    output:
        sra=temp(touch(SRA_RUN + "/{sra_run}/{sra_run}.downloaded"))
    params:
        outdir=SRA_RUN
    log:
        "logs/SRAdownload/prefetch/{sra_run}.log"
    benchmark:
        "logs/benchmarks/SRAdownload/prefetch/{sra_run}.tsv"
    threads: 1
    resources:
        mem_mb=1000,
    conda:
        "../envs/sra_download.yaml"
    shell:
        " mkdir -p {params.outdir} 2> {log} "
        " ; "
        " prefetch "
        " --output-directory {params.outdir} "
        " -X 999999999 "
        " --progress "
        " --log-level info "
        " {wildcards.sra_run} &>> {log} "
        " ; "
        " vdb-validate {params.outdir}/{wildcards.sra_run}/{wildcards.sra_run}.sra &>> {log} "

# def seq_fraction(wildcards):
#     seq_type = SAMPLES_INDEX.loc[SAMPLES_INDEX.iloc[:, 0] == wildcards.sra_run, SAMPLES_INDEX.columns[1]].values[0]
#     return ["_1", "_2"] if seq_type == 'Y' else [""]

PAIRED_END = SAMPLES_TABLE.columns.str.contains("PAIRED").any() or config.get(
    "interleaved_fastqs", False
)

Sra_frac = ["_1", "_2"] if PAIRED_END else [""]

rule extract_run:
    input:
        flag=rules.prefetch.output
    output:
        expand(
            SRA_RUN + "/{{sra_run}}/{{sra_run}}{frac}.fastq.gz",
            frac=Sra_frac
        ),
    params:
        outdir=SRA_RUN + "/{sra_run}",
        sra_file=SRA_RUN + "/{sra_run}/{sra_run}.sra",
    log:
        "logs/SRAdownload/fasterqdump/{sra_run}.log",
    benchmark:
        "logs/benchmarks/SRAdownload/fasterqdump/{sra_run}.tsv"
    threads: 4
    resources:
        mem_mb=1000
    conda:
        "../envs/sra_download.yaml"
    shell:
        " parallel-fastq-dump "
        " --threads {threads} "
        " --gzip --split-files "
        " --outdir {params.outdir} "
        " --tmpdir {resources.tmpdir} "
        " --skip-technical --split-3 "
        " -s {params.sra_file} &>> {log} "
        " ; "
        " rm -f {params.sra_file} 2>> {log} "

# checkpoint extract_run:
#     input:
#         flag=rules.prefetch.output,
#     output:
#         directory(SRA_RUN + "/{sra_run}"),
#     params:
#         sra_file=SRA_RUN + "/{sra_run}/{sra_run}.sra",
#     log:
#         "logs/SRAdownload/fasterqdump/{sra_run}.log",
#     benchmark:
#         "logs/benchmarks/SRAdownload/fasterqdump/{sra_run}.tsv"
#     threads: 4
#     resources:
#         mem_mb=1000
#     conda:
#         "../envs/sra_download.yaml"
#     shell:
#         " parallel-fastq-dump "
#         " --threads {threads} "
#         " --gzip --split-files "
#         " --outdir {output} "
#         " --tmpdir {resources.tmpdir} "
#         " --skip-technical --split-3 "
#         " -s {params.sra_file} &>> {log} "
#         " ; "
#         " rm -f {params.sra_file} 2>> {log} "
#
#
# def get_fastq_files(wildcards):
#     checkpoint_output = checkpoints.extract_run.get(**wildcards).output[0]
#     fastq_files = glob.glob(os.path.join(checkpoint_output,f"{wildcards.sra_run}*.fastq.gz"))
#
#     if not fastq_files:
#         raise ValueError(f"No fastq files found for {wildcards.sra_run}")
#
#     return sorted(fastq_files)