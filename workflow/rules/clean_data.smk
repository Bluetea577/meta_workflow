localrules:
    cleandata_nor,

CLEAN_RUN = os.path.join(WORKDIR, config["path"]["clean_path"])

# def seq_se_or_paired(wildcards):
#     directory = SRA_RUN + {wildcards.sra_run}
#     fq_files = [f for f in os.listdir(directory) if f.endswith('.fastq.gz')]
#     if len(fq_files) == 2:
#         return {
#             'fq1': SRA_RUN + f"/{wildcards.sra_run}/{wildcards.sra_run}_1.fastq.gz",
#             'fq2': SRA_RUN + f"/{wildcards.sra_run}/{wildcards.sra_run}_2.fastq.gz"
#         }
#     else:
#         return{
#             'fq': SRA_RUN + f"/{wildcards.sra_run}/{wildcards.sra_run}.fastq.gz"
#         }

# def seq_fraction(wildcards):
#     directory = SRA_RUN + f"/{wildcards.sra_run}"
#     fq_files = [f for f in os.listdir(directory) if f.endswith('.fastq.gz')]
#     if len(fq_files) == 2:
#         fraction = ["_1", "_2"]
#     else:
#         fraction = [""]
#     return fraction

def kneaddata_input(input):
    Nfiles = len(input)

    if Nfiles == 1:
        out = f"--input {input[0]}"
    else:
        out = f"--input1 {input[0]} --input2 {input[1]} --reorder"

    return out

# ###### ??????
# def find_item(wildcards, input, output):
#     Nfiles = len(input)
#     outdir = os.path.dirname(output[0])
#     if Nfiles != 1:
#         out = f"""
#         cd {outdir}
#         rm -rf !({wildcards.sra_run}_1_kneaddata_paired_1.fastq | {wildcards.sra_run}_1_kneaddata_paired_2.fastq | {wildcards.sra_run}.tsv)
#         """
#     else:
#         out = f"""
#         cd {outdir}
#         rm -rf !({wildcards.sra_run}_kneaddata.fastq | {wildcards.sra_run}.tsv)
#         """
#
#     return out.strip()

def gzip_item(wildcards, input, output):
    Nfiles = len(input)
    outdir = os.path.dirname(output[0])
    if Nfiles == 1:
        out = f"""
        gzip -1 {outdir}/{wildcards.sra_run}_kneaddata.fastq
        """
    else:
        out = f"""
        gzip -1 {outdir}/{wildcards.sra_run}_1_kneaddata_paired_1.fastq
        gzip -1 {outdir}/{wildcards.sra_run}_1_kneaddata_paired_2.fastq
        """

    return out.strip()

def rename_item(wildcards, input, output):
    Nfiles = len(input)
    outdir = os.path.dirname(output[0])
    if Nfiles == 1:
        out = ""
    else:
        out = f"""
        mv {outdir}/{wildcards.sra_run}_1_kneaddata_paired_1.fastq.gz {outdir}/{wildcards.sra_run}_kneaddata_1.fastq.gz
        mv {outdir}/{wildcards.sra_run}_1_kneaddata_paired_2.fastq.gz {outdir}/{wildcards.sra_run}_kneaddata_2.fastq.gz
        """

    return out.strip()

rule kneaddata_process:
    input:
        expand(SRA_RUN + "/{{sra_run}}/{{sra_run}}{frac}.fastq.gz", frac = Sra_frac)
    output:
        temp(expand(CLEAN_RUN + "/{{sra_run}}_temp/{{sra_run}}_kneaddata{frac}.fastq.gz", frac = Sra_frac)),
    params:
        inputs = lambda wc, input: kneaddata_input(input),
        outdir=CLEAN_RUN + "/{sra_run}_temp",
        humanref=config["dbs"]["human"],
        trimopt=config["kneaddata_opt"]["trimmomatic"],
        gzipitem = lambda wc, input, output: gzip_item(wc,input,output),
        renameitem = lambda wc, input, output: rename_item(wc,input,output),
    log:
        "logs/cleandata/kneaddata/{sra_run}.log",
    benchmark:
        "logs/benchmarks/cleandata/kneaddata/{sra_run}.tsv",
    threads: 16
    resources:
        mem_mb=3500*16,
    conda:
        "../envs/kneaddata.yaml"
    shell:
        """
        export PATH=$CONDA_PREFIX/share/trimmomatic-0.39-2:$PATH
        
        kneaddata \
        {params.inputs} \
        -db {params.humanref} \
        -t {threads} \
        --output {params.outdir} \
        --trimmomatic-options {params.trimopt} \
        --remove-intermediate-output 2> {log}
        
        kneaddata_read_count_table \
        --input {params.outdir} \
        --output {params.outdir}/{wildcards.sra_run}.tsv
        
        {params.gzipitem}
        
        {params.renameitem}
        """

def mv_item(wildcards, input, output):
    Nfiles = len(input)
    inputdir = os.path.dirname(input[0])
    outdir = os.path.dirname(output[0])
    if Nfiles == 1:
        out = f"""
        mv {inputdir}/{wildcards.sra_run}_kneaddata.fastq.gz {outdir}
        mv {inputdir}/{wildcards.sra_run}.tsv {outdir}
        """
    else:
        out = f"""
        mv {inputdir}/{wildcards.sra_run}_kneaddata_1.fastq.gz {outdir}
        mv {inputdir}/{wildcards.sra_run}_kneaddata_2.fastq.gz {outdir}
        mv {inputdir}/{wildcards.sra_run}.tsv {outdir}
        """
    return out.strip()

rule cleandata_nor:
    input:
        expand(CLEAN_RUN + "/{{sra_run}}_temp/{{sra_run}}_kneaddata{frac}.fastq.gz",frac=Sra_frac),
    output:
        expand(CLEAN_RUN + "/{{sra_run}}/{{sra_run}}_kneaddata{frac}.fastq.gz",frac=Sra_frac)
    params:
        inputdir=CLEAN_RUN + "/{sra_run}_temp",
        mvitem=lambda wc, input, output: mv_item(wc, input, output)
    threads: 1
    shell:
        """
        {params.mvitem}
        
        rm -rf {params.inputdir}
        """


    # run:
    #     if 'fq1' in input and 'fq2' in input:
    #         shell(
    #             " export PATH=$CONDA_PREFIX/share/trimmomatic-0.39-2:$PATH "
    #             " ; "
    #             " kneaddata "
    #             " --input1 {input.fq1} "
    #             " --input2 {input.fq2} "
    #             " -db {params.humanref} "
    #             " -t {threads} "
    #             " --reorder "
    #             " --output {params.outdir} "
    #             " --trimmomatic-options {params.trimopt} "
    #             " --remove-intermediate-output 2> {log} "
    #             " ; "
    #             " kneaddata_read_count_table "
    #             " --input {params.outdir} "
    #             " --output {params.outdir}/{wildcards.sra_run}.tsv "
    #             " ; "
    #             " find {params.outdir} -type f "
    #             " ! -name '{wildcards.sra_run}_1_kneaddata_paired_1.fastq' "
    #             " ! -name '{wildcards.sra_run}_1_kneaddata_paired_2.fastq' "
    #             " ! -name '{wildcards.sra_run}.tsv' -exec rm -f {} + "
    #             " ; "
    #             " gzip -1 {params.outdir}/{wildcards.sra_run}_1_kneaddata_paired_1.fastq "
    #             " ; "
    #             " gzip -1 {params.outdir}/{wildcards.sra_run}_1_kneaddata_paired_2.fastq "
    #             " ; "
    #             " mv {params.outdir}/{wildcards.sra_run}_1_kneaddata_paired_1.fastq.gz {params.outdir}/{wildcards.sra_run}_kneaddata_1.fastq.gz "
    #             " ; "
    #             " mv {params.outdir}/{wildcards.sra_run}_1_kneaddata_paired_2.fastq.gz {params.outdir}/{wildcards.sra_run}_kneaddata_2.fastq.gz "
    #         )
    #     elif 'fq' in input:
    #         shell(
    #             " export PATH=$CONDA_PREFIX/share/trimmomatic-0.39-2:$PATH "
    #             " ; "
    #             " kneaddata "
    #             " --input {input.fq} "
    #             " -db {params.humanref} "
    #             " -t {threads} "
    #             " --output {params.outdir} "
    #             " --trimmomatic-options {params.trimopt} "
    #             " --remove-intermediate-output 2> {log} "
    #             " ; "
    #             " kneaddata_read_count_table "
    #             " --input {params.outdir} "
    #             " --output {params.outdir}/{wildcards.sra_run}.tsv "
    #             " ; "
    #             " find {params.outdir} -type f "
    #             " ! -name '{wildcards.sra_run}_kneaddata.fastq' "
    #             " ! -name '{wildcards.sra_run}.tsv' -exec rm -f {} + "
    #             " ; "
    #             " gzip -1 {params.outdir}/{wildcards.sra_run}_kneaddata.fastq "
    #         )