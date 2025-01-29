localrules:
    cleandata_nor,

CLEAN_RUN = os.path.join(WORKDIR, config["path"]["clean_path"])

# SRA_RUN = os.path.join(WORKDIR, config["path"]["sra_path"])
#
# PAIRED_END = SAMPLES_TABLE.columns.str.contains("PAIRED").any() or config.get(
#     "interleaved_fastqs", False
# )
#
# Sra_frac = ["_1", "_2"] if PAIRED_END else [""]

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

def hostdb(wildcards):
    base_dir = config["host"]["base_dir"]
    db_params = [f'-db {config["dbs"]["human"]}']

    if 'host' not in SAMPLES_TABLE.columns:
        return " ".join(db_params)

    host = SAMPLES_TABLE.loc[wildcards.sra_run, "host"]
    if pd.isna(host) or host in ["", "na", "none", None]:
        return " ".join(db_params)

    hosts = [h.strip().lower() for h in str(host).split(',')]

    host_paths = {
        # 水产
        'anchovy': config["host"]["animal"]["aquatic"]["anchovy"],
        'carp': config["host"]["animal"]["aquatic"]["carp"],
        'catfish': config["host"]["animal"]["aquatic"]["catfish"],
        'cod': config["host"]["animal"]["aquatic"]["cod"],
        'oyster': config["host"]["animal"]["aquatic"]["oyster"],
        'salmon': config["host"]["animal"]["aquatic"]["salmon"],
        'shrimp': config["host"]["animal"]["aquatic"]["shrimp"],

        # 哺乳动物
        'buffalo': config["host"]["animal"]["mammals"]["buffalo"],
        'cattle': config["host"]["animal"]["mammals"]["cattle"],
        'goat': config["host"]["animal"]["mammals"]["goat"],
        'horse': config["host"]["animal"]["mammals"]["horse"],
        'pig': config["host"]["animal"]["mammals"]["pig"],
        'sheep': config["host"]["animal"]["mammals"]["sheep"],

        # 家禽
        'chicken': config["host"]["animal"]["poultry"]["chicken"],

        # 谷物
        'barley': config["host"]["plant"]["grains"]["barley"],
        'cowpea': config["host"]["plant"]["grains"]["cowpea"],
        'rice': config["host"]["plant"]["grains"]["rice"],
        'sorghum': config["host"]["plant"]["grains"]["sorghum"],
        'soybean': config["host"]["plant"]["grains"]["soybean"],
        'wheat': config["host"]["plant"]["grains"]["wheat"],

        # 香料
        'chili': config["host"]["plant"]["spices"]["chili"],
        'cocoa': config["host"]["plant"]["spices"]["cocoa"],
        'garlic': config["host"]["plant"]["spices"]["garlic"],
        'mustard': config["host"]["plant"]["spices"]["mustard"],
        'tea': config["host"]["plant"]["spices"]["tea"],

        # 蔬菜
        'cabbage': config["host"]["plant"]["vegetables"]["cabbage"],
        'carrot': config["host"]["plant"]["vegetables"]["carrot"],
        'cucumber': config["host"]["plant"]["vegetables"]["cucumber"],
        'grape': config["host"]["plant"]["vegetables"]["grape"],
        'radish': config["host"]["plant"]["vegetables"]["radish"],
        'spinach': config["host"]["plant"]["vegetables"]["spinach"],
        'tomato': config["host"]["plant"]["vegetables"]["tomato"]
    }

    for host in hosts:
        if host in host_paths:
            db_params.append(f'-db {os.path.join(base_dir,host_paths[host])}')
        elif host not in ["", "na", "none", None]:
            logger.warning(f"Unknown host species '{host}' for sample {wildcards.sra_run}")

    return " ".join(db_params)

rule kneaddata_process:
    input:
        fastqs=expand(SRA_RUN + "/{{sra_run}}/{{sra_run}}{frac}.fastq.gz", frac = Sra_frac),
        sra_done=SRA_RUN + "/{sra_run}/.{sra_run}.sra_download.done",
    output:
        temp(expand(CLEAN_RUN + "/{{sra_run}}_temp/{{sra_run}}_kneaddata{frac}.fastq.gz", frac = Sra_frac)),
    params:
        inputs = lambda wc, input: kneaddata_input(input.fastqs),
        outdir=CLEAN_RUN + "/{sra_run}_temp",
        trimopt=config["kneaddata_opt"]["trimmomatic"],
        gzipitem = lambda wc, input, output: gzip_item(wc,input.fastqs,output),
        renameitem = lambda wc, input, output: rename_item(wc,input.fastqs,output),
        dbs=lambda wc: hostdb(wc),
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
        {params.dbs} \
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
        temp(expand(CLEAN_RUN + "/{{sra_run}}/{{sra_run}}_kneaddata{frac}.fastq.gz",frac=Sra_frac)),
        CLEAN_RUN+ "/{sra_run}/{sra_run}.tsv",
        mark=touch(CLEAN_RUN + "/{sra_run}/.{sra_run}.clean_data.done"),  # 添加完成标记
    params:
        inputdir=CLEAN_RUN + "/{sra_run}_temp",
        mvitem=lambda wc, input, output: mv_item(wc, input, output)
    threads: 1
    shell:
        """
        {params.mvitem}
        
        rm -rf {params.inputdir}
        """

localrules:
    upload_clean_data,
    finish_clean_data,

if config.get("upload", False):
    rule upload_clean_data:
        input:
            fastqs = expand(CLEAN_RUN + "/{{sra_run}}/{{sra_run}}_kneaddata{frac}.fastq.gz",frac=Sra_frac),
            stats = CLEAN_RUN + "/{sra_run}/{sra_run}.tsv",
            clean_done = CLEAN_RUN + "/{sra_run}/.{sra_run}.clean_data.done",
        output:
            mark=touch(CLEAN_RUN + "/{sra_run}/.{sra_run}.clean_data.uploaded")
        params:
            remote_dir=config.get("upload_tag", "") + "clean_data/{sra_run}",
            sra_dir=SRA_RUN + "/{sra_run}",
            config_dir=lambda wildcards: f"/tmp/bypy_clean_{wildcards.sra_run}",
        log:
            "logs/cleandata/upload/{sra_run}.log"
        threads: 1
        conda:
            "../envs/baiduyun.yaml"
        resources:
            upload_slots=1,
        shell:
            """  
            mkdir -p {params.config_dir}  
            cp ~/.bypy/bypy.json {params.config_dir}/
            
            bypy --config-dir {params.config_dir} mkdir {params.remote_dir} 2>> {log}  
            
            for f in {input.fastqs}; do  
                bypy --config-dir {params.config_dir} upload "$f" {params.remote_dir}/ 2>> {log}   
            done 
        
            bypy --config-dir {params.config_dir} upload {input.stats} {params.remote_dir}/ 2>> {log}    
            
            rm -rf {params.config_dir}
            
            # rm -rf {params.sra_dir}/*.fastq.gz 2>> {log} 
            """

    rule finish_clean_data_with_upload:
        input:
            upload_marks=expand(CLEAN_RUN + "/{sra_run}/.{sra_run}.clean_data.uploaded",sra_run=IDS)
        output:
            touch(CLEAN_RUN + "/rule_clean_data.done")
else:
    rule finish_clean_data:
        input:
            fastqs = expand(CLEAN_RUN + "/{{sra_run}}/{{sra_run}}_kneaddata{frac}.fastq.gz",frac=Sra_frac),
            stats = expand(CLEAN_RUN + "/{sra_run}/{sra_run}.tsv", sra_run=IDS),
            clean_done = expand(CLEAN_RUN + "/{sra_run}/.{sra_run}.clean_data.done", sra_run=IDS),
        output:
            touch(CLEAN_RUN + "/rule_clean_data.done")


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