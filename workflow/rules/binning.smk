BIN_RUN = os.path.join(WORKDIR, config["path"]["binning_path"])

# METABAT
rule get_metabat_depth:
    input:
        ASSE_RUN + "/megahit/{sra_run}/sequence_alignment/{sra_run}_sort.bam",
    output:
        temp(BIN_RUN + "/metabat/{sra_run}/metabat_depth.txt"),
    log:
        "logs/binning/matabat/{sra_run}_depth.log",
    benchmark:
        "logs/benchmarks/binning/matabat/{sra_run}_depth.tsv",
    conda:
        "../envs/metabat.yaml"
    threads: 16
    resources:
        mem_mb=config["mem"] * 1000
    params:
        minid = lambda wc, input: (
            config["cobinning_readmapping_id"] * 100 if len(input) > 1 else 97
        ),
    shell:
        """
        jgi_summarize_bam_contig_depths \
        --outputDepth {output} \
        {input} &> {log}
        """

def get_metabat_sensitivity():
    if config["metabat"]["sensitivity"] == "sensitive":
        return 500
    else:
        200

rule metabat:
    input:
        depth_file=BIN_RUN + "/metabat/{sra_run}/metabat_depth.txt",
        fasta=ASSE_RUN + "/megahit/{sra_run}/{sra_run}_final_contigs.fasta",
    output:
        clsfile=temp(BIN_RUN + "/metabat/{sra_run}/{sra_run}.tmp"),
    params:
        min_contig_len = config["metabat"]["min_contig_length"],
        sensitivity=get_metabat_sensitivity(),
    log:
        "logs/binning/metabat/{sra_run}_metabat.log",
    benchmark:
        "logs/benchmarks/binning/metabat/{sra_run}_matabat.tsv",
    conda:
        "../envs/metabat.yaml"
    threads: 16
    resources:
        mem_mb=config["mem"] * 1000
    shell:
        """
        metabat2 \
        -i {input.fasta} \
        --abdFile {input.depth_file} \
        --minContig {params.min_contig_len} \
        --numThreads {threads} \
        --maxEdges {params.sensitivity} \
        --saveCls \
        --noBinOut \
        -o {output.clsfile} \
        &> {log}
        """

rule get_unique_cluster_attribution:
    input:
        BIN_RUN + "/metabat/{sra_run}/{sra_run}.tmp"
    output:
        attribution=temp(BIN_RUN + "/metabat/{sra_run}/cluster_attribution.tsv"),
        status=BIN_RUN + "/metabat/{sra_run}/bin_status"
    log:
        "logs/binning/metabat/{sra_run}_cluster_attribution.log",
    run:
        import pandas as pd

        d = pd.read_csv(input[0], index_col=0, header=None, sep="\t").squeeze()

        assert isinstance(d, pd.Series), f"Expected two-column file: {input[0]}"

        old_cluster_ids = [id for id in d.unique() if id != 0]
        map_cluster_ids = {old: f"{wildcards.sra_run}_metabat_{i}"
                           for i, old in enumerate(old_cluster_ids, 1)}

        new_d = d.map(map_cluster_ids)
        new_d.dropna(inplace=True)

        if new_d.empty:
            logger.warning(f"No bins detected in sample {wildcards.sra_run}")
            new_d[wildcards.sra_run] = f"{wildcards.sra_run}_metabat_1"
            with open(output.status, 'w', encoding='utf-8') as f:
                f.write('empty')
        else:
            with open(output.status,'w',encoding='utf-8') as f:
                f.write('valid')

        new_d.to_csv(output[0], sep="\t", header=False)

rule get_bins:
    input:
        cluster_attribution=BIN_RUN + "/metabat/{sra_run}/cluster_attribution.tsv",
        contigs=ASSE_RUN + "/megahit/{sra_run}/{sra_run}_final_contigs.fasta",
        status=BIN_RUN + "/metabat/{sra_run}/bin_status",
    output:
        directory(BIN_RUN + "/metabat/{sra_run}/bin"),
    conda:
        "../envs/get_bins.yaml"
    log:
        "logs/binning/metabat/{sra_run}_generate_bins.log",
    script:
        "../scripts/get_fasta_of_bins.py"

"""
calculate bin quality
"""

rule run_checkm2:
    input:
        fasta_dir=BIN_RUN + "/metabat/{sra_run}/bin",
        db=rules.checkm2_download_db.output,
        status=BIN_RUN + "/metabat/{sra_run}/bin_status",
    output:
        table=BIN_RUN + "/metabat/{sra_run}/checkm2_report.tsv",
        faa=temp(directory(BIN_RUN + "/metabat/{sra_run}/faa")),
    params:
        lowmem=" --lowmem " if config["mem"] < 10 else "",
        dir=lambda wc, output: os.path.join(os.path.dirname(output.table), "checkm2"),
    conda:
        "../envs/checkm2.yaml"
    threads: 16
    log:
        "logs/binning/metabat/{sra_run}_checkm2.log",
    benchmark:
        "logs/benchmarks/binning/metabat/{sra_run}_checkm2.tsv",
    resources:
        time_min = 60 * int(config["runtime"]["default"]),
        mem_mb = config["mem"] * 1000,
    shell:
        """  
        if [ $(cat {input.status}) = "empty" ]; then  
            mkdir -p {output.faa}  
            echo -e "Name\tCompleteness\tContamination" > {output.table}  
            echo -e "{wildcards.sra_run}_metabat_1\t0\t0" >> {output.table}  
        else  
            checkm2 predict \
            --threads {threads} \
            {params.lowmem} \
            --force \
            --allmodels \
            -x .fasta \
            --input {input.fasta_dir} \
            --output-directory {params.dir} \
            &> {log}  
            
            cp {params.dir}/quality_report.tsv {output.table} 2>> {log}  
            mv {params.dir}/protein_files {output.faa} 2>> {log}  
            
            rm -rf {params.dir}  
        fi
        """

rule run_gunc:
    input:
        db=DBDIR + "/gunc_db/gunc_db_progenomes2.1.dmnd",
        fasta_dir=BIN_RUN + "/metabat/{sra_run}/faa",
        status=BIN_RUN + "/metabat/{sra_run}/bin_status",
    output:
        table=BIN_RUN + "/metabat/{sra_run}/gunc_report.tsv",
        folder=temp(directory(BIN_RUN + "/metabat/{sra_run}/gunc")),
    params:
        extension=".faa",
    conda:
        "../envs/gunc.yaml"
    threads: 16
    log:
        "logs/binning/metabat/{sra_run}_gunc.log",
    resources:
        time_min=60 * int(config["runtime"]["default"]),
        mem_mb=config["mem"] * 1000,
    shell:
        """
        if [ $(cat {input.status}) = "empty" ]; then  
            mkdir -p {output.folder}  
            echo -e "genome\tpass.GUNC\twarning" > {output.table}  
            echo -e "{wildcards.sra_run}_metabat_1\tFalse\tEmpty bin" >> {output.table}  
        else  
            mkdir -p {output.folder} 2> {log}  
            
            gunc run \
            --threads {threads} \
            --gene_calls \
            --db_file {input.db} \
            --input_dir {input.fasta_dir} \
            --file_suffix {params.extension} \
            --out_dir {output.folder} &>> {log}  
            
            cp {output.folder}/*.tsv {output.table} 2>> {log}  
        fi
        """

localrules:
    combine_checkm2,
    combine_gunc,

rule combine_checkm2:
    input:
        completeness_files = expand(
            BIN_RUN + "/metabat/{sra_run}/checkm2_report.tsv",
            sra_run=IDS,
        ),
        status_files=expand(
            BIN_RUN + "/metabat/{sra_run}/bin_status",
            sra_run=IDS,
        )
    output:
        bin_table=BIN_RUN + "/metabat/checkm2_quality_report.tsv",
    params:
        samples=IDS,
    log:
        "logs/binning/metabat/combine_checkm2.log",
    script:
        "../scripts/combine_checkm2.py"

rule combine_gunc:
    input:
        gunc_files=expand(
            BIN_RUN + "/metabat/{sra_run}/gunc_report.tsv",
            sra_run=IDS,
        ),
        status_files=expand(
            BIN_RUN + "/metabat/{sra_run}/bin_status",
            sra_run=IDS,
        )
    output:
        bin_table=BIN_RUN + "/metabat/gunc_quality_report.tsv",
    params:
        samples=IDS,
    log:
        "logs/binning/metabat/combine_gunc.log",
    script:
        "../scripts/combine_gunc.py"

localrules:
    get_bin_filenames,

rule get_bin_filenames:
    input:
        dirs=expand(
            BIN_RUN + "/metabat/{sra_run}/bin",
            sra_run=IDS,
        ),
        protein_dirs=expand(
            BIN_RUN + "/metabat/{sra_run}/faa",
            sra_run=IDS,
        ),
    log:
        "logs/binning/metabat/get_bin_filename.log",
    output:
        filenames=BIN_RUN + "/metabat/bins_paths.tsv",
    script:
        "../scripts/bin_filename.py"

localrules:
    all_contigs2bins,


rule all_contigs2bins:
    input:
        tsv_files=expand(
            BIN_RUN + "/metabat/{sra_run}/cluster_attribution.tsv",
            sra_run=IDS,
        ),
        status_files=expand(
            BIN_RUN + "/metabat/{sra_run}/bin_status",
            sra_run=IDS,
        )
    output:
        temp(BIN_RUN + "/metabat/contigs2bins.tsv.gz"),
    run:
        import gzip

        # 获取有效的文件列表
        valid_files = [
            tsv for tsv, status in zip(input.tsv_files,input.status_files)
            if open(status).read().strip() == 'valid'
        ]

        # 合并有效文件
        with gzip.open(output[0],'wb') as f_out:
            for file in valid_files:
                with open(file,'rb') as f_in:
                    f_out.write(f_in.read())

rule calculate_stats:
    input:
        bin_dir=BIN_RUN + "/metabat/{sra_run}/bin",
        status=BIN_RUN + "/metabat/{sra_run}/bin_status",
    output:
        BIN_RUN + "/metabat/{sra_run}/genome_stats.tsv",
    threads: 16
    params:
        extension=".fasta",
    log:
        "logs/binning/metabat/{sra_run}calculate_stats.log",
    script:
        "../scripts/calculate_bin_stats.py"


rule combine_bin_stats:
    input:
        stats=expand(
            BIN_RUN + "/metabat/{sra_run}/genome_stats.tsv",
            sra_run=IDS,
        ),
        status=expand(
            BIN_RUN + "/metabat/{sra_run}/bin_status",
            sra_run=IDS,
        ),
    output:
        BIN_RUN + "/metabat/genome_stats.tsv",
    params:
        samples=IDS,
    log:
        "logs/binning/metabat/combine_stats.log",
    script:
        "../scripts/combine_bin_stats.py"

def quality_filter_bins_input(wildcards):
    "Specify input files for quality_filter_bins rule"

    input_files = dict(
        paths=BIN_RUN + "/metabat/bins_paths.tsv",
        stats=BIN_RUN + "/metabat/genome_stats.tsv",
        quality=BIN_RUN + "/metabat/checkm2_quality_report.tsv",
        gunc=BIN_RUN + "/metabat/gunc_quality_report.tsv",
    )

    # check if gunc is in config file
    filter_chimieric_bins = config.get("filter_chimieric_bins", False)
    assert (
        isinstance(filter_chimieric_bins, bool)
    ), f"filter_chimieric_bins in config file must be a boolean, got {filter_chimieric_bins}"

    if not filter_chimieric_bins:
        del input_files["gunc"]

    return input_files


rule quality_filter_bins:
    input:
        unpack(quality_filter_bins_input),
    output:
        info=BIN_RUN + "/metabat/filtered/filtered_bins_info.tsv",
        paths=BIN_RUN + "/metabat/filtered/filtered_bins_paths.txt",
    threads: 1
    log:
        "logs/binning/metabat/filter_bins.log",
    params:
        filter_criteria=config["genome_filter_criteria"],
    script:
        "../scripts/filter_genomes.py"

localrules:
    upload_bins,
    upload_bin_report,
    finish_binning,

if config.get("upload", False):
    rule upload_bins:
        input:
            bin_dir=BIN_RUN + "/metabat/{sra_run}/bin",
            checkm2_report=BIN_RUN + "/metabat/{sra_run}/checkm2_report.tsv",
            gunc_report=BIN_RUN + "/metabat/{sra_run}/gunc_report.tsv" if config.get("filter_chimieric_bins",False) else [],
            genome_stats=BIN_RUN + "/metabat/{sra_run}/genome_stats.tsv",
            cluster_attribution=BIN_RUN + "/metabat/{sra_run}/cluster_attribution.tsv",
            faa_dir=BIN_RUN + "/metabat/{sra_run}/faa",
            status=BIN_RUN + "/metabat/{sra_run}/bin_status",
        output:
            mark=touch(BIN_RUN + "/metabat/{sra_run}/.{sra_run}.upload.done"),
        params:
            remote_dir=config.get("upload_tag","") + "binning/{sra_run}",
            config_dir=lambda wildcards: f"/tmp/bypy_bins_{wildcards.sra_run}"
        conda:
            "../envs/baiduyun.yaml"
        resources:
            upload_slots=1,
        log:
            "logs/binning/upload/{sra_run}.log"
        shell:
            """  
            mkdir -p {params.config_dir}  
            cp ~/.bypy/bypy.json {params.config_dir}/  
    
            bypy --config-dir {params.config_dir} mkdir {params.remote_dir} 2>> {log}  
            bypy --config-dir {params.config_dir} mkdir {params.remote_dir}/bin 2>> {log}  
            bypy --config-dir {params.config_dir} mkdir {params.remote_dir}/faa 2>> {log}  
    
            bypy --config-dir {params.config_dir} upload {input.bin_dir}/ {params.remote_dir}/bin/ 2>> {log}   
            bypy --config-dir {params.config_dir} upload {input.faa_dir}/ {params.remote_dir}/faa/ 2>> {log}  
            bypy --config-dir {params.config_dir} upload {input.checkm2_report} {params.remote_dir}/ 2>> {log}  
            bypy --config-dir {params.config_dir} upload {input.genome_stats} {params.remote_dir}/ 2>> {log}  
            bypy --config-dir {params.config_dir} upload {input.cluster_attribution} {params.remote_dir}/ 2>> {log}  
            bypy --config-dir {params.config_dir} upload {input.status} {params.remote_dir}/ 2>> {log}  
    
            if [ -f "{input.gunc_report}" ]; then  
                bypy --config-dir {params.config_dir} upload {input.gunc_report} {params.remote_dir}/ 2>> {log}  
            fi  
    
            rm -rf {params.config_dir}  
            """

    rule upload_bin_report:
        input:
            bins_paths=BIN_RUN + "/metabat/bins_paths.tsv",
            checkm2_report=BIN_RUN + "/metabat/checkm2_quality_report.tsv",
            gunc_report=BIN_RUN + "/metabat/gunc_quality_report.tsv" if config.get("filter_chimieric_bins",False) else [],
            genome_stats=BIN_RUN + "/metabat/genome_stats.tsv",
            contigs2bins=BIN_RUN + "/metabat/contigs2bins.tsv.gz",
            filtered_info=BIN_RUN + "/metabat/filtered/filtered_bins_info.tsv",
            filtered_paths=BIN_RUN + "/metabat/filtered/filtered_bins_paths.txt"
        output:
            mark=touch(BIN_RUN + "/.bin_report_upload.done")
        params:
            remote_dir=config.get("upload_tag","") + "binning/report",
            config_dir="/tmp/bypy_bins_report"
        conda:
            "../envs/baiduyun.yaml"
        resources:
            upload_slots=1,
        log:
            "logs/binning/upload_report.log"
        shell:
            """  
            mkdir -p {params.config_dir}  
            cp ~/.bypy/bypy.json {params.config_dir}/  
    
            bypy --config-dir {params.config_dir} mkdir {params.remote_dir} 2>> {log}  
            bypy --config-dir {params.config_dir} mkdir {params.remote_dir}/filtered 2>> {log}  
    
            bypy --config-dir {params.config_dir} upload {input.bins_paths} {params.remote_dir}/ 2>> {log}  
            bypy --config-dir {params.config_dir} upload {input.checkm2_report} {params.remote_dir}/ 2>> {log}  
            bypy --config-dir {params.config_dir} upload {input.genome_stats} {params.remote_dir}/ 2>> {log}  
            bypy --config-dir {params.config_dir} upload {input.contigs2bins} {params.remote_dir}/ 2>> {log}  
            bypy --config-dir {params.config_dir} upload {input.filtered_info} {params.remote_dir}/filtered/ 2>> {log}  
            bypy --config-dir {params.config_dir} upload {input.filtered_paths} {params.remote_dir}/filtered/ 2>> {log}  
    
            if [ -f "{input.gunc_report}" ]; then  
                bypy --config-dir {params.config_dir} upload {input.gunc_report} {params.remote_dir}/ 2>> {log}  
            fi  
    
            rm -rf {params.config_dir}  
            """

    rule finish_binning_with_upload:
        input:
            upload_done=expand(BIN_RUN + "/metabat/{sra_run}/.{sra_run}.upload.done",sra_run=IDS),
            report_upload_done=BIN_RUN + "/.bin_report_upload.done"
        output:
            touch(BIN_RUN + "/rule_binning.done")
else:
    rule finish_binning:
        input:
            bin_dir = expand(BIN_RUN + "/metabat/{sra_run}/bin", sra_run=IDS),
            checkm2_report = expand(BIN_RUN + "/metabat/{sra_run}/checkm2_report.tsv", sra_run=IDS),
            gunc_report = expand(BIN_RUN + "/metabat/{sra_run}/gunc_report.tsv", sra_run=IDS) if config.get("filter_chimieric_bins",False) else [],
            genome_stats = expand(BIN_RUN + "/metabat/{sra_run}/genome_stats.tsv", sra_run=IDS),
            cluster_attribution = expand(BIN_RUN + "/metabat/{sra_run}/cluster_attribution.tsv", sra_run=IDS),
            faa_dir = expand(BIN_RUN + "/metabat/{sra_run}/faa", sra_run=IDS),
            status = expand(BIN_RUN + "/metabat/{sra_run}/bin_status", sra_run=IDS),
            bins_paths=BIN_RUN + "/metabat/bins_paths.tsv",
            checkm2_reports=BIN_RUN + "/metabat/checkm2_quality_report.tsv",
            gunc_reports=BIN_RUN + "/metabat/gunc_quality_report.tsv" if config.get("filter_chimieric_bins",False) else [],
            genome_statses=BIN_RUN + "/metabat/genome_stats.tsv",
            contigs2bins=BIN_RUN + "/metabat/contigs2bins.tsv.gz",
            filtered_info=BIN_RUN + "/metabat/filtered/filtered_bins_info.tsv",
            filtered_paths=BIN_RUN + "/metabat/filtered/filtered_bins_paths.txt",
        output:
            touch(BIN_RUN + "/rule_binning.done")


# rule build_bin_report:
#     input:
#         checkm2_report=BIN_RUN + "metabat/checkm2_quality_report.tsv",
#     output:
#         report=WORKDIR + "reports/binning_report.html",
#     conda:
#         "../envs/report.yaml"
#     log:
#         "logs/binning/report.log",
#     script:
#         "../scripts/bin_report.py"
