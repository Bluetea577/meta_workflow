BIN_RUN = os.path.join(WORKDIR, config["path"]["binning_path"]) # 单独使用增加

rule quality_file:
    input:
        BIN_RUN + "/metabat/filtered/filtered_bins_info" + config.get("samples_batch", "") +".tsv",
    output:
        temp(BIN_RUN + "/metabat/filtered/filtered_bins_info_for_drep" + config.get("samples_batch", "") + ".csv"),
    threads: 1
    run:
        import pandas as pd
        import sys

        def convert_to_drep_format(input_file, output_file):
            df = pd.read_csv(input_file, sep="\t")

            drep_format = pd.DataFrame({
                'genome': df['Bin Id'] + '.fasta',
                'completeness': df['Completeness'],
                'contamination': df['Contamination']
            })

            drep_format.to_csv(output_file, index=False)

        convert_to_drep_format(input[0], output[0])

# use drep
rule run_drep:
    input:
        paths=BIN_RUN + "/metabat/filtered/filtered_bins_paths" + config.get("samples_batch", "") + ".txt",
        quality=BIN_RUN + "/metabat/filtered/filtered_bins_info_for_drep" + config.get("samples_batch", "") + ".csv"
    output:
        workdir=directory(BIN_RUN + "/derep/workdir" + config.get("samples_batch", "")),
        genome_dir=directory(BIN_RUN + "/derep/workdir" + config.get("samples_batch", "") + "/dereplicated_genomes"),
    log:
        "logs/binning/drep/drep" + config.get("samples_batch", "") + ".log"
    params:
        min_comp=config["drep"].get("min_completeness",75),
        max_cont=config["drep"].get("max_contamination",25),
        min_ani=config["drep"].get("primary_ani",0.9),
        final_ani=config["drep"].get("final_ani",0.95),
        min_cov=config["drep"].get("coverage_threshold",0.1),
        extra=config["drep"].get("extra","")
    threads: 16
    conda:
        "../envs/drep.yaml"
    resources:
        mem_mb=config["mem"] * 1000,
        time_min=60 * int(config["runtime"]["default"])
    shell:
        """  
        dRep dereplicate \
        {output.workdir} \
        --genomes {input.paths} \
        --processors {threads} \
        --genomeInfo {input.quality} \
        --completeness {params.min_comp} \
        --contamination {params.max_cont} \
        --P_ani {params.min_ani} \
        --S_ani {params.final_ani} \
        --cov_thresh {params.min_cov} \
        {params.extra} \
        &> {log}
        """

rule process_drep_results:
    input:
        workdir=BIN_RUN + "/derep/workdir" + config.get("samples_batch", ""),
        quality=BIN_RUN + "/metabat/filtered/filtered_bins_info" + config.get("samples_batch", "") +".tsv"
    output:
        merged_info=BIN_RUN + "/derep/bin_info" + config.get("samples_batch", "") + ".tsv",
        bins2species=BIN_RUN + "/derep/bins2species" + config.get("samples_batch", "") + ".tsv",
    log:
        "logs/binning/drep/process_drep_results" + config.get("samples_batch", "") + ".log"
    conda:
        "../envs/drep.yaml"
    script:
        "../scripts/process_drep_results.py"

"""
# skani
rule run_skani:
    input:
        paths=BIN_RUN + "/metabat/filtered/filtered_bins_paths.txt",
    output:
        temp("tmp/derep/distance_matrix.txt"),
    log:
        "logs/binning/drep/skani_calculation.log",
    resources:
        mem_mb=config["mem"] * 1000,
        time_min=60 * config["runtime"]["default"],
    params:
        #preset= "medium", # fast, medium or slow
        min_af=config["genome_dereplication"]["overlap"] * 100,
        extra="",
    threads: config["threads"]
    conda:
        "../envs/skani.yaml"
    shell:
        "skani triangle "
        " {params.extra} "
        " -l {input.paths} "
        " -o {output} "
        " -t {threads} "
        " --sparse --ci "
        " --min-af {params.min_af} "
        " &> {log} "


rule skani_2_parquet:
    input:
        rules.run_skani.output,
    output:
        BIN_RUN + "/derep/genome_similarities.parquet",
    resources:
        mem_mb=config["mem"] * 1000,
        time_min=60 * config["runtime"]["simplejob"],
    log:
        "logs/binning/drep/skani_2_parquet.log",
    threads: 1
    run:
        try:
            skani_column_dtypes = {
                "Ref_file": "category",
                "Query_file": "category",
                "ANI": float,
                "Align_fraction_ref": float,
                "Align_fraction_query": float,
                "ANI_5_percentile": float,
                "ANI_95_percentile": float,
            }  # Ref_name        Query_name

            import pandas as pd

            import pandas as pd

            df = pd.read_table(input[0])

            def simplify_path(path, remove_gz=True):
                path = Path(path)
                name = path.stem
                ext = path.suffix
                if remove_gz & (ext == ".gz"):
                    name = Path(name).stem
                return name

            df = pd.read_table(
                input[0],
                usecols=list(skani_column_dtypes.keys()),
                dtype=skani_column_dtypes,
            )

            df["Ref"] = df.Ref_file.cat.rename_categories(simplify_path)
            df["Query"] = df.Query_file.cat.rename_categories(simplify_path)

            df.to_parquet(output[0])

        except Exception as e:
            import traceback

            with open(log[0], "w") as logfile:
                traceback.print_exc(file=logfile)

            raise e


rule cluster_species:
    input:
        dist=BIN_RUN + "/derep/genome_similarities.parquet",
        bin_info=BIN_RUN + "/metabat/filtered/filtered_bins_info.tsv",
    output:  
        bin_info=BIN_RUN + "/derep/bin_info.tsv",  
        bins2species=BIN_RUN + "/derep/bins2species.tsv",
    params:
        linkage_method="average",
        pre_cluster_threshold=0.925,
        threshold=config["genome_dereplication"]["ANI"],
    conda:
        "../envs/species_clustering.yaml"
    log:
        "logs/binning/drep/species_clustering.log",
    script:
        "../scripts/cluster_species.py"

"""

localrules:
    build_bin_report,

rule build_bin_report:
    input:
        bin_info=BIN_RUN + "/derep/bin_info" + config.get("samples_batch", "") + ".tsv",
        bins2species=BIN_RUN + "/derep/bins2species" + config.get("samples_batch", "") + ".tsv",
    output:
        report=WORKDIR + "reports/bin_report_metabat" + config.get("samples_batch", "") + ".html",
    conda:
        "../envs/report.yaml"
    log:
        "logs/binning/report_binning" + config.get("samples_batch", "") + ".log",
    script:
        "../scripts/bin_report.py"

localrules:
    upload_drep,
    upload_drep_report,
    finish_drep,

if config.get("upload", False):
    rule upload_drep:
        input:
            workdir=BIN_RUN + "/derep/workdir" + config.get("samples_batch", ""),
            bin_info=BIN_RUN + "/derep/bin_info" + config.get("samples_batch", "") + ".tsv",
            bins2species=BIN_RUN + "/derep/bins2species" + config.get("samples_batch", "") + ".tsv",
        output:
            mark=touch(BIN_RUN + "/derep/.drep.upload" + config.get("samples_batch", "") + ".done")
        params:
            remote_dir=config.get("upload_tag","") + "binning/derep",
            config_dir="/tmp/bypy_drep"
        conda:
            "../envs/baiduyun.yaml"
        resources:
            upload_slots=1,
        log:
            "logs/binning/upload/drep" + config.get("samples_batch", "") + ".log"
        shell:
            """  
            mkdir -p {params.config_dir}  
            cp ~/.bypy/bypy.json {params.config_dir}/  
    
            bypy --config-dir {params.config_dir} mkdir {params.remote_dir} 2>> {log}   
            bypy --config-dir {params.config_dir} mkdir {params.remote_dir}/workdir 2>> {log}  
    
            bypy --config-dir {params.config_dir} upload {input.workdir}/ {params.remote_dir}/workdir/ 2>> {log}  
            bypy --config-dir {params.config_dir} upload {input.bin_info} {params.remote_dir}/ 2>> {log}  
            bypy --config-dir {params.config_dir} upload {input.bins2species} {params.remote_dir}/ 2>> {log}  
    
            rm -rf {params.config_dir}  
            """

    rule upload_drep_report:
        input:
            report=WORKDIR + "reports/bin_report_metabat" + config.get("samples_batch", "") + ".html"
        output:
            mark=touch(BIN_RUN + "/.drep_report_upload" + config.get("samples_batch", "") + ".done")
        params:
            remote_dir=config.get("upload_tag","") + "binning/report",
            config_dir="/tmp/bypy_drep_report"
        conda:
            "../envs/baiduyun.yaml"
        resources:
            upload_slots=1,
        log:
            "logs/binning/upload/drep_report" + config.get("samples_batch", "") + ".log"
        shell:
            """  
            mkdir -p {params.config_dir}  
            cp ~/.bypy/bypy.json {params.config_dir}/  
    
            bypy --config-dir {params.config_dir} mkdir {params.remote_dir} 2>> {log}  
            bypy --config-dir {params.config_dir} upload {input.report} {params.remote_dir}/ 2>> {log}  
    
            rm -rf {params.config_dir}  
            """

    rule finish_drep_with_upload:
        input:
            drep_upload=BIN_RUN + "/derep/.drep.upload" + config.get("samples_batch", "") + ".done",
            report_upload=BIN_RUN + "/.drep_report_upload" + config.get("samples_batch", "") + ".done"
        output:
            touch(BIN_RUN + "/rule_drep" + config.get("samples_batch", "") + ".done")
else:
    rule finish_drep:
        input:
            workdir=BIN_RUN + "/derep/workdir" + config.get("samples_batch", ""),
            bin_info=BIN_RUN + "/derep/bin_info" + config.get("samples_batch", "") + ".tsv",
            bins2species=BIN_RUN + "/derep/bins2species" + config.get("samples_batch", "") + ".tsv",
            report=WORKDIR + "reports/bin_report_metabat" + config.get("samples_batch", "") + ".html",
        output:
            touch(BIN_RUN + "/rule_drep" + config.get("samples_batch", "") + ".done")