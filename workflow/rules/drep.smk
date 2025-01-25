rule quality_file:
    input:
        BIN_RUN + "/metabat/filtered/filtered_bins_info.tsv",
    output:
        temp(BIN_RUN + "/metabat/filtered/filtered_bins_info_for_drep.csv"),
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
        paths=BIN_RUN + "/metabat/filtered/filtered_bins_paths.txt",
        quality=BIN_RUN + "/metabat/filtered/filtered_bins_info_for_drep.csv"
    output:
        workdir=directory(BIN_RUN + "/derep/workdir"),
        genome_dir=directory(BIN_RUN + "/derep/workdir/dereplicated_genomes"),
    log:
        "logs/binning/drep/drep.log"
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
        workdir=BIN_RUN + "/derep/workdir",
        quality=BIN_RUN + "/metabat/filtered/filtered_bins_info.tsv"
    output:
        merged_info=BIN_RUN + "/derep/bin_info.tsv",
        bins2species=BIN_RUN + "/derep/bins2species.tsv",
    log:
        "logs/binning/drep/process_drep_results.log"
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

rule build_bin_report:
    input:
        bin_info=BIN_RUN + "/derep/bin_info.tsv",
        bins2species=BIN_RUN + "/derep/bins2species.tsv",
    output:
        report=WORKDIR + "reports/bin_report_metabat.html",
    conda:
        "../envs/report.yaml"
    log:
        "logs/binning/report_binning.log",
    script:
        "../scripts/bin_report.py"

localrules:
    upload_drep,
    upload_drep_report,
    finish_drep,

rule upload_drep:
    input:
        workdir=BIN_RUN + "/derep/workdir",
        bin_info=BIN_RUN + "/derep/bin_info.tsv",
        bins2species=BIN_RUN + "/derep/bins2species.tsv",
    output:
        mark=touch(BIN_RUN + "/derep/.drep.upload.done")
    params:
        remote_dir="binning/derep"
    conda:
        "../envs/baiduyun.yaml"
    log:
        "logs/binning/upload/drep.log"
    shell:
        """
        bypy mkdir {params.remote_dir} 2>> {log} 
        bypy mkdir {params.remote_dir}/workdir 2>> {log}

        bypy upload {input.workdir}/* {params.remote_dir}/workdir/ 2>> {log}
        bypy upload {input.bin_info} {params.remote_dir}/ 2>> {log}
        bypy upload {input.bins2species} {params.remote_dir}/ 2>> {log}
        """

rule upload_drep_report:
    input:
        report=WORKDIR + "reports/bin_report_metabat.html"
    output:
        mark=touch(BIN_RUN + "/.drep_report_upload.done")
    params:
        remote_dir="binning/report"
    conda:
        "../envs/baiduyun.yaml"
    log:
        "logs/binning/upload/drep_report.log"
    shell:
        """
        bypy mkdir {params.remote_dir} 2>> {log}
        bypy upload {input.report} {params.remote_dir}/bin_report_metabat.html 2>> {log}
        """

rule finish_drep:
    input:
        drep_upload=BIN_RUN + "/derep/.drep.upload.done",
        report_upload=BIN_RUN + "/.drep_report_upload.done"
    output:
        touch(BIN_RUN + "/rule_drep.done")