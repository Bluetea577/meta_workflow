GTDB_DIR = os.path.join(WORKDIR,config["path"]["gtdbtk_path"])

"""
rule gtdbtk_db_download:
    output:
        gtdbtkdb=DBDIR + "/gtdb_r220/gtdbtk_data.tar.gz",
    params: 
        outdir=lambda wc,output: os.path.dirname(output.gtdbtkdb),
    shell:
        " cd {params.outdir} "
        " ; "
        " wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz "
"""

localrules:
    gtdbtk_db_preparation

rule gtdbtk_db_preparation:
    input:
        database=DBDIR + "/gtdb_r220/gtdbtk_data.tar.gz"
    output:
        marker=touch(DBDIR + "/gtdb_r220/.db_extract_complete")
    threads: 1
    resources:
        time_min=60 * int(config.get("runtime",{"long": 10})["long"])
    log:
        stdout="logs/download/gtdbtk_untar.log",
        stderr="logs/download/gtdbtk_untar.err"
    params:
        extract_path=DBDIR + "/gtdb_r220"
    shell:
        "tar -xzvf {input.database} -C {params.extract_path} --strip 1 2> {log.stderr} > {log.stdout}"

rule gtdbtk_classify_wf:
    input:
        genome_dir=BIN_RUN + "/derep/workdir/dereplicated_genomes",
        db_ready=DBDIR + "/gtdb_r220/.db_extract_complete",
    output:
        directory=directory(GTDB_DIR),
        marker=touch(GTDB_DIR + "/.gtdbtk_complete"),
    params:
        extension=".fasta",
        gtdb_database_path=DBDIR + "/gtdb_r220",
        mash_mode=config.get("mash_mode", "--skip_ani_screen"),
    threads: 32
    resources:
        mem_mb=3500 * 32,
        time_min=60 * 24
    conda:
        "../envs/gtdbtk.yaml"
    log:
        "logs/gtdbtk/gtdbtk_classify.log"
    benchmark:
        "logs/benchmarks/gtdbtk/gtdbtk_classify.tsv"
    shell:
        """  
        mkdir -p {output.directory}  

        export GTDBTK_DATA_PATH="{params.gtdb_database_path}"  

        gtdbtk classify_wf \
            --genome_dir {input.genome_dir} \
            --out_dir {output.directory} \
            --prefix gtdbtk_result \
            --cpus {threads} \
            --extension {params.extension} \
            {params.mash_mode} \
            &> {log}  
        """

"""
rule gtdbtk_infer:
    input:
        msa=GTDB_DIR + "/gtdbtk_result.bac120.user_msa.fasta",
    output:
        tree=directory(GTDB_DIR + "/phylogeny/infer"),
    threads: 16
    resources:
        mem_mb=3500 * 8,
        time_min=60 * 6
    log:
        "logs/phylogeny/gtdbtk_infer.log"
    conda:
        "../envs/gtdbtk.yaml"
    shell:
        " gtdbtk infer --msa_file {input.msa} --out_dir {output} --cpu {threads} "
"""

rule fasttree_phylogeny:
    input:
        flags=GTDB_DIR + "/.gtdbtk_complete",
    output:
        tree=GTDB_DIR + "/phylogeny/fasttree.nwk",
    params:
        input_msa=GTDB_DIR + "/align/gtdbtk_result.bac120.user_msa.fasta.gz",
        outdir=lambda wc, output: os.path.dirname(output.tree),
    threads: 16
    resources:
        mem_mb=3500 * 8,
        time_min=60 * 6
    log:
        "logs/phylogeny/fasttree.log"
    conda:
        "../envs/fasttree.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        
        zcat {params.input_msa} | \
        FastTreeMP \
            -nt \
            -gamma \
            -gtr \
            > {output.tree} \
            2> {log}
        """

rule upload_gtdbtk:
    """上传GTDB-Tk分析结果"""
    input:
        gtdbtk_dir=GTDB_DIR,
        classify_done=GTDB_DIR + "/.gtdbtk_complete",
        tree=GTDB_DIR + "/phylogeny/fasttree.nwk"
    output:
        mark=touch(GTDB_DIR + "/.upload_complete")
    params:
        remote_dir="taxonomy/gtdbtk"
    conda:
        "../envs/baiduyun.yaml"
    retries: 3
    log:
        "logs/gtdbtk/upload.log"
    shell:
        """  
        # 创建远程目录  
        bypy mkdir {params.remote_dir} 2>> {log}  

        # 上传GTDB-Tk结果  
        bypy upload {input.gtdbtk_dir}/ {params.remote_dir}/ 2>> {log}  
        """

rule finish_gtdbtk:
    """标记GTDB-Tk分析完成"""
    input:
        classify=GTDB_DIR + "/.gtdbtk_complete",
        tree=GTDB_DIR + "/phylogeny/fasttree.nwk",
        upload=GTDB_DIR + "/.upload_complete"
    output:
        touch(GTDB_DIR + "/rule_gtdbtk.done")