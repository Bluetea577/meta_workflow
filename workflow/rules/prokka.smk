PROKKA_DIR = os.path.join(WORKDIR,config["path"]["prokka_path"])

# Prokka注释单个基因组
rule prokka_annotation:
    input:
        genome=BIN_RUN + "/derep/workdir/dereplicated_genomes/{genome}.fasta"
    output:
        gff=PROKKA_DIR + "/{genome}/annotation.gff",
        faa=PROKKA_DIR + "/{genome}/annotation.faa",
        ffn=PROKKA_DIR + "/{genome}/annotation.ffn",
        gbk=PROKKA_DIR + "/{genome}/annotation.gbk",
        tsv=PROKKA_DIR + "/{genome}/annotation.tsv"
    params:
        outdir=PROKKA_DIR + "/{genome}",
        prefix="annotation",
        kingdom="Bacteria",# or "Archaea"
        locustag="{genome}",
        mincontiglen=200,# 最小contig长度
        compliant=True  # 生成NCBI兼容的输出
    threads: 16
    resources:
        mem_mb=3500 * 16,
        time_min=30
    log:
        "logs/prokka/{genome}.log"
    conda:
        "../envs/prokka.yaml"
    shell:
        """  
        mkdir -p {params.outdir}  

        prokka \
            --outdir {params.outdir} \
            --prefix {params.prefix} \
            --kingdom {params.kingdom} \
            --locustag {params.locustag} \
            --mincontiglen {params.mincontiglen} \
            --cpus {threads} \
            --force \
            --compliant \
            {input.genome} \
            &> {log}  
        """

localrules:
    get_derep_genomes,

# 获取derep后的基因组列表
checkpoint get_derep_genomes:
    input:
        derep_dir=BIN_RUN + "/derep/workdir/dereplicated_genomes"
    output:
        genome_list=PROKKA_DIR + "/genome_list.txt"
    shell:
        """  
        mkdir -p $(dirname {output.genome_list})  

        if [ ! -d "{input.derep_dir}" ]; then  
            echo "Error: Dereplicated genomes directory not found" >&2  
            exit 1  
        fi  

        ls {input.derep_dir}/*.fasta | xargs -n 1 basename | sed 's/.fasta//' > {output.genome_list}  
        """


# 获取Prokka输入文件列表的函数
def get_prokka_input(wildcards):
    checkpoint_output = checkpoints.get_derep_genomes.get(**wildcards).output[0]
    with open(checkpoint_output) as f:
        genomes = [line.strip() for line in f]
    return expand(PROKKA_DIR + "/{genome}/annotation.{ext}",
        genome=genomes,
        ext=["faa", "ffn", "tsv"])


# 合并所有Prokka结果
rule combine_prokka_results:
    input:
        annotations=get_prokka_input
    output:
        combined_faa=PROKKA_DIR + "/all_proteins.faa",
        combined_ffn=PROKKA_DIR + "/all_genes.ffn",
        combined_tsv=PROKKA_DIR + "/all_annotations.tsv"
    params:
        prokka_dir=PROKKA_DIR
    log:
        "logs/prokka/combine_results.log"
    shell:
        """  
        mkdir -p $(dirname {output.combined_faa})  

        # 检查输入文件  
        for f in {params.prokka_dir}/*/annotation.{{faa,ffn,tsv}}; do  
            if [ ! -f "$f" ]; then  
                echo "Error: Missing file $f" >&2  
                exit 1  
            fi  
        done  

        # 合并蛋白质序列  
        cat {params.prokka_dir}/*/annotation.faa > {output.combined_faa} 2>> {log}  

        # 合并核酸序列  
        cat {params.prokka_dir}/*/annotation.ffn > {output.combined_ffn} 2>> {log}  

        # 合并注释信息  
        echo -e "Genome\\tLocus_tag\\tFtype\\tLength_bp\\tGene\\tEC_number\\tProduct" > {output.combined_tsv}  
        for f in {params.prokka_dir}/*/annotation.tsv; do  
            genome=$(basename $(dirname $f))  
            tail -n +2 $f | sed "s/^/$genome\\t/" >> {output.combined_tsv}  
        done 2>> {log}  
        """

# 统计注释结果
rule prokka_stats:
    input:
        combined_tsv=PROKKA_DIR + "/all_annotations.tsv"
    output:
        stats=PROKKA_DIR + "/annotation_statistics.txt"
    threads: 2
    resources:
        mem_mb=3500 * 2,
        time_min=30
    log:
        "logs/prokka/stats.log"
    script:
        "../scripts/summarize_prokka.py"

localrules:
    upload_prokka,
    finish_prokka,

rule upload_prokka:
    """上传Prokka注释结果"""
    input:
        combined_faa=PROKKA_DIR + "/all_proteins.faa",
        combined_ffn=PROKKA_DIR + "/all_genes.ffn",
        combined_tsv=PROKKA_DIR + "/all_annotations.tsv",
        stats=PROKKA_DIR + "/annotation_statistics.txt"
    output:
        mark=touch(PROKKA_DIR + "/.upload_complete")
    params:
        remote_dir="annotation/prokka",
        prokka_dir=PROKKA_DIR
    conda:
        config["upload"]
    retries: 3
    log:
        "logs/prokka/upload.log"
    shell:
        """  
        # 创建远程目录  
        bypy mkdir {params.remote_dir} 2>> {log}  

        # 上传合并后的结果文件  
        bypy upload {input.combined_faa} {params.remote_dir}/ 2>> {log}  
        bypy upload {input.combined_ffn} {params.remote_dir}/ 2>> {log}  
        bypy upload {input.combined_tsv} {params.remote_dir}/ 2>> {log}  
        bypy upload {input.stats} {params.remote_dir}/ 2>> {log}  

        # 上传单个基因组的注释结果  
        for genome_dir in {params.prokka_dir}/*/; do  
            if [ -d "$genome_dir" ]; then  
                genome=$(basename $genome_dir)  
                bypy mkdir {params.remote_dir}/$genome 2>> {log}  
                bypy upload $genome_dir/ {params.remote_dir}/$genome/ 2>> {log}  
            fi  
        done  
        """

rule finish_prokka:
    """标记Prokka注释完成"""
    input:
        combined_faa=PROKKA_DIR + "/all_proteins.faa",
        combined_ffn=PROKKA_DIR + "/all_genes.ffn",
        combined_tsv=PROKKA_DIR + "/all_annotations.tsv",
        stats=PROKKA_DIR + "/annotation_statistics.txt",
        upload=PROKKA_DIR + "/.upload_complete"
    output:
        touch(PROKKA_DIR + "/rule_prokka.done")