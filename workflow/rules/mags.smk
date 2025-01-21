localrules:
    calculate_contig_coverage,


rule calculate_contig_coverage:
    """计算单个样本中每个contig的覆盖度"""
    input:
        bam=ASSE_RUN + "/megahit/{sra_run}/sequence_alignment/{sra_run}_sort.bam",
        bai=ASSE_RUN + "/megahit/{sra_run}/sequence_alignment/{sra_run}_sort.bam.bai",
    output:
        coverage=BIN_RUN + "/contigs_coverage/{sra_run}_contig_coverage.tsv",
    conda:
        "../envs/align.yaml"
    benchmark:
        "logs/benchmarks/binning/mags/{sra_run}_calculate_contig_coverage.tsv"
    log:
        "logs/binning/mags/{sra_run}_calculate_contig_coverage.log"
    threads: 1
    resources:
        mem_mb=3000
    shell:
        """
        samtools depth -a {input.bam} | \
        awk 'BEGIN {{print "contig\tlength\tcoverage"}}
            {{sum[$1]+=$3; len[$1]++}}
            END {{for(i in sum) print i "\t" len[i] "\t" sum[i]/len[i]}}' \
            > {output.coverage} 2> {log}
        """

rule calculate_bin_abundance:
    """将contig覆盖度汇总到bin级别"""
    input:
        coverage=BIN_RUN + "/contigs_coverage/{sra_run}_contig_coverage.tsv",
        contig2bin=BIN_RUN + "/metabat/{sra_run}/cluster_attribution.tsv",
    output:
        abundance=BIN_RUN + "/bin_abundances/{sra_run}_bin_abundance.tsv",
    threads: 2
    resources:
        mem_mb=3500*2
    run:
        import pandas as pd

        # 读取contig覆盖度
        coverage_df = pd.read_csv(input.coverage, sep='\t')

        # 读取contig到bin的映射
        contig2bin = {}
        with open(input.contig2bin) as f:
            for line in f:
                contig, bin_name = line.strip().split('\t')
                contig2bin[contig] = bin_name

        # 初始化bin统计字典
        bin_stats = {}

        # 计算bin级别的统计
        for _, row in coverage_df.iterrows():
            contig = row['contig']
            if contig in contig2bin:
                bin_name = contig2bin[contig]
                length = row['length']
                coverage = row['coverage']

                if bin_name not in bin_stats:
                    bin_stats[bin_name] = {
                        'total_length': 0,
                        'weighted_coverage': 0
                    }

                bin_stats[bin_name]['total_length'] += length
                bin_stats[bin_name]['weighted_coverage'] += length * coverage

        # 计算每个bin的平均覆盖度
        bin_abundance = {
            bin_name: stats['weighted_coverage'] / stats['total_length']
            for bin_name, stats in bin_stats.items()
        }

        # 保存结果
        abundance_df = pd.DataFrame({
            'bin': list(bin_abundance.keys()),
            'abundance': list(bin_abundance.values())
        })
        abundance_df.to_csv(output.abundance, sep='\t', index=False)

rule combine_bin_abundances:
    """合并所有样本的bin丰度为一个矩阵"""
    input:
        abundances=expand(BIN_RUN + "/bin_abundances/{sra_run}_bin_abundance.tsv",
            sra_run=IDS),
    output:
        matrix=BIN_RUN + "/bin_abundances/bin_abundance_matrix.tsv",
    threads: 2
    resources:
        mem_mb=3500*2
    run:
        import pandas as pd
        from pathlib import Path

        # 读取所有样本的丰度文件
        dfs = []
        for file in input.abundances:
            sample = Path(file).stem.replace('_bin_abundance', '')
            df = pd.read_csv(file, sep='\t')
            df = df.set_index('bin')
            df.columns = [sample]  # 重命名列为样本名
            dfs.append(df)

        # 合并所有数据框
        abundance_matrix = pd.concat(dfs, axis=1)
        abundance_matrix = abundance_matrix.fillna(0)  # 填充缺失值为0

        # 保存结果
        abundance_matrix.to_csv(output.matrix, sep='\t')

rule map_to_representative_mags:
    """将bin丰度映射到代表性MAGs"""
    input:
        abundance_matrix=BIN_RUN + "/bin_abundances/bin_abundance_matrix.tsv",
        bins2species=BIN_RUN + "/derep/bins2species.tsv",
    output:
        species_abundance=BIN_RUN + "/bin_abundances/species_abundance_matrix.tsv",
    threads: 2
    resources:
        mem_mb=3500*2
    run:
        import pandas as pd

        abundance_df = pd.read_csv(input.abundance_matrix,sep='\t',index_col=0)
        species_df = pd.read_csv(input.bins2species,sep='\t')

        species_mapping = dict(zip(species_df['genome'],species_df['Species']))
        representative_bins = species_df[species_df['Representative']]['genome'].tolist()
        species_abundance = {}

        for rep_bin in representative_bins:
            species_abundance[rep_bin] = {sample: 0 for sample in abundance_df.columns}

        # 计算每个物种的最高丰度
        for bin_name in abundance_df.index:
            if bin_name in species_mapping:
                species = species_mapping[bin_name]
                # 找到对应的代表性bin
                rep_bin = next(rb for rb in representative_bins
                               if species_mapping[rb] == species)

                # 更新每个样本的丰度
                for sample in abundance_df.columns:
                    current_abundance = abundance_df.loc[bin_name, sample]
                    species_abundance[rep_bin][sample] = max(
                        species_abundance[rep_bin][sample],
                        current_abundance
                    )

        # 转换为DataFrame并保存
        species_abundance_df = pd.DataFrame(species_abundance)
        species_abundance_df.to_csv(output.species_abundance,sep='\t')