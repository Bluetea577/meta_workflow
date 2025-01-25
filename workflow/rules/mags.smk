localrules:
    calculate_contig_coverage,

# 计算contig覆盖度和基础统计量
rule calculate_contig_coverage:
    """计算单个样本中每个contig的覆盖度、reads数和TPM"""
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
        # 获取总reads数和mapped reads数  
        total_reads=$(samtools view -c {input.bam})  
        mapped_reads=$(samtools view -F 0x4 -c {input.bam})  

        # 先获取每个contig的实际reads数  
        samtools idxstats {input.bam} | cut -f 1,3 > contig_reads.tmp  

        # 然后计算深度  
        samtools depth -a {input.bam} | \
        awk -v total="$total_reads" -v mapped="$mapped_reads" '  
        BEGIN {{  
            # 读取每个contig的实际reads数  
            while(getline < "contig_reads.tmp") {{  
                contig_reads[$1] = $2  
            }}  
            print "contig\tlength\tcoverage\treads\tTPM\trel_abundance"  
        }}  
        {{  
            depth[$1] += $3  
            bases[$1]++  
        }}  
        END {{  
            # 计算RPK  
            rpk_sum = 0  
            total_length = 0  
            total_bases_covered = 0  
            max_coverage = 0  
            min_coverage = 999999  

            for(i in depth) {{  
                coverage = depth[i]/bases[i]  
                reads = contig_reads[i]  # 使用实际的reads数  
                rpk = (reads * 1000.0)/bases[i]  
                rpk_values[i] = rpk  
                rpk_sum += rpk  

                # 保存数据  
                reads_count[i] = reads  
                coverage_values[i] = coverage  
                length_values[i] = bases[i]  

                # 统计信息  
                total_length += bases[i]  
                if(coverage > max_coverage) max_coverage = coverage  
                if(coverage < min_coverage) min_coverage = coverage  
            }}  

            # 计算TPM  
            scaling_factor = rpk_sum/1000000  

            # 输出结果  
            for(i in depth) {{  
                tpm = rpk_values[i]/scaling_factor  
                rel_abundance = (reads_count[i]/total) * 100  

                print i "\t" length_values[i] "\t" coverage_values[i] "\t" reads_count[i] "\t" tpm "\t" rel_abundance  
            }}  

            # 输出统计信息到log文件  
            printf "=== Alignment Statistics ===\n" > "/dev/stderr"  
            printf "Total reads: %d\n", total > "/dev/stderr"  
            printf "Mapped reads: %d (%.2f%%)\n", mapped, (mapped/total)*100 > "/dev/stderr"  
            printf "Unmapped reads: %d (%.2f%%)\n", total-mapped, ((total-mapped)/total)*100 > "/dev/stderr"  
            printf "\n=== Contig Statistics ===\n" > "/dev/stderr"  
            printf "Total contigs: %d\n", length(depth) > "/dev/stderr"  
            printf "Total length: %d bp\n", total_length > "/dev/stderr"  
            printf "Maximum coverage: %.2f\n", max_coverage > "/dev/stderr"  
            printf "Minimum coverage: %.2f\n", min_coverage > "/dev/stderr"  
            printf "Average coverage: %.2f\n", (rpk_sum/length(depth)) > "/dev/stderr"  
        }}' > {output.coverage} 2> {log}  

        rm contig_reads.tmp
        """

rule calculate_bin_abundance:
    """将contig统计量汇总到bin级别"""
    input:
        coverage=BIN_RUN + "/contigs_coverage/{sra_run}_contig_coverage.tsv",
        contig2bin=BIN_RUN + "/metabat/{sra_run}/cluster_attribution.tsv",
    output:
        abundance=BIN_RUN + "/bin_abundances/{sra_run}_bin_abundance.tsv",
    threads: 2
    resources:
        mem_mb=3500 * 2
    run:
        import pandas as pd

        # 读取contig统计量
        coverage_df = pd.read_csv(input.coverage, sep='\t')

        # 读取contig到bin的映射
        contig2bin = {}
        with open(input.contig2bin) as f:
            for line in f:
                contig, bin_name = line.strip().split('\t')
                contig2bin[contig] = bin_name

        # 初始化bin统计字典
        bin_stats = {}

        # 计算bin级别的统计量
        for _, row in coverage_df.iterrows():
            contig = row['contig']
            if contig in contig2bin:
                bin_name = contig2bin[contig]
                length = row['length']
                coverage = row['coverage']
                reads = row['reads']
                tpm = row['TPM']
                rel_abundance = row['rel_abundance']

                if bin_name not in bin_stats:
                    bin_stats[bin_name] = {
                        'total_length': 0,
                        'weighted_coverage': 0,
                        'total_reads': 0,
                        'weighted_tpm': 0,
                        'total_rel_abundance': 0
                    }

                bin_stats[bin_name]['total_length'] += length
                bin_stats[bin_name]['weighted_coverage'] += length * coverage
                bin_stats[bin_name]['total_reads'] += reads
                bin_stats[bin_name]['weighted_tpm'] += length * tpm
                bin_stats[bin_name]['total_rel_abundance'] += rel_abundance  # 直接累加相对丰度

        # 计算bin的最终统计量
        bin_abundance = []
        for bin_name, stats in bin_stats.items():
            bin_abundance.append({
                'bin': bin_name,
                'length': stats['total_length'],
                'coverage': stats['weighted_coverage'] / stats['total_length'],
                'reads': stats['total_reads'],
                'TPM': stats['weighted_tpm'] / stats['total_length'],
                'rel_abundance': stats['total_rel_abundance']  # bin的相对丰度是其包含的contigs相对丰度之和
            })

        # 保存结果
        abundance_df = pd.DataFrame(bin_abundance)
        # 按相对丰度降序排序
        abundance_df = abundance_df.sort_values('rel_abundance', ascending=False)
        abundance_df.to_csv(output.abundance, sep='\t', index=False)

rule combine_bin_abundances:
    """合并所有样本的bin统计量为矩阵"""
    input:
        abundances=expand(BIN_RUN + "/bin_abundances/{sra_run}_bin_abundance.tsv",
            sra_run=IDS),
    output:
        matrix=BIN_RUN + "/bin_abundances/bin_abundance_matrix.tsv",
        summary=BIN_RUN + "/bin_abundances/bin_abundance_summary.tsv",# 添加一个汇总文件
    threads: 2
    resources:
        mem_mb=3500 * 2
    run:
        import pandas as pd
        from pathlib import Path

        # 读取所有样本的丰度文件
        abundance_dfs = []  # 用于存储相对丰度
        summary_dfs = []  # 用于存储完整统计信息

        for file in input.abundances:
            sample = Path(file).stem.replace('_bin_abundance','')
            df = pd.read_csv(file,sep='\t')

            # 为汇总文件准备数据
            summary_df = df.copy()
            summary_df['sample'] = sample
            summary_dfs.append(summary_df)

            # 为丰度矩阵准备数据
            abundance_df = df[['bin', 'rel_abundance']].copy()
            abundance_df = abundance_df.rename(columns={'rel_abundance': sample})
            abundance_df = abundance_df.set_index('bin')
            abundance_dfs.append(abundance_df)

        # 创建丰度矩阵
        abundance_matrix = pd.concat(abundance_dfs,axis=1)
        abundance_matrix = abundance_matrix.fillna(0)  # 填充缺失值为0

        # 创建汇总表
        summary_matrix = pd.concat(summary_dfs,ignore_index=True)
        summary_matrix = summary_matrix.sort_values(['sample', 'rel_abundance'],
            ascending=[True, False])

        # 保存结果
        abundance_matrix.to_csv(output.matrix,sep='\t')
        summary_matrix.to_csv(output.summary,sep='\t',index=False)

rule map_to_representative_mags:
    """将bin统计量映射到代表性MAGs"""
    input:
        summary=BIN_RUN + "/bin_abundances/bin_abundance_summary.tsv",
        bins2species=BIN_RUN + "/derep/bins2species.tsv",
    output:
        coverage=BIN_RUN + "/bin_abundances/species_coverage_matrix.tsv",
        abundance=BIN_RUN + "/bin_abundances/species_relative_abundance_matrix.tsv",
        tpm=BIN_RUN + "/bin_abundances/species_tpm_matrix.tsv"
    threads: 2
    resources:
        mem_mb=3500 * 2
    run:
        import pandas as pd

        # 读取数据文件
        summary_df = pd.read_csv(input.summary, sep='\t')
        species_df = pd.read_csv(input.bins2species, sep='\t')

        # 获取样本列表
        samples = summary_df['sample'].unique()

        # 创建物种映射字典和代表性bins列表
        species_mapping = dict(zip(species_df['genome'], species_df['Species']))
        representative_bins = species_df[species_df['Representative'].astype(str).str.lower() == 'true']['genome'].tolist()

        # 初始化结果字典
        species_coverage = {sample: {} for sample in samples}
        species_abundance = {sample: {} for sample in samples}
        species_tpm = {sample: {} for sample in samples}

        # 对每个样本进行处理
        for sample in samples:
            # 提取该样本的数据
            sample_data = summary_df[summary_df['sample'] == sample]

            # 初始化所有代表性bins的值为0
            for rep_bin in representative_bins:
                species_coverage[sample][rep_bin] = 0
                species_abundance[sample][rep_bin] = 0
                species_tpm[sample][rep_bin] = 0

            # 计算总reads用于相对丰度
            total_reads = sample_data['reads'].sum()

            # 映射到代表性MAGs
            for _, row in sample_data.iterrows():
                bin_name = row['bin']
                if bin_name in species_mapping:
                    species = species_mapping[bin_name]
                    rep_bin = next(rb for rb in representative_bins
                                 if species_mapping[rb] == species)

                    # 更新覆盖度（取最大值）
                    species_coverage[sample][rep_bin] = max(
                        species_coverage[sample][rep_bin],
                        row['coverage']
                    )

                    # 更新相对丰度（累加reads后再计算）
                    if total_reads > 0:
                        relative_abundance = (row['reads'] / total_reads) * 100
                    else:
                        relative_abundance = 0
                    species_abundance[sample][rep_bin] = max(
                        species_abundance[sample][rep_bin],
                        relative_abundance
                    )

                    # 更新TPM（取最大值）
                    species_tpm[sample][rep_bin] = max(
                        species_tpm[sample][rep_bin],
                        row['TPM']
                    )

        # 转换为DataFrame
        coverage_matrix = pd.DataFrame(species_coverage)
        abundance_matrix = pd.DataFrame(species_abundance)
        tpm_matrix = pd.DataFrame(species_tpm)

        # 确保所有缺失值都是0
        coverage_matrix = coverage_matrix.fillna(0)
        abundance_matrix = abundance_matrix.fillna(0)
        tpm_matrix = tpm_matrix.fillna(0)

        # 保存结果
        coverage_matrix.to_csv(output.coverage, sep='\t')
        abundance_matrix.to_csv(output.abundance, sep='\t')
        tpm_matrix.to_csv(output.tpm, sep='\t')

localrules:
    upload_abundances,
    finish_mags,

rule upload_abundances:
    """上传所有统计结果"""
    input:
        contig_coverage=expand(BIN_RUN + "/contigs_coverage/{sra_run}_contig_coverage.tsv",
            sra_run=IDS),
        bin_abundance=expand(BIN_RUN + "/bin_abundances/{sra_run}_bin_abundance.tsv",
            sra_run=IDS),
        bin_matrix=BIN_RUN + "/bin_abundances/bin_abundance_matrix.tsv",  # 改名以避免重复
        species_coverage=BIN_RUN + "/bin_abundances/species_coverage_matrix.tsv",
        species_abundance=BIN_RUN + "/bin_abundances/species_relative_abundance_matrix.tsv",
        species_tpm=BIN_RUN + "/bin_abundances/species_tpm_matrix.tsv"
    output:
        mark=touch(BIN_RUN + "/.abundance_upload.done")
    params:
        remote_dir="binning/abundance"
    conda:
        config["upload"]
    log:
        "logs/binning/upload/abundance.log"
    shell:
        """  
        # 创建远程目录结构  
        bypy mkdir {params.remote_dir} 2>> {log}  
        bypy mkdir {params.remote_dir}/contigs_coverage 2>> {log}  
        bypy mkdir {params.remote_dir}/bin_abundances 2>> {log}  

        # 上传contig覆盖度文件  
        for f in {input.contig_coverage}; do  
            bypy upload "$f" {params.remote_dir}/contigs_coverage/ 2>> {log}  
        done  

        # 上传bin丰度文件  
        for f in {input.bin_abundance}; do  
            bypy upload "$f" {params.remote_dir}/bin_abundances/ 2>> {log}  
        done  

        # 上传矩阵文件  
        bypy upload {input.bin_matrix} {params.remote_dir}/ 2>> {log}  
        bypy upload {input.species_coverage} {params.remote_dir}/ 2>> {log}  
        bypy upload {input.species_abundance} {params.remote_dir}/ 2>> {log}  
        bypy upload {input.species_tpm} {params.remote_dir}/ 2>> {log}  
        """

rule finish_mags:
    input:
        abundance_upload=BIN_RUN + "/.abundance_upload.done"
    output:
        touch(BIN_RUN + "/rule_mags.done")