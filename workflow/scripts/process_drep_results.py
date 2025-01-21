#! /usr/bin/env python


import sys
import logging, traceback

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logging.error(
        "".join(
            [
                "Uncaught exception: ",
                *traceback.format_exception(exc_type, exc_value, exc_traceback),
            ]
        )
    )


sys.excepthook = handle_exception


import pandas as pd
import os

def process_drep_results(workdir, quality_file):
    """
    处理dRep的去冗余结果

    Args:
        workdir: dRep工作目录路径
        quality_file: 原始质量信息文件路径

    Returns:
        merged_info: 合并后的质量信息
        bins2species: bin到物种的映射关系
    """
    # 读取文件
    cdb = pd.read_csv(os.path.join(workdir, "data_tables", "Cdb.csv"))
    wdb = pd.read_csv(os.path.join(workdir, "data_tables", "Wdb.csv"))
    quality = pd.read_csv(quality_file, sep='\t')

    # 找出代表性基因组（每个cluster中score最高的）
    representatives = wdb.loc[wdb.groupby('cluster')['score'].idxmax(), 'genome'].tolist()

    # 处理聚类结果
    bins2species = pd.DataFrame({
        "Species": cdb["secondary_cluster"],
        "Representative": cdb["genome"].isin(representatives)
    })

    # 设置索引为genome（不带.fasta后缀）
    bins2species.index = cdb["genome"].str.replace('.fasta', '')

    quality.set_index('Bin Id', inplace=True)

    return quality, bins2species


if __name__ == "__main__":
    try:
        merged_info, bins2species = process_drep_results(
            snakemake.input.workdir,
            snakemake.input.quality
        )

        # 保存结果
        merged_info.to_csv(snakemake.output.merged_info, sep='\t')
        bins2species.to_csv(snakemake.output.bins2species, sep='\t')
    except Exception as e:
        print(f"Error processing dRep results: {str(e)}")
        import traceback

        traceback.print_exc()
        raise