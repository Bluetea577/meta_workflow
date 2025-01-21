#! /usr/bin/env python

import os, sys
import logging, traceback

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

logging.captureWarnings(True)


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


# Install exception handler
sys.excepthook = handle_exception

#### Begining
import pandas as pd
import sys

# 读取注释结果
df = pd.read_csv(snakemake.input.combined_tsv, sep='\t')

# 计算统计信息
stats = {
    'Total genomes': df['Genome'].nunique(),
    'Total genes': len(df),
    'Genes with EC numbers': df['EC_number'].notna().sum(),
    'Unique EC numbers': df['EC_number'].dropna().nunique(),
    'Genes with function': df['Product'].notna().sum(),
    'Average genes per genome': len(df) / df['Genome'].nunique()
}

# 写入统计结果
with open(snakemake.output.stats, 'w') as f:
    for key, value in stats.items():
        f.write(f'{key}: {value}\n')