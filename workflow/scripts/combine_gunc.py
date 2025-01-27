#! /usr/bin/env python

import os, sys
import logging, traceback

from snakemake.script import snakemake

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

#### Begining of scripts

import pandas as pd

def read_gunc_report(file_path):
    """读取单个GUNC报告文件,处理空bins的情况"""
    try:
        df = pd.read_csv(file_path, sep='\t')

        # 检查是否为空bins的输出
        if (df.shape[0] == 1 and
                'genome' in df.columns and
                df['genome'].iloc[0].endswith('_metabat_1')):
            return df

        return df

    except Exception as e:
        logging.error(f"Error reading file {file_path}: {str(e)}")
        # 返回空DataFrame
        return pd.DataFrame()


def combine_gunc_reports(gunc_files, status_files, samples, output_file):
    """合并所有GUNC报告"""
    df_list = []

    for sample, gunc_file, status_file in zip(samples,gunc_files, status_files):
        try:
            with open(status_file, 'r') as f:
                status = f.read().strip()

            if status == 'valid':
                df = pd.read_csv(gunc_file,sep='\t')
                df_list.append(df)
                logging.info(f"Including GUNC data for sample {sample}")
            else:
                logging.info(f"Skipping sample {sample} (status: {status})")
        except Exception as e:
            logging.error(f"Failed processing sample {sample}: {str(e)}")
            continue

    if df_list:
        # 合并所有数据
        combined_df = pd.concat(df_list, axis=0, ignore_index=True)
        # 保存结果
        combined_df.to_csv(output_file, sep='\t', index=False)
        logging.info(f"Successfully combined {len(df_list)} reports")
    else:
        # 如果没有有效数据,创建空的输出文件
        pd.DataFrame(columns=['genome', 'pass.GUNC', 'warning']) \
            .to_csv(output_file, sep='\t', index=False)
        logging.warning("No valid data to combine")


if __name__ == "__main__":
    try:
        combine_gunc_reports(
            gunc_files=snakemake.input.gunc_files,
            status_files=snakemake.input.status_files,
            samples=snakemake.params.samples,
            output_file=snakemake.output.bin_table
        )
    except Exception as e:
        logging.error(f"Script failed: {str(e)}")
        traceback.print_exc(file=sys.stderr)
        raise e