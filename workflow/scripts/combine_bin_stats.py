#! /usr/bin/env python


import sys
import logging, traceback

from snakemake.script import snakemake

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

import pandas as pd

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

def _pandas_concat_in_memory(
    input_tables,
    output_table,
    sep,
    index_col,
    axis,
    read_arguments,
    save_arguments,
    concat_arguments,
):

    Tables = [
        pd.read_csv(file, index_col=index_col, sep=sep, **read_arguments)
        for file in input_tables
    ]

    out = pd.concat(Tables, axis=axis, **concat_arguments).sort_index()

    del Tables

    out.to_csv(output_table, sep=sep, **save_arguments)


def _pandas_concat_disck_based(
    input_tables,
    output_table,
    sep,
    index_col,
    read_arguments,
    save_arguments,
    selected_headers=None,
):
    """combine different tables but one after the other in disk based"""


    try:
        from tqdm import tqdm
    except ImportError:
        tqdm = tuple

    if selected_headers is not None:
        try:
            selected_headers = list(selected_headers)
        except Exception as e:
            raise Exception("selected_headers should be a list-like") from e

    else:
        # read all_headers
        selected_headers = set()
        for file in input_tables:
            headers_of_file = pd.read_csv(
                file, index_col=index_col, sep=sep, nrows=2, dtype=str, **read_arguments
            )

            selected_headers.update(list(headers_of_file.columns))

        selected_headers = list(selected_headers)
        logging.info(f"Inferred following list of headers {selected_headers}")

    # parse one file after another

    logging.info("Read an append table by table")
    for file in tqdm(input_tables):
        # read full table
        table = pd.read_csv(
            file, index_col=index_col, sep=sep, dtype=str, **read_arguments
        )
        # set to common header
        table = table.reindex(selected_headers, axis=1)

        if file == input_tables[0]:
            mode = "w"
            print_header = True
        else:
            mode = "a"
            print_header = False

        table.to_csv(
            output_table, sep=sep, mode=mode, header=print_header, **save_arguments
        )


def pandas_concat(
    input_tables,
    output_table,
    sep="\t",
    index_col=0,
    axis=0,
    read_arguments=None,
    save_arguments=None,
    concat_arguments=None,
    disk_based=False,
    selected_headers=None,  # only used in disk based, not passed to usecols
):
    """
    Uses pandas to read,concatenate and save tables using pandas.concat
    """

    if read_arguments is None:
        read_arguments = {}
    if save_arguments is None:
        save_arguments = {}

    if type(input_tables) == str:
        input_tables = [input_tables]

    common_arguments = dict(
        input_tables=input_tables,
        output_table=output_table,
        sep=sep,
        index_col=index_col,
        read_arguments=read_arguments,
        save_arguments=save_arguments,
    )

    if disk_based:
        if concat_arguments is not None:
            raise Exception(
                f"cannot hanndle concat arguments by disck based append, got {concat_arguments}"
            )

        assert axis == 0, "Can only append on axis= 0"

        _pandas_concat_disck_based(
            selected_headers=selected_headers, **common_arguments
        )

    else:
        # in memory concat
        if concat_arguments is None:
            concat_arguments = {}

        if selected_headers is not None:
            raise Exception(
                "argument 'selected_headers' is not used in 'in memory' concat. Use read_arguments=dict(usecols=selected_headers) instead "
            )

        _pandas_concat_in_memory(
            axis=axis, concat_arguments=concat_arguments, **common_arguments
        )

try:
    # 获取有效的统计文件
    valid_stats_files = []
    for stats_file, status_file in zip(snakemake.input.stats, snakemake.input.status):
        with open(status_file, 'r') as f:
            status = f.read().strip()
        if status == 'valid':
            valid_stats_files.append(stats_file)
            logging.info(f"Including stats file: {stats_file}")
        else:
            logging.info(f"Skipping empty bin stats: {stats_file}")

    # 处理统计文件
    if not valid_stats_files:
        empty_stats = pd.DataFrame(columns=[
            "File", "Length", "N_scaffolds", "N50",
            "Length_contigs", "N_contigs", "Ambigious_bases"
        ])
        empty_stats.to_csv(snakemake.output[0], sep="\t", index=False)
        logging.info("No valid bins found, created empty combined stats file")
    else:
        pandas_concat(valid_stats_files, snakemake.output[0])
        logging.info(f"Successfully combined {len(valid_stats_files)} stats files")

except Exception as e:
    logging.error(f"Error: {str(e)}")
    raise