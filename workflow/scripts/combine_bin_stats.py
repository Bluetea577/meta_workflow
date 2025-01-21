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
    import pandas as pd

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

    import pandas as pd

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
        logger.info(f"Inferred following list of headers {selected_headers}")

    # parse one file after another

    logger.info("Read an append table by table")
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

    common_arrguments = dict(
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
            selected_headers=selected_headers, **common_arrguments
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
            axis=axis, concat_arguments=concat_arguments, **common_arrguments
        )

try:

    pandas_concat(snakemake.input, snakemake.output[0])

except Exception as e:
    import traceback

    with open(snakemake.log[0], "w") as logfile:
        traceback.print_exc(file=logfile)

    raise e
