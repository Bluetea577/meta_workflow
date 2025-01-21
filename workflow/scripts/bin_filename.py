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
from pathlib import Path

from snakemake.script import snakemake


def simplify_path(path, remove_gz=True):
    path = Path(path)

    name = path.stem
    ext = path.suffix

    if remove_gz & (ext == ".gz"):
        name = Path(name).stem

    return name

def get_list_of_files(dirs, pattern):
    fasta_files = []

    for dir in dirs:
        dir = Path(dir)
        fasta_files += list(dir.glob(pattern))

    filenames = pd.DataFrame(fasta_files, columns=["Filename"])
    filenames.index = filenames.Filename.apply(simplify_path)
    filenames.index.name = "Bin"

    filenames.sort_index(inplace=True)

    return filenames


fasta_filenames = get_list_of_files(snakemake.input.dirs, "*.f*")
faa_filenames = get_list_of_files(snakemake.input.protein_dirs, "*.faa")

assert all(
    faa_filenames.index == fasta_filenames.index
), "faa index and faa index are nt the same"

faa_filenames.columns = ["Proteins"]

filenames = pd.concat((fasta_filenames, faa_filenames), axis=1)

filenames.to_csv(snakemake.output.filenames, sep="\t")