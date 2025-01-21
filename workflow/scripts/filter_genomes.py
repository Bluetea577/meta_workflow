#! /usr/bin/env python


import sys, os
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


# Install exception handler
sys.excepthook = handle_exception


import pandas as pd
from glob import glob
from numpy import log
from warnings import warn

def load_quality(quality_file):
    Q = pd.read_csv(quality_file, index_col=0, sep="\t")

    # remove extension if present
    if Q.index.str.contains(".fa").all():
        warn("Found fasta extension in index. I remove them")
        Q.index = Q.index.str.split(".fa", expand=True).to_frame()[0]

    # Q.columns = Q.columns.str.lower()

    necessary_columns = ["Completeness", "Contamination"]

    # rename lower and uppercase to necessary_columns
    Q = Q.rename(
        columns={
            fun(s[0]) + s[1:]: s
            for s in necessary_columns
            for fun in (str.lower, str.upper)
        }
    )

    if Q.columns.isin(necessary_columns).sum() != len(necessary_columns):
        raise Exception(
            f"{necessary_columns} should be in the quality table, only got {Q.columns}"
        )

    assert not Q.index.duplicated().any(), f"duplicated indexes in {quality_file}"

    return Q

def process_quality_and_stats(quality_file, stats_file):
    Q = load_quality(quality_file)

    stats = pd.read_csv(stats_file, index_col=0, sep="\t")
    stats["logN50"] = log(stats.N50)

    Q = Q.join(stats.loc[Q.index, stats.columns.difference(Q.columns)])

    return Q

def apply_gunc_filter(Q, gunc_file):
    gunc = pd.read_table(gunc_file, index_col=0)
    gunc = gunc.loc[Q.index]

    bad_genomes = gunc.index[gunc["pass.GUNC"] == False]
    logging.info(f"{len(bad_genomes)} Don't pass gunc filtering")

    Q.drop(bad_genomes, inplace=True)
    return Q

# Load and process quality and stats data
Q = process_quality_and_stats(snakemake.input.quality, snakemake.input.stats)
n_all_bins = Q.shape[0]

# Apply filtering criteria
filter_criteria = snakemake.params["filter_criteria"]
logging.info(f"Filter genomes according to criteria:\n {filter_criteria}")

Q = Q.query(filter_criteria)
logging.info(f"Retain {Q.shape[0]} genomes from {n_all_bins}")

# Apply GUNC filtering if available
if hasattr(snakemake.input, "gunc"):
    Q = apply_gunc_filter(Q, snakemake.input.gunc)
else:
    logging.info(" Don't filter based on gunc")

# Check if any bins remain
if Q.shape[0] == 0:
    logging.error(
        f"No bins passed filtering criteria! Bad luck!. You might want to tweek the filtering criteria. Also check the {snakemake.input.quality}"
    )
    exit(1)

# Save filtered quality data
Q.to_csv(snakemake.output.info, sep="\t")

# Filter and save path information
F = pd.read_table(snakemake.input.paths, index_col=0).squeeze()
F = F.loc[Q.index].iloc[:, 0]
F.to_csv(snakemake.output.paths, index=False, header=False)