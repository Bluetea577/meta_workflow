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

#### Begining of scripts

import pandas as pd

def read_checkm2_output(
    completness_table, quality_score_formula="Completeness - 5*Contamination"
):
    df = pd.read_table(completness_table, index_col=0)

    if not "Completeness" in df.columns:
        # create empty column
        df.insert(0, "Completeness", 0.0)

        # add completeness depending on selected model
        specific = df.Completeness_Model_Used.str.contains("Specific Model")
        df.loc[specific, "Completeness"] = df.loc[specific, "Completeness_Specific"]
        df.loc[~specific, "Completeness"] = df.loc[~specific, "Completeness_General"]

    df.eval(
        "Quality_score = " + quality_score_formula,
        inplace=True,
    )

    df.index.name = "Bin Id"

    return df


def main(samples, completeness_files, status_files, bin_table):
    sample_data = {}
    div = {}

    df_list = []

    for i, (sample, completeness_file, status_file) in enumerate(zip(samples, completeness_files, status_files)):
        try:
            with open(status_file, "r") as f:
                status = f.read().strip()

            if status == 'valid':
                sample_data = read_checkm2_output(completness_table=completeness_files[i])
                sample_data["Sample"] = sample

                df_list.append(sample_data)
                logging.info(f"Including checkm2 data for sample {sample}")
            else:
                logging.info(f"Skipping sample {sample} (status: {status})")

        except Exception as e:
            logging.error(f"Error processing sample {sample}: {str(e)}")

    if df_list:
        df = pd.concat(df_list, axis=0)
        df.to_csv(bin_table, sep="\t")
        logging.info(f"Successfully combined {len(df_list)} valid samples")
    else:
        logging.warning("No valid data to combine")
        # 创建空的输出文件
        pd.DataFrame(columns=["Bin Id", "Completeness", "Contamination", "Quality_score", "Sample"]) \
            .to_csv(bin_table, sep="\t", index=False)


if __name__ == "__main__":
    main(
        samples=snakemake.params.samples,
        completeness_files=snakemake.input.completeness_files,
        status_files=snakemake.input.status_files,
        bin_table=snakemake.output.bin_table,
    )
