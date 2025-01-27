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

from pathlib import Path
from multiprocessing import Pool
import pandas as pd
import numpy as np
import gzip as gz


def read_fasta(filename):
    """read FASTA file with gz"""
    open_func = gz.open if str(filename).endswith('.gz') else open
    mode = 'rt' if str(filename).endswith('.gz') else 'r'

    with open_func(filename, mode) as f:
        header = ''
        sequence = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if header and sequence:
                    yield header, ''.join(sequence)
                header = line[1:]
                sequence = []
            else:
                sequence.append(line.upper())
        if header and sequence:
            yield header, ''.join(sequence)


def get_stats_from_lengths(lengths):
    sorted_lengths = sorted(lengths, reverse=True)
    csum = np.cumsum(sorted_lengths)

    total_length = int(sum(lengths))
    n = len(lengths)

    n2 = int(total_length / 2)

    # get index for cumsum >= N/2
    csumn2 = min(csum[csum >= n2])
    ind = int(np.where(csum == csumn2)[0][0])

    n50 = sorted_lengths[ind]

    return total_length, n, n50


def genome_stats(fasta_file, number_of_n_for_split=10):
    """Get genome stats from a fasta file. Outputs a tuple with:
    name,Length, n_seq,N50
    """

    try:
        name = Path(fasta_file).stem

        scaffold_lengths = []
        contig_lengths = []
        ambigious_bases = 0

        for _,sequence in read_fasta(fasta_file):
            if not set(sequence).issubset({'A', 'C', 'G', 'T', 'N'}):
                raise ValueError(f"Invalid DNA sequence in {fasta_file}")

            ambigious_bases += sequence.count('N')
            scaffold_lengths.append(len(sequence))
            contig_lengths.extend([len(contig) for contig in sequence.split('N' * number_of_n_for_split)])

        length_scaffolds, n_scaffolds, n50 = get_stats_from_lengths(scaffold_lengths)

        length_contigs, n_contigs, _ = get_stats_from_lengths(contig_lengths)

        return {
            "File": name,
            "Length_scaffolds": length_scaffolds,
            "N_scaffolds": n_scaffolds,
            "N50": n50,
            "Length_contigs": length_contigs,
            "N_contigs": n_contigs,
            "Ambigious_bases": ambigious_bases,
        }
    except Exception as e:
        print(f"Error processing {fasta_file}: {e}")
        raise

def get_many_genome_stats(filenames, output_filename, threads=1):
    """Small function to calculate total genome length and N50"""
    pool = Pool(threads)
    results = pool.map(genome_stats, filenames)
    stats = pd.DataFrame(results).rename({"Length_scaffolds": "Length"})
    stats.to_csv(output_filename, sep="\t", index=False)


with open(snakemake.input.status, 'r') as f:
    status = f.read().strip()

if status == 'empty':
    # 如果是空bins，创建空的统计文件
    empty_stats = pd.DataFrame(columns=[
        "File", "Length_scaffolds", "N_scaffolds", "N50",
        "Length_contigs", "N_contigs", "Ambigious_bases"
    ])
    empty_stats.to_csv(snakemake.output[0], sep="\t", index=False)
    logging.info("Empty bin detected, created empty stats file")
else:
    # 处理正常的bins
    filenames = list(Path(snakemake.input.bin_dir).glob("*" + snakemake.params.extension))
    get_many_genome_stats(filenames, snakemake.output[0], snakemake.threads)