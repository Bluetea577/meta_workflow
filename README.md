# Metagenomics Data Processing Workflow

This Snakemake workflow is designed for processing metagenomic sequencing data, handling everything from data download to cleaning and analysis.

## Workflow Structure

<img src="images/rulegraph.png" alt="workflow" width="800"/>

## Workflow usage

1. modify config.yaml
	-	path.user_root
	-	path.meta_path
	-	samples
	-	dbs.human
	-	database_dir
	-	host
	-	upload_tag
	-	kneaddata_opt.trimmomatic
	-	megahit_set
	-	genome_filter_criteria
	-	filter_chimieric_bins
2. modify Snakefile to chose some rules
3. confirm some remove rules