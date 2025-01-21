DBDIR = config["database_dir"]

localrules:
    checkm2_download_db,
    # download_gunc,

rule checkm2_download_db:
    output:
        directory(f"{DBDIR}/CheckM2"),
    conda:
        "../envs/checkm2.yaml"
    threads: 1
    log:
        "logs/download/checkm2.log",
    resources:
        time_min=60 * int(config.get("runtime", {"long": 10})["long"]),
    shell:
        " checkm2 database --download --path {output} "
        " &>> {log}"

"""
rule download_gunc:
    output:
        directory(f"{DBDIR}/gunc_db"),
    conda:
        "../envs/gunc.yaml"
    threads: 1
    resources:
        time_min=60 * int(config.get("runtime", {"default": 5})["default"]),
        mem_mb=config.get("simplejob_mem", 1) * 1000,
    log:
        "logs/downloads/gunc_download.log",
    shell:
        " rm -rf {output} "
        " ; "
        " mkdir -p {output} "
        " ; "
        " gunc download_db {output} &> {log} "
"""
