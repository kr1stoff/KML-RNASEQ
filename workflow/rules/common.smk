from snakemake.utils import validate
import pandas as pd
from pathlib import Path


configfile: str(Path(workflow.basedir).parent.joinpath("config/config.yaml"))


# validate config
validate(config, schema="../schemas/config.schema.yaml")

# sample config
samples = (
    pd.read_csv(config["samples"], sep="\t", dtype={"name": str, "group": str})
    .set_index("name", drop=False)
    .sort_index()
)

validate(samples, schema="../schemas/samples.schema.yaml")


wildcard_constraints:
    sample="|".join(samples["name"]),


def get_copy_fq_input(wildcards):
    return samples.loc[wildcards.sample, ["fq1", "fq2"]]


def get_fastqc_input1(wildcards):
    # * wrappers fastqc 仅支持单 fq 输入
    return samples.loc[wildcards.sample, "fq1"]


def get_fastqc_input2(wildcards):
    return samples.loc[wildcards.sample, "fq2"]
