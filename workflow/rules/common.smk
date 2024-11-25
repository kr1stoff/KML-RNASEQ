from snakemake.utils import validate
import pandas as pd
from pathlib import Path


configfile: str(Path(workflow.basedir).parent.joinpath("config/config.yaml"))


# validate config
validate(config, schema="../schemas/config.schema.yaml")

# sample config
df_sample = (
    pd.read_csv(config["samples"], sep="\t", dtype={"name": str, "group": str})
    .set_index("name", drop=False)
    .sort_index()
)

validate(df_sample, schema="../schemas/samples.schema.yaml")
