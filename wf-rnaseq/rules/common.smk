from snakemake.utils import validate

# validate config
validate(config, schema="../schemas/config.schema.yaml")
