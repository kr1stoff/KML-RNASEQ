from kml_rnaseq import create_snakemake_configfile
from kml_rnaseq import get_sample_names_by_samptab
# from kml_rnaseq import run_snakemake


work_dir = '/data/mengxf/Project/KML250113_RNAseq_pipeline/result/250115'
sample_table = '/data/mengxf/Project/KML250113_RNAseq_pipeline/input/250115.input.tsv'
metadata = '/data/mengxf/Project/KML250113_RNAseq_pipeline/input/250115.metadata.tsv'


def test_create():
    sample_names = get_sample_names_by_samptab(sample_table)
    print(create_snakemake_configfile(sample_names, work_dir, metadata))


# def test_run():
#     run_snakemake(work_dir)
