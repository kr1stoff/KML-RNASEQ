from kml_rnaseq import create_snakemake_configfile
from kml_rnaseq import get_sample_names_by_samptab
# from kml_rnaseq import run_snakemake


work_dir = '/data/mengxf/Project/KML250122_rnaseq_ZiYan/result/250207'
sample_table = '/data/mengxf/Project/KML250122_rnaseq_ZiYan/work/250122/input.1M.tsv'
metadata = '/data/mengxf/Project/KML250122_rnaseq_ZiYan/work/250122/metadata.tsv'


def test_create():
    sample_names = get_sample_names_by_samptab(sample_table)
    print(create_snakemake_configfile(sample_names, work_dir, metadata))


# def test_run():
#     run_snakemake(work_dir)
