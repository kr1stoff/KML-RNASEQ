from kml_rnaseq import create_snakemake_configfile
from kml_rnaseq import get_sample_names_by_samptab
# from kml_rnaseq import run_snakemake

sample_table = '/data/mengxf/Project/KML250122-rnaseq-ZiYan/work/250122/input.tsv'
metadata = '/data/mengxf/Project/KML250122-rnaseq-ZiYan/work/250122/metadata.tsv'
work_dir = '/data/mengxf/Project/KML250122-rnaseq-ZiYan/result/250731'

sample_names = get_sample_names_by_samptab(sample_table)
print(create_snakemake_configfile(sample_names, work_dir, metadata))
