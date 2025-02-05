from kml_rnaseq import prepare_fastq_by_samptab
from kml_rnaseq import get_sample_names_by_samptab

work_dir = '/data/mengxf/Project/KML250122_rnaseqe_ZiYan/result/250123'
sample_table = '/data/mengxf/Project/KML250122_rnaseqe_ZiYan/jupyter/250122/input.1M.tsv'
metadata = '/data/mengxf/Project/KML250122_rnaseqe_ZiYan/jupyter/250122/metadata.tsv'


def test_prepare():
    prepare_fastq_by_samptab(work_dir, sample_table)


def test_get():
    names = get_sample_names_by_samptab(sample_table)
    print(names)
