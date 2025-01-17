from kml_rnaseq import prepare_fastq_by_samptab
from kml_rnaseq import get_sample_names_by_samptab

work_dir = '/data/mengxf/Project/KML250113_RNAseq_pipeline/result/250115'
sample_table = '/data/mengxf/Project/KML250113_RNAseq_pipeline/input/250115.input.tsv'
metadata = '/data/mengxf/Project/KML250113_RNAseq_pipeline/input/250115.metadata.tsv'


def test_prepare():
    prepare_fastq_by_samptab(work_dir, sample_table)


def test_get():
    names = get_sample_names_by_samptab(sample_table)
    print(names)
