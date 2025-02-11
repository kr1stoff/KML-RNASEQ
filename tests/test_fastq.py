from kml_rnaseq import prepare_fastq_by_samptab
from kml_rnaseq import get_sample_names_by_samptab

# * 原始数据
# work_dir = '/data/mengxf/Project/KML250122_rnaseq_ZiYan/result/25021102'
# sample_table = '/data/mengxf/Project/KML250122_rnaseq_ZiYan/work/250122/input.tsv'

# * 抽数据，小样本
work_dir = '/data/mengxf/Project/KML250122_rnaseq_ZiYan/result/250207'
sample_table = '/data/mengxf/Project/KML250122_rnaseq_ZiYan/work/250122/input.1M.tsv'

metadata = '/data/mengxf/Project/KML250122_rnaseq_ZiYan/work/250122/metadata.tsv'


def test_prepare():
    prepare_fastq_by_samptab(work_dir, sample_table)


def test_get():
    names = get_sample_names_by_samptab(sample_table)
    print(names)
