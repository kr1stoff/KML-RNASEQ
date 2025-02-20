import logging
import click
from pathlib import Path
from kml_rnaseq import prepare_fastq_by_samptab
from kml_rnaseq import get_sample_names_by_samptab
from kml_rnaseq import create_snakemake_configfile
from kml_rnaseq import run_snakemake


logging.basicConfig(level=logging.DEBUG,
                    format="%(asctime)s - %(levelname)s - %(filename)s - %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S")


@click.command()
@click.option('--sample_table', '-s', type=click.Path(exists=True), required=True, help='样本信息表.')
@click.option('--metadata', '-m', type=click.Path(exists=True), required=True, help='metadata 元数据表.')
@click.option('--work_dir', '-w', type=str, default='lvisa_result', help='结果生成目录. [default: lvisa_result]')
@click.help_option('-h', '--help')
def main(work_dir, sample_table, metadata):
    """RNA-seq 分析流程"""
    logging.info(f'开始分析!')

    metadata = str(Path(metadata).resolve())
    sample_table = str(Path(sample_table).resolve())
    work_dir = Path(work_dir).resolve()

    # fastq
    sample_names = get_sample_names_by_samptab(sample_table)
    prepare_fastq_by_samptab(work_dir, sample_table)

    # snakemake
    create_snakemake_configfile(sample_names, work_dir, metadata=metadata)
    run_snakemake(work_dir)

    logging.info(f'分析完成!')


if __name__ == '__main__':
    main()
