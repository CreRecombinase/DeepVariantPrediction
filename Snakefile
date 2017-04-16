configfile:
    'config.yaml'
include:
    'modules/variant2input.snakemake',
    'modules/input2score.snakemake',
    'modules/score2performance.snakemake'

rule all:
    input:
        expand('performance/{group}/{method}/report.html', method=config['performance'], group=config['label']),
        expand('score/{group}/result.rds', group=config['label'])
