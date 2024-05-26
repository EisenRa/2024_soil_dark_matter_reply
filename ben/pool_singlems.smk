import polars as pl

bioprojects_and_archives = pl.read_csv('bioproject_joined.tsv', separator='\t')

sample_to_archive = {}
for row in bioprojects_and_archives.rows(named=True):
    sample_to_archive[row['sample_accession']] = row['archive'].split(' ')
sample_names = list(sample_to_archive.keys())

print("Read in {} samples".format(len(sample_to_archive)))

singlem = '/home/woodcrob/git/singlem/bin/singlem'
metapackage = '/work/microbiome/db/singlem/S4.2.2.GTDB_r214.metapackage_20240502.smpkg.zb'

rule all:
    input:
        'snakedata/smf.csv',

rule merge:
    output:
        merged = 'snakedata/merged_archives/{sample}.json',
        done = touch('snakedata/dumped_reads/{sample}.done'),
    params:
        archives = lambda wildcards: ' '.join([f'<(zcat {j})' for j in sample_to_archive[wildcards.sample]])
    conda:
        "singlem-dev"
    shell:
        """
        {singlem} summarise --input-archive-otu-table {params.archives} --output-archive-otu-table {output.merged} --collapse-to-sample-name {wildcards.sample}
        """

rule renew:
    input:
        json = 'snakedata/merged_archives/{sample}.json'
    output:
        json = 'snakedata/renewed_archives/{sample}.json',
        done = touch('snakedata/renewed_archives/{sample}.done'),
        profile = 'snakedata/renewed_archives/{sample}.profile'
    conda:
        "singlem-dev"
    shell:
        """
        {singlem} renew --input-archive-otu-table {input.json} --archive-otu-table {output.json} --metapackage {metapackage} --taxonomic-profile {output.profile} --taxonomic-profile-krona {output.profile}.html
        """

rule merge_renewed_profiles:
    input:
        profiles = expand('snakedata/renewed_archives/{sample}.profile', sample=sample_names)
    output:
        merged = 'snakedata/merged_profiles.csv'
    conda:
        "singlem-dev"
    shell:
        """
        {singlem} summarise --input-taxonomic-profile {input.profiles} --output-taxonomic-profile {output.merged}
        """

rule smf:
    input:
        profile = 'snakedata/merged_profiles.csv',
    output:
        smf = 'snakedata/smf.csv',
    conda:
        "singlem-dev"
    shell:
        """
        {singlem} microbial_fraction --input-profile {input.profile} --output-tsv {output.smf} --metapackage {metapackage} --input-metagenome-sizes num_bases.csv
        """
