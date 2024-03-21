rule get_additional_compatible_sequences_bed:
    input:
        regions="results/genome_comparison/all_numts_updated.csv",
        genomes="config/genomes.tsv"
    output:
        "data/genome_comparison/additional_regions_compatible_seqs/{species}/{assembly_name}/{assembly_accession}/{genome}.bed"
    params:
        genome=lambda wc:wc.genome,
        species=lambda wc:wc.species,
        assembly_name=lambda wc:wc.assembly_name,
        assembly_accession=lambda wc:wc.assembly_accession
    script:
        "../../scripts/assembled_genomes/numt_assembly_comparison/get_compatible_seqs_bed.py"

rule get_additional_compatible_sequences_fasta:
    input:
        bed="data/genome_comparison/additional_regions_compatible_seqs/{species}/{assembly_name}/{assembly_accession}/{genome}.bed",
        genome="resources/genomes/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_genomic.fna.gz"
    output:
        "data/genome_comparison/additional_regions_compatible_seqs/{species}/{assembly_name}/{assembly_accession}/{genome}.fasta"
    shell:
        "bedtools getfasta -fi {input.genome} -bed {input.bed} -fo {output}"

rule lastdb_sequences:
    input:
        "data/genome_comparison/additional_regions_compatible_seqs/{species}/{assembly_name}/{assembly_accession}/{genome}.fasta"
    output:
        "data/genome_comparison/additional_regions_compatible_seqs/{species}/{assembly_name}/{assembly_accession}/{genome}/lastdb/lastdb.prj"
    params:
        db_name="data/genome_comparison/additional_regions_compatible_seqs/{species}/{assembly_name}/{assembly_accession}/{genome}/lastdb/lastdb"
    shell:
        "lastdb {params.db_name} {input}"

rule align_mtc_to_additional_seqs_lastdb:
    input:
        db="data/genome_comparison/additional_regions_compatible_seqs/{species}/{assembly_name}/{assembly_accession}/{genome}/lastdb/lastdb.prj",
        mtc_genome="data/mtc_genomes/{species}/mtc_genome.fasta",
    output:
        maf="data/genome_comparison/additional_regions_compatible_seqs/{species}/{assembly_name}/{assembly_accession}/{genome}/lastal/{genome}.maf",
    params:
        multiplicity=150000
    threads: 4
    shell:
        "lastal $(echo {input.db} | sed 's/\.[^.]*$//') -P {threads} -m{params.multiplicity} {input.mtc_genome} | last-postmask > {output.maf}"

rule additional_seqs_maf_to_blasttab:
    input:
        "data/genome_comparison/additional_regions_compatible_seqs/{species}/{assembly_name}/{assembly_accession}/{genome}/lastal/{genome}.maf"
    output:
        "data/genome_comparison/additional_regions_compatible_seqs/{species}/{assembly_name}/{assembly_accession}/{genome}/lastal/{genome}.blasttab"
    shell:
        "maf-convert blasttab {input} > {output}"

rule add_additional_seqs_to_original_numts:
    input:
        report="resources/genomes/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_assembly_report.tsv",
        mtc_scaffolds="stats/general/mtc_scaffolds.tsv",
        maf="data/genome_comparison/additional_regions_compatible_seqs/{species}/{assembly_name}/{assembly_accession}/{genome}/lastal/{genome}.maf",
        bt="data/genome_comparison/additional_regions_compatible_seqs/{species}/{assembly_name}/{assembly_accession}/{genome}/lastal/{genome}.blasttab",
        numts="results/numt_discovery/numts/{species}/{assembly_name}/{assembly_accession}_{assembly_name}_numts.csv",
        additional_numts_bed="data/genome_comparison/additional_regions_compatible_seqs/{species}/{assembly_name}/{assembly_accession}/{genome}.bed"
    output:
        "data/genome_comparison/additional_regions_compatible_seqs/{species}/{assembly_name}/{assembly_accession}/{genome}/all_numts/{genome}.csv"
    params:
        functions_path=config["functions_path"],
        species=lambda wc: wc.species,
    script:
        "../../scripts/assembled_genomes/numt_assembly_comparison/merge_compatible_seqs.py"

rule updated_stats:
    input:
        numt_paths=expand(
            "data/genome_comparison/additional_regions_compatible_seqs/{species}/{assembly_name}/{assembly_accession}/{genome_id}/all_numts/{genome_id}.csv",
            zip,
            species=genomes[~genomes['Species'].isin(['Bos_taurus','Capra_hircus'])]["Species"],
            assembly_accession=genomes[~genomes['Species'].isin(['Bos_taurus','Capra_hircus'])]["Assembly Accession"],
            assembly_name=genomes[~genomes['Species'].isin(['Bos_taurus','Capra_hircus'])]["Assembly Name"],
            genome_id=genomes[~genomes['Species'].isin(['Bos_taurus','Capra_hircus'])]['Genome Abbreviation']
        ),
        genome_info="config/genomes.tsv",
    params:
        distribution_plots_path_individual="figures/numts_and_additional_regions/numt_size_and_identity_distribution_individual/",
        distribution_plots_path_overall="figures/numts_and_additional_regions/numt_size_and_identity_distribution_overall/",
        functions_path=config["functions_path"],
    output:
        "tables/numts_and_additional_regions/all_genomes_numt_info.csv",
    script:
        "../../scripts/assembled_genomes/numt_assembly_comparison/stats_and_figures_with_additional_regions.py"
