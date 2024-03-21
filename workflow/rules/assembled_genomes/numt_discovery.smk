rule align_circular_mtgenome:
    input:
        indexed_genome="resources/genomes/{species}/{assembly_name}/lastdb_basic/{assembly_accession}_{assembly_name}_genomic.prj",
        mtc_genome="data/mtc_genomes/{species}/mtc_genome_circularized.fasta",
    output:
        maf="data/numt_discovery/lastal/circular_sequence/{species}/{assembly_name}/{assembly_accession}_{assembly_name}.maf",
    params:
        lastal_score_threshold=30,
        multiplicity=150000
    threads: 4
    shell:
        "lastal $(echo {input.indexed_genome} | sed 's/\.[^.]*$//') -e {params.lastal_score_threshold} -P {threads} -m{params.multiplicity} {input.mtc_genome} | last-postmask > {output.maf}"

rule align_linear_mtgenome:
    input:
        indexed_genome="resources/genomes/{species}/{assembly_name}/lastdb_basic/{assembly_accession}_{assembly_name}_genomic.prj",
        mtc_genome="data/mtc_genomes/{species}/mtc_genome.fasta",
    output:
        maf="data/numt_discovery/lastal/linear_sequence/{species}/{assembly_name}/{assembly_accession}_{assembly_name}.maf",
    params:
        lastal_score_threshold=30,
        multiplicity=150000
    threads: 4
    shell:
        "lastal $(echo {input.indexed_genome} | sed 's/\.[^.]*$//') -e {params.lastal_score_threshold} -P {threads} -m{params.multiplicity} {input.mtc_genome} | last-postmask> {output.maf}"


rule convert_maf_to_blasttab:
    input:
        linear_maf="data/numt_discovery/lastal/linear_sequence/{species}/{assembly_name}/{assembly_accession}_{assembly_name}.maf",
        circular_maf="data/numt_discovery/lastal/circular_sequence/{species}/{assembly_name}/{assembly_accession}_{assembly_name}.maf",
    output:
        linear_blasttab="data/numt_discovery/lastal/linear_sequence/{species}/{assembly_name}/{assembly_accession}_{assembly_name}.blasttab",
        circular_blasttab="data/numt_discovery/lastal/circular_sequence/{species}/{assembly_name}/{assembly_accession}_{assembly_name}.blasttab",
    shell:
        "maf-convert blasttab {input.linear_maf} > {output.linear_blasttab} && maf-convert blasttab {input.circular_maf} > {output.circular_blasttab}"


rule replace_in_linear_df_with_alignments_over_edge:
    input:
        circular_maf="data/numt_discovery/lastal/circular_sequence/{species}/{assembly_name}/{assembly_accession}_{assembly_name}.maf",
        circular_BlastTab="data/numt_discovery/lastal/circular_sequence/{species}/{assembly_name}/{assembly_accession}_{assembly_name}.blasttab",
        linear_maf="data/numt_discovery/lastal/linear_sequence/{species}/{assembly_name}/{assembly_accession}_{assembly_name}.maf",
        linear_BlastTab="data/numt_discovery/lastal/linear_sequence/{species}/{assembly_name}/{assembly_accession}_{assembly_name}.blasttab",
    output:
        alignments_to_drop="data/numt_discovery/lastal/dropped_from_linear/{species}/{assembly_name}/{assembly_accession}_{assembly_name}/dropped_from_linear.csv",
        alignments_to_add="data/numt_discovery/lastal/added_to_linear/{species}/{assembly_name}/{assembly_accession}_{assembly_name}/added_to_linear.csv",
        merged_alignments="data/numt_discovery/lastal/linear_and_circular_merged/{species}/{assembly_name}/{assembly_accession}_{assembly_name}/linear_and_circular_merged.csv",
    params:
        functions_path=config["functions_path"],
    script:
        "../../scripts/assembled_genomes/numt_discovery/alignments_over_edge.py"

rule drop_duplicates_and_define_regions_and_names_assembly_numts:
    input:
        input_data="data/numt_discovery/lastal/linear_and_circular_merged/{species}/{assembly_name}/{assembly_accession}_{assembly_name}/linear_and_circular_merged.csv",
        report="resources/genomes/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_assembly_report.tsv",
        mtc_scaffolds="stats/general/mtc_scaffolds.tsv",
    output:
        numts_table="results/numt_discovery/numts/{species}/{assembly_name}/{assembly_accession}_{assembly_name}_numts.csv",
        dropped_putative_numts_close_to_extremity="stats/numt_discovery/dropped_numts_close_to_extremity/{species}/{assembly_name}/{assembly_accession}_{assembly_name}.csv",
        duplicate_numts="stats/numt_discovery/dropped_duplicate_numts/{species}/{assembly_name}/{assembly_accession}_{assembly_name}.csv",
        #duplicate_numts_plots_directory=directory("figures/numt_discovery/dropped_duplicate_numts/{species}/{assembly_name}/{assembly_accession}_{assembly_name}/"),
    params:
        species=lambda wc: wc.species,
        assembly_name=lambda wc: wc.assembly_name,
        functions_path=config["functions_path"],
    script:
        "../../scripts/assembled_genomes/numt_discovery/assembly_numts.py"

rule stats_and_plots_numt_discovery:
    input:
        numt_paths=expand(
            "results/numt_discovery/numts/{species}/{assembly_name}/{assembly_accession}_{assembly_name}_numts.csv",
            zip,
            species=genomes["Species"],
            assembly_accession=genomes["Assembly Accession"],
            assembly_name=genomes["Assembly Name"],
        ),
        genome_info="config/genomes.tsv",
    params:
        region_plots_path="figures/numt_discovery/fragmented_numt_regions/",
        distribution_plots_path_individual="figures/numt_discovery/numt_size_and_identity_distribution_individual/",
        distribution_plots_path_overall="figures/numt_discovery/numt_size_and_identity_distribution_overall/",
        functions_path=config["functions_path"],
    output:
        "tables/numt_discovery/all_genomes_numt_info.csv",
    script:
        "../../scripts/assembled_genomes/numt_discovery/stats_and_figures.py"


rule prepare_numt_region_bed:
    input:
        numts_table="results/numt_discovery/numts/{species}/{assembly_name}/{assembly_accession}_{assembly_name}_numts.csv",
    output:
        bed="data/numt_discovery/region_sequences/{species}/{assembly_name}/{assembly_accession}_{assembly_name}.bed",
    run:
        pd.read_csv(input.numts_table)[
            ["chromosome_code", "region_start", "region_end", "region_id"]
        ].drop_duplicates().to_csv(output.bed, sep="\t", index=False, header=None)
