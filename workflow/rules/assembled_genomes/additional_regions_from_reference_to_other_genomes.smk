rule get_new_reference_numt_coords:
    input:
        numts="results/genome_comparison/all_numts.csv"
    output:
        coords="data/genome_comparison/additional_regions_found_in_reference/coords.txt"
    run:
        all=pd.read_csv(input.numts)
        coords=all[all['REF11_region_id'].str.contains(":")]['REF11_region_id'].tolist()
        with open(output.coords, 'w') as f:
            for i in coords:
                f.write(i)
                f.write('\n')

rule maf_cut_new_numt_coords_ref_to_other:
    input:
        maf="resources/genomes/Sus_scrofa/Sscrofa11.1/last_split_near/{other_species}/{other_assembly_name}/GCF_000003025.6_{other_assembly_accession}.split.sorted.maf",
        regions="data/genome_comparison/additional_regions_found_in_reference/coords.txt",
    output:
        "data/genome_comparison/additional_regions_found_in_reference/regions_extracted_from_maf/Sscrofa11.1/{other_species}/{other_assembly_name}/{other_assembly_accession}_in_GCF_000003025.6.maf"
    shell:
        "touch {output} > {output} &&  for i in $(cat {input.regions}) ; do echo $i && maf-cut ${{i}} {input.maf} >> {output}; done"

rule maf_cut_new_numt_coords_other_to_ref_reverse:
    input:
        maf="resources/genomes/{other_species}/{other_assembly_name}/last_split_near/Sus_scrofa/Sscrofa11.1/{other_assembly_accession}_GCF_000003025.6.split.sorted.maf",
        regions="data/genome_comparison/additional_regions_found_in_reference/coords.txt",
    output:
        "data/genome_comparison/additional_regions_found_in_reference/regions_extracted_from_maf/Sscrofa11.1/{other_species}/{other_assembly_name}/GCF_000003025.6_in_{other_assembly_accession}.maf"
    shell:
        "touch {output} > {output} &&  for i in $(cat {input.regions}) ; do echo $i && maf-cut ${{i}} {input.maf} >> {output}; done"

rule convert_new_regions_maf_to_blasttab:
    input:
        a="data/genome_comparison/additional_regions_found_in_reference/regions_extracted_from_maf/Sscrofa11.1/{other_species}/{other_assembly_name}/{other_assembly_accession}_in_GCF_000003025.6.maf",
        b="data/genome_comparison/additional_regions_found_in_reference/regions_extracted_from_maf/Sscrofa11.1/{other_species}/{other_assembly_name}/GCF_000003025.6_in_{other_assembly_accession}.maf",
    output:
        a="data/genome_comparison/additional_regions_found_in_reference/regions_extracted_from_maf/Sscrofa11.1/{other_species}/{other_assembly_name}/{other_assembly_accession}_in_GCF_000003025.6.blasttab",
        b="data/genome_comparison/additional_regions_found_in_reference/regions_extracted_from_maf/Sscrofa11.1/{other_species}/{other_assembly_name}/GCF_000003025.6_in_{other_assembly_accession}.blasttab",

    shell:
        "maf-convert blasttab {input.a} > {output.a} && maf-convert blasttab {input.b} > {output.b}  "

rule get_dfs_from_maf_bt:
    input:
        main_maf="data/genome_comparison/additional_regions_found_in_reference/regions_extracted_from_maf/Sscrofa11.1/{other_species}/{other_assembly_name}/{other_assembly_accession}_in_GCF_000003025.6.maf",
        main_bt="data/genome_comparison/additional_regions_found_in_reference/regions_extracted_from_maf/Sscrofa11.1/{other_species}/{other_assembly_name}/{other_assembly_accession}_in_GCF_000003025.6.blasttab",
        rev_maf="data/genome_comparison/additional_regions_found_in_reference/regions_extracted_from_maf/Sscrofa11.1/{other_species}/{other_assembly_name}/GCF_000003025.6_in_{other_assembly_accession}.maf",
        rev_bt="data/genome_comparison/additional_regions_found_in_reference/regions_extracted_from_maf/Sscrofa11.1/{other_species}/{other_assembly_name}/GCF_000003025.6_in_{other_assembly_accession}.blasttab",
    output:
        df="data/genome_comparison/additional_regions_found_in_reference/regions_extracted_from_maf/Sscrofa11.1/{other_species}/{other_assembly_name}/{other_assembly_accession}_in_GCF_000003025.6.csv",
        rev_df="data/genome_comparison/additional_regions_found_in_reference/regions_extracted_from_maf/Sscrofa11.1/{other_species}/{other_assembly_name}/GCF_000003025.6_in_{other_assembly_accession}.csv",
    script:
        "../../scripts/assembled_genomes/additional_reference_numts_in_other_genomes/read_maf.py"


rule update_all_numts:
    input:
        genomes='config/genomes.tsv',
        all_numts="results/genome_comparison/all_numts.csv",
        dfs=expand("data/genome_comparison/additional_regions_found_in_reference/regions_extracted_from_maf/Sscrofa11.1/{other_species}/{other_assembly_name}/{other_assembly_accession}_in_GCF_000003025.6.csv",
                zip,
                other_species=genomes[(genomes['Assembly Name']!='Sscrofa11.1')&(~genomes['Species'].isin(['Bos_taurus','Capra_hircus']))]['Species'],
                other_assembly_name=genomes[(genomes['Assembly Name']!='Sscrofa11.1')&(~genomes['Species'].isin(['Bos_taurus','Capra_hircus']))]['Assembly Name'],
                other_assembly_accession=genomes[(genomes['Assembly Name']!='Sscrofa11.1')&(~genomes['Species'].isin(['Bos_taurus','Capra_hircus']))]['Assembly Accession'])
    output:
        df="results/genome_comparison/all_numts_updated.csv"
    script:
        "../../scripts/assembled_genomes/additional_reference_numts_in_other_genomes/parse_dfs.py"
