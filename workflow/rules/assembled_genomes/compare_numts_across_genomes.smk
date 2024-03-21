rule maf_cut_extract_numt_region_from_alignment:
    input:
        maf="resources/genomes/Sus_scrofa/Sscrofa11.1/last_split_near/{other_species}/{other_assembly_name}/GCF_000003025.6_{other_assembly_accession}.split.sorted.maf",
        region_bed="data/numt_discovery/region_sequences/{other_species}/{other_assembly_name}/{other_assembly_accession}_{other_assembly_name}.bed",
    output:
        "data/genome_comparison/regions_extracted_from_maf/Sscrofa11.1/{other_species}/{other_assembly_name}/{other_assembly_accession}_in_GCF_000003025.6.maf",
    shell:
        "touch {output} > {output} && cut -f 1,2,3 {input.region_bed} |sed 's/\t/:/' | sed 's/\t/-/' | for i in $(cat -) ; do echo $i && maf-cut ${{i}} {input.maf} >> {output}; done"


rule maf_cut_extract_reference_regions_from_alignment:
    input:
        maf="resources/genomes/Sus_scrofa/Sscrofa11.1/last_split_near/{other_species}/{other_assembly_name}/GCF_000003025.6_{other_assembly_accession}.split.sorted.maf",
        region_bed_ref="data/numt_discovery/region_sequences/Sus_scrofa/Sscrofa11.1/GCF_000003025.6_Sscrofa11.1.bed",
    output:
        "data/genome_comparison/regions_extracted_from_maf/Sscrofa11.1/{other_species}/{other_assembly_name}/GCF_000003025.6_in_{other_assembly_accession}.maf",
    shell:
        "touch {output} > {output} && cut -f 1,2,3 {input.region_bed_ref} |sed 's/\t/:/' | sed 's/\t/-/' | for i in $(cat -) ; do echo $i && maf-cut ${{i}} {input.maf} >> {output}; done"


rule maf_cut_extract_numt_region_from_reverse_alignment:
    input:
        maf="resources/genomes/{other_species}/{other_assembly_name}/last_split_near/Sus_scrofa/Sscrofa11.1/{other_assembly_accession}_GCF_000003025.6.split.sorted.maf",
        region_bed="data/numt_discovery/region_sequences/{other_species}/{other_assembly_name}/{other_assembly_accession}_{other_assembly_name}.bed",
    output:
        "data/genome_comparison/regions_extracted_from_maf/{other_species}/{other_assembly_name}/{other_assembly_accession}_in_GCF_000003025.6.maf",
    shell:
        "touch {output} > {output} && cut -f 1,2,3 {input.region_bed} |sed 's/\t/:/' | sed 's/\t/-/' | for i in $(cat -) ; do echo $i && maf-cut ${{i}} {input.maf} >> {output}; done"


rule maf_cut_extract_reference_regions_from_reverse_alignment:
    input:
        maf="resources/genomes/{other_species}/{other_assembly_name}/last_split_near/Sus_scrofa/Sscrofa11.1/{other_assembly_accession}_GCF_000003025.6.split.sorted.maf",
        region_bed_ref="data/numt_discovery/region_sequences/Sus_scrofa/Sscrofa11.1/GCF_000003025.6_Sscrofa11.1.bed",
    output:
        "data/genome_comparison/regions_extracted_from_maf/{other_species}/{other_assembly_name}/GCF_000003025.6_in_{other_assembly_accession}.maf",
    shell:
        "touch {output} > {output} && cut -f 1,2,3 {input.region_bed_ref} |sed 's/\t/:/' | sed 's/\t/-/' | for i in $(cat -) ; do echo $i && maf-cut ${{i}} {input.maf} >> {output}; done"


rule convert_cut_regions_maf_to_blasttab:
    input:
        a="data/genome_comparison/regions_extracted_from_maf/Sscrofa11.1/{other_species}/{other_assembly_name}/{other_assembly_accession}_in_GCF_000003025.6.maf",
        b="data/genome_comparison/regions_extracted_from_maf/Sscrofa11.1/{other_species}/{other_assembly_name}/GCF_000003025.6_in_{other_assembly_accession}.maf",
    output:
        a="data/genome_comparison/regions_extracted_from_maf/Sscrofa11.1/{other_species}/{other_assembly_name}/{other_assembly_accession}_in_GCF_000003025.6.blasttab",
        b="data/genome_comparison/regions_extracted_from_maf/Sscrofa11.1/{other_species}/{other_assembly_name}/GCF_000003025.6_in_{other_assembly_accession}.blasttab",
    shell:
        "maf-convert blasttab {input.a} > {output.a} && maf-convert blasttab {input.b} > {output.b}  "


rule convert_cut_regions_reverse_maf_to_blasttab:
    input:
        a="data/genome_comparison/regions_extracted_from_maf/{other_species}/{other_assembly_name}/{other_assembly_accession}_in_GCF_000003025.6.maf",
        b="data/genome_comparison/regions_extracted_from_maf/{other_species}/{other_assembly_name}/GCF_000003025.6_in_{other_assembly_accession}.maf",
    output:
        a="data/genome_comparison/regions_extracted_from_maf/{other_species}/{other_assembly_name}/{other_assembly_accession}_in_GCF_000003025.6.blasttab",
        b="data/genome_comparison/regions_extracted_from_maf/{other_species}/{other_assembly_name}/GCF_000003025.6_in_{other_assembly_accession}.blasttab",
    shell:
        "maf-convert blasttab {input.a} > {output.a} && maf-convert blasttab {input.b} > {output.b}  "


rule other_in_reference:
    input:
        main_numts="results/numt_discovery/numts/{other_species}/{other_assembly_name}/{other_assembly_accession}_{other_assembly_name}_numts.csv",
        other_numts="results/numt_discovery/numts/Sus_scrofa/Sscrofa11.1/GCF_000003025.6_Sscrofa11.1_numts.csv",
        maf="data/genome_comparison/regions_extracted_from_maf/{other_species}/{other_assembly_name}/{other_assembly_accession}_in_GCF_000003025.6.maf",
        blasttab="data/genome_comparison/regions_extracted_from_maf/{other_species}/{other_assembly_name}/{other_assembly_accession}_in_GCF_000003025.6.blasttab",
        main_regions="data/numt_discovery/region_sequences/{other_species}/{other_assembly_name}/{other_assembly_accession}_{other_assembly_name}.bed",
        other_regions="data/numt_discovery/region_sequences/Sus_scrofa/Sscrofa11.1/GCF_000003025.6_Sscrofa11.1.bed",
    output:
        comparison="data/genome_comparison/regions_extracted_from_maf/{other_species}/{other_assembly_name}/{other_assembly_accession}_in_GCF_000003025.6.merged.csv",
    script:
        "../../scripts/assembled_genomes/numt_assembly_comparison/maf_merge_numts.py"


rule reference_in_other:
    input:
        main_numts="results/numt_discovery/numts/Sus_scrofa/Sscrofa11.1/GCF_000003025.6_Sscrofa11.1_numts.csv",
        other_numts="results/numt_discovery/numts/{other_species}/{other_assembly_name}/{other_assembly_accession}_{other_assembly_name}_numts.csv",
        maf="data/genome_comparison/regions_extracted_from_maf/Sscrofa11.1/{other_species}/{other_assembly_name}/GCF_000003025.6_in_{other_assembly_accession}.maf",
        blasttab="data/genome_comparison/regions_extracted_from_maf/Sscrofa11.1/{other_species}/{other_assembly_name}/GCF_000003025.6_in_{other_assembly_accession}.blasttab",
        main_regions="data/numt_discovery/region_sequences/Sus_scrofa/Sscrofa11.1/GCF_000003025.6_Sscrofa11.1.bed",
        other_regions="data/numt_discovery/region_sequences/{other_species}/{other_assembly_name}/{other_assembly_accession}_{other_assembly_name}.bed",
    output:
        comparison="data/genome_comparison/regions_extracted_from_maf/Sscrofa11.1/{other_species}/{other_assembly_name}/GCF_000003025.6_in_{other_assembly_accession}.merged.csv",
    script:
        "../../scripts/assembled_genomes/numt_assembly_comparison/maf_merge_numts.py"


rule other_in_reference_reverse:
    input:
        main_numts="results/numt_discovery/numts/{other_species}/{other_assembly_name}/{other_assembly_accession}_{other_assembly_name}_numts.csv",
        other_numts="results/numt_discovery/numts/Sus_scrofa/Sscrofa11.1/GCF_000003025.6_Sscrofa11.1_numts.csv",
        maf="data/genome_comparison/regions_extracted_from_maf/{other_species}/{other_assembly_name}/GCF_000003025.6_in_{other_assembly_accession}.maf",
        blasttab="data/genome_comparison/regions_extracted_from_maf/{other_species}/{other_assembly_name}/GCF_000003025.6_in_{other_assembly_accession}.blasttab",
        main_regions="data/numt_discovery/region_sequences/{other_species}/{other_assembly_name}/{other_assembly_accession}_{other_assembly_name}.bed",
        other_regions="data/numt_discovery/region_sequences/Sus_scrofa/Sscrofa11.1/GCF_000003025.6_Sscrofa11.1.bed",
    output:
        comparison="data/genome_comparison/regions_extracted_from_maf/{other_species}/{other_assembly_name}/GCF_000003025.6_in_{other_assembly_accession}.merged.csv",
    script:
        "../../scripts/assembled_genomes/numt_assembly_comparison/maf_merge_numts_reverse.py"


rule reference_in_other_reverse:
    input:
        main_numts="results/numt_discovery/numts/Sus_scrofa/Sscrofa11.1/GCF_000003025.6_Sscrofa11.1_numts.csv",
        other_numts="results/numt_discovery/numts/{other_species}/{other_assembly_name}/{other_assembly_accession}_{other_assembly_name}_numts.csv",
        maf="data/genome_comparison/regions_extracted_from_maf/Sscrofa11.1/{other_species}/{other_assembly_name}/{other_assembly_accession}_in_GCF_000003025.6.maf",
        blasttab="data/genome_comparison/regions_extracted_from_maf/Sscrofa11.1/{other_species}/{other_assembly_name}/{other_assembly_accession}_in_GCF_000003025.6.blasttab",
        main_regions="data/numt_discovery/region_sequences/Sus_scrofa/Sscrofa11.1/GCF_000003025.6_Sscrofa11.1.bed",
        other_regions="data/numt_discovery/region_sequences/{other_species}/{other_assembly_name}/{other_assembly_accession}_{other_assembly_name}.bed",
    output:
        comparison="data/genome_comparison/regions_extracted_from_maf/Sscrofa11.1/{other_species}/{other_assembly_name}/{other_assembly_accession}_in_GCF_000003025.6.merged.csv",
    script:
        "../../scripts/assembled_genomes/numt_assembly_comparison/maf_merge_numts_reverse.py"


rule merge_comparisons_forward_reverse:
    input:
        ref_numts="results/numt_discovery/numts/Sus_scrofa/Sscrofa11.1/GCF_000003025.6_Sscrofa11.1_numts.csv",
        other_numts="results/numt_discovery/numts/{other_species}/{other_assembly_name}/{other_assembly_accession}_{other_assembly_name}_numts.csv",
        a="data/genome_comparison/regions_extracted_from_maf/Sscrofa11.1/{other_species}/{other_assembly_name}/{other_assembly_accession}_in_GCF_000003025.6.merged.csv",
        b="data/genome_comparison/regions_extracted_from_maf/Sscrofa11.1/{other_species}/{other_assembly_name}/GCF_000003025.6_in_{other_assembly_accession}.merged.csv",
        c="data/genome_comparison/regions_extracted_from_maf/{other_species}/{other_assembly_name}/GCF_000003025.6_in_{other_assembly_accession}.merged.csv",
        d="data/genome_comparison/regions_extracted_from_maf/{other_species}/{other_assembly_name}/{other_assembly_accession}_in_GCF_000003025.6.merged.csv",
    output:
        "data/genome_comparison/regions_extracted_from_maf/Sscrofa11.1/{other_species}/{other_assembly_name}/both_sides_compared/GCF_000003025.6_vs_{other_assembly_accession}.csv",
    script:
        "../../scripts/assembled_genomes/numt_assembly_comparison/merge_normal_inverse.py"


rule merge_all_comparisons:
    input:
        report="resources/genomes/Sus_scrofa/Sscrofa11.1/ncbi/GCF_000003025.6_Sscrofa11.1_assembly_report.tsv",
        comparisons=expand(
            "data/genome_comparison/regions_extracted_from_maf/Sscrofa11.1/{other_species}/{other_assembly_name}/both_sides_compared/GCF_000003025.6_vs_{other_assembly_accession}.csv",
            zip,
            other_species=genomes[
                (genomes["Assembly Name"] != "Sscrofa11.1")
                & (~genomes["Species"].isin(["Bos_taurus", "Capra_hircus"]))
            ]["Species"],
            other_assembly_name=genomes[
                (genomes["Assembly Name"] != "Sscrofa11.1")
                & (~genomes["Species"].isin(["Bos_taurus", "Capra_hircus"]))
            ]["Assembly Name"],
            other_assembly_accession=genomes[
                (genomes["Assembly Name"] != "Sscrofa11.1")
                & (~genomes["Species"].isin(["Bos_taurus", "Capra_hircus"]))
            ]["Assembly Accession"],
        ),
        ref_numts="results/numt_discovery/numts/Sus_scrofa/Sscrofa11.1/GCF_000003025.6_Sscrofa11.1_numts.csv",
        genomes="config/genomes.tsv",
    output:
        "results/genome_comparison/all_numts.csv",
    script:
        "../../scripts/assembled_genomes/numt_assembly_comparison/merge_all_both_sides.py"


rule get_orthologous_matrix:
    input:
        df="results/genome_comparison/all_numts_updated.csv",
        genomes="config/genomes.tsv",
    output:
        "tables/genome_comparison/orthologues_matrix.csv",
    script:
        "../../scripts/assembled_genomes/numt_assembly_comparison/orthologues_matrix.py"
