rule get_mtc_scaffolds_from_assembly_reports:
    input:
        reports=expand(
            "resources/genomes/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_assembly_report.tsv",
            zip,
            species=genomes["Species"],
            assembly_name=genomes["Assembly Name"],
            assembly_accession=genomes["Assembly Accession"],
        ),
        genomes="config/genomes.tsv",
    output:
        "stats/general/mtc_scaffolds.tsv",
    script:
        "../scripts/general/get_mtc_scaffolds_from_assembly_reports.py"
