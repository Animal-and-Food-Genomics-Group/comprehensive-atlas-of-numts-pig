import pandas as pd

genomes = pd.read_csv(snakemake.input.genomes, sep="\t")

reportlist = []
for tempreport in snakemake.input.reports:
    for index, row in genomes.iterrows():
        if row["Assembly Accession"] in tempreport:
            tempgenome = row
    report = pd.read_csv(tempreport, sep="\t", comment="#")
    report["Assembly Name"] = tempgenome["Assembly Name"]
    report["Assembly Accession"] = tempgenome["Assembly Accession"]
    reportlist.append(report)

reports = pd.concat(reportlist).reset_index(drop=True)

mtc_indexes = []
mtc_indexes += list(
    reports[reports["Sequence-Name"].str.contains("MT", na=False)].index
)
mtc_indexes += list(
    reports[reports["Assigned-Molecule-Location/Type"].str.contains("Mito")].index
)
mtc_indexes += list(reports[reports["Assigned-Molecule"].str.contains("MT")].index)
mtc_indexes += list(reports[reports["Assembly-Unit"] == "non-nuclear"].index)
mtc_scaffolds = reports.loc[mtc_indexes].drop_duplicates()
mtc_scaffolds.to_csv(snakemake.output[0], sep="\t", index=False)
