import pandas as pd
from pybedtools import BedTool
import numpy as np
import seaborn as sns

# ref=pd.read_csv("results/numt_discovery/numts/Sus_scrofa/Sscrofa11.1/GCF_000003025.6_Sscrofa11.1_numts.csv")
ref = pd.read_csv(snakemake.input.ref_numts)
ref = ref[["region_id"]].drop_duplicates()
ref.columns = ["REF11_region_id"]

# genomes=pd.read_table("config/genomes.tsv")
genomes = pd.read_table(snakemake.input.genomes)

dfs = {}
for index, row in genomes.iloc[1:-2].iterrows():
    other_species = row["Species"]
    other_assembly_name = row["Assembly Name"]
    other_assembly_accession = row["Assembly Accession"]
    tempdf = pd.read_csv(
        f"test/genome_comparison/regions_extracted_from_maf/Sscrofa11.1/{other_species}/{other_assembly_name}/both_sides_compared/GCF_000003025.6_vs_{other_assembly_accession}.csv"
    )
    tempdf["other_status"] = (
        tempdf["main_region_status"] + "|" + tempdf["other_region_status"]
    )
    tempdf = tempdf.drop(columns=["main_region_status", "other_region_status"])
    for i in [
        "main_region_start",
        "main_region_end",
        "other_region_start",
        "other_region_end",
    ]:
        tempdf[i] = tempdf[i].astype(int)
    tempdf["main_region_size"] = tempdf["main_region_end"] - tempdf["main_region_start"]
    tempdf["other_region_size"] = (
        tempdf["other_region_end"] - tempdf["other_region_start"]
    )
    tempdf = tempdf.sort_values(by="other_region_size", ascending=False)
    tempdf.columns = [i.replace("main", "REF11") for i in tempdf.columns]
    tempdf.columns = [
        i.replace("other", row["Genome Abbreviation"]) for i in tempdf.columns
    ]
    dfs[row["Genome Abbreviation"]] = tempdf


df = pd.concat(dfs)

df_noref = df[df["REF11_region_id"] == "."]
df = df[df["REF11_region_id"] != "."]


df_ref = []
for name, group in df.groupby(by="REF11_region_id"):
    df_ref.append(group.groupby("REF11_region_id").first())
df_ref = pd.concat(df_ref)
df_ref = df_ref.reset_index()

otherdf = df_noref

otherdf = otherdf.sort_values(
    by=["REF11_chromosome", "REF11_region_start", "REF11_region_end"]
)

otherdf_bed = BedTool.from_dataframe(
    otherdf[["REF11_chromosome", "REF11_region_start", "REF11_region_end"]]
)

otherdf_clustered = otherdf_bed.cluster().to_dataframe()

otherdf["cluster"] = otherdf_clustered["name"].tolist()

merged_otherdf = []

for name, group in otherdf.groupby(by="cluster"):
    merged_otherdf.append(group.groupby("cluster").first())
merged_otherdf = pd.concat(merged_otherdf)
merged_otherdf = merged_otherdf.reset_index(drop=True)
outdf = pd.concat([df_ref, merged_otherdf])

# report=pd.read_table("resources/genomes/Sus_scrofa/Sscrofa11.1/ncbi/GCF_000003025.6_Sscrofa11.1_assembly_report.tsv")
report = pd.read_table(snakemake.input.report)
chromdict = {
    i: j for i, j in zip(report["RefSeq-Accn"], report["Chromosome/Scaffold Name"])
}

for k, v in chromdict.items():
    if "NC_" in k:
        chromdict[k] = "SSC" + v

outdf["REF11_chromosome_name"] = outdf["REF11_chromosome"].apply(lambda x: chromdict[x])

outdf = outdf.sort_values(by=["REF11_chromosome", "REF11_region_start"])

outdf = outdf.reset_index(drop=True)
outdf = outdf.reset_index()
region_names = {}

for name, group in outdf.groupby(by="REF11_chromosome"):
    counter = 0
    for index, row in group.iterrows():
        counter += 1
        region_names[index] = f"{row['REF11_chromosome_name']}_REGION_{counter}"


outdf["region_id"] = outdf["index"].apply(lambda x: region_names[x])

outdf = outdf.drop(columns="index")


statuscols = [i for i in df.columns if "_status" in i]

ref_status = []
for index, row in outdf.iterrows():
    tempstatus = None
    tempstatus_other = None
    for col in statuscols:
        try:
            if "PRIVATE|" in row[col]:
                tempstatus = "PRIVATE"
            elif "COMPATIBLE_SEQUENCE|" in row[col]:
                tempstatus = "COMPATIBLE_SEQUENCE"
            else:
                tempstatus = "CORRESPONDING"
        except:
            pass
        try:
            if "|PRIVATE" in row[col]:
                tempstatus_other = "PRIVATE"
            elif "|CORRESPONDING" in row[col]:
                tempstatus_other = "CORRESPONDING"
            elif "|COMPATIBLE_SEQUENCE" in row[col]:
                tempstatus_other = "COMPATIBLE_SEQUENCE"
            else:
                tempstatus_other = "PORCOLADRO"
        except:
            pass
    ref_status.append(f"{tempstatus}|{tempstatus_other}")

outdf["REF11_status"] = ref_status

newcolumns = ["region_id", "REF11_chromosome_name"] + [
    i for i in outdf.columns if ("REF11" in i) & ("chromosome_name" not in i)
]

for i in outdf.columns:
    if i not in newcolumns:
        newcolumns.append(i)

outdf = outdf[newcolumns]


statuscols = [i for i in outdf.columns if "_status" in i]

ref11_seqs = outdf[outdf["REF11_status"].str.contains("COMPATIBLE_SEQUENCE\|")]
ref11_seqs = ref11_seqs[["REF11_chromosome", "REF11_region_start", "REF11_region_end"]]
outdf["REF11_region_id"] = (
    outdf["REF11_chromosome"].astype(str)
    + ":"
    + outdf["REF11_region_start"].astype(str)
    + "-"
    + outdf["REF11_region_end"].astype(str)
).where(
    outdf["REF11_status"].str.contains("COMPATIBLE_SEQUENCE\|"),
    outdf["REF11_region_id"],
)

regioncols = [i for i in outdf.columns if "_region_id" in i]

genomes = genomes["Genome Abbreviation"].iloc[0:-2]

for i in genomes[1::]:
    outdf[f"{i}_region_id"] = (
        outdf[f"{i}_chromosome"].astype(str)
        + ":"
        + outdf[f"{i}_region_start"].fillna(0).astype(int).astype(str)
        + "-"
        + outdf[f"{i}_region_end"].fillna(0).astype(int).astype(str)
    ).where(
        outdf[f"{i}_status"].str.contains("\|COMPATIBLE_SEQUENCE"),
        outdf[f"{i}_region_id"],
    )


outdf.to_csv(snakemake.output[0], index=False)
