import pandas as pd
import pybedtools
from pybedtools import BedTool
import numpy as np

# genomes=pd.read_table("config/genomes.tsv")
genomes = pd.read_table(snakemake.input.genomes)

genomes = genomes.iloc[0:-2]

# all_full=pd.read_csv("results/genome_comparison/all_numts.csv")
all_full = pd.read_csv(snakemake.input.all_numts)
all = all_full[all_full["REF11_region_id"].str.contains(":")]
all_size = {i: j for i, j in zip(all["REF11_region_id"], all["REF11_region_size"])}

all["REF11_strand"] = "+"
all_bed = all[
    [
        "REF11_chromosome",
        "REF11_region_start",
        "REF11_region_end",
        "REF11_strand",
        "REF11_region_id",
    ]
]
all_bed = BedTool.from_dataframe(all_bed)

dfs = {}
revdfs = {}

for index, row in genomes.iloc[1::].iterrows():
    other_species = row["Species"]
    other_assembly_accession = row["Assembly Accession"]
    other_assembly_name = row["Assembly Name"]
    tempdf = f"data/genome_comparison/additional_regions_found_in_reference/regions_extracted_from_maf/Sscrofa11.1/{other_species}/{other_assembly_name}/{other_assembly_accession}_in_GCF_000003025.6.csv"
    temprevdf = f"data/genome_comparison/additional_regions_found_in_reference/regions_extracted_from_maf/Sscrofa11.1/{other_species}/{other_assembly_name}/GCF_000003025.6_in_{other_assembly_accession}.csv"
    dfs[row["Genome Abbreviation"]] = pd.read_csv(tempdf)
    revdfs[row["Genome Abbreviation"]] = pd.read_csv(temprevdf)


df = pd.concat(dfs)
revdf = pd.concat(revdfs)

df = df.reset_index()
df = df.drop(columns=["level_1"])
df = df.rename(columns={"level_0": "genome"})

revdf = revdf.reset_index()
revdf = revdf.drop(columns=["level_1"])
revdf = revdf.rename(columns={"level_0": "genome"})


df_bed = BedTool.from_dataframe(
    df[["REF11_chromosome", "REF11_start", "REF11_end", "REF11_strand", "genome"]]
)

df_merged = df_bed.intersect(all_bed, wao=True).to_dataframe()

revdf_bed = BedTool.from_dataframe(
    revdf[["REF11_chromosome", "REF11_start", "REF11_end", "REF11_strand", "genome"]]
)
revdf_merged = revdf_bed.intersect(all_bed, wao=True).to_dataframe()


df_merged = df_merged[["chrom", "start", "end", "score", "blockCount", "blockSizes"]]
df_merged.columns = [
    "REF11_chromosome",
    "REF11_start",
    "REF11_end",
    "genome",
    "REF11_region_id",
    "overlap",
]
df = df_merged.merge(df, on=list(df_merged.columns[0:4]))

revdf_merged = revdf_merged[revdf_merged["blockCount"] != "."]
revdf_merged = revdf_merged[
    ["chrom", "start", "end", "score", "blockCount", "blockSizes"]
]
revdf_merged.columns = [
    "REF11_chromosome",
    "REF11_start",
    "REF11_end",
    "genome",
    "REF11_region_id",
    "overlap",
]
revdf = revdf_merged.merge(revdf, on=list(revdf_merged.columns[0:4]))

df["REF11_region_size"] = df["REF11_region_id"].apply(lambda x: all_size[x])
revdf["REF11_region_size"] = revdf["REF11_region_id"].apply(lambda x: all_size[x])

df["aln_by_size"] = df["aln_length"] / df["REF11_region_size"]
revdf["aln_by_size"] = revdf["aln_length"] / revdf["REF11_region_size"]
revdf = revdf[revdf["aln_by_size"] > 0.8]

df = df[df["aln_by_size"] > 0.8]
genome_names = genomes["Genome Abbreviation"].iloc[1::].tolist()
all = all.drop(columns=["REF11_strand"])

newdf = pd.DataFrame(columns=all.columns)

for index, row in all.iterrows():
    newrow = [
        row["region_id"],
        row["REF11_chromosome_name"],
        row["REF11_region_id"],
        row["REF11_chromosome"],
        row["REF11_region_start"],
        row["REF11_region_end"],
        row["REF11_region_size"],
        row["REF11_status"],
    ]
    tempdf = df[df["REF11_region_id"] == row["REF11_region_id"]]
    for genome in genome_names:
        if row[f"{genome}_region_id"] != row[f"{genome}_region_id"]:
            if len(tempdf[tempdf["genome"] == genome]) > 0:
                genome_df = tempdf[tempdf["genome"] == genome].iloc[0]
                newid = f"{genome_df['other_chromosome']}:{genome_df['other_start_forward']}-{genome_df['other_end_forward']}"
                for i in [
                    newid,
                    genome_df["other_chromosome"],
                    genome_df["other_start_forward"],
                    genome_df["other_end_forward"],
                    "COMPATIBLE_SEQUENCE|COMPATIBLE_SEQUENCE",
                    (genome_df["other_end"] - genome_df["other_start"]),
                ]:
                    newrow.append(i)
            else:
                for i in range(0, 6):
                    newrow.append(np.nan)
        else:
            for i in [
                row[f"{genome}_region_id"],
                row[f"{genome}_chromosome"],
                row[f"{genome}_region_start"],
                row[f"{genome}_region_end"],
                row[f"{genome}_status"],
                row[f"{genome}_region_size"],
            ]:
                newrow.append(i)
    newdf.loc[len(newdf)] = newrow

regions = all_full["region_id"]

all_full = all_full[~all_full["REF11_region_id"].str.contains(":")]

newdf = pd.concat([all_full, newdf])

newdf = newdf.set_index("region_id")
newdf = newdf.loc[regions]

newdf = newdf.reset_index()

newdf.to_csv(snakemake.output.df)
