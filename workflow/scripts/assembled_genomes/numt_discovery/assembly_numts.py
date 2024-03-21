import pandas as pd
import sys
import os
import numpy as np
from natsort import natsort_keygen
from operator import attrgetter
from pybedtools import BedTool
from pathlib import Path

sys.path.append(snakemake.params.functions_path)
import functions as custom_functions

df = pd.read_csv(snakemake.input.input_data)
report = pd.read_table(snakemake.input.report)
mtc_ids = pd.read_table(snakemake.input.mtc_scaffolds)
species = snakemake.params.species

# df=pd.read_csv("data/numt_discovery/lastal/linear_and_circular_merged/Sus_scrofa/Sscrofa11.1/GCF_000003025.6_Sscrofa11.1/linear_and_circular_merged.csv")
# report=pd.read_table("resources/genomes/Sus_scrofa/Sscrofa11.1/ncbi/GCF_000003025.6_Sscrofa11.1_assembly_report.tsv", comment='#')
# mtc_ids=pd.read_table("stats/general/mtc_scaffolds.tsv")
mtc_ids_list = mtc_ids["Sequence-ID"].tolist()

# species='Sus_scrofa'
species_specific_chromosome_name = (species[0] + species.split("_")[1][0]).upper() + "C"

df = df[~df["chromosome"].isin(mtc_ids_list)].copy()

chromosome_range = [str(i) for i in range(1, 50)] + ["X", "Y"]

renaming_dict = {}
scaffold_or_chromosome_dict = {}
for i, j in zip(report["Sequence-ID"], report["Chromosome/Scaffold Name"]):
    if j in chromosome_range:
        j = species_specific_chromosome_name + j
        scaffold_or_chromosome_dict[i] = "chromosome"
    else:
        scaffold_or_chromosome_dict[i] = "scaffold"
    renaming_dict[i] = j

df["chromosome_code"] = df["chromosome"]
df["chromosome"] = df["chromosome"].apply(lambda x: renaming_dict[x])
df["chromosome_type"] = df["chromosome_code"].apply(
    lambda x: scaffold_or_chromosome_dict[x]
)

df = df.sort_values(
    by=["chromosome_type", "chromosome", "nuc_start"], key=natsort_keygen()
).reset_index(drop=True)

numt_ids_dict = {}

for name, group in df.groupby(by="chromosome"):
    counter = 0
    for index, row in group.iterrows():
        counter += 1
        numt_ids_dict[index] = f"{row['chromosome']}_NUMT_{counter}"

df["numt_id"] = df.reset_index()["index"].apply(lambda x: numt_ids_dict[x])

df_bed = df[["chromosome", "nuc_start", "nuc_end", "numt_id"]]

regions_list = []
for name, group in df_bed.groupby(by="chromosome"):
    tempbed = BedTool.from_dataframe(group)
    regions_bed = tempbed.cluster(d=20000).to_dataframe()
    regions_bed.columns = ["chrom", "start", "end", "numt_id", "region_id"]
    regions_bed["region_id"] = (
        regions_bed["chrom"]
        + "_"
        + "REGION"
        + "_"
        + regions_bed["region_id"].astype(str)
    )
    regions_list.append(regions_bed)

region_df = pd.concat(regions_list)

region_dict = {i: j for i, j in zip(region_df["numt_id"], region_df["region_id"])}
BedTool.from_dataframe(group).cluster(d=20000)

df["region_id"] = df["numt_id"].apply(lambda x: region_dict[x])


# assign distance between alignments in nuc and mtc

distance_dict_nuc = {}
distance_dict_mtc = {}

for name, group in df.groupby(by="chromosome"):
    group["distance_to_next_numt_nuc"] = group["nuc_start"].shift(-1) - group["nuc_end"]
    group["distance_to_next_numt_mtc"] = group["mtc_start"].shift(-1) - group["mtc_end"]

    for index, row in group.iterrows():
        distance_dict_nuc[row["numt_id"]] = row["distance_to_next_numt_nuc"]
        distance_dict_mtc[row["numt_id"]] = row["distance_to_next_numt_mtc"]

df["distance_to_next_numt_nuc"] = df["numt_id"].apply(lambda x: distance_dict_nuc[x])
df["distance_to_next_numt_mtc"] = df["numt_id"].apply(lambda x: distance_dict_mtc[x])

region_starts = {}
region_ends = {}

for name, group in df.groupby(by="region_id"):
    region_starts[name] = group["nuc_start"].min()
    region_ends[name] = group["nuc_end"].max()

df["region_start"] = df["region_id"].apply(lambda x: region_starts[x])
df["region_end"] = df["region_id"].apply(lambda x: region_ends[x])

df["region_length"] = df["region_end"] - df["region_start"]


# REMOVE ALL NUMT REGIONS WITH START OR END WITHIN 1000BP OF SCAFFOLD/CHROMOSOME EXTREMITY -> flanking shorter than 1000

numts_too_close_to_extremity = df[
    (df["region_start"] < 1000) | ((df["nuc_srcsize"] - df["region_end"]) < 1000)
]
numts_too_close_to_extremity.to_csv(
    snakemake.output.dropped_putative_numts_close_to_extremity
)

df = df[df["region_start"] >= 1000]
df = df[(df["nuc_srcsize"] - df["region_end"]) >= 1000]


# select numts repeating/overlapping with each other and drop them from main df, keeping only one

overlapping_numts = []

for name, group in df.groupby(by="region_id"):
    for index, row in group.iterrows():
        if row["distance_to_next_numt_nuc"] <= 0:
            overlapping_numts.append(index)
            overlapping_numts.append(index + 1)

overlapping_numts_df = df.loc[overlapping_numts].drop_duplicates()
overlapping_numts_df = overlapping_numts_df.sort_index()

significantly_overlapping_regions = df[
    df["region_id"].isin(
        overlapping_numts_df[
            overlapping_numts_df["distance_to_next_numt_nuc"]
            / overlapping_numts_df["nuc_size"]
            <= -0.5
        ]["region_id"]
    )
]

all_numts_to_keep = []
for name, group in significantly_overlapping_regions.groupby(by="region_id"):
    numts_to_keep = set()
    overlapped_numts = set()
    for index, row in group.iterrows():
        if row["distance_to_next_numt_nuc"] / row["nuc_size"] <= -0.5:
            overlapped_numts.add(row["numt_id"])
            overlapped_numts.add(df.loc[index + 1]["numt_id"])
        else:
            numts_to_keep.add(row["numt_id"])
    numts_to_keep = numts_to_keep - overlapped_numts
    new_group = group[group["numt_id"].isin(list(overlapped_numts))]
    biggestnumt = 0
    numt_to_keep = None
    for new_index, new_row in new_group.iterrows():
        if new_row["nuc_size"] > biggestnumt:
            biggestnumt = new_row["nuc_size"]
            numt_to_keep = new_row["numt_id"]
        else:
            continue
    # this is to check that the one of the overlapped identical numts is not already in the numt before
    check_for_last_overlap = group[group["numt_id"] == numt_to_keep].iloc[0]
    for i, r in group[group["numt_id"].isin(list(numts_to_keep))].iterrows():
        if (
            len(
                set(
                    range(
                        check_for_last_overlap["nuc_start"],
                        check_for_last_overlap["nuc_end"] + 1,
                    )
                ).intersection(set(range(r["nuc_start"], r["nuc_end"] + 1)))
            )
            >= check_for_last_overlap["nuc_size"] / 2
        ):
            numt_to_keep = None
        else:
            continue
    numts_to_keep.add(numt_to_keep)
    all_numts_to_keep.extend(list(numts_to_keep))

Path(snakemake.output.duplicate_numts_plots_directory).mkdir(
    parents=True, exist_ok=True
)

for name, group in significantly_overlapping_regions.groupby(by="region_id"):
    custom_functions.plot_region(
        group,
        save_path=os.path.join(
            snakemake.output.duplicate_numts_plots_directory, name + ".png"
        ),
    )

numts_to_drop_from_overlap_regions_df = significantly_overlapping_regions[
    ~significantly_overlapping_regions["numt_id"].isin(all_numts_to_keep)
]
numts_to_drop_from_overlap_regions_df.to_csv(
    snakemake.output.duplicate_numts, index=False
)
numts_to_drop_from_overlap_regions = numts_to_drop_from_overlap_regions_df["numt_id"]

df = df[~df["numt_id"].isin(numts_to_drop_from_overlap_regions)]


numt_ids_dict = {}
df = df.reset_index(drop=True)
for name, group in df.groupby(by="chromosome"):
    counter = 0
    for index, row in group.iterrows():
        counter += 1
        numt_ids_dict[index] = f"{row['chromosome']}_NUMT_{counter}"

df["numt_id"] = df.reset_index()["index"].apply(lambda x: numt_ids_dict[x])

df.to_csv(snakemake.output.numts_table, index=False)
