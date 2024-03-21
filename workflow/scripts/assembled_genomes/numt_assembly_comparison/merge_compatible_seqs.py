import pandas as pd
import sys
from natsort import natsort_keygen
from pybedtools import BedTool
import numpy as np

sys.path.append(snakemake.params.functions_path)
# import workflow.scripts.functions as custom_functions
import functions as custom_functions


mtc_ids = pd.read_table(snakemake.input.mtc_scaffolds)
species = snakemake.params.species

mtc_ids_list = mtc_ids["Sequence-ID"].tolist()
species_specific_chromosome_name = (species[0] + species.split("_")[1][0]).upper() + "C"

report = pd.read_table(snakemake.input.report)

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


maf = snakemake.input.maf
bt = snakemake.input.bt


try:
    df = custom_functions.merge_maf_and_blasttab(maf, bt)
except:
    df = pd.DataFrame()

ref = pd.read_csv(snakemake.input.numts)

if len(df) > 0:

    df["seq_start"] = (
        df["chromosome"].str.split(":").str[1].str.split("-").str[0].astype(int)
    )
    df["seq_end"] = (
        df["chromosome"].str.split(":").str[1].str.split("-").str[1].astype(int)
    )
    df["seq_chr"] = df["chromosome"].str.split(":").str[0]

    df = df.rename(columns={"chromosome": "additional_region_id"})

    df["nuc_start"] = df["seq_start"] + df["nuc_start"]
    df["nuc_end"] = df["seq_start"] + df["nuc_end"]
    df["nuc_size"] = df["nuc_end"] - df["nuc_start"]

    df = df.drop(columns=["seq_start", "seq_end"])
    df = df.rename(columns={"seq_chr": "chromosome_code"})

    df["subject_start"] = df["nuc_start"]
    df["subject_end"] = df["nuc_end"]

    df["chromosome"] = df["chromosome_code"].apply(lambda x: renaming_dict[x])

    df["chromosome_type"] = df["chromosome_code"].apply(
        lambda x: scaffold_or_chromosome_dict[x]
    )

    df = df.sort_values(
        by=["chromosome_type", "chromosome", "nuc_start"], key=natsort_keygen()
    ).reset_index(drop=True)

    chromsizes = {
        i: j for i, j in zip(report["Sequence-ID"], report["Sequence-Length"])
    }
    df["nuc_srcsize"] = df["chromosome_code"].apply(lambda x: chromsizes[x])

    max_numt_dict = {i: 0 for i in df["chromosome"].unique()}

    for name, group in ref.groupby(by="chromosome"):
        max_numt_dict[name] = group["numt_id"].str.split("_").str[-1].astype(int).max()

    numt_ids_dict = {}

    for name, group in df.groupby(by="chromosome"):
        counter = max_numt_dict[name]
        for index, row in group.iterrows():
            counter += 1
            numt_ids_dict[index] = f"{row['chromosome']}_NUMT_{counter}"

    df["numt_id"] = df.reset_index()["index"].apply(lambda x: numt_ids_dict[x])

    df_bed = df[["chromosome", "nuc_start", "nuc_end", "numt_id"]]

    max_region_dict = {i: 0 for i in df["chromosome"].unique()}

    for name, group in ref.groupby(by="chromosome"):
        max_region_dict[name] = (
            group["region_id"].str.split("_").str[-1].astype(int).max()
        )

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
            + (regions_bed["region_id"] + max_region_dict[name]).astype(str)
        )
        regions_list.append(regions_bed)

    region_df = pd.concat(regions_list)

    region_dict = {i: j for i, j in zip(region_df["numt_id"], region_df["region_id"])}

    df["region_id"] = df["numt_id"].apply(lambda x: region_dict[x])

    distance_dict_nuc = {}
    distance_dict_mtc = {}

    for name, group in df.groupby(by="chromosome"):
        group["distance_to_next_numt_nuc"] = (
            group["nuc_start"].shift(-1) - group["nuc_end"]
        )
        group["distance_to_next_numt_mtc"] = (
            group["mtc_start"].shift(-1) - group["mtc_end"]
        )

        for index, row in group.iterrows():
            distance_dict_nuc[row["numt_id"]] = row["distance_to_next_numt_nuc"]
            distance_dict_mtc[row["numt_id"]] = row["distance_to_next_numt_mtc"]

    df["distance_to_next_numt_nuc"] = df["numt_id"].apply(
        lambda x: distance_dict_nuc[x]
    )
    df["distance_to_next_numt_mtc"] = df["numt_id"].apply(
        lambda x: distance_dict_mtc[x]
    )

    region_starts = {}
    region_ends = {}

    for name, group in df.groupby(by="region_id"):
        region_starts[name] = group["nuc_start"].min()
        region_ends[name] = group["nuc_end"].max()

    df["region_start"] = df["region_id"].apply(lambda x: region_starts[x])
    df["region_end"] = df["region_id"].apply(lambda x: region_ends[x])

    df["region_length"] = df["region_end"] - df["region_start"]

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

    numts_to_drop_from_overlap_regions_df = significantly_overlapping_regions[
        ~significantly_overlapping_regions["numt_id"].isin(all_numts_to_keep)
    ]
    numts_to_drop_from_overlap_regions = numts_to_drop_from_overlap_regions_df[
        "numt_id"
    ]

    df = df[~df["numt_id"].isin(numts_to_drop_from_overlap_regions)]

    numt_ids_dict = {}
    df = df.reset_index(drop=True)
    for name, group in df.groupby(by="chromosome"):
        counter = max_numt_dict[name]
        for index, row in group.iterrows():
            counter += 1
            numt_ids_dict[index] = f"{row['chromosome']}_NUMT_{counter}"

    df["numt_id"] = df.reset_index()["index"].apply(lambda x: numt_ids_dict[x])

    df = pd.concat([ref, df])
    df = df.reset_index(drop=True)
    df = df.sort_values(
        by=["chromosome_type", "chromosome", "nuc_start"], key=natsort_keygen()
    ).reset_index(drop=True)

    distance_dict_nuc = {}
    distance_dict_mtc = {}

    for name, group in df.groupby(by="chromosome"):
        group["distance_to_next_numt_nuc"] = (
            group["nuc_start"].shift(-1) - group["nuc_end"]
        )
        group["distance_to_next_numt_mtc"] = (
            group["mtc_start"].shift(-1) - group["mtc_end"]
        )

        for index, row in group.iterrows():
            distance_dict_nuc[row["numt_id"]] = row["distance_to_next_numt_nuc"]
            distance_dict_mtc[row["numt_id"]] = row["distance_to_next_numt_mtc"]

    df["distance_to_next_numt_nuc"] = df["numt_id"].apply(
        lambda x: distance_dict_nuc[x]
    )
    df["distance_to_next_numt_mtc"] = df["numt_id"].apply(
        lambda x: distance_dict_mtc[x]
    )

    additional_numts = pd.read_table(snakemake.input.additional_numts_bed, header=None)

    not_found = additional_numts[~additional_numts[3].isin(df["additional_region_id"])]
    not_found.columns = [
        "chromosome_code",
        "nuc_start",
        "nuc_end",
        "additional_region_id",
    ]
    if len(not_found) > 0:

        for i in df.columns:
            if i not in not_found.columns:
                not_found[i] = np.nan

        not_found = not_found[list(df.columns)]

        not_found["chromosome"] = not_found["chromosome_code"].apply(
            lambda x: renaming_dict[x]
        )
        not_found["nuc_strand"] = 1
        not_found["nuc_size"] = not_found["nuc_end"] - not_found["nuc_start"]
        not_found["mtc_size"] = not_found["nuc_size"]
        not_found["mtc_srcsize"] = df["mtc_srcsize"].min()
        not_found["chromosome_type"] = not_found["chromosome_code"].apply(
            lambda x: scaffold_or_chromosome_dict[x]
        )

        not_found = not_found.reset_index(drop=True)

        max_numt_dict = {i: 0 for i in not_found["chromosome"].unique()}

        for name, group in df.groupby(by="chromosome"):
            max_numt_dict[name] = (
                group["numt_id"].str.split("_").str[-1].astype(int).max()
            )

        numt_ids_dict = {}

        for name, group in not_found.groupby(by="chromosome"):
            counter = max_numt_dict[name]
            for index, row in group.iterrows():
                counter += 1
                numt_ids_dict[index] = f"{row['chromosome']}_NUMT_{counter}"

        not_found["numt_id"] = not_found.reset_index()["index"].apply(
            lambda x: numt_ids_dict[x]
        )

        df_bed = not_found[["chromosome", "nuc_start", "nuc_end", "numt_id"]]
        df_bed = df_bed.sort_values(by=["chromosome", "nuc_start"])

        max_region_dict = {i: 0 for i in df["chromosome"].unique()}

        for name, group in df.groupby(by="chromosome"):
            max_region_dict[name] = (
                group["region_id"].str.split("_").str[-1].astype(int).max()
            )
        for i in not_found["chromosome"].unique():
            max_region_dict.setdefault(i, 0)
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
                + (regions_bed["region_id"] + max_region_dict[name]).astype(str)
            )
            regions_list.append(regions_bed)
        region_df = pd.concat(regions_list)
        region_dict = {
            i: j for i, j in zip(region_df["numt_id"], region_df["region_id"])
        }

        not_found["region_id"] = not_found["numt_id"].apply(lambda x: region_dict[x])

        df = pd.concat([df, not_found])

    df.to_csv(snakemake.output[0], index=False)

else:
    ref.to_csv(snakemake.output[0], index=False)
