import pandas as pd
from Bio import AlignIO
import numpy as np

main = "REF11"


def maf_to_csv(maf_file, score_threshold=None, csv_outpath=None):
    # this function takes as input a .maf file and returns a pandas DataFrame containing all of the file's information
    # maf_file = os.path.abspath(maf_file)
    maf_dict = {
        "score": [],
        "main_chromosome": [],
        "other_chromosome": [],
        "mismap": [],
        "main_start": [],
        "main_strand": [],
        "main_srcsize": [],
        "main_size": [],
        "other_start": [],
        "other_strand": [],
        "other_srcsize": [],
        "other_size": [],
        "main_seq": [],
        "other_seq": [],
    }

    with open(maf_file) as f:
        # parse the maf file to get score, EG2 and E-value
        lines = f.readlines()
        for line in lines:
            if line.startswith("a score="):
                line = line.strip().split()
                line = line[1:]
                maf_dict["score"].append(line[0].split("=")[1])
                maf_dict["mismap"].append(line[1].split("=")[1])

    for multiple_alignment in AlignIO.parse(maf_file, "maf"):
        maf_dict["main_chromosome"].append(multiple_alignment[0].id)
        maf_dict["other_chromosome"].append(multiple_alignment[1].id)

        maf_dict["main_start"].append(multiple_alignment[0].annotations["start"])
        maf_dict["main_strand"].append(multiple_alignment[0].annotations["strand"])
        maf_dict["main_srcsize"].append(multiple_alignment[0].annotations["srcSize"])
        maf_dict["main_size"].append(multiple_alignment[0].annotations["size"])
        maf_dict["main_seq"].append(str(multiple_alignment[0].seq))
        maf_dict["other_start"].append(multiple_alignment[1].annotations["start"])
        maf_dict["other_strand"].append(multiple_alignment[1].annotations["strand"])
        maf_dict["other_srcsize"].append(multiple_alignment[1].annotations["srcSize"])
        maf_dict["other_size"].append(multiple_alignment[1].annotations["size"])
        maf_dict["other_seq"].append(str(multiple_alignment[1].seq))

    maf_dataframe = pd.DataFrame.from_dict(maf_dict)

    if score_threshold is not None:
        maf_dataframe = maf_dataframe.loc[
            maf_dataframe["score"].astype(int) > int(score_threshold)
        ]

    maf_dataframe["other_end"] = (
        maf_dataframe["other_start"] + maf_dataframe["other_size"]
    )
    maf_dataframe["main_end"] = maf_dataframe["main_start"] + maf_dataframe["main_size"]
    maf_dataframe["other_start_forward"] = maf_dataframe["other_start"].where(
        maf_dataframe["other_strand"] == 1,
        (maf_dataframe["other_srcsize"] - maf_dataframe["other_end"]) + 1,
    )
    maf_dataframe["other_end_forward"] = maf_dataframe["other_end"].where(
        maf_dataframe["other_strand"] == 1,
        (maf_dataframe["other_srcsize"] - maf_dataframe["other_start"]),
    )
    # if private numt perfectly, other start and end are identical, and in negative strand conversion gives end=start-1 basically
    maf_dataframe["other_start_forward"] = maf_dataframe["other_start"].where(
        maf_dataframe["other_start"] == maf_dataframe["other_end"],
        maf_dataframe["other_start_forward"],
    )
    maf_dataframe["other_end_forward"] = maf_dataframe["other_end"].where(
        maf_dataframe["other_start"] == maf_dataframe["other_end"],
        maf_dataframe["other_end_forward"],
    )

    if csv_outpath is not None:
        maf_dataframe.to_csv(csv_outpath)

    return maf_dataframe


# df=maf_to_csv("data/genome_comparison/additional_regions_found_in_reference/regions_extracted_from_maf/Sscrofa11.1/Sus_scrofa/Jinhua_pig_v1/GCA_001700295.1_in_GCF_000003025.6.maf")
df = maf_to_csv(snakemake.input.main_maf)

# df_bt=pd.read_table("data/genome_comparison/additional_regions_found_in_reference/regions_extracted_from_maf/Sscrofa11.1/Sus_scrofa/Jinhua_pig_v1/GCA_001700295.1_in_GCF_000003025.6.blasttab", header=None)
df_bt = pd.read_table(snakemake.input.main_bt, header=None)

df_bt.columns = [
    "other_chromosome",
    "main_chromosome",
    "%_identity",
    "aln_length",
    "mismatches",
    "gap_opens",
    "other_start",
    "other_end",
    "main_start",
    "main_end",
]

df_bt = df_bt[["%_identity", "aln_length", "mismatches", "gap_opens"]]
df = pd.concat([df, df_bt], axis=1)

df["main_strand"] = df["main_strand"].replace(1, "+")
df["main_strand"] = df["main_strand"].replace(-1, "-")
df["other_strand"] = df["other_strand"].replace(1, "+")
df["other_strand"] = df["other_strand"].replace(-1, "-")

# rev_df=maf_to_csv("data/genome_comparison/additional_regions_found_in_reference/regions_extracted_from_maf/Sscrofa11.1/Sus_scrofa/Jinhua_pig_v1/GCF_000003025.6_in_GCA_001700295.1.maf")
rev_df = maf_to_csv(snakemake.input.rev_maf)

# rev_df_bt=pd.read_table("data/genome_comparison/additional_regions_found_in_reference/regions_extracted_from_maf/Sscrofa11.1/Sus_scrofa/Jinhua_pig_v1/GCF_000003025.6_in_GCA_001700295.1.blasttab", header=None)
rev_df_bt = pd.read_table(snakemake.input.rev_bt, header=None)

rev_df_bt.columns = [
    "other_chromosome",
    "main_chromosome",
    "%_identity",
    "aln_length",
    "mismatches",
    "gap_opens",
    "other_start",
    "other_end",
    "main_start",
    "main_end",
]

rev_df_bt = rev_df_bt[["%_identity", "aln_length", "mismatches", "gap_opens"]]


rev_df = pd.concat([rev_df, rev_df_bt], axis=1)

rev_df["main_strand"] = rev_df["main_strand"].replace(1, "+")
rev_df["main_strand"] = rev_df["main_strand"].replace(-1, "-")
rev_df["other_strand"] = rev_df["other_strand"].replace(1, "+")
rev_df["other_strand"] = rev_df["other_strand"].replace(-1, "-")

df.columns = [i.replace("main_", f"{main}_") for i in df.columns]
rev_df.columns = [
    i.replace("other_", f"{main}_").replace("main_", f"other_") for i in rev_df.columns
]

df.to_csv(snakemake.output.df, index=False)
rev_df.to_csv(snakemake.output.rev_df, index=False)
