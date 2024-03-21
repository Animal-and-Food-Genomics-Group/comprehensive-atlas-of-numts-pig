import pandas as pd
from Bio import AlignIO
from Bio import SeqIO
import os


def split_multifasta(multifasta, outdir):
    """
    This function splits a multifasta in its individual fasta files
    and saves each fasta in a separate file under the same directory
    """
    with open(multifasta, "r") as f:
        mf = SeqIO.parse(f, "fasta")
        for i in mf:
            fasta_id = i.id
            seq = i.seq
            with open(os.path.join(outdir, fasta_id + ".fasta"), "w") as output:
                output.write(">" + str(id) + "\n" + str(seq) + "\n")
                output.close()
    return ()


def maf_to_csv(maf_file, score_threshold=None, csv_outpath=None):
    """
    This function takes as input a .maf file
    and returns a pandas DataFrame containing all of the file's information
    """
    maf_dict = {
        "score": [],
        "EG2": [],
        "E": [],
        "chromosome": [],
        "nuc_start": [],
        "nuc_strand": [],
        "nuc_srcsize": [],
        "nuc_size": [],
        "mtc_start": [],
        "mtc_strand": [],
        "mtc_srcsize": [],
        "mtc_size": [],
        "nuc_seq": [],
        "mtc_seq": [],
    }

    with open(maf_file) as f:
        # parse the maf file to get score, EG2 and E-value
        lines = f.readlines()
        for line in lines:
            if line.startswith("a score="):
                line = line.strip().split()
                line = line[1:]
                maf_dict["score"].append(line[0].split("=")[1])
                maf_dict["EG2"].append(line[1].split("=")[1])
                maf_dict["E"].append(line[2].split("=")[1])

    for multiple_alignment in AlignIO.parse(maf_file, "maf"):
        maf_dict["chromosome"].append(multiple_alignment[0].id)
        maf_dict["nuc_start"].append(multiple_alignment[0].annotations["start"])
        maf_dict["nuc_strand"].append(multiple_alignment[0].annotations["strand"])
        maf_dict["nuc_srcsize"].append(multiple_alignment[0].annotations["srcSize"])
        maf_dict["nuc_size"].append(multiple_alignment[0].annotations["size"])
        maf_dict["nuc_seq"].append(str(multiple_alignment[0].seq))
        maf_dict["mtc_start"].append(multiple_alignment[1].annotations["start"])
        maf_dict["mtc_strand"].append(multiple_alignment[1].annotations["strand"])
        maf_dict["mtc_srcsize"].append(multiple_alignment[1].annotations["srcSize"])
        maf_dict["mtc_size"].append(multiple_alignment[1].annotations["size"])
        maf_dict["mtc_seq"].append(str(multiple_alignment[1].seq))

    maf_dataframe = pd.DataFrame.from_dict(maf_dict)

    if score_threshold is not None:
        maf_dataframe = maf_dataframe.loc[
            maf_dataframe["score"].astype(int) > int(score_threshold)
        ]

    maf_dataframe["mtc_end"] = maf_dataframe["mtc_start"] + maf_dataframe["mtc_size"]
    maf_dataframe["nuc_end"] = maf_dataframe["nuc_start"] + maf_dataframe["nuc_size"]
    if csv_outpath is not None:
        maf_dataframe.to_csv(csv_outpath)

    return maf_dataframe


def merge_maf_and_blasttab(mafpath, btpath, outfile=None):
    """
    This function takes as input a maf and its blasttab conversion
    and merges them in a single pandas DataFrame
    """
    with open(mafpath, "r+") as f:
        if f.read() == "":
            return ()
    maf_df = maf_to_csv(mafpath)
    columns = [
        "query id",
        "subject_id",
        "%_identity",
        "alignment_length",
        "mismatches",
        "gap_opens",
        "query_start",
        "query_end",
        "subject_start",
        "subject_end",
        "evalue",
        "bit_score",
    ]
    bt_df = pd.read_csv(btpath, comment="#", sep="\t", header=None)
    bt_df.columns = columns
    columns_to_update_bt = ["query_start", "subject_start"]
    for col in columns_to_update_bt:
        bt_df[col] = bt_df[col] - 1
    maf_df["score"] = maf_df["score"].astype(int)
    merged_df = pd.merge(maf_df, bt_df, left_index=True, right_index=True)
    if outfile:
        merged_df.to_csv(outfile, index=False)
    return merged_df


import seaborn as sns
from labellines import labelLine, labelLines
import matplotlib.pyplot as plt
from natsort import natsort_keygen
import pandas as pd
import numpy as np


def plot_region(region_df, save_path=None, consider_strands=True):
    """
    This function takes as input the pandas dataframe containing
    a NUMT region and the coordinates of its NUMT fragments
    and returns the fragment plot, with sequence identity
    """
    name = region_df.iloc[0]["region_id"]
    f, ax = plt.subplots()
    sns.set_style("ticks")
    sns.set_theme(style="whitegrid")
    if consider_strands == True:
        if len(region_df["mtc_strand"].unique()) > 1:
            for index, row in region_df.iterrows():
                if row["mtc_strand"] == 1:
                    g = sns.lineplot(
                        x=[row["nuc_start"], row["nuc_end"]],
                        y=[row["mtc_start"], row["mtc_end"]],
                        linewidth=3,
                        color="blue",
                    )  # , label=str(row['numt_id']))
                else:
                    g = sns.lineplot(
                        x=[row["nuc_start"], row["nuc_end"]],
                        y=[row["mtc_end"], row["mtc_start"]],
                        linewidth=3,
                        color="blue",
                    )  # , label=str(row['numt_id']))
        else:
            for index, row in region_df.iterrows():
                g = sns.lineplot(
                    x=[row["nuc_start"], row["nuc_end"]],
                    y=[row["mtc_start"], row["mtc_end"]],
                    linewidth=3,
                    color="blue",
                )  # , label=str(row['numt_id']))
    else:
        for index, row in region_df.iterrows():
            g = sns.lineplot(
                x=[row["nuc_start"], row["nuc_end"]],
                y=[row["mtc_start"], row["mtc_end"]],
                linewidth=3,
                color="blue",
            )  # , label=str(row['numt_id']))
    # labelLines(plt.gca().get_lines(), zorder=2.5, fontsize=10)
    plt.legend([], [], frameon=False)
    region_start = int(region_df.iloc[0]["region_start"])
    region_end = int(region_df.iloc[0]["region_end"])
    plot_xticks = []
    i = region_start
    while i < region_end:
        plot_xticks.append(i)
        i += 500
    plot_xticks.append(region_end)
    x_labels = []
    for i in plot_xticks:
        x_labels.append(i / 1000000)
    ax.set_xticks(plot_xticks, labels=x_labels)
    ax.set_xticklabels(labels=x_labels)  # , rotation=90)#, size=10)
    for label in ax.get_xticklabels()[1:-1]:
        label.set_visible(False)
    ax.set_xlabel("Nuclear coordinates (Mbp)", fontsize=15)
    ax.set_ylabel("Mitochondrial coordinates (bp)", fontsize=15)
    ax.set_title(name, fontsize=15)
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
    return ax


def plot_region_with_identity(region_df, save_path=None, consider_strands=True):
    """
    This function takes as input the pandas dataframe containing
    a NUMT region and the coordinates of its NUMT fragments
    and returns the fragment plot, with sequence identity
    """

    name = region_df.iloc[0]["region_id"]
    f, ax = plt.subplots()
    sns.set_style("ticks")
    sns.set_theme(style="whitegrid")
    if consider_strands == True:
        if len(region_df["mtc_strand"].unique()) > 1:
            for index, row in region_df.iterrows():
                if row["mtc_strand"] == 1:
                    g = sns.lineplot(
                        x=[row["nuc_start"], row["nuc_end"]],
                        y=[row["mtc_start"], row["mtc_end"]],
                        linewidth=3,
                        color="blue",
                        label=str(round(row["%_identity"], 0)),
                    )
                else:
                    g = sns.lineplot(
                        x=[row["nuc_start"], row["nuc_end"]],
                        y=[row["mtc_end"], row["mtc_start"]],
                        linewidth=3,
                        color="blue",
                        label=str(round(row["%_identity"], 0)),
                    )
        else:
            for index, row in region_df.iterrows():
                g = sns.lineplot(
                    x=[row["nuc_start"], row["nuc_end"]],
                    y=[row["mtc_start"], row["mtc_end"]],
                    linewidth=3,
                    color="blue",
                    label=str(round(row["%_identity"], 0)),
                )
    else:
        for index, row in region_df.iterrows():
            g = sns.lineplot(
                x=[row["nuc_start"], row["nuc_end"]],
                y=[row["mtc_start"], row["mtc_end"]],
                linewidth=3,
                color="blue",
                label=str(round(row["%_identity"], 0)),
            )
    labelLines(plt.gca().get_lines(), zorder=2.5, fontsize=5)
    plt.legend([], [], frameon=False)
    region_start = int(region_df.iloc[0]["region_start"])
    region_end = int(region_df.iloc[0]["region_end"])
    plot_xticks = []
    i = region_start
    while i < region_end:
        plot_xticks.append(i)
        i += 500
    plot_xticks.append(region_end)
    x_labels = []
    for i in plot_xticks:
        x_labels.append(i / 1000000)
    ax.set_xticks(plot_xticks, labels=x_labels)
    ax.set_xticklabels(labels=x_labels)  # , rotation=90)#, size=10)
    for label in ax.get_xticklabels()[1:-1]:
        label.set_visible(False)
    ax.set_xlabel("Nuclear coordinates (Mbp)", fontsize=15)
    ax.set_ylabel("Mitochondrial coordinates (bp)", fontsize=15)
    ax.set_title("Region: " + name, fontsize=15)
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
    return ax


def assign_range_and_plot(
    df,
    range,
    outcol,
    column="nuc_size",
    outpath=None,
    x_axis_label=None,
    y_axis_label=None,
):
    """
    This function splits a NUMT df in size ranges and plots the distribution
    """

    def get_interval(val, intervals):
        if val == np.nan:
            return None
        for interval in intervals:
            if val in interval:
                return str(interval.left) + "-" + str(interval.right)
            else:
                continue
        return None

    plt.figure()
    sns.set_theme()
    sns.set_style()
    sns.set(rc={"figure.figsize": (11.7, 8.27)})

    count_by_range = df[column].groupby(pd.cut(df[column], range)).count()
    intervals = count_by_range.index
    df[outcol] = df[column].apply(lambda x: get_interval(x, intervals))
    df[outcol] = df[outcol]
    rangeorder = df[outcol].sort_values(key=natsort_keygen()).unique()
    ax = sns.countplot(
        data=df, x=outcol, order=rangeorder, color=sns.color_palette()[0]
    )
    ax.bar_label(ax.containers[0])
    plt.xticks(rotation=90)
    if x_axis_label:
        ax.set(xlabel=x_axis_label)
    if y_axis_label:
        ax.set(ylabel=y_axis_label)
    if outpath:
        plt.savefig(outpath, dpi=300, bbox_inches="tight")
    return df


def compute_age_dayama(input_msa, divergence_age):
    """
    This function computes the NUMT estimated age using the method by Dayama et al.
    starting from a multiple sequence alignment and the age upper boundary
    """

    if not os.path.isfile(input_msa):
        return (np.nan, np.nan, np.nan)
    msa = SeqIO.parse(open(input_msa), "fasta")
    msa_dict = {}
    for i in msa:
        tempid = i.id.replace("_R_", "")
        msa_dict[tempid] = i.seq
    if len(msa_dict) < 3:
        return (np.nan, np.nan, np.nan)
    # print(msa_dict)
    modern = msa_dict["modern_sequence"]
    ancestral = msa_dict["ancestral_sequence"]
    numt = msa_dict["numt_sequence"]

    avm = 0
    nvm = 0
    pos_counter = -1
    for i, j, z in zip(modern, ancestral, numt):
        pos_counter += 1
        mod_base = i.upper()
        anc_base = j.upper()
        numt_base = z.upper()
        if "-" in mod_base + anc_base:
            continue
        else:
            if mod_base != anc_base:
                avm += 1
                if numt_base == mod_base:
                    nvm += 1
    if avm == 0:
        return (-1, avm, nvm)
    else:
        amr = nvm / avm
        age = (1 - amr) * float(divergence_age)
        return (age, avm, nvm)
