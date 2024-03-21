import pandas as pd
import os
import shutil


configfile: "config/config.yaml"


genomes = pd.read_table(config["genomes"])
mtc_genomes = pd.read_table(config["mtc_genomes"])
mtc_genomes_for_dayama = pd.read_table(config["mtc_genomes_for_dayama"], header=None)[
    0
].tolist()


wildcard_constraints:
    species="([a-zA-Z_]+)",
    other_species="([a-zA-Z_]+)",
    main_species="([a-zA-Z_]+)",
    assembly_accession="(GC[AF]_[0-9]{9}(\.[0-9])?)",
    other_assembly_accession="(GC[AF]_[0-9]{9}(\.[0-9])?)",
    main_assembly_accession="(GC[AF]_[0-9]{9}(\.[0-9])?)",


try:
    all_regions, all_numts = glob_wildcards(
        "data/phylogeny/sequences/all_sequences/{region}/{numt}.fasta"
    )
except:
    all_regions, all_numts = (None, None)

groups = ["sus", "suinae", "phacocerus", "cebifrons"]


def get_chromosomes_from_report(species, name, accession):
    tempdf = pd.read_table(
        f"resources/genomes/{species}/{name}/ncbi/{accession}_{name}_assembly_report.tsv"
    )
    chroms = tempdf["Sequence-ID"].tolist()
    return chroms


def get_chromosome_and_start_end_from_report(species, name, accession):
    tempdf = pd.read_table(
        f"resources/genomes/{species}/{name}/ncbi/{accession}_{name}_assembly_report.tsv"
    )
    tempdict = {
        chrom: end
        for chrom, end in zip(tempdf["Sequence-ID"], tempdf["Sequence-Length"])
    }
    return tempdict


def get_end_from_chromosome(species, name, accession, chromosome):
    tempdict = get_chromosome_and_start_end_from_report(species, name, accession)
    end = tempdict[chromosome]
    return end
