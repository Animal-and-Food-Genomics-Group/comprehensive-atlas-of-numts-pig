import pandas as pd

# df=pd.read_csv("results/genome_comparison/all_numts.csv")
df = pd.read_csv(snakemake.input.regions)

# genome_info=pd.read_table("config/genomes.tsv")
genome_info = pd.read_table(snakemake.input.genomes)

# i=row['Genome Abbreviation']
# species=row['Species']
# ass_name=row['Assembly Name']
# ass_acc=row['Assembly Accession']

i = snakemake.params.genome
species = snakemake.params.species
ass_name = snakemake.params.assembly_name
ass_acc = snakemake.params.assembly_accession


outdf = df[df[f"{i}_region_id"].fillna("A").str.contains(":")][
    [f"{i}_chromosome", f"{i}_region_start", f"{i}_region_end", f"{i}_region_id"]
]
outdf[f"{i}_region_start"] = outdf[f"{i}_region_start"].astype(int)
outdf[f"{i}_region_end"] = outdf[f"{i}_region_end"].astype(int)
outdf.to_csv(snakemake.output[0], index=False, header=None, sep="\t")
