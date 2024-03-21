from Bio import SeqIO

fasta = snakemake.input[0]

with open(snakemake.output[0], "w") as outfile:
    with open(fasta) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            outfile.write(f">{snakemake.params.reversed_header}")
            outfile.write("\n")
            outfile.write(str(record.seq[::-1]))
