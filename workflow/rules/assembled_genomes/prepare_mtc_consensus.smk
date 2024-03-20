rule concat_sus_scrofa_mtc_genomes:
    input:
        expand(
            "resources/mtc_genomes/{species}/{mtc_genome}.fasta",
            mtc_genome=mtc_genomes[mtc_genomes["Species"] == "Sus_scrofa"][
            "NCBI Accession ID"
            ],
            species="Sus_scrofa",
        ),
    output:
        "data/mtc_genomes/Sus_scrofa/mtc_genomes_concatenated.fasta",
    shell:
        "cat {input} > {output}"


rule align_sus_scrofa_mtc_genomes:
    input:
        "data/mtc_genomes/Sus_scrofa/mtc_genomes_concatenated.fasta",
    output:
        "data/mtc_genomes/Sus_scrofa/mtc_genomes_aligned.fasta",
    threads: 4
    shell:
        " mafft --maxiterate 1000 --globalpair {input} > {output}"


rule consensus_sus_scrofa_mtc_genomes:
    input:
        "data/mtc_genomes/Sus_scrofa/mtc_genomes_aligned.fasta",
    output:
        "data/mtc_genomes/Sus_scrofa/mtc_genome.fasta",
    shell:
        "cons {input} {output}"


rule circularize_consensus_mtc_genome:
    input:
        "data/mtc_genomes/Sus_scrofa/mtc_genome.fasta",
    output:
        "data/mtc_genomes/Sus_scrofa/mtc_genome_circularized.fasta",
    shell:
        "cat {input} > {output} && grep -v '>' {input} >> {output}"


rule reverse_consensus_mtc_genome:
    input:
        "data/mtc_genomes/Sus_scrofa/mtc_genome.fasta",
    output:
        "data/mtc_genomes/Sus_scrofa/mtc_genome_reversed.fasta",
    params:
        reversed_header="reversed_consensus",
    script:
        "../../scripts/assembled_genomes/prepare_mtc_consensus/reverse_fasta_sequence.py"


rule prepare_non_sus_scrofa_mtc_genomes:
    input:
        expand(
            "resources/mtc_genomes/{species}/{mtc_genome}.fasta",
            zip,
            mtc_genome=mtc_genomes[mtc_genomes["Species"] != "Sus_scrofa"][
            "NCBI Accession ID"
            ],
            species=mtc_genomes[mtc_genomes["Species"] != "Sus_scrofa"]["Species"],
        ),
    output:
        expand(
            "data/mtc_genomes/{species}/mtc_genome.fasta",
            zip,
            species=mtc_genomes[mtc_genomes["Species"] != "Sus_scrofa"]["Species"],
        ),
    run:
        for i, o in zip(input, output):
            shutil.copyfile(i, o)


rule circularize_non_scrofa_mtc_genomes:
    input:
        mtgenomes=expand(
            "resources/mtc_genomes/{species}/{mtc_genome}.fasta",
            zip,
            mtc_genome=mtc_genomes[mtc_genomes["Species"] != "Sus_scrofa"][
            "NCBI Accession ID"
            ],
            species=mtc_genomes[mtc_genomes["Species"] != "Sus_scrofa"]["Species"],
        ),
    output:
        circularized_mtgenomes=expand(
            "data/mtc_genomes/{species}/mtc_genome_circularized.fasta",
            zip,
            species=mtc_genomes[mtc_genomes["Species"] != "Sus_scrofa"]["Species"],
        ),
    run:
        for i, o in zip(input.mtgenomes, output.circularized_mtgenomes):
            with open(o, "w") as f:
                lines = []
                for line in open(i, "r").readlines():
                    f.write(line)
                    if line.startswith(">"):
                        pass
                    else:
                        lines.append(line)
                for line in lines:
                    f.write(line)
