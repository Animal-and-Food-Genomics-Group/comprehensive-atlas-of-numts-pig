# comprehensive-atlas-of-numts-pig
[![Snakemake](https://img.shields.io/badge/snakemake-≥7.19.0-brightgreen.svg)](https://snakemake.github.io)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Code style: snakefmt](https://img.shields.io/badge/code%20style-snakefmt-000000.svg)](https://github.com/snakemake/snakefmt)

Repository including workflow and scripts for the manuscript  
"A comprehensive atlas of nuclear sequences of mitochondrial origin (NUMT) inserted into the pig genome"

For information reach out to [Matteo Bolner](https://github.com/matteobolner) at matteo-dot-bolner2-at-unibo-dot-it

## Workflow structure
This workflow was built using Snakemake (https://snakemake.github.io/).
The structure of the workflow is the following:

├── config/  
├── envs/
└── workflow/  
  ├── rules/  
  │   ├── assembled_genomes/  
  │   ├── common.smk  
  │   ├── general.smk  
  │   └── phylogeny/  
  ├── scripts/  
  │   ├── assembled_genomes/  
  │   │   ├── additional_reference_numts_in_other_genomes/  
  │   │   ├── numt_assembly_comparison/  
  │   │   ├── numt_discovery/  
  │   │   └── prepare_mtc_consensus/  
  │   ├── functions.py  
  │   └── general/  
  └── Snakefile  

- config contains information about the assembled nuclear and mitochondrial genomes used in this work
- envs contains the full list of packages installed in the `mamba` environment used for this project
- Snakefile orchestrates the various rules and specifies the final files to obtain
- rules contains all files specifying rules, divided by topic for clarity; each rule has inputs and outputs
- scripts contains all scripts called by the different rules

### NUMT discovery

A summary of all rules called in the NUMT discovery pipeline is shown below:

[numt_discovery_summary]: https://github.com/matteobolner/comprehensive-atlas-of-numts-pig/blob/main/pipeline_graphs/numt_discovery_summary.svg "Test"
![alt text][numt_discovery_summary]

Briefly:
- Mitochondrial genomes for S. scrofa are concatenated, aligned and used to build a consensus sequence
- The consensus sequence for S. scrofa,  or reference sequence for other species, is circularized and aligned to the nuclear genome (see config for the full list)
- Alignments out of the start/end boundary of the linear sequence are replaced with the complete alignments from the circular sequence
- Duplicate alignments are removed
- NUMT regions are defined from clustered alignments within 20kbp of each other
- Alignments are repeated for each nuclear genome

- Tables, statistics and figures are produced from all the resulting datasets

The full Direct Acyclyc Graph (DAG) for this part of the pipeline can be found in the pipeline_graphs folder ([link](pipeline_graphs/numt_discovery_complete.svg))

### NUMT comparison across assembled genomes

A summary of all rules called in the pipeline for NUMT comparison across genomes  is shown below (the figure includes the previously described NUMT discovery pipeline from which this pipeline gets the input):

[genome_comparison_summary]: https://github.com/matteobolner/comprehensive-atlas-of-numts-pig/blob/main/pipeline_graphs/genome_comparison_summary.png "Test"
![alt text][genome_comparison_summary]

Briefly:
- Each genome's NUMT region coordinates are extracted in a BED file  

- The portions of whole genome alignment overlapping the NUMT coordinates are extracted from the alignment of Sscrofa11.1 and the other assembly  

- The same is repeated for the reverse alignment (other assembly vs Sscrofa11.1)

- The two comparisons are merged to obtain a unique set of alignments

- All genome comparisons are merged in a single table of non-redundant NUMT regions

- NUMT regions that were missed in the NUMT discovery pipeline but align to a NUMT region in Sscrofa11.1 or other genomes are annotated

- The original NUMT discovery datasets are updated with these newly added NUMT regions along with tables, statistics and figures

- An orthologous matrix is produced to summarize all comparisons


The full Direct Acyclyc Graph (DAG) for this part of the pipeline can be found in the pipeline_graphs folder ([link](pipeline_graphs/genome_comparison_complete.svg))
