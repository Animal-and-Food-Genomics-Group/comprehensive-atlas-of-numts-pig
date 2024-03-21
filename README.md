# comprehensive-atlas-of-numts-pig
[![Snakemake](https://img.shields.io/badge/snakemake-≥7.19.0-brightgreen.svg)](https://snakemake.github.io)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Code style: snakefmt](https://img.shields.io/badge/code%20style-snakefmt-000000.svg)](https://github.com/snakemake/snakefmt)

Repository including workflow and scripts for the manuscript  
"A comprehensive atlas of nuclear sequences of mitochondrial origin (NUMT) inserted into the pig genome"

For information reach out to [Matteo Bolner](https://github.com/matteobolner) at matteo.bolner2@unibo.it

## Workflow description
This workflow was built using Snakemake (https://snakemake.github.io/).
The structure of the workflow is the following:

├── config/  
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
- Snakefile orchestrates the various rules and specifies the final files to obtain
- rules contains all files specifying rules, divided by topic for clarity; each rule has inputs and outputs
- scripts contains all scripts called by the different rules

A summary of all rules called in the NUMT discovery pipeline is shown below:

[numt_discovery_summary]: https://github.com/matteobolner/comprehensive-atlas-of-numts-pig/blob/main/pipeline_graphs/numt_discovery_summary.svg "Test"
![alt text][numt_discovery_summary]

The full Direct Acyclyc Graph (DAG) for this part of the pipeline can be found in the pipeline_graphs folder ([link](pipeline_graphs/numt_discovery_complete.svg))


A summary of all rules called in the pipeline for NUMT comparison across genomes  is shown below:

[genome_comparison_summary]: https://github.com/matteobolner/comprehensive-atlas-of-numts-pig/blob/main/pipeline_graphs/genome_comparison_summary.png "Test"
![alt text][genome_comparison_summary]

The full Direct Acyclyc Graph (DAG) for this part of the pipeline can be found in the pipeline_graphs folder ([link](pipeline_graphs/genome_comparison_complete.svg))
