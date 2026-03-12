## Phylotranscriptomics allows distinguishing major gene flow events from incomplete lineage sorting in rapidly diversifying mimetic orchids (genus Ophrys)
Script from the paper "Phylotranscriptomics allows distinguishing major gene flow events from incomplete lineage sorting in rapidly diversifying mimetic orchids (genus Ophrys)"

## assembly_pipeline.sh
De Novo transcriptome assembly pipeline from raw paired-end RNA-seq reads. Performs read quality assessment and filtering, transcriptome assembly using multiple tools, redundacy removal and assembly completness evaluation

### Requirements

The pipeline requires the following software to be installed:

- FastQC  
- Trimmomatic
- Trinity
- Trans-Abyss
- rnaSPAdes
- EviGeneR  
- BUSCO

### Usage

Assembly_pipeline.sh sample_R1.fastq.gz sample_R2.fastq.gz

## data_to_twisst.R
Use genomic position and phylogenetic tree to generate the input files required to produce TWISST (https://github.com/simonhmartin/twisst) plots.
Use a modify version of plot_twisst.R (https://github.com/simonhmartin/twisst) called "plot_twisst_modify.R"

### Input

The script requires the following input files:

| Parameter | Description |
|-----------|-------------|
| `--input_pos` | File containing gene alignment genomic window positions |
| `--tree_folder` | Folder containing phylogenetic trees |
| `--tree_pattern` | Pattern identifying tree files |
| `--taxa` | Comma-separated list of taxa used in the TWISST analysis |
| `--prune_tree` | Whether trees should be pruned to the specified taxa (TRUE/FALSE) |
| `--chromo_size` | File containing chromosome or scaffold sizes |
| `--alignment_size` | File containing gene alignment size information |

### Usage

Assembly_pipeline.sh sample_R1.fastq.gz sample_R2.fastq.gz
