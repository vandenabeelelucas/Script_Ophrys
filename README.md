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

### Parameters

| Parameter | Description |
|-----------|-------------|
| `STRANDED` | Indicates whether the RNA-seq reads are strand-specific (`yes`/`no`) |
| `KMER_LIST` | List of k-mer sizes to use for assemblies |
| `THREADS` | Number of CPU threads to use  |
| `MEMORY` | Amount of RAM (GB) to allocate  |
| `R1` | Path to forward reads (first pair)  |
| `R2` | Path to reverse reads (second pair)  |
| `BUSCO_LINEAGE` | BUSCO lineage dataset used for assembly completeness |

### Usage
```
Assembly_pipeline.sh sample_R1.fastq.gz sample_R2.fastq.gz
```
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
| `--taxa` | Comma-separated list of taxa used in the TWISST analysis (H1,H2,H3,outgroup) |
| `--prune_tree` | Whether trees should be pruned to the specified taxa (`TRUE`/`FALSE`) |
| `--chromo_size` | File containing chromosome or scaffold sizes |
| `--alignment_size` | File containing gene alignment size information |

### Usage
```
Rscript data_to_twisst.R \
  --input_pos file.tsv \
  --tree_folder /path/to/trees \
  --tree_pattern treefile \
  --taxa H1,H2,H3,Outgroup \
  --prune_tree TRUE/FALSE \
  --chromo_size file.tsv \
  --alignment_size file.tsv
```
