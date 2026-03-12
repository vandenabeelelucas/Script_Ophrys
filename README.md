# Phylotranscriptomics allows distinguishing major gene flow events from incomplete lineage sorting in rapidly diversifying mimetic orchids (genus Ophrys)
Script from the study "Phylotranscriptomics allows distinguishing major gene flow events from incomplete lineage sorting in rapidly diversifying mimetic orchids (genus Ophrys)" ()

## assembly_pipeline.sh
_De Novo_ transcriptome assembly pipeline from raw paired-end RNA-seq reads. Performs read quality assessment and filtering, transcriptome assembly using multiple tools, redundacy removal and assembly completness evaluation

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
Use genomic position and phylogenetic tree to generate the input files required to produce TWISST (https://github.com/simonhmartin/twisst) plots. Uses a modified version of plot_twisst.R called `plot_twisst_modify.R`.

### Parameters

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

## select-aliMCL-taxonomy.c

Filters sequence files (.ali) based on number of species, number of paralogs per species, and number of clades present.

### Parameters

| Parameter | Description |
|-----------|-------------|
| `min_spec` | Minimum number of species |
| `min_clade` | Minimum number of clades (optional) |
| `max_para` | Maximum number of paralogs (optional) |
| `file_clade` | File containing the definition of clades (optional) |
| `file_needed_clade` | file containing the list of clades that must be present (optional) |
| `ali` | name of file .ali (optional) |

### Usage
```bash
select-aliMCL-taxonomy min_spec=10 min_clade=5 file_clade=Example_file_clade.txt
```
## fasta2complementary.c
Reverse-complements sequences on the complementary strand of a FASTA file, using the longest sequence as reference to determine alignment orientation.

### Requirements

The script requires the following software to be installed:

-BLASTn

### Usage
```bash
fasta2complementary
```

## ali2stat.c

Provides various statistics on an alignment file in .ali format.

| Parameter | Description |
|-----------|-------------|
| `type` | protein or DNA |
| `infospecies` | yes or no (default=no) |
| `ali` | name of file .ali (optional) |

### Usage
```bash
ali2stat type=DNA
```
## split-gene2equal-info.c

Splits an alignment file (.ali) into two files, each containing half of the informative positions, optionally considering gaps in the calculation.

| Parameter | Description |
|-----------|-------------|
| `type` | protein or DNA (optional)|
| `gap` | known or unknown (default=unknown) |
| `ali` | name of file .ali (optional)  |

### Usage
```bash
split-gene2equal-info gap=unknown
```
## root-max-div-taxon.c

Detects and separates from an alignment file (.ali) two distinct groups of paralogous sequences by splitting a corresponding phylogenetic tree (.ali) into two subtrees that maximize taxonomic diversity.

| Parameter | Description |
|-----------|-------------|
| `arb` | name of file .arb |
| `ali` | name of file .ali |
| `criterion0` | minimum number of species in clade 1, minimum number of species in clade 2, percent of internal branches longer than the one used to split the tree, minimum number of species shared by the two subtrees= (e.g. `criterion0=20,15,10,0`) |
| `criterion1` | minimum number of species in clade 1, minimum number of species in clade 2, percent of internal branches longer than the one used to split the tree, minimum number of species shared by the two subtrees= (e.g. `criterion1=20,15,10,15`) |

### Usage
```bash
root-max-div-taxon ali=file.ali arb=corresponding_tree.arb criterion0=20,15,10,0 criterion1=20,15,10,0
```
## detect-problems-arb.c

Check orthology in alignment file (.ali). 

| Parameter | Description |
|-----------|-------------|
| `arb` | name of file .arb |
| `ali` | name of file .ali |
| `clades` |  file containing the species list for each clade |
| `maxdist_sistergroup` | maximum distance to consider sister-group relationship as suspicious (default=0) |
| `mindist_suspicious` | minimum distance to consider an internal branch length as sufficient to support suspicious taxonomy (default=0) |
| `min_nb_spec_split` | minimum number of species in the smallest subclades to allow splitting (default=1) |
| `min_overlap` | number of amino acids overlapping between two subclades to validate a putative paralogy (default=50) |
| `min_long_BL` | minimum relative branch length to consider a branch as too long (default=50) |
| `paralogy_clade` | false or true (to create separate alignments for clades without out- and in-paralogy and for clades with out- or in-paralogy) |

### Usage
```bash
detect-problems-arb ali=file.ali arb=corresponding_tree.arb clades==Example_file_clade.txt paralogy_clade=false
```
### ali2fasta.c
Reformats '.ali' files in .fasta

| Parameter | Description |
|-----------|-------------|
| `ali` | name of file .ali (optional) |
| `keep_alignment` |yes or no (default = yes) |
| `ident` | space or no (default = space) |
| `min_length` |  minimum number of characters to keep a sequence (default=0) (optional)|
