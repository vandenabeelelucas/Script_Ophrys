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

