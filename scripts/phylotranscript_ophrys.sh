STRANDED="no"
KMER_LIST="21,39,59,79,99"
THREADS=2
MEMORY="20"
R1="$1"
R2="$2"
BUSCO_LINEAGE="liliopsida_odb10"
ADAPTERS="path/to/adapters/TruSeq3-PE.fa"
EVIGENE_DIR="/path/to/evigene/evigene/scripts"
TRANSABYSS_DIR="/path/to/transabyss/transabyss-master"

############################################
# STEP 1: READ QUALITY ASSESSMENT
############################################

SAMPLE=$(basename "$R1" _R1.fastq.gz)

echo "Running FastQC"

fastqc "$R1" "$R2"

############################################
# STEP 2: READ TRIMMING
############################################

echo "Running Trimmomatic"
java -jar trimmomatic-0.39.jar PE -threads $THREADS $R1 $R2 trim_$R1 unpaired_$R1 trim_$R2 unpaired_$R2 ILLUMINACLIP:$ADAPTERS:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:75 AVGQUAL:25


############################################
# STEP 3a: TRINITY
############################################

echo "Running Trinity"
OUTPUT_TRINITY=Trinity-${SAMPLE}
if [[ "$STRANDED" == "yes" ]]; then
	Trinity --seqType fq --left trim_$R1 --right trim_$R2 --SS_lib_type RF --CPU $THREADS --max_memory ${MEMORY}G --output $OUTPUT_TRINITY
else 
    	Trinity --seqType fq --left trim_$R1 --right trim_$R2 --CPU $THREADS --max_memory ${MEMORY}G --output $OUTPUT_TRINITY
fi


############################################
# STEP 3b: SPADES
############################################

echo "Running Spades"

OUTPUT_SPADES=Spades-${SAMPLE}
if [[ "$STRANDED" == "yes" ]]; then
	spades.py --rna --ss-rf --memory $MEMORY -t $THREADS -k $KMER_LIST -1 trim_$R1 -2 trim_$R2 -o $OUTPUT_SPADES
else
    	spades.py --rna --memory $MEMORY -t $THREADS -k $KMER_LIST -1 trim_$R1 -2 trim_$R2 -o $OUTPUT_SPADES
    fi


############################################
# STEP 3b: trans-ABySS
############################################
KMER=$(echo "$KMER_LIST" | tr ',' ' ')

echo "Running trans-abyss"

OUTPUT_ABYSS=abyss-${SAMPLE}
for K in $KMER; do
       	if [[ "$STRANDED" == "yes" ]]; then    
       		$TRANSABYSS_DIR/transabyss --SS --pe trim_$R1 trim_$R2 --name ${K}_$OUTPUT_ABYSS --kmer $K --outdir ${K}_$OUTPUT_ABYSS --threads $THREADS
	else
		$TRANSABYSS_DIR/transabyss --pe trim_$R1 trim_$R2 --name ${K}_$OUTPUT_ABYSS --kmer $K --outdir ${K}_$OUTPUT_ABYSS --threads $THREADS
	fi
done


############################################
# STEP 4: ASSEMBLY MERGING
############################################
echo "Running Trans-abyss merge"

for K in $KMER; do
cp ${K}_${OUTPUT_ABYSS}/${K}_${OUTPUT_ABYSS}-final.fa ./
done

cp ${OUTPUT_SPADES}/transcripts.fasta ${OUTPUT_SPADES}.fasta

FILES=""
PREFIXES=""
for K in $KMER; do
        FILES="$FILES ${K}_${OUTPUT_ABYSS}-final.fa"
        PREFIXES="$PREFIXES abyss${K}_"
done

FILES="$FILES ${OUTPUT_SPADES}.fasta ${OUTPUT_TRINITY}.Trinity.fasta"
PREFIXES="$PREFIXES spades_ trinity_"

$TRANSABYSS_DIR/transabyss-merge $FILES --mink $(echo $KMER | awk '{print $1}') --maxk $(echo $KMER | awk '{print $NF}') --threads $THREADS --prefixes $PREFIXES --out merge_${SAMPLE}.fasta


############################################
# STEP 5: EVIDENTIALGENE
############################################

## Sequence with unresolved nucleotide (N) removal
awk '/^>/ {if (seq && seq !~ /N/) print header "\n" seq; header=$0; seq=""; next} {seq=seq $0} END {if (seq !~ /N/) print header "\n" seq}' merge_${SAMPLE}.fasta > noN_${SAMPLE}.fasta 


## Sequence formating

$EVIGENE_DIR/rnaseq/trformat.pl noN_${SAMPLE}.fasta > ok_${SAMPLE}.fasta  

mkdir evigene_${SAMPLE}
mv ok_${SAMPLE}.fasta  evigene_${SAMPLE}
cd evigene_${SAMPLE}

## Evidencegene
if [[ "$STRANDED" == "yes" ]]; then
	$EVIGENE_DIR/prot/tr2aacds4.pl -NCPU $THREADS -MAXMEM $MEMORY -log -strand yes -cdna ok_${SAMPLE}.fasta 
else
	$EVIGENE_DIR/prot/tr2aacds4.pl -NCPU $THREADS -MAXMEM $MEMORY -log -strand no -cdna ok_${SAMPLE}.fasta
fi

############################################
# STEP 6: FINAL ASSEMBLY
############################################

##Copying main evigene results
cp okayset/ok_${SAMPLE}.okay.mrna ../
cd ..


awk '/^>/ {if (seq) print seq; print; seq=""; next} {seq=seq $0} END {if (seq) print seq}' ok_${SAMPLE}.okay.mrna > oneline_${SAMPLE}.fasta


##Extracting main and noclass transcripts

awk '/^>/ {header=$0; getline seq; if (header ~ /main/ || header ~ /noclass/) print header "\n" seq}' oneline_${SAMPLE}.fasta > ${SAMPLE}.fasta


############################################
# STEP 6: ASSEMBLY COMPLETNESS ASSESSMENT
############################################

busco -i ${SAMPLE}.fasta -m tran -l $BUSCO_LINEAGE -c 10 -o busco_${SAMPLE}

