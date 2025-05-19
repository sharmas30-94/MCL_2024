#!/bin/sh
#SBATCH --job-name=Bowtie2
#SBATCH --ntasks=1
#SBATCH --mem=24G
#SBATCH --ntasks-per-node=8
#SBATCH --time=20:00:00 
#SBATCH --mail-type=FAIL
#SBATCH --licenses=common
#SBATCH --mail-user=abouska@unmc.edu
#SBATCH --chdir=/lustre/work/javeediqbal/abouska/CutAndRun/Bowtie2_hg38

 
module load bowtie/2.4 samtools/1.3 picard/2.9 bedtools/2.27


## instructions from https://yezhengstat.github.io/CUTTag_tutorial/#

## Build the bowtie2 reference genome index if needed:
## bowtie2-build path/to/hg38/fasta/hg38.fa /path/to/bowtie2Index/hg38


cores=8

####################  hg19  #############################################################################################
#build=hg19
#ref="/common/javeediqbal/shared/Reference.Files/hg19/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome"
#chromSize="/common/javeediqbal/shared/Reference.Files/ChromosomeSizes/ucsc.hg19.chrom.size.txt"

###################  hg38  #############################################################################################

build=hg38
ref="/common/javeediqbal/shared/Reference.Files/hg38/Homo_sapiens/NCBI/GRCh38/Sequence/Bowtie2/genome"
chromSize="/common/javeediqbal/shared/Reference.Files/ChromosomeSizes/ucsc.hg38.chrom.size.txt"

###################  e Coli  #############################################################################################

#spikeInRef="/common/javeediqbal/shared/Reference.Files/Escherichia_coli_K_12_DH10B/NCBI/2008-03-17/Sequence/Bowtie2Index/genome"

# I originally downloaded Escherichia_coli_K_12_DH10B but then andy sent the instructions and it said to use Escherichia_coli_K_12_MG1655

spikeInRef="/common/javeediqbal/shared/Reference.Files/Escherichia_coli_K_12_MG1655/NCBI/2001-10-15/Sequence/Bowtie2Index/genome"
#####################################################################################################################################

projPath="/lustre/work/javeediqbal/abouska/CutAndRun/Bowtie2_hg38"

r1=${1}
r2=${2}
histName=${3}

mkdir -p ${projPath}/alignment/sam/bowtie2_summary
mkdir -p ${projPath}/alignment/bam
mkdir -p ${projPath}/alignment/bed
mkdir -p ${projPath}/alignment/bedgraph


#Aligning to reference genome with Bowtie2
echo "[`date`] Paired Read aligning for ${r1} and ${r2} to ${build} using with Bowtie2"

# use --local option if Reads needed to be trimmed otherwise --end-to-end

bowtie2 --local --very-sensitive --no-mixed --no-discordant --phred33 \
-I 10 -X 700 -p ${cores} -x ${ref} -1 ${r1} -2 ${r2} \
-S ${projPath}/alignment/sam/${histName}_bowtie2.sam &> ${projPath}/alignment/sam/bowtie2_summary/${histName}_bowtie2.txt


#Aligning to Spike-in Control genome with Bowtie2
echo "[`date`] Paired Read aligning for ${r1} and ${r2} to e coli using with Bowtie2"


bowtie2 --local --very-sensitive --no-overlap --no-dovetail --no-mixed --no-discordant --phred33 \
-I 10 -X 700 -p ${cores} -x ${spikeInRef} -1 ${r1} -2 ${r2} \
-S ${projPath}/alignment/sam/${histName}_bowtie2_spikeIn.sam &> ${projPath}/alignment/sam/bowtie2_summary/${histName}_bowtie2_spikeIn.txt

seqDepthDouble=`samtools view -F 0x04 ${projPath}/alignment/sam/${histName}_bowtie2_spikeIn.sam | wc -l`
seqDepth=$((seqDepthDouble/2))
echo $seqDepth >${projPath}/alignment/sam/bowtie2_summary/${histName}_bowtie2_spikeIn.seqDepth




#### Remove Dup  doing what pipeline said, but I think you should be able to just remove duplicate if you want and used the sorted bam if you don't want dup removed

mkdir -p $projPath/alignment/removeDuplicate/picard_summary

## Sort by coordinate
samtools view -bS ${projPath}/alignment/sam/${histName}_bowtie2.sam | samtools sort -o ${projPath}/alignment/sam/${histName}_bowtie2.sorted.bam &&
samtools index ${projPath}/alignment/sam/${histName}_bowtie2.sorted.bam 

## mark duplicates
picard MarkDuplicates I=${projPath}/alignment/sam/${histName}_bowtie2.sorted.bam O=${projPath}/alignment/removeDuplicate/${histName}_bowtie2.sorted.dupMarked.bam M=${projPath}/alignment/removeDuplicate/picard_summary/${histName}_picard.dupMark.txt

## remove duplicates
picard MarkDuplicates I=${projPath}/alignment/sam/${histName}_bowtie2.sorted.bam O=${projPath}/alignment/removeDuplicate/${histName}_bowtie2.sorted.rmDup.bam REMOVE_DUPLICATES=true METRICS_FILE=$projPath/alignment/removeDuplicate/picard_summary/${histName}_picard.rmDup.txt



#### extract Frag length

mkdir -p $projPath/alignment/sam/fragmentLen

## Extract the 9th column from the alignment sam file which is the fragment length
samtools view -F 0x04 ${projPath}/alignment/sam/${histName}_bowtie2.sorted.bam  | awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | sort | uniq -c | awk -v OFS="\t" '{print $2, $1/2}' >$projPath/alignment/sam/fragmentLen/${histName}_fragmentLen.txt



#####  If you want to finler out min Quality score do this...... Not doing at present #############################
minQualityScore=2
#samtools view -q $minQualityScore ${projPath}/alignment/sam/${histName}_bowtie2.sorted.bam  >${projPath}/alignment/sam/${histName}_bowtie2.qualityScore$minQualityScore.sortd.bam

######################################################################################################################




samtools view -bS -F 0x04 ${projPath}/alignment/sam/${histName}_bowtie2.sorted.bam > ${projPath}/alignment/bam/${histName}_bowtie2.mapped.bam

#0x04 is for read unmapped
#also remove not primary alignment and unmapped use 260

samtools view -bS -F 260 ${projPath}/alignment/sam/${histName}_bowtie2.sorted.bam > ${projPath}/alignment/bam/${histName}_bowtie2.mapped.primary.align.bam

samtools index ${projPath}/alignment/bam/${histName}_bowtie2.mapped.primary.align.bam


## Convert into bed file format
bedtools bamtobed -i ${projPath}/alignment/bam/${histName}_bowtie2.mapped.bam -bedpe >$projPath/alignment/bed/${histName}_bowtie2.bed

## Keep the read pairs that are on the same chromosome and fragment length less than 1000bp.
awk '$1==$4 && $6-$2 < 1000 {print $0}' ${projPath}/alignment/bed/${histName}_bowtie2.bed >$projPath/alignment/bed/${histName}_bowtie2.clean.bed

## Only extract the fragment related columns
cut -f 1,2,6 $projPath/alignment/bed/${histName}_bowtie2.clean.bed | sort -k1,1 -k2,2n -k3,3n  >$projPath/alignment/bed/${histName}_bowtie2.fragments.bed

## We use the mid point of each fragment to infer which 500bp bins does this fragment belong to.
binLen=500
awk -v w=$binLen '{print $1, int(($2 + $3)/(2*w))*w + w/2}' ${projPath}/alignment/bed/${histName}_bowtie2.fragments.bed | sort -k1,1V -k2,2n | uniq -c | awk -v OFS="\t" '{print $2, $3, $1}' |  sort -k1,1V -k2,2n  >${projPath}/alignment/bed/${histName}_bowtie2.fragmentsCount.bin$binLen.bed


#### Spike-in Calibration  

if [[ "$seqDepth" -gt "1" ]]; then
    
    mkdir -p $projPath/alignment/bedgraph

    scale_factor=`echo "10000 / $seqDepth" | bc -l`
    echo "Scaling factor for $histName is: $scale_factor!"
    bedtools genomecov -bg -scale $scale_factor -i $projPath/alignment/bed/${histName}_bowtie2.fragments.bed -g $chromSize > $projPath/alignment/bedgraph/${histName}_bowtie2.fragments.normalized.bedgraph
    

fi



###### unNormalize to Spikein Bedgraph

 mkdir -p $projPath/alignment/bedgraph_noSpikeNorm

     bedtools genomecov -bg -i $projPath/alignment/bed/${histName}_bowtie2.fragments.bed -g $chromSize > $projPath/alignment/bedgraph_noSpikeNorm/${histName}_bowtie2.fragments.noSpikeNorm.bedgraph
    




