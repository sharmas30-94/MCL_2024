#!/bin/sh
#SBATCH --job-name=Bowtie2
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --ntasks-per-node=1
#SBATCH --time=1:30:00 
#SBATCH --mail-type=FAIL
#SBATCH --licenses=common
#SBATCH --mail-user=abouska@unmc.edu
#SBATCH --chdir=/lustre/work/javeediqbal/abouska/CutAndRun/Bowtie2_hg38

 

module load bowtie/2.4 samtools/1.3 picard/2.9  bedtools/2.27 R/4.2


seacr="$COMMON/PROGRAMS/SEACR/SEACR_1.3.sh"

mkdir -p ${projPath}/peakCalling/SEACR

projPath="/lustre/work/javeediqbal/abouska/CutAndRun/Bowtie2_hg38"

histName=${1}
histControl=${2}
	
bedgph=${projPath}/alignment/bedgraph_noSpikeNorm/${histName}_bowtie2.fragments.noSpikeNorm.bedgraph
controlBed=${projPath}/alignment/bedgraph_noSpikeNorm/${histControl}_bowtie2.fragments.noSpikeNorm.bedgraph
outPath=${projPath}/peakCalling/SEACR_noSpikNorm

mkdir -p $outPath

#############  SEAKR Call peaks ############################


$seacr ${bedgph} \
     ${controlBed} \
     non stringent ${outPath}/${histName}_seacr_control.peaks

$seacr ${bedgph} 0.01 non stringent ${outPath}/${histName}_seacr_top0.01.peaks


