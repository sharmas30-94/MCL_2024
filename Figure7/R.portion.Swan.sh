 module purge
module load mamba

conda activate $COMMON/.conda/envs/deseq2


R

library(dplyr)
library(stringr)
library(ggplot2)
library(viridis)
library(GenomicRanges)
library(chromVAR)
library(DESeq2)
library(ggpubr)
library(corrplot)
library(ChIPseqSpikeInFree)
##=== R command ===## 



## Path to the project and histone list

projPath="/lustre/work/javeediqbal/abouska/CutAndRun/Bowtie2_hg38"

sampleList = c("Maver-Ctrl-H3K4Me3","Maver-Ctrl-IgG","Maver-Ctrl-p53-1","Maver-Ctrl-p53-2","Maver-p53-H3K4Me3","Maver-p53-IgG","Maver-p53-p53-1","Maver-p53-p53-2","Mino-Ctrl-H3K4Me3","Mino-Ctrl-IgG","Mino-Ctrl-p53-1","Mino-Ctrl-p53-2","Mino-p53-H3K4Me3","Mino-p53-IgG","Mino-p53-p53-1","Mino-p53-p53-2","Z138-R248Q-p53-1","Z138-R248Q-p53-2","Z138-R273C-H3k4Me3","Z138-R273C-IgG","Z138-R273C-p53-1","Z138-R273C-p53-2","Z138-wt-p53-1","Z138-wt-p53-2")

histList = c("H3K4Me3","p53-1", "p53-2", "IgG")

## Collect the alignment results from the bowtie2 alignment summary files

Replicate<-c("Rep1","Rep1","Rep1","Rep2","Rep1","Rep1","Rep1","Rep2","Rep1","Rep1","Rep1","Rep2","Rep1","Rep1","Rep1","Rep2","Rep1","Rep2","Rep1","Rep1","Rep1","Rep2","Rep1","Rep2")


alignResult = c()
for(hist in sampleList){
  alignRes = read.table(paste0(projPath, "/alignment/sam/bowtie2_summary/", hist, "_bowtie2.txt"), header = FALSE, fill = TRUE)
  alignRate = substr(alignRes$V1[6], 1, nchar(as.character(alignRes$V1[6]))-1)
  histInfo = strsplit(hist, "_")[[1]]
  alignResult = data.frame(Histone = histInfo[1], 
                           SequencingDepth = alignRes$V1[1] %>% as.character %>% as.numeric, 
                           MappedFragNum_hg38 = alignRes$V1[4] %>% as.character %>% as.numeric + alignRes$V1[5] %>% as.character %>% as.numeric, 
                           AlignmentRate_hg38 = alignRate %>% as.numeric)  %>% rbind(alignResult, .)
}

#alignResult$Histone = factor(alignResult$Histone, levels = histList)
alignResult %>% mutate(AlignmentRate_hg38 = paste0(AlignmentRate_hg38, "%"))

alignResult<-cbind(alignResult,Replicate)
write.table(alignResult,file=paste("hg38alignSummary",Sys.Date(),"txt",sep="."),row.names=F,sep="\t",quote=F)


spikeAlign = c()
for(hist in sampleList){
  spikeRes = read.table(paste0(projPath, "/alignment/sam/bowtie2_summary/", hist, "_bowtie2_spikeIn.txt"), header = FALSE, fill = TRUE)
  alignRate = substr(spikeRes$V1[6], 1, nchar(as.character(spikeRes$V1[6]))-1)
  histInfo = strsplit(hist, "_")[[1]]
  spikeAlign = data.frame(Histone = histInfo[1],
                          SequencingDepth = spikeRes$V1[1] %>% as.character %>% as.numeric, 
                          MappedFragNum_spikeIn = spikeRes$V1[4] %>% as.character %>% as.numeric + spikeRes$V1[5] %>% as.character %>% as.numeric, 
                          AlignmentRate_spikeIn = alignRate %>% as.numeric)  %>% rbind(spikeAlign, .)
}
#spikeAlign$Histone = factor(spikeAlign$Histone, levels = histList)
spikeAlign %>% mutate(AlignmentRate_spikeIn = paste0(AlignmentRate_spikeIn, "%"))

spikeAlign <-cbind(spikeAlign ,Replicate)

write.table(spikeAlign,file=paste("SpikeInAlignSummary",Sys.Date(),"txt",sep="."),row.names=F,sep="\t",quote=F)



alignSummary = left_join(alignResult, spikeAlign, by = c("Histone", "Replicate", "SequencingDepth")) %>%
  mutate(AlignmentRate_hg38 = paste0(AlignmentRate_hg38, "%"), 
         AlignmentRate_spikeIn = paste0(AlignmentRate_spikeIn, "%"))
write.table(alignSummary,file=paste("All.alignSummary",Sys.Date(),"txt",sep="."),row.names=F,sep="\t",quote=F)

############# Summary Figures  ###################33

thm1<-theme(axis.text.x = element_text(angle = 90,hjust=1))
fig3A = alignResult %>% ggplot(aes(x = Histone, y = SequencingDepth/1000000, fill = Histone)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("Sequencing Depth per Million") +
    xlab("") + 
    ggtitle("A. Sequencing Depth")+thm1

fig3B = alignResult %>% ggplot(aes(x = Histone, y = MappedFragNum_hg38/1000000, fill = Histone)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("Mapped Fragments per Million") +
    xlab("") +
    ggtitle("B. Alignable Fragment (hg19)")+thm1

fig3C = alignResult %>% ggplot(aes(x = Histone, y = AlignmentRate_hg38, fill = Histone)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("% of Mapped Fragments") +
    xlab("") +
    ggtitle("C. Alignment Rate (hg19)")+thm1

fig3D = spikeAlign %>% ggplot(aes(x = Histone, y = AlignmentRate_spikeIn, fill = Histone)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("Spike-in Alignment Rate") +
    xlab("") +
    ggtitle("D. Alignment Rate (E.coli)")+thm1

pdf(file=paste("Align.Summary",Sys.Date(),"pdf",sep="."),width=6,height=11,onefile=F)
ggarrange(fig3A, fig3B, fig3C, fig3D, ncol = 2, nrow=2, common.legend = TRUE, legend="bottom")
dev.off()

#### now do my Version

summaryInfo<-read.table("Sample.Summary.txt",header=T,sep="\t")

fig3A.1 = summaryInfo %>% ggplot(aes(x =CellLine.IP , y = SequencingDepth/1000000, fill = IP)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("Sequencing Depth per Million") +
    xlab("") + 
    ggtitle("A. Sequencing Depth")+thm1

fig3B.1 = summaryInfo %>% ggplot(aes(x = CellLine.IP, y = MappedFragNum_hg19/1000000, fill = IP)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("Mapped Fragments per Million") +
    xlab("") +
    ggtitle("B. Alignable Fragment (hg19)")+thm1

fig3C.1 = summaryInfo %>% ggplot(aes(x = CellLine.IP, y = AlignmentRate_hg19, fill = IP)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("% of Mapped Fragments") +
    xlab("") +
    ggtitle("C. Alignment Rate (hg19)")+thm1

fig3D.1 = summaryInfo %>% ggplot(aes(x = CellLine.IP, y = AlignmentRate_spikeIn, fill = IP)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("Spike-in Alignment Rate") +
    xlab("") +
    ggtitle("D. Alignment Rate (E.coli)")+thm1

pdf(file=paste("Alignment.Summary.byCell.line.IP",Sys.Date(),"pdf",sep="."),width=6,height=11,onefile=F)
ggarrange(fig3A.1, fig3B.1, fig3C.1, fig3D.1, ncol = 2, nrow=2, common.legend = TRUE, legend="bottom")
dev.off()




######################### Summarize the duplication information from the picard summary outputs.

dupResult = c()
for(hist in sampleList)
{
  dupRes = read.table(paste0(projPath, "/alignment/removeDuplicate/picard_summary/", hist, "_picard.rmDup.txt"), header = TRUE, fill = TRUE)
  
  histInfo = strsplit(hist, "_")[[1]]
  dupResult = data.frame(Histone = histInfo[1], MappedFragNum_hg19 = dupRes$READ_PAIRS_EXAMINED[1] %>% as.character %>% as.numeric, DuplicationRate = dupRes$PERCENT_DUPLICATION[1] %>% as.character %>% as.numeric * 100, EstimatedLibrarySize = dupRes$ESTIMATED_LIBRARY_SIZE[1] %>% as.character %>% as.numeric) %>% mutate(UniqueFragNum = MappedFragNum_hg19 * (1-DuplicationRate/100))  %>% rbind(dupResult, .)
}
dupResult$Histone = factor(dupResult$Histone)
alignDupSummary = left_join(alignSummary, dupResult, by = c("Histone", "MappedFragNum_hg19")) %>% mutate(DuplicationRate = paste0(DuplicationRate, "%"))
alignDupSummary
write.table(alignDupSummary,file=paste("All.alignSummary.WithDup",Sys.Date(),"txt",sep="."),row.names=F,sep="\t",quote=F)



## generate boxplot figure for the  duplication rate
fig4A = dupResult %>% ggplot(aes(x = Histone, y = DuplicationRate, fill = Histone)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("Duplication Rate (*100%)") +
    xlab("") +thm1

fig4B = dupResult %>% ggplot(aes(x = Histone, y = EstimatedLibrarySize, fill = Histone)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("Estimated Library Size") +
    xlab("") +thm1


fig4C = dupResult %>% ggplot(aes(x = Histone, y = UniqueFragNum, fill = Histone)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("# of Unique Fragments") +
    xlab("")+thm1

pdf(file=paste("Dup.Summary",Sys.Date(),"pdf",sep="."),width=8.5,height=5.6,onefile=F)
ggarrange(fig4A, fig4B, fig4C, ncol = 3, common.legend = TRUE, legend="bottom")
dev.off()



### my mod

fig4A.1 = summaryInfo %>% ggplot(aes(x = CellLine.IP, y = DuplicationRate, fill = IP)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("Duplication Rate") +
    xlab("") 

fig4B.1 = summaryInfo %>% ggplot(aes(x = CellLine.IP, y = EstimatedLibrarySize, fill = IP)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("Estimated Library Size") +
    xlab("") 

fig4C.1 = summaryInfo %>% ggplot(aes(x = CellLine.IP, y = UniqueFragNum, fill = IP)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("# of Unique Fragments") +
    xlab("")


pdf(file=paste("Dup.Summary.myVersion",Sys.Date(),"pdf",sep="."),width=8.5,height=5.6,onefile=F)
ggarrange(fig4A.1, fig4B.1, fig4C.1, ncol = 3, common.legend = TRUE, legend="bottom")
dev.off()




#### fragment size info



fragLen = c()
for(hist in sampleList){
  
  histInfo = strsplit(hist, "_")[[1]]
  fragLen = read.table(paste0(projPath, "/alignment/sam/fragmentLen/", hist, "_fragmentLen.txt"), header = FALSE) %>% mutate(fragLen = V1 %>% as.numeric, fragCount = V2 %>% as.numeric, Weight = as.numeric(V2)/sum(as.numeric(V2)), Histone = histInfo[1], Replicate = histInfo[2], sampleInfo = hist) %>% rbind(fragLen, .) 
}
fragLen$sampleInfo = factor(fragLen$sampleInfo)
fragLen$Histone = factor(fragLen$Histone)
## Generate the fragment size density plot (violin plot)
fig5A = fragLen %>% ggplot(aes(x = sampleInfo, y = fragLen, weight = Weight, fill = Histone)) +
    geom_violin(bw = 5) +
    scale_y_continuous(breaks = seq(0, 800, 50)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 20) +
    ggpubr::rotate_x_text(angle = 20) +
    ylab("Fragment Length") +
    xlab("")

fig5B = fragLen %>% ggplot(aes(x = fragLen, y = fragCount, color = Histone, group = sampleInfo, linetype = Replicate)) +
  geom_line(size = 1) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma") +
  theme_bw(base_size = 20) +
  xlab("Fragment Length") +
  ylab("Count") +
  coord_cartesian(xlim = c(0, 500))

pdf(file=paste("FragLengthSummary",Sys.Date(),"pdf",sep="."),width=8.5,height=5.6,onefile=F)

ggarrange(fig5A, fig5B, ncol = 2)

dev.off()



pdf(file=paste("FragLengthSummary.v2",Sys.Date(),"pdf",sep="."),width=8.5,height=10,onefile=F)

ggarrange(fig5A+facet_grid(~IP,scales="free",space="free"), fig5B+facet_grid(~IP,scales="free",space="free"), ncol = 1)

dev.off()
pdf(file=paste("FragLength.line",Sys.Date(),"pdf",sep="."),width=14,height=4.5,onefile=F)

fig5B+facet_grid(~IP,scales="free",space="free")

dev.off()
pdf(file=paste("FragLength.line",Sys.Date(),"pdf",sep="."),width=35,height=4.5,onefile=F)

fig5B+facet_grid(~IP,scales="free",space="free")+theme_Pub

dev.off()

pdf(file=paste("FragLength.line.vertical",Sys.Date(),"pdf",sep="."),width=11,height=11,onefile=F)

fig5B+facet_grid(IP~.,scales="free",space="free")+theme_Pub

dev.off()


reprod = c()
fragCount = NULL
for(hist in sampleList){
  
  if(is.null(fragCount)){
    
    fragCount = read.table(paste0(projPath, "/alignment/bed/", hist, "_bowtie2.fragmentsCount.bin500.bed"), header = FALSE) 
    colnames(fragCount) = c("chrom", "bin", hist)
  
  }else{
    
    fragCountTmp = read.table(paste0(projPath, "/alignment/bed/", hist, "_bowtie2.fragmentsCount.bin500.bed"), header = FALSE)
    colnames(fragCountTmp) = c("chrom", "bin", hist)
    fragCount = full_join(fragCount, fragCountTmp, by = c("chrom", "bin"))
    
  }
}

M = cor(fragCount %>% select(-c("chrom", "bin")) %>% log2(), use = "complete.obs") 

pdf(file=paste("Corplot.fracCount",Sys.Date(),"pdf",sep="."),width=8,height=8.5,onefile=F)
corrplot(M, method = "color", outline = T, addgrid.col = "darkgray", order="hclust", addrect = 3, rect.col = "black", rect.lwd = 3,cl.pos = "b", tl.col = "indianred4", tl.cex = 1, cl.cex = 1, addCoef.col = "black", number.digits = 2, number.cex = 1, col = colorRampPalette(c("midnightblue","white","darkred"))(100))
dev.off()




######  summarize Seacr

peakN = c()
peakWidth = c()
peakType = c("control", "top0.01")
for(hist in sampleList){
  histInfo = strsplit(hist, "_")[[1]]
  if(histInfo[1] != "IgG"){
    for(type in peakType){
      peakInfo = read.table(paste0(projPath, "/peakCalling/SEACR/", hist, "_seacr_", type, ".peaks.stringent.bed"), header = FALSE, fill = TRUE)  %>% mutate(width = abs(V3-V2))
      peakN = data.frame(peakN = nrow(peakInfo), peakType = type, Histone = histInfo[1], Replicate = Replicate) %>% rbind(peakN, .)
      peakWidth = data.frame(width = peakInfo$width, peakType = type, Histone = histInfo[1], Replicate = Replicate)  %>% rbind(peakWidth, .)
    }
  }
}
peakN %>% select(Histone, Replicate, peakType, peakN)


sampleList2<-sampleList[c(1,3:5,7:9,11:13,15:19,21:24)]


peakN = c()
peakWidth = c()
peakType = c("control", "top0.01")
for(hist in sampleList2){
{
    for(type in peakType){
      peakInfo = read.table(paste0(projPath, "/peakCalling/SEACR/", hist, "_seacr_", type, ".peaks.stringent.bed"), header = FALSE, fill = TRUE)  %>% mutate(width = abs(V3-V2))
      peakN = data.frame(peakN = nrow(peakInfo), peakType = type) %>% rbind(peakN, .)
      peakWidth = data.frame(width = peakInfo$width, peakType = type)  %>% rbind(peakWidth, .)
    }
  }
}
peakN %>% select(Histone, peakType, peakN)




#### reproducibility between peaks


histL = c("Mino-Ctrl-p53","Mino-p53-p53","Maver-Ctrl-p53","Maver-p53-p53","Z138-R248Q-p53","Z138-R273C-p53","Z138-wt-p53")

histL = c("Mino-Ctrl-p53","Mino-p53-p53","Maver-Ctrl-p53","Maver-p53-p53","Z138-wt-p53")
repL = c( 1:2)
peakType = c("control", "top0.01")
peakOverlap = c()
for(type in peakType){
  for(hist in histL){
    overlap.gr = GRanges()
    for(rep in repL){
      peakInfo = read.table(paste0(projPath, "/peakCalling/SEACR/", hist, "-", rep, "_seacr_", type, ".peaks.stringent.bed"), header = FALSE, fill = TRUE)
      peakInfo.gr = GRanges(peakInfo$V1, IRanges(start = peakInfo$V2, end = peakInfo$V3), strand = "*")
      if(length(overlap.gr) >0){
        overlap.gr = overlap.gr[findOverlaps(overlap.gr, peakInfo.gr)@from]
      }else{
        overlap.gr = peakInfo.gr
        
      }
    }
    peakOverlap = data.frame(peakReprod = length(overlap.gr), Histone = hist, peakType = type) %>% rbind(peakOverlap, .)
  }
}

peakReprod = left_join(peakN, peakOverlap, by = c("Histone", "peakType")) %>% mutate(peakReprodRate = peakReprod/peakN * 100)
peakReprod %>% select(Histone, Replicate, peakType, peakN, peakReprodNum = peakReprod, peakReprodRate)


files<-read.table("file.for.peak.overlap.chk.txt",header=T,sep="\t")


overlapInfo =NULL

for(i in 1:nrow(files))
{

peakInfo = read.table(files[i,1], header = FALSE, fill = TRUE)
peakInfo2 = read.table(files[i,2], header = FALSE, fill = TRUE)

peakInfo.gr = GRanges(peakInfo$V1, IRanges(start = peakInfo$V2, end = peakInfo$V3), strand = "*")
peakInfo.gr2 = GRanges(peakInfo2$V1, IRanges(start = peakInfo2$V2, end = peakInfo2$V3), strand = "*")

overlap.gr<-try(findOverlaps(peakInfo.gr, peakInfo.gr2))
overlap <-length(overlap.gr)
overlap2<-cbind(files[i,3],overlap)
overlapInfo<-rbind(overlapInfo,overlap2 )

}

overlap.gr = overlap.gr[findOverlaps(peakInfo.gr, peakInfo.gr2)@from]



##=== R command ===## Create the peak×sample  matrix.

## only doing p53 IPs:

histL = c("Mino-Ctrl-p53","Mino-p53-p53","Maver-Ctrl-p53","Maver-p53-p53","Z138-R248Q-p53","Z138-R273C-p53","Z138-wt-p53")
repL = c( 1:2)

mPeak<-NULL
masterPeak <-NULL
mPeak = GRanges()
## overlap with bam file to get count
for(hist in histL){
  for(rep in repL){
    peakRes = try(read.table(paste0(projPath, "/peakCalling/SEACR/", hist, "-", rep, "_seacr_control.peaks.stringent.bed"), header = FALSE, fill = TRUE))
if(!is.null(nrow(peakRes))){
print(hist)
print(rep)
print(nrow(peakRes))
  mPeak =   GRanges(seqnames = peakRes$V1, IRanges(start = peakRes$V2, end = peakRes$V3), strand = "*") %>% append(mPeak, .)}
  }
}
masterPeak = reduce(mPeak)





##=== R command ===##  Get the fragment counts for each peak in the master peak list. 


library(DESeq2)
bamDir = paste0(projPath, "/alignment/bam")
countMat = matrix(NA, length(masterPeak), length(histL)*length(repL))
## overlap with bam file to get count

allNames<-NULL
i = 1
for(hist in histL){
  for(rep in repL){

    bamFile = paste0(bamDir, "/", hist, "-", rep, "_bowtie2.mapped.bam")
name=paste0( hist, "-", rep, "_bowtie2.mapped.bam")
print(name)
allNames<-c(allNames,name)
    fragment_counts <- getCounts(bamFile, masterPeak, paired = TRUE, by_rg = FALSE, format = "bam")
    countMat[, i] = counts(fragment_counts)[,1]
    i = i + 1
  }
}

#colnames(countMat) = paste(rep(histL, 2), rep(repL, each = 2), sep = "-")This code is from the website but is is Not correct

colnames(countMat) = allNames

df1 <- as.data.frame(masterPeak)
write.table(cbind(df1,countMat),"fragmentCounts.for.merged.p53IPs.seacr_control.peaks.stringent.bed.txt",row.names=F,sep="\t",quote=F)




###   Sequencing depth normalization and differential enriched peaks detection

selectR = which(rowSums(countMat) > 5) ## remove low count genes
dataS = countMat[selectR,]
condition = factor(rep(histL, each = length(repL)))
dds = DESeqDataSetFromMatrix(countData = dataS,
                              colData = DataFrame(condition),
                              design = ~ condition)
DDS = DESeq(dds)
normDDS = counts(DDS, normalized = TRUE) ## normalization with respect to the sequencing depth
colnames(normDDS) = paste0(colnames(normDDS), "_norm")

write.table(cbind(df1,normDDS),"SeqDepthNormalizedFragmentCounts.for.merged.p53IPs.seacr_control.peaks.stringent.bed.txt",row.names=F,sep="\t",quote=F)

res = results(DDS, independentFiltering = FALSE, altHypothesis = "greaterAbs")

# Export wanted DDS results
resultsNames(DDS)
write.table(results(DDS, name="condition_Maver.p53.p53_vs_Maver.Ctrl.p53",independentFiltering = FALSE,cooksCutoff=F),"test.txt",row.names=T,sep="\t")


DDS$condition <- relevel(DDS$condition, ref = "Mino-Ctrl-p53")
DDS = DESeq(DDS)
res = results(DDS, independentFiltering = FALSE, altHypothesis = "greaterAbs")
resultsNames(DDS)
write.table(results(DDS, name="condition_Mino.p53.p53_vs_Mino.Ctrl.p53"),"clipboard-500000",row.names=T,sep="\t")



DDS$condition <- relevel(DDS$condition, ref = "Z138-wt-p53")
DDS = DESeq(DDS)
res = results(DDS, independentFiltering = FALSE, altHypothesis = "greaterAbs")
resultsNames(DDS)
write.table(results(DDS, name="condition_Z138.R248Q.p53_vs_Z138.wt.p53",independentFiltering = FALSE,cooksCutoff=F),"temp.txt",row.names=T,sep="\t")
write.table(results(DDS, name="condition_Z138.R273C.p53_vs_Z138.wt.p53",independentFiltering = FALSE,cooksCutoff=F),"temp.txt",row.names=T,sep="\t")


### Now use Cell Lines as Replicates

condition = factor(rep(histL, each = length(repL)))
dds = DESeqDataSetFromMatrix(countData = dataS,
                              colData = DataFrame(condition),
                              design = ~ condition)
DDS = DESeq(dds)
normDDS = counts(DDS, normalized = TRUE) ## normalization with respect to the sequencing depth
colnames(normDDS) = paste0(colnames(normDDS), "_norm")

write.table(cbind(df1,normDDS),"SeqDepthNormalizedFragmentCounts.for.merged.p53IPs.seacr_control.peaks.stringent.bed.txt",row.names=F,sep="\t",quote=F)

res = results(DDS, independentFiltering = FALSE, altHypothesis = "greaterAbs")

# Export wanted DDS results
resultsNames(DDS)
write.table(results(DDS, name="condition_Maver.p53.p53_vs_Maver.Ctrl.p53",independentFiltering = FALSE,cooksCutoff=F),"test.txt",row.names=T,sep="\t")








#############################    ##################

##=== R command ===## Create the peak×sample  matrix.

## only doing p53 IPs:
projPath="/lustre/work/javeediqbal/abouska/CutAndRun/Bowtie2_hg38"

histL = c("Mino-Ctrl-p53","Mino-p53-p53","Maver-Ctrl-p53","Maver-p53-p53","Z138-R248Q-p53","Z138-R273C-p53","Z138-wt-p53")
repL = c( 1:2)

mPeak<-NULL
masterPeak <-NULL
mPeak = GRanges()
## overlap with bam file to get count
for(hist in histL){
  for(rep in repL){
    peakRes = try(read.table(paste0(projPath, "/peakCalling/SEACR/", hist, "-", rep, "_seacr_top0.01.peaks.stringent.bed"), header = FALSE, fill = TRUE))
if(!is.null(nrow(peakRes))){
print(hist)
print(rep)
print(nrow(peakRes))
  mPeak =   GRanges(seqnames = peakRes$V1, IRanges(start = peakRes$V2, end = peakRes$V3), strand = "*") %>% append(mPeak, .)}
  }
}
masterPeak = reduce(mPeak)





##=== R command ===##  Get the fragment counts for each peak in the master peak list. 


library(DESeq2)
bamDir = paste0(projPath, "/alignment/bam")
countMat = matrix(NA, length(masterPeak), length(histL)*length(repL))
## overlap with bam file to get count

allNames<-NULL
i = 1
for(hist in histL){
  for(rep in repL){

    bamFile = paste0(bamDir, "/", hist, "-", rep, "_bowtie2.mapped.bam")
name=paste0( hist, "-", rep, "_bowtie2.mapped.bam")
print(name)
allNames<-c(allNames,name)
    fragment_counts <- getCounts(bamFile, masterPeak, paired = TRUE, by_rg = FALSE, format = "bam")
    countMat[, i] = counts(fragment_counts)[,1]
    i = i + 1
  }
}

#colnames(countMat) = paste(rep(histL, 2), rep(repL, each = 2), sep = "-")This code is from the website but is is Not correct

colnames(countMat) = allNames

df1 <- as.data.frame(masterPeak)
write.table(cbind(df1,countMat),"fragmentCounts.for.merged.p53IPs._seacr_top0.01.peaks.stringent.bed.txt",row.names=F,sep="\t",quote=F)




###   Sequencing depth normalization and differential enriched peaks detection

#selectR = which(rowSums(countMat) > 5) ## remove low count genes
#dataS = countMat[selectR,]
dataS<-countMat
condition = factor(rep(histL, each = length(repL)))
dds = DESeqDataSetFromMatrix(countData = dataS,
                              colData = DataFrame(condition),
                              design = ~ condition)

DDS = DESeq(dds)
normDDS = counts(DDS, normalized = TRUE) ## normalization with respect to the sequencing depth
colnames(normDDS) = paste0(colnames(normDDS), "_norm")

write.table(cbind(df1,normDDS),"SeqDepthNormalizedFragmentCounts.for.merged.p53IPs.seacr_top0.01.peaks.stringent.bed.txt",row.names=F,sep="\t",quote=F)

res = results(DDS, cooksCutoff=FALSE,independentFiltering = F, altHypothesis = "greaterAbs")



# Export wanted DDS results
resultsNames(DDS)
write.table(results(DDS,  cooksCutoff=FALSE,independentFiltering = F,name="condition_Maver.p53.p53_vs_Maver.Ctrl.p53"),"test.txt",row.names=T,sep="\t")


DDS$condition <- relevel(DDS$condition, ref = "Mino-Ctrl-p53")
DDS = DESeq(DDS)
res = results(DDS, cooksCutoff=FALSE,independentFiltering = FALSE, altHypothesis = "greaterAbs")
resultsNames(DDS)
write.table(results(DDS, name="condition_Mino.p53.p53_vs_Mino.Ctrl.p53",cooksCutoff=FALSE,independentFiltering = FALSE),"temp.txt",row.names=F,sep="\t")



DDS$condition <- relevel(DDS$condition, ref = "Z138-wt-p53")
DDS = DESeq(DDS)
res = results(DDS, independentFiltering = FALSE, altHypothesis = "greaterAbs")
resultsNames(DDS)
write.table(results(DDS, name="condition_Z138.R248Q.p53_vs_Z138.wt.p53",cooksCutoff=FALSE,independentFiltering = FALSE),"temp.txt",row.names=T,sep="\t")
write.table(results(DDS, name="condition_Z138.R273C.p53_vs_Z138.wt.p53",cooksCutoff=FALSE,independentFiltering = FALSE),"temp.txt",row.names=T,sep="\t")


DDS$condition <- relevel(DDS$condition, ref = "Z138-R248Q-p53")

DDS = DESeq(DDS)
write.table(results(DDS, name="condition_Z138.R273C.p53_vs_Z138.R248Q.p53",cooksCutoff=FALSE,independentFiltering = FALSE),"temp.txt",row.names=T,sep="\t")

### now vs no Control

dat<-read.table("fragmentCounts.for.merged.p53IPs._seacr_top0.01.peaks.stringent.bed.txt",header=T,sep="\t")
dataS<-dat[,6:19]
condition = factor(rep(histL, each = length(repL)))
dds = DESeqDataSetFromMatrix(countData = dataS,
                              colData = DataFrame(condition),
                              design = ~ condition)
dds$condition <- relevel(dds$condition, ref = "Z138-R248Q-p53")
DDS = DESeq(dds)
res = results(DDS, independentFiltering = FALSE, altHypothesis = "greaterAbs")
resultsNames(DDS)
write.table(results(DDS, name="condition_Z138.R273C.p53_vs_Z138.R248Q.p53"),"temp.txt",row.names=T,sep="\t")




### Now use Cell Lines as Replicates

condition = factor(rep(histL, each = length(repL)))
dds = DESeqDataSetFromMatrix(countData = dataS,
                              colData = DataFrame(condition),
                              design = ~ condition)
DDS = DESeq(dds)
normDDS = counts(DDS, normalized = TRUE) ## normalization with respect to the sequencing depth
colnames(normDDS) = paste0(colnames(normDDS), "_norm")

write.table(cbind(df1,normDDS),"SeqDepthNormalizedFragmentCounts.for.merged.p53IPs.seacr_control.peaks.stringent.bed.txt",row.names=F,sep="\t",quote=F)

res = results(DDS, independentFiltering = FALSE, altHypothesis = "greaterAbs")

# Export wanted DDS results
resultsNames(DDS)
write.table(results(DDS, name="condition_Maver.p53.p53_vs_Maver.Ctrl.p53"),"test.txt",row.names=T,sep="\t")



#### ################     trying on Counts that I normalized to the individual spikeins

###   Sequencing depth normalization and differential enriched peaks detection
setwd( "/lustre/work/javeediqbal/abouska/CutAndRun/Bowtie2_hg38/alignment")
library(DESeq2)

histL = c("Mino-Ctrl-p53","Mino-p53-p53","Maver-Ctrl-p53","Maver-p53-p53","Z138-R248Q-p53","Z138-R273C-p53","Z138-wt-p53")
repL = c( 1:2)

dataS<-read.table("SpikeNorm.fragmentCounts.for.merged.p53IPs.seacr_control.peaks.stringent.bed.txt",header=T,sep="\t")
dataS<-round(dataS,digits=0)

condition = factor(rep(histL, each = length(repL)))

dds = DESeqDataSetFromMatrix(countData = dataS,colData = DataFrame(condition),design = ~ condition)
DDS = DESeq(dds)
normDDS = counts(DDS, normalized = TRUE) ## normalization with respect to the sequencing depth
colnames(normDDS) = paste0(colnames(normDDS), "_norm")

write.table(cbind(df1,normDDS),"ScaleFactorAdj.then.SeqDepthNormalizedFragmentCounts.for.merged.p53IPs.macs2.v.IgG_peak_q0.1_peaks.narrowPeak.txt",row.names=F,sep="\t",quote=F)

res = results(DDS, independentFiltering = FALSE, altHypothesis = "greaterAbs")

######    Export wanted DDS results   ###############

resultsNames(DDS)

a<-results(DDS, name="condition_Maver.p53.p53_vs_Maver.Ctrl.p53")


DDS$condition <- relevel(DDS$condition, ref = "Mino-Ctrl-p53")
DDS = DESeq(DDS)
res = results(DDS, independentFiltering = FALSE, altHypothesis = "greaterAbs")
resultsNames(DDS)

b<-results(DDS,cooksCutoff=FALSE,independentFiltering = F, altHypothesis = "greaterAbs", name="condition_Mino.p53.p53_vs_Mino.Ctrl.p53")


DDS$condition <- relevel(DDS$condition, ref = "Z138-wt-p53")
DDS = DESeq(DDS)
res = results(DDS, independentFiltering = FALSE, altHypothesis = "greaterAbs")
resultsNames(DDS)
c<-results(DDS, cooksCutoff=FALSE,independentFiltering = F, altHypothesis = "greaterAbs",name="condition_Z138.R248Q.p53_vs_Z138.wt.p53")

d<-(results(DDS, cooksCutoff=FALSE,independentFiltering = F, altHypothesis = "greaterAbs",name="condition_Z138.R273C.p53_vs_Z138.wt.p53")



DDS$condition <- relevel(DDS$condition, ref = "Z138-R248Q-p53")
DDS = DESeq(DDS)
res = results(DDS, independentFiltering = FALSE, altHypothesis = "greaterAbs")
resultsNames(DDS)

e<-results(DDS,cooksCutoff=FALSE,independentFiltering = F, altHypothesis = "greaterAbs", name="condition_Z138.R273C.p53_vs_Z138.R248Q.p53")



combo<-cbind(a,b,c,d,e)
write.table(combo,"Deseq2Results.SpikeAdjustedCounts.merged.p53IPs.seacr_control.peaks.stringent.bed.txt",row.names=F,sep="\t")



##############  vs Control  #################################
projPath="/lustre/work/javeediqbal/abouska/CutAndRun/Bowtie2_hg38"

histL = c("Mino-Ctrl-p53","Mino-p53-p53","Maver-Ctrl-p53","Maver-p53-p53","Z138-R248Q-p53","Z138-R273C-p53","Z138-wt-p53")
repL = c( 1:2)
# have to take out Z138-R273C-p53 beause it has no peaks
histL2 = c("Mino-Ctrl-p53","Mino-p53-p53","Maver-Ctrl-p53","Maver-p53-p53","Z138-R248Q-p53","Z138-wt-p53")

mPeak<-NULL
masterPeak <-NULL
mPeak = GRanges()
## overlap with bam file to get count
for(hist in histL2){
  for(rep in repL){
    peakRes = try(read.table(paste0(projPath, "/peakCalling/SEACR/", hist, "-", rep, "_seacr_control.peaks.stringent.bed"), header = FALSE, fill = TRUE))
if(!is.null(nrow(peakRes))){
print(hist)
print(rep)
print(nrow(peakRes))
  mPeak =   GRanges(seqnames = peakRes$V1, IRanges(start = peakRes$V2, end = peakRes$V3), strand = "*") %>% append(mPeak, .)}
  }
}
masterPeak = reduce(mPeak)


bamDir = paste0(projPath, "/alignment/bam")
countMat = matrix(NA, length(masterPeak), length(histL)*length(repL))
## overlap with bam file to get count

allNames<-NULL
i = 1
for(hist in histL){
  for(rep in repL){

    bamFile = paste0(bamDir, "/", hist, "-", rep, "_bowtie2.mapped.bam")
name=paste0( hist, "-", rep, "_bowtie2.mapped.bam")
print(name)
allNames<-c(allNames,name)
    fragment_counts <- getCounts(bamFile, masterPeak, paired = TRUE, by_rg = FALSE, format = "bam")
    countMat[, i] = counts(fragment_counts)[,1]
    i = i + 1
  }
}

#colnames(countMat) = paste(rep(histL, 2), rep(repL, each = 2), sep = "-")This code is from the website but is is Not correct

colnames(countMat) = allNames

df1 <- as.data.frame(masterPeak)
write.table(cbind(df1,countMat),"fragmentCounts.for.merged.p53IPs.seacr_control.peaks.stringent.bed.txt",row.names=F,sep="\t",quote=F)




dataS<-countMat
condition = factor(rep(histL, each = length(repL)))
dds = DESeqDataSetFromMatrix(countData = dataS,
                              colData = DataFrame(condition),
                              design = ~ condition)

DDS = DESeq(dds)
normDDS = counts(DDS, normalized = TRUE) ## normalization with respect to the sequencing depth
colnames(normDDS) = paste0(colnames(normDDS), "_norm")

write.table(cbind(df1,normDDS),"SeqDepthNormalizedFragmentCounts.for.merged.p53IPs.seacr_control.peaks.stringent.bed.txt",row.names=F,sep="\t",quote=F)




conda deactivate



#########    5/9/23



###### doing additional comparisons combining cell lines
library(DESeq2)

countMat<-read.table("fragmentCounts.for.merged.p53IPs._seacr_top0.01.peaks.stringent.bed.txt",header=T,sep="\t",check.names=F)

countMat<-countMat[,6:19]

info<-read.table("Sample.info.table.txt",header=T,sep="\t")

dataS<-countMat

condition<-info$comparison2

dds = DESeqDataSetFromMatrix(countData = dataS,
                              colData = DataFrame(condition),
                              design = ~ condition)
DDS = DESeq(dds)

resultsNames(DDS)


a<-results(DDS,cooksCutoff=FALSE,independentFiltering = F, altHypothesis = "greaterAbs" ,name="condition_MinoOrMaver.p53.p53_vs_MinoOrMaver.Ctrl.p53")
write.table(a,"MinoOrMaver.Deseq2Results._seacr_top0.01.peaks.stringent.bed.txt",row.names=F,sep="\t")


### vs control
countMat<-read.table("fragmentCounts.for.merged.p53IPs.seacr_control.peaks.stringent.bed.txt",header=T,sep="\t",check.names=F)

countMat<-countMat[,6:19]

dataS<-countMat

condition<-info$comparison2

dds = DESeqDataSetFromMatrix(countData = dataS,
                              colData = DataFrame(condition),
                              design = ~ condition)
DDS = DESeq(dds)

resultsNames(DDS)


a<-results(DDS,cooksCutoff=FALSE,independentFiltering = F, altHypothesis = "greaterAbs" ,name="condition_MinoOrMaver.p53.p53_vs_MinoOrMaver.Ctrl.p53")
write.table(a,"MinoOrMaver.Deseq2Results._seacr_control.peaks.stringent.bed.txt",row.names=F,sep="\t")


###  library(DESeq2) version if you want to exclude your counts

# Assuming you have already loaded and prepared your count matrix and metadata

# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = metadata, design = ~ condition)

# Filter out genes with zero counts
dds <- filterByExpr(dds, min_count = 1)

# Perform differential expression analysis
dds <- DESeq(dds)

# Continue with downstream analysis (e.g., results, visualization, etc.)

