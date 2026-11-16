###### Microbial 16S rRNA DADA2 Pipeline
##### Archaeal DADA2 analysis for Seepicola viridiplumi microbiome manuscript

## Rice LN, Louca, S, Cordes EE, Young CM. (2026). Trophic niches and the associated
## microbial community for a newly described methane seep sabellid. Deep-Sea Research
## Part 1


## This code is used to take paired end sequencing reads from the sequencing facility
## and run them through the DADA2 pipeline. 
## Full tutorial can be found here: https://benjjneb.github.io/dada2/tutorial.html
## This pipeline includes steps for filtering and trimming primers, low quality regions, 
## identifying learning errors and sampling inferences, merging paired end reads, generating a sequence table,
## removing chimeras, and assigning taxonomy. 

####### Set Up #######

## Set Working directory:
#setwd() #CHANGE ME to desired directory

#Load packages:
library(dada2)
library(dplyr)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
######## Archaeal Run ########

## Define paths to raw (F and R) fasta files

arch <- "/archrun" #CHANGE ME to the directory containing the fastq files

list.files(arch) # Verify the file list

## Perform string manipulation
archFs <- sort(list.files(arch, pattern="_R1_001.fastq", full.names = TRUE))
archRs <- sort(list.files(arch, pattern = "_R2_001.fastq", full.names = TRUE))

## Extract sample names from file list:
archsample.names <- sapply(strsplit(basename(archFs), basename(archRs)), `[`, 1)

filtFsarc <- file.path(arch, paste0(archsample.names,"_F_filt.fastq.gz"))
filtRsarc <- file.path(arch, paste0(archsample.names,"_R_filt.fastq.gz"))


### Assessing Quality ------
png("arcrun_QplotF.png", width = 1600, height = 1200)
plotQualityProfile(archFs)
dev.off

png("arcrun_QplotR.png", width = 1600, height = 1200)
plotQualityProfile(archRs)
dev.off()


#### Filter and Trim -----
outarc <- filterAndTrim(archFs, filtFsarc, archRs, filtRsarc, truncLen = c(270, 210), maxN = 0, maxEE = c(3,3), truncQ = 2,
                         rm.phix = TRUE, compress = TRUE, trimLeft = 18, multithread = FALSE)
head(outarc)
mean(outarc[,2]/outarc[,1])


#### Estimate Error Rates -----
errFarch <- learnErrors(filtFsarc, multithread = TRUE)

errRarch <- learnErrors(filtRsarc, multithread = TRUE)

png("arch_errF.png", width = 1600, height = 1200)
plotErrors(errFarch, nominalQ = TRUE)
dev.off()

png("arch_errR.png", width = 1600, height = 1200)
plotErrors(errRarch, nominalQ = TRUE)
dev.off()


#### Sample inference -----
derepFsarch <- derepFastq(filtFsarc, verbose = TRUE)

derepRsarch <- derepFastq(filtRsarc, verbose = TRUE)

names(derepFsarch) <- archsample.names
names(derepRsarch) <- archsample.names

dadaarchFs <- dada(derepFsarch, err = errFarch, multithread = TRUE)
dadaarchRs <- dada(derepRsarch, err = errRarch, multithread = TRUE)


#### Merge the Paired Reads ----
mergerarch <- mergePairs(dadaarchFs, derepFsarch, dadaarchRs, derepRsarch, minOverlap = 8,trimOverhang = TRUE, verbose = TRUE)
head(mergerarch[[1]])

## Checking the merged amplicon length distribution ----
arctab <- makeSequenceTable(mergerarch)

table(nchar(getSequences(arctab)))

arctab.filt <- arctab[, nchar(colnames(arctab)) %in% 395:425]

#### Create the Sequence Table -----
seqtabArch <- makeSequenceTable(mergerarch)
dim(seqtabArch)   # Gives you info on the number of Samples and ASVs identified in the archaeal run

# Check sequence lengths
table(nchar(getSequences(seqtabArch)))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
######## Chimera Removal ########

seqtabArch.nonchim <- removeBimeraDenovo(seqtabArch, method = "consensus", multithread = TRUE, verbose = TRUE)

dim(seqtabArch.nonchim)
dim(seqtabArch)
sum(seqtabArch.nonchim)/sum(seqtabArch)
# Gives you the percentage of sequences recovered after chimera removal
## This is a good check point as most reads should be retained (>80%).


saveRDS(seqtabArch.nonchim, "archaeal_seqtab.nonchim_merged_11-10.rds")

# Save Tracking File
getN <- function(x) sum(getUniques(x))
trackarch <- cbind(outarc, sapply(dadaarchFs, getN), sapply(dadaarchRs, getN), sapply(mergerarch, getN), 
                   rowSums(seqtabArch.nonchim))
colnames(trackarch) <- c("input","filtered","denoisedF","denoisedR","merged","nonchim")
rownames(trackarch) <- archsample.names
head(trackarch)
write.csv(trackarch, "archaeal_track_reads_11-10.csv")


##### Assign Taxonomy ----
# CHANGE ME to set file path to point to where the silva database is
silva <- "/silv"
silvaspecies <- "/silvspec" 

taxaarc <- assignTaxonomy(seqtabArch.nonchim, silva, multithread = TRUE)

taxaarc.plus <- addSpecies(taxaarc, silvaspecies, verbose = TRUE)

taxaarc.print <- taxaarc.plus
rownames(taxaarc.print) <- NULL
write.csv(taxaarc.print, "Archaeal_taxonomy_11-10.csv")

saveRDS(taxaarc.plus, "Archaeal_taxa_merged_11-10.rds")







