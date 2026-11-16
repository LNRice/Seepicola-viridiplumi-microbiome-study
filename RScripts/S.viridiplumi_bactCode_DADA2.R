###### Microbial 16S rRNA DADA2 Pipeline
##### DADA2 analysis for Seepicola viridiplumi microbiome manuscript

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
######## Bacterial Run ########

## Define paths to raw (F and R) fasta files

bact <- "/bactrun" #CHANGE ME to the directory containing the fastq files

list.files(bact) # Verify the file list

## Perform string manipulation
bactFs <- sort(list.files(bact, pattern="_R1_001.fastq", full.names = TRUE))
bactRs <- sort(list.files(bact, pattern = "_R2_001.fastq", full.names = TRUE))

## Extract sample names from file list:
bactsample.names <- sapply(strsplit(basename(bactFs), basename(bactRs)), `[`, 1)

filtFsbact <- file.path(bact, paste0(bactsample.names,"_F_filt.fastq.gz"))
filtRsbact <- file.path(bact, paste0(bactsample.names,"_R_filt.fastq.gz"))


### Assessing Quality ------
png("bactrun_QplotF.png", width = 1600, height = 1200)
plotQualityProfile(bactFs)
dev.off

png("bactrun_QplotR.png", width = 1600, height = 1200)
plotQualityProfile(bactRs)
dev.off()


#### Filter and Trim -----
outbact <- filterAndTrim(bactFs, filtFsbact, bactRs, filtRsbact, truncLen = c(270, 232), maxN = 0, maxEE = c(3,3), truncQ = 2,
                         rm.phix = TRUE, compress = TRUE, trimLeft = 17, multithread = FALSE)
head(outbact)
mean(outbact[,2]/outbact[,1])


#### Estimate Error Rates -----
errFbact <- learnErrors(filtFsbact, multithread = TRUE)

errRbact <- learnErrors(filtRsbact, multithread = TRUE)

# Plot the rates
png("bact_errF.png", width = 1600, height = 1200)
plotErrors(errF, nominalQ = TRUE)
dev.off()

png("bact_errR.png", width = 1600, height = 1200)
plotErrors(errR, nominalQ = TRUE)
dev.off()


#### Sample inference -----
derepFsbact <- derepFastq(filtFsbact, verbose = TRUE)

derepRsbact <- derepFastq(filtRsbact, verbose = TRUE)

names(derepFsbact) <- bactsample.names
names(derepRsbact) <- bactsample.names

dadabactFs <- dada(derepFsbact, err = errFbact, multithread = TRUE)
dadabactRs <- dada(derepRsbact, err = errRbact, multithread = TRUE)


#### Merge the Paired Reads ----
mergerbact <- mergePairs(dadabactFs, derepFsbact, dadabactRs, derepRsbact, trimOverhang = TRUE, verbose = TRUE)

head(mergerbact[[1]])


#### Create the Sequence Table -----
seqtabBact <- makeSequenceTable(mergerbact)
dim(seqtabBact)   # Gives you info on the number of Samples and ASVs identified in the bacterial run

# Check sequence lengths
table(nchar(getSequences(seqtabBact)))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
######## Chimera Removal ########

seqtabBact.nonchim <- removeBimeraDenovo(seqtabBact, method = "consensus", multithread = TRUE, verbose = TRUE)

dim(seqtabBact.nonchim)
dim(seqtabBact)
sum(seqtabBact.nonchim)/sum(seqtabBact)
# Gives you the percentage of sequences recovered after chimera removal
## This is a good check point as most reads should be retained (>80%).


saveRDS(seqtabBact.nonchim, "bacterial_seqtab.nonchim_merged.rds")

# Save Tracking File
getN <- function(x) sum(getUniques(x))
trackbact <- cbind(outbact, sapply(dadabactFs, getN), sapply(dadabactRs, getN), sapply(mergerbact, getN), 
                   rowSums(seqtabBact.nonchim))
colnames(trackbact) <- c("input","filtered","denoisedF","denoisedR","merged","nonchim")
rownames(trackbact) <- bactsample.names
head(trackbact)

write.csv(trackbact, "bacterial_track_reads.csv")


##### Assign Taxonomy ----
# Set file path to point to where the silva database is

silva <- "/silv"

silvaspecies <- "/silvspec"

taxabact <- assignTaxonomy(seqtabBact.nonchim, silva, multithread = TRUE)

unname(taxabact)

taxabact.plus <- addSpecies(taxabact, silvaspecies, verbose = TRUE)

taxabact.print <- taxabact.plus  # removing sequence rownames for display only
rownames(taxabact.print) <- NULL

write.csv(taxabact.print, "Bacterial_taxonomy.csv")

saveRDS(taxabact.plus, "Bacterial_taxa_merged.rds")







