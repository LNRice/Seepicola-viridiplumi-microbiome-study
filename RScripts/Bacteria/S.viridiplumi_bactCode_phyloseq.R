###### Microbial 16S rRNA post-DADA2 workflow, creating Phyloseq and data cleaning ########
##### Bacterial analysis for Seepicola viridiplumi microbiome manuscript

## Rice LN, Louca, S, Cordes EE, Young CM. (2026). Trophic niches and the associated
## microbial community for a newly described methane seep sabellid. Deep-Sea Research
## Part 1



## This code is used to take the raw sequence table produced in the DADA2 pipeline,
## create a phyloseq object, and remove Eukaryotic, chloroplast, and mitochondrial sequences.
## The final output from this workflow will be used in all downstream analysis.

####### Set Up #######

## Set Working directory:
#setwd() #CHANGE ME to desired directory

#Load packages:
library(phyloseq)
library(Biostrings)
library(DECIPHER)
library(MiscMetabar)
library(speedyseq)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
######## Creating Phyloseq Object ########

## Open the needed RDS files
bactseqtab <- readRDS("bacterial_seqtab.nonchim_merged.rds")
bactaxa <- readRDS("Bacterial_taxa_merged.rds")

# bringing in sample metadata
sampdat <- read.table(file = "samplemetadata.csv", sep = ",", header = TRUE, row.names = 1)
rownames(sampdat) <- sampdat$Abbreviated.Sample.Name
rownames(bactseqtab) <- sampdat$Abbreviated.Sample.Name

# Make phyloseq obj
bac_data_raw <- phyloseq(otu_table(bactseqtab, taxa_are_rows = FALSE),
                         sample_data(sampdat),
                         tax_table(bactaxa))
bac_data_raw
sample_sums(bac_data_raw)

### Identify similar ASVs and cluster them 
asv_sequences <- colnames(bactseqtab)
names <- rownames(bactseqtab)
dna <- Biostrings::DNAStringSet(asv_sequences)

nproc = 4

## Align the sequence data
algn <- DECIPHER::AlignSeqs(dna, processors = nproc)

## Generate similarity distance matrix
dist <- DECIPHER::DistanceMatrix(algn, processors = nproc)

## Cluster OTUs with 97% genetic similarity
cluster <- DECIPHER::Treeline(myDistMatrix = dist,
                              method = "complete",
                              cutoff = 0.03,
                              type = "clusters",
                              processors = nproc)


bac_data <- merge_taxa_vec(bac_data_raw, group = cluster$cluster, tax_adjust = 2)
bac_data


######## Cleaning up ########
## Removing unwanted sequences from Eukarya, mitochondria, chloroplast, as well as high-level unknowns

table(tax_table(bac_data)[, "Kingdom"], exclude = NULL) # checking groups

bac_data_king <- subset_taxa(bac_data, !is.na(Kingdom) & !Kingdom %in% c("Eukaryota"))

bac_data_phy <- subset_taxa(bac_data_king, !is.na(Phylum) & !Phylum %in% c("Deinoccota","FCPU426","Elusimicrobiota","MBNT15",
                                                                      "Methylomirabilota"))

# Checking the percentage of reads 
sample_sums(bac_data_phy)
(sample_sums(bac_data_phy)/sample_sums(bac_data))
phyclean <- (1-(sample_sums(bac_data_phy)/sample_sums(bac_data)))*100

bac_data_chl <- subset_taxa(bac_data_phy, !is.na(Order) & !Order %in% c("Chloroplast"))

# Checking the percentage of reads
sample_sums(bac_data_chl)
(sample_sums(bac_data_chl)/sample_sums(bac_data))
chlclean <- (1-(sample_sums(bac_data_chl)/sample_sums(bac_data)))*100

bac_data_mit <- subset_taxa(bac_data_chl, !is.na(Family) & !Family %in% c("Mitochondria"))

# Checking the percentage of reads
sample_sums(bac_data_mit)
(sample_sums(bac_data_mit)/sample_sums(bac_data))
mitclean <- (1-(sample_sums(bac_data_mit)/sample_sums(bac_data)))*100


### Track ASV loss
reads <- cbind(phyclean, chlclean, mitclean)
head(reads)

ASVlost <- c(ntaxa(bac_data), ntaxa(bac_data_king), ntaxa(bac_data_phy), ntaxa(bac_data_chl), ntaxa(bac_data_mit))
Readslost <- c(sum(sample_sums(bac_data)), sum(sample_sums(bac_data_king)), sum(sample_sums(bac_data_phy)), sum(sample_sums(bac_data_chl)), 
               sum(sample_sums(bac_data_mit)))

alltrack <- data.frame("TrackASVs" = ASVlost, "TrackReads" = Readslost)
rownames(alltrack) <- c("Raw","Kingdom","Phyla","Chloroplasts","Mitrochondria")
write.csv(alltrack, "bacterial_ps_clearning_trackReads.csv")


######## Normalizing the phyloseq object to account for differences in sequencing depth ########
bac_norm <- norm_tmm(otu_table(bac_data_mit),
                     ref_column = NULL,
                     logratio_trim = 0.3,
                     sum_trim = 0.05,
                     do_weighting = TRUE,
                     Acutoff = -1e+10)

bac_norm
write.csv(bac_norm, "normalized_OTU_table.csv")

bacteria <- phyloseq(otu_table(bac_norm), taxa_are_rows = TRUE, sample_data(sampdat), tax_table(bactaxa))

# Add library size to the data frame
sample_data(bacteria)$libsize <- sample_sums(bacteria)

# Simplify seq names from raw sequences
dna <- Biostrings::DNAStringSet(taxa_names(bacteria))
names(dna) <- taxa_names(bacteria)

bacteria <- merge_phyloseq(bacteria, dna)
taxa_names(bacteria) <- paste0("ASV", seq(ntaxa(bacteria)))
names(dna) <- taxa_names(bacteria)

writeXStringSet(dna, filepath = "allbackASVs4blast.fasta", format = "fasta")

# Save phyloseq object for all downstream analysis
saveRDS(bacteria, "bacteria_final_ps.rds")


melt <- psmelt(bacteria)
write.csv(melt, "Final_bacterial_phyloseq_object.csv")






