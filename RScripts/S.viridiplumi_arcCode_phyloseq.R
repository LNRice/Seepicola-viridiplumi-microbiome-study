###### Microbial 16S rRNA post-DADA2 workflow, creating Phyloseq and data cleaning ########
##### Final analysis for Seepicola viridiplumi microbiome manuscript

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
library(vegan)
library(ranacapa)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
######## Creating Phyloseq Object ########

## Open the needed RDS files
arcseqtab <- readRDS("archaeal_seqtab.nonchim_merged_11-10.rds")
arctaxa <- readRDS("Archaeal_taxa_merged_11-10.rds")


# bringing in sample metadata
sampdat <- read.table(file = "samplemetadata.csv", sep = ",", header = TRUE, row.names = 1)
rownames(sampdat) <- sampdat$Abbreviated.Sample.Name
rownames(arcseqtab) <- sampdat$Abbreviated.Sample.Name

# Making phyloseq object
arc_data_raw <- phyloseq(otu_table(arcseqtab, taxa_are_rows = FALSE), sample_data(sampdat), tax_table(arctaxa))
arc_data_raw

sample_sums(arc_data_raw)
table(tax_table(arc_data_raw)[, "Kingdom"], exclude = NULL)

# Isolating known Archaeal ASVs
is.archaea <- !is.na(arctaxa[, "Kingdom"]) & arctaxa[, "Kingdom"] == "Archaea"
seqtab.archaea <- arcseqtab[, is.archaea, drop = FALSE]
taxa.archaea <- arctaxa[is.archaea, drop = FALSE]
rowSums(seqtab.archaea)

# Make phyloseq obj
arc_data_filt <- phyloseq(otu_table(seqtab.archaea, taxa_are_rows = FALSE),
                         sample_data(sampdat),
                         tax_table(arctaxa))
arc_data_filt
sample_sums(arc_data_filt)


### Identify similar ASVs and cluster them 
asv_sequences <- colnames(seqtab.archaea)
names <- rownames(seqtab.archaea)
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

arc_data <- merge_taxa_vec(arc_data_filt, group = cluster$cluster, tax_adjust = 2)
arc_data


sample_sums(arc_data) 
sum(sample_sums(arc_data))  
asvs_per_sample <- sample_sums(otu_table(arc_data)>0)
print(asvs_per_sample)
percent_retained <- ((sample_sums(arc_data))/(sample_sums(arc_data_raw)))*100
print(percent_retained)
sum(asvs_per_sample)


########## Making rarefaction curves ########
cmin <- min(sample_sums(arc_data))
otutab <- otu_table(arc_data)
class(otutab) <- "matrix"

col_pal <- c("#7e03a8","#f89540","#122001")


tiff("Rarecurve_archaea_final.tiff", width = 2000, height = 2000, res = 300)

rarecurve(otutab, step = 100, col = col_pal, lwd=2, label = T, ylab = "ASV Richness",
          lty = c(2,2,2,2,2,2,1,2,2,2,1,1,1,1,2,2,2,2,2,2))
abline(v = cmin)
legend("bottomright", legend = c("Free-Living","Epibiotic"), col = "black", lty = c(2,1), bty = "n")
legend(x = 230, y = 5, legend = c("Bush Hill", "Brine Pool NR1","Green Canyon 234"), 
       col = c("#7e03a8","#f89540","#122001"), lty = 1, bty = "n", lwd = 2)

dev.off()


######## Organizing data ########
### As some of my samples obtained 0 Archaeal reads after initial merged reads of 1,100 - 1,600 (BPSa91 and BPSa93), this likely indicates Archaea 
## were genuinely absent from these samples

### Will avoid rarefying the data to prevent loosing low read samples

# Add library size to the data frame
sample_data(arc_data)$libsize <- sample_sums(arc_data)

# Simplify seq names from raw sequences
dna <- Biostrings::DNAStringSet(taxa_names(arc_data))
names(dna) <- taxa_names(arc_data)

archaea <- merge_phyloseq(dna, arc_data)
taxa_names(archaea) <- paste0("ASV", seq(ntaxa(archaea)))
names(dna) <- taxa_names(archaea)

writeXStringSet(dna, filepath = "allArcASVs4blast.fasta", format = "fasta")

# Save phyloseq object for all downstream analysis
saveRDS(archaea, "archaea_final_ps_11-11.rds")


melt <- psmelt(archaea)
write.csv(melt, "Final_archaeal_phyloseq_object_11-11.csv")

