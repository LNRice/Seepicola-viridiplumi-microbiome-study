###### Microbial 16S rRNA analysis  ########
##### Final Archaeal analysis for Seepicola viridiplumi microbiome manuscript

## Rice LN, Louca, S, Cordes EE, Young CM. (2026). Trophic niches and the associated
## microbial community for a newly described methane seep sabellid. Deep-Sea Research
## Part 1

## This code is used to analyze the microbial 16s rRNA data on Seepicola viridiplumi radioles and anterior tissues
## for both free-living and epibiotic individuals. It also includes functional diversity assignments, 
## statistical analyses, and figure creation. 


####### Set Up #######

## Set Working directory:
#setwd() #CHANGE ME to desired directory

#Load Packages
library(phyloseq)
library(Biostrings)
library(DECIPHER)
library(microbiomeMarker)
library(vegan)
library(ggplot2)
library(Tax4Fun)
library(file2meco)
library(microeco)
library(car)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
######## Open the needed data files ########

## Load the phyloseq object
archaea <- readRDS("archaea_final_ps_11-11.rds")
samp <- data.frame(sample_data(archaea))
sample_names(archaea) <- samp$Abbreviated.Sample.Name

## Add in the taxonomic database file path
taxafun <- "/taxa" 
blast_bin <- "/taxablast"


######## Examining functional diversity of the Archaeal community ########
## Transform phyloseq object to a microtable format
microarc <- phyloseq2meco(archaea)

transmicroarc <- trans_func$new(microarc)                 # Sets dataframe format
transmicroarc$cal_func(prok_database = "FAPROTAX")        # Maps taxonomy to FAPROTAX dataframe
transmicroarc$res_spe_func[1:50, 1:4]                     # Shows presence/absence of a functional trait in a sample
transmicroarc$cal_func_FR(abundance_weighted = TRUE)      # Calculates percentages for each community
transmicroarc$res_spe_func_perc[1:20, 1:10]               # Shows percentages for each functional trait by sample
transmicroarc$trans_func_FR()                             # Adjusts output into long table format

tiff("Archaeal_functional_diversity.tiff", width = 700, height = 700)
transmicroarc$plot_func_FR()
dev.off()

transmicroarc$cal_tax4fun2(blast_tool_path = NULL, path_to_reference_data = "Tax4Fun2_ReferenceData_v2", path_to_temp_folder = NULL)

# Make dataframe of taxonomic abundances
arcfunc <- as.data.frame(transmicroarc$res_func_FR_trans)


######## Examining Archaeal Diversity ########
name10 <- names(sort(taxa_sums(archaea), decreasing = TRUE))[1:10]
top10 <- transform_sample_counts(archaea, function(otu)otu/sum(otu))
top10 <- prune_taxa(name10, top10)
melt <- psmelt(top10)


arc_pal <- c("#3C3933","#028A75","#840303","#A59A63","#387B2D","#6D8AA0","#905508","#37268D","#8298F7","#97A30B","#6B15BE","#C02F20","antiquewhite",
             "azure4","#EC51E3","#332B2F","#25D90B","#CE8C12","deepskyblue4","#122001","maroon","orange","blue","darkcyan","antiquewhite3")

tiff("Archaeal_top10.tiff", width = 700, height = 700)

sitetop10 <- ggplot(melt, aes(Sample, Abundance, fill=Family)) +
  geom_bar(stat = "identity", position = "fill") + 
  facet_grid(~Site, scales = "free_x", space = "free") +
  labs(x = "Sample", y = "Relative Abundance (%)") +
  theme(axis.text.x = element_text(angle = 90), 
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 0.5),
        panel.background = element_blank()) +
  scale_fill_manual(values = arc_pal) +
  theme(axis.text.x = element_text(size = 10, color = "black", angle = -90),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14, colour = "black"))
sitetop10

dev.off()


######## Alpha Diversity ########
## Plotting alpha diversity measures by site and microhabitat condition

tiff("Alph_divers_bySiteandCondition.tiff", width = 700, height = 700)
p1 <- suppressWarnings(
  plot_richness(archaea, x = "Site", measures = c("Shannon", "Observed"), color = "Condition") +
  geom_boxplot() + scale_color_manual(values = c("grey20","red3")) +
  theme(panel.border = element_rect(fill = "transparent", color = "black", linewidth = 0.5)) +
  theme(panel.background = element_blank(),
        strip.background = element_blank(),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 0.5)))

p1
dev.off()


## Isolate Brine Pool NR1 free-living samples to avoid conflating microhabitat condition across sites
free <- subset_samples(archaea, !Condition %in% c("Epibiotic"))

freemelt <- psmelt(free) # done to check previous step worked

### Compare diversity measures statistically
richness <- suppressWarnings(estimate_richness(free, split = TRUE, measures = c("Shannon","Observed")))

# Adding groups with sample data to the richness data frame
richness$site <- c("Bush Hill","Bush Hill","Bush Hill","Bush Hill","Bush Hill","Bush Hill","Brine Pool","Brine Pool","Brine Pool",
                   "Green Canyon 234","Green Canyon 234","Green Canyon 234","Green Canyon 234","Green Canyon 234","Green Canyon 234")


# Checking assumptions
shapiro.test(richness$Observed)  
shapiro.test(richness$Shannon) 

leveneTest(richness$Observed~richness$site) 
leveneTest(richness$Shannon~richness$site)  


# Comparing diversity metrics across site
kruskal.test(richness$Shannon~richness$site)
kruskal.test(richness$Observed~richness$site) 



## Comparing diversity metrics between microhabitats within Brine Pool NR1
bp <- subset_samples(archaea, !Site %in% c("Green Canyon 234","Bush Hill"))
bpmelt <- psmelt(bp)

bprich <- suppressWarnings(estimate_richness(bp, split = TRUE, measures = c("Shannon","Observed")))
bprich$condition <- c("Epibiotic","Free-Living","Free-Living","Free-Living","Epibiotic","Epibiotic","Epibiotic","Epibiotic")

# Checking Assumptions
shapiro.test(bprich$Observed)  

leveneTest(bprich$Observed~bprich$condition)   


#Comparing diversity metrics within Brine Pool 
wilcox.test(bprich$Observed~bprich$condition) 



######## Not enough ASV abundance to examine Beta Diversity trends ########








