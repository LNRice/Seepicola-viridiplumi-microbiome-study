###### Microbial 16S rRNA analysis  ########
##### Final bacterial analysis for Seepicola viridiplumi microbiome manuscript

## Rice LN, Louca, S, Cordes EE, Young CM. (2026). Trophic niches and the associated
## microbial community for a newly described methane seep sabellid. Deep-Sea Research
## Part 1

## This code is used to analyze the microbial 16s rRNA data on Seepicola viridiplumi radioles and anterior tissues
## for both free-living and epibiotic individuals. It also includes functional diversity assignments, 
## statistical analyses, and figure creation. 


####### Set Up #######

## Set Working directory:
#setwd() #CHANGE ME to desired directory

#Load packages:
library(phyloseq)
library(Biostrings)
library(DECIPHER)
library(microbiomeMarker)
library(vegan)
library(ggplot2)
library(Tax4Fun)
library(RColorBrewer)
library(microeco)
library(file2meco)
library(ranacapa)
library(car)
library(plyr)
library(pairwiseAdonis)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
######## Open the needed data files ########

## The phyloseq object
bac <- readRDS("bacteria_final_ps.rds")

## Extract sample data frame for easier manipulation
sampdat <- data.frame(sample_data(bac))

sample_names(bac) <- sampdat$Abbreviated.Sample.Name

## Add in needed databases; CHANGE ME to needed filepath
taxafun <- "/taxa" 
blast_bin <- "/taxablast" 


######## Examining microbial functional diversity ########

## Transform phyloseq object to a microtable to be used in microeco package
micro <- phyloseq2meco(bac)

transmicro <- trans_func$new(micro)               # Data set format
transmicro$cal_func(prok_database = "FAPROTAX")   # Maps taxonomy to FAPROTAX data frame
transmicro$res_spe_func[1:50, 1:4]                # Shows presence/absence of a functional trait in a sample
transmicro$cal_func_FR(abundance_weighted = TRUE) # Calculates percentages for each community
transmicro$res_spe_func_perc[1:20, 1:10]          # Shows percentages for each functional trait by sample
transmicro$trans_func_FR()                        # Adjusts above output into long table format 

tiff("bacterial_functional_diversity.tiff", width = 1600, height = 1200)
transmicro$plot_func_FR()
dev.off()

transmicro$cal_tax4fun2(blast_tool_path = NULL, path_to_reference_data = "Tax4Fun2_ReferenceData_v2", 
                        path_to_temp_folder = NULL)

# Make data frame of taxonomic abundances
func <- as.data.frame(transmicro$res_func_FR_trans)


######## Examining Microbial Diversity by Taxonomic Family ########
name20 <- names(sort(taxa_sums(bac), decreasing = TRUE))[1:20]
top20 <- transform_sample_counts(bac, function(otu)otu/sum(otu))
top20 <- prune_taxa(name20,top20)
melt20 <- psmelt(top20)  # Enables viewing of phyloseq object as a .csv file
col_pal <- c("#3C3933","#028A75","#840303","#A59A63","#387B2D","#6D8AA0","#905508","#37268D","#8298F7","#97A30B","#6B15BE","#C02F20","antiquewhite",
             "azure4","#EC51E3","#332B2F","#25D90B","#CE8C12","deepskyblue4","#122001","maroon","orange","blue","darkcyan","tan3")

tiff("bacterial_top20.tiff", width = 700, height = 700)

sitetop20 <- ggplot(melt20, aes(Sample, Abundance, fill=Family)) +
  geom_bar(stat = "identity", position = "fill") + 
  facet_grid(~Site, scales = "free_x", space = "free") +
  labs(x = "Sample", y = "Relative Abundance (%)") +
  theme(axis.text.x = element_text(angle = 90), 
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 0.5),
        panel.background = element_blank()) +
  scale_fill_manual(values = col_pal) +
  theme(axis.text.x = element_text(size = 10, color = "black", angle = -90),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14, colour = "black"))
sitetop20

dev.off()

## Percentages by Family 
bac_fam <- tax_glom(bac, "Family")
rel_abund <- transform_sample_counts(bac_fam, function(x) 100 * x / sum(x))
fam <- psmelt(rel_abund)
write.csv(fam, file = "Relative_abundance_by_sample.csv")

# Percentages by family overall
summary <- fam %>% 
  dplyr::group_by(Family) %>%
  dplyr::summarize(Mean_Percent = mean(Abundance, na.rm = TRUE)) %>%
  dplyr::arrange(desc(Mean_Percent))


## Plotting Rarefaction curves, color by site, line type by microhabitat
cmin <- min(sample_sums(bac))
otutab <- otu_table(bac)
class(otutab) <- "matrix"

cols <- c("Bush Hill" = "#7e03a8", "Brine Pool NR1" = "#f89540", "Green Canyon 234" = "#122001")

tiff("Rarecurve_final.tiff", width = 2000, height = 2000, res = 300)

rarecurve(otutab, step = 100, col = cols, lwd=2, label = T, ylab = "ASV Richness",
          lty = c(2,2,2,2,2,2,1,2,2,2,1,1,1,1,2,2,2,2,2,2))
abline(v = cmin)
legend("bottomright", legend = c("Free-Living","Epibiotic"), col = "black", lty = c(2,1), bty = "n")
legend("topright", legend = c("Bush Hill", "Brine Pool NR1","Green Canyon 234"), 
       col = c("#7e03a8","#f89540","#122001"), lty = 1, bty = "n", lwd = 2)

dev.off()


######## Alpha Diversity ########
## Plotting alpha diversity measures by site and microhabitat condition
tiff("Alph_divers_bySiteandCondition.tiff", width = 700, height = 700)
p1 <- plot_richness(bac, x = "Site", measures = c("Shannon", "Observed"), color = "Condition") +
  geom_boxplot() + scale_color_manual(values = c("grey20","red3")) +
  theme(panel.border = element_rect(fill = "transparent", color = "black", linewidth = 0.5)) +
  theme(panel.background = element_blank(),
        strip.background = element_blank(),
        axis.line.y = element_line(color = "black", linewidth = 0.5),
        axis.line.x = element_line(color = "black", linewidth = 0.5),
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 0.5))
  
p1
dev.off()

## Isolate Brine Pool NR1 free-living samples to avoid conflating microhabitat condition across sites
free <- subset_samples(bac, !Condition %in% c("Epibiotic"))

freemelt <- psmelt(free) # done to check previous step worked

### Compare diversity measures statistically
richness <- estimate_richness(free, split = TRUE, measures = c("Shannon","Observed"))

# Adding groups with sample data to the richness data frame
richness$site <- c("Bush Hill","Bush Hill","Bush Hill","Bush Hill","Bush Hill","Bush Hill","Brine Pool","Brine Pool","Brine Pool",
                   "Green Canyon 234","Green Canyon 234","Green Canyon 234","Green Canyon 234","Green Canyon 234","Green Canyon 234")

# Checking assumptions
shapiro.test(richness$Observed) 
shapiro.test(richness$Shannon)  

leveneTest(richness$Observed~richness$site) 
leveneTest(richness$Shannon~richness$site)  

# Comparing diversity metrics across site
summary(aov(richness$Shannon~richness$site))
summary(aov(richness$Observed~richness$site)) 


## Comparing diversity metrics between microhabitats within Brine Pool NR1
bp <- subset_samples(bac, !Site %in% c("Green Canyon 234","Bush Hill"))
bpmelt <- psmelt(bp)

bprich <- estimate_richness(bp, split = TRUE, measures = c("Shannon","Observed"))
bprich$condition <- c("Epibiotic","Free-Living","Free-Living","Free-Living","Epibiotic","Epibiotic","Epibiotic","Epibiotic")

# Checking Assumptions
shapiro.test(bprich$Observed) 
shapiro.test(bprich$Shannon)  

leveneTest(bprich$Observed~bprich$condition)
leveneTest(bprich$Shannon~bprich$condition)  

#Comparing diversity metrics within Brine Pool 
t.test(bprich$Observed~bprich$condition)  
t.test(bprich$Shannon~bprich$condition)  


######## Beta Diverity across sites - PCoA Bray ########
free <- transform_sample_counts(free, function(otu) otu/sum(otu))
ord.pcoa.bray <- ordinate(free, method = "PCoA", distance = "bray")

# Look at plot and get axefree# Look at plot and get axes % variation explained - Axis.1 = 24.8%, Axis.2 = 19.9%
plot_ordination(free, ord.pcoa.bray)
ord.pcoa.bray$values$Eigenvalues

#PCoA get scores
PCoAxes <- ord.pcoa.bray$vectors[,1:2]

# Add in ancillary data
scoresdf <- data.frame(PCoAxes)
scoresdf$ID <- c("BHSa118","BHSa120","BHSa121","BHSa85","BHSa87","BHSa93","BPSa137","BPSa152","BPSa154",
                 "GCSa59","GCSa63","GCSa64","GCSa67","GCSa71","GCSa72")
scoresdf$Site <- c("Bush Hill","Bush Hill","Bush Hill","Bush Hill","Bush Hill","Bush Hill","Brine Pool NR1",
                   "Brine Pool NR1","Brine Pool NR1","Green Canyon 234","Green Canyon 234","Green Canyon 234",
                   "Green Canyon 234","Green Canyon 234","Green Canyon 234")
head(scoresdf)

# Hulls by site
hulls <- ddply(scoresdf, .(Site), function(scoresdf) 
  scoresdf[chull(scoresdf$Axis.1, scoresdf$Axis.2), ])

tiff("PCoA_16sbySite.tiff", width = 700, height = 700)
site <- ggplot(scoresdf, aes(x = Axis.1, y = Axis.2, col = Site)) +
  geom_point(size = 2, aes(colour = Site)) +
  geom_polygon(data = hulls, aes( alpha = 0.1, fill = Site)) + 
  labs(x = "Axis 1 [24.8%]", y = "Axis 2 [19.9%]") +
  scale_fill_manual(values = c("#700875","#0580BB","#1C480F")) +
  scale_color_manual(values = c("#700875","#0580BB","#1C480F")) +
  theme(panel.background = element_blank(), #strip.background = element_blank(),
        axis.line.y = element_line(colour = "black", linewidth = 0.5),
        axis.line.x = element_line(colour = "black", linewidth = 0.5),
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 0.5),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14, colour = "black"))
site

dev.off()

######## Checking for differences in community composition ########
bray_dis <- phyloseq::distance(free, method = "bray")

res_adonis <- adonis2(bray_dis ~ Site, scoresdf)
print(res_adonis)


data.frame(name = rownames(res_adonis), DF = res_adonis$Df,
           SumofSquares = res_adonis$SumOfSqs, R2 = res_adonis$R2,
           f.model = res_adonis$F, p = res_adonis$`Pr(>F)`) %>% write.csv(file = "Site_PERMANOVA_braydis_full.csv")

pairwise.adonis2(bray_dis ~ Site, data = scoresdf)


#check for differences in multivariuate dispersion based on site
beta <- betadisper(bray_dis, scoresdf$Site) # no sig heterogeneity of dispersion
permutest(beta) 


## Redo above stats on data from Brine Pool NR1 to check influence of epibiotic condition
bp <- transform_sample_counts(bp, function(otu) otu/sum(otu))
bp.pcoa.bray <- ordinate(bp, method = "PCoA", distance = "bray")
  
plot_ordination(bp, bp.pcoa.bray)
bp.pcoa.bray$values$Eigenvalues

#PCoA get scores
bpPCoAxes <- bp.pcoa.bray$vectors[,1:2]

# Add in ancillary data
bpdf <- data.frame(bpPCoAxes)
bpdf$ID <- c("BPSa109","BPSa137","BPSa152","BPSa154","BPSa91","BPSa93","BPSa96","BPSa97")
bpdf$Condition <- c("Epibiotic","Free-Living","Free-Living","Free-Living","Epibiotic","Epibiotic","Epibiotic","Epibiotic")
bpdf$libsize <- c("3532","8176","170","1394","943","932","3894","1716")
head(bpdf)

# Hulls by site
bphulls <- ddply(bpdf, .(Condition), function(bpdf) 
  bpdf[chull(bpdf$Axis.1, bpdf$Axis.2), ])

tiff("PCoA_16sbyCondition.tiff", width = 700, height = 700)
condbp <- ggplot(bpdf, aes(x = Axis.1, y = Axis.2, col = Condition)) +
  geom_point(size = 2, aes(colour = Condition)) +
  geom_polygon(data = bphulls, aes( alpha = 0.2, fill = Condition)) + 
  labs(x = "Axis 1 [36.3%]", y = "Axis 2 [26.8%]") +
  scale_fill_manual(values = c("grey20","red3")) +
  scale_color_manual(values = c("grey20","red3")) +
  theme(panel.background = element_blank(), #strip.background = element_blank(),
        axis.line.y = element_line(colour = "black", linewidth = 0.5),
        axis.line.x = element_line(colour = "black", linewidth = 0.5),
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 0.5),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14, colour = "black"))
condbp

dev.off()

######## Checking for differences in community composition ########
bpbray_dis <- phyloseq::distance(bp, method = "bray")
bpres_adonis <- adonis2(bpbray_dis~libsize, data = bpdf)

bpres_adonis <- adonis2(bpbray_dis ~ Condition, data = bpdf)
print(bpres_adonis)


data.frame(name = rownames(bpres_adonis), DF = bpres_adonis$Df,
           SumofSquares = bpres_adonis$SumOfSqs, R2 = bpres_adonis$R2,
           f.model = bpres_adonis$F, p = bpres_adonis$`Pr(>F)`) %>% write.csv(file = "condition_PERMANOVA_braydis_full.csv")


bpbeta <- betadisper(bpbray_dis, bpdf$Condition) 
permutest(bpbeta) 

