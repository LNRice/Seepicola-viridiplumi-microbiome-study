###### Carbon, Nitrogen, and Sulfur Stable Isotope Analysis  ########
##### Statistical analysis and data visualization for Seepicola viridiplumi microbiome manuscript

## Rice LN, Louca, S, Cordes EE, Young CM. (2026). Trophic niches and the associated
## microbial community for a newly described methane seep sabellid. Deep-Sea Research
## Part 1

## This code is used to analyze the stable isotopic signatures from Seepicola viridiplumi radioles and anterior tissues
## for both free-living and epibiotic individuals. It also includes SIBER analyses, 
## statistical analyses, and figure creation. 


######## Set Up ########
## Set Working directory:
#setwd() #CHANGE ME to desired directory

## Load packages
library(readxl)
library(car) 
library(FSA)
library(ggplot2)
library(dplyr)
library(plyr)
library(ggforce)
library(rcompanion)
library(datawizard)
library(glmmTMB)
library(emmeans)
library(DHARMa)
library(tidyr)


## Loading the Needed Data Sheets:
#CHANGE Me to the directory containing the Excel datasheets
alldata <- "/data" 
carbon <- "/carb" 
nitrogen <- "/nit"
sulfur <- "/sulf" 
ave <- "/averages"
avesulf <- "/sulfuraverages"

######## TESTING ASSUMPTIONS ########

## Examining skewness and variances in the isotopic datasets
skewness(carbon$d13C) 
plot(density(carbon$d13C))
leveneTest(carbon$d13C ~ carbon$Collection * carbon$Site) 

skewness(nitrogen$d15N) 
plot(density(nitrogen$d15N))
leveneTest(nitrogen$d15N ~ nitrogen$Collection * nitrogen$Site) 

skewness(sulfur$d34S)
plot(density(sulfur$d34S))
leveneTest(sulfur$d34S ~ sulfur$Collection * sulfur$Site) 

## Using GLMMs to examine each isotope by site, collection, and epibiotic condition

#~~~~~~~~~~~~~~~~~~~~~~~~
######## ISOTOPE 1: Carbon ########
#~~~~~~~~~~~~~~~~~~~~~~~~
##### Model 1 - Using free-living data to examine effect of site and collection

# Reflect carbon data to account for negative values (make it positively right skewed)
c_buffer <- max(carbon$d13C, na.rm = TRUE) +1
carbon$reflected <- c_buffer-carbon$d13C

# Making a sampling batch identifier
carbon$sampling_batch <- paste(carbon$Site, carbon$Collection, sep = "_")

# Subset data to isolate free-living individuals and ensuring all rows with NAs are removed
c_free <- carbon %>%
  filter(Condition == "Free-living") %>%
  drop_na(d13C, Site, Collection, sampling_batch) %>%
  droplevels()

# Fitting the model to allow for dispersion (variance) to change by site
c_model <- glmmTMB(reflected ~ Site + Collection + (1|sampling_batch),
                   dispformula = ~ Site,
                   data = c_free,       
                   family = Gamma(link = "log"))

# Simulating the residuals focusing on the fixed effects structure
c_res <- simulateResiduals(fittedModel = c_model, n = 1000, use.u = FALSE)

# Running a categorical diagnostic plot
plotResiduals(c_res, form = c_free$Site)
testDispersion(c_res)
summary(c_model)

Anova(c_model, type = "III")

## Extract estimated means from the site-specific dispersion model using Satterthwaite degrees of freedom (to accommodate smaller sample sizes)
c_emm <- emmeans(c_model, ~ Site + Collection, ddf = "satterthwaite", type ="response")
c_df <- as.data.frame(c_emm)
  
cpair <- contrast(c_emm, method = "pairwise", adjust = "tukey")
print(cpair)
write.csv(cpair, file = "carbon pairwise comparisions_11-11.csv", row.names = TRUE)


##### Model 2 - Examining effect of epibiosis 
# Identify sites and collections when both microhabitats were sampled
c_counts <- table(carbon$sampling_batch, carbon$Condition)
c_cobatch <- rownames(c_counts)[rowSums(c_counts > 0) > 1]
print(c_cobatch)

c_batch <- carbon %>%
  dplyr::filter(sampling_batch %in% c_cobatch) %>%
  tidyr::drop_na(d13C, Condition, Site, sampling_batch) %>%
  droplevels()

# Run the model
c_epimodel <- glmmTMB(reflected ~ Condition + Site + (1 | sampling_batch),
                      dispformula = ~Site,
                      data = c_batch,
                      family = Gamma(link = "log"))

# Simulating the residuals 
c_epires <- simulateResiduals(fittedModel = c_epimodel, n = 1000, use.u = FALSE)

# Running a categorical diagnostic plot
plotResiduals(c_epires, form = c_batch$Condition)
testDispersion(c_epires)
summary(c_epimodel)

Anova(c_epimodel, type = "III")

## Extract the marginal means and run the pairwise comparisons
c_epiemm <- emmeans(c_epimodel, ~ Condition + Site, ddf = "satterthwaite", type = "response")
print(as.data.frame(c_epiemm))

cepipair <- contrast(c_epiemm, method = "pairwise", adjust = "tukey")
write.csv(cepipair, file = "carbon pairwise epibiotic comparisions_11-11.csv", row.names = TRUE)


#~~~~~~~~~~~~~~~~~~~~~~~~
######## ISOTOPE 2: Nitrogen ########
#~~~~~~~~~~~~~~~~~~~~~~~~
##### Model 1 - Using free-living data to examine effect of site and collection

# Making a sampling batch identifier
nitrogen$sampling_batch <- paste(nitrogen$Site, nitrogen$Collection, sep = "_")

# Subset data to isolate free-living individuals and ensuring all rows with NAs are removed
n_free <- nitrogen %>%
  filter(Condition == "Free-living") %>%
  drop_na(d15N, Site, Collection, sampling_batch) %>%
  droplevels()

# Data values are all positive and normally distributed. No adjustments needed
n_model <- glmmTMB(d15N ~ Site + Collection + (1|sampling_batch), 
                   dispformula = ~ Site,
                   data = nitrogen,
                   family = Gamma(link = "log"))

# Simulate the residuals
n_res <- simulateResiduals(fittedModel = n_model, n = 1000, re.form = NULL)

# Testing for Dispersion
testDispersion(n_res)
plot(n_res)
summary(n_model)

Anova(n_model, type = "III")

## Extract estimated means from the dispersion model, ignoring site level comparisons due to the non-significant site effect
n_emm <- emmeans(n_model, ~ Site + Collection, ddf = "satterthwaite", type = "response")
n_emm_df <- as.data.frame(n_emm)

## Define scaling parameters
n_min <- min(nitrogen$d15N, na.rm = TRUE) - 0.1
n_max <- max(nitrogen$d15N, na.rm = TRUE) + 0.1
n_range <- n_max - n_min

# Looking at the pairwise comparisons from the model
npair <- contrast(n_emm, method = "pairwise", adjust = "tukey")
print(npair)
write.csv(npair, file = "nitrogen pairwise comparisons_11-11.csv", row.names = TRUE)


##### Model 2 - Examining effect of epibiosis 

# Identify sites and collections when both microhabitats were sampled
n_counts <- table(nitrogen$sampling_batch, nitrogen$Condition)
n_cobatch <- rownames(n_counts)[rowSums(n_counts > 0) > 1]
print(n_cobatch)

n_batch <- nitrogen %>%
  dplyr::filter(sampling_batch %in% n_cobatch) %>%
  tidyr::drop_na(d15N, Condition, Site, sampling_batch) %>%
  droplevels()

# Run the model
n_epimodel <- glmmTMB(d15N ~ Condition + Site + (1 | sampling_batch),
                      dispformula = ~Site,
                      data = n_batch,
                      family = Gamma(link = "log"))

# Simulating the residuals 
n_epires <- simulateResiduals(fittedModel = n_epimodel, n = 1000, use.u = FALSE)

# Running a categorical diagnostic plot
plotResiduals(n_epires, form = n_batch$Condition)
testDispersion(n_epires)
summary(n_epimodel)


Anova(n_epimodel, type = "III")

## Extract the marginal means and run the pairwise comparisons
n_epiemm <- emmeans(n_epimodel, ~ Condition + Site, ddf = "satterthwaite", type = "response")
print(as.data.frame(n_epiemm))

nepipair <- contrast(n_epiemm, method = "pairwise", adjust = "tukey")
write.csv(nepipair, file = "nitrogen pairwise epibiotic comparisions_11-11.csv", row.names = TRUE)


#~~~~~~~~~~~~~~~~~~~~~~~~
######## ISOTOPE 3: Sulfur ########
#~~~~~~~~~~~~~~~~~~~~~~~~
##### Model 1 - Using free-living data to examine effect of site and collection

# Min-max normalizing the Sulfur isotopes to accommodate positive and negative data spread
s_min <- min(sulfur$d34S, na.rm= TRUE) - 0.1
s_max <- max(sulfur$d34S, na.rm = TRUE) + 0.1
s_range <- (s_max - s_min)

sulfur$scaled <- (sulfur$d34S - s_min) / (s_max - s_min)

# Making a sampling batch identifier
sulfur$sampling_batch <- paste(sulfur$Site, sulfur$Collection, sep = "_")

# Subsetting data to isolate free-living individuals and ensuring all rows with NAs are removed
s_free <- sulfur %>%
  filter(Condition == "Free-living") %>%
  drop_na(d34S, Site, Collection, sampling_batch) %>%
  droplevels()

# Running the data through a Beta fixed model
s_model <- glmmTMB(scaled ~ Site + Collection + (1|sampling_batch),
                   dispformula = ~ Site,
                   data = sulfur,
                   family = beta_family(link = "logit"))
summary(s_model)


Anova(s_model, type ="III")

# Simulate the residuals
s_res <- simulateResiduals(fittedModel = s_model, n = 1000, re.form = NULL)

# Testing for Dispersion
testDispersion(s_res)
plot(s_res)

# Extract estimated means from my final model
s_emm <- emmeans(s_model, ~ Site + Collection, type = "response")
s_emm_df <- as.data.frame(s_emm) 
  
# Reverse the min-max normalization
s_emm_df <- s_emm_df %>%
  mutate(mean_s = (response * s_range) + s_min,
         lower_s = (asymp.LCL * s_range) + s_min,
         upper_s = (asymp.UCL * s_range) + s_min) %>%
  select(Site, Collection, mean_s, lower_s, upper_s)

# Looking at the pairwise comparisons from the model
spair <- contrast(s_emm, method = "pairwise", adjust = "tukey")
print(spair)
write.csv(spair, file = "sulfur pairwise comparisons_11-11.csv", row.names = TRUE, col.names = TRUE)


##### Model 2 - Examining effect of epibiosis 

# Identify sites and collections when both microhabitats were sampled
s_counts <- table(sulfur$sampling_batch, sulfur$Condition)
s_cobatch <- rownames(s_counts)[rowSums(s_counts > 0) > 1]
print(s_cobatch)

s_batch <- sulfur %>%
  dplyr::filter(sampling_batch %in% s_cobatch) %>%
  tidyr::drop_na(d34S, Condition, Site, sampling_batch) %>%
  droplevels()

# Run the model
s_epimodel <- glmmTMB(scaled ~ Condition + Site + (1 | sampling_batch),
                      dispformula = ~Site,
                      data = s_batch,
                      family = beta_family(link = "logit"))

# Simulating the residuals 
s_epires <- simulateResiduals(fittedModel = s_epimodel, n = 1000, use.u = FALSE)

# Running a categorical diagnostic plot
plotResiduals(s_epires, form = s_batch$Condition)
testDispersion(s_epires)
summary(s_epimodel)

Anova(s_epimodel, type = "III")

## Extract the marginal means and run the pairwise comparisons
s_epiemm <- emmeans(s_epimodel, ~ Condition + Site, ddf = "satterthwaite", type = "response")
print(as.data.frame(s_epiemm))

sepipair <- contrast(s_epiemm, method = "pairwise", adjust = "tukey")
write.csv(sepipair, file = "sulfur pairwise epibiotic comparisions_11-11.csv", row.names = TRUE)


#~~~~~~~~~~~~~~~~~~~~~~~~~
######## Beginning Data visualization of isotope biplots ########
#~~~~~~~~~~~~~~~~~~~~~~~~~

## Carbon and Nitrogen, all data
collcolors <- c("#6B1BD8","#02B6E6","#469425")

tiff("CarbvNit.tiff", width = 700, height = 700)
alldataplot <- ggplot(alldata, aes(x = Car, y = Nit, shape = Collection, color = Site, alpha = Collection)) + 
  geom_point(size = 4) +
  scale_color_manual(values = collcolors) +
  scale_alpha_manual(values = c(1, 0.6, 0.3)) +
  ylim(c(4,10)) + 
  xlim(c(-40,-23))+ 
  theme(panel.background = element_blank(),                              #Removes the background
        plot.margin = unit(x=c(5,13,5,5), units = "mm"),              #Alters margin around plot
        panel.grid = element_blank(),                                    #Eliminates grid lines
        axis.line.x = element_line(colour = "black", linewidth = 0.5),   #Color and thickness of axis lines
        axis.line.y = element_line(colour = "black", linewidth = 0.5), 
        axis.text.x = element_text(size = 14),                           #Adjusts axis tick labeling
        axis.text.y = element_text(size = 14),
        axis.title.x=element_text(colour="black",size=16, vjust=1.3),            #font size for x-axis labels, vjust sets the distance between label and axis
        axis.title.y=element_text(colour="black",margin=margin(0,10,0,0),size=16, angle=90, vjust=1.3),  #font size for y-axis label
        axis.ticks.x=element_line(colour="black", linewidth = 0.5),                      #size of x-axis tick marks
        axis.ticks.y=element_line(colour="black", linewidth = 0.5),                      #size of y-axis tick marks
        strip.background=element_blank()) +
  ylab(expression("δ" ^ 15*"Ν")) + 
  xlab(expression("δ" ^ 13*"C"))
alldataplot

dev.off()

## Carbon and Sulfur all data

tiff("CarbvSulf.tiff", width = 700, height = 700)
allsulfplot <- ggplot(alldata, aes(x = Car, y = Sul, shape = Collection, color = Site, alpha = Collection)) + 
  geom_point(size = 4) +
  scale_color_manual(values = collcolors) +
  scale_alpha_manual(values = c(1, 0.6, 0.3)) +
  scale_shape_manual(values=c(19, 17, 15)) +
  ylim(c(-17,2)) + 
  xlim(c(-40, -23)) +
  theme(panel.background = element_blank(),                             
        plot.margin = unit(x=c(5,13,5,5), units = "mm"),              
        panel.grid = element_blank(),                                    
        axis.line.x = element_line(colour = "black", linewidth = 0.5),   
        axis.line.y = element_line(colour = "black", linewidth = 0.5), 
        axis.text.x = element_text(size = 14),                           
        axis.text.y = element_text(size = 14),
        axis.title.x=element_text(colour="black",size=16, vjust=1.3),            
        axis.title.y=element_text(colour="black",margin=margin(0,10,0,0),size=16, angle=90, vjust=1.3),  
        axis.ticks.x=element_line(colour="black", linewidth = 0.5),                      
        axis.ticks.y=element_line(colour="black", linewidth = 0.5),                      
        strip.background=element_blank()) +
  ylab(expression("δ" ^ 34*"S")) + 
  xlab(expression("δ" ^ 13*"C"))
allsulfplot

dev.off()

## Load confidence intervals Average Carbon and Nitrogen with Error Bars
avecarbCI <- ave$C_95CI
avenitCI <- ave$N_95CI
avesulfCI <- avesulf$S_95CI
carbCI <- avesulf$C_95CI

## Average Carbon and Nitrogen biplot with error bars
tiff("AveCarbvNit.tiff", width = 700, height = 700)
avecarbVnit <-ggplot(ave, aes(x = Carb, y = Nit, shape = Collection, color = Site)) + 
  geom_point(size = 4) + 
  scale_alpha_manual(values = c(1, 0.6, 0.3)) +
  scale_color_manual(values = collcolors) + 
  geom_pointrange(aes(ymin = Nit-avenitCI, ymax =  Nit+avenitCI)) + 
  geom_pointrange(aes(xmin = Carb-avecarbCI, xmax = Carb+avecarbCI)) + 
  ylim(c(4,10)) + 
  xlim(c(-40,-23)) +
  theme(panel.background = element_blank(),                              #Removes the background
        plot.margin = unit(x=c(5,13,5,5), units = "mm"),              #Alters margin around plot
        panel.grid = element_blank(),                                    #Eliminates grid lines
        axis.line.x = element_line(colour = "black", linewidth = 0.5),   #Color and thickness of axis lines
        axis.line.y = element_line(colour = "black", linewidth = 0.5), 
        axis.text.x = element_text(size = 14),                           #Adjusts axis tick labeling
        axis.text.y = element_text(size = 14),
        axis.title.x=element_text(colour="black",size=16, vjust=1.3),            #font size for x-axis labels, vjust sets the distance between label and axis
        axis.title.y=element_text(colour="black",margin=margin(0,10,0,0),size=16, angle=90, vjust=1.3),  #font size for y-axis label
        axis.ticks.x=element_line(colour="black", linewidth = 0.5),                      #size of x-axis tick marks
        axis.ticks.y=element_line(colour="black", linewidth = 0.5),                      #size of y-axis tick marks
        strip.background=element_blank()) +
  ylab(expression("δ" ^ 15*"Ν")) + 
  xlab(expression("δ" ^ 13*"C"))
avecarbVnit

dev.off()


## Average Carbon and Sulfur biplot with error bars
tiff("AveCarbvSulf.tiff", width = 700, height = 700)
avecarbVsulf <-ggplot(avesulf, aes(x = Carb, y = Sulf, shape = Collection, color = Site)) + 
  geom_point(size = 4) + 
  scale_alpha_manual(values = c(0.6, 0.3)) +
  scale_color_manual(values = collcolors) + 
  scale_shape_manual(values=c(17, 15)) +
  geom_pointrange(aes(ymin = Sulf-avesulfCI, ymax =  Sulf+avesulfCI)) + 
  geom_pointrange(aes(xmin = Carb-carbCI, xmax = Carb+carbCI)) + 
  ylim(c(-17,2)) + 
  xlim(c(-40, -23)) +
  theme(panel.background = element_blank(),                              
        plot.margin = unit(x=c(5,13,5,5), units = "mm"),             
        panel.grid = element_blank(),                                   
        axis.line.x = element_line(colour = "black", linewidth = 0.5),   
        axis.line.y = element_line(colour = "black", linewidth = 0.5), 
        axis.text.x = element_text(size = 14),                        
        axis.text.y = element_text(size = 14),
        axis.title.x=element_text(colour="black",size=16, vjust=1.3),            
        axis.title.y=element_text(colour="black",margin=margin(0,10,0,0),size=16, angle=90, vjust=1.3), 
        axis.ticks.x=element_line(colour="black", linewidth = 0.5),                      
        axis.ticks.y=element_line(colour="black", linewidth = 0.5),                      
        strip.background=element_blank()) +
  ylab(expression("δ" ^ 34*"S")) + 
  xlab(expression("δ" ^ 13*"C")) 
avecarbVsulf

dev.off()


