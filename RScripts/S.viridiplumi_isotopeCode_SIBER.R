###### Carbon, Nitrogen, and Sulfur SIBER Analysis  ########
##### Statistical analysis and data visualization for Seepicola viridiplumi microbiome manuscript

## Rice LN, Louca, S, Cordes EE, Young CM. (2026). Trophic niches and the associated
## microbial community for a newly described methane seep sabellid. Deep-Sea Research
## Part 1

## This code is used to analyze the stable isotopic signatures and trophic niche overlap from Seepicola 
## viridiplumi radioles and anterior tissues for both free-living and epibiotic individuals. 
## Epibiotic samples were inconsistent by site and collection time.


######## Set Up ########
## Set Working directory:
#setwd() #CHANGE ME to desired directory

## Load packages
library(readxl)
library(ggplot2)
library(ggforce)
library(rcompanion)
library(SIBER)

## Loading the Needed Data Sheets:
bpcomp <- read.csv("BPisocomp.csv", header = T)       #iso1 = Nitrogen, iso2 = Carbon
bpsulf <- read.csv("bpsulfcomp.csv", header = T)      #iso1 = Sulfur, iso2 = Carbon
bh <- read.csv("BHisocomp.csv", header = T)           #iso1 = Nitrogen, iso2 = Carbon
bhsulf <- read.csv("BHsulfcomp.csv", header = T)      #iso1 = Sulfur, iso2 = Carbon
gc <- read.csv("GCisocomp.csv", header = T)           #iso1 = Nitrogen, iso2 = Carbon, removed 1 sample due to missing Nitrogen measurement


######## Generate plotting arguments to be passed onwards through the visualization functions ########
community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args <- list(n = 100, p.interval = 0.95, lty = 1, lwd = 2)
group.hulls.args <- list(lty = 2, col = "grey20")

#===========================================
######## Examining Nitrogen and Carbon at Brine Pool NR1 ########
#===========================================
##  Make SIBER Object
bpcomp$group <- as.factor(bpcomp$group)
bpcomp$community <- as.factor(bpcomp$community)

bpsiber <- createSiberObject(bpcomp)

## Plot the niche overlap estimates
brine <- plotSiberObject(bpsiber, ax.pad = 4,
                         hulls = F, community.hulls.args,
                         ellipses = T, group.ellipses.args,
                         group.hulls = F, group.hulls.args,
                         bty = "L",
                         iso.order = c(2,1),
                         points.order = c(16,16),
                         xlab = expression("δ" ^ 13*"C"),
                         ylab = expression("δ" ^ 15*"N"),
                         x.limits = c(-44,-22),
                         y.limits = c(2,14))

bpsum <- groupMetricsML(bpsiber)
print(bpsum)        # TA = Total Area, SEA = Standard Ellipse Area, SEAc = Standard Ellipse Area corrected
        
## Calculate Ellipse Overlap by area and proportion
bp <- maxLikOverlap("Brine Pool NR1.Epibiotic", "Brine Pool NR1.Free-living", bpsiber, p = 0.95)
propbp <- as.numeric(bp["overlap"] / (bp["area.1"] + bp["area.2"] - bp["overlap"]) *100)
print(propbp) 


######## Carbon and Sulfur at Brine Pool NR1 ########
##  Make SIBER Object
bpsulf$group <- as.factor(bpsulf$group)
bpsulf$community <- as.factor(bpsulf$community)

bpsulfsiber <- createSiberObject(bpsulf)

## Plot the niche overlap estimates
brinesulf <- plotSiberObject(bpsulfsiber, ax.pad = 4,
                         hulls = F, community.hulls.args,
                         ellipses = T, group.ellipses.args,
                         group.hulls = F, group.hulls.args,
                         bty = "L",
                         iso.order = c(2,1),
                         points.order = c(16,16),
                         xlab = expression("δ" ^ 13*"C"),
                         ylab = expression("δ" ^ 34*"S"),
                         x.limits = c(-44,-22),
                         y.limits = c(-22,4))

bpsulfsum <- groupMetricsML(bpsulfsiber)
print(bpsulfsum)        # TA = Total Area, SEA = Standard Ellipse Area, SEAc = Standard Ellipse Area corrected

## Calculate Ellipse Overlap by area and proportion
sulfbp <- maxLikOverlap("Brine Pool NR1.Epibiotic", "Brine Pool NR1.Free-living", bpsulfsiber, p = 0.95)
propbpsulf <- as.numeric(sulfbp["overlap"] / (sulfbp["area.1"] + sulfbp["area.2"] - sulfbp["overlap"]) *100)
print(propbpsulf) 

######## Carbon and Nitrogen at Bush Hill ########
##  Make SIBER Object
bh$group <- as.factor(bh$group)
bh$community <- as.factor(bh$community)

bhsiber <- createSiberObject(bh)

## Plot the niche overlap estimates
bush <- plotSiberObject(bhsiber, ax.pad = 4,
                             hulls = F, community.hulls.args,
                             ellipses = T, group.ellipses.args,
                             group.hulls = F, group.hulls.args,
                             bty = "L",
                             iso.order = c(2,1),
                             points.order = c(16,16),
                             xlab = expression("δ" ^ 13*"C"),
                             ylab = expression("δ" ^ 15*"N"),
                             x.limits = c(-44,-22),
                             y.limits = c(2,14))

bhsum <- groupMetricsML(bhsiber)
print(bhsum)        # TA = Total Area, SEA = Standard Ellipse Area, SEAc = Standard Ellipse Area corrected

## Calculate Ellipse Overlap by area and proportion
bhprop <- maxLikOverlap("Bush Hill.Epibiotic", "Bush Hill.Free-living", bhsiber, p = 0.95)
bhoverlap <- as.numeric(bhprop["overlap"] / (bhprop["area.1"] + bhprop["area.2"] - bhprop["overlap"]) *100)
print(bhoverlap) 


######## Carbon and Sulfur at Bush Hill ########
##  Make SIBER Object
bhsulf$group <- as.factor(bhsulf$group)
bhsulf$community <- as.factor(bhsulf$community)

bhsulfsiber <- createSiberObject(bhsulf)

## Plot the niche overlap estimates
sulfbush <- plotSiberObject(bhsulfsiber, ax.pad = 4,
                        hulls = F, community.hulls.args,
                        ellipses = T, group.ellipses.args,
                        group.hulls = F, group.hulls.args,
                        bty = "L",
                        iso.order = c(2,1),
                        points.order = c(16,16),
                        xlab = expression("δ" ^ 13*"C"),
                        ylab = expression("δ" ^ 34*"S"),
                        x.limits = c(-44,-22),
                        y.limits = c(-22,4))

bhsulfsum <- groupMetricsML(bhsulfsiber)
print(bhsulfsum)        # TA = Total Area, SEA = Standard Ellipse Area, SEAc = Standard Ellipse Area corrected

## Calculate Ellipse Overlap by area and proportion
bhsulfprop <- maxLikOverlap("Bush Hill.Epibiotic", "Bush Hill.Free-living", bhsulfsiber, p = 0.95)
bhsulfoverlap <- as.numeric(bhsulfprop["overlap"] / (bhsulfprop["area.1"] + bhsulfprop["area.2"] - bhsulfprop["overlap"]) *100)
print(bhsulfoverlap) 


######## Carbon and Nitrogen at Green Canyon 234 ########
##  Make SIBER Object
gc$group <- as.factor(gc$group)
gc$community <- as.factor(gc$community)

gcsiber <- createSiberObject(gc)

## Plot the niche overlap estimates
gciso <- plotSiberObject(gcsiber, ax.pad = 4,
                            hulls = F, community.hulls.args,
                            ellipses = T, group.ellipses.args,
                            group.hulls = F, group.hulls.args,
                            bty = "L",
                            iso.order = c(2,1),
                            points.order = c(16,16),
                            xlab = expression("δ" ^ 13*"C"),
                            ylab = expression("δ" ^ 15*"N"),
                            x.limits = c(-44,-22),
                            y.limits = c(2,14))

gcsum <- groupMetricsML(gcsiber)
print(gcsum)        # TA = Total Area, SEA = Standard Ellipse Area, SEAc = Standard Ellipse Area corrected

## Calculate Ellipse Overlap by area and proportion
gcoverlap <- maxLikOverlap("Green Canyon 234.Epibiotic", "Green Canyon 234.Free-living", gcsiber, p = 0.95)
gcprop <- as.numeric(gcoverlap["overlap"] / (gcoverlap["area.1"] + gcoverlap["area.2"] - gcoverlap["overlap"]) *100)
print(gcprop) 
