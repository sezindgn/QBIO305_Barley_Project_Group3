
# Variance Decomposition ANOVA 

library(ggplot2)    # load the ggplot2 library for plotting 
library(car)        # load the car  library to make a QQplot
library(cowplot)    # allows us to combine multiple plots in one figure
options(scipen=999) # remove scientific notation of very low integers

# load the data
dta.raw <- read.csv("dta---ANOVA_QTLanalysis.csv")

# number of replicates
length(unique(dta.raw$Group.Replicate))

# conduct the raw ANOVA by inserting the correct model
model.raw <- aov(Height ~ Genotype, data = dta.raw)

# extract the residuals
resid.raw <- model.raw$residuals
resid.raw

# plot the histogram
histogram_raw_data = ggplot(data = data.frame(1:length(resid.raw), resid.raw), 
       aes(x = resid.raw)) +
  geom_histogram(aes(y =..density..), color="darkgrey", bins=50) +
  stat_function(fun = dnorm, 
                args = list(mean = mean(resid.raw), 
                            sd = sd(resid.raw)),
                color="darkred")

# plot qq-plot
qqPlot(resid.raw, envelope = F)

# check for homoscedasticity and outliers 
plot(model.raw, which = 1)

# remove outlier(s)
dta.clean <- dta.raw[-c(14,40,182),]

# conduct the ANOVA with the cleaned data
model.clean <- aov(Height ~ Genotype, data = dta.clean)

# extract the residuals of the new ANOVA
resid.clean <- model.clean$residuals

# check normality and homoscedasticity again
histogram_clean_data = ggplot(data = data.frame(1:length(resid.clean), resid.clean), 
       aes(x = resid.clean)) +
  geom_histogram(aes(y =..density..), color="darkgrey", bins=50) +
  stat_function(fun = dnorm, 
                args = list(mean = mean(resid.clean), 
                            sd = sd(resid.clean)),
                color="darkred")
# plot qq-plot
qqplot_clean = qqPlot(resid.clean, envelope = F)

# check for homoscedasticity and outliers 
plot(model.clean, which = 1)

# check the ANOVA table
summary(model.raw)
summary(model.clean)

# we'll use the 'plot_grid' function from the cowplot package to put all histograms and qq- plots,respectively in a figure
plot_grid(histogram_raw_data, histogram_clean_data, labels = c('Histogram:Raw Data', 'Histogram:Clean Data'), label_size = 12,nrow = 2) # to plot two histograms in a figure 

##plot(qqplot_raw, qqplot_clean, labels = c('QQ-Plot:Raw Data', 'QQ-Plot:Clean Data'), label_size = 12,nrow = 2) # to plot 2 qq-plots in a figure

# Create QQ plot for clean residuals
qqplot_clean <- qqPlot(resid.clean, envelope = FALSE, main = "QQ Plot: Clean Residuals")

# Create QQ plot for raw residuals
qqplot_raw <- qqPlot(resid.raw, envelope = FALSE, main = "QQ Plot: Raw Residuals")

# Set up a 1x2 layout for two plots in one row
par(mfrow = c(1, 2))
qqPlot(resid.raw, envelope = F)
qqPlot(resid.clean, envelope = F)
# Reset the graphical parameters to the default
par(mfrow = c(1, 1))

##plot_grid(plot_raw, plot_clean, labels = c('Raw Data', 'Clean Data'), label_size = 12,nrow = 2) # to plot raw and clean data 
par(mfrow = c(1, 2))
plot(model.raw, which = 1,main = "Raw Data")
plot(model.clean, which = 1,main = "Clean Data")
par(mfrow = c(1, 1))

# QTL 

library(qqman)

map <- read.csv("input_map.csv")

# loop through all markers and conduct ANOVAs for each one of them
lis.out <- list()

for(m in 5:ncol(dta.clean)){
  # extract the first two columns (group + spike length) and the m-th column (=marker)
  loop.dta <- dta.clean[,c(2,3,m)]
  # save the name of the current marker
  loop.marker.name <- colnames(loop.dta)[3]
  # in the model in the aov command, we will have to define the independent 
  # variable name. since the marker name changes in each loop, we need to 
  # to change the column names here to have the same marker name in each loop
  colnames(loop.dta) <- c("group","trait", "allele")
  # conduct one-way ANOVA
  loop.aov <- aov(trait ~ allele + group, loop.dta)
  # extract the allele's p-value
  loop.pval.allele <- summary(loop.aov)[[1]][1,5]
  # extract the marker's genetic position from the genetic map
  loop.map <- map[which(map$marker == loop.marker.name),]
  # create the output data frame
  loop.df.out <- data.frame(loop.map,
                            pval.allele=loop.pval.allele)
  # save the output data frame in the list
  lis.out[[m]] <- loop.df.out
}

# combine the loop's output data frames into a single data frame
res.aov <- do.call("rbind", lis.out)

# have a look at the output data frame
head(res.aov)
dim(res.aov)

# calculate the negative logarithm of the Bonferroni corrected significance threshold
sig.threshold.BonfCorrected <- -log10(0.05/nrow(res.aov))

# create a data frame that follows the requirements of the manhattan command of the qqman library
dta.plot <- data.frame(CHR = as.numeric(gsub("H","",res.aov$chr)),
                       BP = res.aov$pos,
                       SNP = 1:nrow(res.aov),
                       pval.allele = res.aov$pval.allele)

# plot genome-wide p-values for markers
manhattan(dta.plot, genomewideline = sig.threshold.BonfCorrected,
          suggestiveline = F, logp=T, p="pval.allele", type="l", 
          lwd=3, ylab="-log10(p-values)", main="Marker effect")


# Genetic Linkage Map 

library(qtl)    # load the qtl library for linkage analysis

#read in linkage map and qualitative phenotypes as denoted in marker annotation
owb <- read.cross("csv", ".", "owb_linkage_map_qualt_phenotypes_WS2324.csv", genotypes=c("a","b"), 
                  alleles=c("a", "b"), crosstype="dh")
#have a look at the complete linkage map
plotMap(owb, main="", show.marker.names=T)

#To map the locus, calculate all pairwise recombination frequencies and LOD scores
awn_length <- tryallpositions(owb, "awn_length", error.prob=0)
#show the best linkage to a marker from each chromosome
summary(awn_length)
#move the locus to the best marker position
owb <- movemarker(owb, "awn_length", "7H", 3.35027)
#update map
plotMap(owb, main="", show.marker.names=T)

