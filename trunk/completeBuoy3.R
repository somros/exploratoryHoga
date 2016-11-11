# 8.11.2016

# General script to for community analysis of the sponges at Buoy 3, Wakatobi Marine Park, Indonesia.
# Starting dataset is sponge counts per quadrats (square meter) from 2005 to 2016 (ongoing). Two years are missing
# (2010, 2012). Sampling desing is as follows: vertical wall on a pristine (or at least not impacted) reef. Three sites,
# A, B, C, replicate of each other. Sites do not recruit from each other. In each site, 5 quadrats.



require(abind)
require(plyr)
require(ggplot2)
require(reshape2)
require(BiodiversityR)
require(ggvegan)

dataAllYears <- read.csv("//Staff/Home/SCIFAC/rovellal/DocumentsRedir/Data/Hoga/buoy3/spongeAbundanceQuadratsCleaned.csv")

# set flags

flagRoutine <- "average" # can be either "sum" to add the 5 quadrats up, or any other string to take average and sd

# load species list

species <- read.csv("//Staff/Home/SCIFAC/rovellal/DocumentsRedir/Data/Hoga/buoy3/speciesKeyCleaned.csv")
species[species=="" | species==0] <- NA # drops empty spaces and zeroes

# substitute description to blank fields in the species column

speciesOrDescription <- list()

for (i in 1:nrow(species)) {
  if (is.na(species[,2][i])==T) {
    speciesOrDescription[[i]] <- species[,3][i]
  } else {
    speciesOrDescription[[i]] <- species[,2][i]
  }
}
speciesOrDescription <- unlist(speciesOrDescription) #

# extract general information about the dataset, to be called later 

years <- levels(factor(substr(names(dataAllYears[-1]),2,3)))
quadrats <- as.numeric(levels(factor(substr(names(dataAllYears), nchar(names(dataAllYears)), nchar(names(dataAllYears))))))
quadrats <- quadrats[is.na(quadrats)==F]# number of quadrats per site
dataColumns <- ncol(dataAllYears)-1 # number of quadrats in all sites in all years, effectively number of columns in the frame
sites <- c("A", "B", "C")


# gets rid of zeros as characters in the datasets

for (i in 1:ncol(dataAllYears)) {
  if (is.numeric(dataAllYears[,i])==T) {
    dataAllYears[,i] <- dataAllYears[,i]
  } else {
    dataAllYears[,i] <- as.numeric(levels(dataAllYears[,i])[dataAllYears[,i]]) # this removes everything that won't fit as numeric
  }
}

dataAllYears[is.na(dataAllYears)] <- 0 # turns NAs to zeroes

# gets rid of the X in front of the column names

colnames(dataAllYears) <- c(names(dataAllYears)[1],
                            substr(names(dataAllYears[,-1]), 2, nchar(names(dataAllYears[,-1]))))


# this data frame is essentially the raw data and starting point for a series of plots and analysis coming up later in the script




# *******************************************************************************************************************************

# goal from now on is taking an average of the quadrats for site means

# set indexing to operate across pools of columns, 5 in our case as we want to take averages for one site across all the quadrats

dataAllYears <- dataAllYears[,-1] # run this!

n <- 1:ncol(dataAllYears)
ind <- data.frame(matrix(c(n, rep(NA, 5 - ncol(dataAllYears)%%5)), byrow=F, nrow = 5))
nonna <- sapply(ind, function(x) all(!is.na(x)))
ind <- ind[, nonna]

# if statement to decide whether we want to add the quadrats or take their average, pointing to initial flag

if (flagRoutine == "sum") {
  electedFunction <- function(i)rowSums(dataAllYears[,i])
} else {
  electedFunction <- function(i) {
    meanColumn <- rowMeans(dataAllYears[,i])
    sdColumn <- apply(dataAllYears[,i], 1, sd)
    cbind(meanColumn, sdColumn) 
  }
}

sitesFrame <- as.data.frame(do.call(cbind, lapply(ind, electedFunction))) 

# reorganize the dataframe with first the means and then the sd

sitesFrame <- cbind(speciesOrDescription,
                    sitesFrame[,grepl("mean", names(sitesFrame))], 
                    sitesFrame[,grepl("sd", names(sitesFrame))])

# prepare vector with the column names ofor the site averages

myNames <- list()

for (i in 1:length(years)) {
  myNames[[i]] <- paste(years[i], sites, sep = "")
}


colnames(sitesFrame) <- c("Species", paste("Mean", unlist(myNames), sep = ""), paste("SD", unlist(myNames), sep = ""))

# what we have here is the data frame of the mean values per site, i.e. the average of all the quadrats has been taken


# ********************************************************************************************************************

# subset the mean abundances per site to the species accounting for 90% of the sponge numbers. 
# List of species may be slightly reduced than if performed on the single quadrats. Also MDS cares about differences between
# sites in time, not quadrats, so this is actually legitimate if we do not pool species in undetermined (in which case it
# does make a relevant difference if we consider quadrats or site means). The only point I'm not comfortable with is that
# MDS+PCA operated on the means per instance kills the uncertainty across quadrats, which kind of belittles the point of having
# replicates in the first place and also is not statistically sound, but then again the ordination on the quadrats makes no sense.
# All in all it's about experimental design here.

sortedMostAbundant <- list()
for (i in 1:((length(sitesFrame)+1)/2)) {
  sortedMostAbundant[[i]] <- sitesFrame[order(sitesFrame[,i], decreasing = T),]
  sortedMostAbundant[[i]] <- sortedMostAbundant[[i]][,c(1,i)]
}

sortedMostAbundant <- sortedMostAbundant[-1]

# can decide that dominant species together account for 90% of the numbers


chainsaw <- function(x, percent) {
  # first it needs to calculate the sum of the numbers
  sumOfNumbers <- sum(x[,2])
  # then calculate n as in n = sumOfNumbers*domPerc / 100
  n = sumOfNumbers*percent/100
  # # then count how many elements of x[,2] are needed to go as big as n, where the count is z
  v <- 0
  z <- 0
  for (i in 1:nrow(x)) {
    if (v < n) {
      v <- v + x[i,2]
      z <- z + 1
    } else {
      v <- v
      z <- z
    }
  }
  # # finally subset the original frame x taking the first z counts, which should add up to the closest possible to
  # # the chosen percentage
  y <- x[1:z,]
  return(y)
  
}

dominantPerc <- lapply(sortedMostAbundant, chainsaw, 90)

speciesSet <- factor(levels(factor(unlist(lapply(dominantPerc, function(x) as.character(x[,1]))))))

# now use it to subset the meanDataFrame, which must be entry for barchart

subsetData <- sitesFrame[sitesFrame[,1] %in% speciesSet,] 

# need to separate means from SD and melt them

means <- subsetData[,1:((ncol(subsetData)+1)/2)]
sd <- subsetData[c(1, (((ncol(subsetData)+1)/2)+1):ncol(subsetData))]
  
meltMean <- melt(means, id.vars = "Species", variable.name = "Instance", value.name = "Mean")
meltSD <- melt(sd, id.vars = "Species", variable.name = "Instance", value.name = "Mean")
meltAll <- cbind(meltMean, meltSD[,ncol(meltSD)])
meltAll$Instance <- as.character(meltAll$Instance)
colnames(meltAll) <- c(names(meltAll)[-length(meltAll)], "SD")
meltAll$Year <- as.numeric(substr(meltAll$Instance, nchar(meltAll$Instance)-2, nchar(meltAll$Instance)-1))
meltAll$Site <- factor(substr(meltAll$Instance, nchar(meltAll$Instance), nchar(meltAll$Instance)))



# 1. barchart for community composition. Rather uninformative with the classification at speceis level

barPlot <- ggplot(data = meltAll, aes(x = Year, y = Mean, fill = Species, group = Site))+
  geom_bar(stat = "identity", position = "stack")+
  scale_x_continuous(breaks = seq(5,16,1),
                     labels = seq(5,16,1),
                     limits = c(4,17))+
  scale_y_continuous(limits = c(0,180),
                     breaks = seq(0,180,20))+
  #scale_fill_manual(values = getPalette(nOfColors))+
  labs(y=expression(paste(Sponge~density~(n~m^2))))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank())+
  theme(plot.title = element_text(size=14, vjust=2))+
  theme(axis.title.x = element_text(size=10,vjust=-0.5),
        axis.title.y = element_text(size=10,vjust=0.5))+
  theme(axis.text.x=element_text(size=10))+
  theme(axis.text.y=element_text(size=10))+
  #theme(legend.position="none")+
  facet_grid(.~Site)
barPlot


# *****************************************************************************************************************


# 2. Single species line plots (not visualising uncertainty)

linesPlot <- ggplot(data = meltAll, aes(x = Year, y = log(Mean+1), color = Species))+
  geom_line(size = 1)+
  geom_point()+
  scale_x_continuous(breaks = seq(5,16,1),
                     labels = seq(5,16,1),
                     limits = c(4,17))+
  scale_y_continuous(limits = c(0,5),
                     breaks = seq(0,5,.5))+
  #scale_fill_manual(values = getPalette(nOfColors))+
  labs(y=expression(paste(log(Sponge~density~(n~m^2)))))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank())+
  theme(plot.title = element_text(size=14, vjust=2))+
  theme(axis.title.x = element_text(size=10,vjust=-0.5),
        axis.title.y = element_text(size=10,vjust=0.5))+
  theme(axis.text.x=element_text(size=10))+
  theme(axis.text.y=element_text(size=10))+
  #theme(legend.position="none")+
  facet_grid(.~Site)
linesPlot




# *****************************************************************************************************************

# 3. Total species line plots with gaps for missing years

# go back to the original data frame (dataAllYears) and transform it. 
# No species reduction, we do not care, we work on totals here

transDatatmp <- as.data.frame(t(dataAllYears))
transDatatmp <- transDatatmp[-1,] # gets rid of the line with the species IDs

# need to turn it back to numeric because it all beccame factor in the meantime

transData <- as.data.frame(apply(transDatatmp, 2, function(x) as.numeric(as.character(x))))

rownames(transData) <- rownames(transDatatmp)


transData$Total <- rowSums(transData)

# need to separate years (my x variable) and quadrats

quadratID <- rownames(transData)
transData$Year <- as.numeric(substr(quadratID, 1, 2))
transData$Site <- substr(quadratID, nchar(quadratID)-1, nchar(quadratID)-1)
transData$Quadrat <- substr(quadratID, nchar(quadratID), nchar(quadratID))  
transData$QuadratSite <- substr(quadratID, nchar(quadratID)-1, nchar(quadratID))

# restrict to small data frame with only the total numbers and the quadrats IDs

smallFrame <- transData[,(ncol(transData)-4):ncol(transData)]

library(RColorBrewer)
par(mar = c(0, 4, 0, 0))
display.brewer.all()
nOfColorsQuad <- 5 # numer of quadrats
getPaletteQuad <- colorRampPalette(brewer.pal(5, "Set1"))

quadratsFate <- ggplot(data = smallFrame, aes(x = Year, y = Total, group = QuadratSite, color = Quadrat))+
  geom_line(aes(color = Quadrat), size = 1)+
  geom_point()+
  scale_color_manual(values = getPaletteQuad(nOfColorsQuad))+
  scale_x_continuous(breaks = seq(5,16,1),
                     labels = seq(5,16,1),
                     limits = c(5,16))+
  scale_y_continuous(limits = c(0,350),
                     breaks = seq(0,350,50))+
  labs(y=expression(paste(Sponge~density~(n~m^2))))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank())+
  theme(plot.title = element_text(size=14, vjust=2))+
  theme(axis.title.x = element_text(size=10,vjust=-0.5),
        axis.title.y = element_text(size=10,vjust=0.5))+
  theme(axis.text.x=element_text(size=10, angle = 60, 
                                 hjust = 1, vjust = .9))+
  theme(axis.text.y=element_text(size=10))+
  facet_grid(. ~ Site, scales = "free")
quadratsFate




# ***************************************************************************************************************

# 4. Growth of total numbers factorized by quadrats
# Uses the same data frame of point 3. 

smallFrame$TotalNextYear <- c(smallFrame$Total[(smallFrame$Year != 5)], 
                              rep(NA, nrow(smallFrame[smallFrame$Year == 5,])))

spongesInTime <- ggplot(data = smallFrame, 
                        aes(x = Total, y = TotalNextYear, group = Quadrat))+
  geom_point(aes(color = Quadrat))+
  scale_color_manual(values = getPaletteQuad(nOfColorsQuad))+
  scale_x_continuous(breaks = seq(0,350,50),
                     labels = seq(0,350,50),
                     limits = c(0,350))+
  scale_y_continuous(limits = c(0,350),
                     breaks = seq(0,350,50))+
  labs(x="D(t)", y="D(t+1)")+
  #geom_smooth(method = "lm")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank())+
  theme(plot.title = element_text(size=14, vjust=2))+
  theme(axis.title.x = element_text(size=10,vjust=-0.5),
        axis.title.y = element_text(size=10,vjust=0.5))+
  theme(axis.text.x=element_text(size=10, vjust = .9))+
  theme(axis.text.y=element_text(size=10))+
  facet_grid(. ~ Site)
spongesInTime




# *******************************************************************************************************




# Combination of nMDS and PCA with PERMANOVA or NPMANOVA, operating from Bray-Curtis matrix

allMeans <- means[,-1]

allMeansTrans <- log10(allMeans+1)


# first we transpose (which is not necessary now but should become at some point I reckon)

transMatrix <- t(allMeansTrans)

colnames(transMatrix) <- speciesOrDescription[as.numeric(rownames(allMeans))]


########################################################
########################################################
########################################################
# EXPORT FOR PRIMER

forPrimer <- t(allMeans)
colnames(forPrimer) <- speciesOrDescription[as.numeric(rownames(allMeans))]
#write.csv(forPrimer, "//Staff/Home/SCIFAC/rovellal/DocumentsRedir/Data/Hoga/buoy3/primerData.csv")

########################################################
########################################################
########################################################
                                         



# now structure is N as rows and species as columns. However, we osu old structure for computation (same ending point and I don't
# have to recode it all up)

grandBC <- list()

for (i in 1:length(allMeansTrans)) { # home-made Bray-Curtis, makes little sense as there's vegdist in vegan doing it but it's solid (validated)
  BC <- vector(mode = "logical", length = length(allMeansTrans))
  for (j in 1:length(allMeansTrans)) {
    reduced <- allMeansTrans[,c(i,j)]
    BCpartial <- apply(reduced, 1, function(x) {
      numBC <- 2*min(x[1], x[2])
      denBC <- x[1]+x[2]
      metrics <- c(numBC, denBC)
    })
    BC[j] <- rowSums(BCpartial)[1]/rowSums(BCpartial)[2]
  }
  grandBC[[i]] <- BC
}

simMatrix <- abind(grandBC, along = 0)
colnames(simMatrix) <- names(allMeansTrans)
rownames(simMatrix) <- names(simMatrix)

dissMatrix <- 1-simMatrix # turn into dissimilarity

fit <- isoMDS(dissMatrix, k = 2) # this is a ***NON-METRIC*** MDS plot based on Brey-Curtis dissimilarity.

distances <- as.data.frame(fit$points)
rownames(distances) <- names(allMeansTrans)
distances$Names <- as.character(rownames(distances))
for (i in 1:nrow(distances)) {
  distances$Year[i] <- paste("20", substr(distances$Names[i], 
                                          (nchar(distances$Names[i])-2), 
                                          (nchar(distances$Names[i])-1)),
                             sep = "")
}
for (i in 1:nrow(distances)) {
  distances$Site[i] <- substr(distances$Names[i], nchar(distances$Names[i]), nchar(distances$Names[i]))
}

# plot

# validated against PRIMER

MDSplot <- ggplot(data = distances, aes(x = V1, y = V2, group = Site, color = Year))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_vline(xintercept = 0, linetype = "dashed")+
  geom_text(aes(label = Year, color = Year), size = 4.5, 
            fontface = "bold", vjust = 0, nudge_y = 0.02)+
  guides(color=FALSE)+
  geom_point(aes(shape = Site, color = Year), size = 2)+
  scale_shape_manual(values = c(0,1,2))+
  labs(x = "MDS1", y = "MDS2")+
  #geom_path(size = 1)+
  scale_x_continuous(limits = c(-.7,.7))+
  scale_y_continuous(limits = c(-.7,.7))+
  theme_bw()+
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank())#+
#facet_grid(.~Site)
MDSplot

# ggsave("//Staff/Home/SCIFAC/rovellal/DocumentsRedir/Wellington/words/latex/hogaExploratory/Pics/MDSplot.pdf",
#        MDSplot, useDingbats = T, width = 6, height = 5)


# *********************************************************************************************************************************

# CAP uses a dedicated package and starts from the species list

# need to build a dummy frame apparently with data

dummyMatrix <- as.data.frame(cbind(distances$Year, distances$Site)) # shouldn't this include the species somehow too?
colnames(dummyMatrix) <- c("Year", "Site")
rownames(dummyMatrix) <- rownames(transMatrix)

# different formulas mean different ordinations, depending on the included variables

CAPall <- capscale(transMatrix ~ Year + Site, dummyMatrix, distance = "bray", add = T)
CAPTime <- capscale(transMatrix ~ Year, dummyMatrix, distance = "bray", add = T)
CAPSpace <- capscale(transMatrix ~ Site, dummyMatrix, distance = "bray", add = T)

# plots for different formulas

# look like those obtained in PRIMER just a bit off

plot(CAPall)
plot(CAPTime)
plot(CAPSpace)

# all the below is just to get fancy with ggplot. Packages vegan and ggvegan will do all the job

components <- fortify(CAPall) # extract the coordinates from the cca objects calculated with capscale 
components <- data.frame(lapply(components, function(x) {
  x <- gsub("Year", "", x)
  x <- gsub("Site", "", x)
  x <- gsub("Mean", "", x)
  return(x)
}))
components$Dim1 <- as.numeric(as.character(components$Dim1))
components$Dim2 <- as.numeric(as.character(components$Dim2))
components$Label <- gsub("spe", "", components$Label)


# plotting region

plotData <- subset(components, components$Score != "biplot" & components$Score != "constraints" &
                     components$Score != "sites")
plotData$ColorKey <- c(rep("Species", nrow(plotData[plotData$Score=="species",])), 
                       rep("Year", 10), rep("Site", 3)) # make this flexible on the iterations
# fix levels
plotData$ColorKey <- factor(plotData$ColorKey, levels = unique(plotData$ColorKey))


CAPplot <- ggplot(data = plotData, aes(x = Dim1, y = Dim2, group = ColorKey))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_vline(xintercept = 0, linetype = "dashed")+
  geom_segment(data = subset(plotData, plotData$Score == "centroids"),
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2), arrow = arrow(length = unit(1/2, 'picas')))+
  geom_text(aes(label = Label, color = ColorKey, size = ColorKey))+
  scale_color_manual(values = c("darkgrey", "red3", "cyan4"))+
  scale_size_discrete(range = c(4,8))+
  theme_bw()+
  labs(x = "CAP1", y = "CAP2")+
  scale_x_continuous(limits = c(-.8,.8))+
  scale_y_continuous(limits = c(-.8,.8))+
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank())
CAPplot

# ggsave("//Staff/Home/SCIFAC/rovellal/DocumentsRedir/Wellington/words/latex/hogaExploratory/Pics/CAPplot.pdf",
#        CAPplot, useDingbats = T, width = 6.5, height = 5)


# PERMANOVA

# validated against PRIMER, same results, but I'm not sure it makes any sense. 
# Lack of replication (because we do not run it over the quadrats)

# The community is different but we don't know why. Less conceptually, how do we ordinate a dataset with no
# environmental covariates?

adonis(transMatrix~Year + Site, data = dummyMatrix, permutations = 999, method = "bray")

# yeah validated with PRIMER and all. I will actually have to give myself that I did a good job with more 
# luck than skill with this routine here, it's all mint.


