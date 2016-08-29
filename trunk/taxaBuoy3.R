# 23.08.2016 script to substitute the order names to the species in the quadrat data
# so basically fuck R, take the original frame in excel and sit down one hour with WORMS
# then use this script to transpose the frame. shoud become a part of other quadrats script
# do order and family for each entry

# 28.08.2016 Best script so far, to be cleaned, but it does what quadrats fate does and it's 
# more complete. needs to be adapted with the means to do what betterBuoy3 does

require(abind)
require(plyr)
require(ggplot2)
require(reshape2)

setwd("/home/somros/Documents/R/exploratoryHoga/input/")


# sets the flags

flagRoutine <- "average" # can be either "sum" to add the 5 quadrats up, or any other string to take average and sd

#############

taxonData <- read.csv("/home/somros/Documents/R/exploratoryHoga/input/taxonData.csv", sep = ",", header = T)

# using the raw data, find out what entries make up for most of the abundance per quadrat per year, i.e. per column
# first step is to get rid of the description column and fill all the gaps with "Undetermined"

#taxonDataNoDesc <- taxonData[,-1]

# set information about the sampling design to be used later for factorization

years <- levels(factor(substr(names(taxonData[-(1:4)]),2,3)))
quadrats <- as.numeric(levels(factor(substr(names(taxonData), 
                                            nchar(names(taxonData)), 
                                            nchar(names(taxonData))))))
quadrats <- quadrats[is.na(quadrats)==F]# number of quadrats per site
dataColumns <- ncol(taxonData)-4 # number of quadrats in all sites in all years, effectivelty number of columns in the frame
sites <- c("A", "B", "C")

# second step is to change all the empties with Undetermined

taxonData$Species <- as.character(taxonData$Species)
taxonData$Species[taxonData$Species == "" | taxonData$Species == "0"] <- "Undetermined"
taxonData$Species <- as.factor(taxonData$Species) # and turn to factor

# reduce the data to the columns and rows of interest

#taxonData <- taxonDataNoDesc[1:127,]

# get rid of spaces

taxonData[taxonData==""] <- NA # drop all that is not a name, to be refined though

# turn all the columns but the first into numeric

taxonData[,4:length(taxonData)] <- apply(taxonData[,4:length(taxonData)], 2, 
                                         function(x) {
  x <- as.numeric(as.character(x))
  return(x)
})

taxonData[is.na(taxonData)] <- 0 # replace NAs with 0 for operations

# remove X from column names

colnames(taxonData) <- c(names(taxonData)[1:4],
                            substr(names(taxonData[,-(1:4)]), 2, nchar(names(taxonData[,-(1:4)]))))




# first sum all the undetermined together. split lapply combine requires less work after than "by"

splitTaxonData <- split(taxonData, taxonData$Species)
splitTaxonSum <- lapply(splitTaxonData, function(x) {
  y <- matrix(colSums(x[,-(1:4)]), nrow = 1)
  z <- data.frame(x[1,1:4], y)
  return(z)
})
taxonSmall <- as.data.frame(abind(splitTaxonSum, along = 1))
colnames(taxonSmall) <- names(taxonData)
rownames(taxonSmall) <- 1:nrow(taxonSmall)
taxonSmall <- data.frame(taxonSmall[,1:4], apply(taxonSmall[,-(1:4)], 2, function(x) as.numeric(as.character(x))))


# write out a .csv to be kept

#write.csv(taxonSmall, "/home/somros/Documents/R/exploratoryHoga/output/spongeNumbersTrophic.csv")



# taxonDataSmall <- by(taxonData[,-(1:4)], taxonData$Species, FUN = colSums)
# taxonFrameSmall <- as.data.frame(do.call(rbind, taxonDataSmall)) 
# 
# # here row names are the species, order column is lost. ripristinate. polish data frame
# 
# taxonFrameSmall$Species <- rownames(taxonFrameSmall)
# newOrders <- unlist(apply(taxonFrameSmall, 1, function (x) {
#   index <- grep(pattern = x[length(x)], taxonData$Species)[1]
#   myOrder <- taxonData$Order[index]
#   myOrder
# }))
# 
# taxonFrameSmall$Order <- newOrders

# following secion is to calculate most abundant species in each column

# make a list where each element is the last 2 columns and 1 column

# listOfQuadrats <- vector(mode = "list", length = ncol(taxonFrameSmall)-2)
# for (i in 1:(length(taxonFrameSmall)-2)) {
#   listOfQuadrats[[i]] <- data.frame(taxonFrameSmall[,i], taxonFrameSmall[,(length(taxonFrameSmall)-1)],
#   taxonFrameSmall[,length(taxonFrameSmall)])
#   colnames(listOfQuadrats[[i]]) <- c("Number", "Species", "Order")
# }
# 
# # now apply function to sort them all and take the species making up for the max abundance (if lucky I get rid of one species,
# # now that the undetermined are lumped together)
# 
# sortedQuadrats <- lapply(listOfQuadrats, function(x) {
#   x <- x[order(x[,1], decreasing = T),]
#   x
# })
# 
# # now cut. as you are at it you should also calculate what the proportion of the undetermined is
# 
# chainsawRebuild <- function(x) {
#   # first it needs to calculate the sum of the numbers
#   sumOfNumbers <- sum(x[,1])
#   # then calculate n as in n = sumOfNumbers*domPerc / 100
#   n = sumOfNumbers*90/100
#   # # then count how many elements of x[,1] are needed to go as big as n, where the count is z
#   v <- 0
#   z <- 0
#   for (i in 1:nrow(x)) {
#     if (v < n) {
#       v <- v + x[i,1]
#       z <- z + 1
#     } else {
#       v <- v
#       z <- z
#     }
#   }
#   # # finally subset the original frame x taking the first z counts, which should add up to the closest possible to
#   # # the chosen percentage
#   y <- x[1:z,]
#   return(y)
#   
# }
# 
# 
# dominantSpecies <- lapply(sortedQuadrats, chainsawRebuild)
# 
# # there's an issue with some orders, that become NA for whatever reason. First we need to know if all of this is worth it
# 
# # bind together all the dataframes to isolate the species and compare with the untrimmed community
# 
# allFrames <- do.call(rbind, dominantSpecies)
# levels(factor(allFrames[,2]))



# now sum according to the factor, which is the trophic group

trophicSums <- by(taxonSmall[,-(1:4)], taxonSmall$Feeding, FUN = colSums)
trophicMode <- as.data.frame(do.call(rbind, trophicSums))
total <- colSums(trophicMode)
trophicMode <- rbind(trophicMode, total)
rownames(trophicMode) <- c("h", "p", "u", "tot")

# remove X from column names

colnames(trophicMode) <- names(taxonData[-(1:4)])






######################################################################

# brach for average within sites

# code to take average of the columns (recycle from betterBuoy3)
# this routine loses the single quadrats information, therefore it must be a second and
# separate branch of the main trunk

n <- 1:ncol(trophicMode)
ind <- data.frame(matrix(c(n, rep(NA, 5 - ncol(trophicMode)%%5)), byrow=F, nrow = 5))
nonna <- sapply(ind, function(x) all(!is.na(x)))
ind <- ind[, nonna]

# if statement to decide whether we want to add the quadrats or take their average

if (flagRoutine == "sum") {
  electedFunction <- function(i)rowSums(trophicMode[,i])
} else {
  electedFunction <- function(i) {
    meanColumn <- rowMeans(trophicMode[,i])
    varColumn <- apply(trophicMode[,i], 1, sd)
    cbind(meanColumn, (varColumn)^2) # to get the variance, used for the propagation
  }
}

sitesFrame <- as.data.frame(do.call(cbind, lapply(ind, electedFunction))) # hahah nailed it
sitesFrame <- cbind(sitesFrame[,grepl("mean", names(sitesFrame))], 
                    sitesFrame[,grepl("V", names(sitesFrame))])

# skip the grand mean, mean values per site are enough, too much information is lost looking
# at the entire area (although after all they are replicates of the same habitat so it will
# be necessary)

# need to rename the columns for melting (and to understand)

# create array of sites and years

sitesYearsList <- list()
for(i in 1:length(years)) {
  sitesYearsList[[i]] <- paste(sites, years[i], sep = "")
}
sitesYears <- unlist(sitesYearsList)

# apply as column names

colnames(sitesFrame) <- c(paste("mean", sitesYears, sep = ""),
                          paste("sd", sitesYears, sep = ""))

# first transpose then melt, goal is same plot as quadrat fate but site-wise (bars and shit)

transSite <- as.data.frame(t(sitesFrame))

# add relevant column for factorization

siteID <- rownames(transSite)
transSite$Year <- as.numeric(substr(siteID, nchar(siteID)-1, nchar(siteID)))
transSite$Site <- substr(siteID, nchar(siteID)-2, nchar(siteID)-2)
transSite$Stat <- substr(siteID, 1, 2)
#transSite$QuadratSite <- substr(siteID, nchar(siteID)-1, nchar(siteID))

# separate into two frames and melt

transList <- split(transSite, transSite$Stat)
meltTransList <- lapply(transList, function(x) {
  moltenFrame <- melt(x, id.vars = c("Year", "Site", "Stat"), variable.name = "Group",
                      value.name = "Number")
  moltenFrame
}
)
meltTrophicMeans <- as.data.frame(cbind(meltTransList[[1]], 
                                        meltTransList[[2]][,length(meltTransList[[2]])]))
colnames(meltTrophicMeans) <- c(names(meltTrophicMeans[-length(meltTrophicMeans)]),
                                "Var")
meltTrophicMeans$Number <- as.numeric(meltTrophicMeans$Number)
meltTrophicMeans$Var <- as.numeric(meltTrophicMeans$Var)

myDodge <- position_dodge(width = .2)



library(RColorBrewer)
par(mar = c(0, 4, 0, 0))
display.brewer.all()
nOfColorsQuad <- 5 # numer of quadrats
getPaletteQuad <- colorRampPalette(brewer.pal(5, "Set1"))


# plotter

meanPlot <- ggplot(data = meltTrophicMeans, aes(x = Year, y = Number, color = Group))+
  geom_line(position = myDodge)+
  geom_point(position = myDodge)+
  geom_errorbar(aes(ymax = Number + sqrt(Var),
                ymin = Number - sqrt(Var)),
                position = myDodge)+
  scale_color_manual(values = getPaletteQuad(nOfColorsQuad))+
  scale_x_continuous(breaks = seq(5,16,1),
                     labels = seq(5,16,1))+
  # scale_y_continuous(limits = c(0,350),
  #                    breaks = seq(0,350,50))+
  labs(y=expression(paste(Sponge~density~(n~m^2))))+
  theme_bw()+
  # theme(panel.grid.minor = element_blank(), 
  #       panel.grid.major = element_blank())+
  theme(plot.title = element_text(size=14, vjust=2))+
  theme(axis.title.x = element_text(size=10,vjust=-0.5),
        axis.title.y = element_text(size=10,vjust=0.5))+
  theme(axis.text.x=element_text(size=10, angle = 60, 
                                 hjust = 1, vjust = .9))+
  theme(axis.text.y=element_text(size=10))+
  facet_wrap(~Site, nrow = 1) 
meanPlot

ggsave("/home/somros/Documents/R/exploratoryHoga/output/pics/meanSites.pdf", meanPlot,
        width=6.5, height=4, useDingbats=T)


#######################################################################









# branch for fate of quadrats separately

transposeFrame <- as.data.frame(t(trophicMode))

#transposeFrame$Total <- rowSums(transposeFrame)

# need to separate years (my x variable) and quadrats

quadratID <- rownames(transposeFrame)
transposeFrame$Year <- as.numeric(substr(quadratID, 1, 2))
transposeFrame$Site <- substr(quadratID, nchar(quadratID)-1, nchar(quadratID)-1)
transposeFrame$Quadrat <- substr(quadratID, nchar(quadratID), nchar(quadratID))  
transposeFrame$QuadratSite <- substr(quadratID, nchar(quadratID)-1, nchar(quadratID))


# must melt this 

meltMyData <- melt(transposeFrame, id.vars = c("Year", "Site", "Quadrat", "QuadratSite"),
             variable.name = "Group", value.name = "Number")


# plotting region

library(RColorBrewer)
par(mar = c(0, 4, 0, 0))
display.brewer.all()
nOfColorsQuad <- 5 # numer of quadrats
getPaletteQuad <- colorRampPalette(brewer.pal(5, "Set1"))

# 

quadratsFate <- ggplot(data = meltMyData[meltMyData$Group == "p",],
                       aes(x = Year, y = Number, group = Quadrat, color = Quadrat))+
  geom_line()+
  geom_point()+
  scale_color_manual(values = getPaletteQuad(nOfColorsQuad))+
  scale_x_continuous(breaks = seq(5,16,1),
                     labels = seq(5,16,1),
                     limits = c(5,16))+
  # scale_y_continuous(limits = c(0,350),
  #                    breaks = seq(0,350,50))+
  labs(y=expression(paste(Sponge~density~(n~m^2))))+
  theme_bw()+
  # theme(panel.grid.minor = element_blank(), 
  #       panel.grid.major = element_blank())+
  theme(plot.title = element_text(size=14, vjust=2))+
  theme(axis.title.x = element_text(size=10,vjust=-0.5),
        axis.title.y = element_text(size=10,vjust=0.5))+
  theme(axis.text.x=element_text(size=10, angle = 60, 
                                 hjust = 1, vjust = .9))+
  theme(axis.text.y=element_text(size=10))+
  facet_grid(.~ Site, scales = "free")
quadratsFate

ggsave("/home/somros/Documents/R/exploratoryHoga/output/pics/quadratsFatePhot.pdf", quadratsFate,
        width=6.5, height=3, useDingbats=T)

# next thing to explore is the variability within and between sites

# geom_stat uses the predict routine using by default the LOESS method (LOcal regrESSion)
# which is a point-wise regression method. Note that the line is NOT an average. The 
# shaded area is 95% point-wise confidence interval

quadratsFateSmooth <- ggplot(data = meltMyData[meltMyData$Order== "Total",], aes(x = Year, y = Number))+
  geom_smooth()+
  geom_point(aes(shape = Quadrat))+
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
quadratsFateSmooth

# linear Nt Nt+1 plots

# trash data from 2010, 11, 12 to use only adjacent years

# calculate the number of previous years for all the years, also the not adjacent pairs, then trash them

dummyList <- list()
for (gr in 1:4) {
  dummyList[[gr]] <- c(rep(NA, nrow(transposeFrame[transposeFrame$Year == 5,])),
                       transposeFrame[,gr][(transposeFrame$Year != 16)])
}

prevYrFrame <- abind(dummyList, along = 2)
pairFrame <- as.data.frame(cbind(transposeFrame, prevYrFrame))
colnames(pairFrame) <- c(names(transposeFrame), "hPrevYr", "pPrevYr", "uPrevYr", "totPrevYr")
truePairs <- pairFrame[pairFrame$Year != 5 & pairFrame$Year != 11 & pairFrame$Year != 13,]

# create df for greys in the plots

greys <- truePairs[,c(1,2,3,4,9,10,11,12)]

library(RColorBrewer)
par(mar = c(0, 4, 0, 0))
display.brewer.all()
nOfColorsQuad <- 7 # numer of quadrats
getPaletteQuad <- colorRampPalette(brewer.pal(7, "Set1"))

# now plot this shit please

pairsPlot <- ggplot(data = truePairs,
                    aes(x = pPrevYr, y = p, color = factor(Year)))+
  geom_abline(intercept = 0, slope = 1, color = "grey90") +
  geom_point(data = greys, aes(x = pPrevYr, y = p), color = "grey90")+
  geom_point()+
  scale_colour_manual(values=getPaletteQuad(nOfColorsQuad))+
  theme_bw()+
  facet_grid(~ Site)
pairsPlot

ggsave("/home/somros/Documents/R/exploratoryHoga/output/pics/pairsPlotP.pdf", pairsPlot,
       width=6.5, height=2.3, useDingbats=T)
