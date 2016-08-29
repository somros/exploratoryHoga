# 09.08.2016 Alberto.Rovellini@vuw.ac.nz
# this is a script to read the opwall monitoring data

# 11.08.2016 Script is a bit better now. Next thing to do is to order the factors of the location
# for the plotting, but the renaming of the empty benthic types is working now and the code is
# better commented

# rename in a way that all object in the project have unique names please

require(XLConnect)
require(abind)
require(plyr)
require(ggplot2)

setwd("/home/somros/Documents/Data/Hoga/MonitoringProgram")

# set flags (later)

listOfSheets <- list.files("/home/somros/Documents/Data/Hoga/MonitoringProgram",
                           pattern = "Benthic", recursive = T)

setSiteName <- c("B3", "S1")
setDepthName <- c("F", "C", "S")
setReplicate <- 1:3
namesOfSheets <- as.vector(sapply(setSiteName,
                               function(x) {c(paste(x, as.vector(sapply(setDepthName,
                                                                        function(x) {c(paste(".", x, setReplicate, sep = ""))})), 
                                                    sep = ""))}))
namesOfLocations <- substr(namesOfSheets, 1, nchar(namesOfSheets)-1)
namesOfLocations <- levels(factor(namesOfLocations, levels = unique(namesOfLocations)))

summaryData <- vector(mode = "list", length = length(namesOfSheets))

monitBentReader <- function(spreadsheet) {
  book <- loadWorkbook(spreadsheet)
  bookSheets <- getSheets(book)
  for (i in bookSheets) {
    summaryData[[i]] <- readWorksheet(book, i, startRow = 1, endRow = 202,
                                      startCol = 1, endCol = 4)
  }
  return(summaryData)
}

listOfBenthicTransectsTmp <- lapply(listOfSheets, monitBentReader)
# for some reason the function read 18 empty entries, need to investigate
listOfBenthicTransects <- lapply(listOfBenthicTransectsTmp, function(x) x[19:27]) # list of 2 lists of 9 dataframe each. each level 1 list is one location, each level 2 a transect

# lump all lists into one single list

benthicTransects <- listOfBenthicTransects[[1]]
for (i in 2:length(listOfBenthicTransects)) {
  benthicTransects <- append(benthicTransects, listOfBenthicTransects[[i]])
} 

##########################################

# section to correct all the misspelled entries, which are a lot. need to come up with a function, won't be 
# easy. Limit to the Benthic.Type for now. This requires some analysis of the typos

benthicTypeLevels <- levels(factor(unlist(lapply(benthicTransects, function(x) levels(factor(x$Benthic.Type))))))
benthicTypeLevels

# correction routine, manual specification of the faulty entries

benthicTransectsCorrect <- lapply(benthicTransects, function(x) {
  x[x=="algae"] <- "Algae"
  x[x=="rubble"] <- "Rubble"
  x[x=="sand"] <- "Sand"
  x[x=="silt"] <- "Silt"
  x[x=="sponge"] <- "Sponge"
  x[x=="Sponge "] <- "Sponge"
  return(x)
})

newLevels <- levels(factor(unlist(lapply(benthicTransectsCorrect, function(x) levels(factor(x$Benthic.Type))))))


#************************************************************************************************#

# routine to substitute the missing entries in Morphology with Benthic.Type and 
# of Further.Info with Morphology. Aim is to have entries for all the levels of details.

# for reasons of the function it's necessary to get rid of the NAs, turn them to 0

benthicTransectsCorrect <- lapply(benthicTransectsCorrect, function(x) {
  x[is.na(x)] <- "0"
  return(x)}
)

# function to substitute the relevant entries

groupsRewriter <- function(frameTransect) {
  frameTransect1 <- within(frameTransect, Morphology[Morphology=="0"] <- Benthic.Type[Morphology=="0"])
  frameTransect2 <- within(frameTransect1, Further.Info[Further.Info=="0"] <- Morphology[Further.Info=="0"])
  return(frameTransect2)
}

benthicTransectsComplete <- lapply(benthicTransectsCorrect, groupsRewriter) # yep

################################################################################################

# routine to add column with coarse benthic type

levels(factor(benthicTransectsComplete[[1]]$Benthic.Type))

# function: take the benthic type out as vector, replace with coarse entry with similar function as above,
# add column to the data frame with the new benthic type. applied to all frames in the list
# recode the rest with a dynamic column identifier instead of benthic type to pick the desired aggregation
# method



# routine to calculate the percentage cover from the tape points per category. Column IDs should be
# specified outside the function, at the beginning of the script as flags. Flags will have to become 
# function arguments if the set of scripts has to be turned into a package at some point, which it should

percentCoverCalc <- function(frameReplicate) {
  subsetFrame <- frameReplicate[,1:2]
  nOfPoints <- as.data.frame(colSums(table(subsetFrame)))
  nOfPoints$Benthic.Type <- rownames(nOfPoints)
  nOfPoints$Percentage.Cover <- nOfPoints[,1]*100/sum(nOfPoints[,1])
  rownames(nOfPoints) <- 1:nrow(nOfPoints)
  nOfPoints <- nOfPoints[,c(2,1,3)] # reorders the columns with type, % cover and points
  colnames(nOfPoints) <- c("Type", "Points", "Cover")
  return(nOfPoints)
}

pointsAndCover <- lapply(benthicTransectsComplete, percentCoverCalc)

# need to have all frames with the same levels. they have to be the sum of all the available levels

# first build a dummy frame with all the levels and 0 as entries for cover and points
dummyFrame <- as.data.frame(cbind(newLevels, rep(0, length(newLevels)), rep(0, length(newLevels))))
colnames(dummyFrame) <- names(pointsAndCover[[1]]) # rename the columns for consistency

# then append the dummy data frame at the end of each transect data frame. A bit of reorganization,
# column renaming and class manipulation is in the routine too to keep things as smooth as possible

pointsAndCoverComplete <- lapply(pointsAndCover, function (x) {
  completeLevels <- as.data.frame(rbind(x, dummyFrame)) # append dummy frame
  colnames(dummyFrame) <- names(dummyFrame) # rename columns as step above changes the names
  completeLevels$Points <- as.numeric(as.character(completeLevels$Points)) # redefine class to numeric
  completeLevels$Cover <- as.numeric(as.character(completeLevels$Cover)) # redefine class to numeric
  
  # split-apply-combine routine follows
  
  oneEntryList <- split(completeLevels, completeLevels$Type) # splits on the benthic type
  oneEntryListAgg <- lapply(oneEntryList, function(y) {
    z <- c(levels(factor(y[,1])), colSums(y[,2:3])) # sums points and cover per benthic type
    return(z)}
  )
  oneEntry <- as.data.frame(abind(oneEntryListAgg, along = 0)) # recombines all of it into one frame
  
  # polishing follows
  
  colnames(oneEntry) <- names(dummyFrame) 
  oneEntry$Points <- as.numeric(as.character(oneEntry$Points))
  oneEntry$Cover <- as.numeric(as.character(oneEntry$Cover))
  return(oneEntry)
})

# adds a column with the corresponding replicate name to each dataframe

for (i in 1:length(pointsAndCoverComplete)) {
  pointsAndCoverComplete[[i]]$Replicate <- rep(names(pointsAndCoverComplete[i]), 
                                               nrow(pointsAndCoverComplete[[i]]))
  pointsAndCoverComplete[[i]]$Location <- substr(pointsAndCoverComplete[[i]]$Replicate,
                                                 1, nchar(pointsAndCoverComplete[[i]]$Replicate)-1)
}

# now the transect frames are correctly named, organized and the %cover is calculated.
# all the dataframes are merged together by row, renaming and reclassing follows

benthicData <- as.data.frame(abind(pointsAndCoverComplete, along = 1))
benthicData$Points <- as.numeric(as.character(benthicData$Points))
benthicData$Cover <- as.numeric(as.character(benthicData$Cover))
rownames(benthicData) <- 1:nrow(benthicData)

# now can do summary statistics over the replicate, mean and sd or whatever else. 
# recursive split-apply-combine on the benthic type and the location

typeSplit <- split(benthicData, benthicData$Type)
locationSplit <- lapply(typeSplit, function(x) split(x, x$Location))
meanAndSdList <- lapply(locationSplit, function(x) {
  lapply(x, function(y) {
    means <- mean(y$Cover)
    stdev <- sd(y$Cover)
    meanSd <- data.frame(levels(factor(y$Type)), means, stdev, levels(factor(y$Location)))
    return(meanSd)
  })
})

# nested lists on the second dimension are remerged first on the benthic type factor...

meanAndSdTmp <- lapply(meanAndSdList, function(x) as.data.frame(abind(x, along = 1)))

# ... and then on the location factor. Usual renaming and reclassing follows

meanAndSd <- as.data.frame(abind(meanAndSdTmp, along = 1))
rownames(meanAndSd) <- 1:nrow(meanAndSd)
colnames(meanAndSd) <- c("Type", "Mean", "Sd", "Location")
meanAndSd$Mean <- as.numeric(as.character(meanAndSd$Mean))
meanAndSd$Sd <- as.numeric(as.character(meanAndSd$Sd))

# order factors of stations for plot: two steps required:
# first order the data frame according to the indices of the namesOfLocations object

meanAndSd <- meanAndSd[order(match(meanAndSd$Location, namesOfLocations)),] # easy peasy

# then assign unique values to the levels of the Location column to keep the order for the plot

meanAndSd$Location <- factor(meanAndSd$Location, levels = unique(meanAndSd$Location))

# abiotic types can be lumped into one single type. However, to do that I'd wait to see other datasets

# create column for coarse benthic type

levels(meanAndSd$Type)

coarse <- c("Abiotic", "Algae", "Ascidian", "Abiotic", "Abiotic", "Hard coral",
            "Other", "Other", "Abiotic", "Abiotic", "Abiotic", "Abiotic", "Soft coral",
            "Sponge", "Seagrass", "Unknown", "Abiotic")
meanAndSd$Coarse <- as.factor(rep(coarse, length(levels(meanAndSd$Location))))

# need to add the means or do this early on, not working this way

# plot

library(RColorBrewer)
par(mar = c(0, 4, 0, 0))
#display.brewer.all()
nOfColors <- length(levels(meanAndSd$Coarse))
getPalette <- colorRampPalette(brewer.pal(9, "BrBG"))
#myPalette <- doublePalette[seq(3,length(doublePalette),1)]

benthicMonitoring <- ggplot(meanAndSd, aes(x=Location, y=Mean, fill=Coarse))+
  geom_bar(stat = "identity", width = .7)+
  # geom_errorbar(data = buoy3Sampela, 
  #               aes(ymax = Mean + StdErr, ymin = Mean - StdErr),
  #               width = .7)+
  #scale_fill_grey(start = 0, end = 0.95)+
  scale_fill_manual(values = getPalette(nOfColors))+
  labs(y = "Average % cover")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank())+
  theme(plot.title = element_text(size=14, vjust=2))+
  theme(axis.title.x = element_text(size=10,vjust=-0.5),
        axis.title.y = element_text(size=10,vjust=0.5))+
  theme(axis.text.x=element_text(size=10, angle = 45, 
                                 hjust = 1, vjust = .9))+
  theme(axis.text.y=element_text(size=10))
benthicMonitoring

ggsave("/home/somros/Documents/R/exploratoryHoga/output/pics/benthicMonitoring.pdf", benthicMonitoring,
       width=10, height=4, useDingbats=T)
