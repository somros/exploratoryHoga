# 09.08.2016 Alberto.Rovellini@vuw.ac.nz
# this is a script to read the opwall monitoring data

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
setLocationName <- c("F", "C", "S")
setReplicate <- 1:3
namesOfSheets <- as.vector(sapply(setSiteName,
                               function(x) {c(paste(x, as.vector(sapply(setLocationName,
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
                                      startCol = 1, endCol = 8)
  }
  return(summaryData)
}

listOfBenthicTransectsTmp <- lapply(listOfSheets, monitBentReader)
# for some reason the function read 18 empty entries, need to investigate
listOfBenthicTransects <- lapply(listOfBenthicTransectsTmp, function(x) x[19:27]) # list of 2 lists of 9 dataframe each. each level 1 list is one location, each level 2 a transect

# might want to have a single list instead

benthicTransects <- c(listOfBenthicTransects[[1]], listOfBenthicTransects[[2]]) # merge into 1 list
benthicTransects <- lapply(benthicTransects, function(x) x[-c(5:8)]) # get rid of useless columns, need to specify this in the reading routine instead

##########################################

# section to correct all the misspelled entries, which are a lot. need to come up with a function, won't be 
# easy. Limit to the Benthic.Type now, I need to go on.

benthicTypeLevels <- levels(factor(unlist(lapply(benthicTransects, function(x) levels(factor(x$Benthic.Type))))))
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


# # criteria must be composite though, as Morphology and Further.Info are restricted to corals
# # one solution is to substitute all 0 with the most detailed info available
# 
# # not working, and also not efficient if it was, please learn how to work in R
# 
# criteriaRefiner <- function(frame) {
#   for (i in 1:nrow(frame)) {
#     if (frame$Morphology[i] == "0") {
#       frame$Morphology[i] <- frame$Benthic.Type[i]
#     }
#     if (frame$Further.Info[i] == "0") {
#       frame$Further.Info[i] <- frame$Morphology[i]
#     }
#   }
#   return(frame)
# }
# 
# # testFrame <- benthicTransects[[2]]
# # for (i in 1:nrow(testFrame)) {
# #   if (testFrame$Morphology[i] == "0") {
# #     testFrame$Morphology[i] <- testFrame$Benthic.Type[i]
# #   }
# #   if (testFrame$Further.Info[i] == "0") {
# #     testFrame$Further.Info[i] <- testFrame$Morphology[i]
# #   }
# # }
# 
# newBT <- lapply(benthicTransects, criteriaRefiner)
# 
# 
# # DIO PORCO SCANNATO IN CROCE COI CHIODI NELLE MANI E LA MADONNA CHE RIDE
# 
# 
# # i must keep going I have no time 


# to be applied to all frames in the list

percentCoverCalc <- function(frameReplicate) {
  subsetFrame <- frameReplicate[,1:2]
  nOfPoints <- as.data.frame(colSums(table(subsetFrame)))
  nOfPoints$Benthic.Type <- rownames(nOfPoints)
  nOfPoints$Percentage.Cover <- nOfPoints[,1]*100/sum(nOfPoints[,1])
  rownames(nOfPoints) <- 1:nrow(nOfPoints)
  nOfPoints <- nOfPoints[,c(2,1,3)]
  colnames(nOfPoints) <- c("Type", "Points", "Cover")
  return(nOfPoints)
}

pointsAndCover <- lapply(benthicTransectsCorrect, percentCoverCalc)

# need to have all frames with the same levels. they have to be the sum of all the available levels

# now it is a list of data frames of different length and different levels. Empty levels must be
# added as new lines to each frame, with 0 percentage (don't even want to think about the nightmare of
# going through this for separate years and probably different observers and compilers)

dummyFrame <- as.data.frame(cbind(newLevels, rep(0, length(newLevels)), rep(0, length(newLevels))))
colnames(dummyFrame) <- names(pointsAndCover[[1]])
dummyFrame

# aggregate the levels

pointsAndCoverComplete <- lapply(pointsAndCover, function (x) {
  completeLevels <- as.data.frame(rbind(x, dummyFrame))
  colnames(dummyFrame) <- names(dummyFrame)
  completeLevels$Points <- as.numeric(as.character(completeLevels$Points))
  completeLevels$Cover <- as.numeric(as.character(completeLevels$Cover))
  oneEntryList <- split(completeLevels, completeLevels$Type)
  oneEntryListAgg <- lapply(oneEntryList, function(y) {
    z <- c(levels(factor(y[,1])), colSums(y[,2:3]))
    return(z)}
  )
  oneEntry <- as.data.frame(abind(oneEntryListAgg, along = 0))
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
# the script so far has not done anything

benthicData <- as.data.frame(abind(pointsAndCoverComplete, along = 1))
benthicData$Points <- as.numeric(as.character(benthicData$Points))
benthicData$Cover <- as.numeric(as.character(benthicData$Cover))
rownames(benthicData) <- 1:nrow(benthicData)

# now can do summary statistics over the replicate

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
meanAndSdTmp <- lapply(meanAndSdList, function(x) as.data.frame(abind(x, along = 1)))
meanAndSd <- as.data.frame(abind(meanAndSdTmp, along = 1))
rownames(meanAndSd) <- 1:nrow(meanAndSd)
colnames(meanAndSd) <- c("Type", "Mean", "Sd", "Location")
meanAndSd$Mean <- as.numeric(as.character(meanAndSd$Mean))
meanAndSd$Sd <- as.numeric(as.character(meanAndSd$Sd))
# order factors of stations for plot


# plot

library(RColorBrewer)
par(mar = c(0, 4, 0, 0))
display.brewer.all()
nOfColors <- length(levels(meanAndSd$Type))
getPalette <- colorRampPalette(brewer.pal(11, "BrBG"))
#myPalette <- doublePalette[seq(3,length(doublePalette),1)]

benthicMonitoring <- ggplot(meanAndSd, aes(x=Location, y=Mean, fill=Type))+
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
