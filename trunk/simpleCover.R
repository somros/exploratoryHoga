# 29.07.2016
# Alberto.Rovellini@vuw.ac.nz
# script to average the benthic % cover from all the stations in the Wakatobi from 2014 and plot
# it. At the moment only a simple barchart for the cover, can come up with something better though.
# next step is to include some formal statistic to tell the sites apart (although this must have been
# done already). At the moment, two separate twin scripts according to the desired degree of detail of the
# cover, but it must become a method of the same script

# edit 02.08.2016

# script has a flag to be set when running it to decide how to read the .xlsx files, whether
# with the detailed cover or not. next step is to set a flag to subset which datasets
# of the main directory must be read. ideally, in the future this should become a function
# connected to another script where the functions of interest are called.

require(XLConnect)
require(abind)
require(ggplot2)

setwd("/home/somros/Documents/R/exploratoryHoga/input/CPC_HOLLY")

# set flags and other settings for the script

flagIsComplete <- F # set TRUE or FALSE depending on the detail of the desired % cover
myLocations <- c("BUOY3", "SAMPELA") # type array of filenames with the locations of interest. 
#default is NULL and reads all the data. some flexibility allowed, location name is enough but
# it must be spelled right

myDepths <- c("_5", 
              "10")#, 
              #"15") # sets depths of interest, please mind the "_" before 5

listOfSheets <- list.files("/home/somros/Documents/R/exploratoryHoga/input/CPC_HOLLY",
                           pattern = ".xls", recursive = T)

# subsets the list of files to be read into the script basing on the myLocation

locationsOfInterest <- list() # initiate list
for (i in myLocations) {
  locationsOfInterest[[i]] <- subset(listOfSheets, grepl(i, listOfSheets))
}
locationsOfInterest <- unlist(locationsOfInterest)

# does the same for the depth

depthOfInterest <- list()
for (i in myDepths) {
  depthOfInterest[[i]] <- subset(locationsOfInterest, grepl(i, locationsOfInterest))
}

dataOfInterest <- unlist(depthOfInterest)

# if the dataOfInterest is NULL because of myLocation being NULL, the scripts loops over the whole 
# list of data

if (is.null(dataOfInterest)) {
  listOfSheets <- listOfSheets
} else {
  listOfSheets <- dataOfInterest
}


# reads the 2014 data sheets depending on the selected flag

spreadsheetReader <- function(spreadsheet, isComplete) {
  book <- loadWorkbook(spreadsheet)
  if (isComplete==T) {
    summaryData <- readWorksheet(book, "Data Summary", startRow = 27, endRow = 65,
                                 startCol = 1, endCol = 2)
  } else {
  summaryData <- readWorksheet(book, "Data Summary", startRow = 10, endRow = 23,
                               startCol = 1, endCol = 2)
  }
  summaryData$Location <- rep(substr(spreadsheet, 1, nchar(spreadsheet)-5), nrow(summaryData))
  colnames(summaryData) <- c("Category", "% cover", "Location")
  summaryData <- summaryData[complete.cases((summaryData)),]
  return(summaryData)
}

listOfSets <- lapply(listOfSheets, spreadsheetReader, flagIsComplete)

# now all the datasets are stored in a list. replicates from same depth of same site must be averaged
# rearrange list to list of one frame per depth per site with the three replicates as columns


# create list with the location names, i.e. site and depth but no replicate

listOfLocations <- levels(factor(substr(listOfSheets, 1, nchar(listOfSheets)-7)))

# create a list with the indices of each replicate in the original listOfSets

listOfIndices <- list()
for (i in listOfLocations) {
  listOfIndices[[i]] <- grep(i, listOfSets)
}

# creates a temporary list where the three replicates per site are stored in a nested sublist

tempList <- rep(list(list()), length(listOfLocations))

for (j in 1:length(listOfIndices)) {
  for (i in listOfIndices[[j]]) {
    tempList[[j]][[i]] <- listOfSets[[i]]
  }
}

# unites the lists of 3 replicates per site in one data frame per site, without dumping any column yet

sortedList <- lapply(tempList, function(x) abind(x, along=2)) # sweet

# create function to rearrange the new data frames in neat frames with three data columns, the 
# key to the coverage and the name of the site and depth

columnSorter <- function(x) {
  x <- as.data.frame(x)
  dummy <- as.character(x[,length(x)])
  x <- x[,c(1,2,5,8)]
  x$Location <- rep(substr(dummy, 1, nchar(dummy)-2))
  return(x)
}

# applies the above function to the sorted list

sites <- lapply(sortedList, columnSorter)

# now averages must be taken and standard error calculated

meanSDtool <- function(frame) {
  # defines the data columns
  dataColumns <- frame[,c(2,3,4)]
  # defines relevant functions:
  # 1_ standard error
  standardError <- function(x) sd(x)/sqrt(length(x))
  # 2_ turn all data columns to numeric
  toNumeric <- function(x) as.numeric(as.character(x))
  
  # turn data columns to numeric
  dataColumns <- apply(dataColumns, 2, toNumeric)
  mean <- rowMeans(dataColumns)
  se <- apply(dataColumns, 1, standardError)
  meanFrame <- data.frame(frame[,length(frame)], frame[,1], mean, se)
  colnames(meanFrame) <- c("Location", "Category", "Mean", "StdErr")
  return(meanFrame)
}

means <- lapply(sites, meanSDtool)

# good. now it makes sense to put it all together in one frame for ggplot. worth to go back to
# the start and reorder the factors though

tmpMeanFrame <- abind(means, along=1)
tmpMeanFrame <- gsub("_5", "_05", tmpMeanFrame)

meansData <- data.frame(as.factor(tmpMeanFrame[,1]),
                        as.factor(tmpMeanFrame[,2]),
                        as.numeric(tmpMeanFrame[,3]),
                        as.numeric(tmpMeanFrame[,4]))
colnames(meansData) <- c("Location", "Category", "Mean", "StdErr")

meansData <- meansData[order(meansData$Location),]


# fix levels of factor as per usual. how about learning to use attach someday?

meansData$Category <- factor(as.character(meansData$Category), levels = unique(meansData$Category))

# set factors of the location for the faceting of the plot
charFac <- as.character(meansData$Location)
meansData$facetFactors <- factor(substr(charFac, 
                              nchar(charFac)-1, 
                              nchar(charFac)))
meansData$facetFactors <- factor(paste(meansData$facetFactors, "m", sep = " "))
# plot, will need to become automated as function to plot specified locations to be entered at 
# the beginning. then again if I want to automate I should do it on reading a subset of the
# data and plot it all afterwards, instead of subsetting the whole thing after dragging it
# around without needing it for the whole script
# 

p <- ggplot(data=meansData, aes(x = Location, y = Mean, fill = Category))+
  geom_bar(stat = "identity", width = .7)+
  # geom_errorbar(data = buoy3Sampela, 
  #               aes(ymax = Mean + StdErr, ymin = Mean - StdErr),
  #               width = .7)+
  scale_fill_grey(start = 0, end = 0.95)+
  labs(y = "Average % cover")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank())+
  theme(plot.title = element_text(size=14, vjust=2))+
  theme(axis.title.x = element_text(size=10,vjust=-0.5),
        axis.title.y = element_text(size=10,vjust=0.5))+
  theme(axis.text.x=element_text(size=10, angle = 45, vjust = 0.5))+
  theme(axis.text.y=element_text(size=10))+
  facet_grid(. ~ facetFactors, scales = "free")
p

ggsave("/home/somros/Documents/R/exploratoryHoga/output/pics/benthicCover.pdf", p,
       width=8, height=6, useDingbats=T)






