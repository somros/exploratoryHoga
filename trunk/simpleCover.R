# script to extract the percentage cover data from the 2014 datasets for Buoy3 and Sampela (in fact
# ideally for all the stations) and then plot them and compare them. If I were to do a proper job, this
# should be skimming through the crude folder and extract the data I need from each excel file.
# worth a try but don't go nuts, this is just exploratory.

require(XLConnect)
require(abind)
require(ggplot2)

setwd("/home/somros/Documents/R/exploratoryHoga/input/CPC_HOLLY")

listOfSheets <- list.files("/home/somros/Documents/R/exploratoryHoga/input/CPC_HOLLY",
                           pattern = ".xls", recursive = T)


spreadsheetReader <- function(spreadsheet) {
  book <- loadWorkbook(spreadsheet)
  summaryData <- readWorksheet(book, "Data Summary", startRow = 10, endRow = 24,
                               startCol = 1, endCol = 2)
  summaryData$Location <- rep(substr(spreadsheet, 1, nchar(spreadsheet)-5), nrow(summaryData))
  colnames(summaryData) <- c("Category", "% cover", "Location")
  summaryData <- summaryData[complete.cases((summaryData)),]
  return(summaryData)
}

listOfSets <- lapply(listOfSheets, spreadsheetReader)

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

# substitute _5 with _05


# subset to buoy 3 only

buoy3 <- subset(meansData, grepl("BUOY3", meansData$Location))

# and buoy3 vs sampela at common depths only

buoy3Sampela <- subset(meansData, grepl("BUOY3", meansData$Location) & meansData$Location!="BUOY3_15" | 
                         grepl("SAMPELA", meansData$Location))

p <- ggplot(data=buoy3Sampela, aes(x = Location, y = Mean, fill = Category))+
  geom_bar(stat = "identity", position = "dodge", width = .7)+
  geom_errorbar(data = buoy3Sampela, 
                aes(ymax = Mean + StdErr, ymin = Mean - StdErr),
                position = "dodge", width = .7)+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank())+
  theme(plot.title = element_text(size=14, vjust=2))+
  theme(axis.title.x = element_text(size=10,vjust=-0.5),
        axis.title.y = element_text(size=10,vjust=0.5))+
  theme(axis.text.x=element_text(size=10))+
  theme(axis.text.y=element_text(size=10))
p

ggsave("/home/somros/Documents/R/exploratoryHoga/output/pics/benthicCover.pdf", p,
       width=8, height=6, useDingbats=T)






