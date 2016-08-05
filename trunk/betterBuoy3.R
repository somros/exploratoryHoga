# 27/07/2016 script to open and explore the abundance data from Buoy3
# Alberto.Rovellini@vuw.ac.nz

# 02.08.2016 script reads from .csv because .xlsx is very messy. Goal now is to analyse a bit better the available
# data, such as the most abundant species and so on, limit the analysis to those species.

require(abind)
require(plyr)
require(ggplot2)
require(reshape2)

dataAllYears <- read.csv("/home/somros/Documents/R/exploratoryHoga/input/spongeAbundanceQuadrats.csv")

# set flags

flagRoutine <- "average" # can be either "sum" to add the 5 quadrats up, or any other string to take average and sd

# initiate information about the dataset

species <- read.csv("/home/somros/Documents/R/exploratoryHoga/input/speciesKey.csv")
species[species=="" | species==0] <- NA # drop all that is not a name, to be refined though

# option 1: keep only the species (or at least what is close)

identified <- species[,-ncol(species)]
identified <- identified[complete.cases(identified),]

# option 2: use species when available and description otherwise

speciesOrDescription <- list()

for (i in 1:nrow(species)) {
  if (is.na(species[,2][i])==T) {
    speciesOrDescription[[i]] <- species[,3][i]
  } else {
    speciesOrDescription[[i]] <- species[,2][i]
  }
}
speciesOrDescription <- unlist(speciesOrDescription) #

# follows general information about the dataset, such as number of years, quadrats, sites and so on

years <- levels(factor(substr(names(dataAllYears[-1]),2,3)))
quadrats <- as.numeric(levels(factor(substr(names(dataAllYears), nchar(names(dataAllYears)), nchar(names(dataAllYears))))))
quadrats <- quadrats[is.na(quadrats)==F]# number of quadrats per site
dataColumns <- ncol(dataAllYears)-1 # number of quadrats in all sites in all years, effectivelty number of columns in the frame
sites <- c("A", "B", "C")

# each data column is a quadrat. There are five quadrats per site. first letter of the header is X because 
# read.csv() (and read.table() for that matter) do not allow a column name to start wit a number. 
# first 2 digits after X are the year (05 to 16). S is site, and second to last letter (A, B or C) is the site code.
# final digit is the quadrat within the site.


# all the columns need to be turned to numeric. some typos in the excel spreadsheet might cause the columns to be
# factors or characters. as the method removes everything that cannot be coerced to numeric, typos such as "19'"
# will also be dropped, hence it's recommended to take a look at the data first.
# in an ideal world I will adapt the script to skim through typos the right way

# cut the whole thing to the first 125 rows, as below that it gets confused

dataAllYears <- dataAllYears[1:124,]

for (i in 1:ncol(dataAllYears)) {
  if (is.numeric(dataAllYears[,i])==T) {
    dataAllYears[,i] <- dataAllYears[,i]
  } else {
    dataAllYears[,i] <- as.numeric(levels(dataAllYears[,i])[dataAllYears[,i]]) # this removes everything that won't fit as numeric
  }
}

dataAllYears[is.na(dataAllYears)] <- 0 # turns NAs to zeroes

# gets rid of the X in front of the column names because it bothers me. but then again it puts it back

colnames(dataAllYears) <- c(names(dataAllYears)[1],
                                  substr(names(dataAllYears[,-1]), 2, nchar(names(dataAllYears[,-1]))))


# yet another strategy, STOLEN from SO

dataAllYears <- dataAllYears[,-1]

# please understand and rename

n <- 1:ncol(dataAllYears)
ind <- data.frame(matrix(c(n, rep(NA, 5 - ncol(dataAllYears)%%5)), byrow=F, nrow = 5))
nonna <- sapply(ind, function(x) all(!is.na(x)))
ind <- ind[, nonna]

# if statement to decide whether we want to add the quadrats or take their average

if (flagRoutine == "sum") {
  electedFunction <- function(i)rowSums(dataAllYears[,i])
} else {
  electedFunction <- function(i) {
    meanColumn <- rowMeans(dataAllYears[,i])
    varColumn <- apply(dataAllYears[,i], 1, sd)
    cbind(meanColumn, (varColumn)^2) # to get the variance, used for the propagation
  }
}

sitesFrame <- as.data.frame(do.call(cbind, lapply(ind, electedFunction))) # hahah nailed it
sitesFrame$Species <- speciesOrDescription

# reorganize the dataframe with first the means and then the sd

sitesFrame <- cbind(sitesFrame[,grepl("mean", names(sitesFrame))], 
                    sitesFrame[,grepl("V", names(sitesFrame))])


################################################################################################################

# need to know the most abundant species per column, for loop because it's easier and I'm short on time
# sd is not absolutely necessary here

sortedMostAbundant <- list()
for (i in 1:length(sitesFrame)) {
  sortedMostAbundant[[i]] <- sitesFrame[order(sitesFrame[,i], decreasing = T),]
  sortedMostAbundant[[i]] <- sortedMostAbundant[[i]][,c(length(sortedMostAbundant[[i]]),i)]
}

sortedMostAbundant <- sortedMostAbundant[-length(sortedMostAbundant)]

# can decide that dominant species together account for 90% of the numbers

chainsaw <- function(x) {
  # first it needs to calculate the sum of the numbers
  sumOfNumbers <- sum(x[,2])
  # then calculate n as in n = sumOfNumbers*domPerc / 100
  n = sumOfNumbers*90/100
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


dominantPerc <- lapply(sortedMostAbundant, chainsaw)


# from here the most abundant species can be extracted. I hoped the community to be a bit less diverse tbh,
# in some cases about 20 species account for 90% of the diversity (which is not even that big of a percentage to consider)

###########################################################################################################

# resuming back from the list of data frames per site (listOfSites)
# need to average 3 by 3. Same strategy as above I reckon

# here I must calculate the grand mean and the grand sd. ominous
# for the purpose, I split the original frame into 2 though


sitesFrameMeans <- sitesFrame[,grepl("mean", names(sitesFrame))]
sitesFrameVar <- sitesFrame[,grepl("V", names(sitesFrame))]


m <- 1:ncol(sitesFrameMeans)
indYears <- data.frame(matrix(c(m, rep(NA, 3 - ncol(sitesFrameMeans)%%3)), byrow=F, nrow = 3))
noNAYears <- sapply(indYears, function(x) all(!is.na(x)))
indYears <- indYears[, noNAYears]

yearsFrame <- as.data.frame(do.call(cbind, lapply(indYears, function(i) rowMeans(sitesFrameMeans[,i])))) # it seems to be working although I don't know how
#yearsFrame$Species <- speciesOrDescription
#colnames(yearsFrame) <- c(years, "Species")

# define function to calculate the combined sd according to
# GSD = sqrt((ESS+TGSS)/N-1)
# ESS <- (var1*dof1) + (var2*dof2) + ... + (varN*dofN) # where n is the group, i.e. the site A B C, i.e. the column
# TGSS <- (µ1-GM)^2 * n1 + (µ2-GM)^2 * n2 + ... + (µN-GM)^2 * nN

yearsFrameCheat <- cbind(yearsFrame, yearsFrame, yearsFrame)
columnIndex <- c(1,11,21,2,12,22,3,13,23,4,14,24,5,15,25,6,16,26,7,17,27,8,18,28,9,19,29,10,20,30) # please don't

yearsFrameCheat <- yearsFrameCheat[,columnIndex] 

GSDcalculator <- function(x) {
  ESS <- rowSums(sitesFrameVar[,x]*(length(quadrats)-1)) # so far so good, don't chack again
  TGSS <- rowSums((sitesFrameMeans[,x]-yearsFrameCheat[,x])^2*length(quadrats)) # faulted SON OF A BITCH WILL MESS THE COLUMNS
  GSD <- sqrt((ESS+TGSS)/(length(years)-1))
  return(GSD)  
}

GSDframe <- as.data.frame(do.call(cbind, lapply(indYears, GSDcalculator))) 

# now put species, means and sd together

yearsFrame$Species <- speciesOrDescription
colnames(yearsFrame) <- c(years, "Species")

GSDframe$Species <- speciesOrDescription
colnames(GSDframe) <- c(years, "Species")

meltYears <- melt(yearsFrame, id.vars = "Species", variable.name = "Year", value.name = "Mean" )
meltYears$Year <- as.numeric(as.character(meltYears$Year))
meltGSD <- melt(GSDframe, id.vars = "Species", variable.name = "Year", value.name = "GSD" )

meltComplete <- cbind(meltYears, meltGSD$GSD)

# area chart

area <- ggplot(data = meltComplete, aes(x = Year, y = Mean, group = Species, fill = Species))+
  geom_area()+
  guides(fill = F)+
  scale_x_continuous(breaks = seq(5,16,1),
                     labels = seq(5,16,1),
                     limits = c(5,16))+
  scale_y_continuous(limits = c(0,200),
                     breaks = seq(0,200,50))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank())+
  theme(plot.title = element_text(size=14, vjust=2))+
  theme(axis.title.x = element_text(size=10,vjust=-0.5),
        axis.title.y = element_text(size=10,vjust=0.5))+
  theme(axis.text.x=element_text(size=10))+
  theme(axis.text.y=element_text(size=10))
area

# barchart

bar <- ggplot(data = meltComplete, aes(x = Year, y = Mean, group = Species, fill = Species))+
  geom_bar(stat = "identity")+
  guides(fill = F)
bar


# all the possible entries are plotted. script needs to be adapted to include only the desired species, i.e. to
# restrict the analysis to either the most abundant or the complete cases. Keep in mind that this is only the number




??melt



