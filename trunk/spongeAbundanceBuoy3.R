# 27/07/2016 script to open and explore the abundance data from Buoy3
# Alberto.Rovellini@vuw.ac.nz

require(abind)
require(plyr)
require(ggplot2)
require(reshape2)

data <- read.csv("/home/somros/Documents/R/exploratoryHoga/input/spongeAbundanceQuadrats.csv")

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

years <- as.numeric(levels(factor(substr(names(data),2,3))))
years <- years[is.na(years)==F] # getting rid of the NA and being left with only the years
quadrats <- as.numeric(levels(factor(substr(names(data), nchar(names(data)), nchar(names(data))))))
quadrats <- quadrats[is.na(quadrats)==F]# number of quadrats per site
dataColumns <- ncol(data)-1 # number of quadrats in all sites in all years, effectivelty number of columns in the frame


# each data column is a quadrat. There are five quadrats per site. first letter of the header is X because 
# read.csv() (and read.table() for that matter) do not allow a column name to start wit a number. 
# first 2 digits after X are the year (05 to 16). S is site, and second to last letter (A, B or C) is the site code.
# final digit is the quadrat within the site.

# a way to go is to add the quadrats and average between sites.


# all the columns need to be turned to numeric. some typos in the excel spreadsheet might cause the columns to be
# factors or characters. as the method removes everything that cannot be coerced to numeric, typos such as "19'"
# will also be dropped, hence it's recommended to take a look at the data first.
# in an ideal world I will adapt the script to skim through typos the right way

# cut the whole thing to the first 125 rows, as below that it gets confused

data <- data[1:124,]

for (i in 1:ncol(data)) {
  if (is.numeric(data[,i])==T) {
    data[,i] <- data[,i]
  } else {
    data[,i] <- as.numeric(levels(data[,i])[data[,i]]) # this removes everything that won't fit as numeric
  }
}

data[is.na(data)] <- 0 # turns NAs to zeroes

# the following is completely unacceptable but I need to keep going
# objective of the following code is to create a list of data frames. each element of the list is opne frame for
# one of the sites at one of the years. column with the species id is not necessary

indices <- list() # initiate empty list for the loop

for (i in 1:(dataColumns/quadrats)) {
  indices[[i]] <- (1:5)+5*i-4
}

# next is to create the actual list of frames

listColumns <- list()

for (i in 1:length(indices)) {
  listColumns[[i]] <- data.frame(data[,indices[[i]]])
}

quadratsSums <- lapply(listColumns, rowSums) # here the quadrats for one site are added up. can also be averaged

siteFrame <- abind(quadratsSums, along=2) # back to a data frame structure, now each column is one site

# now sites need to be averaged three by three. same notion as before, now mean and sd though (to become se)

indicesSites <- list()

for (i in 1:(dataColumns/quadrats/length(sites))) {
  indicesSites[[i]] <- (1:3)+3*i-3
}

listMeans <- list()

for (i in 1:length(indicesSites)) {
  listMeans[[i]] <- data.frame(siteFrame[,indicesSites[[i]]])
}

# and then average. all of this was probably supposed to be done in one line

siteMeans <- lapply(listMeans, rowMeans)
finalFrame <- as.data.frame(abind(siteMeans, along=2))
colnames(finalFrame) <- years
head(finalFrame)

# following option 1, I subset the finalFrame with only the species keys

restrictedFrame <- finalFrame[identified[,1],] # yep
restrictedFrame$Species <- identified$Species.name
meltedFrame <- melt(restrictedFrame, 
                    variable.name = "Year",
                    value.name = "Number", id.vars = "Species")
meltedFrame$Year <- as.numeric(as.character(meltedFrame$Year))

# plot as area chart the abundance through time. 2 years (2010 and 2012) are missing from the data.
# also a calculation of the uncertainty is required.

p <- ggplot(data = meltedFrame, aes(x = Year, y = Number, group = Species))+
  geom_area(aes(fill = Species))
p

ggsave("/home/somros/Documents/R/exploratoryHoga/output/pics/spongeAbundanceBuoy3.pdf", p,
       width=12, height=7, useDingbats=T)

# species need to be aggregated until functional group level, then analysis can be rerun.
