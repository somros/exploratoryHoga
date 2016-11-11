# Sandbox script for a mixture of things, collecting MDS plots, quadrats fate and possibly also the growth
# model at some point. Most recent script, should be of reference for current analysis.

require(abind)
require(MASS)
require(reshape2)
require(ggplot2)


setwd("//Staff/Home/SCIFAC/rovellal/DocumentsRedir/Data/Hoga/buoy3")
myData <- read.csv("spongeAbundanceQuadratsCleaned.csv")

# turn to numeric immediately

#myData <- myData[c(1:124),] keep only for the uncleaned set here
myData <- apply(myData, 2, as.numeric)
myData[is.na(myData)] <- 0 # turns NAs in 0, as necessary for Bray-Curtis
myData <- as.data.frame(myData)

# definition of the species

species <- read.csv("//Staff/Home/SCIFAC/rovellal/DocumentsRedir/Data/Hoga/buoy3/speciesKeyCleaned.csv")
species[species=="" | species==0] <- NA # drop all that is not a name, to be refined though

speciesOrDescription <- list()

for (i in 1:nrow(species)) {
  if (is.na(species[,2][i])==T) {
    speciesOrDescription[[i]] <- species[,3][i]
  } else {
    speciesOrDescription[[i]] <- species[,2][i]
  }
}
speciesOrDescription <- unlist(speciesOrDescription)


# some steps to break down the frame in a list of frames, one per year per site.

nOfQuadrats <- 5

myData$Species <- speciesOrDescription

# TODO: introduce routine to cut data frames to the species accounting for 90% of the abundance. must be applied
# to each quadrat and then all the species considered make the subset key for the entire data frame

# we do it with apply. jk lol use a loop apply is slow af anyway

# perc <- 90
# 
# species <- vector("list", length = ncol(myData)-1)
# 
# for (k in 2:ncol(myData)) {
#   dataOrd <- myData[order(myData[,k], decreasing = T),]
#   # does the math: sum all the entries in one column to get the total
#   totNum <- sum(dataOrd[,k])
#   n <- totNum*perc/100
#   ntmp <- 0
#   spy <- 0
#   for (i in 1:nrow(dataOrd)) {
#     if (ntmp < n) {
#       ntmp <- ntmp + dataOrd[i,k]
#       spy <- i
#     } else {
#       ntmp <- ntmp
#       spy <- spy
#     }
#   }
#   species[[k]] <- dataOrd[c(1:spy),1]
# }
# 
# # NOTE: using the names of the taxonomic entries instead of numerical indexing means that the levels of the former are fewer
# # than those of the latter, as the data frame has duplicate entries in the indentification column
# 
# species <- factor(levels(factor(unlist(lapply(species, as.character)))))

# now subset original data

# subsetData <- myData[myData[,1] %in% species,] # there we go




# if for whatever reason we want to use the entire set of species for the analysis, we can get rid of the above
# section (this should become a flag). 
# the above routine may be better off if used on the means. That would rule out all the outliers from the
# single quadrats. The number of species is however half the original amount.


dataQuad <- myData[,-1] # remove species for convenience


# create vector of indeces for consequent loop (breakage)

indeces <- rep(1:30, each = 5)

dataQuadN <- rbind(indeces, dataQuad)

# split into list depending on the index

newList <- vector(mode = "list", length = length(indeces)/nOfQuadrats)

for (i in 1:length(newList)) {
  subsetFrame <- dataQuadN[1,]==i
  newList[[i]] <- dataQuadN[,subsetFrame]
}


# get rid of the index and make row means (sd later)

newListMeans <- lapply(newList, function(x) {
  y <- x[-1,] # gets rid of the "indeces" row used for the splitting
  z <- rowMeans(y)
  frame <- as.data.frame(z[1:nrow(dataQuad)])
  colnames(frame) <- paste("Mean", substr(names(x[1,])[1], 2,3), substr(names(x[1,])[1], 5, 5), sep = "")
  return(frame)
})


allMeansTMP <- as.data.frame(abind(newListMeans, along = 2))

allMeans <- cbind(speciesOrDescription, allMeansTMP)







####################  testing phase  ######################

# add the subsetting routine here and see how many species we manage to cut

perc <- 90

speciesSet <- vector("list", length = ncol(myData)-1)

for (k in 2:ncol(allMeans)) {
  dataOrd <- allMeans[order(allMeans[,k], decreasing = T),]
  # does the math: sum all the entries in one column to get the total
  totNum <- sum(dataOrd[,k])
  n <- totNum*perc/100
  ntmp <- 0
  spy <- 0
  for (i in 1:nrow(dataOrd)) {
    if (ntmp < n) {
      ntmp <- ntmp + dataOrd[i,k]
      spy <- i
    } else {
      ntmp <- ntmp
      spy <- spy
    }
  }
  speciesSet[[k]] <- dataOrd[c(1:spy),1]
}

# NOTE: using the names of the taxonomic entries instead of numerical indexing means that the levels of the former are fewer
# than those of the latter, as the data frame has duplicate entries in the indentification column

speciesSet <- factor(levels(factor(unlist(lapply(speciesSet, as.character)))))

#now subset original data

subsetData <- allMeans[allMeans[,1] %in% speciesSet,]  # there we go


# comment: this sucks because the original identification is approximate and full of duplicates stemming from typos.
# this must be redone with a decent input file, that requires some work. Also all the descriptions must be pooled
# as undetermined. this removes the original "species or description" routine, as it will be done in excel. IN any case
# there don't seem to be too many duplicates, it's more about a huge quantity of approximate descriptions.

# way out is calling undetermined all that is not a species AT THIS STAGE, because if done at the beginning then Undetermined
# dominate the community and we miss the whole point of getting rid of not abundant species. easiest way is excel

write.csv(subsetData$speciesOrDescription, "subsetOfSpecies.csv")

string <- read.csv("//Staff/Home/SCIFAC/rovellal/DocumentsRedir/Data/Hoga/buoy3/speciesString.csv", header = F)[,1]

# get the levels

finalSpecies <- levels(factor(string))

# to do now is to substitute the first column with the appropriate entry

subsetData$speciesOrDescription <- as.character(subsetData$speciesOrDescription)

for (i in 1:nrow(subsetData)) {
  if (subsetData[i,1] %in% finalSpecies) {
    subsetData[i,1] <- subsetData[i,1]
  } else {
    subsetData[i,1] <- "Undetermined"
  }
}

subsetData$speciesOrDescription <- as.factor(subsetData$speciesOrDescription)

# ad hoc correction for pseudocertaina (which is even spelled wrong actually)

subsetData[subsetData$speciesOrDescription == "Pseudocertaina sp. ?",1] <- "Pseudocertaina sp."
subsetData$speciesOrDescription <- droplevels(subsetData$speciesOrDescription)

# now need to split and sum. is this legit with the means though? there isn't really any other solution anyway, as
# we need to choose between this "inaccurate" method and keeping a lot of species in the routine that are not abundant

# split-apply-combine routine

splitFrame <- split(subsetData, subsetData$speciesOrDescription)
splitFrameSums <- lapply(splitFrame, function(x) {
  y <- x[,-1]
  z <- colSums(y)
  return(c(as.character(x[1,1]), z)) # turns all into character, annoyingly enough
}
)
frameSums <- abind(splitFrameSums, along = 0)

frameSums <- cbind.data.frame(factor(frameSums[,1]), apply(frameSums[,-1], 2, as.numeric))

# after all the workarounds it is still full of typos, but especially the undetermined make up for a big chunk of the community.
# MDS is going to be very different and fundamentally wrong because of the heavy pooling. not recommended for MDS!!

#########################################################################################



# from here on the script must take different paths. you did it wrong for the quadrats fate though.

# getting messy. plot averages now. Transposition should be enough. In all of this we completely ignore the 
# uncertainty and this is not OKAY

transDatatmp <- as.data.frame(t(frameSums))
transDatatmp <- transDatatmp[-1,]
transData <- as.data.frame(apply(transDatatmp, 2, function(x) as.numeric(as.character(x))))

rownames(transData) <- rownames(transDatatmp)

transData$Total <- rowSums(transData) 

ID <- rownames(transData)
transData$Year <- as.numeric(substr(ID, 5, 6))
transData$Site <- substr(ID, nchar(ID), nchar(ID))

# need to melt here

meltData <- melt(transData, id.vars = c("Year", "Site"))

meanPlot <- ggplot(data = meltData[meltData$variable != "Total",], aes(x = Year, y = value, color = variable))+
  geom_line()+
  geom_point()+
  #scale_y_continuous(trans = "log10")+
  theme_bw()+
  facet_grid(.~Site)
meanPlot

ggsave("//Staff/Home/SCIFAC/rovellal/DocumentsRedir/Wellington/words/latex/hogaExploratory/Pics/meanPlot.pdf",
       meanPlot, useDingbats = T) 







###### ATTENTION!!! ######

# MDS routine is extremely biased by the pooling of undetermined species. I strongly advise against using the pooled taxonomy
# for this analysis, as the information is WRONG.
# mds ROUTINE MUST RESUME THE MEAN DATAFRAME FROM THE BEGINNING

# transform in log(x+1). we can debate wether we want to have this done here or earlier. as a rule of thumb,
# the earlier the better. Need to investigate the sensitivity to this of the MDS plot though, prepare 2 versions

allMeansTrans <- log10(allMeansTMP+1)


# now get a list out of this, each element of a list is a time series with the means for a site. One element per site
# i.e. 3 elements

sites <- levels(factor(unlist(lapply(names(allMeansTrans), function(x) substr(x, nchar(x), nchar(x))))))

sitesList <- vector(mode = "list", length = length(sites))


for(i in sites) {
  sitesList[[match(i, sites)]] <- allMeansTrans[,grepl(i, names(allMeansTrans))]
}

# now the similarity index must be calculated

matrixList <- list()

for (s in 1:length(sitesList)) {
  
  S <- sitesList[[s]]
  
  grandBC <- list()
  
  for (i in 1:length(S)) {
    BC <- vector(mode = "logical", length = length(S))
    for (j in 1:length(S)) {
      reduced <- S[,c(i,j)]
      BCpartial <- apply(reduced, 1, function(x) {
        numBC <- 2*min(x[1], x[2])
        denBC <- x[1]+x[2]
        metrics <- c(numBC, denBC)
      })
      BC[j] <- rowSums(BCpartial)[1]/rowSums(BCpartial)[2]
    }
    grandBC[[i]] <- BC
  } # minchia che culo
  
  kekes <- abind(grandBC, along = 0)
  colnames(kekes) <- names(S)
  rownames(kekes) <- names(kekes)
  matrixList[[s]] <- kekes
}


# now MDS plots how the hell do they work

# turn into bloody dissimilarity

dissList <- lapply(matrixList, function(x) 1-x)

fitList <- lapply(dissList, function(x) {
  fit <- isoMDS(x, k = 2)
  y <- fit$points[,1]
  z <- fit$points[,2]
  yzFrame <- data.frame(y,z)
  yzFrame$Names <- names(x[1,])
  return(yzFrame)
}
)


plot(fitList[[3]][,1], fitList[[3]][,2], xlab="Coordinate 1", ylab="Coordinate 2", 
     main="Nonmetric MDS", type="n")
text(fitList[[3]][,1], fitList[[3]][,2], labels = names(dissList[[1]][1,]), cex=.7)

# please turn this to ggplot for the sake of decency

allFrames <- as.data.frame(abind(fitList, along = 1))
allFrames$y <- as.numeric(as.character(allFrames$y))
allFrames$z <- as.numeric(as.character(allFrames$z))
allFrames$Names <- as.character(allFrames$Names) # for substr
allFrames$Site <- rep(sites, each = 10)
# add year for colour key

allFrames$Year <- rep(NA, nrow(allFrames))

for (i in 1:nrow(allFrames)) {
  allFrames$Year[i] <- paste("20", substr(allFrames$Names[i], 
                                          (nchar(allFrames$Names[i])-2), 
                                          (nchar(allFrames$Names[i])-1)),
                             sep = "")
}

# plot

MDSplotYears <- ggplot(data = allFrames, aes(x = y, y = z, group = Site, color = Year))+
  geom_text(aes(label = Year, color = Year))+
  theme_bw()+
  facet_grid(.~Site)
MDSplotYears



# Note: MDS plot here is different because it's calulated on the entire community, and not on the subsetted
# number of taxonomic unites forming the 90% of the abundance. This is because in this script the subsetting
# routine is carried out on the means, in order to reduce (supposedly) the number of taxonomic units to
# a smaller number. Should be flagged and transformed in binary routine depending on the desired analysis