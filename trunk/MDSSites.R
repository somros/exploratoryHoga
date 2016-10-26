# MDS plot per site instead of per year

require(abind)
require(MASS)
require(ggplot2)

setwd("//Staff/Home/SCIFAC/rovellal/DocumentsRedir/Data/Hoga/buoy3")
myData <- read.csv("spongeAbundanceQuadrats.csv")

# turn to numeric immediately

myData <- myData[c(1:127),]
myData <- apply(myData, 2, as.numeric)
myData[is.na(myData)] <- 0 # turns NAs in 0, as necessary for Bray-Curtis


# some steps to break down the frame in a list of frames, one per year per site.

nOfQuadrats <- 5



# TODO: introduce routine to cut data frames to the species accounting for 90% of the abundance. must be applied
# to each quadrat and then all the species considered make the subset key for the entire data frame

# we do it with apply. jk lol use a loop apply is slow af anyway

perc <- 90

species <- vector("list", length = ncol(myData)-1)

for (k in 2:ncol(myData)) {
  dataOrd <- myData[order(myData[,k], decreasing = T),]
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
  species[[k]] <- dataOrd[c(1:spy),1]
}

species <- factor(levels(factor(unlist(species))))

# now subset original data

subsetData <- myData[myData[,1] %in% species,] # there we go


# if for whatever reason we want to use the entire set of species for the analysis, we can get rid of the above
# section (this should become a flag). 
# the above routine may be better off if used on the means. That would rule out all the outliers from the
# single quadrats. The number of species is however half the original amount.


dataQuad <- subsetData[,-1] # remove species for convenience


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
  y <- x[-1,]
  z <- rowMeans(y)
  frame <- as.data.frame(z[1:nrow(dataQuad)])
  colnames(frame) <- paste("Mean", substr(names(x[1,])[1], 2,3), substr(names(x[1,])[1], 5, 5), sep = "")
  return(frame)
})


allMeans <- as.data.frame(abind(newListMeans, along = 2))

# transform in log(x+1). we can debate wether we want to have this done here or earlier. as a rule of thumb,
# the earlier the better. Need to investigate the sensitivity to this of the MDS plot though, prepare 2 versions


allMeansTrans <- log10(allMeans+1)

# from now on as per in other script but using years instead


years <- levels(factor(unlist(lapply(names(allMeansTrans), function(x) {
  substr(x, nchar(x)-2, nchar(x)-1)}))))

yearsList <- vector(mode = "list", length = length(years))


for(i in years) {
  yearsList[[match(i, years)]] <- allMeansTrans[,grepl(i, names(allMeansTrans))]
}

# now the similarity index must be calculated

matrixList <- list()

for (y in 1:length(yearsList)) {
  
  Y <- yearsList[[y]]
  
  grandBC <- list()
  
  for (i in 1:length(Y)) {
    BC <- vector(mode = "logical", length = length(Y))
    for (j in 1:length(Y)) {
      reduced <- Y[,c(i,j)]
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
  colnames(kekes) <- names(Y)
  rownames(kekes) <- names(kekes)
  matrixList[[y]] <- kekes
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


# basically the other way round than the other scrpt. might make them methods to switch wiht flags

allFrames <- as.data.frame(abind(fitList, along = 1))
allFrames$y <- as.numeric(as.character(allFrames$y))
allFrames$z <- as.numeric(as.character(allFrames$z))
allFrames$Names <- as.character(allFrames$Names) # for substr
allFrames$Year <- rep(years, each = 3)
# add year for colour key

allFrames$Site <- rep(NA, nrow(allFrames))

for (i in 1:nrow(allFrames)) {
  allFrames$Site[i] <- substr(allFrames$Names[i], 
                                          (nchar(allFrames$Names[i])), 
                                          (nchar(allFrames$Names[i])))
}

# plot

MDSplotSites <- ggplot(data = allFrames, aes(x = y, y = z, group = Year, color = Site))+
  geom_text(aes(label = Site, color = Site))+
  theme_bw()+
  facet_wrap(~Year, nrow = 2)
MDSplotSites


# save this heresy in the latex folder

ggsave("//Staff/Home/SCIFAC/rovellal/DocumentsRedir/Wellington/words/latex/hogaExploratory/Pics/MDSSites.pdf",
       MDSplotSites, useDingbats = T) 
