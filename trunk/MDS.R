# 17/10/2016 script for MDS analysis on the Buoy 3 data. ONe of tw twin scripts for the moment, one
# for across year comparison of the same site and one for across sites comparison of the same year. 


require(abind)
require(MASS)
require(ggplot2)
require(BiodiversityR)


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
  fit <- isoMDS(x, k = 2) # in case of metric MDS just swap isoMDS() with cmdscale() and remove $points below (flag this)
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
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_text(aes(label = Year, color = Year), size = 4.5, fontface = "bold")+
  geom_path(size = 1)+
  scale_x_continuous(limits = c(-.25,.25))+
  scale_y_continuous(limits = c(-.25,.25))+
  theme_bw()+
  facet_grid(.~Site)
MDSplotYears

# ggsave("//Staff/Home/SCIFAC/rovellal/DocumentsRedir/Wellington/words/latex/hogaExploratory/Pics/MDSYears.pdf",
#        MDSplotYears, useDingbats = T) 




# CAP uses a dedicated package and starts from the species list
#


