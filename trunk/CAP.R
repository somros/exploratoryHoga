# 27.10.2016 script for  CAP analysis. Ordination of sites depending on their composition, as decided by the abundance of each species,
# correlated with environmental variables (namely time and space in this case). Environmental varibales are ANOVA-like factors and
# not continuous dimensions. 

# 8.11.2016 contains the MDS analysis as well. use this for any communty analysis and keep the sandbox for the rest

# TODO: play with T data and other plausible environmental variables if you want to go on a nice journal


require(abind)
require(MASS)
require(ggplot2)
require(BiodiversityR)
require(ggvegan)


setwd("//Staff/Home/SCIFAC/rovellal/DocumentsRedir/Data/Hoga/buoy3")
myData <- read.csv("spongeAbundanceQuadrats.csv")

# turn to numeric immediately

myData <- myData[c(1:124),]
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


# first we transpose (which is not necessary now but should become at some point I reckon)

transMatrix <- t(allMeansTrans)

# now structure is N as rows and species as columns. However, we osu old structure for computation (same ending point and I don't
# have to recode it all up)

grandBC <- list()

for (i in 1:length(allMeansTrans)) { # home-made Bray-Curtis, makes little sense as there's vegdist in vegan doing it but it's solid (validated)
  BC <- vector(mode = "logical", length = length(allMeansTrans))
  for (j in 1:length(allMeansTrans)) {
    reduced <- allMeansTrans[,c(i,j)]
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
colnames(kekes) <- names(allMeansTrans)
rownames(kekes) <- names(kekes)


# now metric MDS plots, or PCO

# turn intodissimilarity

dissMatrix <- 1-kekes

fit <- isoMDS(dissMatrix, k = 2) # this is a ***NON-METRIC*** MDS plot based on Brey-Curtis dissimilarity.

distances <- as.data.frame(fit$points)
rownames(distances) <- names(allMeansTrans)
distances$Names <- as.character(rownames(distances))
for (i in 1:nrow(distances)) {
  distances$Year[i] <- paste("20", substr(distances$Names[i], 
                                          (nchar(distances$Names[i])-2), 
                                          (nchar(distances$Names[i])-1)),
                             sep = "")
}
for (i in 1:nrow(distances)) {
  distances$Site[i] <- substr(distances$Names[i], nchar(distances$Names[i]), nchar(distances$Names[i]))
}

# plot

MDSplot <- ggplot(data = distances, aes(x = V1, y = V2, group = Site, color = Year))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_vline(xintercept = 0, linetype = "dashed")+
  geom_text(aes(label = Year, color = Year), size = 4.5, 
            fontface = "bold", vjust = 0, nudge_y = 0.02)+
  guides(color=FALSE)+
  geom_point(aes(shape = Site, color = Year), size = 2)+
  scale_shape_manual(values = c(0,1,2))+
  labs(x = "MDS1", y = "MDS2")+
  #geom_path(size = 1)+
  scale_x_continuous(limits = c(-.7,.7))+
  scale_y_continuous(limits = c(-.7,.7))+
  theme_bw()+
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank())#+
  #facet_grid(.~Site)
MDSplot

ggsave("//Staff/Home/SCIFAC/rovellal/DocumentsRedir/Wellington/words/latex/hogaExploratory/Pics/MDSplot.pdf",
               MDSplot, useDingbats = T, width = 6, height = 5)

# calculating the dissimilarity matrix among all instances (site*year) we see an influence of space in clustering.
# temporal variation and return to original conditions is a bit less evident but still well visible considering the 
# sites one at a time. good plot, now need to go on and end up with the CAP for the species!!!


# CAP uses a dedicated package and starts from the species list

# need to build a dummy frame apparently with data

dummyMatrix <- as.data.frame(cbind(distances$Year, distances$Site))
colnames(dummyMatrix) <- c("Year", "Site")
rownames(dummyMatrix) <- rownames(transMatrix)

# different formulas mean different ordinations, depending on the included variables

CAPall <- capscale(transMatrix ~ Year + Site, dummyMatrix, distance = "bray", add = T)
CAPTime <- capscale(transMatrix ~ Year, dummyMatrix, distance = "bray", add = T)
CAPSpace <- capscale(transMatrix ~ Site, dummyMatrix, distance = "bray", add = T)

# plots for different formulas

plot(CAPall)
plot(CAPTime)
plot(CAPSpace)

# all the below is just to get fancy with ggplot. Packages vegan and ggvegan will do all the job

components <- fortify(CAPall) # extract the coordinates from the cca objects calculated with capscale 
components <- data.frame(lapply(components, function(x) {
  x <- gsub("Year", "", x)
  x <- gsub("Site", "", x)
  x <- gsub("Mean", "", x)
  return(x)
}))
components$Dim1 <- as.numeric(as.character(components$Dim1))
components$Dim2 <- as.numeric(as.character(components$Dim2))
components$Label <- gsub("spe", "", components$Label)


# plotting region

plotData <- subset(components, components$Score != "biplot" & components$Score != "constraints" &
                    components$Score != "sites")
plotData$ColorKey <- c(rep("species", 75), rep("year", 10), rep("site", 3)) # make this flexible on the iterations
# fix levels
plotData$ColorKey <- factor(plotData$ColorKey, levels = unique(plotData$ColorKey))


CAPplot <- ggplot(data = plotData, aes(x = Dim1, y = Dim2, group = ColorKey))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_vline(xintercept = 0, linetype = "dashed")+
  geom_segment(data = subset(plotData, plotData$Score == "centroids"),
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2), arrow = arrow(length = unit(1/2, 'picas')))+
  geom_text(aes(label = Label, color = ColorKey, size = ColorKey))+
  scale_color_manual(values = c("darkgrey", "red3", "cyan4"))+
  scale_size_discrete(range = c(4,8))+
  theme_bw()+
  labs(x = "CAP1", y = "CAP2")+
  scale_x_continuous(limits = c(-.8,.8))+
  scale_y_continuous(limits = c(-.8,.8))+
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank())
CAPplot

ggsave("//Staff/Home/SCIFAC/rovellal/DocumentsRedir/Wellington/words/latex/hogaExploratory/Pics/CAPplot.pdf",
       CAPplot, useDingbats = T, width = 6.5, height = 5)
