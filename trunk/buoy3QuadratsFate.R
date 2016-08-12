# 12.08.2016 Alberto.Rovellini@vuw.ac.nz 
# script to follow the fate of sponges within a quadrat

require(abind)
require(plyr)
require(ggplot2)
require(reshape2)

dataAllYears <- read.csv("/home/somros/Documents/R/exploratoryHoga/input/spongeAbundanceQuadrats.csv")

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

# gets rid of zeros as characters in the datasets

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


# literally about transposing the dataframe here, nothing more. maybe melt?

transData <- as.data.frame(t(dataAllYears))
transData <- transData[-1,] # gets rid of the line with the species IDs
transData$Total <- rowSums(transData)

# need to separate years (my x variable) and quadrats

quadratID <- rownames(transData)
transData$Year <- as.numeric(substr(quadratID, 1, 2))
transData$Site <- substr(quadratID, nchar(quadratID)-1, nchar(quadratID)-1)
transData$Quadrat <- substr(quadratID, nchar(quadratID), nchar(quadratID))  
transData$QuadratSite <- substr(quadratID, nchar(quadratID)-1, nchar(quadratID))

# restrict to small data frame with only the total numbers and the quadrats IDs

smallFrame <- transData[,125:129]

library(RColorBrewer)
par(mar = c(0, 4, 0, 0))
display.brewer.all()
nOfColorsQuad <- 5 # numer of quadrats
getPaletteQuad <- colorRampPalette(brewer.pal(5, "Set1"))

quadratsFate <- ggplot(data = smallFrame, aes(x = Year, y = Total, group = QuadratSite, color = Quadrat))+
  geom_line(aes(color = Quadrat), size = .6)+
  geom_point()+
  scale_color_manual(values = getPaletteQuad(nOfColorsQuad))+
  scale_x_continuous(breaks = seq(5,16,1),
                     labels = seq(5,16,1),
                     limits = c(5,16))+
  scale_y_continuous(limits = c(0,350),
                     breaks = seq(0,350,50))+
  labs(y=expression(paste(Sponge~density~(n~m^2))))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank())+
  theme(plot.title = element_text(size=14, vjust=2))+
  theme(axis.title.x = element_text(size=10,vjust=-0.5),
        axis.title.y = element_text(size=10,vjust=0.5))+
  theme(axis.text.x=element_text(size=10, angle = 60, 
                                 hjust = 1, vjust = .9))+
  theme(axis.text.y=element_text(size=10))+
  facet_grid(. ~ Site, scales = "free")
quadratsFate

ggsave("/home/somros/Documents/R/exploratoryHoga/output/pics/quadratsFate.pdf", quadratsFate,
       width=7, height=4, useDingbats=T)

# next thing to explore is the variability within and between sites



