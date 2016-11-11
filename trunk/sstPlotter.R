# 12.08.2016 Alberto.Rovellini@vuw.ac.nz
# script to plot the SST data from Coral Reef Watch (NOAA), from 2001 to 2016.
# selected station is Wakatobi (see http://coralreefwatch.noaa.gov/satellite/vs/docs/list_vs_group_latlon_201103.php)
# for information about the dataset, see "description.txt", or also http://coralreefwatch.noaa.gov/satellite/data_nrt/descriptions/timeseries.txt

require(ggplot2)
require(abind)
require(plyr)
require(lubridate)

dataSet <- read.table("//Staff/Home/SCIFAC/rovellal/DocumentsRedir/Data/Hoga/buoy3/wakatobi_SST.txt", 
                      header = T, sep = "", dec = ".")
days <- seq(as.Date("2000/11/28"), as.Date("2016/08/18"), "days")
julianDays <- yday(days)


# create a new column in the data set with the correct format of the days to match with progs
myDays <- list()
for (i in 1:nrow(dataSet)) {
myDays[[i]] <- paste(dataSet$BYYY[i],
                if (nchar(dataSet$BM[i]) == 2) {
                  dataSet$BM[i]
                } else {
                  paste("0", dataSet$BM[i], sep = "")
                }, 
                if (nchar(dataSet$BD[i]) == 2) {
                  dataSet$BD[i]
                } else {
                  paste("0", dataSet$BD[i], sep = "")
                } , sep = "-")
}

myDays <- unlist(myDays)

# now need to associate the progressive number of the days to my dataSet
# could use match from column of data frame to vector of dates to see the position apparently

indexDay <- match(myDays, as.character(days))

# now extract from the julianDays vector the elements indexed according to indexDay

julians <- julianDays[indexDay]
trueDates <- as.Date(days[indexDay])

dataSet$Julian <- julians
dataSet$Date <- trueDates


#####################################################
# non-averaged plot section

myYears <- 2005:2016
months <- c("J","F","M","A","M","J","J","A","S","O","N","D")
yearsOfInterest <- dataSet[dataSet$BYYY %in% myYears,]

allTemp <- ggplot(data = yearsOfInterest, 
                  aes(x=Julian, y=SST, group = BYYY))+
  geom_line()+
  scale_x_continuous(breaks = seq(15, 365, 30),
                     labels = 1:12)+
  labs(y = "SST (°C)", x = "Month")+
  theme_bw()+
  theme(panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=14, vjust=2))+
  theme(axis.title.x = element_text(size=10,vjust=-0.5),
        axis.title.y = element_text(size=10,vjust=0.5))+
  theme(axis.text.x=element_text(size=10))+
  theme(axis.text.y=element_text(size=10))+
  facet_wrap(~ BYYY)
allTemp

# plot anomalies
# anomalies are calculated as recorded SST - daily climatology. for incormation about the calculation
# of the latter, refer to http://coralreefwatch.noaa.gov/satellite/methodology/methodology.php#clim

anomalies <- ggplot(data = yearsOfInterest,
                    aes(x = Julian, y = SSTANOM))+
  geom_line()+
  geom_hline(yintercept=0, col="red")+
  scale_x_continuous(breaks = seq(15, 365, 30),
                     labels = months)+ # approximation to visualise
  facet_wrap(~ BYYY)
anomalies
  
# plot proper time series


ts <- ggplot(data = yearsOfInterest,
             aes(x = Date, y = SST))+
  geom_line()+
  scale_y_continuous(limits = c(15, 35))
ts


# high resolution, sampled every 3 days. can take monthly means and look for trends during time maybe
# (don't expect to see any actually)
# might actually want to use a oneliner instead of the split-apply-combine
# data from Hoga come from the same season each year, but it makes no sense to compare with
# mean values for that time because there is certainly a delay, assuming that there is an
# effect at all. So how do we detect it if it's there? Where to look?

# guess I have to work with monthly averages at the very least

monthListTmp <- split(dataSet, list(factor(dataSet$BM), factor(dataSet$BYYY)))
monthListTmp1 <- monthListTmp[-(1:10)]
monthList <- monthListTmp1[(length(monthListTmp1)-4:length(monthListTmp1))] # what the fuck
monthlyMeans <- lapply(monthList, function(x) {
  meanSST <- mean(x$SST)
  sdSST <- sd(x$SST)
  frameOfMeans <- data.frame(levels(factor(x$BYYY)), levels(factor(x$BM)),
                             meanSST, sdSST)
  return(frameOfMeans)
})

timeSeries <- as.data.frame(abind(monthlyMeans, along = 0))
colnames(timeSeries) <- c("Year", "Month", "Mean", "SD")
timeSeries <- as.data.frame(apply(timeSeries, 2, function(x) as.numeric(as.character(x))))

# set years of interest

myYears <- 2005:2016

meanSSTYears <- ggplot(data = timeSeries[timeSeries$Year %in% myYears,], 
               aes(x = Month, y = Mean))+
  geom_line(colour = "red")+
  scale_x_continuous(limits = c(1,12),
                     breaks = 1:12)+
  labs(y = "Mean monthly temperature (° C)")+
  theme_bw()+
  theme(panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=14, vjust=2))+
  theme(axis.title.x = element_text(size=10,vjust=-0.5),
        axis.title.y = element_text(size=10,vjust=0.5))+
  theme(axis.text.x=element_text(size=10, hjust = 1, vjust = .9))+
  theme(axis.text.y=element_text(size=10))+
  facet_wrap( ~ Year, ncol = 2)
meanSSTYears

# ggsave("/home/somros/Documents/R/exploratoryHoga/output/pics/meanSSTYears.pdf",
#        meanSSTYears, width=6, height=9, useDingbats=T)


# plot the monthly means over time

monthlyMeanInTime <- ggplot(data = timeSeries[timeSeries$Year %in% myYears,],
                            aes(x = Year, y = Mean))+
  geom_line()+
  facet_wrap(~ Month)
monthlyMeanInTime


# get climatology for Dec-Jan-Feb-Mar. 

myYears <- 2004:2016
myMonths <- c(12,1:3)

dataSST <- dataSet[dataSet$BYYY %in% myYears,]
dataSST <- dataSST[dataSST$BM %in%  myMonths,]

# do it in excel, I'll come back to this later if the results encourage it

# write.csv(dataSST, "//Staff/Home/SCIFAC/rovellal/DocumentsRedir/Data/Hoga/buoy3/SST.csv")
dataSST <- read.csv("//Staff/Home/SCIFAC/rovellal/DocumentsRedir/Data/Hoga/buoy3/SST.csv")

dataSplit <- split(dataSST, dataSST$Clim)
dataLapply <- lapply(dataSplit, function(x) mean(x$SST))
dataCombined <- cbind(myYears[-1], unlist(dataLapply))
dataCombined <- dataCombined[-c(6,8),]
write.csv(dataCombined, "//Staff/Home/SCIFAC/rovellal/DocumentsRedir/Data/Hoga/buoy3/cheat.csv")



# SOI

dataSOI <- read.table("//Staff/Home/SCIFAC/rovellal/DocumentsRedir/Data/Hoga/buoy3/SOI.txt", header = T, sep = ",")
dataSOI$Year <- substr(dataSOI$Date, 1, nchar(dataSOI$Date)-2)
dataSOI$Month <- substr(dataSOI$Date, nchar(dataSOI$Date)-1, nchar(dataSOI$Date))
mySOI <- dataSOI[dataSOI$Year %in% myYears,]
mySOI$Month <- as.numeric(as.character(mySOI$Month))
mySOI <- mySOI[mySOI$Month %in% myMonths,]
mySOI <- mySOI[-(1:3),]
mySOI$Index <- rep(1:(length(levels(factor(mySOI$Year)))-1), each = 4)
splitSOI <- split(mySOI, mySOI$Index)
lapplySOI <- lapply(splitSOI, function(x) mean(x$Value))
meanSOI <- cbind(rep(myYears[-1], each = 3), rep(unlist(lapplySOI), each = 3))
meanSOI <- meanSOI[meanSOI[,1] != 2010 & meanSOI[,1] != 2012,]
write.csv(meanSOI, "//Staff/Home/SCIFAC/rovellal/DocumentsRedir/Data/Hoga/buoy3/cheatSOI.csv")






