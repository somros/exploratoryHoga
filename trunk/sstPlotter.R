# 12.08.2016 Alberto.Rovellini@vuw.ac.nz
# script to plot the SST data from Coral Reef Watch (NOAA), from 2001 to 2016.
# selected station is Wakatobi (see http://coralreefwatch.noaa.gov/satellite/vs/docs/list_vs_group_latlon_201103.php)
# for information about the dataset, see "description.txt", or also http://coralreefwatch.noaa.gov/satellite/data_nrt/descriptions/timeseries.txt

setwd("/home/somros/Documents/Data/Hoga/Temperature")

dataSet <- read.table("wakatobi.txt", header = T, sep = "", dec = ".")

# high resolution, sampled every 3 days. can take monthly means and look for trends during time maybe
# (don't expect to see any actually)
# might actually want to use a oneliner instead of the split-apply-combine
# data from Hoga come from the same season each year, but it makes no sense to compare with
# mean values for that time because there is certainly a delay, assuming that there is an
# effect at all. So how do we detect it if it's there? Where to look?



