##############################################################
##  Clean CPC 3007 data and create summary statistics
##
##  David Morley, d.morley@imperial.ac.uk
##  Imperial College London
##  Jul-2016
##############################################################

require(lubridate)
require(ggplot2)

# 1)  Data must be downloaded from CPC using Aerosol Instrumnet Manager Software
# 2)  Input must be in CSV format with 1st Column as time in seconds, 
#     subsequent columns are individual measurements (usually 30min) and have a 
#     header for the site name (e.g. X1, X2...) see object 'd1' in the code below
# 3)  Input data must also be supplied with a key linking by site name (X1, X2...)
#     see object 'h1' in the code below
# 4)  Output creates a cleaned series and summary statistics


######################################
## Cleaning functions
######################################

filterCPCv1 <- function(vectorTOfilter) {
  tempdat <- vectorTOfilter[-1] 
  validSize <- tempdat > 100
  Factor <- ratioV3(vectorTOfilter, 10)
  validFactor <- !Factor[,2]
  keepIndex <- validSize & validFactor
  tempdat[which(keepIndex == FALSE)] <- NA
  return(tempdat)
}

# This function written by Ming Tsai of Swiss TPH
ratioV3=function(v,x){ #(vector v, factor x)
  vectorLen<-length(v)
  fact<-x
  tempdat0<-v[-vectorLen] 
  tempdat1<-v[-1]       
  ratio<-tempdat1/tempdat0
  ratio1<-as.data.frame(cbind(ratio,ratio>fact|ratio<1/fact)) #the length of this vector is 1 less than the input
  cumsum1<-cumsum(ratio1[,2]) #cumulative sum.
  rle1<-rle(cumsum1) #run length encoding--Compute the lengths and values of runs of equal values
  tossDF<-as.data.frame(cbind(rle1$values,rle1$lengths)) #dataframe with values and their lengths
  names(tossDF)<-c('values','lengths')
  tossDF<-tossDF[-1,] #remove value of '0' from being tossed.
  index<-which(tossDF$lengths==1) #if a value occurs only once, that means that it is surrounded on both sides by values that are a factor of >10 different.
  valsToRm<-tossDF$values[index] #using the index, pick out the cumulative sums to be removed.
  to_remove<-is.element(cumsum1,valsToRm)                   
  ratio1<-cbind(ratio1,cumsum1,to_remove)
  names(ratio1)[2]<-paste('GTorLTfact_',fact,sep='')
  return(ratio1)
}       

######################################
## Read and clean raw CPC data
######################################

# Read in CPC UFP 1-second data
d1 <- read.csv("I:\\..\\CPCUFP_norwich_1.csv")

# These are data outputted from Aerosol Instrument Manager
# 1st Column is time from start in seconds
# Remaining X columns are site numbers (161 sites in Norwich)
head(d1)
# TIME X1 X2 X3 X4 X5 X6 X7 X8 X9 X10 X11 X12 X13 X14 X15 X16 X17 X18 X19 X20 X21 X22 X23 X24 
# 1    1  0  0  0  0  0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   
# 2    2  0  0  0  0  0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   
# 3    3  0  0  0  0  0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   
# 4    4  0  0  0  0  0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   
# 5    5  0  0  1  0  0  0  0  0  0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   
# 6    6  0  0  0  0  0  0  0  0  0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   

# Read UFP time and date info
h1 <- read.csv("I:\\..\\CPCUFP_norwich_1_key.csv")

head(h1)
# SITE       DATE VISIT STARTTIME ENDTIME                                   NOTES
# 1    1 13/08/2014     1     13:24   13:54                                        
# 2    2 13/08/2014     1     12:42   13:12   battery ran out on noise, so shorter
# 3    3 21/08/2014     1     11:48   12:18                                        
# 4    4 21/08/2014     1     13:18   13:48   Actually outside 178 Dereham, not 153
# 5    5 21/08/2014     1     13:59   14:29   Bonfire nearby? Especially first 10mins
# 6    6 13/08/2014     1     12:04   12:34   

# Clean CPC sites
# Removed values are left as NA
d1.c <- data.frame(apply(d1[2:ncol(d1)], 2, filterCPCv1))


# Plot a particular site's series
site = "X2"
ggplot(data = data.frame(seconds = seq_along(d1.c[, site]), UFP = d1.c[, site]), 
       aes(x = seconds, y = UFP)) + geom_point()


# Create summary statistics
cpc <- data.frame(Site = h1$SITE)
cpc$mean <- sapply(d1.c, mean, na.rm = TRUE)
cpc$median <- sapply(d1.c, median, na.rm = TRUE)
cpc$stdev <- sapply(d1.c, sd, na.rm = TRUE)
quant <- apply(d1.c, 2, quantile, probs = c(0.05, 0.25, 0.75, 0.95),  na.rm = TRUE)
cpc$P05 <- quant[1, ]
cpc$P25 <- quant[2, ]
cpc$P75 <- quant[3, ]
cpc$P95 <- quant[4, ]
cpc$date <- h1$DATE
cpc$time <- h1$STARTTIME
cpc$posix <- strptime(paste(cpc[,"date"], cpc[,"time"], sep = ":"), format = "%d/%m/%Y:%H:%M")
head(cpc)

# Site     mean median     stdev     P05     P25      P75     P95       date  time               posix
# 1    1 4995.725 4699.0 1834.7784 3285.40 3986.00  5523.00  7106.8 13/08/2014 13:24 2014-08-13 13:24:00
# 2    2 8361.276 8030.0 1751.6371 6158.70 6719.00  9946.50 11112.2 13/08/2014 12:42 2014-08-13 12:42:00
# 3    3 9963.710 8277.5 5703.0031 4829.50 6479.00 11894.25 19580.5 21/08/2014 11:48 2014-08-21 11:48:00
# 4    4 7943.057 6757.0 4541.0931 5146.30 5846.50  8022.50 14499.6 21/08/2014 13:18 2014-08-21 13:18:00
# 5    5 8650.912 7356.5 3092.0647 6029.45 6740.25  9559.50 15444.3 21/08/2014 13:59 2014-08-21 13:59:00
# 6    6 6094.898 6105.0  931.1926 5018.55 5497.00  6510.25  6994.9 13/08/2014 12:04 2014-08-13 12:04:00


## Save Summary and cleaned series
write.csv(cpc, "I:\\..cpc_summary.csv")
write.csv(d1.c, "I:\\..cpc_cleaned.csv")





