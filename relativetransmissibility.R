## Code adapted from that used by Ben Swallow in Panovska-Griffiths et al., 2022

#Ensures that all necessary packages are installed

if(!require(tidyverse)){install.packages("tidyverse");require(tidyverse)}
if(!require(mgcv)){install.packages("mgcv");require(mgcv)}
if(!require(emmeans)){install.packages("emmeans");require(emmeans)}
if(!require(parallel)){install.packages("parallel");require(parallel)}
if(!require(ggplot2)){install.packages("ggplot2");require(ggplot2)}
if(!require(broom)){install.packages("broom");require(broom)}
if(!require(rgdal)){install.packages("rgdal");require(rgdal)}
if(!require(dplyr)){install.packages("dplyr");require(dplyr)}
# Calls relevant packages that are already in the library
library("magrittr")

# Reads cases data, downloaded from COG-UK Consortium: https://covid19.sanger.ac.uk/lineages/raw?latitude=52.001820&longitude=-2.155859&zoom=4.713940&show=A%2CBA.5
data <- read.delim("/Data/lineages_by_ltla_and_week_4_12.tsv")
# Reads longitudes and latitudes for Local Authority Districts, downloaded from: https://geoportal.statistics.gov.uk/datasets/ons::local-authority-districts-december-2020-uk-buc/about
bounds<-read.csv("/Data/Local_Authority_Districts_(December_2020)_UK_BUC.csv", stringsAsFactors=TRUE)
data$day<-as.numeric(as.Date(data$WeekEndDate)-as.Date(data$WeekEndDate[1])+1)
# Set lineages of interest
d2 <- data[data$Lineage=="Variant 1"|data$Lineage=="Variant 2",]
# Converts data from numeric to factor
d2=d2%>%mutate(LTLA = as.factor(LTLA))
d2=d2%>%mutate(Lineage = as.factor(Lineage))
# Adds latitude and longitude for each LTLA to d2
d2$lat <- bounds$LAT[match(d2$LTLA,bounds$LAD20CD)]
d2$lon <- bounds$LONG[match(d2$LTLA,bounds$LAD20CD)]
# Sorts grouping of d2
d2<-d2%>%group_by(day, Lineage)

## Model S
# A HGAM is fitted with the smoother estimated separately for each variant but with the same wiggliness in each case.


cov_modS <- bam(Count ~ t2(lat, lon, day, Lineage, bs=c("tp", "tp", "tp", "re"),
                           k=c(10, 10, 10, 6), m=2, full=TRUE),
                data=d2, method="fREML",discrete=T, family="nb",nthreads=3)


## Model I
# A HGAM is fitted with the smoother estimated separately for each variant and with different wiggliness.  

cov_modI <- bam(Count ~ Lineage + te(lat, lon, day, by=Lineage,
                                     bs=c("tp", "tp", "tp"), k=c(10, 10, 10), m=2),
                data=d2, method="fREML", family="nb",nthreads=3)


# checks the basis sizes are well chosen
gam.check(cov_modS)
gam.check(cov_modI)
# 
k.check(cov_modS)
k.check(cov_modI)

summary(cov_modS)
summary(cov_modI)

# Creates a cluster of 3 nodes for parallel computation
cl <- makeCluster(3)
## Model G
# Fits a GAM with a global smoother fixed across all variants
cov_modG <- bam(Count ~ t2(lat, lon, day, bs=c("tp", "tp", "tp"), k=c(10, 10, 10)),
                data=d2, method="fREML",discrete=T, family="nb",nthreads=3)


## Model GS
# Adds a HGAM similar to model S but with an additional global smoother
cov_modGS <- bam(Count ~ t2(lat, lon, day, bs=c("tp", "tp", "tp"),
                            k=c(10, 10), m=2) +
                   t2(lat, lon, day, Lineage, bs=c("tp", "tp", "tp", "re"),
                      k=c(10, 10, 10, 6), m=2, full=TRUE),
                 data=d2, method="fREML", family="nb")


# Creates a cluster of 3 nodes for parallel computation
cl <- parallel::makeCluster(3)
## Model GI
# Fits a HGAM similar to model I but with an additional global smoother
cov_modGI <- bam(Count ~ Lineage+
                   t2(lat, lon, day, bs=c("tp", "tp", "tp"), k=c(10, 10, 10), m=2) +
                   te(lat, lon, day, by=Lineage, bs= c("tp", "tp", "tp"),
                      k=c(10, 10, 10), m=1),
                 data=d2, method="fREML", family="nb",cluster=cl)


# checks the basis sizes are well chosen
gam.check(cov_modG)
gam.check(cov_modGS)
gam.check(cov_modGI)
# 
k.check(cov_modG)
k.check(cov_modGS)
k.check(cov_modGI)

# Check fit uses Akaike's Information Criterion
AIC<-AIC(cov_modS,cov_modI,cov_modG,cov_modGS,cov_modGI)

## Analysis of Best-fit model

# (this code is adapted from Davies et al., 2021, Science)
# Function to change difference in exponential growth rate to a relative transmissibility value. See page 15.
M.from.delta_r = function (delta_r, g=5.5) {
  delta_R = exp(delta_r*g)
  return( delta_R )
}

Rt.from.r = function(r, g=4.7, sigma=2.9) {
  k <- (sigma / g)^2
  Rt <- (1 + k * r * g)^(1 / k)
  return(Rt) }

M.from.delta_r_df = function (df, g1=5.5, g2=3.6,
                              coln=c("M1","M1.LCL","M1.UCL",
                                     "M2","M2.LCL","M2.UCL")) {
  df_num = df[,which(unlist(lapply(df, is.numeric))), drop=F]
  df_nonnum = df[,which(!unlist(lapply(df, is.numeric))), drop=F]
  df_out1 = apply(df_num, 2, function (delta_r) M.from.delta_r(delta_r, g1))
  if (class(df_out1)[1]=="numeric") df_out1=as.data.frame(t(df_out1), check.names=F)
  df_out2 = apply(df_num, 2, function (delta_r) M.from.delta_r(delta_r, g2))
  if (class(df_out2)[1]=="numeric") df_out2=as.data.frame(t(df_out2), check.names=F)
  df_out = data.frame(df_out1, df_out2, check.names=F)
  if (!is.null(coln)) colnames(df_out) = coln
  return( data.frame(df_nonnum, df_num, df_out, check.names=F) )
}

# Specify specific latitudes/longitudes of all LTLAs
ull <- na.omit(unique(d2[,c("lat","lon")]))

output <- NULL

# For each LTLA calculate relative multiplicative growth rate factor
for(i in 1:nrow(ull)){
  
  # Get trends during time period of interest of the exponential growth rate, enter time period of interest by replacing X's in format yyyy-mm-dd
  mfit_emtrends <- emtrends(cov_modGI, ~ Lineage*lat*lon, "day", mode="latent",lat=list(lat=unlist(ull[i,1]),lon=unlist(ull[i,2]),day=seq(as.numeric(as.Date("202X-XX-XX",format="%Y-%m-%d")-as.Date(d2$WeekEndDate[1])),as.numeric(as.Date("202X-XX-XX",format="%Y-%m-%d")-as.Date(d2$WeekEndDate[1])), by=2)))
  # Calculate difference between trends in each LTLA
  mfit_contrasts = data.frame(as.data.frame(contrast(mfit_emtrends, method="trt.vs.ctrl", ref=1, reverse=TRUE, adjust="sidak")),
                              as.data.frame(confint(contrast(mfit_emtrends, method="trt.vs.ctrl", ref=1, reverse=TRUE, adjust="sidak")))[,c("lower.CL","upper.CL")])
  colnames(mfit_contrasts)[which(colnames(mfit_contrasts) %in% c("estimate","lower.CL","upper.CL"))] =
    c("delta_r","delta_r.lower.CL","delta_r.upper.CL")
  
  # Convert the differences to multiplicative advantages
  mfit_contrasts = data.frame(mfit_contrasts[,c("contrast","SE","t.ratio","p.value")],
                              M.from.delta_r_df(mfit_contrasts[,c("delta_r","delta_r.lower.CL","delta_r.upper.CL")]))
  
  
  output <- rbind(output,mfit_contrasts)
  
}

# Output estimate of multiplicative advantage  assuming 5.5 day generation time
out2comp <- matrix(output$M1,nrow=1,byrow=T)


latlon <- expand.grid(unlist(ull[,1]),unlist(ull[,2]))
latlon <- latlon[order(latlon[,2]),]

llind <- sapply(1:311,function(i)which(latlon[,1]==unlist(ull[i,1])&latlon[,2]==unlist(ull[i,2])))


## Plots

# Read in England shape file; downloaded from https://borders.ukdataservice.ac.uk/easy_download_data.html?data=England_ct_2011
my_spdf <- readOGR(
  dsn="/Data/england", layer="england_ct_2011",
  verbose=FALSE
)

summary(my_spdf)

# Transform projection into latitude/longitude values to match the LTLA data
my_spdf <- spTransform(my_spdf, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

spdf_fortified <- tidy(my_spdf)


# Plot shparefile to get the outline of the map of England
names(latlon) <- c("lat","lon")
plt = ggplot() +
  geom_polygon(data = spdf_fortified, aes( x = long, y = lat, group = group), fill="#69b3a2", color="white")


#Add points for LTLAs where relative transmissibility of Variant 2 vs Variant 1 is greater than 1 i.e. there is an increase in transmissibility
plt3 = plt +
  geom_point(data = latlon[llind,][which(out2comp[1,]>1),], aes(x = lon, y = lat, size = out2comp[1,which(out2comp[1,]>1)],col= out2comp[2,which(out2comp[1,]>1)]),
             shape = 20) +
  scale_fill_gradient(name="Multiplicative Rt",low = "yellow", high = "red", na.value = NA,aesthetics = "colour")  +
  labs(title="Variant 2 vs Variant 1")+ guides(size="none") +
  theme_void()

#Save plot
ggsave('Variant2vVariant1.pdf', plt3) 

#Save information on the best fit model
save(cov_modGI, file="Variant2vVariant1GI.rdata")
#Save AIC values for the models
save(AIC, file="Variant2vVariant1AIC.rdata")
#Save relative transmissibilities in each LTLA
save(out2comp, file="Variant2vVariant1out2comp.rdata")
