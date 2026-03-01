

##############

## Script file for running HGAMS and calculating relative transmissibility

##############



####Authored by B. Swallow

packageurl <- "http://cran.r-project.org/src/contrib/Archive/emmeans/emmeans_1.7.5.tar.gz"

if(!require(tidyverse)){install.packages("tidyverse");require(tidyverse)}
if(!require(mgcv)){install.packages("mgcv");require(mgcv)}
if(!require(emmeans)){install.packages(packageurl, repos=NULL,type="source");require(emmeans)}
if(!require(parallel)){install.packages("parallel");require(parallel)}
if(!require(ggplot2)){install.packages("ggplot2");require(ggplot2)}
if(!require(broom)){install.packages("broom");require(broom)}
if(!require(rgdal)){install.packages("rgdal");require(rgdal)}
if(!require(dplyr)){install.packages("dplyr");require(dplyr)}
install.packages('sf')
library("dplyr")
library("magrittr")
data <- read_delim("Data/lineages_by_ltla_and_week_4_12.tsv")
bounds<-read.csv("Data/Local_Authority_Districts_(December_2020)_UK_BUC.csv", stringsAsFactors=TRUE)
data$day<-as.numeric(as.Date(data$WeekEndDate)-as.Date(data$WeekEndDate[1])+1)
d2 <- data[data$Lineage=="B.1.1.7"|data$Lineage=="B.1.617.2",]
d2=d2%>%mutate(LTLA = as.factor(LTLA))
d2=d2%>%mutate(Lineage = as.factor(Lineage))
d2$lat <- bounds$LAT[match(d2$LTLA,bounds$LAD20CD)]
d2$lon <- bounds$LONG[match(d2$LTLA,bounds$LAD20CD)]
d2<-d2%>%group_by(day, Lineage)



#### NB: Can also use parallel versions using below and replacing "gam" function with "bam"
# specifying "nthreads=3" as an extra argument in the gam functions below
# Further speedups can be obtained for some models using 'discrete=T' argument with 'method=ffREML",discrete=T,'
# See help(bam) for more details
#require(parallel)  
#nc <- 2   ## cluster size, set for example portability
#if (detectCores()>1) { ## no point otherwise
#  cl <- makeCluster(nc) 
#  ## could also use makeForkCluster, but read warnings first!
#} else cl <- NULL

#

### Model S

cov_modS <- bam(Count ~ t2(lat, lon, day, Lineage, bs=c("tp", "tp", "tp", "re"),
                           k=c(10, 10, 10, 6), m=2, full=TRUE),
                data=d2, method="fREML",discrete=T, family="nb",nthreads=3)

#### Model I
  
cov_modI <- bam(Count ~ Lineage + te(lat, lon, day, by=Lineage,
                                       bs=c("tp", "tp", "tp"), k=c(10, 10, 10), m=2),
                  data=d2, method="fREML", family="nb",nthreads=3)


gam.check(cov_modS)
gam.check(cov_modI)
# 
k.check(cov_modS)
k.check(cov_modI)

summary(cov_modS)
summary(cov_modI)

#if(cdarg[1]=="G"){
  cl <- makeCluster(3)
  ### Model G
  cov_modG <- bam(Count ~ t2(lat, lon, day, bs=c("tp", "tp", "tp"), k=c(10, 10, 10)),
                  data=d2, method="fREML",discrete=T, family="nb",nthreads=3)
#}

#if(cdarg[1]=="GS"){
  ### Model GS
  
  cov_modGS <- bam(Count ~ t2(lat, lon, day, bs=c("tp", "tp", "tp"),
                              k=c(10, 10), m=2) +
                     t2(lat, lon, day, Lineage, bs=c("tp", "tp", "tp", "re"),
                        k=c(10, 10, 10, 6), m=2, full=TRUE),
                   data=d2, method="fREML", family="nb")
#}

#if(cdarg[1]=="GI"){
  cl <- parallel::makeCluster(3)
  ### Model GI
  cov_modGI <- bam(Count ~ Lineage+
                     t2(lat, lon, day, bs=c("tp", "tp", "tp"), k=c(10, 10, 10), m=2) +
                     te(lat, lon, day, by=Lineage, bs= c("tp", "tp", "tp"),
                        k=c(10, 10, 10), m=1),
                   data=d2, method="fREML", family="nb",cluster=cl)
#}

#save(list=ls(),file=paste0("./Oct_22_res/Omicron_GI.Rdata"),compress=T)


#gam.check(cov_modG)
#gam.check(cov_modGS)
gam.check(cov_modGI)
# 
#k.check(cov_modG)
#k.check(cov_modGS)
k.check(cov_modGI)

AIC(cov_modS,cov_modI,cov_modG,cov_modGS,cov_modGI)



# ### Calculate (code adapted from Davies et al., 2021, Science)
# 
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

# Specify specific lat/longs of all LTLAs
ull <- na.omit(unique(d2[,c("lat","lon")]))

output <- NULL
g_values <- c(Delta=5.67, Alpha=5.87)
#g_values <- c(Delta=4.7, Alpha=5.5)

# For each LTLA calculate relative multiplicative growth rate factor and
for(i in 1:nrow(ull)){
  
  #Get trends for B.1.617.2 v alpha
  mfit_emtrends2 <- emtrends(cov_modGI, ~ Lineage*lat*lon, "day", mode="latent",at=list(lat=unlist(ull[i,1]),lon=unlist(ull[i,2]),day=seq(as.numeric(as.Date("2021-05-05",format="%Y-%m-%d")-as.Date(d2$WeekEndDate[1])),as.numeric(as.Date("2021-07-17",format="%Y-%m-%d")-as.Date(d2$WeekEndDate[1])), by=2)))
  my_contrast <- list("Delta_vs_Alpha" = c(
    Delta = g_values["Delta"],
    Alpha = -g_values["Alpha"]
  ))
  
  mfit_contrasts2 <- data.frame(
    as.data.frame(contrast(mfit_emtrends2, method=my_contrast)),
    as.data.frame(confint(contrast(mfit_emtrends2, method=my_contrast)))[,c("lower.CL","upper.CL")]
  )
  
  colnames(mfit_contrasts2)[which(colnames(mfit_contrasts2) %in% c("estimate","lower.CL","upper.CL"))] =
    c("delta_rg","delta_rg.lower.CL","delta_rg.upper.CL")

  
  mfit_contrasts = data.frame(mfit_contrasts2[,c("contrast","SE","t.ratio","p.value")],
                              exp(mfit_contrasts2[,c("delta_rg","delta_rg.lower.CL","delta_rg.upper.CL")]))
  
  
  output <- rbind(output,mfit_contrasts)
  
}

# Output estimate of mulp. rate assuming above generation times
out2comp <- matrix(output$delta_rg,nrow=1,byrow=T)


latlon <- expand.grid(unlist(ull[,1]),unlist(ull[,2]))
latlon <- latlon[order(latlon[,2]),]

llind <- sapply(1:311,function(i)which(latlon[,1]==unlist(ull[i,1])&latlon[,2]==unlist(ull[i,2])))


#### PLOTS #####
library(sf)
# Read in England shape file; downloaded from https://borders.ukdataservice.ac.uk/easy_download_data.html?data=England_ct_2011
my_spdf <- st_read(
  dsn="Data/england/england_ct_2011.shp", quiet=TRUE
)

summary(my_spdf)

#Transform projection into lat/long to match LTLA data
my_spdf <- st_transform(my_spdf, crs = 4326)

spdf_fortified <- tidy(my_spdf)


# Plot shparefile
names(latlon) <- c("lat","lon")

plt <- ggplot(my_spdf) +
  geom_sf(fill = "#69b3a2", color = "white") +
  coord_sf() +
  theme_minimal()

#Add points for LTLAs where B.1.617.2 vs Alpha is greater than 0.5
plt3 = plt +
  geom_point(data = latlon[llind,][which(out2comp[1,]>0.5),], aes(x = lon, y = lat, size = out2comp[1,which(out2comp[1,]>0.5)],col= out2comp[1,which(out2comp[1,]>0.5)]),
             shape = 20) +
  scale_fill_gradient(name="Multiplicative Rt",low = "yellow", high = "red", na.value = NA,aesthetics = "colour")  +
  labs(title="B.1.617.2 vs Alpha")+ guides(size="none") +
  theme_void()

ggsave('B.1.617.2valphavaryingtimes.pdf',plt3)

save(cov_modGI, file="B.1.617.2valphaGI.rdata")
save(AIC, file="B.1.617.2valphaAIC.rdata")
save(out2comp, file="B.1.617.2valphaout2comp.rdata")