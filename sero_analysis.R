

##############

## Script file for running HGAMS and calculating relative transmissibility

##############

source('sero_data_clean.R')

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
library("dplyr")
library("magrittr")
data <- read_delim("Data/lineages_by_ltla_and_week.tsv")
bounds<-read.csv("Data/Local_Authority_Districts_(December_2020)_UK_BUC.csv", stringsAsFactors=TRUE)
data$day<-as.numeric(as.Date(data$WeekEndDate)-as.Date(data$WeekEndDate[1])+1)
# d2 <- data[data$Lineage=="AY.4"|data$Lineage=="BA.1.1",]
data['Lineage'][data['Lineage']=='B.1.617.2'|data['Lineage']=='AY.4']<-'delta'
d2 <- data[data$Lineage=="B.1.1.7"|data$Lineage=="delta"|data$Lineage=="BA.1.1"|data$Lineage=="B.1.177",]
d2=d2%>%mutate(LTLA = as.factor(LTLA))
d2=d2%>%mutate(Lineage = as.factor(Lineage))
d2$lat <- bounds$LAT[match(d2$LTLA,bounds$LAD20CD)]
d2$lon <- bounds$LONG[match(d2$LTLA,bounds$LAD20CD)]
d2 <- d2%>%group_by(WeekEndDate,day,Lineage)%>%summarise(total=sum(Count))
d2 <- ungroup(d2)




by <- join_by(closest(WeekEndDate >= `179 ng/ml antibody level, England`))
d2 <- left_join(d2, datadownload, by)%>%drop_na()
d2 <- d2%>%
  mutate(sero=rowMeans(select(d2, starts_with("Age ")), na.rm = TRUE))%>%
  select(!starts_with("Age "))


d2<-d2%>%group_by(day, Lineage)

cov_modGI <- gam(total ~ Lineage+s(day)+s(sero)+
                   t2(sero, day, bs=c("tp", "tp"),by=Lineage, k=c(10, 10), m=2),
                 data=d2, method="fREML", family="nb",nthreads = 3)

cov_modGI2 <- gam(total ~ Lineage+s(day)+s(sero),
                 data=d2, method="fREML", family="nb",nthreads = 3)

cov_modGI3 <- gam(total ~ Lineage+s(sero),
                  data=d2, method="fREML", family="nb",nthreads = 3)

library(marginaleffects)
library(ggplot2)

theme_set(theme_sjplot())



plot_predictions(cov_modGI, by = c("sero", "Lineage"))+xlab("Seropositivity")+ylab("Positive tests")



plot(cov_modGI)

