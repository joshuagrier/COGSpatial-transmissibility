if(!require(ape)){
    install.packages('ape');
    library(ape)
}


load('...') #change to XXXout2comp.Rdata file, where XXX is variants to compare


#Load relevant data for lat/longs
data <- read_delim("lineages_by_ltla_and_week.tsv")
bounds<-read.csv("Local_Authority_Districts_(December_2020)_UK_BUC.csv", stringsAsFactors=TRUE)
data$day<-as.numeric(as.Date(data$WeekEndDate)-as.Date(data$WeekEndDate[1])+1)

#Change next three lines to relevant lineages
data['Lineage'][data['Lineage']=='B.1.617.2'|data['Lineage']=='AY.4']<-'delta'
data['Lineage'][data['Lineage']=='BA.1.1'|data['Lineage']=='BA.2']<-'omicron'
d2 <- data[data$Lineage=="delta"|data$Lineage=="omicron",]

d2 <- d2%>%group_by(WeekEndDate, LTLA, day, Lineage) %>% summarise_at(vars(Count), sum) %>%
  ungroup()
d2=d2%>%mutate(LTLA = as.factor(LTLA))
d2=d2%>%mutate(Lineage = as.factor(Lineage))
d2$lat <- bounds$LAT[match(d2$LTLA,bounds$LAD20CD)]
d2$lon <- bounds$LONG[match(d2$LTLA,bounds$LAD20CD)]
d2<-d2%>%group_by(day, Lineage)

ull <- na.omit(unique(d2[,c("lat","lon")]))
latlon <- expand.grid(unlist(ull[,1]),unlist(ull[,2]))
latlon <- latlon[order(latlon[,2]),]

#change below depending on which row of ot2comp is appropriate
moridat <- cbind(latlon[llind,],out2comp[1,]) 
names(moridat)[3]<-"ratio"


#Calcualte distance matrix using Euclidian distance
xx=as.matrix(dist(latlon[llind,1:2]))

xx <- 1/xx
diag(xx) <- 0

Moran.I(moridat[,3], xx)

