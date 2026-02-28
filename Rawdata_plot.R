if(!require(tidyverse)){install.packages("tidyverse");require(tidyverse)}
if(!require(mgcv)){install.packages("mgcv");require(mgcv)}
if(!require(emmeans)){install.packages(packageurl, repos=NULL,type="source");require(emmeans)}
if(!require(parallel)){install.packages("parallel");require(parallel)}
if(!require(ggplot2)){install.packages("ggplot2");require(ggplot2)}
if(!require(broom)){install.packages("broom");require(broom)}
if(!require(rgdal)){install.packages("rgdal");require(rgdal)}
if(!require(dplyr)){install.packages("dplyr");require(dplyr)}
data <- read_delim("lineages_by_ltla_and_week.tsv")
#data$day<-as.numeric(as.Date(data$WeekEndDate)-as.Date(data$WeekEndDate[1])+1)
data['Lineage'][data['Lineage']=='B.1.617.2'|data['Lineage']=='AY.4']<-'Delta'
data['Lineage'][data['Lineage']=='BA.1.1'|data['Lineage']=='BA.2']<-'Omicron'
data['Lineage'][data['Lineage']=='B.1.177']<-'Wild'
data['Lineage'][data['Lineage']=='B.1.1.7']<-'Alpha'
d2 <- data[data$Lineage=="Delta"|data$Lineage=="Omicron"|data$Lineage=="Alpha"|data$Lineage=="Wild",]
d2 <- d2%>%group_by(WeekEndDate, LTLA, Lineage) %>% summarise_at(vars(Count), sum) %>%
  ungroup()
d2=d2%>%mutate(LTLA = as.factor(LTLA))
d2=d2%>%mutate(Lineage = as.factor(Lineage))
d2<-d2%>%group_by(WeekEndDate, Lineage)%>%summarise(Count=sum(Count))

p1 <- ggplot(d2)+geom_line(aes(x=WeekEndDate, y=Count, colour=Lineage),linewidth = 1.)
ggsave(filename = "Rawdata.pdf",plot=p1)
