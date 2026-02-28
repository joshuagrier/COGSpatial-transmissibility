library(readxl)
datadownload <- read_excel("~/Nexus365/Joshua Grier - Results/datadownload.xlsx",sheet = 2)



datadownload[,1] <- as.Date(sapply(1:dim(datadownload)[1],function(x)sub(" to .*", "", datadownload[x,1])),
                            format="%d %B %Y")

write.csv(datadownload,file="seropositivity.csv")
