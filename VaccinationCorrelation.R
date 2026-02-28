# Get names of LTLAs which we have latitude and longitude for
ltlas <- na.omit(unique(d2[,c("lat", "lon", "LTLA")]))
vaccinations <- read.csv("LTLAvaccinations.csv")
# set the week of interest for alpha vs delta weeks 1 - 10, omicron vs delta 37 - 45
# B.1.617.2 vs Alpha 1 - 10, BA.1.1 vs AY.4 36 - 44
for (k in range(44, 44)){
# find the vaccination data for each number of doses
D1 = NULL
D2 = NULL
D3 = NULL
vacweek = k
# iterate over each LTLA of interest
i=1
for (j in ltlas$LTLA){
  D1[i] = vaccination[[j]][["D1"]][vacweek]
  D2[i] = vaccination[[j]][["D2"]][vacweek]
  D3[i] = vaccination[[j]][["D3"]][vacweek]
  i = i+1
}

# format and combine the data
newltla <- as.list(ltlas[['LTLA']])
newout2comp <- t(out2comp)
jointlist <- combineLists(list(newltla, newout2comp, D1, D2, D3))

col_names = list('LTLA','Transmissibility','D1','D2','D3')
colnames(jointlist) <- c(col_names)

input <- data.frame(jointlist)

x <- as.numeric(input[,2])
y1 <- as.numeric(input[,3])
y2 <- as.numeric(input[,4])
y3 <- as.numeric(input[,5])
# model1 <- lm(Transmissibility ~ D1, data = input)
print('Dose 1')
print(cor.test(x,y1))
print('Dose 2')
print(cor.test(x,y2))
print('Dose 3')
print(cor.test(x,y3))
# model2 <- lm(Transmissibility ~ D2, data = input)
# summary(model2)
# model3 <- lm(Transmissibility ~ D3, data = input)
# summary(model3)

}


## function for combining lists
combineLists <- function(manyLists){
  library(plyr)
  newLists <- list()
  for(ixList in 1:length(manyLists)){
    tmpList <- lapply(manyLists[[ixList]], paste, sep = "", collapse = ", ")
    tmpVec  <- as.character(tmpList)
    newLists[[ixList]] <- tmpVec
  }
  newDF   <- t(ldply(newLists))
  return(newDF)
}
