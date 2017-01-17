# Load data from csv file and split into development and internal validation
# datasets

DeBois <- function(Ht, Wt){
  0.007184 * Ht^0.725 * Wt^0.425
}

thedata <- read.csv("GFR_data_clean.csv")

thedata <- thedata[thedata$Age >= 18,]

thedata_full <- thedata

set.seed(1234)
val_index <- sample(1:2471, 2471/5)
thedata_val <- thedata[val_index,]
thedata_dev <- thedata[-val_index,]



# Load external validation dataset
library(xlsx)
GlasgowData <- read.xlsx("SeminomaGFRdataGlasgow.xls",sheetIndex = 1)

# selecting relavant variables
CleanGlasgowData <- GlasgowData[,c(3,4,6,7,14,26)]

# remove rows with missing values, match names to original dataset, add SufA 
# variable 
CleanGlasgowData <- na.omit(CleanGlasgowData)
CleanGlasgowData$Height <- CleanGlasgowData$Height*100
names(CleanGlasgowData) <- c("Age", "Sex", "Ht", "Wt", "Creat", "GFR")
CleanGlasgowData$SufA <- DeBois(CleanGlasgowData$Ht, CleanGlasgowData$Wt)

