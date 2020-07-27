
# Correlation -------------------------------------------------------------
source("complete.R")
corr <- function(directory, treshold = 0){
  all <- complete(directory)
  fsample <- all[all$NOBS>treshold, ]
  cor <- vector("numeric")
  for (i in fsample$ID) {
    if (i<10){
      x_i <- paste("./", directory, "/", "00", i, ".csv", sep = "")
      y_i <- read.csv(x_i)
    } else if (i<100){
      x_i <- paste("./", directory, "/", "0", i, ".csv", sep = "")
      y_i <- read.csv(x_i)
      
    } else {
      x_i <- paste("./", directory, "/", i, ".csv", sep = "")
      y_i <- read.csv(x_i)
      
    }
    cor[i] <- cor(y_i$sulfate,y_i$nitrate, use="complete.obs")
  }
  cor [complete.cases(cor)]
  
}

'cr <- corr("specdata", 150)
head(cr)
summary(cr)

cr <- corr("specdata", 400)
head(cr)
summary(cr)
cr <- corr("specdata", 5000)
summary(cr)
length(cr)
cr <- corr("specdata")
summary(cr)
length(cr)'