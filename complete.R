
# Complete ----------------------------------------------------------------

complete <- function(directory, id=1:332) {
  data <- data.frame(ID=vector("numeric"), NOBS=vector("numeric" ) )
  for (i in id) {
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
    
    z <- complete.cases(y_i)
    nrows <- nrow(y_i[z, ])
    data[i,1] <- i
    data[i,2] <- nrows
  }
  
  data[complete.cases(data), ]
}
'
complete("specdata", 1)
complete("specdata", c(2, 4, 8, 10, 12))
complete("specdata", 30:25)
complete("specdata", 3) '
