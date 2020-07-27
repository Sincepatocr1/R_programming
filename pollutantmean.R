
# pollutantmean -----------------------------------------------------------

pollutantmean <- function(directory, pollutant, id = 1:332 ) {
  lid <- length(id)
  
  if (id[1]<10){
    x_i <- paste("./", directory, "/", "00", id[1], ".csv", sep = "")
    y_i <- read.csv(x_i)
    z_i <- y_i 
      
  } else if (id[1]<100){
    x_i <- paste("./", directory, "/", "0", id[1], ".csv", sep = "")
    y_i <- read.csv(x_i)
    z_i <- y_i
      
  } else {
    x_i <- paste("./", directory, "/", id[1], ".csv", sep = "")
    y_i <- read.csv(x_i)
    z_i <- y_i
  }
  
  if (lid>1) {
    
    for (i in id[2]:id[lid]) {
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
    
    z_i <- rbind(z_i, y_i)
    }
  }

  mean(z_i[[pollutant]], na.rm = T)
}
'
pollutantmean("specdata", "sulfate", 1:10)
pollutantmean("specdata", "nitrate", 70:72)
pollutantmean("specdata", "nitrate", 23)'

