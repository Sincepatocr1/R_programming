rankall<- function(outcome="character", num="best") {
  out <- read.csv("outcome-of-care-measures.csv", colClasses = "character")
  if (!outcome %in% c("heart attack", "heart failure", "pneumonia")) {stop("invalid outcome")}
  outstate <- split(out, out$State)
  results <- data.frame(hospital = vector("character", length = length(outstate)), state = vector("character", length = length(outstate)))
   if (outcome=="heart attack") {
    col <- 11
  } else if (outcome=="heart failure") {
    col <- 17
  } else {
    col <- 23
  }
  values <- sort(unique(out$State))
  j <- 1
  for (i in values) {
    outstateuse <- outstate[[i]]
    outstateuse[, col] <- as.numeric(outstateuse[, col])
    new <- outstateuse[, c(2, col)]
    good <- complete.cases(new)
    new1 <- new[good, ]
    new1 <- new1[order(new1[, 2], new1[, 1] ),]
    last <- nrow(new1)
    
    if (num=="best") {
      results[j, 1] <- new1[[c(1,1)]]
      results[j, 2] <- i
    }
    else if (num=="worst") {
      
      results[j, 1] <- new1[[c(1,last)]]
      results[j, 2] <- i
      
    }
    else if (num>last) {
      results[j, 1] <- NA
      results[j, 2] <- i
      
    }
    else {
      
      results[j, 1] <- new1[[c(1,num)]]
      results[j, 2] <- i
    }
    j <- j+1
  }
  results
}
head(rankall("heart attack", 20), 10)
tail(rankall("pneumonia", "worst"), 3)
tail(rankall("heart failure"), 10)