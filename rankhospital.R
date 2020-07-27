# Assiignment 3 -----------------------------------------------------------

rankhospital <- function(state="character", outcome="character", num="best") {
  out <- read.csv("outcome-of-care-measures.csv", colClasses = "character")
  if (!state %in% out$State) {stop("invalid state")}
  else if (!outcome %in% c("heart attack", "heart failure", "pneumonia")) {stop("invalid outcome")}
  outstate <- split(out, out$State)
  outstateuse <- outstate[[state]]
  
  if (outcome=="heart attack") {
    col <- 11
  } else if (outcome=="heart failure") {
    col <- 17
  } else {
    col <- 23
  }
  outstateuse[, col] <- as.numeric(outstateuse[, col])
  new <- outstateuse[, c(2, col)]
  good <- complete.cases(new)
  new1 <- new[good, ]
  new1 <- new1[order(new1[, 2], new1[, 1] ),]
  
  last <- nrow(new1)
  
  if (num=="best") {
    new1[[c(1,1)]]
  }
  else if (num=="worst") {
    
    new1[[c(1, last)]]
  }
  else if (num>last) {
    
    return(NA)
  }
  else {
    
    new1[[c(1, num)]]
  }
}


# rankhospital("TX", "heart failure", 4)
# rankhospital("MD", "heart attack", "worst")
# rankhospital("MN", "heart attack", 5000)