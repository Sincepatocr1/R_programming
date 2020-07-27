
# Assignment 3 ------------------------------------------------------------
best <- function(state="character", outcome="character" ) {
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
  new1[[c(1,1)]]

}

#Test
# best("TX", "heart attack")
# best("TX", "heart failure")
# best("MD", "heart attack")
# best("MD", "pneumonia")
# best("BB", "heart attack")
# best("NY", "hert attack")
