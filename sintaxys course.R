getwd()
dir()
myfunction <- function() {
  x <- rnorm(100)
  mean(x)
}

ls()

#rm("myfuction")

x <- 1:50
class(x)
as.numeric(x)
as.logical(x)
as.character(x)
as.complex(x)

y <- vector("character", length = 3)

rm("x")

x <- list(3, "hola", FALSE, 4 - 2i)
x

m <- matrix(nrow = 3, ncol = 4)
m
dim(m) #dimesion is a vector itself with two elements
attributes(m)

mat <- matrix(1:10, nrow = 5, ncol = 2) #it puts elements by columns
mat

#createa matrix by vectors, giving a dimension

rm("y")
y <- 1:6
dim(y) <- c(2, 3)
y

#Matrix by cbind, rbind

rm("x", "y")

x <- 1:4
y <- 8:11

cbind(x, y)
matrix <- rbind(x, y)
matrix

#factor (categorical variable)

rm("x", "y")

x <- factor(c("no", "yes", "no", "no", "yes"))
x
unclass(x)

x <-  factor(c("no", "yes", "no", "no", "yes"), levels = c("yes", "no"))
levels(x)

#Missing values

x <- c(4, 2, 3, NA, NaN)
is.na(x)
is.nan(x)

#Data frame
x <- data.frame(logic = c(T, F, T, T), ID = 1:4)
x
nrow(x)
ncol(x)

#names attribute

x <- 1:4
names(x) <- c("Hola", "Fabiola", "Richard", "Goku")
x
y <- list(a = 1, b = 2, c = "Hola")
y
mat <- matrix(1:6, nrow = 2, ncol = 3) #filling column by column
mat
dimnames(mat) <- list(c("a", "b"), c("a", "b", "c"))
mat

#reading larga data sets

#in<-read.table("datatablename.txt", nrows=50)
#clases<-sapply(in, class)
#all<-read.table("datatablename.txt", colClasses = clases)  #colClases="numeric", etc.
help("read.table")

#using dput
y <- data.frame(a = 10:19, b = 1:10)
dput(y)
dput(y, file = "y.R")
newy = dget("y.R")
newy

#dump - multiple objects
dump(c("newy", "matrix"), file = "data.R")
rm(matrix, newy)
source("data.R")

#Conection outside world

str(file) #describe the function

con <- file("name.txt", "r")
data <- read.csv(con)
close(con)

#is the same as data<-read.csv("name.txt"), but is usefull for reading some lines, e.g.

#con<-gzfile("words.gz")
#x<-readLines(con, 10)
#x

#save<-url("https://www.caces.gob.ec")
#x<-readLines(save)
#x


#Subsetting basics
rm(list = ls())
x <- c("a", "b", "c", "d", "a")
x[1]
x[1:3]
x[x > "a"]
u <- x > "a"
u
x[u]

#Subseting list
x <- list(var1 = 1:4,
          var2 = "Hola",
          var3 = T)
x[1]
x[[1]]
x$var1
x["var2"] #list
x[["var2"]] #element
x[c(1, 3)]

namex <- "var1"
x[namex]

#List to extract an interger sequence

x <- list(a = list(1, 2, 3), b = c(1.5, 6.7))
x[[c(1, 3)]] #equivalent to
x[[1]][[3]]
x[[c(2, 2)]]


##Subsetting matrices

xx <- matrix(1:8, 2, 4)
xx
xx[2, 4]
xx[1,]
xx[, 2]
xx[1:2, 3:4]

#We can put off the drop dimention option

xx[2, 4, drop = F]
xx[, 2, drop = F]


#Partial matching

x <- list(ostia = c(1, 2, 3, 7, 0))
x$o
x[["o", exact = F]]
y <- c(1, 2, 4, 5, NA, NaN)
bad <- is.na(y)
bad
y[!bad]
x <- c("a", NA, "b", NA, "c", "d")
good <- complete.cases(x, y)
x[good]

#complete cases in a data frame
y[1:6, ] #show a number of cases
bien <- complete.cases(y)
y[bien,][1:7,]

#vectorized operations

x <- 1:5
y <- 5:9
x + y
x * y
y > 2
y == 2
y >= 7
x / y

#matrix operations

x <- matrix(1:4, 2, 2)
y <- matrix(rep(10, 4), 2, 2)
#element by element
x * y
x / y
x + y
x %*% y #real matriz multiplication

#Quiz

x <- c(3, 5, 1, 10, 12, 6)
x[x <= 5] <- 0
x[x %in% 1:5] <- 0
x[x < 6] <- 0

dir()

#Data
x <- read.csv("quiz1_data/hw1_data.csv")
x
x[1:2,]
x[152:153, ]
x[47, 1]
misozone <- is.na(x[, 1])
is.na(x$Ozone)
x[!misozone, 1]
length(misozone[misozone == T])

mean(x[!misozone, 1])
mean(x$Ozone, na.rm = T)

y <- x[which(x$Ozone > 31 & x$Temp > 90), ]
mean(y$Solar.R)
y <- subset(x, x$Month == 6, )
mean(y$Temp)
help("max")
max(x[which(x$Month == 5), "Ozone"], na.rm = T)
nrow(x)
tail(x, 2) #extract two last rows

restriction <- c(x$Ozone > 31 & x$Temp > 90)
mean(x[restriction, "Solar.R"], na.rm = T)



###Control Structurres: If, else

if (x > 3) {
  y <- 10
} else {
  y <- 0
}
y <- if (x < 3) {
  10
} else {
  0
}

if (x > 3) {
  y <- 10
} else if (x < 50) {
  y <- 0
} else {
  y <- NA
  
}

#for loop

for (i in 1:10) {
  print(i)
}

x <- c("a", "b", "c", "d")

#Equivalent
for (i in 1:4) {
  print(x[i])
}

for (i in seq_along(x))
  #create a vector equal to the length of the vector
{
  print(x[i])
}


for (letters in x) {
  print(letters)
}

for (i in 1:4)
  print(x[i])

x <- matrix(1:20, 5, 4)

for (i in seq_len(nrow(x)))
  #desired length
{
  for (j in seq_len(ncol(x))) {
    print(x[i, j])
  }
}


##While

count <- 0
while (count < 10) {
  print(count)
  count <- count + 1
}

z <- 5
while (z >= 3 && z <= 10) {
  print(z)
  coin <- rbinom(1, 1, 0.5)
  
  if (coin == 1) {
    ##randomm walk
    z <- z + 1
    
  } else {
    z <- z - 1
  }
}

##CONTROL STRUCTURES

#Repeat example (it may not break - not converge)
x0 <- 1
tol <- 1e-8

repeat {
  x1 <- computeEstimate()  #the function does not exist
  
  if (abs(x1 - x0) < tol) {
    break
  } else{
    x0 <- x1
  }
  
}

##next: skip

for (i in 1:100) {
  if (i <= 20) {
    next
  }
  return(i) #exit the function and return a given value
}

##FUNCTIONS

suma2 <- function(x, y) {
  x + y
}
suma2(5, 10)

sobre <- function(x, y)
  x[x > y]
X <- 1:20
sobre(X, 16)

#equivalent

sobre <- function(x, y) {
  use <- x > y
  x[use]
}

sobre(X, 15)

colummeans <- function(x, removeNA = T) {
  ncol <- ncol(x)
  means <- numeric(ncol)
  for (i in 1:ncol) {
    means[i] = mean(x[, i], na.rm = removeNA)
  }
  means
}
colummeans(airquality)
colummeans(airquality, F)

#Argument matching by position or by name, the following are wquivalent:

mydata <- rnorm(100)
sd(mydata)
sd(x = mydata)
sd(x = mydata, na.rm = F)
sd(na.rm = F, x = mydata)
sd(na.rm = F, mydata)

y <- rnorm(100)

args(lm)

#lm(data = mydata, y ~ x, model=F, 1:100)
#lm(y ~ x, mydata, 1:100, model=F )

f <- function(a,
              b = 1,
              c = 2,
              d = NULL)
  #is the same NULL or nothing
  #Lazy evaluation
  f <- function(a, b) {
    a ^ 2
  }
f(2)

f <- function(a, b) {
  print(a)
  print(b)
}
f(45) #error only after print(a)

#"..." argument indicated a variable number of arguments that are usually passed on to other functions, is often used
#when extending another function and you don't want to copy the entire argument list of the original function

myplot <- function(x, y, type = "l", ...) {
  plot(x, y, type = type, ...)
}

formals(myplot) #tp know the formal arguments of a function
#Also used in generic functions  so that extra arguments can be passed for methods
mean

#Also when the number of arguments passed to the function can not be kown in advance
args(paste)
args(cat) #concaten arguments

#Arguments after "..." must be named explicitly
paste("a", "b", sep = ":")
paste("a", "b", se = ":") # no partial matching


search() #R look first at global enviroment (my enviroment). If I define a new funtion with the same name of
#an existing function R took first my function.

#When load a package with library, R puts it in the Second position of the search()


##SCOPING RULES R lexical or static scoping, others dynamic scoping

f <- function(x, y) {
  x + y / z  #z is a free variable, lexical scoping in R means that the value of free variable is searched in the enviroment
  #where  the function was defined
}
##Environment is a collection of (symbol, value) pairs "x" symbol, "3.14" its value
##Every environment have a parent environment "Children"
##Empty enviroment = environment without a parent
##A function + an environment = a closure or function clousure
##To know the value of a free variable first search in the environment where the function was defined and then , to the parents
##until base package /after is the empty environment)

#Example

make.power <- function(n) {
  pow <- function(x) {
    x ^ n
  }
  pow #retunrs a function
}

cube <- make.power(3) #function x^3
square <- make.power(2) # function x^2

cube(3)
square(3)

#Ypou can see what is in a function enviaronment with ls

ls(environment(cube))
get("n", environment(cube))
get("pow", environment(cube))

##Lexical vs dynamic scoping

y <- 10

f <- function(x) {
  y <- 2
  y ^ 2 + g(x)
}

g <- function(x) {
  x * y
}
f(3) #in lexical y=10 in function g *global environment, in dynamic the value of y is 2
##### because value is looked up in the environment for wich the function was called (calling environment)

#f(3) is 2^2+3*10 lexical
#f(3) is 2^2+3*2 dynamic
g(3)

##Optimization> functions in R> optim, nlm, optimize


#constructor function - Negatuve likelihood

make.NegLokLik <- function(data, fixed = c(F, F)) {
  params <- fixed
  function(p) {
    params[!fixed] <- p
    mu <- params[1]
    sigma <- params[2]
    a <- -0.5 * length(data) * log(2 * pi * sigma ^ 2)
    b <- -0.5 * sum((data - mu) ^ 2) / (sigma ^ 2)
    - (a + b)
  }
}

fixed = c(1, F)
params <- fixed
params[!fixed]

rm(fixed, params)

##simulate some normal random variables, mean 1, sd 2
set.seed(1)
normals <- rnorm(100, 1, 2)

##create a negative likelihood function>

nLL <- make.NegLokLik(normals)
nLL

ls(environment(nLL))

##Estimating parameters

args(optim)
optim(c(mu = 0, sigma = 1), nLL)$par
help(optim)

#fixing sigma=2
nLL <- make.NegLokLik(normals, c(F, 2))
optimize(nLL, c(-1, 3))$minimum

#fixing mu=1
nLL <- make.NegLokLik(normals, c(1, F))
optimize(nLL, c(1e-6, 10))$minimum



#Plotting the likelihood

nLL <- make.NegLokLik(normals, c(1, F))
x <- seq(1.7, 1.9, len = 100)
y <- sapply(x, nLL) #apply x to the function nLL
plot(x, exp(-(y - min(y))), type = "l")


nLL <- make.NegLokLik(normals, c(F, 2))
x <- seq(0.5, 1.5, len = 100)
y <- sapply(x, nLL) #apply x to the function nLL
plot(x, exp(-(y - min(y))), type = "l")



# Dates and times ---------------------------------------------------------

x <- as.Date("1970-01-01")
x 
unclass(x)

x <- as.Date("1970-01-02")
unclass(x)

x <- Sys.time()
x# Already in "POSIXct" format
unclass(x)
p <- as.POSIXlt(x)
names(unclass(p))

p$zone

datestring <- c("Enero 10, 1958 10:50", "Julio 23, 1961 12:15")

x <- strptime(datestring, "%B %d, %Y %H:%M")
x
class(x)
?strptime
x <- as.Date("1968-01-02")
y <- strptime("Enero 10, 1958 10:50", "%B %d, %Y %H:%M")

x <- as.POSIXlt(x)
x-y ##Also can compare (== or <=)

x <- as.Date("1968-01-02"); y <- as.Date("1968-01-16")
x-y

x <- as.POSIXct("2020-05-10 01:00:00")
y <- as.POSIXct("2020-05-10 07:00:00", tz = "GMT")

x-y

x <- 1:10
if (x>5){
  x <- 0
}

rm(list = ls())

f <- function(x) {
  g <- function(y) {
    y + z
  }
  z <- 4
  x + g(x)
}

z <- 10
f(3)


# Quiz 2 ------------------------------------------------------------------

source("pollutantmean.R")
source("complete.R")
source("corr.R")

pollutantmean("specdata", "sulfate", 1:10)
pollutantmean("specdata", "nitrate", 70:72)
pollutantmean("specdata", "sulfate", 34)
pollutantmean("specdata", "nitrate")

cc <- complete("specdata", c(6, 10, 20, 34, 100, 200, 310))
print(cc$NOBS)

cc <- complete("specdata", 54)
print(cc$NOBS)

RNGversion("3.5.1")  
set.seed(42)
cc <- complete("specdata", 332:1)
?sample
use <- sample(332, 10)
print(cc[use, "NOBS"])

cr <- corr("specdata")                
cr <- sort(cr)   
RNGversion("3.5.1")
set.seed(868)                
out <- round(cr[sample(length(cr), 5)], 4)
print(out)

cr <- corr("specdata", 129)                
cr <- sort(cr)                
n <- length(cr)    
RNGversion("3.5.1")
set.seed(197)                
out <- c(n, round(cr[sample(n, 5)], 4))
print(out)

cr <- corr("specdata", 2000)                
n <- length(cr)                
cr <- corr("specdata", 1000)                
cr <- sort(cr)
print(c(n, round(cr, 4)))


# Lappy, sapply, apply, tapply, mapply, split -----------------------------
rm(list = ls())
#Lapply
?rnorm
x <- list(a=1:5, b=rnorm(10))
lapply(x, mean) #always retunr a list
x <- list(a=1:5, b=rnorm(10), c=rnorm(15, 3, 3), d=rnorm(20, sd=10))
lapply(x, mean)
?runif #generate uniform random variables
x <- 1:4
lapply(x, runif) # apply runif to the sequence (takes 1, 2, 3, 4 as number of random variables that has to generate)
lapply(x, runif, min=0, max=10)

x <- list(a=matrix(1:4, 2, 2), b=matrix(1:8, 4, 2))
x
#anonymus functions
lapply(x, function(elt) elt[,1]) #elt just exist in the context of lapply thatt takes the first column

x <- list(a=1:5, b=rnorm(10), c=rnorm(15, 3, 3), d=rnorm(20, sd=10))
sapply(x, mean) #not always retunr a list
mean(x) #not work with lits

# apply - often used to evaluate a function (an anonymus function) over the margins of an array
## apply a function to the rows or columns of a matrix
## It can be used with general arrays - e.g. taking the avergae of an array of matrices
## is not faster than writing a loop, but it works in one line

x <- matrix(rnorm(200), 20, 10)
x
apply(x, 2, mean) ## 2 represent the margin (dimension) of the array that has to take
apply(x, 1, sum) ## summ returns the sum all the values
# The following are faster:
# rowSums=apply(x, 1, sum)
# rowMeans=apply(x, 1, mean)
# colSums=apply(x, 2, sum)
# colMeans=apply(x, 2, mean)

apply(x, 1, quantile, probs = c(0.25, 0.75))

a <- array(rnorm(2*2*10), c(2,2,10))
a
?array
apply(a, c(1,2), mean)
rowMeans(a, dims = 2)

# mapply(function, ...) - miltivarate apply pf sort which applies a function in paralell over a set of arguments
str(mapply)
list(rep(1,4), rep(2,3), rep(3,2), rep(4,1))
#Instead
mapply(rep, 1:4, 4:1) #apply function to mulptiple sets of arguments

noise <- function(n, mean, sd){
  rnorm(n, mean, sd)
}

noise(5, 1, 2)

noise(1:5, 1:5, 2) # it does not work well
mapply(noise, 1:5, 1:5, 2) # it works with a set of arguments
#it is the same as
list(noise(1,1,2), noise(2,2,2), noise(3,3,2), noise(4,4,2), noise(5,5,2))

#Tapply - apply functions over subsets of a vector

x <- c(rnorm(10), runif(10), rnorm(10, 1))
?gl
f <- gl(3,10) #generate 3 levels and repeat each level 10 times
f
tapply(x, f, mean) #simplfication by default
tapply(x, f, mean, simplify = F)
tapply(x, f, range)

#Split - takes a vector or other objects and splits it into groups determined by a factor or list of factors
split(x, f) #also we can put a logical vector indicates wheter empty factor levels should be dropped
## Always retunr a list

#It is common to see with lappy
lapply(split(x, f), mean)
sapply(split(x,f), mean)
tapply(x, f, mean)

library(datasets)
head(airquality)
#View(airquality)
airquality$Month
s <- split(airquality, airquality$Month)
#View(s)
lapply(s, function(x) colMeans(x[, c("Ozone", "Solar.R", "Wind")])) #anonymus function
sapply(s, function(x) colMeans(x[, c("Ozone", "Solar.R", "Wind")])) #anonymus function
sapply(s, function(x) colMeans(x[, c("Ozone", "Solar.R", "Wind")], na.rm = T)) #anonymus function


#Spliting in more than one level (for expample if you have more than one factor variable: gender and race)

x <- rnorm(10)
f1 <- gl(5, 2)
f2 <- gl(2, 5)

f1
f2
interaction(f1, f2)

#str display, compactly, the structure of an R object
str(split(x, list(f1, f2)))
str(split(x, list(f1, f2), drop = T))


# Debugging tools ---------------------------------------------------------
## Indication that something is not rigth
## mesaage - generic notification/diagnostic messagge produced by the funtions messagge
## warning - indication that smething is wrong but not neccesarily fatal, execution of the function conitues
## error - fatal problem was occurred - exectution stops
## condition - a generic concept for indicating tat something unexpected can occur
rm(list = ls())
log(-1)
printmessahe <- function(x){
  if (x>0){
    print("x is greater than cero")
  } else {
    print("x is lees tha or equal to cero")  
  }
  invisible(x)
}
printmessahe(1)
printmessahe(NA)

printmessahe <- function(x){
  if (is.na(x))
    print("x is a missing value")
  else if (x>0)
    print("x is greater than cero")
  else 
    print("x is lees tha or equal to cero")  
  
  invisible(x) # not have funtion return when is not assigned
}

printmessahe(NA)

##If something is wrong:
## See Input, how did you call the function?
## What were you expecting?
## What did you get?
##  How does what you get different forom what do you expect?
## Were your expectations correct in the first place?
## Can you reproduce the problem (exactly)?

###Debugging - Basic Tools

# tracebakk say how many function you call and where the error is
mean(x)
traceback()
lm(y ~ x)
traceback() # error occurs in 7 level

# debug - you can go or excecute a funtion line by line
debug(lm)
lm(y ~ x)
# bowser - same debug but not necesary at the beggining, you can call browser anywhere in your code
# trace - allows to insert debuggin without actually edditing the function itself
# recover - chande the error behaviour - not return the error in console, stops right where the error occurs

#options(error = stop)
?options
options(default)
read.csv("sisisisi")


# Quiz 3 ------------------------------------------------------------------

library(datasets)
data(iris)
?iris
s <- split(iris, iris$Species)
View(s)
?round
round(mean(s$virginica[, "Sepal.Length"]), digits = 0)
apply(iris, 2, mean)
apply(iris, 1, mean)
rowMeans(iris[, 1:4])
apply(iris[, 1:4], 2, mean)
apply(iris[, 1:4], 1, mean)

library(datasets)
data(mtcars)
?mtcars
sapply(split(mtcars$mpg, mtcars$cyl), mean)
split(mtcars$mpg, mtcars$cyl)
mtcars
tapply(mtcars$mpg, mtcars$cyl, mean)
with(mtcars, tapply(mpg, cyl, mean))

t <- sapply(split(mtcars$hp, mtcars$cyl), mean)
t
round(abs(t["4"]-t["8"]), digits = 0)


# Weeek 4 -----------------------------------------------------------------

# str compactly display the structure of an R object - diagnostic and alternative of summary
str(str)
str(lm)
x <- rnorm(100, 0, 1)
summary(x)
str(x)
f <- gl(40, 10)
str(f)
summary(f)
#View(f)
library(datasets)
str(airquality)
s <- split(airquality, airquality$Month)
str(s)

###Simulation - Probability distributions
## rnorm - random normal variables with a given mean and sd
## dnorm - evaluate the normal probability density (with a giiven mean/sd) at a point (or vector of points)
## pnorm - evaluate the cumuative distribution function for a normal distribution
## rpois - generate random Poisson variates with a given rate

# d - density, r - random, p - cumulative, q - quantile

x <- rnorm(10, 20, 2)
x
summary(x)

# set.seed ensures reproducibility - because they are pseudo random variables
set.seed(1)
rnorm(5)
rnorm(5)

rpois(10, 1) # mean is the rate

ppois(2,2) ##P(x<= 2) if rate is 2
ppois(4,2) ##P(x<= 4) if rate is 2

## Simulating a Linear Model

## Suppose y=b0+b1x+e, where e~N(0, 2^2), x~N(0, 1^2), b0=0.5, b1=2

set.seed(20)
x <- rnorm(100)
e <- rnorm(100, 0, 2)
y <- 0.5+2*x+e
summary(y)
plot(x,y)

#What if x is binary? - like gender, control vs treatment
?rbinom
set.seed(10)
x <- rbinom(100, 1, 0.5) #n=1, p=0.5 (probability of one is 0.5)
e <- rnorm(100, 0, 2)
y <- 0.5+2*x+e
summary(y)
plot(x,y)

#Poisson model where:
#Y~Poisson(u)
#log(u)=b0+b1x
#b0=0.5 b1=0.3

set.seed(1)
x <- rnorm(100)
logmu <- 0.5+0.3*x #error term is not like the other case
y <- rpois(100, exp(logmu))
summary(y)
plot(x,y)

#Random sampling

#sample function draws randomly from a specific set of (scalar) objects allowing to sample from
# arbitrary distributions

set.seed(1)
sample(1:10, 4) #sample of vector 1:10, 4 elements, without replacement
sample(letters, 5)

sample(1:10) #permutation
sample(1:10, replace = T) ##sample with replacement

### R profiler to know why is taking a lot of time, and suggest strategies to fix the problem

#system.time - not profiler - takes an arbitrary expresion as input and returns the amount of
#time taken to evaluate the expression. Computes time in seconds needed to execute expression (if an error occurs
# gives times until error occured)

#return an object of class proc_time
# user time - time charged to the CPU for this expression - time computer experience
# elapsed time - "Wall clock" time - amount of time that you experience

#Elapsed time can be greater if CPU spend time waiting (CPU is been used by other background actions)
#Elapsed time can be smaller if machine has multicores
#  - Multi -threaded BLAS library (vecLib/Accelerate, ATLAS, ACML- AMD Machines, MKL - Intel machines) (linear algebra type of library)
# - Parallel procesing via parallel package

## Elapsed time > user time
system.time(readLines("http://www.jhsph.edu"))

## Elapsed time < user time

hilbert <- function(n) {
  
  i <- 1:n
  1/outer(i-1, i, "+")
  
}

x <- hilbert(100)
system.time(svd(x))
?svd
system.time({
  n <- 1000
  r <- numeric(1000)
  for (u in 1:n) {
    x <- rnorm(u)
    r[u] <- mean(x)
  }
})


##Rprof - starts the profiler in R
## summaryRprof() summarizes the output of Rprof() (otherwise is not readable)
# not use system.time and Rprof togheter
# Rprof keeps track of function call stack and regularly sampled intervals and tabulates
# how much time is spend in each function

# default sampling interval is 0.02 seconds

#y <- rnorm(100)
#x <- rnorm(100)

#lm(y~x)
#sample.interval=1000

#summaryRprof / two methods for nomalizing data by.total, by.self (first divides the time spend by each function
# by tthe totla run time, second do the same but first substracts out time spents in functions above in the call stack)

# by.total how many times that function appears in the calls - but of the top level function
# by.self tells you the time of the help functions that are called from top lvele function
# $by.total
# $by.self
# $sample.interval
# $sampling.time

#rm(list = ls())
#Rprof()
#y <- rnorm(100)
#x <- rnorm(100)
#lm(y~x)

#sample.interval=100
#summaryRprof()

set.seed(1)
rpois(5, 2)

set.seed(10)
x <- rep(0:1, each = 5)
e <- rnorm(10, 0, 20)
y <- 0.5 + 2 * x + e

rm(list = ls())

library(datasets)
Rprof()
y <- rnorm(10)
x1 <- rnorm(10)
x2 <- rnorm(10)
fit <- lm(y ~ x1 + x2)
Rprof(NULL)
summaryRprof()