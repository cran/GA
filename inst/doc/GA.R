## ----setup, include=FALSE------------------------------------------------
library(knitr)
opts_chunk$set(fig.align = "center", 
               out.width = "80%",
               fig.width = 6, fig.height = 5,
               dev.args=list(pointsize=10),
               par = TRUE, # needed for setting hook 
               collapse = TRUE, # collapse input & ouput code in chunks
               warning = FALSE)

knit_hooks$set(par = function(before, options, envir)
  { if(before && options$fig.show != "none") 
       par(family = "sans", mar=c(4.1,4.1,1.1,1.1), mgp=c(3,1,0), tcl=-0.5)
})

## ---- message = FALSE, echo=1--------------------------------------------
library(GA)
cat(GA:::GAStartupMessage(), sep="")

## ------------------------------------------------------------------------
f <- function(x)  (x^2+x)*cos(x)
lbound <- -10; ubound <- 10
curve(f, from = lbound, to = ubound, n = 1000)

GA <- ga(type = "real-valued", fitness = f, lower = c(th = lbound), upper = ubound)
summary(GA)
plot(GA)

curve(f, from = lbound, to = ubound, n = 1000)
points(GA@solution, GA@fitnessValue, col = 2, pch = 19)

## ------------------------------------------------------------------------
Rastrigin <- function(x1, x2)
{
  20 + x1^2 + x2^2 - 10*(cos(2*pi*x1) + cos(2*pi*x2))
}

x1 <- x2 <- seq(-5.12, 5.12, by = 0.1)
f <- outer(x1, x2, Rastrigin)
persp3D(x1, x2, f, theta = 50, phi = 20, col.palette = bl2gr.colors)
filled.contour(x1, x2, f, color.palette = bl2gr.colors)

## ------------------------------------------------------------------------
GA <- ga(type = "real-valued", 
         fitness =  function(x) -Rastrigin(x[1], x[2]),
         lower = c(-5.12, -5.12), upper = c(5.12, 5.12), 
         popSize = 50, maxiter = 1000, run = 100)
summary(GA)
plot(GA)

## ------------------------------------------------------------------------
filled.contour(x1, x2, f, color.palette = bl2gr.colors, 
  plot.axes = { axis(1); axis(2); 
                points(GA@solution[,1], GA@solution[,2], 
                       pch = 3, cex = 2, col = "white", lwd = 2) }
)

## ---- eval=FALSE---------------------------------------------------------
#  monitor <- function(obj)
#  {
#    contour(x1, x2, f, drawlabels = FALSE, col = grey(0.5))
#    title(paste("iteration =", obj@iter), font.main = 1)
#    points(obj@population, pch = 20, col = 2)
#    Sys.sleep(0.2)
#  }
#  
#  GA <- ga(type = "real-valued",
#           fitness =  function(x) -Rastrigin(x[1], x[2]),
#           lower = c(-5.12, -5.12), upper = c(5.12, 5.12),
#           popSize = 50, maxiter = 100,
#           monitor = monitor)

## ------------------------------------------------------------------------
suggestedSol <- matrix(c(0.2,1.5,-1.5,0.5), nrow = 2, ncol = 2, byrow = TRUE)
GA1 <- ga(type = "real-valued", 
          fitness =  function(x) -Rastrigin(x[1], x[2]),
          lower = c(-5.12, -5.12), upper = c(5.12, 5.12), 
          suggestions = suggestedSol,
          popSize = 50, maxiter = 1)
head(GA1@population)

## ------------------------------------------------------------------------
GA <- ga(type = "real-valued", 
         fitness =  function(x) -Rastrigin(x[1], x[2]),
         lower = c(-5.12, -5.12), upper = c(5.12, 5.12), 
         suggestions = suggestedSol,
         popSize = 50, maxiter = 100)
summary(GA)

## ------------------------------------------------------------------------
f <- function(x)
  { 100 * (x[1]^2 - x[2])^2 + (1 - x[1])^2 }

c1 <- function(x) 
  { x[1]*x[2] + x[1] - x[2] + 1.5 }

c2 <- function(x) 
  { 10 - x[1]*x[2] }

## ------------------------------------------------------------------------
ngrid <- 250
x1 <- seq(0, 1, length = ngrid)
x2 <- seq(0, 13, length = ngrid)
x12 <- expand.grid(x1, x2)
col <- adjustcolor(bl2gr.colors(4)[2:3], alpha = 0.2)
plot(x1, x2, type = "n", xaxs = "i", yaxs = "i")
image(x1, x2, matrix(ifelse(apply(x12, 1, c1) <= 0, 0, NA), ngrid, ngrid), 
      col = col[1], add = TRUE)
image(x1, x2, matrix(ifelse(apply(x12, 1, c2) <= 0, 0, NA), ngrid, ngrid), 
      col = col[2], add = TRUE)
contour(x1, x2, matrix(apply(x12, 1, f), ngrid, ngrid), 
        nlevels = 21, add = TRUE)

## ------------------------------------------------------------------------
x <- c(0.8122, 12.3104)
f(x)

## ------------------------------------------------------------------------
c1(x)
c2(x)

## ------------------------------------------------------------------------
fitness <- function(x) 
{ 
  f <- -f(x)                         # we need to maximise -f(x)
  pen <- sqrt(.Machine$double.xmax)  # penalty term
  penalty1 <- max(c1(x),0)*pen       # penalisation for 1st inequality constraint
  penalty2 <- max(c2(x),0)*pen       # penalisation for 2nd inequality constraint
  f - penalty1 - penalty2            # fitness function value
}

## ------------------------------------------------------------------------
GA <- ga("real-valued", fitness = fitness, 
         lower = c(0,0), upper = c(1,13), 
         # selection = GA:::gareal_lsSelection_R,
         maxiter = 1000, run = 200, seed = 123)
summary(GA)

fitness(GA@solution)
f(GA@solution)
c1(GA@solution)
c2(GA@solution)

## ------------------------------------------------------------------------
plot(x1, x2, type = "n", xaxs = "i", yaxs = "i")
image(x1, x2, matrix(ifelse(apply(x12, 1, c1) <= 0, 0, NA), ngrid, ngrid), 
      col = col[1], add = TRUE)
image(x1, x2, matrix(ifelse(apply(x12, 1, c2) <= 0, 0, NA), ngrid, ngrid), 
      col = col[2], add = TRUE)
contour(x1, x2, matrix(apply(x12, 1, f), ngrid, ngrid), 
        nlevels = 21, add = TRUE)
points(GA@solution[1], GA@solution[2], col = "dodgerblue3", pch = 3)  # GA solution

## ------------------------------------------------------------------------
AQL   <- 0.01; alpha <- 0.05
LTPD  <- 0.06; beta  <- 0.10
plot(0, 0, type="n", xlim=c(0,0.2), ylim=c(0,1), bty="l", xaxs="i", yaxs="i", 
     ylab="Prob. of acceptance", xlab=expression(p))
lines(c(0,AQL), rep(1-alpha,2), lty=2, col="grey")
lines(rep(AQL,2), c(1-alpha,0), lty=2, col="grey")
lines(c(0,LTPD), rep(beta,2), lty=2, col="grey")
lines(rep(LTPD,2), c(beta,0), lty=2, col="grey")
points(c(AQL, LTPD), c(1-alpha, beta), pch=16)
text(AQL, 1-alpha, labels=expression(paste("(", AQL, ", ", 1-alpha, ")")), pos=4)
text(LTPD, beta, labels=expression(paste("(", LTPD, ", ", beta, ")")), pos=4)

## ------------------------------------------------------------------------
decode1 <- function(x)
{ 
  x <- gray2binary(x)
  n <- binary2decimal(x[1:l1])
  c <- min(n, binary2decimal(x[(l1+1):(l1+l2)]))
  out <- structure(c(n,c), names = c("n", "c"))
  return(out)
}

fitness1 <- function(x) 
{ 
  par <- decode1(x)
  n <- par[1]  # sample size
  c <- par[2]  # acceptance number
  Pa1 <- pbinom(c, n, AQL)
  Pa2 <- pbinom(c, n, LTPD)
  Loss <- (Pa1-(1-alpha))^2 + (Pa2-beta)^2
  -Loss
}

n  <- 2:200                  # range of values to search
b1 <- decimal2binary(max(n)) # max number of bits requires
l1 <- length(b1)             # length of bits needed for encoding
c  <- 0:20                   # range of values to search
b2 <- decimal2binary(max(c)) # max number of bits requires
l2 <- length(b2)             # length of bits needed for encoding

GA1 <- ga(type = "binary", fitness = fitness1, 
          nBits = l1+l2, 
          popSize = 100, maxiter = 1000, run = 100)
summary(GA1)
decode1(GA1@solution)

## ------------------------------------------------------------------------
plot(0,0,type="n", xlim=c(0,0.2), ylim=c(0,1), bty="l", xaxs="i", yaxs="i", 
     ylab=expression(P[a]), xlab=expression(p))
lines(c(0,AQL), rep(1-alpha,2), lty=2, col="grey")
lines(rep(AQL,2), c(1-alpha,0), lty=2, col="grey")
lines(c(0,LTPD), rep(beta,2), lty=2, col="grey")
lines(rep(LTPD,2), c(beta,0), lty=2, col="grey")
points(c(AQL, LTPD), c(1-alpha, beta), pch=16)
text(AQL, 1-alpha, labels=expression(paste("(", AQL, ", ", 1-alpha, ")")), pos=4)
text(LTPD, beta, labels=expression(paste("(", LTPD, ", ", beta, ")")), pos=4)
n <- 87; c <- 2
p <- seq(0, 0.2, by = 0.001)
Pa <- pbinom(2, 87, p)
lines(p, Pa, col = 2)

## ------------------------------------------------------------------------
decode2 <- function(x)
{ 
  n <- floor(x[1])         # sample size
  c <- min(n, floor(x[2])) # acceptance number
  out <- structure(c(n,c), names = c("n", "c"))
  return(out)
}

fitness2 <- function(x) 
{ 
  x <- decode2(x)
  n <- x[1]  # sample size
  c <- x[2]  # acceptance number
  Pa1 <- pbinom(c, n, AQL)
  Pa2 <- pbinom(c, n, LTPD)
  Loss <- (Pa1-(1-alpha))^2 + (Pa2-beta)^2
  return(-Loss)
}

GA2 <- ga(type = "real-valued", fitness = fitness2, 
          lower = c(2,0), upper = c(200,20)+1,
          popSize = 100, maxiter = 1000, run = 100)
summary(GA2)
t(apply(GA2@solution, 1, decode2))

## ---- eval=FALSE, echo=-(1:2)--------------------------------------------
#  set.seed(20181111)
#  options(digits = 4)
#  nrep <- 100
#  systime <- loss <- niter <- matrix(as.double(NA), nrow = nrep, ncol = 2,
#                                     dimnames = list(NULL, c("Binary", "Real-valued")))
#  for(i in 1:nrep)
#  {
#    t <- system.time(GA1 <- ga(type = "binary", fitness = fitness1,
#                               nBits = l1+l2, monitor = FALSE,
#                               popSize = 100, maxiter = 1000, run = 100))
#    systime[i,1] <- t[3]
#    loss[i,1]    <- -GA1@fitnessValue
#    niter[i,1]   <- GA1@iter
#    #
#    t <- system.time(GA2 <- ga(type = "real-valued", fitness = fitness2,
#                               lower = c(2,0), upper = c(200,20)+1,
#                               monitor = FALSE,
#                               popSize = 100, maxiter = 1000, run = 100))
#    systime[i,2] <- t[3]
#    loss[i,2]    <- -GA2@fitnessValue
#    niter[i,2]   <- GA2@iter
#  }
#  
#  describe <- function(x) c(Mean = mean(x), sd = sd(x), quantile(x))
#  
#  t(apply(systime, 2, describe))
#  #               Mean      sd    0%   25%    50%    75%  100%
#  # Binary      0.6902 0.20688 0.421 0.553 0.6340 0.7455 1.463
#  # Real-valued 0.3251 0.07551 0.252 0.275 0.2995 0.3470 0.665
#  
#  t(apply(loss, 2, describe))*1000
#  #                Mean     sd      0%     25%     50%     75%   100%
#  # Binary      0.09382 0.1919 0.05049 0.05049 0.05049 0.05049 1.5386
#  # Real-valued 0.09600 0.1551 0.05049 0.05049 0.05049 0.05049 0.6193
#  
#  t(apply(niter, 2, describe))
#  #              Mean    sd  0% 25%   50%   75% 100%
#  # Binary      160.8 48.31 100 129 146.0 172.2  337
#  # Real-valued 122.5 27.99 100 104 110.5 130.0  246

## ------------------------------------------------------------------------
GA <- ga(type = "real-valued", 
         fitness =  function(x) -Rastrigin(x[1], x[2]),
         lower = c(-5.12, -5.12), upper = c(5.12, 5.12), 
         popSize = 50, maxiter = 1000, run = 100,
         optim = TRUE)
summary(GA)
plot(GA)

## ---- eval=FALSE---------------------------------------------------------
#  library(GA)
#  fitness <- function(x)
#  {
#    Sys.sleep(0.01)
#    x*runif(1)
#  }
#  
#  library(rbenchmark)
#  out <- benchmark(
#    GA1 = ga(type = "real-valued",
#             fitness = fitness, lower = 0, upper = 1,
#             popSize = 50, maxiter = 100, monitor = FALSE,
#             seed = 12345),
#    GA2 = ga(type = "real-valued",
#             fitness = fitness, lower = 0, upper = 1,
#             popSize = 50, maxiter = 100, monitor = FALSE,
#             seed = 12345, parallel = TRUE),
#    GA3 = ga(type = "real-valued",
#             fitness = fitness, lower = 0, upper = 1,
#             popSize = 50, maxiter = 100, monitor = FALSE,
#             seed = 12345, parallel = 2),
#    GA4 = ga(type = "real-valued",
#             fitness = fitness, lower = 0, upper = 1,
#             popSize = 50, maxiter = 100, monitor = FALSE,
#             seed = 12345, parallel = "snow"),
#    columns = c("test", "replications", "elapsed", "relative"),
#    order = "test",
#    replications = 10)
#  out$average <- with(out, average <- elapsed/replications)
#  out[,c(1:3,5,4)]
#  ##   test replications elapsed average relative
#  ## 1  GA1           10 565.075 56.5075    3.975
#  ## 2  GA2           10 142.174 14.2174    1.000
#  ## 3  GA3           10 263.285 26.3285    1.852
#  ## 4  GA4           10 155.777 15.5777    1.096

## ---- eval=FALSE---------------------------------------------------------
#  library(doParallel)
#  workers <- rep(c("141.250.100.1", "141.250.105.3"), each = 8)
#  cl <- makeCluster(workers, type = "PSOCK")
#  registerDoParallel(cl)

## ---- eval=FALSE---------------------------------------------------------
#  clusterExport(cl, varlist = c("x", "fun"))
#  clusterCall(cl, library, package = "mclust", character.only = TRUE)

## ---- eval=FALSE---------------------------------------------------------
#  GA5 <- ga(type = "real-valued",
#            fitness = fitness, lower = 0, upper = 1,
#            popSize = 50, maxiter = 100, monitor = FALSE,
#            seed = 12345, parallel = cl)

## ---- eval=FALSE---------------------------------------------------------
#  stopCluster(cl)

## ---- echo=FALSE---------------------------------------------------------
# run not in parallel because it is not allowed in CRAN checks
GA <- gaisl(type = "real-valued", 
            fitness =  function(x) -Rastrigin(x[1], x[2]),
            lower = c(-5.12, -5.12), upper = c(5.12, 5.12), 
            popSize = 100, 
            maxiter = 1000, run = 100, 
            numIslands = 4, 
            migrationRate = 0.2, 
            migrationInterval = 50,
            parallel = FALSE)

## ---- eval=FALSE---------------------------------------------------------
#  GA <- gaisl(type = "real-valued",
#              fitness =  function(x) -Rastrigin(x[1], x[2]),
#              lower = c(-5.12, -5.12), upper = c(5.12, 5.12),
#              popSize = 100,
#              maxiter = 1000, run = 100,
#              numIslands = 4,
#              migrationRate = 0.2,
#              migrationInterval = 50)

## ------------------------------------------------------------------------
summary(GA)
plot(GA, log = "x")

## ---- eval = FALSE-------------------------------------------------------
#  data(fat, package = "UsingR")
#  mod <- lm(body.fat.siri ~ age + weight + height + neck + chest + abdomen +
#            hip + thigh + knee + ankle + bicep + forearm + wrist, data = fat)
#  summary(mod)
#  x <- model.matrix(mod)[,-1]
#  y <- model.response(mod$model)
#  
#  fitness <- function(string)
#  {
#    mod <- lm(y ~ x[,string==1])
#    -BIC(mod)
#  }
#  
#  library(memoise)
#  mfitness <- memoise(fitness)
#  is.memoised(fitness)
#  ## [1] FALSE
#  is.memoised(mfitness)
#  ## [1] TRUE
#  
#  library(rbenchmark)
#  tab <- benchmark(
#    GA1 = ga("binary", fitness = fitness, nBits = ncol(x),
#             popSize = 100, maxiter = 100, seed = 1, monitor = FALSE),
#    GA2 = ga("binary", fitness = mfitness, nBits = ncol(x),
#             popSize = 100, maxiter = 100, seed = 1, monitor = FALSE),
#    columns = c("test", "replications", "elapsed", "relative"),
#    replications = 10)
#  tab$average <- with(tab, elapsed/replications)
#  tab
#  ##   test replications elapsed relative average
#  ## 1  GA1           10  59.071    5.673  5.9071
#  ## 2  GA2           10  10.413    1.000  1.0413
#  
#  # to clear cache use
#  forget(mfitness)

## ------------------------------------------------------------------------
sessionInfo()

