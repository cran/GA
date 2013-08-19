##############################################################################
#                                                                            #
#                        GENETIC ALGORITHMS in R                             #
#                                                                            #
##############################################################################

ga <- function(type = c("binary", "real-valued", "permutation"), 
               fitness, ...,
               min, max, nBits,
               population = gaControl(type)$population,
               selection = gaControl(type)$selection,
               crossover = gaControl(type)$crossover, 
               mutation = gaControl(type)$mutation,
               popSize = 50, 
               pcrossover = 0.8, 
               pmutation = 0.1, 
               elitism = base::max(1, round(popSize*0.05)), 
               monitor = gaMonitor,
               maxiter = 100,
               run = maxiter,
               maxfitness = -Inf,
               names = NULL,
               suggestions, 
               keepBest = FALSE,
               parallel = FALSE,
               seed) 
{

  call <- match.call()
  
  type <- match.arg(type)
  if(!is.function(population)) population <- get(population)
  if(!is.function(selection))  selection  <- get(selection)
  if(!is.function(crossover))  crossover  <- get(crossover)
  if(!is.function(mutation))   mutation   <- get(mutation)
  
  if(missing(fitness))
    { stop("A fitness function must be provided") }
  if(!is.function(fitness)) 
    { stop("A fitness function must be provided") }
  if(popSize < 10) 
    { stop("The population size must be at least 10.") }
  if(maxiter < 1) 
    { stop("The maximum number of iterations must be at least 1.") }
  if(elitism > popSize) 
    { stop("The elitism must be larger that population size.") }
  if(pcrossover < 0 | pcrossover > 1)
    { stop("Probability of crossover must be between 0 and 1.") }
  if(pmutation < 0 | pmutation > 1)
    { stop("Probability of mutation must be between 0 and 1.") }
  if(missing(min) & missing(max) & missing(nBits))
    { stop("A min and max range of values (for 'real-valued' or 'permutation' GA) or nBits (for 'binary' GA) must be provided!") }
  
  switch(type, 
         "binary"      = { nBits <- as.vector(nBits)[1]
                           min <- max <- NA
                           nvars <- nBits 
                         },
         "real-valued" = { min <- as.vector(min)
                           max <- as.vector(max)
                           nBits <- NA
                           if(length(min) != length(max))
                             { stop("min and max must be vector of the same length!") }
                           nvars <- length(max) 
                         },
         "permutation" = { min <- as.vector(min)[1]
                           max <- as.vector(max)[1]
                           nBits <- NA
                           nvars <- length(seq(min,max)) 
                         }
        )

  if(missing(suggestions))
    { suggestions <- matrix(nrow = 0, ncol = nvars) }
  else
    { if(is.vector(suggestions)) 
        { if(nvars > 1) suggestions <- matrix(suggestions, nrow = 1)
          else          suggestions <- matrix(suggestions, ncol = 1) }
      else
        { suggestions <- as.matrix(suggestions) }
      if(nvars != ncol(suggestions))
        stop("Provided suggestions (ncol) matrix do not match number of variables of the problem!")
    }
  
  # Start parallel computing (if needed)
  parallel <- if(is.logical(parallel)) 
                { if(parallel) gaParallel(parallel) else FALSE }
              else { gaParallel(parallel) }
  
#   ## parallel computing #####################################################
#   # set default parallel functionality depending on system OS
#   parallelType <- if(.Platform$OS.type == "windows") 
#                     "snow" else "multicore"
#   if(is.character(parallel))
#     { parallelType <- parallel
#       parallel <- TRUE }
#   # start "parallel backend" if needed
#   if(parallel)
#     { # load package doParallel and all its dependencies
#       suppressPackageStartupMessages(availPkgs <- require("doParallel"))
#       if(!availPkgs)
#         stop("Required packages for parallel computing not available!")
#       if(is.numeric(parallel)) 
#         { # get number of cores from parallel if numeric
#           numCores <- as.integer(parallel) }
#       else 
#         { # get the current number of cores available
#           numCores <- detectCores() }
#       if(parallelType == "snow")
#         { # snow functionality on Unix-like systems & Windows
#           cl <- makeCluster(numCores)
#           # export environment
#           clusterExport(cl, varlist = ls(envir = parent.frame(), 
#                                          all.names = TRUE),
#                         envir = parent.frame()
#                         # envir = parent.env(environment())
#                        )
#           # load current packages in workers
#           pkgs <- .packages()
#           lapply(pkgs, function(pkg) 
#                  clusterCall(cl, library, package = pkg, 
#                              character.only = TRUE))
#           #
#           registerDoParallel(cl) 
#       }
#       else 
#         { # multicore functionality on Unix-like systems
#           registerDoParallel(cores = numCores) }
#     }
#   ###########################################################################

  bestEval <- rep(NA, maxiter)
  meanEval <- rep(NA, maxiter)
  bestSol <- if(keepBest) vector(mode = "list", length = maxiter)
             else         list()
  Fitness <- rep(NA, popSize)

  object <- new("ga", 
                call = call, 
                type = type,
                min = min, 
                max = max, 
                nBits = nBits, 
                names = if(is.null(names)) character() else names,
                popSize = popSize,
                iter = 0, 
                run = 1, 
                maxiter = maxiter,
                suggestions = suggestions,
                population = matrix(), 
                elitism = elitism, 
                pcrossover = pcrossover, 
                pmutation = pmutation,
                fitness = Fitness, 
                best = bestEval, 
                mean = meanEval,
                bestSol = bestSol)

  if(!missing(seed)) set.seed(seed)

  # generate beginning population
  Pop <- matrix(NA, nrow = popSize, ncol = nvars)
  ng <- min(nrow(suggestions), popSize)
  if(ng > 0) # use suggestion if provided
    { Pop[1:ng,] <- suggestions }
  # fill the rest with a random population
  if(popSize > ng)
    { Pop[(ng+1):popSize,] <- population(object)[1:(popSize-ng),] }
  object@population <- Pop

  # start iterations
  for(iter in 1:maxiter)
     {
      # evalute fitness function (if needed) 
      if(!parallel)
        { for(i in 1:popSize) 
             if(is.na(Fitness[i]))
               { Fitness[i] <- fitness(Pop[i,], ...) } 
        }
      else
        { Fitness <- foreach(i = 1:popSize, .combine = "c") %dopar% 
                     { if(is.na(Fitness[i])) fitness(Pop[i,], ...) 
                       else                  Fitness[i] }
        }
      bestEval[iter] <- max(Fitness, na.rm = TRUE)
      meanEval[iter] <- mean(Fitness, na.rm = TRUE)
      
      # update object
      object@iter <- iter
      object@population <- Pop
      object@fitness <- Fitness
      object@best <- bestEval
      object@mean <- meanEval
      
      if(keepBest) 
        object@bestSol[[iter]] <- unique(Pop[Fitness == bestEval[iter],,drop=FALSE])

      if(is.function(monitor)) 
        { monitor(object) }

      # check stopping criteria
      if(iter > 1)
        { if(bestEval[iter] > bestEval[iter-1]) 
               object@run <- 1 
          else 
               object@run <- object@run + 1 
        }
      if(object@run >= run) break  
      if(maxfitness > max(Fitness, na.rm = TRUE)) break
      if(object@iter == maxiter) break  

      # PopNew <- matrix(NA, nrow = popSize, ncol = nvars)
      # fitnessNew <- rep(NA, popSize)
      ord <- order(Fitness, decreasing = TRUE)
      PopSorted <- Pop[ord,,drop=FALSE]
      FitnessSorted <- Fitness[ord]
        
      # selection
      if(is.function(selection))
        { sel <- selection(object)
          Pop <- sel$population
          Fitness <- sel$fitness
        }
      else
        { sel <- sample(1:popSize, size = popSize, replace = TRUE)
          Pop <- object@population[sel,]
          Fitness <- object@fitness[sel]
        }
      object@population <- Pop
      object@fitness <- Fitness
    
      # crossover
      if(is.function(crossover) & pcrossover > 0)
        { nmating <- floor(popSize/2)
          mating <- matrix(sample(1:(2*nmating), size = (2*nmating)), ncol = 2)
          for(i in 1:nmating)
             { if(pcrossover > runif(1))
                 { parents <- mating[i,]
                   Crossover <- crossover(object, parents)
                   Pop[parents,] <- Crossover$children
                   Fitness[parents] <- Crossover$fitness
                 }
             }
          object@population <- Pop
          object@fitness <- Fitness
        }

      # mutation
      if(is.function(mutation) & pmutation > 0)
        { for(i in 1:popSize) 
             { if(pmutation > runif(1)) 
                 { Mutation <- mutation(object, i)
                   Pop[i,] <- Mutation
                   Fitness[i] <- NA
                 }
             }
          object@population <- Pop
          object@fitness <- Fitness
        }

      # elitism
      if(elitism > 0) 
        { ord <- order(object@fitness, na.last = TRUE)
          u <- which(!duplicated(PopSorted, margin = 1))
          Pop[ord[1:elitism],] <- PopSorted[u[1:elitism],]
          Fitness[ord[1:elitism]] <- FitnessSorted[u[1:elitism]]
          object@population <- Pop
          object@fitness <- Fitness
        } 
  
  }
      
  # in case of premature convergence remove NA from fitness evalutations
  object@best <- object@best[!is.na(object@best)]
  object@mean <- object@mean[!is.na(object@mean)]
  
  # get solution(s)
  object@fitnessValue <- max(object@fitness, na.rm = TRUE)
  valueAt <- which(object@fitness == object@fitnessValue)
  solution <- object@population[valueAt,,drop=FALSE]
  if(nrow(solution) > 1)
    { # find unique solutions to precision given by default tolerance
      eps <- gaControl("eps")
      solution <- unique(round(solution/eps)*eps, margin = 1)
    }
  colnames(solution) <- parNames(object)
  object@solution <- solution
  if(keepBest)
    object@bestSol <- object@bestSol[!sapply(object@bestSol, is.null)]  
  
  # return an object of class 'ga'
  return(object)
}

setClassUnion("numericOrNA", members = c("numeric", "logical"))

setClass(Class = "ga", 
         representation(call = "language",
                        type = "character",
                        min = "numericOrNA", 
                        max = "numericOrNA", 
                        nBits = "numericOrNA", 
                        names = "character",
                        popSize = "numeric",
                        iter = "numeric", 
                        run = "numeric", 
                        maxiter = "numeric",
                        suggestions = "matrix",
                        population = "matrix",
                        elitism = "numeric", 
                        pcrossover = "numeric", 
                        pmutation = "numeric",
                        fitness = "numericOrNA",
                        best = "numericOrNA", 
                        mean = "numericOrNA",
                        bestSol = "list",
                        fitnessValue = "numeric",
                        solution = "matrix"
                      ),
         package = "GA" 
) 

setMethod("print", "ga", function(x, ...) str(x))

setMethod("show", "ga",
function(object)
 { cat("An object of class \"ga\"\n")
   cat("\nCall:\n", deparse(object@call), "\n\n",sep="")
   cat("Available slots:\n")
   print(slotNames(object))
   # cat("  type:", object@type, "\n")
   # cat("  popSize:", object@popSize, "\n")
   # cat("  pcrossover:", object@pcrossover, "\n")
   # cat("  pmutation:", object@pmutation, "\n")
   # cat("  elitism:", object@elitism, "\n")
   # cat("  iter:", object@iter, "\n")
}) 


setMethod("summary", "ga", 
function(object, ...)
{

  nvars <- ncol(object@population)
  varnames <- parNames(object)
  domain <- NULL
  if(object@type == "real-valued")
    { domain <- rbind(object@min, object@max)
      rownames(domain) <- c("Min", "Max")
      if(ncol(domain) == nvars) 
         colnames(domain) <- varnames
    }
  suggestions <- NULL
  if(nrow(object@suggestions) > 0) 
    { suggestions <- object@suggestions
      dimnames(suggestions) <- list(1:nrow(suggestions), varnames) 
    }
  
  out <- list(type = object@type,
              popSize = object@popSize,
              maxiter = object@maxiter,
              elistism = object@elitism,
              pcrossover = object@pcrossover,
              pmutation = object@pmutation,
              domain = domain,
              suggestions = suggestions,
              iter = object@iter,
              fitness = object@fitnessValue,
              solution = object@solution)  
  class(out) <- "summary.ga"
  return(out)
})

print.summary.ga <- function(x, digits = getOption("digits"), ...)
{
  cat("+-----------------------------------+\n")
  cat("|         Genetic Algorithm         |\n")
  cat("+-----------------------------------+\n\n")
  cat("GA settings: \n")
  cat(paste("Type                  = ", x$type, "\n"))
  cat(paste("Population size       = ", x$popSize, "\n"))
  cat(paste("Number of generations = ", x$maxiter, "\n"))
  cat(paste("Elitism               = ", x$elitism, "\n"))
  cat(paste("Crossover probability = ", format(x$pcrossover, digits = digits), "\n"))
  cat(paste("Mutation probability  = ", format(x$pmutation, digits = digits), "\n"))
  
  if(x$type == "real-valued")
    { cat(paste("Search domain \n"))
      print(x$domain, digits = digits)
    }

  if(!is.null(x$suggestions))
    { cat(paste("Suggestions", "\n"))
      print(x$suggestions, digits = digits)
    }

  cat("\nGA results: \n")
  cat(paste("Iterations             =", format(x$iter, digits = digits), "\n"))
  cat(paste("Fitness function value =", format(x$fitness, digits = digits), "\n"))
  if(nrow(x$solution) > 1) 
    { cat(paste("Solutions              = \n")) }
  else
    { cat(paste("Solution               = \n")) }
  print(x$solution, digits = digits)

  invisible()
}


setMethod("plot", "ga", 
function(x, y, ylim, cex.points = 0.8, col = c("green3", "dodgerblue3"), pch = c(20, 17), lty = c(1,2), grid = graphics:::grid, ...)
{
  object <- x  # Argh.  Really want to use 'object' anyway
  is.final <- !(any(is.na(object@mean)))
  iters <- if(is.final) 1:object@iter else 1:object@maxiter
  plot(iters, object@best, type = "n",
       ylim = if(missing(ylim)) range(object@mean, object@best, na.rm = TRUE)
              else ylim, 
       xlab = "Generation", ylab = "Fitness value")
  if(is.final & is.function(grid)) 
    { grid() }
  points(iters, object@best,
         type = ifelse(is.final, "b", "p"),
         pch = pch[1], lty = lty[1], col = col[1], cex = cex.points)
  points(iters, object@mean,
         type = ifelse(is.final, "b", "p"),
         pch = pch[2], lty = lty[2], col = col[2], cex = cex.points)
  if(is.final)
    { legend("bottomright", legend = c("Best", "Mean"), 
             col = col, pch = pch, lty = lty, inset = 0.01) }
  else
    { title(paste("Iteration", object@iter), font.main = 1) }
    
  out <- list(iter = iters,  
              fitnessBest = object@best, 
              fitnessMean = object@mean)
  invisible(out)
})

# questa non funziona quando installa il pacchetto con NAMESPACE
setGeneric(name = "parNames", 
           def = function(object, ...) { standardGeneric("parNames") }
          )

setMethod("parNames", "ga",
function(object, ...)
{ 
  names <- object@names
  nvars <- ncol(object@population)
  if(length(names) == 0)
    { names <- paste("x", 1:nvars, sep = "") }
  return(names)
})
# per ora uso questo ma si dovrebbe ripristinare il metodo sopra:
# gaParNames <- function(object, ...)
# { 
#   names <- object@names
#   nvars <- ncol(object@population)
#   if(length(names) == 0)
#     { names <- paste("x", 1:nvars, sep = "") }
#   return(names)
# }

gaMonitor <- function(object, digits = getOption("digits"), ...)
{ 
  cat(paste("Iter =", object@iter, 
            " | Mean =", format(object@mean[object@iter], digits = digits), 
            " | Best =", format(object@best[object@iter], digits = digits), "\n"))
}
