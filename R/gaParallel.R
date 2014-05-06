gaParallel <- function(parallel = TRUE, ...)
{
# Start parallel computing on GA package
  
  # load package doParallel and all its dependencies (i.e., foreach,
  # iterators)
  suppressPackageStartupMessages(availPkgs <- require("doParallel"))
  if(!availPkgs)
    stop("Required packages (doParallel, foreach and iterators) for parallel computation not available!")

  # set default parallel functionality depending on system OS:
  # - snow functionality on Windows OS
  # - multicore functionality on Unix-like systems (Unix/Linux & Mac OSX)
  parallelType <- if(.Platform$OS.type == "windows") 
                    "snow" else "multicore"

  # get the current number of cores available
  numCores <- detectCores()

  # set parameters for parallelization
  if(is.logical(parallel))
    { NULL }
  else if(is.numeric(parallel))
    { numCores <- as.integer(parallel)
      parallel <- TRUE }
  else if(is.character(parallel))
    { parallelType <- parallel
      parallel <- TRUE 
    }
  else parallel <- FALSE
  
  attr(parallel, "type") <- parallelType
  attr(parallel, "cores") <- numCores

  # start "parallel backend" if needed
  if(parallel)
  { 
    if(parallelType == "snow")
      { 
        # snow functionality on Unix-like systems & Windows
        cl <- makeCluster(numCores, type = "PSOCK")
        attr(parallel, "cluster") <- cl
        # export parent environment
        varlist <- ls(envir = parent.frame(), all.names = TRUE)
        varlist <- varlist[varlist != "..."]
        clusterExport(cl, varlist = varlist,
                          envir = parent.frame()
                          # envir = parent.env(environment())
                     )
        # export global environment (workspace)
        clusterExport(cl, varlist = ls(envir = globalenv(), all.names = TRUE),
                          envir = globalenv())
        # load current packages in workers
        pkgs <- .packages()
        lapply(pkgs, function(pkg) 
               clusterCall(cl, library, package = pkg, character.only = TRUE))
        #
        registerDoParallel(cl, cores = numCores)
      }
      else if(parallelType == "multicore")
        { # multicore functionality on Unix-like systems
          cl <- makeCluster(numCores, type = "FORK")
          registerDoParallel(cl, cores = numCores) 
          attr(parallel, "cluster") <- cl
        }
      else 
        { stop("Only 'snow' and 'multicore' clusters allowed!") }
  }

  return(parallel)
}