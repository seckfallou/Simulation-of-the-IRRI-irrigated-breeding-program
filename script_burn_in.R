# load packages -----------------------------------------------------------
rm(list = ls())
library(BreedingSchemeLanguage)
library(lme4)
library(magrittr)

source('Extraction_function.R')
source('burn_in.R')


m <- Sys.getenv("SGE_TASK_ID")
#m <- 1

set.seed(m)

# defineSpecies parameters ------------------------------------------------
nMrkr <- 5000
nInd <- 300 #number of founders
sharing <-  "none"
nCore <- 1
nSim <- 1
nChr <- 12
lengthChr <- 120
effPopSize <- 300
propDomi <- 0
#nEpiLoci <- 0

#out <- list()
res <- list()
a <- 0
fracGxEAdd = 0.4
errVar <- list(err_1 = c(8, 8, 8), err_2 = c(16, 16, 16))
Inter <- list(int_1 = c(6, 8), int_2 = c(24, 32))
param <- list(
  QTL = c(),
  EPISTASIS = c(),
  ERROR = errVar,
  INTERACTION = Inter
)
cat("--------------------------------------------------------------- \n")
print(errVar)
print(Inter)

time1 <- Sys.time()
for (QTL in c(250, 500)) {
  cat("--------------------------------------------------------------- \n")
  print(paste("Number of QTL:", QTL))
  param[["QTL"]] <- append(param[["QTL"]], QTL)
  for (EPI in c(0, 2)) {
    cat("--------------------------------------------------------------- \n")
    print(paste("Number of epistatic locus:", EPI))
    if (length(param[["QTL"]]) == 2) {
      param[["EPISTASIS"]] <- append(param[["EPISTASIS"]], EPI)
    }
    if (exists("simEnv")) {
      rm(list = names(simEnv), envir = simEnv)
      rm(simEnv)
      print("Delete previous environment")
    }
    cat("------------------------------- \n")
    print("Create new environment")
    simEnv <- defineSpecies(
      loadData = NULL,
      nSim = nSim,
      nCore = nCore,
      saveDataFileName = paste("save" , m,  sep = ""),
      nChr = nChr,
      lengthChr = lengthChr,
      effPopSize = effPopSize,
      nMarkers = nMrkr,
      nQTL = QTL,
      propDomi = propDomi,
      nEpiLoci = EPI
    )
    for (ERR in names(errVar)) {
      cat("------------------------------- \n")
      print(paste("Level of error variance:", ERR))
      for (INT in names(Inter)) {
        cat("------------------------------- \n")
        print(paste("Level of interaction:", INT))
        if (!is.null(simEnv$sims[[1]]$genoRec)) {
          rm(list = names(simEnv), envir = simEnv)
          rm(simEnv)
          print("Load previous environment")
          cat("\n")
          simEnv <- defineSpecies(loadData = paste("save" , m,  sep = ""))
        }
        a <- a + 1
        res[[a]] <- burnin(
          ErrVars = errVar[[ERR]],
          interact = Inter[[INT]],
          fracGxEAdd = 0.4,
          locCorr = NULL,
          pop = FALSE,
          verbose = TRUE
        )
        
        for (n in names(res[[a]])) {
          res[[a]] <- add_param(
            res[[a]],
            n,
            nQTL = QTL,
            nEPI = EPI,
            ERR = ERR,
            INT = INT
          )
        }
      }
    }
  }
  
  if (length(param[["QTL"]]) == 2) {
    out <- list(parameters = param, results = res)
    rm(res, param)
  }
}

difftime(Sys.time(), time1)
pryr::mem_used()
readr::write_rds(out, paste("res", m, ".rds", sep = ""))
