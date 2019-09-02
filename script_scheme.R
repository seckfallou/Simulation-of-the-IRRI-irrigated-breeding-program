rm(list = ls())
library(BreedingSchemeLanguage)
library(magrittr)
library(lme4)
library(synbreed)
library(regress)
source("add_functions.R")
source("scheme_function.R")
source("burnin_phase.R")

#m <- Sys.getenv("SGE_TASK_ID")
m <- 1

set.seed(m)

# fixed parameters --------------------------------------------------------
effPopSize1 <- 30
effPopSize2 <- 300
nQTL <- 50#200
nChr <- 12
lengthChr <- 120
nMarkers <- 250#5000
propDomi <- 0
nEpiLoci <- 0#1
nInd <- 200#300
gVariance <- 2
plotTypeErrVars <- c(reduced = 8,
                     standard = 8)
locCorrelations <- NULL
fracGxEAdd <- 0.4


#out <- list()
res <- list()
a <- 0
Inter <- list(int_1 = c(2, 3), int_2 = c(4, 6))
param <- list(
  BURNIN=c(),
  EPISTASIS = c(),
  INTERACTION = Inter
  )

time1 <- Sys.time()
for (B in c(FALSE, TRUE)) {
  cat("--------------------------------------------------------------- \n")
  print(paste("Burn-in phase:", B))
  param[["BURNIN"]] <- append(param[["BURNIN"]], B)
  for (EPI in c(0, 0.02)) {
    cat("--------------------------------------------------------------- \n")
    print(paste("Number of epistatic locus:", EPI))
    if (length(param[["BURNIN"]]) == 2) {
      param[["EPISTASIS"]] <- append(param[["EPISTASIS"]], EPI)
    }
    if (exists("simEnv")) {
      rm(list = names(simEnv), envir = simEnv)
      rm(simEnv)
      print("Delete previous environment")
    }
    cat("------------------------------- \n")
  
    if (B == F) {
      print("Create new environment - effPopSize1")
      simEnv <- defineSpecies(
        nSim = 1,
        nCore = 1,
        saveDataFileName = paste("save" , m,  sep = ""),
        nChr = nChr,
        lengthChr = lengthChr,
        effPopSize = effPopSize1,
        nMarkers = nMarkers,
        nQTL = nQTL,
        propDomi = propDomi,
        nEpiLoci = EPI
      )
    }else{
      print("Create new environment - effPopSize2")
      simEnv <- defineSpecies(
        nSim = 1,
        nCore = 1,
        saveDataFileName = paste("save" , m,  sep = ""),
        nChr = nChr,
        lengthChr = lengthChr,
        effPopSize = effPopSize2,
        nMarkers = nMarkers,
        nQTL = nQTL,
        propDomi = propDomi,
        nEpiLoci = EPI
      )
    }
    
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
        res[[a]] <- scheme(burnin = B,
                           EBV = F,
                           interact = Inter[[INT]]
                          )
        for (n in names(res[[a]])) {
          res[[a]] <- add_param3(
            res[[a]],
            n,
            BURNIN = B,
            nEPI = EPI,
            INT = INT
          )
        }
    }
  }
  }
  if (length(param[["BURNIN"]]) == 2) {
    out <- list(parameters = param, results = res)
    rm(res, param)
  }

difftime(Sys.time(), time1)
pryr::mem_used()
readr::write_rds(out, paste("res", m, ".rds", sep = ""))
