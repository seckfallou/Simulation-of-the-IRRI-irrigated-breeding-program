# Data extraction in the environment  -------------------------------------
dataExtrac <- function(x) {
  simEnv[["sims"]][[1]][[x]]
}

# Showing last generated population ---------------------------------------
currentPop <- function() {
  max(dataExtrac("genoRec")$popID)
}

# Current population actually and size ------------------------------------
popSize <- function() {
  table(dataExtrac("genoRec")$popID)
}

# Base population ---------------------------------------------------------
basePop <- function() {
  table(dataExtrac("genoRec")$basePopID)
}

# Genotype data extraction for a given base population ID -----------------
genoData <- function(popID) {
  mrkPos <- simEnv[["sims"]][[1]]$mapData$markerPos#Marker popsition
  geno <- dataExtrac("geno")#haplotypes
  genoR <- dataExtrac("genoRec")
  if (any(popID == genoR$basePopID) == TRUE) {
    #Checking if it is subset pop
    nPop <- sum(genoR$basePopID == popID)#population size
    gid.diplo <- genoR$GID[genoR$basePopID == popID]
    score <-
      geno[gid.diplo * 2 - 1, mrkPos] + geno[gid.diplo * 2, mrkPos]
  } else{
    #pop is a subset from select() function
    nPop <- sum(genoR$popID == popID)#population size
    gid.diplo <- genoR$GID[genoR$popID == popID]
    score <-
      geno[gid.diplo * 2 - 1, mrkPos] + geno[gid.diplo * 2, mrkPos]
  }
  score <- score / 2 + 1 #Transform score to snp count
  colnames(score) <- paste("snp", 1:length(mrkPos), sep = ".")
  rownames(score) <- gid.diplo
  return(geno.Pop = score)
}

# Phenotypic data extraction  ---------------------------------------------
phenoData <- function(popID = NULL) {
  if (is.null(popID)) {
    stop("the population ID does not indicate")
  }
  pheno <-
    dataExtrac("phenoRec")#Record the phenotypic data in the dataframe pheno
  genoR <- dataExtrac("genoRec")
  if (any(popID == genoR$basePopID) == TRUE) {
    gid.pheno <-
      genoR$GID[genoR$basePopID == popID]#Extract the gid of the target population
  } else{
    #Phenotype for the corresponding subset population
    gid.pheno <- genoR$GID[genoR$popID == popID]
  }
  pheno.now <-
    pheno[pheno$phenoGID %in% gid.pheno,]#Extract the phenotypic data corresponding to the indicate population
  return(pheno.now)
}

# True Genetic value  -----------------------------------------------------
gvalData <- function(popID) {
  genoR <- dataExtrac("genoRec")
  gValue <-
    data.frame(GID = genoR$GID, gValue = dataExtrac("gValue"))
  if (any(popID == genoR$basePopID) == TRUE) {
    #Whole population
    gid.now <- genoR$GID[genoR$basePopID == popID]
  } else{
    #Subset population
    gid.now <- genoR$GID[genoR$popID == popID]
  }
  gVal.now <- gValue[gid.now,]
  return(genoVal.now = gVal.now)
}

# Allele frequencies calculation ------------------------------------------
allFreq <- function(geno) {
  n0 = apply(geno == 0, MARGIN = 2, FUN = sum) # total individuals not carrying the referrent allele
  n1 = apply(geno == 1, MARGIN = 2, FUN = sum)# total individuals carrying 1 copy of the referrent allele
  n2 = apply(geno == 2, MARGIN = 2, FUN = sum)# total individuals carrying 2 copies of the referrent allele
  N = n0 + n1 + n2 # equql to the sample size if there are not missign data
  p = (2 * n2 + n1) / (2 * N)
  q = 1 - p
  maf = pmin(p, q) #minor allele frequency (MAF) the second most frequent allele value
  return(list(
    maf = maf,
    p = p,
    q = q,
    sizePop = unique(N)
  ))
}

# #Add parameters to the burn-in output -----------------------------------
add_param <- function(list, obj, nQTL, nEPI, ERR, INT) {
  list[[obj]]$nQTL <- QTL
  list[[obj]]$nEPI <- EPI
  list[[obj]]$Error <- ERR
  list[[obj]]$Interaction <- INT
  return(list)
}

add_param2 <- function(list, obj, nEPI, INT) {
  list[[obj]]$nEPI <- EPI
  list[[obj]]$Interaction <- INT
  return(list)
}

add_param3 <- function(list, obj, BURNIN, nEPI, INT) {
  list[[obj]]$BURNIN <- B
  list[[obj]]$nEPI <- EPI
  list[[obj]]$Interaction <- INT
  return(list)
}

# add_ymean ---------------------------------------------------------------
# function to add a year mean to phenotypes (agronomic trend)
add_ymean <- function(sEnv = SimEnv,
                      year,
                      popID = NULL,
                      YrMn) {
  gRec <- with(simEnv, sims[[1]]$genoRec)
  pRec <- with(simEnv, sims[[1]]$phenoRec)
  if (is.null(popID)) {
    popID <- max(gRec$popID)
  }
  gidadd <- gRec[which(gRec$popID == popID), 'GID']
  ixmodify <- which(pRec$year == year & pRec$phenoGID %in% gidadd)
  pRec[ixmodify, 'pValue'] <- pRec[ixmodify, 'pValue'] + YrMn
  assign('precJR', value = pRec, envir = simEnv)
  with(simEnv, sims[[1]]$phenoRec <- precJR)
}

# modified predictValue function ------------------------------------------

predictValue_func <- function(bsl,
                              popID,
                              trainingPopID,
                              locations,
                              years,
                              sharingInfo) {
  if (is.null(locations))
    locations <- unique(bsl$phenoRec$loc)
  if (is.null(years))
    years <- unique(bsl$phenoRec$year)
  phenoRec <- subset(bsl$phenoRec,
                     subset = bsl$phenoRec$loc %in%
                       locations & bsl$phenoRec$year %in% years)
  if (sharingInfo == "none") {
    mt1ObsPerGID <- length(unique(phenoRec$phenoGID)) <
      nrow(phenoRec)
    if (mt1ObsPerGID) {
      fmla <- "pValue ~ (1|phenoGID)"
      if (length(unique(phenoRec$year)) > 1)
        fmla <- paste(fmla, "+ year + (1|phenoGID:year)")
      if (length(unique(phenoRec$loc)) > 1)
        fmla <- paste(fmla, "+ loc + (1|phenoGID:loc)")
      phenoRec$phenoGID <- factor(phenoRec$phenoGID)
      phenoRec$loc <- factor(phenoRec$loc)
      phenoRec$year <- factor(phenoRec$year)
      fitIID <- lmer(
        formula = as.formula(fmla),
        data = phenoRec,
        weights = 1 / phenoRec$error
      )
      predict <- ranef(fitIID)$phenoGID
      predGID <- as.numeric(rownames(predict))
    }
    else {
      predict <- phenoRec$pValue
      predGID <- as.numeric(phenoRec$phenoGID)
    }
  }
  if (sharingInfo == "parents") {
    #get predicted values of F1s based on parents
    gRec <- bsl$genoRec
    prRec <- bsl$predRec
    f1genoix <- which(gRec[, 'popID'] == popID)
    gRec <- gRec[f1genoix, ]
    p1blp <- prRec[match(gRec[, 2], prRec[, 1]), 3]
    p2blp <- prRec[match(gRec[, 3], prRec[, 1]), 3]
    predict <- (p1blp + p2blp) / 2 #predictions of the genos in gRec
    if (length(unique(predict)) < length(predict)) {
      predict[-match(unique(predict), predict)] <- -Inf
    }
    predGID <- gRec[, 1]
  }
  if (is.null(bsl$predRec)) {
    predNo <- 1
  }
  else {
    predNo <- max(bsl$predRec$predNo) + 1
  }
  toAdd <- data.frame(predGID, predNo, predict)
  colnames(toAdd) <- c("predGID", "predNo", "predict")
  bsl$predRec <- rbind(bsl$predRec, toAdd)
  bsl$selCriterion <- list(popID = popID, criterion = "pred")
  if (exists("totalCost", bsl))
    bsl$totalCost <- bsl$totalCost + bsl$costs$predCost
  return(bsl)
}


# Burnin Function ---------------------------------------------------------





# -------------------------------------------------------------------------

#' Calculate an additive relationship matrix
#'
#' \code{pedigreeToAmatrix} returns an additive relationship matrix from a
#'  pedigree specified in three columns. The first column has to be the row
#'  number and sire and dam columns refer directly to rows of the pedigree.
#'
#' \code{pedigreeToAmatrix} has some functionality useful for plants.  It can
#'  handle a pedigree of doubled haploid lines. Individuals can be
#'  self-pollinated for an arbitrary number of generations.
#'
#' @param pedColumns A data.frame with four columns. The first column
#'  has to be the row number and sire and dam columns refer directly to rows
#'  of the pedigree. Parents of founders need to be set to 0 (ZERO). The row
#'  of a child has to be after (i.e. a higher row number) that of its parents.
#'  If an individual has one known and one unknown parent, set the unknown
#'  parent to 0. The fourth column indicates: If negative, the individual is
#'  a DH: the F1 is created, generates a gamete, which is then doubled.
#'  If positive, the number of generations an individual was self-pollinated
#'  after it's F1 ancestor was created (can be 0 if the individual is the F1).
#'
#' @param aMatIn A square matrix that contains the additive relationship
#'  matrix between individuals at the beginning of the pedigree. If given,
#'  the function saves time by starting calculations after those individuals
#'  This aMatIn functionality is NOT compatible with calculating A inverse
#'
#' @return A matrix, \code{aMat}, the additive relationship matrix
#'
calcAmatrix <- function(pedColumns, aMatIn = NULL) {
  calcAmatRow <-
    function(pedRec) {
      # Function to process one row of pedColumns
      prog <- pedRec[1]
      sire <-
        max(pedRec[2], 0) # Founder population has negative parents
      dam <- max(pedRec[3], 0)
      progM1 <- prog - 1
      if (sire) {
        sireRow <- aMat[sire, 1:progM1]
      } else{
        sireRow <- integer(progM1)
      }
      if (dam) {
        damRow <- aMat[dam, 1:progM1]
      } else{
        damRow <- integer(progM1)
      }
      aMat[prog, 1:progM1] <<- (sireRow + damRow) / 2
      aMat[1:progM1, prog] <<- (sireRow + damRow) / 2
      if (pedRec[4] < 0) {
        aSelf <- 2
      } else{
        if (sire > 0 & dam > 0) {
          aSelf <- 1 + aMat[sire, dam] / 2
        } else{
          aSelf <- 1
        }
        if (pedRec[4] > 0) {
          # Number generations individual was selfed
          for (i in 1:pedRec[4])
            aSelf <- 1 + aSelf / 2
        }
      }
      aMat[prog, prog] <<- aSelf
    }#END calcAmatRow
  
  # calculate A here
  nInd <- nrow(pedColumns)
  aMat <- matrix(0, nInd, nInd)
  if (is.null(aMatIn)) {
    # the very first individual in the pedigree has to be a founder
    initSelf <- 1
    if (pedColumns[1, 4] < 0)
      initSelf <- 2
    if (pedColumns[1, 4] > 0)
      initSelf <- 2 - 2 ^ (-pedColumns[1, 4])
    aMatIn <- as.matrix(initSelf)
  }
  nIn <- nrow(aMatIn)
  aMat[1:nIn, 1:nIn] <- aMatIn
  start <- nIn + 1
  apply(pedColumns[start:nInd, ], 1, calcAmatRow)
  rownames(aMat) <- colnames(aMat) <- pedColumns[, 1]
  return(aMat)
}
