burnin <- function(ErrVars = c(4, 3, 1),
                   interact = c(1.25, 1.75),
                   fracGxEAdd = 0.4,
                   locCorr = NULL,
                   pop = TRUE,
                   verbose = TRUE,
                   graph = FALSE,
                   pdf_save = FALSE) {
  # defineVariances parameters ----------------------------------------------
  gV <- 2#1
  defineVariances(
    gVariance = gV,
    gByLocVar = interact[1],
    gByYearVar = interact[2],
    locCorrelations = locCorr,
    plotTypeErrVars = c(
      Headrow = ErrVars[1],
      Preliminary = ErrVars[2],
      Advanced = ErrVars[3]
    ),
    fracGxEAdd = fracGxEAdd #nLoci-->nMarkers + nQTL * (nEpiLoci + 1) * 2
  )
  
  # Breeding parameters -----------------------------------------------------
  # Number of Locations - Plote type and number of Replication --------------
  loc1 <- 1
  loc2 <- 1
  loc3 <- 1:6
  # number of year to perform the breeding (taken as season)
  yr <- 1:7
  pTyp1 <- "Headrow"
  pTyp2 <- "Preliminary"
  pTyp3 <- "Advanced"
  rep1 <- 1
  rep2 <- 1
  rep3 <- 2#4
  
  # Length of a breeding cycle ----------------------------------------------
  Yr1 <- 1:8#fix the founders
  lSSD <- 5 #SSD
  nCycl <- 6#1
  
  # Choice parents and number of lines to make ------------------------------
  nF1 <- 100 # 100 best F1 among the nCros realized
  nLT <- 50 * nF1 # Number of lines in fisrt testing stage
  
  # Selection intensity in each testing stage -------------------------------
  #Line Testing--------> 800 selected for PYT
  nPYT <- nLT * 0.16
  #Preliminary yield trial------> 240 selected for AYT
  nAYT <- nPYT * 0.3
  #Advanced yield trial-------> 80 selected
  nparent <- 100
  
  
  # First simulation --------------------------------------------------------
  
  # Initiation de la selection a partir des haplotypes ancestraux simules.
  # On suppose demarrer sur une population de base de 300 individus
  ## initialize and selfing founders to get lines
  if (verbose == T) {
    print("Create founders")
  }
  initializePopulation(nInd = nInd)
  if (verbose == T) {
    print('number of monomorphic SNPs')
    print(length(which(allFreq(
      geno = genoData(popID = currentPop(simEnv))
    )$maf == 0)))
  }
  
  pop_select <- currentPop(simEnv)
  
  Cycl <- "Cycle_0"
  
  if (verbose == T) {
    print("Lines fixation via selfing")
  }
  for (i in Yr1) {
    selfFertilize(nProgeny = nInd)
  }
  pop_select <- append(pop_select, currentPop(simEnv))
  if (verbose == T) {
    print('number of monomorphic SNPs')
    print(length(which(allFreq(
      geno = genoData(popID = currentPop(simEnv))
    )$maf == 0)))
  }
  # Create objects to store information (maf, tbv, heritability)
  MAF <- matrix(NA, nrow = (nCycl + 1), ncol = nMrkr)
  MAF[1,] <- allFreq(genoData(popID = 0))$maf
  #MAF[2,] <- allFreq(genoData(popID = currentPop(simEnv)))$maf
  T_BV <- data.frame(matrix(NA, nrow = nCycl, ncol = 2))
  colnames(T_BV) <- c("Cycle", "T_BV")
  T_BV[1, 1] <- Cycl
  T_BV[1, 2] <- mean(gvalData(0)$gValue)
  T_BV[1, 2] <- T_BV[1, 2] - mean(gvalData(0)$gValue)
  
  EBV <- data.frame(matrix(NA, nrow = nCycl, ncol = 2))
  colnames(EBV) <- c("Cycle", "ebv")
  
  # Correlation of phenotypic and genetic values ----------------------------
  Corr <-
    matrix(NA,
           nrow = 4,
           ncol = nCycl,
           dimnames = list(c("Cycle", "LT", "PYT", "AYT")))
  
  Corr <- as.data.frame(t(Corr))
  
  # Genetic variance --------------------------------------------------------
  gvar <-
    matrix(NA,
           nrow = 4,
           ncol = nCycl,
           dimnames = list(c("Cycle", "LT", "PYT", "AYT")))
  
  gvar <- as.data.frame(t(gvar))
  
  # Heritabily for each trial stage over the breeding cycle-------------
  h2 <-
    matrix(NA,
           nrow = 4,
           ncol = nCycl,
           dimnames = list(c("Cycle", "LT", "PYT", "AYT")))
  h2 <- as.data.frame(t(h2))
  
  
  
  H2 <-
    matrix(NA,
           nrow = 3,
           ncol = nCycl,
           dimnames = list(c("Cycle", "PYT", "AYT")))
  H2 <- as.data.frame(t(H2))
  
  
  
  if (pdf_save == T) {
    pdf(file = "Correlation of phenotypic and genetic values.pdf")
  }
  par(mfrow = c(nCycl, 3), mar = c(2, 4, 2, 2))
  icp <- NULL
  era_pop <- NULL
  # Start cycle -------------------------------------------------------------
  
  for (i in 1:nCycl) {
    if (verbose == T) {
      cat('\n')
      print(paste('Cycle: ', i, sep = ''))
      cat('\n')
      print('number of monomorphic SNPs')
      print(length(which(allFreq(
        geno = genoData(popID = currentPop(simEnv))
      )$maf == 0)))
    }
    #Parental selection and crossing
    if (i == 1) {
      #Until phenotypic data are available the 50 parents are selected randomly
      select(nSelect = nparent, random = T)
    } else{
      select(popID = era_pop[i - 1] - 1,
             nSelect = nparent,
             random = F)
    }
    if (pop == T) {
      print(popSize(simEnv))
    }
    cross(
      popID = currentPop(simEnv),
      nProgeny = nF1,
      equalContribution = F,
      notWithinFam = T
    )
    if (pop == T) {
      print(popSize(simEnv))
    }
    if (verbose == T) {
      print('SSD')
    }
    # Line fixation rapid generation (SSD) ------------------------------------
    self <- 1:lSSD
    for (s in self) {
      selfFertilize(nProgeny = nLT)
      if (pop == T) {
        print(popSize(simEnv))
      }
    }
    # 1st testing and advancement/create F7 by Selfing F6 ---------------------
    phenotype(
      popID = currentPop(simEnv),
      locations = loc1,
      nRep = rep1,
      plotType = pTyp1,
      years = yr[1]
    )
    
    gval <- gvalData(currentPop(simEnv))$gValue
    pval <- phenoData(currentPop(simEnv))$pValue
    Corr[i, 2] <- cor(gval, pval)
    gvar[i, 2] <- var(gval)
    
    if (verbose == T) {
      print("phenotyping - headrow")
      print('Correlation = ')
      print(cor(gval, pval))
      print('h² =')
      print(var(gval) / var(pval))
    }
    if (graph == T) {
      plot(
        gval,
        pval,
        xlab = 'genetic value',
        ylab = 'phenotypic value',
        main = c('Headrow \n' ,
                 paste("Corr = " , round(
                   cor(gval, pval$pval), 2
                 ))),
        cex.main = 0.7
      )
    }
    
    #Heritability
    h2[i, 2] <-
      var(gval) / var(pval)
    #Advancement
    select(nSelect = nPYT)
    selfFertilize(nProgeny = nPYT)#F7 generation
    if (pop == T) {
      print(popSize(simEnv))
    }
    
    # 2nd testing and advancement/create F8 by Selfing F7 ---------------------
    phenotype(
      popID = currentPop(simEnv),
      locations = loc2,
      nRep = rep2,
      plotType = pTyp2,
      years = yr[2:3]
    )
    
    gval <- gvalData(currentPop(simEnv))$gValue
    pheno <- phenoData(currentPop(simEnv))
    pval <- pheno %>%
      dplyr::group_by((phenoGID)) %>%
      dplyr::summarise(pval = mean(pValue))
    Corr[i, 3] <- cor(gval, pval$pval)
    gvar[i, 3] <- var(gval)
    
    if (verbose == T) {
      print("phenotyping - OYT")
      print('Correlation = ')
      print(cor(gval, pval$pval))
      print('h² =')
      print(var(gval) / var(pval$pval))
    }
    if (graph == T) {
      plot(
        gval,
        pval$pval,
        xlab = 'genetic value',
        ylab = 'phenotypic value',
        main = c('OYT \n' ,
                 paste("Corr = " , round(
                   cor(gval, pval$pval), 2
                 ))),
        cex.main = 0.7
      )
    }
    #Heritability of the PYT stage
    pheno$phenoGID <- as.factor(pheno$phenoGID)
    pheno$year <- as.factor(pheno$year)
    
    mod_pyt <- lmer(pValue ~ 1 + year +
                      (1 | phenoGID),
                    data = pheno,
                    REML = TRUE)
    VarComp_PYT <- as.data.frame(VarCorr(mod_pyt))[, -c(2, 3)]
    var_tot <-
      VarComp_PYT[VarComp_PYT$grp == "phenoGID", 'vcov'] +
      VarComp_PYT[VarComp_PYT$grp == "Residual", 'vcov'] / 2
    H2[i, 2] <-
      VarComp_PYT[VarComp_PYT$grp == "phenoGID", 'vcov'] / var_tot
    h2[i, 3] <- var(gval) / var_tot
    
    select(nSelect = nAYT)
    selfFertilize(nProgeny = nAYT)#F8 generation
    if (pop == T) {
      print(popSize(simEnv))
    }
    pop_select <-  append(pop_select, currentPop(simEnv))
    Cycl <- append(Cycl, paste("Cycle", i, sep = "_"))
    
    # 3th testing and Last lines selected the F9 ------------------------------
    phenotype(
      popID = currentPop(simEnv),
      locations = loc3,
      nRep = rep3,
      plotType = pTyp3,
      years = yr[4:5]
    )
    
    gval <- gvalData(currentPop(simEnv))$gValue
    pheno <- phenoData(pop_select[i+2])
    pval <- pheno %>%
      dplyr::group_by((phenoGID)) %>%
      dplyr::summarise(pval = mean(pValue))
    Corr[i, 4] <- cor(gval, pval$pval)
    gvar[i, 4] <- var(gval)
    if (verbose == T) {
      print("phenotyping - AYT")
      print('Correlation = ')
      print(cor(gval, pval$pval))
      print('h² =')
      print(var(gval) / var(pval$pval))
      
    }
    
    if (graph == T) {
      plot(
        gval,
        pval$pval,
        xlab = 'genetic value',
        ylab = 'phenotypic value',
        main = c('AYT \n' ,
                 paste("Corr = " , round(
                   cor(gval, pval$pval), 2
                 ))),
        cex.main = 0.7
      )
    }
    
    #Heritability of the AYT stage
    pheno$phenoGID <- as.factor(pheno$phenoGID)
    pheno$loc <- as.factor(pheno$loc)
    pheno$year <- as.factor(pheno$year)
    mod_ayt <-
      lmer(
        pValue ~ 1 +
          year +
          loc +
          (1 | phenoGID) +
          (1 | phenoGID:year) +
          (1 | phenoGID:loc),
        data = pheno,
        REML = TRUE
      )
    VarComp_AYT <- as.data.frame(VarCorr(mod_ayt))[, -c(2, 3)]
    var_tot <-
      VarComp_AYT[VarComp_AYT$grp == "phenoGID", 'vcov'] +
      VarComp_AYT[VarComp_AYT$grp == "phenoGID:loc", 'vcov'] / length(loc3) +
      VarComp_AYT[VarComp_AYT$grp == "phenoGID:year", 'vcov'] / 2 +
      VarComp_AYT[VarComp_AYT$grp == "Residual", 'vcov'] / (length(loc3) * 2)
    H2[i, 3] <-
      VarComp_AYT[VarComp_AYT$grp == "phenoGID", 'vcov'] / var_tot
    h2[i, 4] <-
      var(gval) / var_tot
    EBV[i , 1] <- i
    EBV[i , 2] <- mean(ranef(mod_ayt)$phenoGID[[1]])
    
    # -------------------------------------------------------------------------
    
    
    MAF[i + 1,] <-
      allFreq(genoData(popID = currentPop(simEnv)))$maf
    T_BV[i + 1, 1] <- Cycl[i + 1]#currentPop(simEnv)
    
    if (mean(gvalData(0)$gValue) < 0) {
      T_BV[i + 1, 2] <- mean(gvalData(currentPop(simEnv))$gValue)
    } else{
      T_BV[i + 1, 2] <-
        mean(gvalData(currentPop(simEnv))$gValue) - mean(gvalData(0)$gValue)
    }
    
    select(nSelect = 2 , random = F)
    era_pop <- append(era_pop, currentPop(simEnv))
    
    if (i == nCycl) {
      #The last evaluated population taken as Breeding germplasm
      icp <- currentPop(simEnv) - 1
      #-------------------------------------------------------------------
      gVal_PYT <- gvalData(icp - 2)[, 2]
      #corresponding to the PYT level
      pheno_PYT <- phenoData(icp - 2)
      pheno_PYT$phenoGID <- as.factor(pheno_PYT$phenoGID)
      pheno_PYT$year <- as.factor(pheno_PYT$year)
      
      mod_pyt <- lmer(pValue ~ 1 + year +
                        (1 | phenoGID),
                      data = pheno_PYT,
                      REML = TRUE)
      
      
      Blup_PYT <- ranef(mod_pyt)$phenoGID[[1]]
      gVal_PYT <- data.frame(T_BV = gVal_PYT, Blup = Blup_PYT)
      
      #-------------------------------------------------------------------
      
      gVal_AYT <- gvalData(icp)[, 2]
      pheno_AYT <- phenoData(icp) #corresponding to the AYT level
      pheno_AYT$phenoGID <- as.factor(pheno_AYT$phenoGID)
      pheno_AYT$loc <- as.factor(pheno_AYT$loc)
      pheno_AYT$year <- as.factor(pheno_AYT$year)
      mod_ayt <-
        lmer(
          pValue ~ 1 +
            year +
            loc +
            (1 | phenoGID) +
            (1 | phenoGID:year) +
            (1 | phenoGID:loc),
          data = pheno_AYT,
          REML = TRUE
        )
      
      
      Blup_AYT <- ranef(mod_ayt)$phenoGID[[1]]
      gVal_AYT <- data.frame(T_BV = gVal_AYT, Blup = Blup_AYT)
      
      AYT_loc <- phenoData(icp) %>%
        tidyr::spread(loc, pValue, sep = ".") %>%
        dplyr::group_by(phenoGID) %>%
        dplyr::summarise(
          mloc1 = mean(loc.1),
          mloc2 = mean(loc.2),
          mloc3 = mean(loc.3),
          mloc4 = mean(loc.4),
          mloc5 = mean(loc.5),
          mloc6 = mean(loc.6)
        )
      corr_loc <- as.data.frame(cor(AYT_loc[, -1]))
      
      AYT_year <- phenoData(icp) %>%
        tidyr::spread(year, pValue, sep = ".") %>%
        dplyr::group_by(phenoGID) %>%
        dplyr::summarise(mYr4 = mean(year.4),
                         mYr5 = mean(year.5))
      
      corr_year <- as.data.frame(cor(AYT_year[, -1]))
      
    }
  }
  
  #ERA
  phenotype(
    popID = era_pop,
    locations = loc3,
    nRep = rep3,
    plotType = pTyp3,
    years = yr[6:7]
  )
  era_pheno <- phenoData(era_pop)
  era_pheno <- era_pheno[era_pheno$year %in% c(6, 7),]
  Cycle <- rep(rep(rep(1:nCycl, each = 2), nCycl), 2)
  era_pheno <- cbind(Cycle, era_pheno) %>%
    tibble::as_tibble()
  ERA <- era_pheno %>% dplyr::group_by(Cycle) %>%
    dplyr::summarise(pheno = mean(pValue))
  ERA$Cycle %<>%  as.numeric()
  era_gg <- lm(pheno ~ Cycle , data = ERA)
  GG <- coefficients(era_gg)[[2]]
  ERA$GGain <- GG
  ERA$PercentGain <- round(GG / ERA$pheno[1], 2) * 100
  #EBV
  EBV <- EBV %>%
    tibble::as_tibble()
  EBV$Cycle %<>%  as.numeric()
  ebv_gg <- lm(ebv ~ Cycle , data = EBV)
  GG <- coefficients(ebv_gg)[[2]]
  EBV$GGain <- GG
  EBV$PercentGain <- round(GG / EBV$ebv[1], 2) * 100
  
  
  if (pdf_save == T) {
    dev.off()
  }
  rownames(MAF) <- Cycl#Whole cycles
  Corr[, 1] <- Cycl[-1]
  gvar[, 1] <- Cycl[-1]
  h2[, 1] <- Cycl[-1]
  H2[, 1] <- Cycl[-1]
  MAF <- as.data.frame(t(MAF))
  return(
    list(
      MAF = MAF,
      T_BV = T_BV,
      gVal_PYT_f = gVal_PYT,
      gVal_AYT_f = gVal_AYT,
      Corr = Corr,
      gvar = gvar,
      h2 = h2,
      H2 = H2,
      corr_loc = corr_loc,
      corr_year = corr_year,
      EBV = EBV,
      ERA = ERA
    )
  )
}
