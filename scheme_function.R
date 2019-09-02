scheme <- function(burnin = FALSE,
                   interact = c(16, 12),
                   print_pop = FALSE,
                   verbose = TRUE,
                   graph = FALSE,
                   pdf_save = FALSE,
                   sharingInfo = "none",
                   blup_lst = FALSE,
                   EBV = FALSE) {
  # Scheme parameters
  nCycl <- 6# burn-in scheme
  cyc_dur <- 10 # seasons (2 per years)
  nyears <- 2 * cyc_dur #6 * cyc_dur
  nparents <- 15#50
  nCross <- (nparents * (nparents - 1)) / 2 # number of crossing
  nF1 <- 100#
  nRGA <- 4
  ntest1 <- nF1 * 10#80
  ntest2 <- 200#2000
  ntest3 <- 50
  within_fam <- FALSE
  max_per_fam2 <- 1
  max_per_fam3 <- 1
  #Non genetic trend
  agroTrend <- FALSE
  agro_trend <- 0.2
  # Parent are selected from the current cycle only
  yearbyyear <- TRUE
  # Parent are selected by "family"
  sel_family <- FALSE
  # can the parent be reuse if yearbyyear is false
  reuse_parent <- TRUE
  loc1 <- 1
  loc2 <- 1
  loc3 <- 1:6
  rep1 <- 1
  rep2 <- 1
  rep3 <- 2
  plottype1 <- "reduced"
  plottype2 <- "standard"
  plottype3 <- "standard"
  
  #Store information from burn-in scheme
 genoFounders <- as.data.frame(matrix(NA, nrow = nInd, ncol = nMarkers))
  B_MAF <- matrix(NA, nrow = (nCycl + 1), ncol = nMarkers)
  B_MAF <- as.data.frame(t(B_MAF))
  B_gvar <-
    matrix(NA,
           nrow = 2,
           ncol = nCycl,
           dimnames = list(c("Year", "gVar")))
  B_gvar <- as.data.frame(t(B_gvar))
  
  if (burnin == T) {
    if (verbose == T) {
      cat("------------------------------------------- \n")
      print("Start of the burn-in phase")
    }
    ## ----defineVariances-----------------------------------------------------
    defineVariances(
      gVariance = gVariance,
      locCorrelations = locCorrelations,
      gByLocVar = interact[1],
      gByYearVar = interact[2],
      fracGxEAdd = fracGxEAdd,
      plotTypeErrVars = plotTypeErrVars
    )
    # Length of a breeding cycle ----------------------------------------------
    nCycl <- 6#1
    # Select parents
    #nF1 <- 100
    nLT <- 50 * nF1
    nPYT <- 1000
    nAYT <- 200
    #nparent <- 30
    yr <- 1:5
    ## initialize and selfing founders to get lines
    
    if (verbose == T) {
      print("Create founders")
    }
    initializePopulation(nInd = nInd)
    if (verbose == T) {
      print("Lines fixation via selfing")
    }
    for (s in 1:8) {
      selfFertilize(nProgeny = nInd)
    }
    
    B_MAF <- matrix(NA, nrow = (nCycl + 1), ncol = nMarkers)
    B_MAF[1, ] <- allFreq(genoData(popID = 0))$maf
   
    for (t in 1:nCycl) {
      #Start cycle
      if (verbose == T) {
        print(paste('Cycle :', t , sep = " "))
      }
      #Parental selection and crossing
      if (t == 1) {
        #Until phenotypic data are available the 50 parents are selected randomly
        select(nSelect = nparents, random = T)
      } else{
        select(popID = currentPop(),
               nSelect = nparents,
               random = F)
      }
      if (print_pop == T) {
        print(popSize())
      }
      cross(
        popID = currentPop(),
        nProgeny = nF1,
        equalContribution = F,
        notWithinFam = T
      )
      if (print_pop == T) {
        print(popSize())
      }
      if (verbose == T) {
        print('SSD')
      }
      # Line fixation rapid generation (SSD) ------------------------------------
      for (s in 1:5) {
        selfFertilize(nProgeny = nLT)
        if (print_pop == T) {
          print(popSize())
        }
      }
      # 1st testing and advancement/create F7 by Selfing F6 ---------------------
      if (verbose == T) {
        print('LST')
      }
      phenotype(
        popID = currentPop(),
        locations = loc1,
        nRep = rep1,
        plotType = plottype1,
        years = yr[1]
      )
      
      # 2nd testing and advancement/create F8 by Selfing F7 ---------------------
      if (verbose == T) {
        print('PYT')
      }
      #Advancement
      select(nSelect = nPYT)
      selfFertilize(nProgeny = nPYT)#F7 generation
      
      phenotype(
        popID = currentPop(),
        locations = loc2,
        nRep = rep2,
        plotType = plottype2,
        years = yr[2:3]
      )
      
      select(nSelect = nAYT)
      selfFertilize(nProgeny = nAYT)#F8 generation
      
      # 3th testing and Last lines selected the F9 ------------------------------
      if (verbose == T) {
        print('AYT')
      }
      phenotype(
        popID = currentPop(),
        locations = loc3,
        nRep = rep3,
        plotType = plottype3,
        years = yr[4:5]
      )
      
      
      B_MAF[t + 1,] <-
        allFreq(genoData(popID = currentPop()))$maf
      
      gval <- gvalData(currentPop())$gValue
      B_gvar[t, 1] <- t
      B_gvar[t, 2] <- var(gval)
      
      
      if (t == nCycl) {
        B_MAF <- as.data.frame(t(B_MAF))
        colnames(B_MAF) <- 0:nCycl
        genoFounders <- genoData(popID = currentPop())
        
        #The last evaluated population taken as Breeding germplasm
        if (verbose == T) {
          print('End of the burnin phase')
        }
      }
    }
    popn_founders <- currentPop()
    old.phenoRecord <- simEnv[["sims"]][[1]][["phenoRec"]]
    simEnv[["sims"]][[1]][["phenoRec"]] <- NULL
  } else{
    ## ----defineVariances-----------------------------------------------------
    defineVariances(
      gVariance = gVariance,
      locCorrelations = locCorrelations,
      gByLocVar = interact[1],
      gByYearVar = interact[2],
      fracGxEAdd = fracGxEAdd,
      plotTypeErrVars = plotTypeErrVars
    )
    ## ----initializePopulation------------------------------------------------
    initializePopulation(nInd = nInd)
    ## ----selfFertilize_founders----------------------------------------------
    for (s in 1:8) {
      selfFertilize(nProgeny = nInd)
    }
    popn_founders <- currentPop()
  }
  # Loop on years -----------------------------------------------------------
  popn_vec <- c()
  activity <- c()
  seasons <- c()
  year <- c()
  Cycle <- c()
  popn_select <- c()
  popn_evaluated <- c()
  popn_era <- c()
  #era <- NULL
  sel_diff <- c()
  exp_r <- c()
  Yr <- c()
  if (print_pop == T) {
    print(popSize())
  }
  
  if (verbose == T) {
    cat("------------------------------------------- \n")
    print("Start breeding cycles")
  }
  #i <- 1
  # i = i + 1
  for (i in seq(1, nyears, 2)) {
    if (i == 1) {
      Yr <- i
    } else{
      Yr <- c(Yr, (i - dplyr::last(Yr)))
    }
    if (verbose == T) {
      cat("\n ------------------------------------- \n")
      print(paste("Year: ", dplyr::last(Yr)))
      T1 <-  Sys.time()
      print(T1)
    }
    # Selections of parents ---------------------------------------------------
    if (i < (cyc_dur + 1)) {
      # sample parents from base until tested lines are avaliable
      if (verbose == T) {
        print("Select parents from founders population")
      }
      select(nSelect = nparents,
             popID = popn_founders,
             random = T)
      exp_r <- append(exp_r, 0)
      if (print_pop == T) {
        print(popSize())
      }
    } else {
      if (verbose == T) {
        print("Select parents from the previous cycle")
      }
      # determine which populations are available as selection candidates
      tab <- data.frame(seasons, activity, popn_vec, Cycle)
      prevsea <- tab[tab[, 1] %in% c(i - 2, i - 1),]#evaluated pop
      prev_test <-
        prevsea[grep('test', prevsea[, 'activity']),] %>%
        dplyr::filter(activity != 'test1')
      popn_test <-
        unique(prev_test$popn_vec) #popn_era[length(Yr)-5]
      popn_evaluated <- c(popn_evaluated, popn_test)
      # predict performance of the possible parents
      nyrincl <- as.numeric(ntest2 > 0) + as.numeric(ntest3 > 0)
      predictValue(popID = popn_test,
                   sharingInfo = sharingInfo,
                   years = c((i - nyrincl):(i - 1)))
      
      if (yearbyyear) {
        if (verbose == T) {
          print("Based the current cycle only")
        }
        # if selecting only based on previous year data
        
        if (sel_family) {
          select(nSelect = maxselPerfam,
                 popID = popn_test,
                 type = "WithinFamily")
          select(nSelect = nparents,
                 popID = currentPop(),
                 type = "Mass")
        } else {
          select(nSelect = nparents,
                 popID = popn_test)
        }
      } else {
        print("Based the all previous cycles")
        # if selecting based on all previous seasons data
        if (reuse_parent) {
          popn_evaluated <- append(popn_evaluated, popn_select)
        }
        predictValue(
          popID = popn_evaluated,
          sharingInfo = sharingInfo,
          years = c(1:c(i - 1))
        )
        
        if (sel_family) {
          select(nSelect = maxselPerfam,
                 popID = popn_evaluated,
                 
                 type = "WithinFamily")
          select(nSelect = nparents,
                 popID = currentPop(),
                 type = "Mass")
        } else {
          select(nSelect = nparents,
                 popID = popn_evaluated)
        }
      }
      print(Sys.time() - T1)
      if (print_pop == T) {
        print(popSize())
      }
      ## -----Calculate expected response (mean of initial population)
      gRec <- dataExtrac("genoRec")
      preds <- dataExtrac("predRec")
      preds <- preds[which(preds$predNo == max(preds$predNo)),]
      gidsel <- gRec[which(gRec$popID == currentPop()), "GID"]
      mnSel <- mean(preds[match(gidsel, preds$predGID), 3])
      mnTot <- mean(preds[, 3])
      sel_diff <-
        append(sel_diff, mnSel - mnTot)# diff. de selection
    }#END of the selection of parents
    
    popn_select <- append(popn_select, currentPop())
    season <- i
    # save information about the selected individuals in a table
    prec <- dataExtrac("genoRec")
    ind_sel <-
      which(prec$popID == currentPop())# GID of selected individuals
    selec_df <-
      data.frame(
        Selected = prec[ind_sel, "GID"],
        Basepop = prec[ind_sel, "basePopID"],
        Cycle = i,
        Season = season
      )
    if (i == 1) {
      selec_dfs <- selec_df
    } else {
      selec_dfs <- rbind(selec_dfs, selec_df)
    }
    
    # Crosses -----------------------------------------------------------------
    if (verbose == T) {
      print("Crosses")
    }
    cross(
      nProgeny = nCross,
      popID = currentPop(),
      equalContribution = T,
      notWithinFam = T
    )
    if (print_pop == T) {
      print(popSize())
    }
    popn_vec <- c(popn_vec, currentPop())
    popn_cross <- currentPop()
    seasons <- append(seasons, season)
    activity <- append(activity, "crossing")
    
    if (i == 1) {
      cycl <- rep(i, 12)
      Cycle <- cycl
    } else{
      cycl <- rep((i - unique(tail(Cycle, 1))), 12)
      Cycle <- append(Cycle, cycl)
    }
    
    # Eliminate some F1s ------------------------------------------------------
    if (verbose == T) {
      print("F1")
    }
    select(nSelect = nF1,
           random = T,
           popID = popn_cross)
    if (print_pop == T) {
      print(popSize())
    }
    popn_vec <- c(popn_vec, currentPop())
    activity <- append(activity, "F1")
    season <- season + 1
    seasons <- append(seasons, season)
    popn_cross <- currentPop()
    
    # RGA ---------------------------------------------------------------------
    if (verbose == T) {
      print("RGA")
    }
    for (k in 1:nRGA) {
      selfFertilize(nProgeny = ntest1)
      if (print_pop == T) {
        print(popSize())
      }
      popn_vec <- c(popn_vec, currentPop())
      activity <- append(activity, "selfing")
      if (k %% 2 != 0) {
        season <- season + 1
      }
      seasons <- append(seasons, season)
    }
    # Testing - first stage ---------------------------------------------------
    if (verbose == T) {
      print("Testing - first stage")
    }
    season <- season + 1
    seasons <- append(seasons, season)
    popn_vec <- c(popn_vec, currentPop())
    activity <- append(activity, "test1")
    phenotype(
      popID = currentPop(),
      nRep = rep1,
      plotType = plottype1,
      locations = loc1,
      years = season
    )
    if (agroTrend == T) {
      add_ymean(
        year = season,
        popID = NULL,
        YrMn = c(agro_trend * season)
      )
    }
    if (blup_lst) {
      # blup for selection among current candidates
      if (ntest2 + ntest3 != 0) {
        predictValue(popID = c(currentPop()),
                     sharingInfo = "none")
      }
    }
    if (verbose == T) {
      gval <- gvalData(currentPop())$gValue
      pval <- phenoData(currentPop())$pValue
      print("phenotyping - test1")
      print('Correlation = ')
      print(cor(gval, pval))
      print('h² =')
      print(var(gval) / var(pval))
      if (graph == T) {
        plot(
          gval,
          pval,
          xlab = 'genetic value',
          ylab = 'phenotypic value',
          main = c('test 1 \n' ,
                   paste("Corr = " , round(
                     cor(gval, pval), 2
                   ))),
          cex.main = 0.7
        )
      }
    }
    #popnobase <- currentPop()
    
    # Testing - second stage --------------------------------------------------
    if (ntest2 > 0) {
      if (verbose == T) {
        print("Testing - second stage")
      }
      
      if (within_fam) {
        select(
          popID = currentPop(),
          nSelect = max_per_fam2,
          random = F,
          type = "WithinFamily"
        )
        if (print_pop == T) {
          print(popSize())
        }
      }
      select(
        nSelect = ntest2,
        type = "Mass",
        random = F,
        popID = currentPop()
      )
      selfFertilize(nProgeny = ntest2)
      if (print_pop == T) {
        print(popSize())
      }
      
      popn_vec <- c(popn_vec, rep(currentPop(), 2))
      season <- season + 2
      years <- c(season, season + 1)
      seasons <- append(seasons, years)
      activity <- append(activity, c("test2", "test2"))
      phenotype(
        nRep = rep2,
        plotType = plottype2,
        locations = loc2,
        years = years
      ) # advanced lines 1
      if (agroTrend == T) {
        add_ymean(
          year = season,
          popID = NULL,
          YrMn = c(agro_trend * season)
        )
      }
      if (verbose == T) {
        gval <- gvalData(currentPop())$gValue
        pheno <- phenoData(currentPop())
        pval <- pheno %>%
          dplyr::group_by((phenoGID)) %>%
          dplyr::summarise(pval = mean(pValue))
        print("phenotyping - test2")
        print('Correlation = ')
        print(cor(gval, pval$pval))
        print('h² =')
        print(var(gval) / var(pval$pval))
        if (graph == T) {
          plot(
            gval,
            pval$pval,
            xlab = 'genetic value',
            ylab = 'phenotypic value',
            main = c('test 2 \n' ,
                     paste("Corr = " , round(
                       cor(gval, pval$pval), 2
                     ))),
            cex.main = 0.7
          )
        }
      }
      # blup for selection among current candidates
      predictValue (popID = currentPop(),
                    sharingInfo = sharingInfo)
    }
    # Testing - third stage ---------------------------------------------------
    if (ntest3 > 0) {
      if (verbose == T) {
        print("Testing - third stage")
      }
      if (within_fam) {
        select(nSelect = max_per_fam3,
               type = "WithinFamily")
        if (print_pop == T) {
          print(popSize())
        }
      }
      select(nSelect = ntest3,
             type = "Mass",
             popID = currentPop())
      selfFertilize(nProgeny = ntest3)
      if (print_pop == T) {
        print(popSize())
      }
      popn_vec <- c(popn_vec, rep(currentPop(), 2))
      season <- season + 2
      years <- c(season, season + 1)
      seasons <- append(seasons, years)
      activity <- append(activity, c("test3", "test3"))
      
      phenotype(
        nRep = rep3,
        plotType = plottype3,
        locations = loc3,
        years = years
      ) # advanced lines 2
      
      if (agroTrend == T) {
        add_ymean(
          year = season,
          popID = NULL,
          YrMn = c(agro_trend * season)
        )
      }
      if (verbose == T) {
        gval <- gvalData(currentPop())$gValue
        pheno <- phenoData(currentPop())
        pval <- pheno %>%
          dplyr::group_by((phenoGID)) %>%
          dplyr::summarise(pval = mean(pValue))
        print("phenotyping - test3")
        print('Correlation = ')
        print(cor(gval, pval$pval))
        print('h² =')
        print(var(gval) / var(pval$pval))
        if (graph == T) {
          plot(
            gval,
            pval$pval,
            xlab = 'genetic value',
            ylab = 'phenotypic value',
            main = c('test3 \n' ,
                     paste("Corr = " , round(
                       cor(gval, pval$pval), 2
                     ))),
            cex.main = 0.7
          )
        }
      }
      #Variaty registration for era trial
      select(nSelect = 2,
             type = "Mass",
             popID = currentPop())
      era <- currentPop()
      popn_era <- c(popn_era, era)
      popn_vec <- c(popn_vec, era)
      seasons <- append(seasons, season + 1)
      activity <- append(activity, "release")
    }
    
    if (i == max(seq(1, nyears, 2))) {
      if (verbose == T) {
        cat("------------------------------------------- \n")
        print("ERA trial")
      }
      # popn_vec <- c(popn_vec, rep(era,4))
      # seasons <- append(seasons, era_season)
      # activity <- append(activity, rep("ERA",4))
      # Cycle <- append(Cycle, rep(max(Cycle)+1, 4)
      
      era_season <- (max(years) + 1):(max(years) + 4)
      phenotype(
        popID = popn_era,
        nRep = rep3,
        plotType = plottype3,
        locations = loc3,
        years = era_season
      )
    }
    if (i == max(seq(1, nyears, 2))) {
      if (verbose == T) {
        cat("------------------------------------------- \n")
        print("End of breeding program")
        cat("------------------------------------------- \n")
        print("Result Analysis ... ")
        print("Please wait a few moments")
      }
    }
  } # End of breeding program
  #-------------------Get important information tables to be used later
  dt <- with(simEnv, sims)
  hdat <- dt[[1]]$phenoRec
  gdat <- dt[[1]]$genoRec
  
  #------------------Preapare the results table---------------------
  #Add the pop id
  pid <- gdat[match(hdat[, 1], gdat[, 1]), 'popID']
  bpid <- gdat[match(hdat[, 1], gdat[, 1]), 'basePopID']
  hdat <- cbind(hdat, pid, bpid)
  
  #Breeding activities
  tab_breed <- data.frame(seasons, activity, popn_vec, Cycle)
  
  #Get the gid year
  ugid <- unique(hdat$phenoGID)
  gidyr <- c()
  for (i in 1:length(ugid)) {
    gidyr <-
      append(gidyr, min(hdat[which(ugid[i] == hdat$phenoGID), 'year']))
  }
  gidyr <- gidyr[match(hdat$phenoGID, ugid)]
  hdat <- cbind(hdat, gidyr)
  
  #add the true breeding value
  gv <- dt[[1]]$gValue
  tbv <- gv[match(hdat$phenoGID, gdat$GID)]
  hdat <- cbind(hdat, tbv)
  
  
  # True breeding Value -----------------------------------------------------
  T_BV <-
    data.frame(matrix(NA, nrow = length(unique(Cycle)) + 1, ncol = 2))
  colnames(T_BV) <- c("Year", "T_BV")
  T_BV[1, 1] <- 0
  T_BV[1, 2] <-
    mean(gvalData(popn_founders)$gValue) - mean(gvalData(popn_founders)$gValue)
  for (b in 1:length(unique(Cycle))) {
    T_BV[b + 1, 1] <- b
    pop <-
      unique(
        tab_breed %>%
          dplyr::filter(Cycle == b & activity == 'test3') %>%
          dplyr::select(popn_vec)
      )$popn_vec
    T_BV[b + 1, 2] <-
      mean(gvalData(pop)$gValue) - mean(gvalData(popn_founders)$gValue)
  }
  
  # Heritability and Genetic Variance ---------------------------------------
  
  h2 <-
    matrix(
      NA,
      nrow = 3,
      ncol = length(unique(Cycle)),
      dimnames = list(c("Year", "OYT", "AYT"))
    )
  h2 <- as.data.frame(t(h2))
  
  gvar <-
    matrix(
      NA,
      nrow = 2,
      ncol = length(unique(Cycle)),
      dimnames = list(c("Year", "gVar"))
    )
  gvar <- as.data.frame(t(gvar))
  
  for (h in 1:nrow(h2)) {
    h2[h, 1] <- h
    gvar[h, 1] <- h
    #Second testing stage
    OYT <- unique(
      tab_breed %>%
        dplyr::filter(Cycle == h &
                        activity == 'test2') %>%
        dplyr::select(popn_vec)
    )$popn_vec
    
    pheno <- hdat %>%
      dplyr::filter(bpid == OYT) %>%
      dplyr::select(phenoGID, loc, year, pValue)
    
    pheno$phenoGID <- as.factor(pheno$phenoGID)
    pheno$year <- as.factor(pheno$year)
    
    mod_oyt <- lmer(pValue ~ 1 + year +
                      (1 | phenoGID),
                    data = pheno,
                    REML = TRUE)
    VarCompOYT <- as.data.frame(VarCorr(mod_oyt))[, -c(2, 3)]
    var_tot <-
      VarCompOYT[VarCompOYT$grp == "phenoGID", 'vcov'] +
      VarCompOYT[VarCompOYT$grp == "Residual", 'vcov'] / 2
    
    h2[h, 2] <- round(var(gvalData(OYT)$gValue) / var_tot, 2)
    #Third testing stage
    AYT <- unique(
      tab_breed %>%
        dplyr::filter(Cycle == h &
                        activity == 'test3') %>%
        dplyr::select(popn_vec)
    )$popn_vec
    gVar <- var(gvalData(AYT)$gValue)
    
    pheno <- hdat %>%
      dplyr::filter(bpid == AYT) %>%
      dplyr::select(phenoGID, loc, year, pValue)
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
    h2[h, 3] <-
      round(gVar / var_tot, 2)
    gvar[h, 2] <- gVar
    
  }
  
  #  Create objects to store information maf -------------------------------
  
  MAF <- matrix(NA, nrow = (length(unique(Cycle)) + 1), ncol = nMarkers)
  MAF[1, ] <- allFreq(genoData(popID = popn_founders))$maf
  
  for (m in 1:length(unique(Cycle))) {
    pop <- 
      unique(tab_breed %>% 
               dplyr::filter(Cycle==m & activity=='test3') %>% 
               dplyr::select(popn_vec))$popn_vec
    MAF[m+1 , ] <- allFreq(genoData(popID = pop))$maf
    
    if (m == length(unique(Cycle))) {
      MAF <- as.data.frame(t(MAF))
      colnames(MAF) <- 0:length(unique(Cycle))
    }
  }
  
  # Genetic gain - ERA Trial ------------------------------------------------
  
  Year <- rep(1:length(unique(Cycle)), each = 2)
  
  era_pheno <- hdat[hdat$year %in% era_season, ]
  
  ERA <- era_pheno %>%
    dplyr::group_by(phenoGID) %>%
    dplyr::summarise(pheno = mean(pValue))
  ERA <- cbind(Year, ERA)
  ERA <- ERA %>%
    dplyr::group_by(Year) %>%
    dplyr::summarise(pheno = mean(pheno)) %>%
    tibble::as_tibble()
  ERA$Year %<>%  as.numeric()
  era_gg <- lm(pheno ~ Year , data = ERA)
  GG <- coefficients(era_gg)[[2]]
  ERA$GGain <- GG
  ERA$PercentGain <- round(GG / ERA$pheno[1], 2) * 100
  
  # Genetic fain - EBV ------------------------------------------------------
  
  #Calculate the additive relionship matrix using pedigree records
  if (EBV == T) {
    pedigree <- gdat %>%
      dplyr::select(1:4)
    aMatrx <- calcAmatrix(pedigree)
    EBV <-
      matrix(
        NA,
        nrow = 2,
        ncol = length(unique(Cycle)),
        dimnames = list(c("Year", "ebv"))
      )
    EBV <- as.data.frame(t(EBV))
    
    for (b in 1:nrow(EBV)) {
      EBV[b, 1] <- b
      POP <- unique(
        tab_breed %>%
          dplyr::filter(Cycle == b &
                          activity == 'test3') %>%
          dplyr::select(popn_vec)
      )$popn_vec
      pheno <- hdat %>%
        dplyr::filter(bpid == POP) %>%
        dplyr::select(phenoGID, loc, year, pValue)
      pheno$phenoGID <- as.factor(pheno$phenoGID)
      pheno$loc <- as.factor(pheno$loc)
      pheno$year <- as.factor(pheno$year)
      
      mod1 <-
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
      blup <- ranef(mod1)[['phenoGID']]
      names(blup) <- 'phenoGID'
      indx <- rownames(blup)
      Amat <- aMatrx[indx, indx]
      mod2 <- regress(blup[, 1] ~ 1,  ~ Amat)
      ebv <- BLUP(mod2)
      EBV[b, 2] <- mean(ebv$Mean)
      
      EBV$Year %<>%  as.numeric()
      ebv_gg <- lm(ebv ~ Year , data = EBV)
      GG <- coefficients(ebv_gg)[[2]]
      EBV$GGain <- GG
      EBV$PercentGain <- round(GG / EBV$ebv[1], 2) * 100
    }
    
  }
  return(list(
    MAF = MAF,
    T_BV = T_BV,
    h2 = h2,
    gvar = gvar,
    ERA = ERA,
    B_genoFounders = genoFounders,
    B_MAF = B_MAF,
    B_gvar = B_gvar,
    EBV = EBV
  ))
}