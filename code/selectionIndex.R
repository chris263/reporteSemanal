

cleaning_data <- function(inFile, crit){
  # inFile <- altaPP
  # crit <- "yes"

  if(crit=="yes"){
    listaZero <- as.data.frame(tapply(inFile$yield, inFile$genotipo, mean))
    colnames(listaZero) <- "yield"
    listaZero <- listaZero %>% filter(yield != 0)
    listaZero$Nomes <- rownames(listaZero)

    # Vetor com os genotipos que permanecem
    vZero <- listaZero$Nomes

    # Filtrando apenas para os nomes que devem permanecer
    preOut <- inFile[which(inFile$genotipo%in%vZero),]
    # Removendo dados errados, baixa nota e baixíssimo Yield
    # Lista de plots a serem removidos
    preOut2 <- preOut %>% filter(Nota <= 3, yield < 6)
    plotRM <- preOut2$barcode
    outFile <<- preOut[-which(preOut$barcode%in%plotRM),]
  }else{
    preOut2 <- inFile %>% filter(Nota <= 3, yield < 6)
    plotRM <- preOut2$barcode
    outFile <<- inFile[-which(inFile$barcode%in%plotRM),]
  }


}

local_blup <- function(inLocal, crt4, inDFF, inType){
  # inLocal = posDfLocal$local
  # crt4="no"
  # inDFF <- myDFB
  # inType = "BU"

  inDFF$Nota <- as.numeric(inDFF$Nota)
  resumoLoc <- data.frame(local = c(1:length(inLocal)),H2N=rep(0,length(inLocal)),
                           H2Y=rep(0,length(inLocal)),Acuracia=rep(0,length(inLocal)))


  if(crt4 == "no"){resumoLoc[,-3]}
  library(lme4)
  BLUEf <- c()


  for(i in 1:length(inLocal)){
    # i=1
    BLUPN <- c()
    BLUEY <- c()

    if(inType == "BU"){
      selectedLoc <- inDFF%>%filter(BU==inLocal[i])
    }else if(inType == "trial"){
      selectedLoc <- inDFF%>%filter(local==inLocal[i])
    }

    selectedLoc$Nota <- as.numeric(selectedLoc$Nota)

    testeNota <- sum(selectedLoc$Nota, na.rm = T)

    if(testeNota == 0 | is.na(testeNota) == T){cat("Local", i," ", inLocal[i], "sem notas.\n");next}
    # selectedLoc <- myDF%>%filter(local=="21WNBPEYG6J1BT122021")

    testeY <- sum(selectedLoc$yield)

    if (is.na(testeY) == TRUE | testeY == 0){
      crt5 = "no"
      # cat("Aqui ok ", i, "\n")
    }else{
      crt5 = "yes"
      # cat(crt5, i, "\n")
      }



    # selectedLoc <- selectedLoc %>% filter(genotipo %in% genFora$Nomes)
    nReps <- max(as.numeric(selectedLoc$rep))

    testeRep <- is.na(nReps)
    if(testeRep == T){ selectedLoc$Rep <- 1}

    # colnames(selectedLoc)[4]<-"Nota"
    testeVAR <- var(selectedLoc$Nota, na.rm=T)
    checkVar <- is.na(testeVAR)
    if(checkVar==T){next}

    # cat("The numer is", i, " ",nrow(selectedLoc),"\n")

    options(error = NULL)
    if (testeVAR == 0){
      # Guardando os dados
      resumoLoc$local[i] <- inLocal[i]
      resumoLoc$H2N[i] <- 0
      genotipos <- unique(selectedLoc$genotipo)
      BLUPN <- data.frame(blueN=rep(selectedLoc$Nota[1],length(genotipos)),
                                   local=rep(unique(selectedLoc$local),length(genotipos)),
                                   genotipo=genotipos)
    }else{

      selectedLoc$rep <- as.factor(selectedLoc$rep)
      if(length(unique(selectedLoc$rep)) == 1){
        try(modeloLocalN <- lmer(Nota~1+(1|genotipo), data=selectedLoc), silent = T)
      }else{
        try(modeloLocalN <- lmer(Nota~rep+(1|genotipo), data=selectedLoc), silent=T)
      }
      errorCheck <- exists("modeloLocalN")
      if(errorCheck == F){ next }


      # unique(selectedLoc$genotipo)
      # selectedLoc$Nota
      #
      # Capturando as variancias para calcular H2
      variance = as.data.frame(VarCorr(modeloLocalN))
      gvar = variance [1,"vcov"]
      resvar = variance [2, "vcov"]

      if(gvar == 0 ){
        H2 = testeVAR/ (testeVAR + resvar)
      }else{
        H2 = gvar/ (gvar + resvar)
      }


      if(crt5 == "yes"){
        selectedLoc$yield <- as.numeric(selectedLoc$yield)
        selectedY <- selectedLoc %>% filter(yield>0)


        if(testeY>0){
          try(modeloLocalY <- lmer(yield~Rep+(1|genotipo), data=selectedLoc))
          teste2<-exists("modeloLocaly")
          if(teste2==F){
            cat("Sem analise de yield...","\n")
            next
          }
          varianceY = as.data.frame(VarCorr(modeloLocalY))
          gvarY = varianceY [1,"vcov"]
          resvarY = varianceY [2, "vcov"]
          H2Y = gvarY/ (gvarY + resvarY)

          # Armazenando os dados de H2
          resumoLoc$local[i] <- inLocal[i]
          resumoLoc$H2N[i] <- H2
          resumoLoc$H2Y[i] <- H2Y
          # resumoLoc$Acuracia[i] <- (H2+H2Y)/2

          # Extraindo os blups para Yield
          adjY = coef(modeloLocalY)$genotipo
          adjY = adjY[1]
          adjY =cbind(adjY,inLocal[i],rownames(adjY))

          colnames(adjY)<-c("blueY","local","genotipo")

          locaisNomeY <- c(1:nrow(adjY))
          BLUEY <- rbind(BLUEY,adjY)

        }else{ # caso não tenha yield
          H2Y = H2

          adjY = coef(modeloLocalN)$genotipo
          adjY[1] = NA
          adjY = cbind(adjY[1], inLocal[i],rownames(adjY))
        }

      }else{ # caso crt5 = no
        resumoLoc$local[i] <- inLocal[i]
        resumoLoc$H2N[i] <- H2
        # resumoLoc$Acuracia[i] <- H2
      }


      # Extraindo os blups para Nota
      adj = coef(modeloLocalN)$genotipo
      adj = adj[1]
      adj =cbind(adj,inLocal[i],rownames(adj))

      colnames(adj)<-c("blueN","local","genotipo")

      locaisNome <- c(1:nrow(adj))
      BLUPN <- rbind(BLUPN,adj)
      BLUPN <- BLUPN %>% dplyr::select(genotipo,local,blueN)
    }
    if(crt5 == "yes"){
      preBLUE <- left_join(BLUPN, BLUEY, by=c("genotipo"="genotipo")) %>% dplyr::select(genotipo,local.x,blueN,blueY)
      colnames(preBLUE) <- c("genotipo", "local", "blueN","blueY")
      BLUEf <- rbind(BLUEf,preBLUE)
      colnames(BLUEf) <- c("genotipo", "local", "blueN","blueY")
    }else{
      BLUPN$blueY <- NA
      BLUEf <- rbind(BLUEf,BLUPN)

    }

  }

  if(crt4 == "yes"){
    BLUEf <- BLUEf %>% dplyr::select(genotipo, local, blueN, blueY)
    colnames(BLUEf) <- c("genotipo", "local", "blueN","blueY")
  }else{
    # BLUEf <- BLUPN
    if(is.null(BLUEf)==T){
      print("No BLUE data...")
      # break
    }else{
      try(BLUEf <- BLUEf %>% dplyr::select(genotipo,local,blueN))
    }

  }

  resumoLoc <<- resumoLoc
  return(BLUEf)

}

calc_blups <- function(inModel,cr){

  # inModel <- modelo1
  # cr <- "no"

  if(is.null(inModel$U$`u:genotipo`$cNota)){
    inBlupNota <- as.data.frame(inModel$U$`u:genotipo`$blueN)
  }else{
    inBlupNota <- as.data.frame(inModel$U$`u:genotipo`$cNota)
  }

  inBlupNota$genotipo <- rownames(inBlupNota)
  colnames(inBlupNota) <- c("fNota","genotipo")
  inBlupNota <- inBlupNota %>% arrange(genotipo)

  if(cr == "yes"){
    inBlupYield <- as.data.frame(inModel$U$`u:genotipo`$cYield)
    inBlupYield$genotipo <- rownames(inBlupYield)
    colnames(inBlupYield) <- c("fYield","genotipo")
    inBlupYield <- inBlupYield %>% arrange(genotipo)
    # inBlupYield <- inBlupYield %>% filter(fYield>0)
    fBLUP <- left_join(inBlupYield,inBlupNota,by=c("genotipo"="genotipo"))
    fBlup <<- fBLUP %>% relocate(genotipo, .before = fYield)
  }else{
    fBLUP <<- inBlupNota
  }
}

  defaultW <- getOption("warn")
  options(warn = -1)

  efeitoNY <- function(effDF){
    # effDF <- altaPP
    aa <- as.numeric(effDF$Nota)
    bb <- as.numeric(effDF$yield)
    effDF <- altaP_DF

  lineaR <- lm(yield~Nota, data=effDF)

  # y = bx+a
  a <- coef(lineaR)[1]
  b <- coef(lineaR)[2]

  # Se aumentar a nota de 3 para quantro, qual o efeito em yield?
  pYield1 <- b*3+a #nota 3
  pYield2 <- b*4+a #nota 4

  resDiff<- round(pYield2 - pYield1,2)
  names(resDiff)<-c("yield")
  return(resDiff)
}



selectionIndex <- function(inModelo, inWgt){
  # inModelo <- modeloConjunto
  # inWgt <- pesos

  # Specify the additive variance and correlation: 1 on the diagonal
  conjuntaV <- summary(inModelo)$varcomp
  addVar <- c(Yield=conjuntaV$VarComp[1],CS=conjuntaV$VarComp[3])
  calcCor <- conjuntaV$VarComp[2]/(sqrt(conjuntaV$VarComp[1])*sqrt(conjuntaV$VarComp[3]))
  addCor <-  matrix(c(1, calcCor, calcCor, 1), nrow=2)
  addSD <- diag(sqrt(addVar)) # Convert cor to SD
  addCov <- addSD %*% addCor %*% addSD

  # Specify the error correlation and calculate error covariance
  errVar <- c(Yield=conjuntaV$VarCompSE[1],CS=conjuntaV$VarCompSE[3])
  calcErr <- conjuntaV$VarCompSE[2]/(sqrt(conjuntaV$VarCompSE[1])*sqrt(conjuntaV$VarCompSE[3]))
  errCor <- matrix(c(1, calcErr, calcErr, 1), nrow=2)
  errSD <- diag(sqrt(errVar))
  errCov <- errSD %*% errCor %*% errSD

  phenCov <- addCov + errCov

  #Indice: b=inv(P)Av
  indexSE <- solve(phenCov)%*%addCov%*%inWgt

  rownames(indexSE) <- c("Yield","CS") # Ordem dos pesos
  indexSE[2] <- (-1)*indexSE[2] # Invertido porque queremos a menor nota

  # Acuracia do Index
  preAcc <- t(indexSE)%*%addCor%*%inWgt
  varH <- t(inWgt)%*%addCor%*%inWgt
  accIndex <<- sqrt(preAcc+varH)

  return(indexSE)

}


ranking_data <- function(inBlup,inDFF,cr2,corr2, indexSel){
  # inBlup<-finalBLUP
  # inDFF<-altaP_DF
  # cr2<-duasT
  # corr2 <- corrl
  # indexSel <- pesos

  # Organizando o Blup para Nota
  if(length(unique(inDFF$local))>1){
    modeloAllNota <- lmer(Nota ~  local + (1| genotipo), data = inDFF)
  }else{
    modeloAllNota <- lmer(Nota ~  1 + (1| genotipo), data = inDFF)
  }

  BlupNota <- coef(modeloAllNota)$genotipo
  BlupNota <- BlupNota[1]
  finalNota <- cbind(rownames(BlupNota),BlupNota$`(Intercept)`)
  colnames(finalNota) <- c("genotipo","blupNota")
  finalNota <- as.data.frame(finalNota)

  if(cr2=="yes"){ # Se for analisar Yield em conjunto
    matrixBlup <- matrix(c(as.numeric(inBlup$fYield),as.numeric(inBlup$fNota)), ncol = 2)

    inBlup$Indice <-  matrixBlup %*% indexSel

    if(length(unique(inDFF$local))==1){
      modeloAllYield <- lmer(yield ~ 1 + (1| genotipo), data = inDFF)
    }else{
      modeloAllYield <- lmer(yield ~ local + (1| genotipo), data = inDFF)
    }

    BlupYield <- coef(modeloAllYield)$genotipo
    BlupYield <- BlupYield[1]
    finalYield <- cbind(rownames(BlupYield),BlupYield$`(Intercept)`)
    colnames(finalYield) <- c("genotipo","blupYield")
    finalYield <- as.data.frame(finalYield)
    finalData <- left_join(finalNota,finalYield, by=c("genotipo"="genotipo"))
    preResult <- left_join(finalData,inBlup,by=c("genotipo"="genotipo")) %>% dplyr::select(genotipo,blupNota, blupYield,Indice)
    # preResult <- na.omit(preResult)
    preResult$blupNota <- as.numeric(preResult$blupNota)
    preResult$blupYield <- as.numeric(preResult$blupYield)
    preResult$Indice <- as.numeric(preResult$Indice)
    if(corr2 < (-0.95)){
      preResult$Indice <- (-1)*preResult$Indice
    }
    preResult <- preResult %>% arrange(desc(Indice))
  }else{
    finalData <- finalNota
    preResult <- left_join(finalData,inBlup,by=c("genotipo"="genotipo"))
    preResult$blupNota <- as.numeric(preResult$blupNota)
    preResult <- preResult %>% arrange(blupNota)
  }

  # preResult$Indice <- -1*preResult$Indice
  return(preResult)
}

qualityGeno <- function(entrada, retd, cr3, crt6){
  # entrada <- dadosGenoLoc
  # retd <- resultados
  # cr3 <- duasT
  # crt6 <- minLoc

  f10 <- entrada[1,] # 1 linha adicionada apenas para criar o df.
  f10$nLocais <- as.numeric(0)
  f10$genoQC <- as.numeric(0)

  f0 <- unique(entrada$genotipo)

  for(i in 1:length(f0)){
    # i=2
    nomeGenotipo <- f0[i]
    f1 <- entrada %>% filter(genotipo == nomeGenotipo) %>% filter(Rep ==1)
    nLinhas <- length(unique(f1$local))

    f2 <- f1[1,]
    f2$nLocais <- as.numeric(nLinhas)
    f2$genoQC <- as.numeric(0)

    # Capturando as medias de acuracia
    f3 <- entrada %>% filter(genotipo == nomeGenotipo)

    f3m <- tapply(f3$Acuracia, f3$genotipo,mean, na.omit=T)

    if(nLinhas > 0){
    f2$genoQC <- round(f3m[[1]]*(1-(1/2^(nLinhas))),3)
    }else{
      f2$genoQC[1] <- "NA"
    }

    # if(nLinhas > crt6){
    #   f2$genoQC[1] <- round(sum(as.numeric(f1$Acuracia))/nLinhas,2)
    # }else if (nLinhas > 0 && nLinhas <= crt6){
    #   accT <- f1$Acuracia
    #   accTF <- accT[1:nLinhas]
    #   f2$genoQC[1] <- round(sum(as.numeric(accTF))/crt6,2)
    # }else if (entrada$nLocais[i] == 0){
    #   f2$genoQC[1] <- "NA"
    # }

    f10 <- rbind(f10,f2)

  }

  f10 <- f10[-1,] # Retira a primeira linha que adicionei apenas para criar o df.

  preFINAL <- left_join(f10,retd,by=c("genotipo"="genotipo"))

  if(cr3 == "yes"){
    resFINAL <- data.frame(preFINAL$genotipo, as.numeric(preFINAL$blupNota), as.numeric(preFINAL$blupYield),
                           preFINAL$Indice,preFINAL$nLocais, preFINAL$genoQC)
    colnames(resFINAL) <- c("genotipo","blupNota","blupYield","Indice","nLocais","Acuracia")
    resFINAL$blupNota <- floor(resFINAL$blupNota)
    resFINAL$blupYield <- round(resFINAL$blupYield,2)
    resFINAL <- resFINAL %>% arrange(desc(Indice))
  }else{
    resFINAL <- data.frame(preFINAL$genotipo, as.numeric(preFINAL$blupNota), blupYield=0,
                           Indice=0,preFINAL$nLocais, preFINAL$genoQC)
    colnames(resFINAL) <- c("genotipo","blupNota","blupYield","Indice","nLocais","Acuracia")
    resFINAL$blupNota <- round(resFINAL$blupNota,0)
    resFINAL <- resFINAL %>% arrange(blupNota)
  }
  resFINAL$Class <- NA
  for(i in 1:nrow(resFINAL)){
    if(is.na(resFINAL$blupNota[i])){
    }else{
      if(resFINAL$blupNota[i] < 3){
        resFINAL$Class[i] <- "T"
      }else if(resFINAL$blupNota[i]>=3 && resFINAL$blupNota[i]<5){
        resFINAL$Class[i] <- "MT"
      }else if(resFINAL$blupNota[i]>=5 && resFINAL$blupNota[i]<=6){
        resFINAL$Class[i]<- "MS"
      }else if(resFINAL$blupNota[i]>6){
        resFINAL$Class[i] <- "S"
      }
    }

  }

  # resFINAL <<- na.omit(resFINAL)
  return(resFINAL)
}

calc_blupB <- function(baixaP_DF){
  # Calculando BLUP para baixa pressão
  baixaP_DF$yield <- as.numeric(baixaP_DF$yield)
  baixaP_DF$local <- as.factor(baixaP_DF$local)
  baixaP_DF$Rep <- as.factor(baixaP_DF$Rep)

  modeloB <- lme4::lmer(yield ~ local + Rep + (1|genotipo), data = baixaP_DF)

  blupB <- coef(modeloB)$genotipo
  blupB <- blupB[1]
  finalB <- cbind(rownames(blupB),blupB$`(Intercept)`)
  colnames(finalB) <- c("genotipo","blupB")
  finalB <- as.data.frame(finalB)
  return(finalB)
}
