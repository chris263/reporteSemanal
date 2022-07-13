
classLocations <- function(preAll, checks){
  # preAll <- preBLUEs_BU
  # checks <- checksF

  # checksF$Esperado[checksF$genotipo=="K9606VIP3"] <- 5
  # colnames(checks) <- c("genotipo","Check","Trait","Esperado","Tipo","TPP","valor" )

  phenoAll <- c()
  preAll <- preAll%>%mutate(Esperado = 0)
  preAll <- preAll%>%mutate(Check = 0)
  preAll <- preAll%>%mutate(Tipo="Obs")


  # myCheck <- unique(checks$genotipo)
  preCheck <- preAll %>% filter(genotipo %in% checks$genotipo)
  myCheck <- unique(preCheck$genotipo)

  if(length(myCheck)==0){stop("Nao foram encontrados checks :-( ")}


  # Separando os que nao sao checks
  nonCheck <- preAll%>%filter(!genotipo %in% myCheck)

  # Loop para separar por check
  for(i in 1:length(myCheck)){
    # i=1
    cat("Genotipo ", myCheck[i], i, "\n")

    checkPheno <- preAll %>% filter(genotipo == myCheck[i])

    if(nrow(checkPheno)>0){
      notaCheck <- as.numeric(unique(checks %>% filter(genotipo == myCheck[i]) %>% dplyr::select(Esperado)))
      checkPheno$Esperado <- notaCheck
      checkPheno$Check <- 1
      if(notaCheck < 3){
        checkPheno$Tipo <- "T"
      }else if(notaCheck>=3 && notaCheck<5){
        checkPheno$Tipo <- "MT"
      }else if(notaCheck>=5 && notaCheck<7){
        checkPheno$Tipo <- "MS"
      }else{
        checkPheno$Tipo <- "S"
      }
      phenoAll <- rbind(phenoAll,checkPheno)
    }else{
      cat("Sem check ",i,"\n")
    }



  }


  # Adicionando os demais materiais que não são checks
  phenoAll <- rbind(phenoAll,nonCheck)

  # colnames(phenoAll)[1] <- "blueN"
  phenoAll$blueN <- as.numeric(phenoAll$blueN)
  phenoAll$Esperado <- as.numeric(phenoAll$Esperado)

  # Pegando apenas os locais que tem checks
  comChecks <- as.data.frame(tapply(phenoAll$Check,phenoAll$local,sum))
  colnames(comChecks)<-c("blups")
  comChecks <- comChecks %>% filter(comChecks$blups != 0)

  # Vetor com os nomes dos locais que possuem checks
  nomesLocais <- rownames(comChecks)

  if(length(nomesLocais) ==0){
    cat("Sem classificação...","\n")
    break
    }

  # Criando um dataframe para armazenar o resultado
  classFinal <- data.frame("Local"=as.character(nomesLocais),
                           "classe"= as.character("classe"),
                           "Probabilidade"= as.numeric(0),
                           "Controles"=as.character("controles"),
                           "Desvio"=as.numeric(0))

  # Classificando os locais

  # match("TRUE", str_detect(nomesLocais,"BU07"))

  for(i in 1:length(nomesLocais)){
    # i=1
    # cat("Local: ",i,nomesLocais[i],"\n")
    # Filtrando por local

    selectedLocal <- phenoAll %>% filter(local == nomesLocais[i])
    testSeason <- unique(stri_sub(selectedLocal$local, 3,4))
    selChecks <- checks %>% filter(Season %in% testSeason, TPP == unique(selectedLocal$TPP))
    checksTPP <- selectedLocal %>% filter(genotipo %in% selChecks$genotipo)

    classFinal$Local[i] <- nomesLocais[i]

    #Verificando se possui checks para as 4 classes
    tipos <- as.data.frame(tapply(checks$Esperado,checks$Tipo,mean))
    tipos$tipo <- rownames(tipos)
    colnames(tipos) <- c("notas","tipos")
    ansTipos <- c("NA","NA","NA","NA")
    ansTipos[1] <- "T"  %in% tipos$tipos
    ansTipos[2] <- "MT" %in% tipos$tipos
    ansTipos[3] <- "MS" %in% tipos$tipos
    ansTipos[4] <- "S"  %in% tipos$tipos

    # Vetor com os checks
    # checks <- selectedLocal %>% filter(Check == 1, Esperado == 4 | Esperado == 5 )
    # ChecksTPP <- selectedLocal %>% filter(Check == 1, )
    chiDF <- data.frame(tapply(checksTPP$blueN,checksTPP$genotipo,mean)) %>% data.frame(tapply(checksTPP$Esperado,checksTPP$genotipo,mean))
    colnames(chiDF) <- c("Obs","Esp")

    # Calculando chi quadrado
    chiDF$X2 <- (chiDF$Obs - chiDF$Esp)^2/chiDF$Esp
    if(nrow(chiDF==1) | nrow(chiDF==0)){
      GL=1
    }else{
      GL=nrow(chiDF)-1
    }
    classFinal$Probabilidade[i] <- 1-pchisq(sum(chiDF$X2),GL)

    # Soma das diferenças para ajudar na classificação
    chiDF$Diff <- chiDF$Obs - chiDF$Esp
    classFinal$Desvio[i] <- sum(chiDF$Diff) # para auxiliar na classificação

    # Armazenando os tipos de checks do local
    classCheck = ""
    for(z in 1:length(tipos$tipos)){
     # z=1
      classCheck1 <- tipos$tipos[z]
      classCheck <- paste0(classCheck1," ",classCheck)
    }

    classCheck <- na.omit(classCheck)
    classFinal$Controles[i] <- classCheck



  }
  classFinal$TPP <- unique(selectedLocal$TPP)
  classFinal <- classFinal %>% arrange(desc(Desvio))
  classFinal <- na.omit(classFinal)
  return(classFinal)
}


finalClass<- function(enterTable, qualC, limA, limB){
  # enterTable = novaTabela
  # qualC=qControl
  # limA=limiteAlta
  # limB=limiteBaixa

  # enterTable <- na.omit(enterTable)
  outTabela <- enterTable

  for(i in 1:nrow(outTabela)){
    # i=20
    # print(i)
    somaDiff <- outTabela$Desvio[i]
    qc <- outTabela$Acuracia[i]

    if(is.na(qc)==T){qc=1}

    if( somaDiff > limA && qc > qualC){

      outTabela$classe[i] <- "Alta"

    }else if(somaDiff > limB && somaDiff < limA){

      outTabela$classe[i] <- "Media"

    }else if (somaDiff < limB) {

      outTabela$classe[i] <- "Baixa"

    }else if(somaDiff > limA && qc < qualC){

      outTabela$classe[i] <- "Low QC"

    }
  }

  # nomesRef <- distinct(BLUEs,local,stage)
  # outTabela2 <- left_join(outTabela,nomesRef, by=c("Local"="local"))
  # outTabela2 <- outTabela2 %>% select(Local,stage,classe,Probabilidade,
  #                                     Controles,Desvio,H2N,H2Y,Acuracia)
  # outTabela2$Probabilidade <- round(outTabela2$Probabilidade,3)

 return(outTabela)

}

joint_sommer <- function(inBLUE,locP){
  # inBLUE = altaPPBU
  # locP = "yes"


  library(sommer)
  inBLUE$BUGeno <- paste0(inBLUE$BU,"_",inBLUE$genotipo)
  inBLUE$Nota <- as.numeric(inBLUE$Nota)

  if(locP=="yes"){
    modeloSommer <- mmer(Nota ~ BU,
                           random= ~ vsr(genotipo, Gtc=unsm(1)) + vsr(BUGeno, Gtc=unsm(1)),
                           rcov= ~ vsr(units, Gtc=unsm(1)),
                           data=inBLUE, verbose = FALSE)
  }else{
    modeloConjunto <- mmer(Nota ~ 1,
                           random= ~ vsr(genotipo, Gtc=unsm(1)) + vsr(localGeno, Gtc=unsm(1)),
                           rcov= ~ vsr(units, Gtc=unsm(1)),
                           # tolparinv = 1e100000,
                           data=inBLUE, verbose = FALSE)
  }


  teste <- as.data.frame(modeloConjunto$U$`u:genotipo`$Nota)
  testeS <- as.data.frame(modeloSommer$U$`u:genotipo`$Nota)
  testeS$genotipo <- rownames(testeS)

  # randef(modeloSommer)
  # summary(modeloSommer)
  #
  # modeloSommer$U$`u:genotipo`

  return(modeloConjunto)
}

joint_blup <- function(inData, locP){
  # inData = altaPPBU
  # locP = "yes"

  inData$Nota <- as.numeric(inData$Nota)

  if(locP == "yes"){
      modeloAllNota <- lmer(Nota ~  BU + (1| genotipo), data = inData)
  }else{
    modeloAllNota <- lmer(Nota ~  1 + (1| genotipo), data = inData)
  }

  BlupNota <- coef(modeloAllNota)$genotipo
  BlupNota <- BlupNota[1]
  finalNota <- cbind(rownames(BlupNota),BlupNota$`(Intercept)`)
  colnames(finalNota) <- c("genotipo","blupNota")
  finalNota <- as.data.frame(finalNota)

  finalData <- finalNota %>% arrange(desc(blupNota))
  preResult <- left_join(finalData,BLUEs,by=c("genotipo"="genotipo"))
  preResult$blupNota <- as.numeric(preResult$blupNota)
  preResult <- preResult %>% arrange(blupNota)

  Result <- preResult %>% distinct(genotipo, blupNota) %>% dplyr::select(genotipo,blupNota)

  return(Result)

}

#
# testeAOV <- aov(Nota ~  year + season + BU + genotipo, data = altaPPBU)
# anova(testeAOV)








