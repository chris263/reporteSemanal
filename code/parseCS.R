

library(stringr)
read_fun <- function(inLista){
  inLista = filenames
  finalData <- data.frame(TPP=NA,
                          local=NA,
                          genotipo=NA,
                          plotID=NA,
                          rep=NA,
                          COSTR=NA,
                          HELMR=NA,
                          GRLSR=NA,
                          LFSPR=NA,
                          SCLBR=NA,
                          SRSTR=NA,
                          DLLFR=NA,
                          PRMDN=NA,
                          yieldGSM=NA,
                          range=NA,
                          row=NA,
                          BU=NA,
                          stg=NA,
                          season=NA,
                          barcode=NA,
                          year=NA
  )

  posFinal <- finalData

  for(i in 1:length(inLista)){
    # i=2
    inFile <- inLista[i]
    # inFile <- gsub(pattern = "../",replacement = "",x = inFile,fixed = T)

    checkError <- FALSE
    tryCatch( { readData <- read.csv(file=inFile, header=T, as.is=1, sep=",") }
              , warning = function(w) { cat("Arquivo \"",inFile,"\" nao encontrado.\n"); checkError <- TRUE })

    if(checkError == T){stop("Arquivo nao encontrado, checar o local do arquivo.")}

    readData$TPP <- NA
    nomesCol <- colnames(readData)
    if(length(nomesCol)==1){ readData <- read.csv(file=inFile, header=T, as.is=1, sep=";")}

    selNames <- c("TPP","EXTID","ABBRC","PLTID","REPNO","COSTR","HELMR","GRLSR","LFSPR","SCLBR","SRSTR","DLLFR","PRMDN","PSS.YGSMN","STRNG","STROW","PSS.BARCD")
    checkNames <- match("TRUE",str_detect(nomesCol,paste(selNames[2:5],collapse = "|")))
    if(is.na(checkNames)==T){stop("Stop! Esta faltando umas das colunas principais, checar arquivo.")}

    newData <- readData

    count = 1
    for(i in 1:length(selNames)){
      # i=1
      # tem que rodar o loop completo para renomear todas as colunas corretamente
      posName <- match("TRUE",str_detect(nomesCol,selNames[i]))
      if(is.na(posName==T)){
        cat("Esta faltando a coluna",selNames[i],"no arquivo", inFile ,"adicionando como NA.\n")
        newData[,length(nomesCol)+count] <- NA
        colnames(newData)[length(nomesCol)+count] <- selNames[i]
        count=count+1
      }else{
        colnames(newData)[posName] <- selNames[i]
      }

    }
    newData <- newData %>% select(all_of(selNames))

    colnames(newData) <- c("TPP","local","genotipo","plotID","rep","COSTR","HELMR",
                           "GRLSR","LFSPR","SCLBR","SRSTR","DLLFR","PRMDN","yieldGSM",
                           "range","row","barcode")


    # capturando as TPPs
    newData$TPP <- stri_sub(newData$local, -4,-3)
    newData$TPP <- gsub("BL","TPP09",newData$TPP)
    newData$TPP <- gsub("BU","TPP11",newData$TPP)
    newData$TPP <- gsub("BT","TPP10",newData$TPP)
    newData$TPP <- gsub("BC","TPP12",newData$TPP)


    newData$BU <- NA
    newData$stg <- NA
    newData$season <- NA
    newData$year <- NA

    finalData <- rbind(finalData,newData)

  }
  finalData <- finalData[-1,]

  todosLoc <- unique(finalData$local)
  for (i in 1:length(todosLoc)) {
    # i=325
    filterLoc <- finalData %>% filter(local %in% todosLoc[i])
    filterLoc$season <- substring(filterLoc$local,3,4)
    filterLoc$stg <- str_sub(filterLoc$local,-7,-7)
    filterLoc$stg <- gsub(8,6,filterLoc$stg)
    filterLoc$stg <- gsub("L",4,filterLoc$stg)
    checkStage <- unique(filterLoc$stg)

    if(grepl("\\D",checkStage) == T){
      cat("Alterando stage local ", todosLoc[i],"\n")
      filterLoc$stg <- str_sub(filterLoc$local,-6,-6)
    }

    filterLoc$year <- str_sub(filterLoc$local, 1,2)
    filterLoc$BU <- paste0(filterLoc$year,filterLoc$season,filterLoc$stg,str_sub(filterLoc$local,-4,-1))

    posFinal <- rbind(posFinal, filterLoc)
  }

  posFinal <- posFinal %>% filter(stg %in% c(4,5,6))

  posFinal <- posFinal[-1,]
  posFinal <- tidyr::gather(posFinal, trait, Nota, COSTR:yieldGSM)
  posFinal$Nota <- as.numeric(posFinal$Nota)
  posFinal$year <- as.numeric(posFinal$year)
  posFinal$stg <- as.numeric(posFinal$stg)

  # finalData$Nota <- as.numeric(finalData$Nota)
  posFinal <- posFinal %>% filter(is.na(Nota) == FALSE)
  posFinal <- distinct(posFinal, barcode, trait,rep,.keep_all = T)

  colnames(posFinal)[1] <- "TPP"
  colnames(posFinal)[8] <- "BU"
  colnames(posFinal)[14] <- "Nota"

  return(posFinal)

}

transCHEC <- function(inCHEC){
  sinN <<- distinct(sinN, genotipo,Novo)
  inCHEC$Esperado <- as.numeric(inCHEC$Esperado)
  errorMessage=FALSE
  tryCatch(tidyr::gather(inCHEC,TPP,valor,7:10),
           error = function(e)
             errorMessage <<- TRUE)

  if(errorMessage==T){
    inCHEC <- inCHEC %>% mutate(TPP09=1,
                                    TPP10=1,
                                    TPP11=1,
                                    TPP12=1)

    inCHEC <- tidyr::gather(inCHEC,TPP,valor,7:10)
  }else{
    inCHEC <- tidyr::gather(inCHEC,TPP,valor,7:10)
  }
}


sinoFun <- function(myData, sinIN){
  # sinIN <- sinN
  # myData <- checksCS

  colnames(myData) <- tolower(colnames(myData))
  colnames(sinIN) <- c("genotipo", "novo")

  sinIN$novo <- gsub(" ","",sinIN$novo)
  sinIN$novo <- gsub("\\.","",sinIN$novo)
  sinIN$genotipo <- gsub(" ","",sinIN$genotipo)
  sinIN$genotipo <- gsub("\\\\","_",sinIN$genotipo)
  sinIN$genotipo <- gsub("\\.","",sinIN$genotipo)
  sinIN$genotipo <- toupper(sinIN$genotipo)

  myData$genotipo <- gsub(" ","", myData$genotipo)
  myData$genotipo <- gsub("\\\\","_", myData$genotipo)
  myData$genotipo <- gsub("\\.","",myData$genotipo)
  myData$genotipo <- toupper(myData$genotipo)

  if(ncol(myData)>9){
    myData$local <- gsub(" ","",myData$local)
    myData$local <- toupper(myData$local)
  }


  for( i in 1:nrow(sinIN)){
    # i=30
    # print(i)
    myOld <- sinIN$genotipo[i]
    myNew <- sinIN %>% filter(genotipo == myOld) %>% dplyr::select(novo)

    if(nrow(myNew)>1){
      cat("Genotipo ", myOld," possui mais de 1 sinonimo","\n")
      next
    }
    getPosition <- grep(myOld,myData$genotipo)
    myData$genotipo <- replace(myData$genotipo, getPosition, myNew)
  }
  myData$genotipo <- gsub("\\_","\\\\", myData$genotipo)

  if(ncol(myData) > ncol(checksCS)){
    colnames(myData)[1] <- "TPP"
    colnames(myData)[8] <- "BU"
    colnames(myData)[14] <- "Nota"
  }

  return(myData)
}



count_fun_teste <- function(inData, inYR, inSE, inGE){

  # inData = countDF
  # inYR = c(21,22)
  # inSE = c("WN")
  # inGE = c("SZD5031VIP3")


  countDF <- aggregate(Nota ~ genotipo + BU + trait + year + season, data = inData, FUN = length)
  dataPointsTable<- count_fun(tableTraits,allData,"local") # arquivo parseCS.R


  if(inGE == "all"){
    inData <- inData %>% filter(year %in% inYR, season %in% inSE)
  }else{
    inData <- inData %>% filter(year %in% inYR, season %in% inSE, genotipo %in% inGE)
  }

  inDF1<- data.frame(tapply(inData$BU, list(inData$BU,inData$trait),length))

  inDF2<- data.frame(tapply(inData$Nota, inData$trait, length))
  inDF2$trait <- rownames(inDF2)

  inDF3 <- data.frame(colSums(!is.na(inDF1)))
  inDF3$trait <- rownames(inDF3)

  outDF <- left_join(inDF2, inDF3, by=c("trait"="trait"))
  colnames(outDF) <- c("dataPoints","trait","Locais")
  outDF <-outDF %>% relocate(trait, .before = dataPoints)

  return(outDF)
}


count_fun <- function(inTable, inData, inType){
  # inTable = tableTraits
  # inData = allData_filter
  # inType <- "genotipos"
  # Funcao antiga para contar

  if(inType == "local"){
    for(i in 1:length(inTable$codigo)){
      # i=1
      dataPoint <- inData %>% filter(trait == inTable$codigo[i]) %>%
        select(local,genotipo,BU,barcode,rep,trait,Nota)
      dataPoint <- na.omit(dataPoint)
      # dataPoint <- distinct(dataPoint,barcode,trait, .keep_all = TRUE)
      trt <- nrow(dataPoint)
      inTable$dataPoints[i] <- trt
      inTable$BU[i] <- length(unique(dataPoint$BU))
      inTable$experimentos[i] <- length(unique(dataPoint$local))
    }
    # Adicionando uma linha com os totais:
    inTable[nrow(inTable)+1,] <- c("Total","Total",sum(inTable$dataPoints),
                                   sum(inTable$BU),sum(inTable$experimentos))

  }else{
    inGenos <- unique(inData$genotipo)
    inTable$dataPoints <- NA
    inTable$BU <-NA
    inTable$experimentos <- NA
    inTable$genotipo <- NA
    inTable$TPP <- myTPP
    inTable1 = inTable
    inTable2 = inTable
    cc=1

    for(j in 1:length(inGenos)){
      if(j>nrow(inTable1)){
        inTable <- rbind(inTable,inTable1)
      }
      for(i in 1:length(inTable1$codigo)){
        cat(cc, i,"\n")
        dataPoint <- inData %>% filter(genotipo == inGenos[j], trait ==inTable$codigo[i]) %>%
          select(local,genotipo,BU,barcode,rep,trait,Nota)
        dataPoint <- na.omit(dataPoint)
        if(nrow(dataPoint)==0){next}
        # dataPoint <- distinct(dataPoint,barcode,trait, .keep_all = TRUE)
        trt <- nrow(dataPoint)
        inTable$dataPoints[cc] <- trt
        inTable$BU[cc] <- length(unique(dataPoint$BU))
        inTable$experimentos[cc] <- length(unique(dataPoint$local))
        inTable$genotipo[cc] <- inGenos[j]
        cc=cc+1
      }
      inTable2 <- rbind(inTable2,inTable)
      j=j+1
    }

    inTable <- inTable2
    # Adicionando uma linha com os totais:
    inTable[nrow(inTable)+1,] <- c("Total","Total",sum(inTable$dataPoints),
                                   sum(inTable$BU),sum(inTable$experimentos),"-",",")
  }

  return(inTable)
}


dfFUN <- function( inData, crit){
  # inData = allData_filter
  # crit = "Yes"

  crit = tolower(crit)
  if(crit == "yes"){
    outData <- inData %>% filter(trait %in% selTrait, BU %in% locHy) %>%
      select(BU,local,barcode,stg,genotipo,rep,trait,Nota,TPP, year, season)
  }else{
    outData <- inData %>% filter(trait %in% filterTrait$codigo[1]) %>%
      select(BU,local,barcode,stg,genotipo,rep,trait,Nota,TPP, year, season)

  }

  return(outData)

}
