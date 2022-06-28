

library(stringr)
read_fun <- function(inLista){
  # inLista = filenames
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
                          range=NA,
                          row=NA,
                          BU=NA,
                          stg=NA,
                          season=NA,
                          barcode=NA,
                          year=NA
  )

  for(i in 1:length(inLista)){
    # i=1
    inFile <- inLista[i]
    # inFile <- gsub(pattern = "../",replacement = "",x = inFile,fixed = T)

    checkError <- FALSE
    tryCatch( { readData <- read.csv(file=inFile, header=T, as.is=1, sep=",") }
              , warning = function(w) { cat("Arquivo \"",inFile,"\" nao encontrado.\n"); checkError <- TRUE })

    if(checkError == T){stop("Arquivo nao encontrado, checar o local do arquivo.")}

    readData$TPP <- NA
    nomesCol <- colnames(readData)
    if(length(nomesCol)==1){ readData <- read.csv(file=inFile, header=T, as.is=1, sep=";")}

    selNames <- c("TPP","EXTID","ABBRC","PLTID","REPNO","COSTR","HELMR","GRLSR","LFSPR","SCLBR","SRSTR","DLLFR","STRNG","STROW","PSS.BARCD")
    checkNames <- match("TRUE",str_detect(nomesCol,paste(selNames[2:5],collapse = "|")))
    if(is.na(checkNames)==T){stop("Stop! Esta faltando umas das colunas principais, checar arquivo.")}

    newData <- readData

    count = 1
    for(i in 1:length(selNames)){
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
                           "GRLSR","LFSPR","SCLBR","SRSTR","DLLFR",
                           "range","row","barcode")


    # capturando as TPPs
    newData$TPP <- stri_sub(newData$local, -4,-3)
    newData$TPP <- gsub("BL","TPP09",newData$TPP)
    newData$TPP <- gsub("BU","TPP11",newData$TPP)
    newData$TPP <- gsub("BT","TPP10",newData$TPP)
    newData$TPP <- gsub("BC","TPP12",newData$TPP)


    newData$BU <- stri_sub(newData$local, -4) # capturando as BUs

    newData$stg <- as.numeric(stri_sub(newData$local, -7,-7))
    newData$stg <- gsub(8,6,newData$stg) # transformando estagio 8 em 6
    newData$stg <- gsub("L",4,newData$stg)

    newData$season <- substring(newData$local,3,4)
    newData$year <- substring(newData$local,1,2)
    newData$year <- as.numeric(newData$year)

    finalData <- rbind(finalData,newData)

  }
  finalData <- finalData[-1,]
  finalData <- tidyr::gather(finalData, trait, Nota, COSTR:DLLFR)
  finalData$Nota <- gsub(",", "", finalData$Nota)
  # finalData$Nota <- as.numeric(finalData$Nota)
  finalData <- distinct(finalData, barcode, trait,rep,.keep_all = T)

  return(finalData)

}

#filtrando para locais que possuem dados de corn stunt
filter_locais <- function(inLoc, stg, yrs, crt5, minEnter, maxEnter,inSea){
  # inLoc <- myDF
  # stg <- Stage
  # yrs <- as.numeric(years)
  # crt5 <- filterYear
  # minEnter <- 50
  # maxEnter <- 200
  # inSea <- "Winter"

  inLoc <- inLoc %>% filter(Season==inSea) %>% filter(Year==selYear) %>% filter(Stage == stg)

  unique(inLoc$local)

  entLoc <- as.data.frame(tapply(inLoc$genotipo, inLoc$local,length))
  entLoc <- entLoc %>% mutate(local = rownames(entLoc))
  colnames(entLoc)[1] <- "Entradas"
  entLoc <- entLoc %>% filter(Entradas > minEnter) %>% filter(Entradas < maxEnter) %>% arrange(Entradas)

  vLoc <- entLoc$local

  y1 <- max(yrs)
  y2 <- yrs[2]
  y3 <- min(yrs)

  if(crt5 == "yes"){
    stg2 <- stg-1
    inLoc$Stage <- as.numeric(inLoc$Stage)
    outLoc1 <- inLoc[which(inLoc$local%in%vLoc),] %>% filter(Year== y1)
    outLoc2 <- inLoc[which(inLoc$local%in%vLoc),] %>% filter(Year== y2)
    outLoc3 <- inLoc[which(inLoc$local%in%vLoc),] %>% filter(Year== y3)

    outLoc <- rbind(outLoc1, outLoc2, outLoc3)
  }else{
    outLoc <- inLoc %>% filter(local %in% vLoc)
  }

  return(outLoc)
}

sinoFun <- function(myData, sinIN){
  # sinIN <- sinN
  # myData <- checksF

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

  if(ncol(myData)>8){
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
  return(myData)
}


parsing_data <- function(stageData, nomes_entrada){
  # Criando um vetor para guardar o numero da coluna
  index_nomes <- c(1:length(nomes_selecionados))

  for (i in 1:length(nomes_selecionados)){
     index_nomes[i] <- which( colnames(stageData)==nomes_entrada[i] )
     colnames(stageData)[index_nomes[i]] <- novos_nomes[i]
   }

  stageData$genotipo <- gsub(" ", "", stageData$genotipo, fixed = TRUE)
  stageData$local <- gsub(" ", "", stageData$local, fixed = TRUE)

  stageData$genotipo <- toupper(stageData$genotipo)
  stageData$local <- toupper(stageData$local)

  stageOut <- stageData %>% select(year, local, plotName, Rep, genotipo, Nota, yield)
  filter_locais(stageOut)
}


count_fun <- function(inTable,inData, inType){
  # inTable = tableTraits
  # inData = myCountData
  # inType <- "genotipos"

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
  crit = "Yes"

  crit = tolower(crit)
  if(crit == "yes"){
    outData <- inData %>% filter(trait %in% filterTrait$codigo[1], stg %in% myStag, local %in% locHy) %>%
      select(BU,local,barcode,stg,genotipo,rep,trait,Nota,TPP, year, season)
  }else{
    outData <- inData %>% filter(trait %in% filterTrait$codigo[1], stg %in% myStag) %>%
      select(BU,local,barcode,stg,genotipo,rep,trait,Nota,TPP, year, season)

  }
  return(outData)

}
