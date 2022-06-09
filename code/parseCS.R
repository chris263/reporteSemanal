

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
  require(stringr)
  # sinIN <- sinN
  # myData <- myDF

  # sinIN$genotipo <- gsub(" ","",sinIN$genotipo)
  sinIN$Novo <- gsub(" ","",sinIN$Novo)
  sinIN$Novo <- gsub("\\\\","_",sinIN$Novo)
  sinIN$genotipo <- gsub(" ","",sinIN$genotipo)
  sinIN$genotipo <- gsub("\\\\","_",sinIN$genotipo)
  myData$genotipo <- gsub(" ","", myData$genotipo)
  myData$genotipo <- gsub("\\\\","_", myData$genotipo)

  for( i in 1:nrow(sinIN)){
    # i=1
    # print(i)
    myOld <- sinIN$genotipo[i]
    myNew <- sinIN %>% filter(genotipo == myOld) %>% dplyr::select(Novo)
    pp <-str_c("\\b", myOld, "\\b", collapse="|")
    myData$genotipo<-stringr::str_replace_all(myData$genotipo, pp, myNew$Novo)

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
