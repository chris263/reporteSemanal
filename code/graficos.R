
#removing warnings
options(warn = -1)

# Script para gerar os graficos
require(grid)
require(hrbrthemes)
require(viridis)
require(ggplot2)

# install.packages("geobr")
# library(brazilmaps)
# require(geobr)
library(sf)
# library(maptools)
# library(leaflet)

# theme_set(theme_bw())

# mapa <- brazilmaps::get_brmap("State")

# mapa <- read_state(showProgress = FALSE)

# class(mapa)


makeG0 <- function(inGraph0, fonte, item){
  inGraph0 %>%
    ggplot( aes(x=get(fonte), y=get(item), fill=local)) +
    geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
    theme_ipsum() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    ggtitle(unique(inGraph0$local)) +
    ylim(0,10) +
    xlab(fonte) +
    ylab(item)

}

makeH0 <- function(inGrapH0, termo){
  H0 <- inGrapH0 %>%
    ggplot( aes(x=get(termo))) +
    geom_histogram( binwidth=3, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
    # ggtitle("") +
    # theme_ipsum() +
    theme(
      plot.title = element_text(size=15)
    ) + xlab(termo)

  return(H0)
}


makeG1 <- function(inGraph1){
  g.mid<-ggplot(inGraph1,aes(x=1,y=genotipo))+geom_text(aes(label=genotipo))+
    geom_segment(aes(x=0.94,xend=0.96,yend=genotipo))+
    geom_segment(aes(x=1.04,xend=1.065,yend=genotipo))+
    ggtitle("")+
    ylab(NULL)+
    scale_x_continuous(expand=c(0,0),limits=c(0.94,1.065))+
    theme(axis.title=element_blank(),
          panel.grid=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.background=element_blank(),
          axis.text.x=element_text(color=NA),
          axis.ticks.x=element_line(color=NA),
          plot.margin = unit(c(1,-1,1,-1), "mm"))

  g1 <- ggplot(data = inGraph1, aes(x = genotipo, y = blupNota)) +
    geom_bar(stat = "identity") + ggtitle("Nota Enfezamento") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin = unit(c(1,-1,1,0), "mm")) +
    scale_y_reverse() + coord_flip()

  g2 <- ggplot(data = inGraph1, aes(x = genotipo, y = blupYield)) +xlab(NULL)+
    geom_bar(stat = "identity") + ggtitle("Yield") +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          plot.margin = unit(c(1,0,1,-1), "mm")) +
    coord_flip()
  library(gridExtra)
  gg1 <- ggplot_gtable(ggplot_build(g1))
  gg2 <- ggplot_gtable(ggplot_build(g2))
  gg.mid <- ggplot_gtable(ggplot_build(g.mid))

  grid.arrange(gg1,gg.mid,gg2,ncol=3,widths=c(4/9,2/9,4/9))
}

makeG2 <- function(inGraph2){
  ggplot(inGraph2, aes(x = reorder(genotipo,-Indice))) +
    xlab("Genotype")+
    geom_col(aes( y = Indice, fill="redfill")) +
    geom_text(aes(y = Indice, label = round(Indice,2)), fontface = "bold", vjust = 1.4, color = "black", size = 3) +
    geom_line(aes(y = blupNota * 1.5, group = 1, color = 'blackline')) +
    geom_text(aes(y = blupNota * 1.5, label = round(blupNota, 0)), vjust = 1.4, color = "black", size = 3) +
    scale_y_continuous(sec.axis = sec_axis(trans = ~ . / 1.5)) +
    scale_fill_manual('', labels = 'Indice', values = "99FF66") +
    scale_color_manual('', labels = 'Nota', values = 'black') +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

}

makeG3 <- function(inGraph3){

  inGraph3$Tnota <- 10-inGraph3$blupNota

  bp <-  reshape2::melt(inGraph3, id.vars=c("genotipo"))

  bpf <- bp %>% filter(variable != "Acuracia(%)",
                       variable != "nLocais",
                       variable != "Indice",
                       variable != "blupNota")

  ggplot(bpf, aes(x=reorder(genotipo,-value, sum),y=value, fill = variable)) +
    geom_bar(stat = "identity")+
    scale_fill_viridis(discrete = T) +
    ggtitle("Agrupando Yield e Inverso da Nota CS") +
    xlab("Genotipos")+
    # theme_ipsum()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

makeG4 <- function(inGraph4){
  ggplot(inGraph4, aes(x=reorder(genotipo,blupNota, sum),y=blupNota, fill = genotipo)) +
    geom_bar(stat = "identity")+
    scale_fill_viridis(discrete = T) +
    # ggtitle("Agrupando Yield e Inverso da Nota CS") +
    xlab("Genotipos")+
    # theme_ipsum()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

makeG5 <- function(ctl, hb1, inGraph5){

  ctlHB <- append(ctl,hb1)
  inGraph5$genotipo <- as.character(inGraph5$genotipo)

  cd <- inGraph5[1,]
  for( i in 1: length(ctlHB)){
    cddf <-filter(resultadoFINAL,genotipo == ctlHB[i])
    cd <- rbind(cd,cddf)
  }
  cd <- cd[-1,]

  cdf <-  reshape2::melt(cd, id.vars=c("genotipo"))

  cdf <- cdf %>% filter(variable == "blupNota")

  ggplot(cdf %>% arrange(variable, desc(value)) %>%
           mutate(genotipo=factor(genotipo, levels=genotipo)),
         aes(x=genotipo,y=value, fill = variable)) +
    geom_bar(stat="identity", show.legend=FALSE) +
    geom_text(aes(label=round(value,2), y=0.5*value), colour="white", size=3) +
    facet_grid(. ~ variable, scales="free_x", space="free_x") +
    scale_y_continuous(limits=c(-0.005, 1.05*max(cdf$value)), expand=c(0,0)) +
    theme_classic() +
    theme(panel.spacing=unit(0,"pt"),
          panel.border=element_rect(colour="grey50", fill=NA),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

}

makeG6 <- function(inGraph6){
  ggplot(inGraph6, aes(x = genotipo)) +
    geom_col(aes( y = blupYield, fill="redfill")) +
    geom_text(aes(y =blupYield, label = blupYield), fontface = "bold", vjust = 1.4, color = "black", size = 3) +
    geom_line(aes(y = blupNota * 1.5, group = 1, color = 'blackline')) +
    geom_text(aes(y = blupNota * 1.5, label = round(blupNota, 0)), vjust = 1.4, color = "black", size = 3) +
    scale_y_continuous(sec.axis = sec_axis(trans = ~ . / 1.5)) +
    scale_fill_manual('', labels = 'Yield', values = "darkgreen") +
    scale_color_manual('', labels = 'Nota', values = 'black') +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

}

makeG7 <- function(inGraph7){
  inGraph7 <- altaPP
  head(inGraph7)

  colnames(inGraph7)[3]<-"Nota"
  colnames(inGraph7)[4]<-"yield"
  inGraph7 <-inGraph7 %>% filter(inGraph7$yield>2000)

  ggplot(inGraph7,aes(x=Nota,y=yield)) +
    geom_point() +
    scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
    geom_smooth(method="lm")

}


makeG8 <- function(controls, candis, accHB){

  # controls <- controles
  # candis <- candidato
  # accHB <- accH1


  # Filtrando controles
  card1 <- filter(resultadoFINAL, genotipo %in% controls) %>% arrange(Indice)
  card1$correct <- card1$Indice + (-1)*min(card1$Indice)
  card1$Percent <- card1$correct/max(card1$correct)*100

  # Filtrando candidato
  cand1 <- filter(resultadoFINAL, genotipo %in% candis)
  cand1$correct <- cand1$Indice + (-1)*min(card1$Indice)
  cand1$Percent <- cand1$correct/max(card1$correct)*100
  cH <- cand1$Percent

  cardT <<- rbind(card1,cand1)

  # brks <- card1$Percent
  # nomesG <- card1$genotipo
  # nomeH <- cand1$genotipo

  # Padrao
  gg.gauge <- function(pos,breaks=card1$Percent) {
    require(ggplot2)
    get.poly <- function(a,b,r1=0.5,r2=1.0) {
      th.start <- pi*(1-a/100)
      th.end   <- pi*(1-b/100)
      th       <- seq(th.start,th.end,length=100)
      x        <- c(r1*cos(th),rev(r2*cos(th)))
      y        <- c(r1*sin(th),rev(r2*sin(th)))
      return(data.frame(x,y))
    }
    ggplot()+
      ggtitle(paste0("Certeza: ",accHB, "%"))+
      geom_polygon(data=get.poly(breaks[1],breaks[2]),aes(x,y),fill="#EB5C5F")+
      geom_polygon(data=get.poly(breaks[2],breaks[3]),aes(x,y),fill="#FA9594")+
      geom_polygon(data=get.poly(breaks[3],breaks[4]),aes(x,y),fill="#FBF2EA")+
      geom_polygon(data=get.poly(breaks[4],breaks[5]),aes(x,y),fill="#F2E599")+
      geom_polygon(data=get.poly(breaks[5],breaks[6]),aes(x,y),fill="#ECD969")+
      geom_polygon(data=get.poly(breaks[6],breaks[7]),aes(x,y),fill="#99D388")+
      geom_polygon(data=get.poly(breaks[7],breaks[8]),aes(x,y),fill="forestgreen")+
      geom_polygon(data=get.poly(pos-1,pos+1,0.2),aes(x,y))+
      geom_text(data=as.data.frame(breaks), size=2, fontface="bold", vjust=0,
                aes(x=1.1*cos(pi*(1-breaks/100)),y=1.1*sin(pi*(1-breaks/100)),label=card1$genotipo))+
      annotate("text",x=0,y=0,label=cand1$genotipo,vjust=0,size=4,fontface="bold")+
      coord_fixed()+
      theme_bw()+
      theme(axis.text=element_blank(),
            axis.title=element_blank(),
            axis.ticks=element_blank(),
            panel.grid=element_blank(),
            panel.border=element_blank())
  }
  gg.gauge(cH,breaks=card1$Percent)

}
#round(maxAcc,2)
makeS1 <- function(inS1){
  ggplot(inS1, aes(x = qualidade)) +
    geom_col(aes( y = avgAcc, fill="redfill")) +
    ylab("Acuracia")+
    geom_text(aes(y = avgAcc, label = round(avgAcc,2)), fontface = "bold", vjust = 1.4, color = "black", size = 3) +
    geom_line(aes(y = maxAcc * 1.5, group = 1, color = 'blackline')) +
    geom_text(aes(y = maxAcc * 1.5, label = maxAcc), vjust = 1.4, color = "black", size = 3) +
    geom_line(aes(y = nLocaisAlta), group =2, color = "darkred") +
    geom_text(aes(y = nLocaisAlta * 1.5, label = nLocaisAlta), vjust = 1.4, color = "red", size = 3) +
    scale_y_continuous(sec.axis = sec_axis(trans = ~ . / 1.5)) +
    scale_fill_manual("", labels = "Qualidade Media", values = "#82D1F1") +
    scale_color_manual("",values = c("Qualidade Maxima"="black","Numero Locais"="red")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

}

makeG9 <- function(inDados,grupos){
  # inDados = myGraph
  # grupos = c("MS","MT","Candidato")

  inTPP <- unique(inDados$TPP)

  cand1 <- data.frame(genotipo=myHB,Check=1,Trait = filterTrait$nome[1],Esperado=NA,Tipo="Candidato",TPP=inTPP)
  checksS <- checksF %>% filter(Tipo %in% grupos, TPP == inTPP)
  checksS <- rbind(checksS, cand1)

  inDados <- inDados %>% filter(genotipo %in% checksS$genotipo)
  try(inDados <- left_join(inDados,checksS,x.by="genotipo", y.by="genotipo"),silent = T)
  inDados$Nota <- as.numeric(inDados$Nota)

  ggplot(inDados, aes(x=genotipo, y=Nota, fill=Tipo)) +
    geom_boxplot() +
    ggtitle(paste0("Local:", unique(inDados$BU)))+#," - PLC", unique(inDados$stg))) +
    scale_y_continuous(breaks=seq(0.0, 10, 1), limits=c(0, 10))+
    geom_hline(yintercept=2, linetype="dashed", color = "blue", size=1.2) +
    geom_hline(yintercept=4, linetype="dashed", color = "orange", size=1.2) +
    geom_hline(yintercept=6, linetype="dashed", color = "red", size=1.2) +

    coord_flip()
}

makeHM <- function(inHM, locs){
  # inHM = myDF
  # locs = locAlta

  myDF2 <- inHM %>% filter(local %in% locs)
  data <- data.frame(X=myDF2$range,Y=myDF2$row,Nota=myDF2$Nota, genotipo=myDF2$genotipo)

  ggplot(data, aes(X, Y, fill= Nota)) +
    geom_tile() +
    scale_fill_viridis(discrete=FALSE) +
    theme_ipsum()

}


makePie <- function(inTable, inStage){
  # inTable = tableCount
  # inStage = 6

  inTable$dataPoints <- as.numeric(inTable$dataPoints)
  inTable <- inTable[!inTable$nome %in% "Total", ]

  # Add label position
  inTable <- inTable %>%
    arrange(desc(nome)) %>%
    mutate(lab.ypos = cumsum(dataPoints) - 0.5*dataPoints)


  mycols <- c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF",
              "#99D1E9", "#9391C9","#EEB4A7")

  ggplot(inTable, aes(x = "", y = dataPoints, fill = nome)) +
    ggtitle(inStage) +
    geom_bar(width = 1, stat = "identity", color = "white") +
    coord_polar("y", start = 0)+
    geom_text(aes(y = lab.ypos, label = dataPoints), color = "white")+
    scale_fill_manual(values = mycols) +
    theme_void()


}



makeBar <- function(){
  library(ggplot2)

  # create dummy data
  data <- data.frame(
    name=letters[1:5],
    value=sample(seq(4,15),5),
    sd=c(1,0.2,3,2,4)
  )



  # Most basic error bar
  ggplot(dataFinal) +
    geom_bar( aes(x=reorder(genotipo,desc(blupNota)), y=blupNota), stat="identity", fill="skyblue", alpha=0.7) +
    geom_errorbar( aes(x=genotipo, ymin=blupNota-sdErr, ymax=blupNota+sdErr), width=0.4, colour="orange", alpha=0.9, size=1.3)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

makePlots <- function(inLocDF, inStg){
  myColor <- c("#51BD80","#574C8F","#FEA287")
  inLocDF %>% mutate(local=fct_reorder(local, desc(Entradas))) %>%
    ggplot( aes(x=local, y=Entradas)) +
    ggtitle(paste0("Stage: ",inStg)) +
    geom_bar(stat="identity", fill=myColor[inStg-3], alpha=.6, width=.4) + coord_flip() +
    xlab("Locais") + ylab("Numero de Entradas") +
    theme_bw()
}




