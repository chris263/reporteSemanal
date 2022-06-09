#Calculando a estabilidade


# Preparando o arquivo
locN <- unique(dadosGenoLoc$local)
ccs3 <- dadosGenoLoc %>% dplyr::select(local, genotipo, Nota) %>% filter(local=="BU24")
ccs3<-as.data.frame(tapply(ccs3$Nota, ccs3$genotipo, mean, na.rm=T))
colnames(ccs3)[1]<-"Delete"
for (i in 1:length(locN)){
  ccs1 <- dadosGenoLoc %>% dplyr::select(local, genotipo, yield) %>% filter(local==locN[i])
  ccs2<-as.data.frame(tapply(ccs1$Nota, ccs1$genotipo, mean, na.rm=T))
  colnames(ccs2)[1] <- locN[i]
  ccs3 <- cbind(ccs3,ccs2)
}
ccs4 <- as.data.frame(ccs3[,-1], na.rm=T)
ccs4 <- as.matrix(round(ccs4,2))

library(Bilinear)
AMMIfit <- bilinear(x = ccs4, verbose=F)
AMMIplot(AMMIfit)
AMMIplot(AMMIfit, PC=2)
