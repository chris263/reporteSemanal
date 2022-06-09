a <- c(1:10)
b <- c(11:20)
v <- c(10,-5)
A <- matrix(a,b,nrow = 10, ncol=2, byrow = FALSE)

xy <- data.frame(a,b)
xy <- as.matrix(xy)

R <- xy %*% v

lineaR <- lm(yield~Nota, data=altaP_DFF)
plot(yield~Nota, data = altaP_DFF)
abline(lineaR)

# y = bx+a
a <- coef(lineaR)[1]
b <- coef(lineaR)[2]

pYield1 <- b*3+a
pYield2 <- b*4+a
resDiff<- round(pYield1 - pYield2,2)
names(resDiff)<-c("yield")

require(ggplot2)
ggplot(altaP_DFF,aes(x=Nota,y=yield)) +
  geom_point() +
  geom_smooth(method="lm")

teste <- altaP_DFF %>% filter(Nota <= 2, yield < 8)
# Capturando dados errados, baixa nota e baixíssimo Yield
dfOut <- inFile %>% filter(Nota < 2, yield < 2)
listaOut <- dfOut$genotipo


ccs <- append(altaP_DFF$genotipo,altaP_DFF$local)


# Teste de modelos

if(duasT == "yes"){
  modeloConjuntoA <- mmer(cbind(yield, Nota) ~ local,
                          random= ~ vs(genotipo, Gtc=unsm(2)) + vs(localGeno, Gtc=unsm(2)),
                          rcov= ~ vs(units, Gtc=unsm(2)),
                          data=altaP_DFF, verbose = FALSE)
}else{
  modeloConjunto <- mmer(cNota ~ local,
                         random =~ genotipo+vs(ds(local),genotipo),
                         rcov =~ vs(ds(local), units),
                         data = altaP_DFF, verbose = FALSE)
}


modeloConjuntoA$Vi


summary(modeloConjunto)

modeloConjunto$U$genotipo$Nota
modeloConjunto$fitted

modeloConjuntoA$U$`u:genotipo`$Nota <- as.data.frame(modeloConjuntoA$U$`u:genotipo`$Nota)
rownames(modeloConjuntoA$U$`u:genotipo`$Nota) <- gsub("id","",rownames(modeloConjuntoA$U$`u:genotipo`$Nota))

tabela <- modeloConjuntoA$U$`u:genotipo`$Nota


cc1 <- summary(modeloConjuntoA)

cc2 <- as.data.frame(cc1$betas)
cc3<-cc2[2,3]
tabelaPQP <- tabela[,1]+cc3


# Segunda tentativa
modeloConjuntoB <- mmer(cbind(yield, Nota) ~ genotipo + local,
                        random= ~ vs(genotipo, Gtc=unsm(2)) + vs(localGeno, Gtc=unsm(2)),
                        rcov= ~ vs(units, Gtc=unsm(2)),
                        data=altaP_DFF, verbose = FALSE)

modB <- summary(modeloConjuntoB)

dfG1 <- as.data.frame(modeloConjuntoB$U$`u:genotipo`$Nota)
dfG1$nomes <- rownames(dfG1)
nomesG1 <- dfG1$nomes




locG <-  as.data.frame(modB$betas)
locG1 <- locG %>% filter(Trait == "Nota")
nomesEx <- locG1$Effect
interN <- nomesG1[!nomesG1%in%nomesEx]
inter <- locG1$Estimate[1]
locG1$Adj <- locG1$Estimate + inter
locG1$Adj[1] <- inter
locG1$Effect <- gsub("genotipo","",locG1$Effect)
locG1$Effect[1] <- interN
# locG1$Estimate[2:nrow(locG1)] <- locG1$Estimate[2:nrow(locG1)]+locG1$Estimate[1]

locGF <- locG1[1:length(nomesG1),]
locFF <- locG1[(length(nomesG1)+1):nrow(locG1),]

conjuntoC <- mmer(cbind(yield,Nota)~local,
                  random= ~genotipo + vs(ds(localGeno),genotipo),
                  rcov= ~ vs(ds(local),units),
                  data=altaP_DFF, verbose = FALSE)

conjuntoD <- mmer(cbind(yield,Nota)~local,
                  random= ~genotipo + vs(ds(local),genotipo),
                  rcov= ~ vs(ds(genotipo),units),
                  data=altaP_DFF, verbose = FALSE)

 anova(modeloConjuntoA, conjuntoC)

unsBLUP(modeloConjuntoA$U[1:6])

modeloConjuntoA$U$`u:genotipo`$Nota
mm <- summary(modeloConjuntoA)
mm$varcomp

p0 <- predict.mmer(object=modeloConjuntoA, classify = "Nota")
