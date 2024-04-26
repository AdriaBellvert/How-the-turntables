library(rgl)
library(gMOIP)
library(cxhull)
library(plyr)
library(dplyr)
library(dismo)
library(gam)

load("data_for_climatic_range.bin")
model <- data_for_climatic_range$model
predictors <- data_for_climatic_range$predictors
ll <- data_for_climatic_range[[3]]

bkg <- list()
is <- c("g", "p", "h")

for (i in 1:length(is)) {
  
  rc <- crop(predictors, ll[[i]])	
  
  background <- randomPoints(mask=rc, n=5000)
  background<-data.frame(background)
  background.variables<-data.frame(extract(predictors, background))
  background<-cbind(background, background.variables)
  bkg[[i]] <- background
  
}

names(bkg)<- is

sil_bio12 <-model[["silvatica_P"]][["data"]][["bio12"]][1:9]
sil_bio5 <- model[["silvatica_P"]][["data"]][["bio5"]][1:9]
sil_elev <- model[["silvatica_P"]][["data"]][["elev"]][1:9]
sil_p <- cbind(sil_bio12,sil_bio5,sil_elev)
sil_p <- as.data.frame(sil_p)
sil_p$sp <- "silvatica_P"

sil_bio12 <-model[["silvatica_G"]][["data"]][["bio12"]][1:27]
sil_bio5 <- model[["silvatica_G"]][["data"]][["bio5"]][1:27]
sil_elev <- model[["silvatica_G"]][["data"]][["elev"]][1:27]
sil_g <- cbind(sil_bio12,sil_bio5,sil_elev)
sil_g <- as.data.frame(sil_g)
sil_g$sp <- "silvatica_G"

sil_bio12 <-model[["silvatica_H"]][["data"]][["bio12"]][1:7]
sil_bio5 <- model[["silvatica_H"]][["data"]][["bio5"]][1:7]
sil_elev <- model[["silvatica_H"]][["data"]][["elev"]][1:7]
sil_h <- cbind(sil_bio12,sil_bio5,sil_elev)
sil_h <- as.data.frame(sil_h)
sil_h$sp <- "silvatica_H"


sil <- rbind(sil_g,sil_h,sil_p)




gg12<- data.frame(bkg[["g"]][["bio12"]], "g")
hh12<- data.frame(bkg[["h"]][["bio12"]], "h")
pp12<- data.frame(bkg[["p"]][["bio12"]], "p")
colnames(gg12) <- c("sil_bio12", "sp")
colnames(hh12) <- c("sil_bio12", "sp")
colnames(pp12) <- c("sil_bio12", "sp")
T12 <- rbind(gg12,hh12,pp12)

gg12<- data.frame(bkg[["g"]][["bio5"]], "g")
hh12<- data.frame(bkg[["h"]][["bio5"]], "h")
pp12<- data.frame(bkg[["p"]][["bio5"]], "p")
colnames(gg12) <- c("sil_bio5", "sp")
colnames(hh12) <- c("sil_bio5", "sp")
colnames(pp12) <- c("sil_bio5", "sp")
T5 <- rbind(gg12,hh12,pp12)

gg12<- data.frame(bkg[["g"]][["elev"]], "g")
hh12<- data.frame(bkg[["h"]][["elev"]], "h")
pp12<- data.frame(bkg[["p"]][["elev"]], "p")
colnames(gg12) <- c("elev", "sp")
colnames(hh12) <- c("elev", "sp")
colnames(pp12) <- c("elev", "sp")
Te <- rbind(gg12,hh12,pp12)

TT <- cbind(T12$sil_bio12, T5$sil_bio5, Te$elev)
TT<- as.data.frame(TT)
TT$sp <- T12$sp
colnames(TT) <- colnames(sil)
silTT <- rbind(sil, TT)
target <- c("g", "silvatica_G", "h", "silvatica_H", "p", "silvatica_P")
silTT <- silTT %>% arrange(factor(sp, levels = target))
silTT <- na.omit(silTT)


sil_bio12 <-model[["gomerensis_G"]][["data"]][["bio12"]][1:41]
sil_bio5 <- model[["gomerensis_G"]][["data"]][["bio5"]][1:41]
sil_elev <- model[["gomerensis_G"]][["data"]][["elev"]][1:41]
sil_g <- cbind(sil_bio12,sil_bio5,sil_elev)
sil_g <- as.data.frame(sil_g)
sil_g$sp <- "gomerensis_G"

sil_bio12 <-model[["gomerensis_H"]][["data"]][["bio12"]][1:24]
sil_bio5 <- model[["gomerensis_H"]][["data"]][["bio5"]][1:24]
sil_elev <- model[["gomerensis_H"]][["data"]][["elev"]][1:24]
sil_h <- cbind(sil_bio12,sil_bio5,sil_elev)
sil_h <- as.data.frame(sil_h)
sil_h$sp <- "gomerensis_H"


gom <- rbind(sil_g,sil_h)
gomTT <- rbind(gom, TT)
target <- c("g", "gomerensis_G", "h", "gomerensis_H")
gomTT <- gomTT %>% arrange(factor(sp, levels = target))



sil_bio12 <-model[["calderensis_P"]][["data"]][["bio12"]][1:9]
sil_bio5 <- model[["calderensis_P"]][["data"]][["bio5"]][1:9]
sil_elev <- model[["calderensis_P"]][["data"]][["elev"]][1:9]
sil_p <- cbind(sil_bio12,sil_bio5,sil_elev)
sil_p <- as.data.frame(sil_p)
sil_p$sp <- "calderensis_P"

sil_bio12 <-model[["calderensis_G"]][["data"]][["bio12"]][1:15]
sil_bio5 <- model[["calderensis_G"]][["data"]][["bio5"]][1:15]
sil_elev <- model[["calderensis_G"]][["data"]][["elev"]][1:15]
sil_g <- cbind(sil_bio12,sil_bio5,sil_elev)
sil_g <- as.data.frame(sil_g)
sil_g$sp <- "calderensis_G"




cal <- rbind(sil_g,sil_p)
calTT <- rbind(cal, TT)
target <- c("g", "calderensis_G", "p", "calderensis_P")
calTT <- calTT %>% arrange(factor(sp, levels = target))

spTT <- rbind(sil, cal, gom, TT)

sptti <- as.data.frame(scale(spTT[,1:3]))
sptti$V4 <- spTT$sp
spTT <- sptti
colnames(spTT) <- c("sil_bio12", "sil_bio5",  "sil_elev",  "sp")



sppi <- unique(spTT$sp)
hul_ar2 <- matrix(nrow = 1, ncol = 10)
for (i in 1:length(sppi)) {
  subdata <- spTT[spTT$sp%in%sppi[i],]
  subdata = subset(spTT, sp == sppi[i])
  trait <- subdata[,-4]
  vertices <- na.omit(as.matrix(trait))
  vertices <- unique(vertices[,c(1:3)])
  vrtx <- convexHull(vertices)
  vrtx <- vrtx$pts
  vrtx <- vrtx[vrtx$vtx%in%"TRUE",]
  vrtx2 <- as.matrix(vrtx[,1:3])
  hull <- cxhull(vrtx2, triangulate = TRUE)
  hul_ar2[,i] <- sum(sapply(hull$facets, `[[`, "volume"))
  
}

colnames(hul_ar2) <- sppi

#z-test

res_cal <- prop.test(x = c(hul_ar2[,"calderensis_P"]/hul_ar2[,"p"]*100,
                           hul_ar2[,"calderensis_G"]/hul_ar2[,"g"]*100), n = c(100, 100), alternative = "less")
res_cal


res_gom <- prop.test(x = c(hul_ar2[,"gomerensis_H"]/hul_ar2[,"h"]*100,
                           hul_ar2[,"gomerensis_G"]/hul_ar2[,"g"]*100), n = c(100, 100), alternative = "less")
res_gom



res_sil_p <- prop.test(x = c(hul_ar2[,"silvatica_P"]/hul_ar2[,"p"]*100,
                             hul_ar2[,"silvatica_G"]/hul_ar2[,"g"]*100), n = c(100, 100), alternative = "less")
res_sil_p


res_sil_h <- prop.test(x = c(hul_ar2[,"silvatica_H"]/hul_ar2[,"h"]*100,
                             hul_ar2[,"silvatica_G"]/hul_ar2[,"g"]*100), n = c(100, 100), alternative = "less")
res_sil_h

