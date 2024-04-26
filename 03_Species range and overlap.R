rm(list = ls())
#dev.off()

setwd("C:/Users/Adri/OneDrive - Universitat de Barcelona/Niche/Laura Pollock/PresentClimate-4AdriaMiquel/4.1.presentClimatologies/100")

library(plyr)
library(dplyr)
library(dismo)
library(gam)

# Load data
data <- read.table("Data_sp_local.txt", header = T, sep = "\t")
spp <- as.vector(unique(data$SpeciesName))

# Coordinate transformation
datacoord <- data[,c(9,8)]
colnames(datacoord) <- c("X", "Y")
coordinates(datacoord) <- c("X", "Y")
proj4string(datacoord) <- CRS("+proj=longlat +datum=WGS84") 
res <- spTransform(datacoord, CRS("+proj=utm +zone=28 ellps=WGS84"))
cord <- as.data.frame(res@coords)

data$LocLongitude <- cord$coords.x1
data$LocLatitude <- cord$coords.x2
data$SpeciesName <- as.factor(data$SpeciesName)

#Load and process raster data
bio12 <- raster("CHELSA_CanaryIslands_bio_12_1979_2013.tif")
elev <- raster("CanaryIslands_DEM_100m.tif")
bio5 <- raster("CHELSA_CanaryIslands_bio_05_1979_2013.tif")

namesb <- c("bio12",	"elev",	"bio5")

df <- list(bio12,	elev,	bio5)
hg <- extent(188197,303808, 3059779, 3202078)
for (i in 1:length(namesb)) {
  rc <- crop(df[[i]], hg)
  df[[i]] <- rc
  
}

vbrick <- brick(df)
names(vbrick) <- c("bio12",	"elev",	"bio5")

predictors <- stack(vbrick[[1:3]])

variables.tabla<-as.data.frame(predictors)


vbrick <- brick(df[[1]], df[[2]], df[[3]])
names(vbrick) <- c("bio12",	"elev",	"bio5")
predictors <- stack(vbrick[[1:3]])
names(predictors) <- c("bio12",	"elev",	"bio5")


# Load island specific extents for La Gomera (g), La Palma (p), El Hierro (h)
# and crop raster layers for each island
g <- extent(255126, 304349, 3090274, 3130417)
p <- extent(188197, 246326, 3143136, 3202078)
h <- extent(188197, 221170, 3059779, 3086736)

ll <- list(g,p,h)
names(ll) <- c("g","p","h")


predictors_island <- list()

for (i in 1:length(ll)) {
  predictors_island[[i]] <- crop(predictors, ll[[i]])	
}
names(predictors_island) <- c("g", "p", "h")


# Initialize lists and matrix
model_predict <- list()
model <- list()
coefficient_sps <- list()
variable_response <- list()

ni_ov <- matrix(nrow = 1000, ncol = 4)
colnames(ni_ov) <- c("sgG", "sgH", "scG", "scP")
art <- matrix(nrow = 1000, ncol = 14)
colnames(art) <- c("sg0", "sh0", "sp0", "sg1", "sh1", "sp1", "gg0", "gh0", "gg1", "gh1", "cg0", "cp0", "cg1", "cp1")


#Link each species/group to their islands
spisl <- as.data.frame(matrix(c (spp, c("p", "h", "g", "h", "g", "p", "g")), ncol=2))


#GAM model 
for (k in 1:1000) {
  
for (i in 1:length(spp)) {
  datas <- unique(data[,c(8,9,12)])
  presencia_m <- datas[datas$SpeciesName%in%spp[i],]
  presencia_m$SpeciesName <- droplevels(presencia_m$SpeciesName)
  presencia_m$SpeciesName <- 1
  
  
  x<- as.character(spisl[i,2])
  rc <- crop(predictors, ll[[x]])	
  
  background <- randomPoints(mask=rc, n=500)
  
  background<-data.frame(background)
  background.variables<-data.frame(extract(rc, background))
  background<-cbind(background, background.variables)
  background$SpeciesName<- 0 
  pre<-presencia_m[,-3]
  pre<- pre[,c("LocLongitude", "LocLatitude")]
  presencia.variables<-data.frame(extract(rc, pre))
  presenciaT<-cbind(presencia_m, presencia.variables)
  colnames(background) <- c("LocLongitude", "LocLatitude","bio12",	"elev",	"bio5","SpeciesName")
  
  presencia.background<-rbind(presenciaT, background)
  formula.gam<-as.formula(paste("SpeciesName ~ s(", paste(names(rc), collapse=") + s("), ")", collapse=""))
  m.gam.temp<-gam(formula.gam, family=binomial(link=logit), data=presencia.background, select=TRUE)
  
  
  m.gam.temp.mapa<-predict(rc, m.gam.temp, type="response")
  names(m.gam.temp.mapa)<-spp[i]
  
  
  model[[i]] <- m.gam.temp
  model_predict[i] <- m.gam.temp.mapa
  su <- summary(m.gam.temp)
  coefficient_sps[[i]] <- su
  variable_response[[i]] <- m.gam.temp
  
}

  names(model_predict) <- spp

#Calculate niche overlap for each species pair
ni_ov[k,1] <- nicheOverlap(model_predict[[3]], model_predict[[5]]) #???in gomera
ni_ov[k,2] <- nicheOverlap(model_predict[[2]], model_predict[[4]]) #in hierro
ni_ov[k,3] <- nicheOverlap(model_predict[[3]], model_predict[[7]]) #in gomera
ni_ov[k,4] <- nicheOverlap(model_predict[[1]], model_predict[[6]])  #in palma



#binarize model predictions and account the cells that the island is present (1) or absent (0)
#for each species/island

model_predict_bin <- list()

for (j in 1:length(model_predict)) {
  bin <- model_predict[[j]]
  bin[bin@data@values > 0.05] <- 1
  bin[bin@data@values <= 0.05] <- 0
  model_predict_bin[[j]] <- bin
  names(model_predict_bin[[j]]) <- names(model_predict[[j]])
  
}

names(model_predict_bin) <- names(model_predict)


sil_p <- table(model_predict_bin[["silvatica_P"]]@data@values)
sil_g <- table(model_predict_bin[["silvatica_G"]]@data@values)
sil_h <- table(model_predict_bin[["silvatica_H"]]@data@values)
gom_h <- table(model_predict_bin[["gomerensis_H"]]@data@values)
gom_g <- table(model_predict_bin[["gomerensis_G"]]@data@values)
cal_P <- table(model_predict_bin[["calderensis_P"]]@data@values)
cal_G <- table(model_predict_bin[["calderensis_G"]]@data@values)


#Calculate how many times the species is present or absent in each island
art[k,1] <- sil_g[1]
art[k,2] <- sil_h[1] 
art[k,3] <- sil_p[1] 
art[k,7] <- gom_g[1]
art[k,8] <- gom_h[1]
art[k,11] <- cal_P[1]
art[k,12] <- cal_G[1] 

art[k,4] <- sil_g[2] 
art[k,5] <- sil_h[2] 
art[k,6] <- sil_p[2] 
art[k,9] <- gom_g[2]
art[k,10] <- gom_h[2]
art[k,13] <- cal_G[2]
art[k,14] <- cal_P[2] 

}





names(model) <- spp
names(model_predict) <- spp
names(coefficient_sps) <- spp
names(variable_response) <- spp


data_for_climatic_range <- list(model = model,
                                predictors = predictors,
                                ll <- ll 
                                )

save(data_for_climatic_range, file = "data_for_climatic_range.bin")
load("data_for_climatic_range.bin")
