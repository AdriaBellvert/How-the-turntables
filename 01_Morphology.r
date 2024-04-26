

library(geomorph)
library(factoextra)
library(plyr)

#Load morphometric data
load("Data_morph.bin")


#Data preparation for trajectory analyses
pca.gr <- gm.prcomp(Data_morph$GMM_data)

sil <- "sil_g"

Data_morph$Classifier$spec_type <- Data_morph$Classifier$sp_island
silv <- Data_morph$Classifier[Data_morph$Classifier$sp_island%in%sil,]
silv$spec_type <- revalue(silv$spec_type, c("sil_g" = "sil1"))
new_spec <- rbind(Data_morph$Classifier,silv)
trj <- rbind(Data_morph$Classifier,silv)

trj$col <- trj$spec_type
trj$col <- revalue(trj$spec_type, c("sil_g" = "blue", "sil_h" = "cyan", "sil_p" = "cornflowerblue",
                                    "cal_g" = "red", "cal_p" = "orange", "gom_h" = "darkolivegreen3", "gom_g" = "darkgreen",
                                    "sil1" = "blue"))


trj$spec_type <- revalue(trj$spec_type, c("sil_g" = "sil2", "sil_h" = "sil2", "sil_p" = "sil1",
                                                    "cal_g" = "cal", "cal_p" = "cal", "gom_h" = "gom", "gom_g" = "gom"))

col <- as.character(trj$col)


data_sil <- Data_morph$GMM_data[,,Data_morph$Classifier$sp_island%in%sil]
data_new <- arrayspecs(rbind(two.d.array(Data_morph$GMM_data), 
                             two.d.array(data_sil)),
                       51, 2)
new_spec$is_type <- droplevels(new_spec$is_type)

data_new_m <- two.d.array(data_new)


#Trajectory analyses
rr <- rrpp.data.frame(sh = data_new_m, spps = trj$spec_type,
                      isl = trj$is_type)

fit <- lm.rrpp(sh ~ spps * isl, data = rr)
TA <- trajectory.analysis(fit, groups = trj$spec_type, traj.pts = trj$is_type)


TP <- plot(TA, pch = 22, bg = col,
           cex = 0.7, col = "gray", xlim = c(-0.15,0.1),
           bty ="L", las = 1)


add.trajectories(TP, traj.pch = 21, traj.bg = 1:4)
legend("topleft", levels(new_spec$col), pch =  c(22), 
       pt.bg = c("red" , "orange" , "darkgreen" , "darkolivegreen3" , "blue" ,   "cornflowerblue",  "cyan" ))
legend("topright", levels(new_spec$is_type), pch =  c(21), pt.bg = c(3,2))
legend("topleft", levels(new_spec$is_type), pch =  c(21), pt.bg = ifelse(new_spec$is_type=="High", 3, 4))


summary(TA)
summary(TA, attribute = "TC", angle.type = "deg") 

