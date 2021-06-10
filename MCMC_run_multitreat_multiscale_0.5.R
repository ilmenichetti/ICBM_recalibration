

library(readODS)
library(rjags)
library(R2jags)
library(coda)
library(modeest)
library(viridis)
library(RColorBrewer)
library(imputeTS)
library(modeest)
library(np)

Ultuna_treat_palette<-rev(viridis(16)[2:16])
# ALPHA TO COLOR - cosmetics
# function to add an alpha value to a colour
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2,
        function(x)
          rgb(x[1], x[2], x[3], alpha=alpha))
}

Ultuna_treat_palette_alpha<-add.alpha(Ultuna_treat_palette,0.5)

Ultuna_crop_palette<-brewer.pal(5,"Dark2")
  
# loading Thomas input calculation for plotting them against the calculated inputs
Thomas_2011_inputs<-read_ods("../Data/FRAME56/Thomas_Kätterer_2011_yield_corr_c_och_e_input_eq_depth.ods", sheet=4)



############### ############### ############### ############### ############### 
############### Generating the model object until 1999
############### ############### ############### ############### ############### 

# reading the input data object
#calib_data_Ultuna<- readRDS(file = "Ultuna_input.rds")
calib_data_Ultuna<- readRDS(file = "Ultuna_input_1999.rds") #for excluding the maize years
for(i in 1:length(names(calib_data_Ultuna))){
  name<-names(calib_data_Ultuna)[i]
  Output = paste("../ICBM_recalibration_GLUE/data/",name, ".csv",  sep = "")
  write.csv(calib_data_Ultuna[i], Output)
  
}


yields_list<-calib_data_Ultuna$Yields_cereals_Ultuna
yields_list[calib_data_Ultuna$Yields_cereals_Ultuna>0]<-"cereals"
yields_list[calib_data_Ultuna$Yields_root_crops_Ultuna>0]<-"root_crops"
yields_list[calib_data_Ultuna$Yields_oilseeds_Ultuna>0]<-"oilseeds"
yields_list[calib_data_Ultuna$Yields_maize_Ultuna>0]<-"maize"

write.csv(yields_list[2,], "../ICBM_recalibration_GLUE/data/crop_list.csv", row.names = F)


yields_tot<-calib_data_Ultuna$Yields_cereals_Ultuna
yields_tot[calib_data_Ultuna$Yields_cereals_Ultuna>0]<-calib_data_Ultuna$Yields_cereals_Ultuna[calib_data_Ultuna$Yields_cereals_Ultuna>0]
yields_tot[calib_data_Ultuna$Yields_root_crops_Ultuna>0]<-calib_data_Ultuna$Yields_root_crops_Ultuna[calib_data_Ultuna$Yields_root_crops_Ultuna>0]
yields_tot[calib_data_Ultuna$Yields_oilseeds_Ultuna>0]<-calib_data_Ultuna$Yields_oilseeds_Ultuna[calib_data_Ultuna$Yields_oilseeds_Ultuna>0]
yields_tot[calib_data_Ultuna$Yields_maize_Ultuna>0]<-calib_data_Ultuna$Yields_maize_Ultuna[calib_data_Ultuna$Yields_maize_Ultuna>0]

write.csv(yields_tot, "../ICBM_recalibration_GLUE/data/yields_tot.csv", row.names = F)


#############################################
##### MCMC PART #############################
#############################################

## normalize re_clim on treatment G in Ultuna
#re_normalization<-rowMeans(aggregate(re_Ultuna_long, by=list(Ultuna_treat), FUN=mean)[,-1])[7]

#INIT RATIO PRIORS
# mean_re_Ultuna= mean(as.numeric(re_Ultuna_long[21,]), na.rm=T)
# mean_I_Ultuna= mean(Ultuna_yields_timeseries_long*10^-3, na.rm=T)
# YSS_Ultuna<-mean_I_Ultuna*0.7/(0.8*mean_re_Ultuna)
# OSS_Ultuna<-0.15*(mean_I_Ultuna*0.7/(0.0085*mean_re_Ultuna))
# Init_prior_Ultuna<-1-YSS_Ultuna/OSS_Ultuna


N.ADAPTS=100
N.RUNS=300
N.BURNIN=100
sampling.nr=25
n.chains=50

### prior for k1
# Xie, Yajun. 2020. “A Meta-Analysis of Critique of Litterbag Method Used in Examining Decomposition of Leaf Litters.” Journal of Soils and Sediments 20 (4): 1881–86. https://doi.org/10.1007/s11368-020-02572-9.
#digitized from plot 2 fig 3
Xie_data_k1<-c(
  0.4389654645547525,
  0.5036980104731434,
  0.5418680346882365,
  0.5475953903352611,
  0.7484834314129486,
  0.697937318560429,
  0.6740416981123342,
  0.5938448289432867,
  0.717837217162621,
  0.895107440006297,
  1.2372810571300257)
mean(Xie_data_k1)
sd(Xie_data_k1)
range(Xie_data_k1)
quantile(Xie_data_k1, c(0.05, 0.95))


#N=newyears
###SELECT SPECIFIC TREATMENTS
#selection_treatments<-c(1:8,10,11)
#selection_treatments<-c(1:15)[-c(9,13,15)]
#selection_treatments<-c(1)
selection_treatments<-c(1:15)
LETTERS[selection_treatments]

dim(calib_data_Ultuna$SOC_Ultuna)
dim(calib_data_Ultuna$re_Ultuna)

for(i in 1:(length(calib_data_Ultuna)-1)){
  if(length(dim(calib_data_Ultuna[[i]]))>0){
  calib_data_Ultuna[[i]]<-calib_data_Ultuna[[i]][selection_treatments,]
  }else {
    calib_data_Ultuna[[i]]<-calib_data_Ultuna[[i]][selection_treatments]
  }
}

calib_data_Ultuna$J_Ultuna<-length(selection_treatments)
calib_data_Ultuna$N_Ultuna<-as.numeric(calib_data_Ultuna$N_Ultuna[1])

calib_data_Ultuna$SOC_Ultuna<-as.data.frame(calib_data_Ultuna$SOC_Ultuna)
#calib_data_Ultuna$error_SOC_Ultuna<-as.data.frame(calib_data_Ultuna$error_SOC_Ultuna)


######### FILTERING THE DATA #########
FILTER=T
threshold_rate=0.75 #ratio of the total variance 


outliers<-list()
png("Outliers.png", height=4500, width=4000, res=320)
par(mfrow=c(5,3))
time_evc_single<-seq(from=1956, length.out=length(as.numeric(calib_data_Ultuna$SOC_Ultuna[i,])))
for(i in 1:dim(calib_data_Ultuna$SOC_Ultuna)[1]){
  #rgr_model<-lm(as.numeric(calib_data_Ultuna$SOC_Ultuna[i,])~poly(time_evc_single, 3))
  ann_data<-data.frame(as.numeric(calib_data_Ultuna$SOC_Ultuna[i,]),time_evc_single)
  names(ann_data)<-c("SOC","Time")
  rgr_model<-loess(SOC~Time, data=ann_data , span=0.85)
  pred<-predict(rgr_model, newdata=(time_evc_single))
  plot(as.numeric(calib_data_Ultuna$SOC_Ultuna[i,])~time_evc_single, main=paste("OUtlier detection, treatment ", LETTERS[i]), xlab="time", ylab="SOC t ha-1")
  lines(time_evc_single, pred, lty=2)
  threshold=max(abs(as.numeric(calib_data_Ultuna$SOC_Ultuna[i,])-pred), na.rm=T)*threshold_rate
  outliers[[i]]<-which(abs(as.numeric(calib_data_Ultuna$SOC_Ultuna[i,])-pred)>threshold)
  points(as.numeric(calib_data_Ultuna$SOC_Ultuna[i,])[outliers[[i]]]~time_evc_single[outliers[[i]]], col="red", pch=16)
  }
dev.off()
  
if(FILTER==T){
  for(i in 1:dim(calib_data_Ultuna$SOC_Ultuna)[1]){
      calib_data_Ultuna$SOC_Ultuna[1,][outliers[[1]]]<-NA
  }
}
#########




parameter_list<-c("k1_ult",
                  "k2_ult",
                  "h_R_ult",
                  "h_S_ult",
                  "h_FYM_ult",
                  "h_PEA_ult",
                  "h_SAW_ult",
                  "h_SLU_ult",
                  "stubbles_ratio_Ultuna",
                  "stubbles_ratio_Ultuna_maize",
                  "SR_cereals_ult","SR_root_crops_ult","SR_oilseeds_ult", "SR_maize_ult",
                  "Tot_Ultuna",
                  "Y_tot_Ultuna",
                  "O_Ultuna",
                  "I_R_Ultuna",
                  "I_S_Ultuna",
                  "fluxR_Ultuna",
                  "fluxS_Ultuna",
                  "Y_R_Ultuna",
                  "Init_ratio_Ultuna",
                  "alpha",
                  "alpha_maize",
                  "C_percent",
                  "R_Ultuna",
                  "exudates_coeff",
                  "flux_sum_Ultuna",
                  "Y_FYM_Ultuna",
                  "Y_GM_Ultuna",
                  "Y_PEA_Ultuna",
                  "Y_SAW_Ultuna",
                  "Y_SLU_Ultuna",
                  "Y_STR_Ultuna",
                  "Y_S_Ultuna",
                  "I_R_cereals_Ultuna",
                  "I_R_root_crops_Ultuna",
                  "I_R_oilseeds_Ultuna",
                  "I_R_maize_Ultuna")

  model.ICBM.matrix_Ultuna<- jags.model('./JAGS_ICBM_4.R',
                                        data=calib_data_Ultuna,
                                        n.chains = n.chains,
                                        n.adapt = N.ADAPTS)
  
  update(model.ICBM.matrix_Ultuna, n.iter=N.RUNS, n.thin=10, n.burnin=N.BURNIN)
  mcmc.array.ICBM.Ultuna<-(jags.samples(model.ICBM.matrix_Ultuna,parameter_list, sampling.nr))

  
###Ultuna
#Plot SOC and Input (roots and shoots) simulation
mcmc.list.SOC.Ultuna<-as.mcmc.list(mcmc.array.ICBM.Ultuna$Tot_Ultuna, chains=F)
mcmc.unlist.SOC.Ultuna<-mcmc.list.SOC.Ultuna[[1]]
for(i in 1:n.chains){
  mcmc.unlist.SOC.Ultuna<-rbind(mcmc.unlist.SOC.Ultuna, mcmc.list.SOC.Ultuna[[i]])
}

mcmc.list.Ir.Ultuna<-as.mcmc.list(mcmc.array.ICBM.Ultuna$I_R_Ultuna, chains=F)
mcmc.unlist.Ir.Ultuna<-(mcmc.list.Ir.Ultuna[[1]])
for(i in 1:n.chains){
  mcmc.unlist.Ir.Ultuna<-rbind(mcmc.unlist.Ir.Ultuna, mcmc.list.Ir.Ultuna[[i]])
}

mcmc.list.Is.Ultuna<-as.mcmc.list(mcmc.array.ICBM.Ultuna$I_S_Ultuna, chains=F)
mcmc.unlist.Is.Ultuna<-(mcmc.list.Is.Ultuna[[1]])
for(i in 1:n.chains){
  mcmc.unlist.Is.Ultuna<-rbind(mcmc.unlist.Is.Ultuna, mcmc.list.Is.Ultuna[[i]])
}

mcmc.list.fluxR.Ultuna<-as.mcmc.list(mcmc.array.ICBM.Ultuna$I_S_Ultuna, chains=F)
mcmc.unlist.fluxR.Ultuna<-(mcmc.list.fluxR.Ultuna[[1]])
for(i in 1:n.chains){
  mcmc.unlist.fluxR.Ultuna<-rbind(mcmc.unlist.fluxR.Ultuna, mcmc.list.fluxR.Ultuna[[i]])
}

mcmc.list.fluxS.Ultuna<-as.mcmc.list(mcmc.array.ICBM.Ultuna$I_S_Ultuna, chains=F)
mcmc.unlist.fluxS.Ultuna<-(mcmc.list.fluxS.Ultuna[[1]])
for(i in 1:n.chains){
  mcmc.unlist.fluxS.Ultuna<-rbind(mcmc.unlist.fluxS.Ultuna, mcmc.list.fluxS.Ultuna[[i]])
}

mcmc.list.fluxSUM.Ultuna<-as.mcmc.list(mcmc.array.ICBM.Ultuna$I_S_Ultuna, chains=F)
mcmc.unlist.fluxSUM.Ultuna<-(mcmc.list.fluxSUM.Ultuna[[1]])
for(i in 1:n.chains){
  mcmc.unlist.fluxSUM.Ultuna<-rbind(mcmc.unlist.fluxSUM.Ultuna, mcmc.list.fluxSUM.Ultuna[[i]])
}



#extract the array of the realizations*treatments*years
#dim<-dim(aggregate(Ultuna_SOC_timeseries_long, by=list(Ultuna_treat), FUN=mean)[,-1])
dim<-dim(calib_data_Ultuna$SOC_Ultuna)
dim[1]<-calib_data_Ultuna$J_Ultuna
Ultuna_prediction_array<-array( ,dim=c((sampling.nr*n.chains),dim))
str(Ultuna_prediction_array)
length_treat<-dim[1]
for(j in 1:(sampling.nr*n.chains)){
  for(i in 1:(dim[2])){
    length_treat_add<-length_treat*i
    Ultuna_prediction_array[j,,i]<-mcmc.unlist.SOC.Ultuna[j,(1-length_treat+length_treat_add):length_treat_add]
  }
}
#plot(Ultuna_prediction_array[20,2,], type="l")


Ultuna_inputS_array<-array( ,dim=c((sampling.nr*n.chains),dim))
for(j in 1:(sampling.nr*n.chains)){
  length_treat<-dim[1]
  for(i in 1:(dim[2])){
    length_treat_add<-length_treat*i
    Ultuna_inputS_array[j,,i]<-mcmc.unlist.Is.Ultuna[j,(1-length_treat+length_treat_add):length_treat_add]
  }
}

Ultuna_inputR_array<-array( ,dim=c((sampling.nr*n.chains),dim))
for(j in 1:(sampling.nr*n.chains)){
  length_treat<-dim[1]
  for(i in 1:(dim[2])){
    length_treat_add<-length_treat*i
    Ultuna_inputR_array[j,,i]<-mcmc.unlist.Ir.Ultuna[j,(1-length_treat+length_treat_add):length_treat_add]
  }
}

Ultuna_fluxS_array<-array( ,dim=c((sampling.nr*n.chains),dim))
for(j in 1:(sampling.nr*n.chains)){
  length_treat<-dim[1]
  for(i in 1:(dim[2])){
    length_treat_add<-length_treat*i
    Ultuna_fluxS_array[j,,i]<-mcmc.unlist.fluxS.Ultuna[j,(1-length_treat+length_treat_add):length_treat_add]
  }
}

Ultuna_fluxR_array<-array( ,dim=c((sampling.nr*n.chains),dim))
for(j in 1:(sampling.nr*n.chains)){
  length_treat<-dim[1]
  for(i in 1:(dim[2])){
    length_treat_add<-length_treat*i
    Ultuna_fluxR_array[j,,i]<-mcmc.unlist.fluxR.Ultuna[j,(1-length_treat+length_treat_add):length_treat_add]
  }
}

Ultuna_fluxSUM_array<-array( ,dim=c((sampling.nr*n.chains),dim))
for(j in 1:(sampling.nr*n.chains)){
  length_treat<-dim[1]
  for(i in 1:(dim[2])){
    length_treat_add<-length_treat*i
    Ultuna_fluxSUM_array[j,,i]<-mcmc.unlist.fluxSUM.Ultuna[j,(1-length_treat+length_treat_add):length_treat_add]
  }
}

dim(mcmc.unlist.SOC.Ultuna)

str(Ultuna_prediction_array)

plot(Ultuna_prediction_array[2,2,], type="l")
lines(Ultuna_prediction_array[3,2,], type="l")



#mean, min and max matrices by plot
Ultuna_mean_predictions<-(colMeans(Ultuna_prediction_array))
Ultuna_max_predictions<-apply(Ultuna_prediction_array, MARGIN=c(2,3), max)
Ultuna_min_predictions<-apply(Ultuna_prediction_array, MARGIN=c(2,3), min)
rownames(Ultuna_mean_predictions)<-as.factor(rownames(calib_data_Ultuna$SOC_Ultuna))[1:calib_data_Ultuna$J_Ultuna]
rownames(Ultuna_max_predictions)<-as.factor(rownames(calib_data_Ultuna$SOC_Ultuna))[1:calib_data_Ultuna$J_Ultuna]
rownames(Ultuna_min_predictions)<-as.factor(rownames(calib_data_Ultuna$SOC_Ultuna))[1:calib_data_Ultuna$J_Ultuna]

Ultuna_mean_Ir<-(colMeans(Ultuna_inputR_array))
Ultuna_max_Ir<-apply(Ultuna_inputR_array, MARGIN=c(2,3), max)
Ultuna_min_Ir<-apply(Ultuna_inputR_array, MARGIN=c(2,3), min)
rownames(Ultuna_mean_Ir)<-as.factor(rownames(calib_data_Ultuna$SOC_Ultuna))[1:calib_data_Ultuna$J_Ultuna]
rownames(Ultuna_max_Ir)<-as.factor(rownames(calib_data_Ultuna$SOC_Ultuna))[1:calib_data_Ultuna$J_Ultuna]
rownames(Ultuna_min_Ir)<-as.factor(rownames(calib_data_Ultuna$SOC_Ultuna))[1:calib_data_Ultuna$J_Ultuna]

Ultuna_mean_Is<-(colMeans(Ultuna_inputS_array))
Ultuna_max_Is<-apply(Ultuna_inputS_array, MARGIN=c(2,3), max)
Ultuna_min_Is<-apply(Ultuna_inputS_array, MARGIN=c(2,3), min)
rownames(Ultuna_mean_Is)<-as.factor(rownames(calib_data_Ultuna$SOC_Ultuna))[1:calib_data_Ultuna$J_Ultuna]
rownames(Ultuna_max_Is)<-as.factor(rownames(calib_data_Ultuna$SOC_Ultuna))[1:calib_data_Ultuna$J_Ultuna]
rownames(Ultuna_min_Is)<-as.factor(rownames(calib_data_Ultuna$SOC_Ultuna))[1:calib_data_Ultuna$J_Ultuna]

Ultuna_mean_fluxS<-(colMeans(Ultuna_fluxS_array))
Ultuna_max_fluxS<-apply(Ultuna_fluxS_array, MARGIN=c(2,3), max)
Ultuna_min_fluxS<-apply(Ultuna_fluxS_array, MARGIN=c(2,3), min)
rownames(Ultuna_mean_fluxS)<-as.factor(rownames(calib_data_Ultuna$SOC_Ultuna))[1:calib_data_Ultuna$J_Ultuna]
rownames(Ultuna_max_fluxS)<-as.factor(rownames(calib_data_Ultuna$SOC_Ultuna))[1:calib_data_Ultuna$J_Ultuna]
rownames(Ultuna_min_fluxS)<-as.factor(rownames(calib_data_Ultuna$SOC_Ultuna))[1:calib_data_Ultuna$J_Ultuna]

Ultuna_mean_fluxR<-(colMeans(Ultuna_fluxR_array))
Ultuna_max_fluxR<-apply(Ultuna_fluxR_array, MARGIN=c(2,3), max)
Ultuna_min_fluxR<-apply(Ultuna_fluxR_array, MARGIN=c(2,3), min)
rownames(Ultuna_mean_fluxR)<-as.factor(rownames(calib_data_Ultuna$SOC_Ultuna))[1:calib_data_Ultuna$J_Ultuna]
rownames(Ultuna_max_fluxR)<-as.factor(rownames(calib_data_Ultuna$SOC_Ultuna))[1:calib_data_Ultuna$J_Ultuna]
rownames(Ultuna_min_fluxR)<-as.factor(rownames(calib_data_Ultuna$SOC_Ultuna))[1:calib_data_Ultuna$J_Ultuna]

Ultuna_mean_fluxSUM<-(colMeans(Ultuna_fluxSUM_array))
Ultuna_max_fluxSUM<-apply(Ultuna_fluxSUM_array, MARGIN=c(2,3), max)
Ultuna_min_fluxSUM<-apply(Ultuna_fluxSUM_array, MARGIN=c(2,3), min)
rownames(Ultuna_mean_fluxSUM)<-as.factor(rownames(calib_data_Ultuna$SOC_Ultuna))[1:calib_data_Ultuna$J_Ultuna]
rownames(Ultuna_max_fluxSUM)<-as.factor(rownames(calib_data_Ultuna$SOC_Ultuna))[1:calib_data_Ultuna$J_Ultuna]
rownames(Ultuna_min_fluxSUM)<-as.factor(rownames(calib_data_Ultuna$SOC_Ultuna))[1:calib_data_Ultuna$J_Ultuna]



#mean, min and max matrices by treatment
#Ultuna_treat_vec<-LETTERS[yields_Ultuna[yields_Ultuna$year==1996,]$treat]
#Ultuna_treat_vec_old<-LETTERS[yields_Ultuna[yields_Ultuna$year==1956,]$treat]
Ultuna_treat_vec<-as.factor(rownames(calib_data_Ultuna$SOC_Ultuna))
Ultuna_mean_predictions_bytreat<-mat.or.vec(calib_data_Ultuna$J_Ultuna,dim[2])
Ultuna_max_predictions_bytreat<-mat.or.vec(calib_data_Ultuna$J_Ultuna,dim[2])
Ultuna_min_predictions_bytreat<-mat.or.vec(calib_data_Ultuna$J_Ultuna,dim[2])
Ultuna_mean_Ir_bytreat<-mat.or.vec(calib_data_Ultuna$J_Ultuna,dim[2])
Ultuna_max_Ir_bytreat<-mat.or.vec(calib_data_Ultuna$J_Ultuna,dim[2])
Ultuna_min_Ir_bytreat<-mat.or.vec(calib_data_Ultuna$J_Ultuna,dim[2])
Ultuna_mean_Is_bytreat<-mat.or.vec(calib_data_Ultuna$J_Ultuna,dim[2])
Ultuna_max_Is_bytreat<-mat.or.vec(calib_data_Ultuna$J_Ultuna,dim[2])
Ultuna_min_Is_bytreat<-mat.or.vec(calib_data_Ultuna$J_Ultuna,dim[2])
Ultuna_mean_fluxS_bytreat<-mat.or.vec(calib_data_Ultuna$J_Ultuna,dim[2])
Ultuna_max_fluxS_bytreat<-mat.or.vec(calib_data_Ultuna$J_Ultuna,dim[2])
Ultuna_min_fluxS_bytreat<-mat.or.vec(calib_data_Ultuna$J_Ultuna,dim[2])
Ultuna_mean_fluxR_bytreat<-mat.or.vec(calib_data_Ultuna$J_Ultuna,dim[2])
Ultuna_max_fluxR_bytreat<-mat.or.vec(calib_data_Ultuna$J_Ultuna,dim[2])
Ultuna_min_fluxR_bytreat<-mat.or.vec(calib_data_Ultuna$J_Ultuna,dim[2])
Ultuna_mean_fluxSUM_bytreat<-mat.or.vec(calib_data_Ultuna$J_Ultuna,dim[2])
Ultuna_max_fluxSUM_bytreat<-mat.or.vec(calib_data_Ultuna$J_Ultuna,dim[2])
Ultuna_min_fluxSUM_bytreat<-mat.or.vec(calib_data_Ultuna$J_Ultuna,dim[2])
Ultuna_measured_bytreat<-mat.or.vec(calib_data_Ultuna$J_Ultuna,dim[2])
for(i in 1:calib_data_Ultuna$J_Ultuna){
  Ultuna_mean_predictions_bytreat[i,]<-Ultuna_mean_predictions[i,]
  Ultuna_min_predictions_bytreat[i,]<-Ultuna_min_predictions[i,]
  Ultuna_max_predictions_bytreat[i,]<-Ultuna_max_predictions[i,]
  Ultuna_mean_Ir_bytreat[i,]<-Ultuna_mean_Ir[i,]
  Ultuna_min_Ir_bytreat[i,]<-Ultuna_min_Ir[i,]
  Ultuna_max_Ir_bytreat[i,]<-Ultuna_max_Ir[i,]
  Ultuna_mean_Is_bytreat[i,]<-Ultuna_mean_Is[i,]
  Ultuna_min_Is_bytreat[i,]<-Ultuna_min_Is[i,]
  Ultuna_max_Is_bytreat[i,]<-Ultuna_max_Is[i,]
  Ultuna_mean_fluxS_bytreat[i,]<-Ultuna_mean_fluxS[i,]
  Ultuna_min_fluxS_bytreat[i,]<-Ultuna_min_fluxS[i,]
  Ultuna_max_fluxS_bytreat[i,]<-Ultuna_max_fluxS[i,]
  Ultuna_mean_fluxR_bytreat[i,]<-Ultuna_mean_fluxR[i,]
  Ultuna_min_fluxR_bytreat[i,]<-Ultuna_min_fluxR[i,]
  Ultuna_max_fluxR_bytreat[i,]<-Ultuna_max_fluxR[i,]
  Ultuna_mean_fluxSUM_bytreat[i,]<-Ultuna_mean_fluxSUM[i,]
  Ultuna_min_fluxSUM_bytreat[i,]<-Ultuna_min_fluxSUM[i,]
  Ultuna_max_fluxSUM_bytreat[i,]<-Ultuna_max_fluxSUM[i,]
  Ultuna_measured_bytreat[i,]<-as.numeric(calib_data_Ultuna$SOC_Ultuna[i,])
}


write.csv(t(Ultuna_mean_Ir_bytreat), file="Ultuna_mean_Ir_bytreat.csv")
write.csv(t(Ultuna_mean_Is_bytreat), file="Ultuna_mean_Is_bytreat.csv")

Ultuna_yields_table<-calib_data_Ultuna$Yields_cereals_Ultuna+calib_data_Ultuna$Yields_root_crops_Ultuna+calib_data_Ultuna$Yields_oilseeds_Ultuna+calib_data_Ultuna$Yields_maize_Ultuna
Ultuna_Ir_table<-as.data.frame(Ultuna_mean_Ir_bytreat)
Ultuna_Is_table<-as.data.frame(Ultuna_mean_Is_bytreat)
Ultuna_Pred_table<-as.data.frame(Ultuna_max_predictions_bytreat)

Ultuna_Is_table<-rbind(Ultuna_Is_table, seq(from=1956, to=(1956+71)))
Ultuna_Ir_table<-rbind(Ultuna_Ir_table, seq(from=1956, to=(1956+71)))
Ultuna_yields_table<-rbind(Ultuna_yields_table, seq(from=1956, to=(1956+71)))
Ultuna_Pred_table<-rbind(Ultuna_Pred_table, seq(from=1956, to=(1956+71)))

rownames(Ultuna_Ir_table)[16]<-"Year"
rownames(Ultuna_yields_table)[16]<-"Year"
rownames(Ultuna_Is_table)[16]<-"Year"
rownames(Ultuna_Pred_table)[16]<-"Year"


write_ods(Ultuna_yields_table, "Ultuna_input_check.ods", sheet="Yields", row_names = T)
write_ods(Ultuna_Ir_table, "Ultuna_input_check.ods", sheet="Calculated root inputs (C,ton ha-1)", append=T, row_names = T)
write_ods(Ultuna_Is_table, "Ultuna_input_check.ods", sheet="Calculated shoot inputs (C,ton ha-1)", append=T, row_names = T)
write_ods(Ultuna_Pred_table, "Ultuna_input_check.ods", sheet="Predicted SOC stocks (C,ton ha-1)", append=T, row_names = T)

Ultuna_Thomas_Ir_table<-data.frame(mat.or.vec(15,54))
for(i in 1:calib_data_Ultuna$J_Ultuna){
  Ultuna_Thomas_Ir_table[i,]<-as.numeric(Thomas_2011_inputs[Thomas_2011_inputs$treat==letters[i],]$`CR+CE70%`[1:54])*0.001
}

colnames(Ultuna_Thomas_Ir_table)<-colnames(calib_data_Ultuna$SOC_Ultuna)[1:54]
rownames(Ultuna_Thomas_Ir_table)<-LETTERS[1:calib_data_Ultuna$J_Ultuna]

write_ods(Ultuna_Thomas_Ir_table, "Ultuna_input_check.ods", sheet="Thomas calculations IR (C,ton ha-1)", append=T, row_names = T)



png("ICBM_predictions_Ultuna_specific.png", height=4500, width=4000, res=320)
par(mfrow=c(5,3))
for(i in 1:calib_data_Ultuna$J_Ultuna){
  yearseq<-seq(from=colnames(calib_data_Ultuna$SOC_Ultuna)[1], to=as.numeric(colnames(calib_data_Ultuna$SOC_Ultuna)[1])+dim[2]-1)
  lastyear<-dim(Ultuna_mean_predictions_bytreat)[2]-1
  plot(yearseq, Ultuna_mean_predictions_bytreat[i,], ylim=c(20,125), type="l",  ylab=expression(paste("C (t ha" ^-1, ")")), xlab="Years", main=paste("Treatment", LETTERS[selection_treatments[i]]), col=Ultuna_treat_palette[i])
  polygon(c(yearseq,rev(yearseq)),
          c(Ultuna_max_predictions_bytreat[i,],rev(Ultuna_min_predictions_bytreat[i,])),
          col=Ultuna_treat_palette_alpha[i],border=Ultuna_treat_palette[i])
  lines(yearseq, Ultuna_mean_predictions_bytreat[i,], col="black", lty=1, lwd=1)
  #lines(yearseq[2:62], na_interpolation(Ultuna_measured_bytreat[i,1:61]), col="firebrick", lty=2, lwd=2)
  #lines(yearseq, na_interpolation(Ultuna_measured_bytreat[i,1:53]), col=add.alpha("firebrick",0.2), lty=2, lwd=2)
  points(yearseq, (Ultuna_measured_bytreat[i,]), col="firebrick", lty=2, lwd=2)
}
dev.off()



png("ICBM_predictions_Ultuna_specific_freescale.png", height=4500, width=4000, res=320)
par(mfrow=c(5,3))
for(i in 1:calib_data_Ultuna$J_Ultuna){
  yearseq<-seq(from=colnames(calib_data_Ultuna$SOC_Ultuna)[1], to=as.numeric(colnames(calib_data_Ultuna$SOC_Ultuna)[1])+dim[2]-1)
  lastyear<-dim(Ultuna_mean_predictions_bytreat)[2]-1
  plot(yearseq, Ultuna_mean_predictions_bytreat[i,], ylim=c(min(Ultuna_mean_predictions_bytreat[i,])*0.3,max(Ultuna_mean_predictions_bytreat[i,])*1.6), type="l",  ylab=expression(paste("C (t ha" ^-1, ")")), xlab="Years", main=paste("Treatment", LETTERS[selection_treatments[i]]), col=Ultuna_treat_palette[i])
  polygon(c(yearseq,rev(yearseq)),
          c(Ultuna_max_predictions_bytreat[i,],rev(Ultuna_min_predictions_bytreat[i,])),
          col=Ultuna_treat_palette_alpha[i],border=Ultuna_treat_palette[i])
  lines(yearseq, Ultuna_mean_predictions_bytreat[i,], col="black", lty=1, lwd=1)
  #lines(yearseq[2:62], na_interpolation(Ultuna_measured_bytreat[i,1:61]), col="firebrick", lty=2, lwd=2)
  #lines(yearseq, na_interpolation(Ultuna_measured_bytreat[i,1:53]), col=add.alpha("firebrick",0.2), lty=2, lwd=2)
  points(yearseq, (Ultuna_measured_bytreat[i,]), col="firebrick", lty=2, lwd=2)
}
dev.off()



png("ICBM_Ir_Ultuna_specific.png", height=4500, width=4000, res=320)
par(mfrow=c(5,3))
for(i in 1:calib_data_Ultuna$J_Ultuna){
  yearseq<-seq(from=colnames(calib_data_Ultuna$SOC_Ultuna)[1], to=as.numeric(colnames(calib_data_Ultuna$SOC_Ultuna)[1])+dim[2]-1)
  lastyear<-dim(Ultuna_mean_Ir_bytreat)[2]-1
  plot(yearseq, Ultuna_mean_Ir_bytreat[i,], ylim=c(0,2), type="l",  ylab=expression(paste("C (t ha" ^-1, ")")), xlab="Years", main=paste("Treatment", LETTERS[selection_treatments[i]]), col=Ultuna_treat_palette[i])
  polygon(c(yearseq,rev(yearseq)),
          c(Ultuna_max_Ir_bytreat[i,],rev(Ultuna_min_Ir_bytreat[i,])),
          col=Ultuna_treat_palette_alpha[i],border=Ultuna_treat_palette[i])
  lines(yearseq, Ultuna_mean_Ir_bytreat[i,], col="black", lty=1, lwd=1)
}
dev.off()

png("ICBM_Is_Ultuna_specific.png", height=4500, width=4000, res=320)
par(mfrow=c(5,3))
for(i in 1:calib_data_Ultuna$J_Ultuna){
  yearseq<-seq(from=colnames(calib_data_Ultuna$SOC_Ultuna)[1], to=as.numeric(colnames(calib_data_Ultuna$SOC_Ultuna)[1])+dim[2]-1)
  lastyear<-dim(Ultuna_mean_Is_bytreat)[2]-1
  plot(yearseq, Ultuna_mean_Is_bytreat[i,], ylim=c(0,2), type="l",  ylab=expression(paste("C (t ha" ^-1, ")")), xlab="Years", main=paste("Treatment", LETTERS[selection_treatments[i]]), col=Ultuna_treat_palette[i])
  polygon(c(yearseq,rev(yearseq)),
          c(Ultuna_max_Is_bytreat[i,],rev(Ultuna_min_Is_bytreat[i,])),
          col=Ultuna_treat_palette_alpha[i],border=Ultuna_treat_palette[i])
  lines(yearseq, Ultuna_mean_Is_bytreat[i,], col="black", lty=1, lwd=1)
}
dev.off()

png("ICBM_fluxS_Ultuna_specific.png", height=4500, width=4000, res=320)
par(mfrow=c(5,3))
for(i in 1:calib_data_Ultuna$J_Ultuna){
  yearseq<-seq(from=colnames(calib_data_Ultuna$SOC_Ultuna)[1], to=as.numeric(colnames(calib_data_Ultuna$SOC_Ultuna)[1])+dim[2]-1)
  lastyear<-dim(Ultuna_mean_fluxS_bytreat)[2]-1
  plot(yearseq, Ultuna_mean_fluxS_bytreat[i,], ylim=c(0,2), type="l",  ylab=expression(paste("C (t ha" ^-1, ")")), xlab="Years", main=paste("Treatment", LETTERS[selection_treatments[i]]), col=Ultuna_treat_palette[i])
  polygon(c(yearseq,rev(yearseq)),
          c(Ultuna_max_fluxS_bytreat[i,],rev(Ultuna_min_fluxS_bytreat[i,])),
          col=Ultuna_treat_palette_alpha[i],border=Ultuna_treat_palette[i])
  lines(yearseq, Ultuna_mean_fluxS_bytreat[i,], col="black", lty=1, lwd=1)
}
dev.off()

png("ICBM_fluxR_Ultuna_specific.png", height=4500, width=4000, res=320)
par(mfrow=c(5,3))
for(i in 1:calib_data_Ultuna$J_Ultuna){
  yearseq<-seq(from=colnames(calib_data_Ultuna$SOC_Ultuna)[1], to=as.numeric(colnames(calib_data_Ultuna$SOC_Ultuna)[1])+dim[2]-1)
  lastyear<-dim(Ultuna_mean_fluxR_bytreat)[2]-1
  plot(yearseq, Ultuna_mean_fluxR_bytreat[i,], ylim=c(0.4,-0.5), type="l",  ylab=expression(paste("C (t ha" ^-1, ")")), xlab="Years", main=paste("Treatment", LETTERS[selection_treatments[i]]), col=Ultuna_treat_palette[i])
  polygon(c(yearseq,rev(yearseq)),
          c(Ultuna_max_fluxR_bytreat[i,],rev(Ultuna_min_fluxR_bytreat[i,])),
          col=Ultuna_treat_palette_alpha[i],border=Ultuna_treat_palette[i])
  lines(yearseq, Ultuna_mean_fluxR_bytreat[i,], col="black", lty=1, lwd=1)
}
dev.off()

png("ICBM_fluxSUM_Ultuna_specific.png", height=4500, width=4000, res=320)
par(mfrow=c(5,3))
for(i in 1:calib_data_Ultuna$J_Ultuna){
  yearseq<-seq(from=colnames(calib_data_Ultuna$SOC_Ultuna)[1], to=as.numeric(colnames(calib_data_Ultuna$SOC_Ultuna)[1])+dim[2]-1)
  lastyear<-dim(Ultuna_mean_fluxSUM_bytreat)[2]-1
  plot(yearseq, Ultuna_mean_fluxSUM_bytreat[i,], ylim=c(0.4,-0.5), type="l",  ylab=expression(paste("C (t ha" ^-1, ")")), xlab="Years", main=paste("Treatment", LETTERS[selection_treatments[i]]), col=Ultuna_treat_palette[i])
  polygon(c(yearseq,rev(yearseq)),
          c(Ultuna_max_fluxSUM_bytreat[i,],rev(Ultuna_min_fluxSUM_bytreat[i,])),
          col=Ultuna_treat_palette_alpha[i],border=Ultuna_treat_palette[i])
  lines(yearseq, Ultuna_mean_fluxSUM_bytreat[i,], col="black", lty=1, lwd=1)
}
dev.off()









############### ############### ############### ############### ############### 
############### Generating the model object until 2019
############### ############### ############### ############### ############### 

# reading the input data object
calib_data_Ultuna_long<- readRDS(file = "Ultuna_input.rds")
for(i in 1:length(names(calib_data_Ultuna_long))){
  name<-names(calib_data_Ultuna_long)[i]
  Output = paste("../ICBM_recalibration_GLUE/data/",name, "_long.csv",  sep = "")
  write.csv(calib_data_Ultuna_long[i], Output)
  
}


yields_list_long<-calib_data_Ultuna_long$Yields_cereals_Ultuna
yields_list_long[calib_data_Ultuna_long$Yields_cereals_Ultuna>0]<-"cereals"
yields_list_long[calib_data_Ultuna_long$Yields_root_crops_Ultuna>0]<-"root_crops"
yields_list_long[calib_data_Ultuna_long$Yields_oilseeds_Ultuna>0]<-"oilseeds"
yields_list_long[calib_data_Ultuna_long$Yields_maize_Ultuna>0]<-"maize"

write.csv(yields_list_long[2,], "../ICBM_recalibration_GLUE/data/crop_list_long.csv", row.names = F)


yields_tot_long<-calib_data_Ultuna_long$Yields_cereals_Ultuna
yields_tot_long[calib_data_Ultuna_long$Yields_cereals_Ultuna>0]<-calib_data_Ultuna_long$Yields_cereals_Ultuna[calib_data_Ultuna_long$Yields_cereals_Ultuna>0]
yields_tot_long[calib_data_Ultuna_long$Yields_root_crops_Ultuna>0]<-calib_data_Ultuna_long$Yields_root_crops_Ultuna[calib_data_Ultuna_long$Yields_root_crops_Ultuna>0]
yields_tot_long[calib_data_Ultuna_long$Yields_oilseeds_Ultuna>0]<-calib_data_Ultuna_long$Yields_oilseeds_Ultuna[calib_data_Ultuna_long$Yields_oilseeds_Ultuna>0]
yields_tot_long[calib_data_Ultuna_long$Yields_maize_Ultuna>0]<-calib_data_Ultuna_long$Yields_maize_Ultuna[calib_data_Ultuna_long$Yields_maize_Ultuna>0]

write.csv(yields_tot_long, "../ICBM_recalibration_GLUE/data/yields_tot_long.csv", row.names = F)


##### MCMC PART #############################

for(i in 1:(length(calib_data_Ultuna_long)-1)){
  if(length(dim(calib_data_Ultuna_long[[i]]))>0){
    calib_data_Ultuna_long[[i]]<-calib_data_Ultuna_long[[i]][selection_treatments,]
  }else {
    calib_data_Ultuna_long[[i]]<-calib_data_Ultuna_long[[i]][selection_treatments]
  }
}

calib_data_Ultuna_long$J_Ultuna<-length(selection_treatments)
calib_data_Ultuna_long$N_Ultuna<-as.numeric(calib_data_Ultuna_long$N_Ultuna[1])

calib_data_Ultuna_long$SOC_Ultuna<-as.data.frame(calib_data_Ultuna_long$SOC_Ultuna)


######### FILTERING THE DATA ######### (PARAMETERS FOR FILTERING DEFINED ABOVE IN THE PRE 1999 CALIBRATION)
outliers<-list()
png("Outliers_long.png", height=4500, width=4000, res=320)
par(mfrow=c(5,3))
time_evc_single<-seq(from=1956, length.out=length(as.numeric(calib_data_Ultuna_long$SOC_Ultuna[i,])))
for(i in 1:dim(calib_data_Ultuna_long$SOC_Ultuna)[1]){
  ann_data<-data.frame(as.numeric(calib_data_Ultuna_long$SOC_Ultuna[i,]),time_evc_single)
  names(ann_data)<-c("SOC","Time")
  rgr_model<-loess(SOC~Time, data=ann_data , span=0.85)
  pred<-predict(rgr_model, newdata=(time_evc_single))
  plot(as.numeric(calib_data_Ultuna_long$SOC_Ultuna[i,])~time_evc_single, main=paste("OUtlier detection, treatment ", LETTERS[i]), xlab="time", ylab="SOC t ha-1")
  lines(time_evc_single, pred, lty=2)
  threshold=max(abs(as.numeric(calib_data_Ultuna_long$SOC_Ultuna[i,])-pred), na.rm=T)*threshold_rate
  outliers[[i]]<-which(abs(as.numeric(calib_data_Ultuna_long$SOC_Ultuna[i,])-pred)>threshold)
  points(as.numeric(calib_data_Ultuna_long$SOC_Ultuna[i,])[outliers[[i]]]~time_evc_single[outliers[[i]]], col="red", pch=16)
}
dev.off()

if(FILTER==T){
  for(i in 1:dim(calib_data_Ultuna_long$SOC_Ultuna)[1]){
    calib_data_Ultuna_long$SOC_Ultuna[1,][outliers[[1]]]<-NA
  }
}
#########



model.ICBM.matrix_Ultuna_long<- jags.model('./JAGS_ICBM_4.R',
                                      data=calib_data_Ultuna_long,
                                      n.chains = n.chains,
                                      n.adapt = N.ADAPTS)

update(model.ICBM.matrix_Ultuna_long, n.iter=N.RUNS, n.thin=10, n.burnin=N.BURNIN)
mcmc.array.ICBM.Ultuna_long<-(jags.samples(model.ICBM.matrix_Ultuna_long,parameter_list, sampling.nr))


###Ultuna
#Plot SOC and Input (roots and shoots) simulation
mcmc.list.SOC.Ultuna_long<-as.mcmc.list(mcmc.array.ICBM.Ultuna_long$Tot_Ultuna, chains=F)
mcmc.unlist.SOC.Ultuna_long<-mcmc.list.SOC.Ultuna_long[[1]]
for(i in 1:n.chains){
  mcmc.unlist.SOC.Ultuna_long<-rbind(mcmc.unlist.SOC.Ultuna_long, mcmc.list.SOC.Ultuna_long[[i]])
}

mcmc.list.Ir.Ultuna_long<-as.mcmc.list(mcmc.array.ICBM.Ultuna_long$I_R_Ultuna, chains=F)
mcmc.unlist.Ir.Ultuna_long<-(mcmc.list.Ir.Ultuna_long[[1]])
for(i in 1:n.chains){
  mcmc.unlist.Ir.Ultuna_long<-rbind(mcmc.unlist.Ir.Ultuna_long, mcmc.list.Ir.Ultuna_long[[i]])
}

mcmc.list.Is.Ultuna_long<-as.mcmc.list(mcmc.array.ICBM.Ultuna_long$I_S_Ultuna, chains=F)
mcmc.unlist.Is.Ultuna_long<-(mcmc.list.Is.Ultuna_long[[1]])
for(i in 1:n.chains){
  mcmc.unlist.Is.Ultuna_long<-rbind(mcmc.unlist.Is.Ultuna_long, mcmc.list.Is.Ultuna_long[[i]])
}

mcmc.list.fluxR.Ultuna_long<-as.mcmc.list(mcmc.array.ICBM.Ultuna_long$I_S_Ultuna, chains=F)
mcmc.unlist.fluxR.Ultuna_long<-(mcmc.list.fluxR.Ultuna_long[[1]])
for(i in 1:n.chains){
  mcmc.unlist.fluxR.Ultuna_long<-rbind(mcmc.unlist.fluxR.Ultuna_long, mcmc.list.fluxR.Ultuna_long[[i]])
}

mcmc.list.fluxS.Ultuna_long<-as.mcmc.list(mcmc.array.ICBM.Ultuna_long$I_S_Ultuna, chains=F)
mcmc.unlist.fluxS.Ultuna_long<-(mcmc.list.fluxS.Ultuna_long[[1]])
for(i in 1:n.chains){
  mcmc.unlist.fluxS.Ultuna_long<-rbind(mcmc.unlist.fluxS.Ultuna_long, mcmc.list.fluxS.Ultuna_long[[i]])
}

mcmc.list.fluxSUM.Ultuna_long<-as.mcmc.list(mcmc.array.ICBM.Ultuna_long$I_S_Ultuna, chains=F)
mcmc.unlist.fluxSUM.Ultuna_long<-(mcmc.list.fluxSUM.Ultuna_long[[1]])
for(i in 1:n.chains){
  mcmc.unlist.fluxSUM.Ultuna_long<-rbind(mcmc.unlist.fluxSUM.Ultuna_long, mcmc.list.fluxSUM.Ultuna_long[[i]])
}



#extract the array of the realizations*treatments*years
dim<-dim(calib_data_Ultuna_long$SOC_Ultuna)
dim[1]<-calib_data_Ultuna_long$J_Ultuna
Ultuna_prediction_array_long<-array( ,dim=c((sampling.nr*n.chains),dim))
str(Ultuna_prediction_array_long)
length_treat<-dim[1]
for(j in 1:(sampling.nr*n.chains)){
  for(i in 1:(dim[2])){
    length_treat_add<-length_treat*i
    Ultuna_prediction_array_long[j,,i]<-mcmc.unlist.SOC.Ultuna_long[j,(1-length_treat+length_treat_add):length_treat_add]
  }
}
plot(Ultuna_prediction_array_long[20,2,], type="l")


Ultuna_inputS_array_long<-array( ,dim=c((sampling.nr*n.chains),dim))
for(j in 1:(sampling.nr*n.chains)){
  length_treat<-dim[1]
  for(i in 1:(dim[2])){
    length_treat_add<-length_treat*i
    Ultuna_inputS_array_long[j,,i]<-mcmc.unlist.Is.Ultuna_long[j,(1-length_treat+length_treat_add):length_treat_add]
  }
}

Ultuna_inputR_array_long<-array( ,dim=c((sampling.nr*n.chains),dim))
for(j in 1:(sampling.nr*n.chains)){
  length_treat<-dim[1]
  for(i in 1:(dim[2])){
    length_treat_add<-length_treat*i
    Ultuna_inputR_array_long[j,,i]<-mcmc.unlist.Ir.Ultuna_long[j,(1-length_treat+length_treat_add):length_treat_add]
  }
}

Ultuna_fluxS_array_long<-array( ,dim=c((sampling.nr*n.chains),dim))
for(j in 1:(sampling.nr*n.chains)){
  length_treat<-dim[1]
  for(i in 1:(dim[2])){
    length_treat_add<-length_treat*i
    Ultuna_fluxS_array_long[j,,i]<-mcmc.unlist.fluxS.Ultuna_long[j,(1-length_treat+length_treat_add):length_treat_add]
  }
}

Ultuna_fluxR_array_long<-array( ,dim=c((sampling.nr*n.chains),dim))
for(j in 1:(sampling.nr*n.chains)){
  length_treat<-dim[1]
  for(i in 1:(dim[2])){
    length_treat_add<-length_treat*i
    Ultuna_fluxR_array_long[j,,i]<-mcmc.unlist.fluxR.Ultuna_long[j,(1-length_treat+length_treat_add):length_treat_add]
  }
}

Ultuna_fluxSUM_array_long<-array( ,dim=c((sampling.nr*n.chains),dim))
for(j in 1:(sampling.nr*n.chains)){
  length_treat<-dim[1]
  for(i in 1:(dim[2])){
    length_treat_add<-length_treat*i
    Ultuna_fluxSUM_array_long[j,,i]<-mcmc.unlist.fluxSUM.Ultuna_long[j,(1-length_treat+length_treat_add):length_treat_add]
  }
}

dim(mcmc.unlist.SOC.Ultuna)

str(Ultuna_prediction_array)

plot(Ultuna_prediction_array[2,2,], type="l")
lines(Ultuna_prediction_array[3,2,], type="l")




#mean, min and max matrices by plot
Ultuna_mean_predictions_long<-(colMeans(Ultuna_prediction_array_long))
Ultuna_max_predictions_long<-apply(Ultuna_prediction_array_long, MARGIN=c(2,3), max)
Ultuna_min_predictions_long<-apply(Ultuna_prediction_array_long, MARGIN=c(2,3), min)
rownames(Ultuna_mean_predictions_long)<-as.factor(rownames(calib_data_Ultuna_long$SOC_Ultuna))[1:calib_data_Ultuna$J_Ultuna]
rownames(Ultuna_max_predictions_long)<-as.factor(rownames(calib_data_Ultuna_long$SOC_Ultuna))[1:calib_data_Ultuna$J_Ultuna]
rownames(Ultuna_min_predictions_long)<-as.factor(rownames(calib_data_Ultuna_long$SOC_Ultuna))[1:calib_data_Ultuna$J_Ultuna]

Ultuna_mean_Ir_long<-(colMeans(Ultuna_inputR_array_long))
Ultuna_max_Ir_long<-apply(Ultuna_inputR_array_long, MARGIN=c(2,3), max)
Ultuna_min_Ir_long<-apply(Ultuna_inputR_array_long, MARGIN=c(2,3), min)
rownames(Ultuna_mean_Ir_long)<-as.factor(rownames(calib_data_Ultuna_long$SOC_Ultuna))[1:calib_data_Ultuna_long$J_Ultuna]
rownames(Ultuna_max_Ir_long)<-as.factor(rownames(calib_data_Ultuna_long$SOC_Ultuna))[1:calib_data_Ultuna_long$J_Ultuna]
rownames(Ultuna_min_Ir_long)<-as.factor(rownames(calib_data_Ultuna_long$SOC_Ultuna))[1:calib_data_Ultuna_long$J_Ultuna]

Ultuna_mean_Is_long<-(colMeans(Ultuna_inputS_array_long))
Ultuna_max_Is_long<-apply(Ultuna_inputS_array_long, MARGIN=c(2,3), max)
Ultuna_min_Is_long<-apply(Ultuna_inputS_array_long, MARGIN=c(2,3), min)
rownames(Ultuna_mean_Is_long)<-as.factor(rownames(calib_data_Ultuna_long$SOC_Ultuna))[1:calib_data_Ultuna_long$J_Ultuna]
rownames(Ultuna_max_Is_long)<-as.factor(rownames(calib_data_Ultuna_long$SOC_Ultuna))[1:calib_data_Ultuna_long$J_Ultuna]
rownames(Ultuna_min_Is_long)<-as.factor(rownames(calib_data_Ultuna_long$SOC_Ultuna))[1:calib_data_Ultuna_long$J_Ultuna]

Ultuna_mean_fluxS_long<-(colMeans(Ultuna_fluxS_array_long))
Ultuna_max_fluxS_long<-apply(Ultuna_fluxS_array_long, MARGIN=c(2,3), max)
Ultuna_min_fluxS_long<-apply(Ultuna_fluxS_array_long, MARGIN=c(2,3), min)
rownames(Ultuna_mean_fluxS_long)<-as.factor(rownames(calib_data_Ultuna_long$SOC_Ultuna))[1:calib_data_Ultuna_long$J_Ultuna]
rownames(Ultuna_max_fluxS_long)<-as.factor(rownames(calib_data_Ultuna_long$SOC_Ultuna))[1:calib_data_Ultuna_long$J_Ultuna]
rownames(Ultuna_min_fluxS_long)<-as.factor(rownames(calib_data_Ultuna_long$SOC_Ultuna))[1:calib_data_Ultuna_long$J_Ultuna]

Ultuna_mean_fluxR_long<-(colMeans(Ultuna_fluxR_array_long))
Ultuna_max_fluxR_long<-apply(Ultuna_fluxR_array_long, MARGIN=c(2,3), max)
Ultuna_min_fluxR_long<-apply(Ultuna_fluxR_array_long, MARGIN=c(2,3), min)
rownames(Ultuna_mean_fluxR_long)<-as.factor(rownames(calib_data_Ultuna_long$SOC_Ultuna))[1:calib_data_Ultuna_long$J_Ultuna]
rownames(Ultuna_max_fluxR_long)<-as.factor(rownames(calib_data_Ultuna_long$SOC_Ultuna))[1:calib_data_Ultuna_long$J_Ultuna]
rownames(Ultuna_min_fluxR_long)<-as.factor(rownames(calib_data_Ultuna_long$SOC_Ultuna))[1:calib_data_Ultuna_long$J_Ultuna]

Ultuna_mean_fluxSUM_long<-(colMeans(Ultuna_fluxSUM_array_long))
Ultuna_max_fluxSUM_long<-apply(Ultuna_fluxSUM_array_long, MARGIN=c(2,3), max)
Ultuna_min_fluxSUM_long<-apply(Ultuna_fluxSUM_array_long, MARGIN=c(2,3), min)
rownames(Ultuna_mean_fluxSUM_long)<-as.factor(rownames(calib_data_Ultuna_long$SOC_Ultuna))[1:calib_data_Ultuna_long$J_Ultuna]
rownames(Ultuna_max_fluxSUM_long)<-as.factor(rownames(calib_data_Ultuna_long$SOC_Ultuna))[1:calib_data_Ultuna_long$J_Ultuna]
rownames(Ultuna_min_fluxSUM_long)<-as.factor(rownames(calib_data_Ultuna_long$SOC_Ultuna))[1:calib_data_Ultuna_long$J_Ultuna]



#mean, min and max matrices by treatment
#Ultuna_treat_vec<-as.factor(rownames(calib_data_Ultuna_long$SOC_Ultuna))
Ultuna_mean_predictions_bytreat_long<-mat.or.vec(calib_data_Ultuna_long$J_Ultuna,dim[2])
Ultuna_max_predictions_bytreat_long<-mat.or.vec(calib_data_Ultuna_long$J_Ultuna,dim[2])
Ultuna_min_predictions_bytreat_long<-mat.or.vec(calib_data_Ultuna_long$J_Ultuna,dim[2])
Ultuna_mean_Ir_bytreat_long<-mat.or.vec(calib_data_Ultuna_long$J_Ultuna,dim[2])
Ultuna_max_Ir_bytreat_long<-mat.or.vec(calib_data_Ultuna_long$J_Ultuna,dim[2])
Ultuna_min_Ir_bytreat_long<-mat.or.vec(calib_data_Ultuna_long$J_Ultuna,dim[2])
Ultuna_mean_Is_bytreat_long<-mat.or.vec(calib_data_Ultuna_long$J_Ultuna,dim[2])
Ultuna_max_Is_bytreat_long<-mat.or.vec(calib_data_Ultuna_long$J_Ultuna,dim[2])
Ultuna_min_Is_bytreat_long<-mat.or.vec(calib_data_Ultuna_long$J_Ultuna,dim[2])
Ultuna_mean_fluxS_bytreat_long<-mat.or.vec(calib_data_Ultuna_long$J_Ultuna,dim[2])
Ultuna_max_fluxS_bytreat_long<-mat.or.vec(calib_data_Ultuna_long$J_Ultuna,dim[2])
Ultuna_min_fluxS_bytreat_long<-mat.or.vec(calib_data_Ultuna_long$J_Ultuna,dim[2])
Ultuna_mean_fluxR_bytreat_long<-mat.or.vec(calib_data_Ultuna_long$J_Ultuna,dim[2])
Ultuna_max_fluxR_bytreat_long<-mat.or.vec(calib_data_Ultuna_long$J_Ultuna,dim[2])
Ultuna_min_fluxR_bytreat_long<-mat.or.vec(calib_data_Ultuna_long$J_Ultuna,dim[2])
Ultuna_mean_fluxSUM_bytreat_long<-mat.or.vec(calib_data_Ultuna_long$J_Ultuna,dim[2])
Ultuna_max_fluxSUM_bytreat_long<-mat.or.vec(calib_data_Ultuna_long$J_Ultuna,dim[2])
Ultuna_min_fluxSUM_bytreat_long<-mat.or.vec(calib_data_Ultuna_long$J_Ultuna,dim[2])
Ultuna_measured_bytreat_long<-mat.or.vec(calib_data_Ultuna_long$J_Ultuna,dim[2])
for(i in 1:calib_data_Ultuna_long$J_Ultuna){
  Ultuna_mean_predictions_bytreat_long[i,]<-Ultuna_mean_predictions_long[i,]
  Ultuna_min_predictions_bytreat_long[i,]<-Ultuna_min_predictions_long[i,]
  Ultuna_max_predictions_bytreat_long[i,]<-Ultuna_max_predictions_long[i,]
  Ultuna_mean_Ir_bytreat_long[i,]<-Ultuna_mean_Ir_long[i,]
  Ultuna_min_Ir_bytreat_long[i,]<-Ultuna_min_Ir_long[i,]
  Ultuna_max_Ir_bytreat_long[i,]<-Ultuna_max_Ir_long[i,]
  Ultuna_mean_Is_bytreat_long[i,]<-Ultuna_mean_Is_long[i,]
  Ultuna_min_Is_bytreat_long[i,]<-Ultuna_min_Is_long[i,]
  Ultuna_max_Is_bytreat_long[i,]<-Ultuna_max_Is_long[i,]
  Ultuna_mean_fluxS_bytreat_long[i,]<-Ultuna_mean_fluxS_long[i,]
  Ultuna_min_fluxS_bytreat_long[i,]<-Ultuna_min_fluxS_long[i,]
  Ultuna_max_fluxS_bytreat_long[i,]<-Ultuna_max_fluxS_long[i,]
  Ultuna_mean_fluxR_bytreat_long[i,]<-Ultuna_mean_fluxR_long[i,]
  Ultuna_min_fluxR_bytreat_long[i,]<-Ultuna_min_fluxR_long[i,]
  Ultuna_max_fluxR_bytreat_long[i,]<-Ultuna_max_fluxR_long[i,]
  Ultuna_mean_fluxSUM_bytreat_long[i,]<-Ultuna_mean_fluxSUM_long[i,]
  Ultuna_min_fluxSUM_bytreat_long[i,]<-Ultuna_min_fluxSUM_long[i,]
  Ultuna_max_fluxSUM_bytreat_long[i,]<-Ultuna_max_fluxSUM_long[i,]
  Ultuna_measured_bytreat_long[i,]<-as.numeric(calib_data_Ultuna_long$SOC_Ultuna[i,])
}


write.csv(t(Ultuna_mean_Ir_bytreat_long), file="Ultuna_mean_Ir_bytreat_long.csv")
write.csv(t(Ultuna_mean_Is_bytreat_long), file="Ultuna_mean_Is_bytreat_long.csv")

Ultuna_yields_table_long<-calib_data_Ultuna_long$Yields_cereals_Ultuna+calib_data_Ultuna_long$Yields_root_crops_Ultuna+calib_data_Ultuna_long$Yields_oilseeds_Ultuna+calib_data_Ultuna_long$Yields_maize_Ultuna
Ultuna_Ir_table_long<-as.data.frame(Ultuna_mean_Ir_bytreat_long)
Ultuna_Is_table_long<-as.data.frame(Ultuna_mean_Is_bytreat_long)
Ultuna_Pred_table_long<-as.data.frame(Ultuna_max_predictions_bytreat_long)

Ultuna_Is_table_long<-rbind(Ultuna_Is_table_long, seq(from=1956, to=(1956+71)))
Ultuna_Ir_table_long<-rbind(Ultuna_Ir_table_long, seq(from=1956, to=(1956+71)))
Ultuna_yields_table_long<-rbind(Ultuna_yields_table_long, seq(from=1956, to=(1956+71)))
Ultuna_Pred_table_long<-rbind(Ultuna_Pred_table_long, seq(from=1956, to=(1956+71)))

rownames(Ultuna_Ir_table_long)[16]<-"Year"
rownames(Ultuna_yields_table_long)[16]<-"Year"
rownames(Ultuna_Is_table_long)[16]<-"Year"
rownames(Ultuna_Pred_table_long)[16]<-"Year"


write_ods(Ultuna_yields_table_long, "Ultuna_input_check_long.ods", sheet="Yields", row_names = T)
write_ods(Ultuna_Ir_table_long, "Ultuna_input_check_long.ods", sheet="Calculated root inputs (C,ton ha-1)", append=T, row_names = T)
write_ods(Ultuna_Is_table_long, "Ultuna_input_check_long.ods", sheet="Calculated shoot inputs (C,ton ha-1)", append=T, row_names = T)
write_ods(Ultuna_Pred_table_long, "Ultuna_input_check_long.ods", sheet="Predicted SOC stocks (C,ton ha-1)", append=T, row_names = T)




png("ICBM_predictions_Ultuna_specific_long.png", height=4500, width=4000, res=320)
par(mfrow=c(5,3))
for(i in 1:calib_data_Ultuna_long$J_Ultuna){
  yearseq<-seq(from=colnames(calib_data_Ultuna_long$SOC_Ultuna)[1], to=as.numeric(colnames(calib_data_Ultuna_long$SOC_Ultuna)[1])+dim[2]-1)
  lastyear<-dim(Ultuna_mean_predictions_bytreat_long)[2]-1
  plot(yearseq, Ultuna_mean_predictions_bytreat_long[i,], ylim=c(20,125), type="l",  ylab=expression(paste("C (t ha" ^-1, ")")), xlab="Years", main=paste("Treatment", LETTERS[selection_treatments[i]]), col=Ultuna_treat_palette[i])
  polygon(c(yearseq,rev(yearseq)),
          c(Ultuna_max_predictions_bytreat_long[i,],rev(Ultuna_min_predictions_bytreat_long[i,])),
          col=Ultuna_treat_palette_alpha[i],border=Ultuna_treat_palette[i])
  lines(yearseq, Ultuna_mean_predictions_bytreat_long[i,], col="black", lty=1, lwd=1)
  points(yearseq, (Ultuna_measured_bytreat_long[i,]), col="firebrick", lty=2, lwd=2)
}
dev.off()



png("ICBM_predictions_Ultuna_specific_freescale_long.png", height=4500, width=4000, res=320)
par(mfrow=c(5,3))
for(i in 1:calib_data_Ultuna_long$J_Ultuna){
  yearseq<-seq(from=colnames(calib_data_Ultuna_long$SOC_Ultuna)[1], to=as.numeric(colnames(calib_data_Ultuna_long$SOC_Ultuna)[1])+dim[2]-1)
  lastyear<-dim(Ultuna_mean_predictions_bytreat_long)[2]-1
  plot(yearseq, Ultuna_mean_predictions_bytreat_long[i,], ylim=c(min(Ultuna_mean_predictions_bytreat_long[i,])*0.3,max(Ultuna_mean_predictions_bytreat_long[i,])*1.6), type="l",  ylab=expression(paste("C (t ha" ^-1, ")")), xlab="Years", main=paste("Treatment", LETTERS[selection_treatments[i]]), col=Ultuna_treat_palette[i])
  polygon(c(yearseq,rev(yearseq)),
          c(Ultuna_max_predictions_bytreat_long[i,],rev(Ultuna_min_predictions_bytreat_long[i,])),
          col=Ultuna_treat_palette_alpha[i],border=Ultuna_treat_palette[i])
  lines(yearseq, Ultuna_mean_predictions_bytreat_long[i,], col="black", lty=1, lwd=1)
  points(yearseq, (Ultuna_measured_bytreat_long[i,]), col="firebrick", lty=2, lwd=2)
}
dev.off()



png("ICBM_Ir_Ultuna_specific_long.png", height=4500, width=4000, res=320)
par(mfrow=c(5,3))
for(i in 1:calib_data_Ultuna$J_Ultuna){
  yearseq<-seq(from=colnames(calib_data_Ultuna_long$SOC_Ultuna)[1], to=as.numeric(colnames(calib_data_Ultuna_long$SOC_Ultuna)[1])+dim[2]-1)
  lastyear<-dim(Ultuna_mean_Ir_bytreat_long)[2]-1
  plot(yearseq, Ultuna_mean_Ir_bytreat_long[i,], ylim=c(0,2), type="l",  ylab=expression(paste("C (t ha" ^-1, ")")), xlab="Years", main=paste("Treatment", LETTERS[selection_treatments[i]]), col=Ultuna_treat_palette[i])
  polygon(c(yearseq,rev(yearseq)),
          c(Ultuna_max_Ir_bytreat_long[i,],rev(Ultuna_min_Ir_bytreat_long[i,])),
          col=Ultuna_treat_palette_alpha[i],border=Ultuna_treat_palette[i])
  lines(yearseq, Ultuna_mean_Ir_bytreat_long[i,], col="black", lty=1, lwd=1)
}
dev.off()

png("ICBM_Is_Ultuna_specific_long.png", height=4500, width=4000, res=320)
par(mfrow=c(5,3))
for(i in 1:calib_data_Ultuna_long$J_Ultuna){
  yearseq<-seq(from=colnames(calib_data_Ultuna_long$SOC_Ultuna)[1], to=as.numeric(colnames(calib_data_Ultuna$SOC_Ultuna)[1])+dim[2]-1)
  lastyear<-dim(Ultuna_mean_Is_bytreat)[2]-1
  plot(yearseq, Ultuna_mean_Is_bytreat_long[i,], ylim=c(0,2), type="l",  ylab=expression(paste("C (t ha" ^-1, ")")), xlab="Years", main=paste("Treatment", LETTERS[selection_treatments[i]]), col=Ultuna_treat_palette[i])
  polygon(c(yearseq,rev(yearseq)),
          c(Ultuna_max_Is_bytreat_long[i,],rev(Ultuna_min_Is_bytreat_long[i,])),
          col=Ultuna_treat_palette_alpha[i],border=Ultuna_treat_palette[i])
  lines(yearseq, Ultuna_mean_Is_bytreat_long[i,], col="black", lty=1, lwd=1)
}
dev.off()

png("ICBM_fluxS_Ultuna_specific_long.png", height=4500, width=4000, res=320)
par(mfrow=c(5,3))
for(i in 1:calib_data_Ultuna_long$J_Ultuna){
  yearseq<-seq(from=colnames(calib_data_Ultuna_long$SOC_Ultuna)[1], to=as.numeric(colnames(calib_data_Ultuna_long$SOC_Ultuna)[1])+dim[2]-1)
  lastyear<-dim(Ultuna_mean_fluxS_bytreat_long)[2]-1
  plot(yearseq, Ultuna_mean_fluxS_bytreat_long[i,], ylim=c(0,2), type="l",  ylab=expression(paste("C (t ha" ^-1, ")")), xlab="Years", main=paste("Treatment", LETTERS[selection_treatments[i]]), col=Ultuna_treat_palette[i])
  polygon(c(yearseq,rev(yearseq)),
          c(Ultuna_max_fluxS_bytreat_long[i,],rev(Ultuna_min_fluxS_bytreat_long[i,])),
          col=Ultuna_treat_palette_alpha[i],border=Ultuna_treat_palette[i])
  lines(yearseq, Ultuna_mean_fluxS_bytreat_long[i,], col="black", lty=1, lwd=1)
}
dev.off()

png("ICBM_fluxR_Ultuna_specific_long.png", height=4500, width=4000, res=320)
par(mfrow=c(5,3))
for(i in 1:calib_data_Ultuna$J_Ultuna){
  yearseq<-seq(from=colnames(calib_data_Ultuna_long$SOC_Ultuna)[1], to=as.numeric(colnames(calib_data_Ultuna_long$SOC_Ultuna)[1])+dim[2]-1)
  lastyear<-dim(Ultuna_mean_fluxR_bytreat_long)[2]-1
  plot(yearseq, Ultuna_mean_fluxR_bytreat_long[i,], ylim=c(0.4,-0.5), type="l",  ylab=expression(paste("C (t ha" ^-1, ")")), xlab="Years", main=paste("Treatment", LETTERS[selection_treatments[i]]), col=Ultuna_treat_palette[i])
  polygon(c(yearseq,rev(yearseq)),
          c(Ultuna_max_fluxR_bytreat_long[i,],rev(Ultuna_min_fluxR_bytreat_long[i,])),
          col=Ultuna_treat_palette_alpha[i],border=Ultuna_treat_palette[i])
  lines(yearseq, Ultuna_mean_fluxR_bytreat_long[i,], col="black", lty=1, lwd=1)
}
dev.off()

png("ICBM_fluxSUM_Ultuna_specific.png", height=4500, width=4000, res=320)
par(mfrow=c(5,3))
for(i in 1:calib_data_Ultuna_long$J_Ultuna){
  yearseq<-seq(from=colnames(calib_data_Ultuna_long$SOC_Ultuna)[1], to=as.numeric(colnames(calib_data_Ultuna_long$SOC_Ultuna)[1])+dim[2]-1)
  lastyear<-dim(Ultuna_mean_fluxSUM_bytreat_long)[2]-1
  plot(yearseq, Ultuna_mean_fluxSUM_bytreat_long[i,], ylim=c(0.4,-0.5), type="l",  ylab=expression(paste("C (t ha" ^-1, ")")), xlab="Years", main=paste("Treatment", LETTERS[selection_treatments[i]]), col=Ultuna_treat_palette[i])
  polygon(c(yearseq,rev(yearseq)),
          c(Ultuna_max_fluxSUM_bytreat_long[i,],rev(Ultuna_min_fluxSUM_bytreat_long[i,])),
          col=Ultuna_treat_palette_alpha[i],border=Ultuna_treat_palette[i])
  lines(yearseq, Ultuna_mean_fluxSUM_bytreat_long[i,], col="black", lty=1, lwd=1)
}
dev.off()











############### ############### ############### ############### ############### 
#Input calculation comparison between Thomas and me, 1999
############### ############### ############### ############### ############### 


Ultuna_Thomas_Ir_table<-data.frame(mat.or.vec(15,54))
for(i in 1:calib_data_Ultuna$J_Ultuna){
  Ultuna_Thomas_Ir_table[i,]<-as.numeric(Thomas_2011_inputs[Thomas_2011_inputs$treat==letters[i],]$`CR+CE70%`[1:54])*0.001
}

colnames(Ultuna_Thomas_Ir_table)<-colnames(calib_data_Ultuna$SOC_Ultuna)[1:54]
rownames(Ultuna_Thomas_Ir_table)<-LETTERS[1:calib_data_Ultuna$J_Ultuna]

write_ods(Ultuna_Thomas_Ir_table, "Ultuna_input_check.ods", sheet="Thomas calculations IR (C,ton ha-1)", append=T, row_names = T)


png("input_comparison.png", width=4000, height=4000, res=300)
par(mfrow=c(4,4))
for(i in 1:calib_data_Ultuna$J_Ultuna){
  plot(Ultuna_mean_Ir_bytreat[i,][1:54], type="l", col="darkorange", lwd=2, ylab="C Mg Ha-1", xlab="Year", main=paste("treatment", letters[i]), ylim=c(0,3))
  lines(as.numeric(Thomas_2011_inputs[Thomas_2011_inputs$treat==letters[i],]$`CR+CE70%`[1:54])*0.001, col="darkgreen", lwd=2)
  if(i ==1){legend("topleft", c("me", "Thomas"), col=c("darkorange", "darkgreen"), bty="n", lty=c(1), cex=1.5)}
  }
dev.off()

png("input_difference_me_vs_Thomas.png", width=4000, height=4000, res=300)
par(mfrow=c(4,4))
for(i in 1:calib_data_Ultuna$J_Ultuna){
  barplot(Ultuna_mean_Ir_bytreat[i,][1:54]-as.numeric(Thomas_2011_inputs[Thomas_2011_inputs$treat==letters[i],]$`CR+CE70%`[1:54])*0.001, ylab="C Mg Ha-1", xlab="Year", main=paste("treatment", letters[i]), ylim=c(-2,2))
  text(9, 1.7, paste("mean : ", round(mean(Ultuna_mean_Ir_bytreat[i,][1:54]-as.numeric(Thomas_2011_inputs[Thomas_2011_inputs$treat==letters[i],]$`CR+CE70%`[1:54])*0.001, na.rm=T),3)))
}
dev.off()


png("input_difference_me_vs_Thomas_procent.png", width=4000, height=4000, res=300)
par(mfrow=c(4,4))
for(i in 1:calib_data_Ultuna$J_Ultuna){
  barplot(100*(Ultuna_mean_Ir_bytreat[i,][1:54]-as.numeric(Thomas_2011_inputs[Thomas_2011_inputs$treat==letters[i],]$`CR+CE70%`[1:54])*0.001)/mean(Ultuna_mean_Ir_bytreat[i,][1:54]), ylab="% difference", xlab="Year", main=paste("treatment", letters[i]), ylim=c(-200,200))
  text(18, 185, paste("mean difference : ", round(mean(Ultuna_mean_Ir_bytreat[i,][1:54]-as.numeric(Thomas_2011_inputs[Thomas_2011_inputs$treat==letters[i],]$`CR+CE70%`[1:54])*0.001, na.rm=T)/mean(Ultuna_mean_Ir_bytreat[i,][1:54]),2)*100, "%"))
}
dev.off()





png("ICBM_predictions_Ultuna_specific.png", height=4500, width=4000, res=320)
par(mfrow=c(5,3))
for(i in 1:calib_data_Ultuna$J_Ultuna){
  yearseq<-seq(from=colnames(calib_data_Ultuna$SOC_Ultuna)[1], to=as.numeric(colnames(calib_data_Ultuna$SOC_Ultuna)[1])+dim[2]-1)
  lastyear<-dim(Ultuna_mean_predictions_bytreat)[2]-1
  plot(yearseq, Ultuna_mean_predictions_bytreat[i,], ylim=c(20,125), type="l",  ylab=expression(paste("C (t ha" ^-1, ")")), xlab="Years", main=paste("Treatment", LETTERS[selection_treatments[i]]), col=Ultuna_treat_palette[i])
  polygon(c(yearseq,rev(yearseq)),
          c(Ultuna_max_predictions_bytreat[i,],rev(Ultuna_min_predictions_bytreat[i,])),
          col=Ultuna_treat_palette_alpha[i],border=Ultuna_treat_palette[i])
  lines(yearseq, Ultuna_mean_predictions_bytreat[i,], col="black", lty=1, lwd=1)
  #lines(yearseq[2:62], na_interpolation(Ultuna_measured_bytreat[i,1:61]), col="firebrick", lty=2, lwd=2)
  #lines(yearseq, na_interpolation(Ultuna_measured_bytreat[i,1:53]), col=add.alpha("firebrick",0.2), lty=2, lwd=2)
  points(yearseq, (Ultuna_measured_bytreat[i,]), col="firebrick", lty=2, lwd=2)
}
dev.off()



png("ICBM_predictions_Ultuna_specific_freescale.png", height=4500, width=4000, res=320)
par(mfrow=c(5,3))
for(i in 1:calib_data_Ultuna$J_Ultuna){
  yearseq<-seq(from=colnames(calib_data_Ultuna$SOC_Ultuna)[1], to=as.numeric(colnames(calib_data_Ultuna$SOC_Ultuna)[1])+dim[2]-1)
  lastyear<-dim(Ultuna_mean_predictions_bytreat)[2]-1
  plot(yearseq, Ultuna_mean_predictions_bytreat[i,], ylim=c(min(Ultuna_mean_predictions_bytreat[i,])*0.3,max(Ultuna_mean_predictions_bytreat[i,])*1.6), type="l",  ylab=expression(paste("C (t ha" ^-1, ")")), xlab="Years", main=paste("Treatment", LETTERS[selection_treatments[i]]), col=Ultuna_treat_palette[i])
  polygon(c(yearseq,rev(yearseq)),
          c(Ultuna_max_predictions_bytreat[i,],rev(Ultuna_min_predictions_bytreat[i,])),
          col=Ultuna_treat_palette_alpha[i],border=Ultuna_treat_palette[i])
  lines(yearseq, Ultuna_mean_predictions_bytreat[i,], col="black", lty=1, lwd=1)
  #lines(yearseq[2:62], na_interpolation(Ultuna_measured_bytreat[i,1:61]), col="firebrick", lty=2, lwd=2)
  #lines(yearseq, na_interpolation(Ultuna_measured_bytreat[i,1:53]), col=add.alpha("firebrick",0.2), lty=2, lwd=2)
  points(yearseq, (Ultuna_measured_bytreat[i,]), col="firebrick", lty=2, lwd=2)
}
dev.off()



png("ICBM_Ir_Ultuna_specific.png", height=4500, width=4000, res=320)
par(mfrow=c(5,3))
for(i in 1:calib_data_Ultuna$J_Ultuna){
  yearseq<-seq(from=colnames(calib_data_Ultuna$SOC_Ultuna)[1], to=as.numeric(colnames(calib_data_Ultuna$SOC_Ultuna)[1])+dim[2]-1)
  lastyear<-dim(Ultuna_mean_Ir_bytreat)[2]-1
  plot(yearseq, Ultuna_mean_Ir_bytreat[i,], ylim=c(0,2), type="l",  ylab=expression(paste("C (t ha" ^-1, ")")), xlab="Years", main=paste("Treatment", LETTERS[selection_treatments[i]]), col=Ultuna_treat_palette[i])
  polygon(c(yearseq,rev(yearseq)),
          c(Ultuna_max_Ir_bytreat[i,],rev(Ultuna_min_Ir_bytreat[i,])),
          col=Ultuna_treat_palette_alpha[i],border=Ultuna_treat_palette[i])
  lines(yearseq, Ultuna_mean_Ir_bytreat[i,], col="black", lty=1, lwd=1)
}
dev.off()

png("ICBM_Is_Ultuna_specific.png", height=4500, width=4000, res=320)
par(mfrow=c(5,3))
for(i in 1:calib_data_Ultuna$J_Ultuna){
  yearseq<-seq(from=colnames(calib_data_Ultuna$SOC_Ultuna)[1], to=as.numeric(colnames(calib_data_Ultuna$SOC_Ultuna)[1])+dim[2]-1)
  lastyear<-dim(Ultuna_mean_Is_bytreat)[2]-1
  plot(yearseq, Ultuna_mean_Is_bytreat[i,], ylim=c(0,2), type="l",  ylab=expression(paste("C (t ha" ^-1, ")")), xlab="Years", main=paste("Treatment", LETTERS[selection_treatments[i]]), col=Ultuna_treat_palette[i])
  polygon(c(yearseq,rev(yearseq)),
          c(Ultuna_max_Is_bytreat[i,],rev(Ultuna_min_Is_bytreat[i,])),
          col=Ultuna_treat_palette_alpha[i],border=Ultuna_treat_palette[i])
  lines(yearseq, Ultuna_mean_Is_bytreat[i,], col="black", lty=1, lwd=1)
}
dev.off()

png("ICBM_fluxS_Ultuna_specific.png", height=4500, width=4000, res=320)
par(mfrow=c(5,3))
for(i in 1:calib_data_Ultuna$J_Ultuna){
  yearseq<-seq(from=colnames(calib_data_Ultuna$SOC_Ultuna)[1], to=as.numeric(colnames(calib_data_Ultuna$SOC_Ultuna)[1])+dim[2]-1)
  lastyear<-dim(Ultuna_mean_fluxS_bytreat)[2]-1
  plot(yearseq, Ultuna_mean_fluxS_bytreat[i,], ylim=c(0,2), type="l",  ylab=expression(paste("C (t ha" ^-1, ")")), xlab="Years", main=paste("Treatment", LETTERS[selection_treatments[i]]), col=Ultuna_treat_palette[i])
  polygon(c(yearseq,rev(yearseq)),
          c(Ultuna_max_fluxS_bytreat[i,],rev(Ultuna_min_fluxS_bytreat[i,])),
          col=Ultuna_treat_palette_alpha[i],border=Ultuna_treat_palette[i])
  lines(yearseq, Ultuna_mean_fluxS_bytreat[i,], col="black", lty=1, lwd=1)
}
dev.off()

png("ICBM_fluxR_Ultuna_specific.png", height=4500, width=4000, res=320)
par(mfrow=c(5,3))
for(i in 1:calib_data_Ultuna$J_Ultuna){
  yearseq<-seq(from=colnames(calib_data_Ultuna$SOC_Ultuna)[1], to=as.numeric(colnames(calib_data_Ultuna$SOC_Ultuna)[1])+dim[2]-1)
  lastyear<-dim(Ultuna_mean_fluxR_bytreat)[2]-1
  plot(yearseq, Ultuna_mean_fluxR_bytreat[i,], ylim=c(0.4,-0.5), type="l",  ylab=expression(paste("C (t ha" ^-1, ")")), xlab="Years", main=paste("Treatment", LETTERS[selection_treatments[i]]), col=Ultuna_treat_palette[i])
  polygon(c(yearseq,rev(yearseq)),
          c(Ultuna_max_fluxR_bytreat[i,],rev(Ultuna_min_fluxR_bytreat[i,])),
          col=Ultuna_treat_palette_alpha[i],border=Ultuna_treat_palette[i])
  lines(yearseq, Ultuna_mean_fluxR_bytreat[i,], col="black", lty=1, lwd=1)
}
dev.off()

png("ICBM_fluxSUM_Ultuna_specific.png", height=4500, width=4000, res=320)
par(mfrow=c(5,3))
for(i in 1:calib_data_Ultuna$J_Ultuna){
  yearseq<-seq(from=colnames(calib_data_Ultuna$SOC_Ultuna)[1], to=as.numeric(colnames(calib_data_Ultuna$SOC_Ultuna)[1])+dim[2]-1)
  lastyear<-dim(Ultuna_mean_fluxSUM_bytreat)[2]-1
  plot(yearseq, Ultuna_mean_fluxSUM_bytreat[i,], ylim=c(0.4,-0.5), type="l",  ylab=expression(paste("C (t ha" ^-1, ")")), xlab="Years", main=paste("Treatment", LETTERS[selection_treatments[i]]), col=Ultuna_treat_palette[i])
  polygon(c(yearseq,rev(yearseq)),
          c(Ultuna_max_fluxSUM_bytreat[i,],rev(Ultuna_min_fluxSUM_bytreat[i,])),
          col=Ultuna_treat_palette_alpha[i],border=Ultuna_treat_palette[i])
  lines(yearseq, Ultuna_mean_fluxSUM_bytreat[i,], col="black", lty=1, lwd=1)
}
dev.off()




#Plot Y and O simulation
mcmc.list.Y.Ultuna<-as.mcmc.list(mcmc.array.ICBM.Ultuna$Y_tot_Ultuna, chains=T)
mcmc.unlist.Y.Ultuna<-mcmc.list.Y.Ultuna[[1]]
dim(mcmc.unlist.Y.Ultuna)

#extract the array of the realizations*treatments*years
Ultuna_prediction_array<-array( ,dim=c(sampling.nr,dim))
str(Ultuna_prediction_array)
for(j in 1:sampling.nr){
  length_treat<-dim[1]
  for(i in 1:(dim[2])){
    length_treat_add<-length_treat*i
    Ultuna_prediction_array[j,,i]<-mcmc.unlist.Y.Ultuna[j,(1-length_treat+length_treat_add):length_treat_add]
  }
}
#mean, min and max matrices by plot
Ultuna_mean_predictions<-(colMeans(Ultuna_prediction_array))
Ultuna_max_predictions<-colmax<-apply(Ultuna_prediction_array, MARGIN=c(2,3), max)
Ultuna_min_predictions<-colmax<-apply(Ultuna_prediction_array, MARGIN=c(2,3), min)
rownames(Ultuna_mean_predictions)<-LETTERS[1:calib_data_Ultuna$J_Ultuna]
rownames(Ultuna_max_predictions)<-LETTERS[1:calib_data_Ultuna$J_Ultuna]
rownames(Ultuna_min_predictions)<-LETTERS[1:calib_data_Ultuna$J_Ultuna]

#mean, min and max matrices by treatment
Ultuna_mean_predictions_bytreat_Y<-mat.or.vec(calib_data_Ultuna$J_Ultuna,dim[2])
Ultuna_max_predictions_bytreat_Y<-mat.or.vec(calib_data_Ultuna$J_Ultuna,dim[2])
Ultuna_min_predictions_bytreat_Y<-mat.or.vec(calib_data_Ultuna$J_Ultuna,dim[2])
for(i in 1:calib_data_Ultuna$J_Ultuna){
  Ultuna_mean_predictions_bytreat_Y[i,]<-Ultuna_mean_predictions[i,]
  Ultuna_max_predictions_bytreat_Y[i,]<-Ultuna_min_predictions[i,]
  Ultuna_min_predictions_bytreat_Y[i,]<-Ultuna_max_predictions[i,]
}

mcmc.list.O.Ultuna<-as.mcmc.list(mcmc.array.ICBM.Ultuna$O_Ultuna, chains=T)
mcmc.unlist.O.Ultuna<-mcmc.list.O.Ultuna[[1]]

#extract the array of the realizations*treatments*years
Ultuna_prediction_array<-array( ,dim=c(sampling.nr,dim))
str(Ultuna_prediction_array)
for(j in 1:sampling.nr){
  length_treat<-dim[1]
  for(i in 1:(dim[2])){
    length_treat_add<-length_treat*i
    Ultuna_prediction_array[j,,i]<-mcmc.unlist.O.Ultuna[j,(1-length_treat+length_treat_add):length_treat_add]
  }
}
#mean, min and max matrices by plot
Ultuna_mean_predictions<-(colMeans(Ultuna_prediction_array))
Ultuna_max_predictions<-colmax<-apply(Ultuna_prediction_array, MARGIN=c(2,3), max)
Ultuna_min_predictions<-colmax<-apply(Ultuna_prediction_array, MARGIN=c(2,3), min)
rownames(Ultuna_mean_predictions)<-LETTERS[1:calib_data_Ultuna$J_Ultuna]
rownames(Ultuna_max_predictions)<-LETTERS[1:calib_data_Ultuna$J_Ultuna]
rownames(Ultuna_min_predictions)<-LETTERS[1:calib_data_Ultuna$J_Ultuna]

#mean, min and max matrices by treatment
Ultuna_mean_predictions_bytreat_O<-mat.or.vec(calib_data_Ultuna$J_Ultuna,dim[2])
Ultuna_max_predictions_bytreat_O<-mat.or.vec(calib_data_Ultuna$J_Ultuna,dim[2])
Ultuna_min_predictions_bytreat_O<-mat.or.vec(calib_data_Ultuna$J_Ultuna,dim[2])
for(i in 1:calib_data_Ultuna$J_Ultuna){
  Ultuna_mean_predictions_bytreat_O[i,]<-Ultuna_mean_predictions[i,]
  Ultuna_max_predictions_bytreat_O[i,]<-Ultuna_min_predictions[i,]
  Ultuna_min_predictions_bytreat_O[i,]<-Ultuna_max_predictions[i,]
}

png("ICBM_predictions_Ultuna_O_specific.png", height=4500, width=4000, res=320)
par(mfrow=c(5,3))
for(i in 1:calib_data_Ultuna$J_Ultuna){
  yearseq<-seq(from=colnames(calib_data_Ultuna$SOC_Ultuna)[1], to=as.numeric(colnames(calib_data_Ultuna$SOC_Ultuna)[1])+dim[2]-1)
  lastyear<-dim(Ultuna_mean_predictions_bytreat_O)[2]-1
  plot(yearseq, Ultuna_mean_predictions_bytreat_O[i,], ylim=c(20,125), type="l",  ylab=expression(paste("C (t ha" ^-1, ")")), xlab="Years", main=paste("Treatment", LETTERS[selection_treatments[i]]), col=Ultuna_treat_palette[i])
  polygon(c(yearseq,rev(yearseq)),
          c(Ultuna_min_predictions_bytreat_O[i,],rev(Ultuna_max_predictions_bytreat_O[i,])),
          col=Ultuna_treat_palette_alpha[i],border=Ultuna_treat_palette[i])
  lines(yearseq, Ultuna_mean_predictions_bytreat_O[i,], col="black", lty=1, lwd=1)
}
dev.off()

png("ICBM_predictions_Ultuna_Y_specific.png", height=4500, width=4000, res=320)
par(mfrow=c(5,3))
for(i in 1:calib_data_Ultuna$J_Ultuna){
  yearseq<-seq(from=colnames(calib_data_Ultuna$SOC_Ultuna)[1], to=as.numeric(colnames(calib_data_Ultuna$SOC_Ultuna)[1])+dim[2]-1)
  lastyear<-dim(Ultuna_mean_predictions_bytreat_Y)[2]-1
  plot(yearseq, Ultuna_mean_predictions_bytreat_Y[i,], ylim=c(0,6), type="l",  ylab=expression(paste("C (t ha" ^-1, ")")), xlab="Years", main=paste("Treatment", LETTERS[selection_treatments[i]]), col=Ultuna_treat_palette[i])
  polygon(c(yearseq,rev(yearseq)),
          c(Ultuna_min_predictions_bytreat_Y[i,],rev(Ultuna_max_predictions_bytreat_Y[i,])),
          col=Ultuna_treat_palette_alpha[i],border=Ultuna_treat_palette[i])
  lines(yearseq, Ultuna_mean_predictions_bytreat_Y[i,], col="black", lty=1, lwd=1)
}
dev.off()





#Plot Y and O simulation
mcmc.list.Y.Ultuna<-as.mcmc.list(mcmc.array.ICBM.Ultuna$I_R_oilseeds_Ultuna, chains=T)
mcmc.unlist.Y.Ultuna<-mcmc.list.Y.Ultuna[[1]]
dim(mcmc.unlist.Y.Ultuna)

#extract the array of the realizations*treatments*years
Ultuna_prediction_array<-array( ,dim=c(sampling.nr,dim))
str(Ultuna_prediction_array)
for(j in 1:sampling.nr){
  length_treat<-dim[1]
  for(i in 1:(dim[2])){
    length_treat_add<-length_treat*i
    Ultuna_prediction_array[j,,i]<-mcmc.unlist.Y.Ultuna[j,(1-length_treat+length_treat_add):length_treat_add]
  }
}
#mean, min and max matrices by plot
Ultuna_mean_predictions<-(colMeans(Ultuna_prediction_array))
Ultuna_max_predictions<-colmax<-apply(Ultuna_prediction_array, MARGIN=c(2,3), max)
Ultuna_min_predictions<-colmax<-apply(Ultuna_prediction_array, MARGIN=c(2,3), min)
rownames(Ultuna_mean_predictions)<-LETTERS[1:calib_data_Ultuna$J_Ultuna]
rownames(Ultuna_max_predictions)<-LETTERS[1:calib_data_Ultuna$J_Ultuna]
rownames(Ultuna_min_predictions)<-LETTERS[1:calib_data_Ultuna$J_Ultuna]

#mean, min and max matrices by treatment
Ultuna_mean_predictions_bytreat_Y<-mat.or.vec(calib_data_Ultuna$J_Ultuna,dim[2])
Ultuna_max_predictions_bytreat_Y<-mat.or.vec(calib_data_Ultuna$J_Ultuna,dim[2])
Ultuna_min_predictions_bytreat_Y<-mat.or.vec(calib_data_Ultuna$J_Ultuna,dim[2])
for(i in 1:calib_data_Ultuna$J_Ultuna){
  Ultuna_mean_predictions_bytreat_Y[i,]<-Ultuna_mean_predictions[i,]
  Ultuna_max_predictions_bytreat_Y[i,]<-Ultuna_min_predictions[i,]
  Ultuna_min_predictions_bytreat_Y[i,]<-Ultuna_max_predictions[i,]
}


#png("ICBM_predictions_Ultuna_Y_TEST_specific.png", height=4500, width=4000, res=320)
i=1
yearseq<-seq(from=colnames(calib_data_Ultuna$SOC_Ultuna)[1], to=as.numeric(colnames(calib_data_Ultuna$SOC_Ultuna)[1])+dim[2]-1)
  lastyear<-dim(Ultuna_mean_predictions_bytreat_Y)[2]-1
  plot(yearseq, Ultuna_mean_predictions_bytreat_Y[i,], ylim=c(0,1), type="l",  ylab=expression(paste("C (t ha" ^-1, ")")), xlab="Years", main=paste("Treatment", LETTERS[selection_treatments[i]]), col=Ultuna_treat_palette[i])
  polygon(c(yearseq,rev(yearseq)),
          c(Ultuna_min_predictions_bytreat_Y[i,],rev(Ultuna_max_predictions_bytreat_Y[i,])),
          col=Ultuna_treat_palette_alpha[i],border=Ultuna_treat_palette[i])
  lines(yearseq, Ultuna_mean_predictions_bytreat_Y[i,], col="black", lty=1, lwd=1)












###Plot parameter posteriors

# parameters_palette_multi<-c("darkorange","cadetblue", 
#                             "cornsilk4","darkgreen","darkorchid",
#                             "lightpink4","goldenrod3","burlywood4","blueviolet","darkcyan",
#                             brewer.pal(6, "BuPu")[2:6], "red")

colors = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
set.seed(25)
parameters_palette_multi=sample(colors, 15)
pie(rep(1,15), col=parameters_palette_multi)

parameters_palette_multi_alpha<-add.alpha(parameters_palette_multi,0.7)


parameter_names<-c(expression(paste("k"[1])),
                   expression(paste("k"[2])),
                   expression(paste("h"[R])),
                   expression(paste("h"[S])),
                   expression(paste("h"[FYM])),
                   expression(paste("h"[PEAT])),
                   expression(paste("h"[SAW])),
                   expression(paste("h"[SLU])),
                   "Stubbles ratio",
                   "Stubbles ratio (maize)",
                   "S:R cereals","S:R root crops","S:R oilseeds", "S:R maize",
                   "O/Y ratio")



#priors
library(truncnorm)
##Parameters

#error_h<-0.1
limits_h<-1
limits_k<-1
limit_SR<-1

#k1  <- runif(60000,0.8-0.8*limits_k, 0.8+0.8*limits_k)
#k2  <- runif(60000, min=0.00605-0.00605*limits_k, max=0.00605+0.00605*limits_k)
# k1  <- runif(60000, min=0, max=1)
# k2  <- runif(60000, min=0, max=0.03)

k1  <- rnorm(60000, mean=0.65, sd=0.65*0.1)
k2  <- rnorm(60000, mean=0.0085, sd=0.0085*0.1)



h_S     <- runif(60000, min=0.125-0.125*limits_h, max=0.125+0.125*limits_h)
h_R     <- runif(60000, min=0.35-0.35*limits_h, max=0.35+0.35*limits_h)
h_FYM   <- runif(60000,  min=0.27-0.27*limits_h, max=0.27+0.27*limits_h)
h_PEA   <- runif(60000,  min=0.59-0.59*limits_h, max=0.59+0.59*limits_h)
h_SAW   <- runif(60000,  min=0.25-0.25*limits_h, max=0.25+0.25*limits_h)
h_SLU   <- runif(60000,  min=0.41-0.41*limits_h, max=0.41+0.41*limits_h)

stubbles_ratio_Ultuna <- runif(60000, 0.01,0.08)
stubbles_ratio_Ultuna_maize <-  runif(60000, 0.01,0.08)
# Init_ratio_Ultuna <- rtruncnorm(60000, mean=0.9291667, sd=0.0001, a=0.9, b=0.95)
# Init_ratio_Lanna <- rtruncnorm(60000, mean=0.9291667, sd=0.0001, a=0.9, b=0.95)
exudates_coeff <- runif(60000, min=1.65-(1.65*0.1), max=1.65+(1.65*0.1))
plot(density(exudates_coeff))


#root/shoot ratios priors
#error_SR<-0.25

SR_cereals    <- runif(60000, 3.6,27.9)
SR_root_crops <- runif(60000, 29.49853-29.49853*limit_SR,29.49853+29.49853*limit_SR)
SR_oilseeds   <- runif(60000, 8-8*limit_SR,8+8*limit_SR)
SR_maize   <- runif(60000, 6.25-6.25*limit_SR,30)

Init_ratio_Ultuna   <- runif(60000,  min=0.8, max=0.98)

C_percent <- runif(60000,min=0.40, max=0.51)
alpha  <- runif(60000,  min=0, max=5) #specification of the prior
inert  <- runif(60000,  min=0, max=15) #specification of the prior


##Parameters
confidence=c(0.025,0.975)
parameter_values<-mat.or.vec(length(parameter_list[1:28]),4)
colnames(parameter_values)<-c("mode","mean","min","max")
rownames(parameter_values)<-parameter_list[1:28]

priors_list<-list(k1, #                           1
                  k2, #                           2
                 #k2_org, #                       3
                  h_R, #                          4
                  h_S, #                          5
                  h_FYM, #                        6
                  h_PEA, #                        7
                  h_SAW, #                        8
                  h_SLU, #                        9
                  #exudates_coeff, #              10
                  stubbles_ratio_Ultuna, #        10
                  stubbles_ratio_Ultuna_maize, #  11
                  SR_cereals, #                   12
                  SR_root_crops, #                13
                  SR_oilseeds, #                  14
                  SR_maize,#                      15
                 Init_ratio_Ultuna) #             16

priors_list_d<-list()
for(i in 1:length(priors_list)){priors_list_d[[i]]<-density(priors_list[[i]])}
chain_list_multisite<-list()
MCMC_list_multisite<-list()
density_list_multisite<-list()
params_mode_multisite<-c()
params_meam_multisite<-c()
for(i in 1:(length(parameter_list[1:28]))){
  chain_list_multisite[[i]]<-as.mcmc.list(eval(parse(text=paste("mcmc.array.ICBM.Ultuna","$",parameter_list[i], sep=""))), chains=T)
  density_list_multisite[[i]]<-density(unlist(chain_list_multisite[[i]]), adjust=2)
  if(length(unique(unlist(chain_list_multisite[[i]])))>1){
    params_mode_multisite[i]<-mlv(unlist(chain_list_multisite[[i]]), method = "Parzen")
  }else{
    params_mode_multisite[i]<-mlv(unlist(chain_list_multisite[[i]]))
    
  }
  params_meam_multisite[i]<-mean(unlist(chain_list_multisite[[i]]))
  parameter_values[i,1]<-round(params_mode_multisite[i],6)
  parameter_values[i,2]<-round(params_meam_multisite[i],6)
  parameter_values[i,3]<-round(quantile(as.matrix(chain_list_multisite[[i]][1]), confidence)[1],6)
  parameter_values[i,4]<-round(quantile(as.matrix(chain_list_multisite[[i]][1]), confidence)[2],6)
  #plot(chain_list_multisite[[i]])
}



parameter_values_mat<-cbind(rownames(parameter_values), round(parameter_values,3))
colnames(parameter_values_mat)[1]<-"parameter"
parameter_values_mat[,1]<-gsub("_", "", parameter_values_mat[,1])
rownames(parameter_values_mat) <- c()


palette_Ult_Lan<-add.alpha(c("Firebrick", "Darkgreen"),0.5)

write.csv(parameter_values_mat, file="parameter_values.csv")

#traceplots to control
for (i in 1:14) {
  jpeg(paste("./traceplots/ICBM_posteriors_traceplots", i, '.jpeg', sep = ''), height=1500, width=3600, res=300)
  plot(chain_list_multisite[[i]])
  dev.off()
}

i=23
jpeg(paste("./traceplots/ICBM_posteriors_traceplots", i, '.jpeg', sep = ''), height=1500, width=3600, res=300)
plot(chain_list_multisite[[i]])
dev.off()

i=26
jpeg(paste("./traceplots/ICBM_posteriors_traceplots", i, '.jpeg', sep = ''), height=1500, width=3600, res=300)
plot(chain_list_multisite[[i]])
dev.off()

i=27
jpeg(paste("./traceplots/ICBM_posteriors_traceplots", i, '.jpeg', sep = ''), height=1500, width=3600, res=300)
plot(chain_list_multisite[[i]])
dev.off()

i=28
jpeg(paste("./traceplots/ICBM_posteriors_traceplots", i, '.jpeg', sep = ''), height=1500, width=3600, res=300)
plot(chain_list_multisite[[i]])
dev.off()

max(unlist(chain_list_multisite[[i]]))

i=27
jpeg(paste("./traceplots/ICBM_posteriors_traceplots", i, '.jpeg', sep = ''), height=1500, width=3600, res=300)
plot(chain_list_multisite[[i]])
dev.off()




jpeg("ICBM_posteriors_multisite_specific.jpg", height=2600, width=4000, res=300)
par(mfrow=c(3,5))

i=1
range<-range(c(density_list_multisite[[i]]$x, priors_list_d[[i]]$x, priors_list_d[[i]]$x))
  mean<-mean(c(density_list_multisite[[i]]$x, priors_list_d[[i]]$x))
  rangey<-range(c(density_list_multisite[[i]]$y, priors_list_d[[i]]$y, priors_list_d[[i]]$y))
  meany<-mean(c(density_list_multisite[[i]]$y, priors_list_d[[i]]$y))

  plot(density_list_multisite[[i]], main=parameter_names[i],xlim=c(range[1]-(range[1]*0.015), range[2]+(range[2]*0.015)), ylim=c(0, rangey[2]+(rangey[2]*0.2)))
  polygon(density_list_multisite[[i]], col=add.alpha(parameters_palette_multi[i],0.8))
  polygon(priors_list_d[[i]], col=add.alpha("lightgrey", 0.6), lty=2)
  if(i==1){legend("topright", c("Ultuna", "Prior"), bty="n", pch=c(16),  col=c(add.alpha(parameters_palette_multi[i],0.8), add.alpha("lightgrey", 0.6)))}
  abline(v=params_mode_multisite[i], col="red", lty=3)
  legend("topleft", paste("value=",round(params_mode_multisite[i],2)), bty="n")
  

i=2
range<-range(c(density_list_multisite[[i]]$x, priors_list_d[[i]]$x, priors_list_d[[i]]$x))
mean<-mean(c(density_list_multisite[[i]]$x, priors_list_d[[i]]$x))
rangey<-range(c(density_list_multisite[[i]]$y, priors_list_d[[i]]$y, priors_list_d[[i]]$y))
meany<-mean(c(density_list_multisite[[i]]$y, priors_list_d[[i]]$y))

plot(density_list_multisite[[i]], main=parameter_names[i],xlim=c(range[1]-(range[1]*0.015), range[2]+(range[2]*0.015)), ylim=c(0, rangey[2]+(rangey[2]*0.2)))
polygon(density_list_multisite[[i]], col=add.alpha(parameters_palette_multi[i],0.8))
polygon(priors_list_d[[i]], col=add.alpha("lightgrey", 0.6), lty=2)
if(i==1){legend("topright", c("Ultuna", "Prior"), bty="n", pch=c(16),  col=c(add.alpha(parameters_palette_multi[i],0.8), add.alpha("lightgrey", 0.6)))}
abline(v=params_mode_multisite[i], col="red", lty=3)
legend("topleft", paste("value=",round(params_mode_multisite[i],3)), bty="n")


#H_r
i=3
parameter_list[c(i, i+6)]
range<-range(c(density_list_multisite[[i]]$x, priors_list_d[[i]]$x, priors_list_d[[i]]$x))
mean<-mean(c(density_list_multisite[[i]]$x, priors_list_d[[i]]$x))
rangey<-range(c(density_list_multisite[[i]]$y, priors_list_d[[i]]$y, priors_list_d[[i]]$y))
meany<-mean(c(density_list_multisite[[i]]$y, priors_list_d[[i]]$y))
plot(density_list_multisite[[i]], main=parameter_names[i],xlim=c(range[1]-(range[1]*0.015), range[2]+(range[2]*0.015)), ylim=c(0, rangey[2]+(rangey[2]*0.2)), col=NA)
polygon(density_list_multisite[[i]], col=add.alpha(parameters_palette_multi[i],0.8))
polygon(priors_list_d[[i]], col=add.alpha("lightgrey", 0.6), lty=2)
abline(v=params_mode_multisite[i], col="red", lty=3)
legend("topleft", paste("value=",round(params_mode_multisite[i],2)), bty="n")

#H_s
i=4
parameter_list[c(i, i+6)]
range<-range(c(density_list_multisite[[i]]$x, priors_list_d[[i]]$x, priors_list_d[[i]]$x))
mean<-mean(c(density_list_multisite[[i]]$x, priors_list_d[[i]]$x))
rangey<-range(c(density_list_multisite[[i]]$y, priors_list_d[[i]]$y,priors_list_d[[i]]$y))
meany<-mean(c(density_list_multisite[[i]]$y, priors_list_d[[i]]$y))
plot(density_list_multisite[[i]], main=parameter_names[i],xlim=c(range[1]-(range[1]*0.015), range[2]+(range[2]*0.015)), ylim=c(0, rangey[2]+(rangey[2]*0.2)), col=NA)
polygon(density_list_multisite[[i]], col=add.alpha(parameters_palette_multi[i],0.8))
polygon(priors_list_d[[i]], col=add.alpha("lightgrey", 0.6), lty=2)
abline(v=params_mode_multisite[i], col="red", lty=3)
legend("topleft", paste("value=",round(params_mode_multisite[i],2)), bty="n")


#H_FYM
i=5
parameter_list[i]
range<-range(c(density_list_multisite[[i]]$x, priors_list_d[[i]]$x,priors_list_d[[i]]$x))
mean<-mean(c(density_list_multisite[[i]]$x, priors_list_d[[i]]$x))
rangey<-range(c(density_list_multisite[[i]]$y, priors_list_d[[i]]$y, priors_list_d[[i]]$y))
meany<-mean(c(density_list_multisite[[i]]$y, priors_list_d[[i]]$y))
plot(density_list_multisite[[i]], main=parameter_names[i],xlim=c(range[1]-(range[1]*0.015), range[2]+(range[2]*0.015)), ylim=c(0, rangey[2]+(rangey[2]*0.2)), col=NA)
polygon(density_list_multisite[[i]], col=add.alpha(parameters_palette_multi[i],0.8))
polygon(priors_list_d[[i]], col=add.alpha("lightgrey", 0.6), lty=2)
abline(v=params_mode_multisite[i], col="red", lty=3)
legend("topleft", paste("value=",round(params_mode_multisite[i],2)), bty="n")


#H_PEA
i=6
parameter_list[c(i)]
range<-range(c(density_list_multisite[[i]]$x, priors_list_d[[i]]$x))
mean<-mean(c(density_list_multisite[[i]]$x, priors_list_d[[i]]$x))
rangey<-range(c(density_list_multisite[[i]]$y, priors_list_d[[i]]$y))
meany<-mean(c(density_list_multisite[[i]]$y, priors_list_d[[i]]$y))
plot(density_list_multisite[[i]], main=parameter_names[i],xlim=c(range[1]-(range[1]*0.015), range[2]+(range[2]*0.015)), ylim=c(0, rangey[2]+(rangey[2]*0.2)), col=NA)
polygon(density_list_multisite[[i]], col=add.alpha(parameters_palette_multi[i],0.8)) #col=add.alpha("lightgrey", 0.6))
polygon(priors_list_d[[i]], col=add.alpha("lightgrey", 0.6), lty=2)
abline(v=params_mode_multisite[i], col="red", lty=3)
legend("topleft", paste("value=",round(params_mode_multisite[i],2)), bty="n")


#H_SAW
i=7
parameter_list[c(i, i+6)]
range<-range(c(density_list_multisite[[i]]$x, priors_list_d[[i]]$x))
mean<-mean(c(density_list_multisite[[i]]$x, priors_list_d[[i]]$x))
rangey<-range(c(density_list_multisite[[i]]$y, priors_list_d[[i]]$y))
meany<-mean(c(density_list_multisite[[i]]$y, priors_list_d[[i]]$y))
plot(density_list_multisite[[i]], main=parameter_names[i],xlim=c(range[1]-(range[1]*0.015), range[2]+(range[2]*0.015)), ylim=c(0, rangey[2]+(rangey[2]*0.2)), col=NA)
polygon(density_list_multisite[[i]], col=add.alpha(parameters_palette_multi[i],0.8))
polygon(priors_list_d[[i]], col=add.alpha("lightgrey", 0.6), lty=2)
abline(v=params_mode_multisite[i], col="red", lty=3)
legend("topleft", paste("value=",round(params_mode_multisite[i],2)), bty="n")

#H_SLU
i=8
parameter_list[c(i, i+4)]
range<-range(c(density_list_multisite[[i]]$x, priors_list_d[[i]]$x,priors_list_d[[i]]$x))
mean<-mean(c(density_list_multisite[[i]]$x, priors_list_d[[i]]$x))
rangey<-range(c(density_list_multisite[[i]]$y, priors_list_d[[i]]$y,priors_list_d[[i]]$y))
meany<-mean(c(density_list_multisite[[i]]$y, priors_list_d[[i]]$y))
plot(density_list_multisite[[i]], main=parameter_names[i],xlim=c(range[1]-(range[1]*0.015), range[2]+(range[2]*0.015)), ylim=c(0, rangey[2]+(rangey[2]*0.2)), col=NA)
polygon(density_list_multisite[[i]], col=add.alpha(parameters_palette_multi[i],0.8))
polygon(priors_list_d[[i]], col=add.alpha("lightgrey", 0.6), lty=2)
abline(v=params_mode_multisite[i], col="red", lty=3)
legend("topleft", paste("value=",round(params_mode_multisite[i],2)), bty="n")


#stubbles ratio
i=9
range<-range(c(density_list_multisite[[9]]$x, priors_list_d[[9]]$x))
mean<-mean(c(density_list_multisite[[i]]$x, priors_list_d[[9]]$x))
rangey<-range(c(density_list_multisite[[9]]$y, priors_list_d[[9]]$y))
plot(density_list_multisite[[9]], main=parameter_names[9],xlim=c(range[1]-(range[1]*0.015), range[2]+(range[2]*0.015)), ylim=c(0, rangey[2]+(rangey[2]*0.2)))
polygon(density_list_multisite[[9]], col=add.alpha(parameters_palette_multi[i],0.8))
polygon(priors_list_d[[9]], col=add.alpha("lightgrey", 0.6), lty=2)
abline(v=params_mode_multisite[i], col="red", lty=3)
legend("topleft", paste("value=",round(params_mode_multisite[i],2)), bty="n")


#stubbles ratio maize
i=10
range<-range(c(density_list_multisite[[10]]$x, priors_list_d[[10]]$x))
mean<-mean(c(density_list_multisite[[i]]$x, priors_list_d[[10]]$x))
rangey<-range(c(density_list_multisite[[10]]$y, priors_list_d[[10]]$y))
plot(density_list_multisite[[10]], main=parameter_names[10],xlim=c(range[1]-(range[1]*0.015), range[2]+(range[2]*0.015)), ylim=c(0, rangey[2]+(rangey[2]*0.2)))
polygon(density_list_multisite[[10]], col=add.alpha(parameters_palette_multi[i],0.8))
polygon(priors_list_d[[10]], col=add.alpha("lightgrey", 0.6), lty=2)
abline(v=params_mode_multisite[i], col="red", lty=3)
legend("topleft", paste("value=",round(params_mode_multisite[i],2)), bty="n")


for(i in 11:14){
  range<-range(c(density_list_multisite[[i]]$x, priors_list_d[[i]]$x))
  mean<-mean(c(density_list_multisite[[i]]$x, priors_list_d[[i]]$x))
  rangey<-range(c(density_list_multisite[[i]]$y, priors_list_d[[i]]$y))
  meany<-mean(c(density_list_multisite[[i]]$y, priors_list_d[[i]]$y))
  plot(density_list_multisite[[i]], main=parameter_names[i],xlim=c(range[1]-(range[1]*0.015), range[2]+(range[2]*0.015)), ylim=c(0, rangey[2]+(rangey[2]*0.2)))
  polygon(density_list_multisite[[i]], col=add.alpha(parameters_palette_multi[i],0.8))
  polygon(priors_list_d[[i]], col=add.alpha("lightgrey", 0.6), lty=2)
  abline(v=params_mode_multisite[i], col="red", lty=3)
  legend("topleft", paste("value=",round(params_mode_multisite[i],1)), bty="n")
}


#init
i=23
range<-range(c(density_list_multisite[[i]]$x, priors_list_d[[15]]$x))
mean<-mean(c(density_list_multisite[[i]]$x, priors_list_d[[15]]$x))
rangey<-range(c(density_list_multisite[[i]]$y, priors_list_d[[15]]$y))
plot(density_list_multisite[[i]], main=parameter_names[15],xlim=c(range[1]-(range[1]*0.015), range[2]+(range[2]*0.015)), ylim=c(0, rangey[2]+(rangey[2]*0.2)))
polygon(density_list_multisite[[i]], col=add.alpha(parameters_palette_multi[15],0.8))
polygon(priors_list_d[[15]], col=add.alpha("lightgrey", 0.6), lty=2)
abline(v=params_mode_multisite[i], col="red", lty=3)
legend("topleft", paste("value=",round(params_mode_multisite[i],2)), bty="n")


dev.off()




jpeg("ICBM_posteriors_multisite_specific_allometry.jpg", height=1200, width=1200, res=300)
i=24
range<-range(c(density_list_multisite[[i]]$x, density(alpha)$x))
mean<-mean(c(density_list_multisite[[i]]$x,  density(alpha)$x))
rangey<-range(c(density_list_multisite[[i]]$y,  density(alpha)$y))
plot(density_list_multisite[[i]], main=expression(alpha),xlim=c(range[1]-(range[1]*0.015), range[2]+(range[2]*0.015)), ylim=c(0, rangey[2]+(rangey[2]*0.2)), xlab=expression("Allometric intercept ( t ha"^-1~")"))
polygon(density_list_multisite[[i]], col=add.alpha("darkorange",0.6))
polygon(density_list_multisite[[1+i]], col=add.alpha("chartreuse2",0.6))
polygon(density(alpha), col=add.alpha("lightgrey", 0.6), lty=2)
legend("topright", c("maize", "other crops"), bty="n", pch=16, col=add.alpha(c("chartreuse2","darkorange"), 0.6))
abline(v=params_mode_multisite[i], col="red", lty=3)
abline(v=mlv(density_list_multisite[[i+1]]$x, method="Parzen"), col="red", lty=3)
legend("topleft", c(paste("value others=",round(params_mode_multisite[i],2)),
                    paste("value maize=",round(mlv(density_list_multisite[[i+1]]$x, method="Parzen"),2))), bty="n")

dev.off()





jpeg("ICBM_posteriors_multisite_specific_Cprocent.jpg", height=1200, width=1200, res=300)
i=26
range<-range(c(density_list_multisite[[i]]$x, density(C_percent)$x))
mean<-mean(c(density_list_multisite[[i]]$x,  density(C_percent)$x))
rangey<-range(c(density_list_multisite[[i]]$y,  density(C_percent)$y))
plot(density_list_multisite[[i]], main="C% content of SOC",xlim=c(range[1]-(range[1]*0.015), range[2]+(range[2]*0.015)), ylim=c(0, rangey[2]+(rangey[2]*0.2)))
polygon(density_list_multisite[[i]], col=add.alpha("lightblue",0.6))
polygon(density(C_percent), col=add.alpha("lightgrey", 0.6), lty=2)
abline(v=params_mode_multisite[i], col="red", lty=3)
legend("topleft", paste("value=",round(params_mode_multisite[i],2)), bty="n")
dev.off()




jpeg("ICBM_posteriors_multisite_specific_Inert.jpg", height=1200, width=1200, res=300)
i=27
range<-range(c(density_list_multisite[[i]]$x, density(inert)$x))
mean<-mean(c(density_list_multisite[[i]]$x,  density(inert)$x))
rangey<-range(c(density_list_multisite[[i]]$y,  density(inert)$y))
plot(density_list_multisite[[i]], main="Inert pool",xlim=c(range[1]-(range[1]*0.015), range[2]+(range[2]*0.015)), ylim=c(0, rangey[2]+(rangey[2]*0.2)))
polygon(density_list_multisite[[i]], col=add.alpha("lightblue",0.6))
polygon(density(inert), col=add.alpha("lightgrey", 0.6), lty=2)
abline(v=params_mode_multisite[i], col="red", lty=3)
legend("topleft", paste("value=",round(params_mode_multisite[i],2)), bty="n")
dev.off()


jpeg("ICBM_posteriors_multisite_specific_exudates.jpg", height=1200, width=1200, res=300)
i=28
range<-range(c(density_list_multisite[[i]]$x, density(exudates_coeff)$x))
mean<-mean(c(density_list_multisite[[i]]$x,  density(exudates_coeff)$x))
rangey<-range(c(density_list_multisite[[i]]$y,  density(exudates_coeff)$y))
plot(density_list_multisite[[i]], main="Exudate coefficient",xlim=c(range[1]-(range[1]*0.015), range[2]+(range[2]*0.015)), ylim=c(0, rangey[2]+(rangey[2]*0.2)))
polygon(density_list_multisite[[i]], col=add.alpha("lightblue",0.6))
polygon(density(exudates_coeff), col=add.alpha("lightgrey", 0.6), lty=2)
abline(v=params_mode_multisite[i], col="red", lty=3)
legend("topleft", paste("value=",round(params_mode_multisite[i],2)), bty="n")
dev.off()












################correlations


library(MASS)

z <- kde2d(x=unlist(chain_list_multisite[[3]]), y=unlist(chain_list_multisite[[11]]), n = 50)

png("hr_vs_RS_dottyplot.png", width = 2000, height=2000, res=300)
plot(x=unlist(chain_list_multisite[[3]]), y=unlist(chain_list_multisite[[11]]), pch = 19,xlab=parameter_names[3],ylab=parameter_names[11])
dev.off()

png("hr_vs_RS.png", width = 3000, height=3000, res=300)
filled.contour(z,xlab=parameter_names[3],ylab=parameter_names[11])
dev.off()

library(plotly)
x=unlist(chain_list_multisite[[3]])
y=unlist(chain_list_multisite[[11]])
cor(x,y)
freqz <- with(data.frame(x,y), MASS::kde2d(x, y, n = 50))
with(freqz, plot_ly(x = x, y = y, z = z, type = "surface")%>%
       layout(scene = list(xaxis = list(title = "h_r"), yaxis = list(title = "S:R cereals")))) 

# 
 library(ggplot2)
 df<-data.frame(x,y)
 png("hr_vs_RS_ggplot.png", width = 2000, height=2000, res=300)
 ggplot(df, aes(x = x, y = y)) +
   geom_point() +
   geom_density_2d_filled(alpha = 0.8, show.legend = FALSE, contour_var= "ndensity", bins=15) +
   geom_density_2d(colour = "black", contour_var= "ndensity", bins=15)+
   labs(x = parameter_names[3], y = parameter_names[11])+
   theme_classic()

 dev.off()

 
 x=unlist(chain_list_multisite[[2]])
 y=unlist(chain_list_multisite[[11]])/unlist(chain_list_multisite[[3]])
 freqz <- with(data.frame(x,y), MASS::kde2d(x, y, n = 50))
 with(freqz, plot_ly(x = x, y = y, z = z, type = "surface")%>%
        layout(scene = list(xaxis = list(title = "h_r"), yaxis = list(title = "hr/S:R cereals")))) 
 
 library(ggplot2)
 df<-data.frame(x,y)
 png("k2_vs_RS-hr_ratio_ggplot.png", width = 2000, height=2000, res=300)
 ggplot(df, aes(x = x, y = y)) +
   geom_point() +
   geom_density_2d_filled(alpha = 0.8, show.legend = FALSE, contour_var= "ndensity", bins=15) +
   geom_density_2d(colour = "black", contour_var= "ndensity", bins=15)+
   labs(x = parameter_names[2], y = expression(paste("h"[R],"/ S:R cereals")))+
   theme_classic()
 
 dev.off()
 
 
 
 
 

 #residuals
dim(calib_data_Ultuna$SOC_Ultuna)
dim(Ultuna_mean_predictions_bytreat)
residuals<-calib_data_Ultuna$SOC_Ultuna-Ultuna_mean_predictions_bytreat
num_residuals<-length(na.omit(as.vector(unlist(residuals))))
RMSE<-sqrt(((sum(residuals^2, na.rm=T)))/num_residuals)


residuals_type<-calib_data_Ultuna$Yields_cereals_Ultuna
residuals_type[!calib_data_Ultuna$Yields_cereals_Ultuna==0]="Cereals"
residuals_type[!calib_data_Ultuna$Yields_root_crops_Ultuna==0]="Roots"
residuals_type[!calib_data_Ultuna$Yields_oilseeds_Ultuna==0]="Oilseeds"
residuals_type[!calib_data_Ultuna$Yields_maize_Ultuna==0]="Maize"


png("residuals.png", width = 5000, height=4000, res=300)
par(mfrow=c(5,3))
i=1
plot(as.numeric(residuals[i,]), axes=FALSE, ylab="mean residuals", xlab="", main=paste("Treatment", LETTERS[selection_treatments[i]]),  ylim=range(residuals, na.rm=T), col=Ultuna_crop_palette[as.numeric(as.factor(residuals_type[1,]))],  pch=as.numeric(as.factor(residuals_type[1,])))
axis(1, at=seq(1, length(residuals_type)), labels= colnames(residuals[i,]), las=2, cex.axis=0.7)
axis(2) #default way
abline(h=0, lty=2)
legend("topleft", levels(as.factor(unlist(residuals_type))), pch=seq(1,5), col = Ultuna_crop_palette, bty="n")
box()
for(i in 2:calib_data_Ultuna$J_Ultuna){
plot(as.numeric(residuals[i,]), axes=FALSE, ylab="mean residuals", xlab="", main=paste("Treatment", LETTERS[selection_treatments[i]]),  ylim=range(residuals, na.rm=T), col=Ultuna_crop_palette[as.numeric(as.factor(residuals_type[2,]))],  pch=as.numeric(as.factor(residuals_type[2,])))
axis(1, at=seq(1, length(residuals_type)), labels= colnames(residuals_type[i,]), las=2, cex.axis=0.7)
axis(2) #default way
abline(h=0, lty=2)
legend("topleft", levels(as.factor(unlist(residuals_type))), pch=seq(1,5), col = Ultuna_crop_palette, bty="n")
box()
}
dev.off()


png("residuals_freescale.png", width = 5000, height=4000, res=300)
par(mfrow=c(5,3))
for(i in 1:calib_data_Ultuna$J_Ultuna){
  plot(as.numeric(residuals[i,]), axes=FALSE, ylab="mean residuals", xlab="", main=paste("Treatment", LETTERS[selection_treatments[i]]),  pch=as.numeric(as.factor(residuals_type[1,])), col=Ultuna_crop_palette[as.numeric(as.factor(residuals_type[1,]))],  )
  axis(1, at=seq(1, length(residuals_type)), labels= colnames(residuals[i,]), las=2, cex.axis=0.5)
  legend("topleft", levels(as.factor(unlist(residuals_type))), pch=seq(1,5), col = Ultuna_crop_palette, bty="n")
  axis(2) #default way
  abline(h=0, lty=2)
  box()
}
dev.off()

png("residuals_mean.png", width = 2000, height=2000, res=300)
barplot(rowMeans(residuals, na.rm=T)[1:calib_data_Ultuna$J_Ultuna])
dev.off()


#function to find the p-value of a regression
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}



png("residuals_scatterplots.png", width = 4000, height=2000, res=300)
par(mfrow=c(1,2))
yields<-(calib_data_Ultuna$Yields_cereals_Ultuna[1:calib_data_Ultuna$J,]+calib_data_Ultuna$Yields_root_crops_Ultuna[1:calib_data_Ultuna$J,]+calib_data_Ultuna$Yields_oilseeds_Ultuna[1:calib_data_Ultuna$J,]+calib_data_Ultuna$Yields_maize_Ultuna[1:calib_data_Ultuna$J,])
plot(unlist(yields),unlist(residuals[1:calib_data_Ultuna$J,]), main="Residuals vs yields", 
     col=Ultuna_crop_palette[as.numeric(as.factor(unlist(residuals_type[1:calib_data_Ultuna$J,])))], pch=as.numeric(as.factor(unlist(residuals_type[1:calib_data_Ultuna$J,]))), xlab="Yields", ylab="Residuals")
text(unlist(yields),unlist(residuals[1:calib_data_Ultuna$J,])+max(unlist(residuals[1:calib_data_Ultuna$J,]), na.rm = T)*0.03, rep(rownames(yields), dim(yields)[2]), col=Ultuna_crop_palette[as.numeric(as.factor(unlist(residuals_type[1:calib_data_Ultuna$J,])))], cex=0.7)
legend("bottomright", levels(as.factor(unlist(residuals_type[1:calib_data_Ultuna$J,]))), pch=seq(1,5), col = Ultuna_crop_palette, bty="n")
abline(h=0, lty=2)

regr_1st_term<-(unlist(yields))
regr_2nd_term<-(unlist(residuals[1:calib_data_Ultuna$J,]))
regr<-lm(as.vector(regr_2nd_term)~regr_1st_term )
#abline(regr, lty=1, col="darkgrey")
#legend("topleft", c(paste("R2=",round(summary(regr)$r.squared,2)),
#                    paste("p=",round(lmp(regr),10))), bty="n")


plot(unlist(calib_data_Ultuna$re_Ultuna[1:calib_data_Ultuna$J,]),unlist(residuals[1:calib_data_Ultuna$J,]), main="Residuals vs reclim", col=Ultuna_crop_palette[as.numeric(as.factor(unlist(residuals_type[1:calib_data_Ultuna$J,])))], pch=as.numeric(as.factor(unlist(residuals_type[1:calib_data_Ultuna$J,]))), xlab="reclim", ylab="Residuals")
text(unlist(calib_data_Ultuna$re_Ultuna[1:calib_data_Ultuna$J,]),unlist(residuals[1:calib_data_Ultuna$J,])+max(unlist(residuals[1:calib_data_Ultuna$J,]), na.rm = T)*0.03, rep(rownames(yields), dim(yields)[2]), col=Ultuna_crop_palette[as.numeric(as.factor(unlist(residuals_type[1:calib_data_Ultuna$J,])))], cex=0.7)
abline(h=0, lty=2)

dev.off()


source("ICBM_single.R")
posteriors
rownames(parameter_values)

mean_mode=1 #1 = mode, 2= mean

A_treat_sim<-ICBM_single(
  BF=T,
  k1=parameter_values[1,mean_mode],
  k2=parameter_values[2,mean_mode],
  h_r=parameter_values[3,mean_mode],
  h_s=parameter_values[4,mean_mode],
  h_fym=parameter_values[5,mean_mode],
  h_pea=parameter_values[6,mean_mode],
  h_saw=parameter_values[7,mean_mode],
  h_slu=parameter_values[8,mean_mode],
  stubbles_others=parameter_values[9,mean_mode],
  stubbles_maize=parameter_values[10,mean_mode],
  S_R_cereals=parameter_values[11,mean_mode],
  S_R_roots=parameter_values[12,mean_mode],
  S_R_oilseeds=parameter_values[13,mean_mode],
  S_R_maize=parameter_values[14,mean_mode],
  OY_ratio=parameter_values[23,mean_mode],
  alpha_others=parameter_values[24,mean_mode],
  alpha_maize=parameter_values[25,mean_mode],
  C_percent=parameter_values[26,mean_mode], 
  Yields_cereals=as.numeric(calib_data_Ultuna$Yields_cereals_Ultuna[1,]),
  Yields_root_crops=as.numeric(calib_data_Ultuna$Yields_root_crops_Ultuna[1,]),
  Yields_oilseeds=as.numeric(calib_data_Ultuna$Yields_oilseeds_Ultuna[1,]),
  Yields_maize=as.numeric(calib_data_Ultuna$Yields_maize_Ultuna[1,]),
  I_fym=as.numeric(calib_data_Ultuna$I_FYM_Ultuna[1,]),
  I_gm=as.numeric(calib_data_Ultuna$I_GM_Ultuna[1,]),
  I_pea=as.numeric(calib_data_Ultuna$I_PEA_Ultuna[1,]),
  I_saw=as.numeric(calib_data_Ultuna$I_SAW_Ultuna[1,]),
  I_slu=as.numeric(calib_data_Ultuna$I_SLU_Ultuna[1,]),
  I_str=as.numeric(calib_data_Ultuna$I_STR_Ultuna[1,]),
  SOC_init=as.numeric(calib_data_Ultuna$SOC_Ultuna[1,1]),
  re=as.numeric(calib_data_Ultuna$re_Ultuna[1,]))

plot(A_treat_sim, type="l", xlab="Year", ylab="SOC t ha-1")
points(calib_data_Ultuna$SOC_Ultuna[1,])


write.csv(A_treat_sim, "tmp.csv")


Ultuna_mean_Ir_bytreat_DF<-as.data.frame(Ultuna_mean_Ir_bytreat)
Ultuna_min_Ir_bytreat_DF<-as.data.frame(Ultuna_min_Ir_bytreat)
Ultuna_max_Ir_bytreat_DF<-as.data.frame(Ultuna_max_Ir_bytreat)
Ultuna_mean_Is_bytreat_DF<-as.data.frame(Ultuna_mean_Is_bytreat)
Ultuna_min_Is_bytreat_DF<-as.data.frame(Ultuna_min_Is_bytreat)
Ultuna_max_Is_bytreat_DF<-as.data.frame(Ultuna_max_Is_bytreat)


colnames(Ultuna_mean_Ir_bytreat_DF)<-yearseq
rownames(Ultuna_mean_Ir_bytreat_DF)<-rownames(Ultuna_mean_Ir)
colnames(Ultuna_min_Ir_bytreat_DF)<-yearseq
rownames(Ultuna_min_Ir_bytreat_DF)<-rownames(Ultuna_min_Ir)
colnames(Ultuna_max_Ir_bytreat_DF)<-yearseq
rownames(Ultuna_max_Ir_bytreat_DF)<-rownames(Ultuna_max_Ir)

colnames(Ultuna_mean_Is_bytreat_DF)<-yearseq
rownames(Ultuna_mean_Is_bytreat_DF)<-rownames(Ultuna_mean_Is)
colnames(Ultuna_min_Is_bytreat_DF)<-yearseq
rownames(Ultuna_min_Is_bytreat_DF)<-rownames(Ultuna_min_Is)
colnames(Ultuna_max_Is_bytreat_DF)<-yearseq
rownames(Ultuna_max_Is_bytreat_DF)<-rownames(Ultuna_max_Is)


write_ods(Ultuna_mean_Ir_bytreat_DF, "inputs_estiamtes_for_Carlos.ods", sheet = "mean_Ir")
write_ods(Ultuna_min_Ir_bytreat_DF, "inputs_estiamtes_for_Carlos.ods", sheet = "min_Ir", append=T)
write_ods(Ultuna_max_Ir_bytreat_DF, "inputs_estiamtes_for_Carlos.ods", sheet = "max_Ir", append=T)
write_ods(Ultuna_mean_Is_bytreat_DF, "inputs_estiamtes_for_Carlos.ods", sheet = "mean_Is", append=T)
write_ods(Ultuna_min_Is_bytreat_DF, "inputs_estiamtes_for_Carlos.ods", sheet = "min_Is", append=T)
write_ods(Ultuna_max_Is_bytreat_DF, "inputs_estiamtes_for_Carlos.ods", sheet = "max_Is", append=T)

