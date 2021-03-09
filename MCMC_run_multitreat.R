

library(readODS)
library(rjags)
library(coda)
library(modeest)
library(viridis)
library(RColorBrewer)
library(imputeTS)

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


# reading the input data object
#calib_data_Ultuna<- readRDS(file = "Ultuna_input.rds")
calib_data_Ultuna<- readRDS(file = "Ultuna_input_1999.rds") #for excluding the maize years



#############################################
##### MCMC PART #############################
#############################################

## normalize re_clim on treatment G in Ultuna
#re_normalization<-rowMeans(aggregate(re_Ultuna_long, by=list(Ultuna_treat), FUN=mean)[,-1])[7]



#INIT RATIO PRIORS
mean_re_Ultuna= mean(as.numeric(re_Ultuna_long[21,]), na.rm=T)
mean_I_Ultuna= mean(Ultuna_yields_timeseries_long*10^-3, na.rm=T)
YSS_Ultuna<-mean_I_Ultuna*0.7/(0.8*mean_re_Ultuna)
OSS_Ultuna<-0.15*(mean_I_Ultuna*0.7/(0.0085*mean_re_Ultuna))
Init_prior_Ultuna<-1-YSS_Ultuna/OSS_Ultuna


library(modeest)

N.ADAPTS=300
N.RUNS=5000
N.BURNIN=200
sampling.nr=100


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


#N=newyears


model.ICBM.matrix_Ultuna<- jags.model('./JAGS_ICBM_3.1_Ultuna_multitreat.R',
                               data=calib_data_Ultuna,
                               n.chains = 1,
                               n.adapt = N.ADAPTS)


update(model.ICBM.matrix_Ultuna, n.iter=N.RUNS, n.thin=10, n.burnin=N.BURNIN)



parameter_list<-c("k1_ult",
                  "k2_ult",
                  "k2_ult_or",
                  "h_R_ult",
                  "h_S_ult",
                  "h_FYM_ult",
                  "h_PEA_ult",
                  "h_SAW_ult",
                  "h_SLU_ult",
                  "exudates_coeff",
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
                  "Y_R_Ultuna")

mcmc.array.ICBM.Ultuna<-(jags.samples(model.ICBM.matrix_Ultuna,parameter_list, sampling.nr))


###Ultuna
#Plot SOC and Input (roots and shoots) simulation
mcmc.list.SOC.Ultuna<-as.mcmc.list(mcmc.array.ICBM.Ultuna$Tot_Ultuna, chains=F)
mcmc.unlist.SOC.Ultuna<-mcmc.list.SOC.Ultuna[[1]]

mcmc.list.Ir.Ultuna<-as.mcmc.list(mcmc.array.ICBM.Ultuna$I_R_Ultuna, chains=F)
mcmc.unlist.Ir.Ultuna<-mcmc.list.Ir.Ultuna[[1]]

mcmc.list.Is.Ultuna<-as.mcmc.list(mcmc.array.ICBM.Ultuna$I_S_Ultuna, chains=F)
mcmc.unlist.Is.Ultuna<-mcmc.list.Is.Ultuna[[1]]

mcmc.list.fluxR.Ultuna<-as.mcmc.list(mcmc.array.ICBM.Ultuna$I_S_Ultuna, chains=F)
mcmc.unlist.fluxR.Ultuna<-mcmc.list.fluxR.Ultuna[[1]]

mcmc.list.fluxS.Ultuna<-as.mcmc.list(mcmc.array.ICBM.Ultuna$I_S_Ultuna, chains=F)
mcmc.unlist.fluxS.Ultuna<-mcmc.list.fluxS.Ultuna[[1]]


#newyears=10

#extract the array of the realizations*treatments*years
#dim<-dim(aggregate(Ultuna_SOC_timeseries_long, by=list(Ultuna_treat), FUN=mean)[,-1])
dim<-dim(calib_data_Ultuna$SOC_Ultuna)
#dim[2]<-dim[2]+N
Ultuna_prediction_array<-array( ,dim=c(sampling.nr,dim))
str(Ultuna_prediction_array)
for(j in 1:sampling.nr){
  length_treat<-dim[1]
  for(i in 1:(dim[2])){
    length_treat_add<-length_treat*i
    Ultuna_prediction_array[j,,i]<-mcmc.unlist.SOC.Ultuna[j,(1-length_treat+length_treat_add):length_treat_add]
  }
}


Ultuna_inputS_array<-array( ,dim=c(sampling.nr,dim))
for(j in 1:sampling.nr){
  length_treat<-dim[1]
  for(i in 1:(dim[2])){
    length_treat_add<-length_treat*i
    Ultuna_inputS_array[j,,i]<-mcmc.unlist.Is.Ultuna[j,(1-length_treat+length_treat_add):length_treat_add]
  }
}

Ultuna_inputR_array<-array( ,dim=c(sampling.nr,dim))
for(j in 1:sampling.nr){
  length_treat<-dim[1]
  for(i in 1:(dim[2])){
    length_treat_add<-length_treat*i
    Ultuna_inputR_array[j,,i]<-mcmc.unlist.Ir.Ultuna[j,(1-length_treat+length_treat_add):length_treat_add]
  }
}

Ultuna_fluxS_array<-array( ,dim=c(sampling.nr,dim))
for(j in 1:sampling.nr){
  length_treat<-dim[1]
  for(i in 1:(dim[2])){
    length_treat_add<-length_treat*i
    Ultuna_fluxS_array[j,,i]<-mcmc.unlist.fluxS.Ultuna[j,(1-length_treat+length_treat_add):length_treat_add]
  }
}

Ultuna_fluxR_array<-array( ,dim=c(sampling.nr,dim))
for(j in 1:sampling.nr){
  length_treat<-dim[1]
  for(i in 1:(dim[2])){
    length_treat_add<-length_treat*i
    Ultuna_fluxR_array[j,,i]<-mcmc.unlist.fluxR.Ultuna[j,(1-length_treat+length_treat_add):length_treat_add]
  }
}


#mean, min and max matrices by plot
Ultuna_mean_predictions<-(colMeans(Ultuna_prediction_array))
Ultuna_max_predictions<-apply(Ultuna_prediction_array, MARGIN=c(2,3), max)
Ultuna_min_predictions<-apply(Ultuna_prediction_array, MARGIN=c(2,3), min)
rownames(Ultuna_mean_predictions)<-as.factor(rownames(calib_data_Ultuna$SOC_Ultuna))
rownames(Ultuna_max_predictions)<-as.factor(rownames(calib_data_Ultuna$SOC_Ultuna))
rownames(Ultuna_min_predictions)<-as.factor(rownames(calib_data_Ultuna$SOC_Ultuna))

Ultuna_mean_Ir<-(colMeans(Ultuna_inputR_array))
Ultuna_max_Ir<-apply(Ultuna_inputR_array, MARGIN=c(2,3), max)
Ultuna_min_Ir<-apply(Ultuna_inputR_array, MARGIN=c(2,3), min)
rownames(Ultuna_mean_Ir)<-as.factor(rownames(calib_data_Ultuna$SOC_Ultuna))
rownames(Ultuna_max_Ir)<-as.factor(rownames(calib_data_Ultuna$SOC_Ultuna))
rownames(Ultuna_min_Ir)<-as.factor(rownames(calib_data_Ultuna$SOC_Ultuna))

Ultuna_mean_Is<-(colMeans(Ultuna_inputS_array))
Ultuna_max_Is<-apply(Ultuna_inputS_array, MARGIN=c(2,3), max)
Ultuna_min_Is<-apply(Ultuna_inputS_array, MARGIN=c(2,3), min)
rownames(Ultuna_mean_Is)<-as.factor(rownames(calib_data_Ultuna$SOC_Ultuna))
rownames(Ultuna_max_Is)<-as.factor(rownames(calib_data_Ultuna$SOC_Ultuna))
rownames(Ultuna_min_Is)<-as.factor(rownames(calib_data_Ultuna$SOC_Ultuna))

Ultuna_mean_fluxS<-(colMeans(Ultuna_fluxS_array))
Ultuna_max_fluxS<-apply(Ultuna_fluxS_array, MARGIN=c(2,3), max)
Ultuna_min_fluxS<-apply(Ultuna_fluxS_array, MARGIN=c(2,3), min)
rownames(Ultuna_mean_fluxS)<-as.factor(rownames(calib_data_Ultuna$SOC_Ultuna))
rownames(Ultuna_max_fluxS)<-as.factor(rownames(calib_data_Ultuna$SOC_Ultuna))
rownames(Ultuna_min_fluxS)<-as.factor(rownames(calib_data_Ultuna$SOC_Ultuna))

Ultuna_mean_fluxR<-(colMeans(Ultuna_fluxR_array))
Ultuna_max_fluxR<-apply(Ultuna_fluxR_array, MARGIN=c(2,3), max)
Ultuna_min_fluxR<-apply(Ultuna_fluxR_array, MARGIN=c(2,3), min)
rownames(Ultuna_mean_fluxR)<-as.factor(rownames(calib_data_Ultuna$SOC_Ultuna))
rownames(Ultuna_max_fluxR)<-as.factor(rownames(calib_data_Ultuna$SOC_Ultuna))
rownames(Ultuna_min_fluxR)<-as.factor(rownames(calib_data_Ultuna$SOC_Ultuna))


#mean, min and max matrices by treatment
#Ultuna_treat_vec<-LETTERS[yields_Ultuna[yields_Ultuna$year==1996,]$treat]
#Ultuna_treat_vec_old<-LETTERS[yields_Ultuna[yields_Ultuna$year==1956,]$treat]
Ultuna_treat_vec<-as.factor(rownames(calib_data_Ultuna$SOC_Ultuna))
Ultuna_mean_predictions_bytreat<-mat.or.vec(15,dim[2])
Ultuna_max_predictions_bytreat<-mat.or.vec(15,dim[2])
Ultuna_min_predictions_bytreat<-mat.or.vec(15,dim[2])
Ultuna_mean_Ir_bytreat<-mat.or.vec(15,dim[2])
Ultuna_max_Ir_bytreat<-mat.or.vec(15,dim[2])
Ultuna_min_Ir_bytreat<-mat.or.vec(15,dim[2])
Ultuna_mean_Is_bytreat<-mat.or.vec(15,dim[2])
Ultuna_max_Is_bytreat<-mat.or.vec(15,dim[2])
Ultuna_min_Is_bytreat<-mat.or.vec(15,dim[2])
Ultuna_mean_fluxS_bytreat<-mat.or.vec(15,dim[2])
Ultuna_max_fluxS_bytreat<-mat.or.vec(15,dim[2])
Ultuna_min_fluxS_bytreat<-mat.or.vec(15,dim[2])
Ultuna_mean_fluxR_bytreat<-mat.or.vec(15,dim[2])
Ultuna_max_fluxR_bytreat<-mat.or.vec(15,dim[2])
Ultuna_min_fluxR_bytreat<-mat.or.vec(15,dim[2])
Ultuna_measured_bytreat<-mat.or.vec(15,dim[2])
for(i in 1:15){
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
  Ultuna_measured_bytreat[i,]<-(calib_data_Ultuna$SOC_Ultuna[i,])
}


write.csv(t(Ultuna_mean_Ir_bytreat), file="Ultuna_mean_Ir_bytreat.csv")
write.csv(t(Ultuna_mean_Is_bytreat), file="Ultuna_mean_Is_bytreat.csv")

Ultuna_yields_table<-calib_data_Ultuna$Yields_cereals_Ultuna+calib_data_Ultuna$Yields_root_crops_Ultuna+calib_data_Ultuna$Yields_oilseeds_Ultuna+calib_data_Ultuna$Yields_maize_Ultuna
Ultuna_Ir_table<-as.data.frame(Ultuna_mean_Ir_bytreat)
Ultuna_Is_table<-as.data.frame(Ultuna_mean_Is_bytreat)
Ultuna_Pred_table<-as.data.frame(Ultuna_max_predictions_bytreat)

#colnames(Ultuna_yields_table)<-yields_Ultuna[yields_Ultuna$plot==3,]$crop
#rownames(Ultuna_yields_table)<-LETTERS[1:15]
# colnames(Ultuna_Ir_table)<-yields_Ultuna[yields_Ultuna$plot==3,]$crop
# rownames(Ultuna_Ir_table)<-LETTERS[1:15]
# colnames(Ultuna_Is_table)<-yields_Ultuna[yields_Ultuna$plot==3,]$crop
# rownames(Ultuna_Is_table)<-LETTERS[1:15]
# colnames(Ultuna_Pred_table)<-yields_Ultuna[yields_Ultuna$plot==3,]$crop
# rownames(Ultuna_Pred_table)<-LETTERS[1:15]

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
for(i in 1:15){
  Ultuna_Thomas_Ir_table[i,]<-as.numeric(Thomas_2011_inputs[Thomas_2011_inputs$treat==letters[i],]$`CR+CE70%`[1:54])*0.001
}

colnames(Ultuna_Thomas_Ir_table)<-colnames(calib_data_Ultuna$SOC_Ultuna)[1:54]
rownames(Ultuna_Thomas_Ir_table)<-LETTERS[1:15]

write_ods(Ultuna_Thomas_Ir_table, "Ultuna_input_check.ods", sheet="Thomas calculations IR (C,ton ha-1)", append=T, row_names = T)


png("input_comparison.png", width=4000, height=4000, res=300)
par(mfrow=c(4,4))
for(i in 1:15){
  plot(Ultuna_mean_Ir_bytreat[i,][1:54], type="l", col="darkorange", lwd=2, ylab="C Mg Ha-1", xlab="Year", main=paste("treatment", letters[i]), ylim=c(0,3))
  lines(as.numeric(Thomas_2011_inputs[Thomas_2011_inputs$treat==letters[i],]$`CR+CE70%`[1:54])*0.001, col="darkgreen", lwd=2)
  if(i ==1){legend("topleft", c("me", "Thomas"), col=c("darkorange", "darkgreen"), bty="n", lty=c(1), cex=1.5)}
  }
dev.off()

png("input_difference_me_vs_Thomas.png", width=4000, height=4000, res=300)
par(mfrow=c(4,4))
for(i in 1:15){
  barplot(Ultuna_mean_Ir_bytreat[i,][1:54]-as.numeric(Thomas_2011_inputs[Thomas_2011_inputs$treat==letters[i],]$`CR+CE70%`[1:54])*0.001, ylab="C Mg Ha-1", xlab="Year", main=paste("treatment", letters[i]), ylim=c(-2,2))
  text(9, 1.7, paste("mean : ", round(mean(Ultuna_mean_Ir_bytreat[i,][1:54]-as.numeric(Thomas_2011_inputs[Thomas_2011_inputs$treat==letters[i],]$`CR+CE70%`[1:54])*0.001, na.rm=T),3)))
}
dev.off()


png("input_difference_me_vs_Thomas_procent.png", width=4000, height=4000, res=300)
par(mfrow=c(4,4))
for(i in 1:15){
  barplot(100*(Ultuna_mean_Ir_bytreat[i,][1:54]-as.numeric(Thomas_2011_inputs[Thomas_2011_inputs$treat==letters[i],]$`CR+CE70%`[1:54])*0.001)/mean(Ultuna_mean_Ir_bytreat[i,][1:54]), ylab="% difference", xlab="Year", main=paste("treatment", letters[i]), ylim=c(-200,200))
  text(18, 185, paste("mean difference : ", round(mean(Ultuna_mean_Ir_bytreat[i,][1:54]-as.numeric(Thomas_2011_inputs[Thomas_2011_inputs$treat==letters[i],]$`CR+CE70%`[1:54])*0.001, na.rm=T)/mean(Ultuna_mean_Ir_bytreat[i,][1:54]),2)*100, "%"))
}
dev.off()





png("ICBM_predictions_Ultuna_specific.png", height=4500, width=4000, res=320)
par(mfrow=c(5,3))
for(i in 1:15){
  yearseq<-seq(from=colnames(calib_data_Ultuna$SOC_Ultuna)[1], to=as.numeric(colnames(calib_data_Ultuna$SOC_Ultuna)[1])+dim[2]-1)
  lastyear<-dim(Ultuna_mean_predictions_bytreat)[2]-1
  plot(yearseq, Ultuna_mean_predictions_bytreat[i,], ylim=c(20,125), type="l",  ylab=expression(paste("C (t ha" ^-1, ")")), xlab="Years", main=paste("Treatment", LETTERS[i]), col=Ultuna_treat_palette[i])
  polygon(c(yearseq,rev(yearseq)),
          c(Ultuna_max_predictions_bytreat[i,],rev(Ultuna_min_predictions_bytreat[i,])),
          col=Ultuna_treat_palette_alpha[i],border=Ultuna_treat_palette[i])
  lines(yearseq, Ultuna_mean_predictions_bytreat[i,], col="black", lty=1, lwd=1)
  #lines(yearseq[2:62], na_interpolation(Ultuna_measured_bytreat[i,1:61]), col="firebrick", lty=2, lwd=2)
  lines(yearseq[2:54], na_interpolation(Ultuna_measured_bytreat[i,1:53]), col=add.alpha("firebrick",0.2), lty=2, lwd=2)
  points(yearseq[2:54], (Ultuna_measured_bytreat[i,1:53]), col="firebrick", lty=2, lwd=2)
}
dev.off()



png("ICBM_Ir_Ultuna_specific.png", height=4500, width=4000, res=320)
par(mfrow=c(5,3))
for(i in 1:15){
  yearseq<-seq(from=colnames(calib_data_Ultuna$SOC_Ultuna)[1], to=as.numeric(colnames(calib_data_Ultuna$SOC_Ultuna)[1])+dim[2]-1)
  lastyear<-dim(Ultuna_mean_Ir_bytreat)[2]-1
  plot(yearseq, Ultuna_mean_Ir_bytreat[i,], ylim=c(0,2), type="l",  ylab=expression(paste("C (t ha" ^-1, ")")), xlab="Years", main=paste("Treatment", LETTERS[i]), col=Ultuna_treat_palette[i])
  polygon(c(yearseq,rev(yearseq)),
          c(Ultuna_max_Ir_bytreat[i,],rev(Ultuna_min_Ir_bytreat[i,])),
          col=Ultuna_treat_palette_alpha[i],border=Ultuna_treat_palette[i])
  lines(yearseq, Ultuna_mean_Ir_bytreat[i,], col="black", lty=1, lwd=1)
}
dev.off()

png("ICBM_Is_Ultuna_specific.png", height=4500, width=4000, res=320)
par(mfrow=c(5,3))
for(i in 1:15){
  yearseq<-seq(from=colnames(calib_data_Ultuna$SOC_Ultuna)[1], to=as.numeric(colnames(calib_data_Ultuna$SOC_Ultuna)[1])+dim[2]-1)
  lastyear<-dim(Ultuna_mean_Is_bytreat)[2]-1
  plot(yearseq, Ultuna_mean_Is_bytreat[i,], ylim=c(0,2), type="l",  ylab=expression(paste("C (t ha" ^-1, ")")), xlab="Years", main=paste("Treatment", LETTERS[i]), col=Ultuna_treat_palette[i])
  polygon(c(yearseq,rev(yearseq)),
          c(Ultuna_max_Is_bytreat[i,],rev(Ultuna_min_Is_bytreat[i,])),
          col=Ultuna_treat_palette_alpha[i],border=Ultuna_treat_palette[i])
  lines(yearseq, Ultuna_mean_Is_bytreat[i,], col="black", lty=1, lwd=1)
}
dev.off()

png("ICBM_fluxS_Ultuna_specific.png", height=4500, width=4000, res=320)
par(mfrow=c(5,3))
for(i in 1:15){
  yearseq<-seq(from=colnames(calib_data_Ultuna$SOC_Ultuna)[1], to=as.numeric(colnames(calib_data_Ultuna$SOC_Ultuna)[1])+dim[2]-1)
  lastyear<-dim(Ultuna_mean_fluxS_bytreat)[2]-1
  plot(yearseq, Ultuna_mean_fluxS_bytreat[i,], ylim=c(0,2), type="l",  ylab=expression(paste("C (t ha" ^-1, ")")), xlab="Years", main=paste("Treatment", LETTERS[i]), col=Ultuna_treat_palette[i])
  polygon(c(yearseq,rev(yearseq)),
          c(Ultuna_max_fluxS_bytreat[i,],rev(Ultuna_min_fluxS_bytreat[i,])),
          col=Ultuna_treat_palette_alpha[i],border=Ultuna_treat_palette[i])
  lines(yearseq, Ultuna_mean_fluxS_bytreat[i,], col="black", lty=1, lwd=1)
}
dev.off()

png("ICBM_fluxR_Ultuna_specific.png", height=4500, width=4000, res=320)
par(mfrow=c(5,3))
for(i in 1:15){
  yearseq<-seq(from=colnames(calib_data_Ultuna$SOC_Ultuna)[1], to=as.numeric(colnames(calib_data_Ultuna$SOC_Ultuna)[1])+dim[2]-1)
  lastyear<-dim(Ultuna_mean_fluxR_bytreat)[2]-1
  plot(yearseq, Ultuna_mean_fluxR_bytreat[i,], ylim=c(0.4,-0.5), type="l",  ylab=expression(paste("C (t ha" ^-1, ")")), xlab="Years", main=paste("Treatment", LETTERS[i]), col=Ultuna_treat_palette[i])
  polygon(c(yearseq,rev(yearseq)),
          c(Ultuna_max_fluxR_bytreat[i,],rev(Ultuna_min_fluxR_bytreat[i,])),
          col=Ultuna_treat_palette_alpha[i],border=Ultuna_treat_palette[i])
  lines(yearseq, Ultuna_mean_fluxR_bytreat[i,], col="black", lty=1, lwd=1)
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
rownames(Ultuna_mean_predictions)<-LETTERS[1:15]
rownames(Ultuna_max_predictions)<-LETTERS[1:15]
rownames(Ultuna_min_predictions)<-LETTERS[1:15]

#mean, min and max matrices by treatment
Ultuna_mean_predictions_bytreat_Y<-mat.or.vec(15,dim[2])
Ultuna_max_predictions_bytreat_Y<-mat.or.vec(15,dim[2])
Ultuna_min_predictions_bytreat_Y<-mat.or.vec(15,dim[2])
for(i in 1:15){
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
rownames(Ultuna_mean_predictions)<-LETTERS[1:15]
rownames(Ultuna_max_predictions)<-LETTERS[1:15]
rownames(Ultuna_min_predictions)<-LETTERS[1:15]

#mean, min and max matrices by treatment
Ultuna_mean_predictions_bytreat_O<-mat.or.vec(15,dim[2])
Ultuna_max_predictions_bytreat_O<-mat.or.vec(15,dim[2])
Ultuna_min_predictions_bytreat_O<-mat.or.vec(15,dim[2])
for(i in 1:15){
  Ultuna_mean_predictions_bytreat_O[i,]<-Ultuna_mean_predictions[i,]
  Ultuna_max_predictions_bytreat_O[i,]<-Ultuna_min_predictions[i,]
  Ultuna_min_predictions_bytreat_O[i,]<-Ultuna_max_predictions[i,]
}

png("ICBM_predictions_Ultuna_O_specific.png", height=4500, width=4000, res=320)
par(mfrow=c(5,3))
for(i in 1:15){
  yearseq<-seq(from=colnames(calib_data_Ultuna$SOC_Ultuna)[1], to=as.numeric(colnames(calib_data_Ultuna$SOC_Ultuna)[1])+dim[2]-1)
  lastyear<-dim(Ultuna_mean_predictions_bytreat_O)[2]-1
  plot(yearseq, Ultuna_mean_predictions_bytreat_O[i,], ylim=c(20,125), type="l",  ylab=expression(paste("C (t ha" ^-1, ")")), xlab="Years", main=paste("Treatment", LETTERS[i]), col=Ultuna_treat_palette[i])
  polygon(c(yearseq,rev(yearseq)),
          c(Ultuna_min_predictions_bytreat_O[i,],rev(Ultuna_max_predictions_bytreat_O[i,])),
          col=Ultuna_treat_palette_alpha[i],border=Ultuna_treat_palette[i])
  lines(yearseq, Ultuna_mean_predictions_bytreat_O[i,], col="black", lty=1, lwd=1)
}
dev.off()

png("ICBM_predictions_Ultuna_Y_specific.png", height=4500, width=4000, res=320)
par(mfrow=c(5,3))
for(i in 1:15){
  yearseq<-seq(from=colnames(calib_data_Ultuna$SOC_Ultuna)[1], to=as.numeric(colnames(calib_data_Ultuna$SOC_Ultuna)[1])+dim[2]-1)
  lastyear<-dim(Ultuna_mean_predictions_bytreat_Y)[2]-1
  plot(yearseq, Ultuna_mean_predictions_bytreat_Y[i,], ylim=c(0,6), type="l",  ylab=expression(paste("C (t ha" ^-1, ")")), xlab="Years", main=paste("Treatment", LETTERS[i]), col=Ultuna_treat_palette[i])
  polygon(c(yearseq,rev(yearseq)),
          c(Ultuna_min_predictions_bytreat_Y[i,],rev(Ultuna_max_predictions_bytreat_Y[i,])),
          col=Ultuna_treat_palette_alpha[i],border=Ultuna_treat_palette[i])
  lines(yearseq, Ultuna_mean_predictions_bytreat_Y[i,], col="black", lty=1, lwd=1)
}
dev.off()





###Plot parameter posteriors

parameters_palette_multi<-c("darkorange","cadetblue", "deepskyblue",
                            "cornsilk4","darkgreen","chocolate4",
                            "lightpink4","lightgoldenrod3","burlywood4","blueviolet","darkcyan",
                            brewer.pal(6, "BuPu")[2:6])
parameters_palette_multi_alpha<-add.alpha(parameters_palette_multi,0.7)


parameter_names<-c(expression(paste("k"[1])),
                   expression(paste("k"[2])),
                   expression(paste("k"[2],"_org")),
                   expression(paste("h"[R])),
                   expression(paste("h"[S])),
                   expression(paste("h"[FYM])),
                   expression(paste("h"[PEAT])),
                   expression(paste("h"[SAW])),
                   expression(paste("h"[SLU])),
                   "Exudates coefficient",
                   "Stubbles ratio",
                   "Stubbles ratio (maize)",
                   "S:R cereals","S:R root crops","S:R oilseeds", "S:R maize")

library(modeest)


#priors
library(truncnorm)
##Parameters

error_h<-0.1
limits_h<-0.3

#k1  <- rtruncnorm(5000, mean=0.8, sd=(0.8*0.05), a=0, b=3)
#k1  <- runif(5000, min=0.78, max=1)
k1  <- rnorm(5000,mean=0.6906054, sd=0.2225509)
k2    <- rtruncnorm(5000, mean=0.00605, sd=(0.00605*error_h), a=0.00605-0.00605*0.5, b=0.00605+0.00605*0.5)
k2_org    <- rtruncnorm(5000, mean=0.00605, sd=(0.00605*error_h), a=0.00605-0.00605*0.5, b=0.00605+0.00605*0.5)

h_S     <- rtruncnorm(5000, mean=0.15, sd=(0.15*error_h), a=0.15-0.15*limits_h, b=0.15+0.15*limits_h)
h_R     <- rtruncnorm(5000, mean=0.35, sd=(0.35*error_h), a=0.35-0.35*limits_h, b=0.35+0.35*limits_h)
h_FYM   <- rtruncnorm(5000, mean=0.27, sd=(0.27*error_h), a=0.27-limits_h, b=0.27+limits_h)
h_PEA   <- rtruncnorm(5000, mean=0.59, sd=(0.59*error_h), a=0.59-limits_h, b=0.59+limits_h)
h_SAW   <- rtruncnorm(5000, mean=0.25, sd=(0.25*error_h), a=0.25-limits_h, b=0.25+limits_h)
h_SLU   <- rtruncnorm(5000, mean=0.41, sd=(0.41*error_h), a=0.41-limits_h, b=0.41+limits_h)

stubbles_ratio_Ultuna <- rtruncnorm(5000, mean=0.04, sd=0.01, a=0.001, b=0.08)
stubbles_ratio_Ultuna_maize <- rtruncnorm(5000, mean=0.04, sd=0.01, a=0.001, b=0.08)

# Init_ratio_Ultuna <- rtruncnorm(5000, mean=0.9291667, sd=0.0001, a=0.9, b=0.95)
# Init_ratio_Lanna <- rtruncnorm(5000, mean=0.9291667, sd=0.0001, a=0.9, b=0.95)
exudates_coeff <- rtruncnorm(5000, mean=1.65, sd=1.65*0.1, a=1.65*0.95, b=1.65*1.05)

#root/shoot ratios priors
error_SR<-0.25
#SR_cereals    <- runif(5000, min=11-11*error_SR, max=11+11*error_SR)
#SR_root_crops <- runif(5000, min=30-30*error_SR, max=30+30*error_SR)
#SR_oilseeds   <- runif(5000, min=8-8*error_SR, max=8+8*error_SR)
#SR_maize   <- runif(5000, min=6.25-6.25*error_SR, max=6.25+6.25*error_SR)

SR_cereals    <- rtruncnorm(5000, mean=11, sd=(11*error_SR), a=11-+11*error_SR, b=11+11*error_SR)
SR_root_crops <- rtruncnorm(5000, mean=30, sd=30*error_SR, a=30-30*error_SR, b=30+30*error_SR)
SR_oilseeds   <- rtruncnorm(5000, mean=8, sd=8*error_SR, a=8-8*error_SR, b=8+8*error_SR)
SR_maize   <- rtruncnorm(5000, mean=6.25, sd=6.25*error_SR, a=6.25-6.25*error_SR, b=6.25+6.25*error_SR)


##Parameters
confidence=c(0.025,0.975)
parameter_values<-mat.or.vec(length(parameter_list[1:29]),4)
colnames(parameter_values)<-c("mode","mean","min","max")
rownames(parameter_values)<-parameter_list[1:29]

priors_list<-list(k1, #                           1
                  k2, #                           2
                  k2_org, #                       3
                  h_R, #                          4
                  h_S, #                          5
                  h_FYM, #                        6
                  h_PEA, #                        7
                  h_SAW, #                        8
                  h_SLU, #                        9
                  exudates_coeff, #               10
                  stubbles_ratio_Ultuna, #        11
                  stubbles_ratio_Ultuna_maize, #  12
                  SR_cereals, #                   13
                  SR_root_crops, #                14
                  SR_oilseeds, #                  15
                  SR_maize) #                     16

priors_list_d<-list()
for(i in 1:length(priors_list)){priors_list_d[[i]]<-density(priors_list[[i]])}
chain_list_multisite<-list()
MCMC_list_multisite<-list()
density_list_multisite<-list()
params_mode_multisite<-c()
params_meam_multisite<-c()
for(i in 1:(length(parameter_list[1:21]))){
  chain_list_multisite[[i]]<-as.mcmc.list(eval(parse(text=paste("mcmc.array.ICBM.Ultuna","$",parameter_list[i], sep=""))), chains=T)
  density_list_multisite[[i]]<-density(as.matrix(chain_list_multisite[[i]][1]))
  params_mode_multisite[i]<-mlv(as.matrix(chain_list_multisite[[i]][1]), method="Parzen")
  params_meam_multisite[i]<-mean(as.matrix(chain_list_multisite[[i]][1]))
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







png("ICBM_posteriors_multisite_specific.png", height=3000, width=4000, res=300)
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

i=2
range<-range(c(density_list_multisite[[i]]$x, priors_list_d[[i]]$x, priors_list_d[[i]]$x), density_list_multisite[[i+1]]$x)
mean<-mean(c(density_list_multisite[[i]]$x, priors_list_d[[i]]$x))
rangey<-range(c(density_list_multisite[[i]]$y, priors_list_d[[i]]$y, priors_list_d[[i]]$y))
meany<-mean(c(density_list_multisite[[i]]$y, priors_list_d[[i]]$y))

plot(density_list_multisite[[i]], main=parameter_names[i],xlim=c(range[1]-(range[1]*0.015), range[2]+(range[2]*0.015)), ylim=c(0, rangey[2]+(rangey[2]*0.2)))
polygon(density_list_multisite[[i]], col=add.alpha(parameters_palette_multi[i],0.8))
polygon(priors_list_d[[i]], col=add.alpha("lightgrey", 0.6), lty=2)
polygon(density_list_multisite[[i+1]], col=add.alpha(parameters_palette_multi[i+1],0.8))
polygon(priors_list_d[[i+1]], col=add.alpha("lightgrey", 0.6), lty=2)
legend("topleft", c("others", "peat/sludge"), bty="n", pch=c(16),  col=c(add.alpha(parameters_palette_multi[c(i, i+1)],0.8), add.alpha("lightgrey", 0.6)))



#H_r
i=4
parameter_list[c(i, i+6)]
range<-range(c(density_list_multisite[[i]]$x, priors_list_d[[i]]$x, priors_list_d[[i]]$x))
mean<-mean(c(density_list_multisite[[i]]$x, priors_list_d[[i]]$x))
rangey<-range(c(density_list_multisite[[i]]$y, priors_list_d[[i]]$y, priors_list_d[[i]]$y))
meany<-mean(c(density_list_multisite[[i]]$y, priors_list_d[[i]]$y))
plot(density_list_multisite[[i]], main=parameter_names[i],xlim=c(range[1]-(range[1]*0.015), range[2]+(range[2]*0.015)), ylim=c(0, rangey[2]+(rangey[2]*0.2)), col=NA)
polygon(density_list_multisite[[i]], col=add.alpha(parameters_palette_multi[i],0.8))
polygon(priors_list_d[[i]], col=add.alpha("lightgrey", 0.6), lty=2)


#H_s
i=5
parameter_list[c(i, i+6)]
range<-range(c(density_list_multisite[[i]]$x, priors_list_d[[i]]$x, priors_list_d[[i]]$x))
mean<-mean(c(density_list_multisite[[i]]$x, priors_list_d[[i]]$x))
rangey<-range(c(density_list_multisite[[i]]$y, priors_list_d[[i]]$y,priors_list_d[[i]]$y))
meany<-mean(c(density_list_multisite[[i]]$y, priors_list_d[[i]]$y))
plot(density_list_multisite[[i]], main=parameter_names[i],xlim=c(range[1]-(range[1]*0.015), range[2]+(range[2]*0.015)), ylim=c(0, rangey[2]+(rangey[2]*0.2)), col=NA)
polygon(density_list_multisite[[i]], col=add.alpha(parameters_palette_multi[i],0.8))
polygon(priors_list_d[[i]], col=add.alpha("lightgrey", 0.6), lty=2)


#H_FYM
i=6
parameter_list[c(i, i+6)]
range<-range(c(density_list_multisite[[i]]$x, priors_list_d[[i]]$x,priors_list_d[[i]]$x))
mean<-mean(c(density_list_multisite[[i]]$x, priors_list_d[[i]]$x))
rangey<-range(c(density_list_multisite[[i]]$y, priors_list_d[[i]]$y, priors_list_d[[i]]$y))
meany<-mean(c(density_list_multisite[[i]]$y, priors_list_d[[i]]$y))
plot(density_list_multisite[[i]], main=parameter_names[i],xlim=c(range[1]-(range[1]*0.015), range[2]+(range[2]*0.015)), ylim=c(0, rangey[2]+(rangey[2]*0.2)), col=NA)
polygon(density_list_multisite[[i]], col=add.alpha(parameters_palette_multi[i],0.8))
polygon(priors_list_d[[i]], col=add.alpha("lightgrey", 0.6), lty=2)


#H_PEA
i=7
parameter_list[c(i, i+6)]
range<-range(c(density_list_multisite[[i]]$x, priors_list_d[[i]]$x))
mean<-mean(c(density_list_multisite[[i]]$x, priors_list_d[[i]]$x))
rangey<-range(c(density_list_multisite[[i]]$y, priors_list_d[[i]]$y))
meany<-mean(c(density_list_multisite[[i]]$y, priors_list_d[[i]]$y))
plot(density_list_multisite[[i]], main=parameter_names[i],xlim=c(range[1]-(range[1]*0.015), range[2]+(range[2]*0.015)), ylim=c(0, rangey[2]+(rangey[2]*0.2)), col=NA)
polygon(density_list_multisite[[i]], col=add.alpha(parameters_palette_multi[i],0.8))
polygon(priors_list_d[[i]], col=add.alpha("lightgrey", 0.6), lty=2)

#H_SAW
i=8
parameter_list[c(i, i+6)]
range<-range(c(density_list_multisite[[i]]$x, priors_list_d[[i]]$x))
mean<-mean(c(density_list_multisite[[i]]$x, priors_list_d[[i]]$x))
rangey<-range(c(density_list_multisite[[i]]$y, priors_list_d[[i]]$y))
meany<-mean(c(density_list_multisite[[i]]$y, priors_list_d[[i]]$y))
plot(density_list_multisite[[i]], main=parameter_names[i],xlim=c(range[1]-(range[1]*0.015), range[2]+(range[2]*0.015)), ylim=c(0, rangey[2]+(rangey[2]*0.2)), col=NA)
polygon(density_list_multisite[[i]], col=add.alpha(parameters_palette_multi[i],0.8))
polygon(priors_list_d[[i]], col=add.alpha("lightgrey", 0.6), lty=2)

#H_SLU
i=9
parameter_list[c(i, i+4)]
range<-range(c(density_list_multisite[[i]]$x, priors_list_d[[i]]$x,priors_list_d[[i]]$x))
mean<-mean(c(density_list_multisite[[i]]$x, priors_list_d[[i]]$x))
rangey<-range(c(density_list_multisite[[i]]$y, priors_list_d[[i]]$y,priors_list_d[[i]]$y))
meany<-mean(c(density_list_multisite[[i]]$y, priors_list_d[[i]]$y))
plot(density_list_multisite[[i]], main=parameter_names[i],xlim=c(range[1]-(range[1]*0.015), range[2]+(range[2]*0.015)), ylim=c(0, rangey[2]+(rangey[2]*0.2)), col=NA)
polygon(density_list_multisite[[i]], col=add.alpha(parameters_palette_multi[i],0.8))
polygon(priors_list_d[[i]], col=add.alpha("lightgrey", 0.6), lty=2)


#ex. coefficient
range<-range(c(density_list_multisite[[10]]$x, priors_list_d[[10]]$x))
rangey<-range(c(density_list_multisite[[10]]$y, priors_list_d[[10]]$y))
plot(density_list_multisite[[10]], main=parameter_names[9],xlim=c(range[1]-(range[1]*0.015), range[2]+(range[2]*0.015)), ylim=c(0, rangey[2]+(rangey[2]*0.2)))
polygon(density_list_multisite[[10]], col=add.alpha(parameters_palette_multi[i],0.8))
polygon(priors_list_d[[10]], col=add.alpha("lightgrey", 0.6), lty=2)

#stubbles ratio
range<-range(c(density_list_multisite[[11]]$x, priors_list_d[[11]]$x))
mean<-mean(c(density_list_multisite[[i]]$x, priors_list_d[[11]]$x))
rangey<-range(c(density_list_multisite[[11]]$y, priors_list_d[[11]]$y))
plot(density_list_multisite[[11]], main=parameter_names[10],xlim=c(range[1]-(range[1]*0.015), range[2]+(range[2]*0.015)), ylim=c(0, rangey[2]+(rangey[2]*0.2)))
polygon(density_list_multisite[[11]], col=add.alpha(parameters_palette_multi[i],0.8))
polygon(priors_list_d[[11]], col=add.alpha("lightgrey", 0.1), lty=2)

#stubbles ratio maize
range<-range(c(density_list_multisite[[12]]$x, priors_list_d[[12]]$x))
mean<-mean(c(density_list_multisite[[i]]$x, priors_list_d[[12]]$x))
rangey<-range(c(density_list_multisite[[12]]$y, priors_list_d[[12]]$y))
plot(density_list_multisite[[12]], main=parameter_names[11],xlim=c(range[1]-(range[1]*0.015), range[2]+(range[2]*0.015)), ylim=c(0, rangey[2]+(rangey[2]*0.2)))
polygon(density_list_multisite[[12]], col=add.alpha(parameters_palette_multi[i],0.8))
polygon(priors_list_d[[12]], col=add.alpha("lightgrey", 0.1), lty=2)


for(i in 13:16){
  range<-range(c(density_list_multisite[[i]]$x, priors_list_d[[i]]$x))
  mean<-mean(c(density_list_multisite[[i]]$x, priors_list_d[[i]]$x))
  rangey<-range(c(density_list_multisite[[i]]$y, priors_list_d[[i]]$y))
  meany<-mean(c(density_list_multisite[[i]]$y, priors_list_d[[i]]$y))
  plot(density_list_multisite[[i]], main=parameter_names[i],xlim=c(range[1]-(range[1]*0.015), range[2]+(range[2]*0.015)), ylim=c(0, rangey[2]+(rangey[2]*0.2)))
  polygon(density_list_multisite[[i]], col=add.alpha(parameters_palette_multi[i],0.8))
  polygon(priors_list_d[[i]], col=add.alpha("lightgrey", 0.6), lty=2)
}


dev.off()















png("ICBM_posteriors_multisite_specific_portrait.png", height=4000, width=3000, res=330)
par(mfrow=c(5,3))

for(i in 1:2){
  range<-range(c(density_list_multisite[[i]]$x, priors_list_d[[i]]$x, priors_list_d[[i]]$x))
  mean<-mean(c(density_list_multisite[[i]]$x, priors_list_d[[i]]$x))
  rangey<-range(c(density_list_multisite[[i]]$y, priors_list_d[[i]]$y, priors_list_d[[i]]$y))
  meany<-mean(c(density_list_multisite[[i]]$y, priors_list_d[[i]]$y))

  plot(density_list_multisite[[i]], main=parameter_names[i],xlim=c(range[1]-(range[1]*0.015), range[2]+(range[2]*0.015)), ylim=c(0, rangey[2]+(rangey[2]*0.2)))
  polygon(density_list_multisite[[i]], col=add.alpha(parameters_palette_multi[i],0.8))
  polygon(priors_list_d[[i]], col=add.alpha("lightgrey", 0.6), lty=2)
  if(i==1){legend("topright", c("Ultuna", "Prior"), bty="n", pch=c(21),  col=c(add.alpha(parameters_palette_multi[i],0.8), add.alpha("lightgrey", 0.6)))}
}


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


#H_FYM
i=5
parameter_list[c(i, i+6)]
range<-range(c(density_list_multisite[[i]]$x, priors_list_d[[i]]$x,priors_list_d[[i]]$x))
mean<-mean(c(density_list_multisite[[i]]$x, priors_list_d[[i]]$x))
rangey<-range(c(density_list_multisite[[i]]$y, priors_list_d[[i]]$y, priors_list_d[[i]]$y))
meany<-mean(c(density_list_multisite[[i]]$y, priors_list_d[[i]]$y))
plot(density_list_multisite[[i]], main=parameter_names[i],xlim=c(range[1]-(range[1]*0.015), range[2]+(range[2]*0.015)), ylim=c(0, rangey[2]+(rangey[2]*0.2)), col=NA)
polygon(density_list_multisite[[i]], col=add.alpha(parameters_palette_multi[i],0.8))
polygon(priors_list_d[[i]], col=add.alpha("lightgrey", 0.6), lty=2)


#H_PEA
i=6
parameter_list[c(i, i+6)]
range<-range(c(density_list_multisite[[i]]$x, priors_list_d[[i]]$x))
mean<-mean(c(density_list_multisite[[i]]$x, priors_list_d[[i]]$x))
rangey<-range(c(density_list_multisite[[i]]$y, priors_list_d[[i]]$y))
meany<-mean(c(density_list_multisite[[i]]$y, priors_list_d[[i]]$y))
plot(density_list_multisite[[i]], main=parameter_names[i],xlim=c(range[1]-(range[1]*0.015), range[2]+(range[2]*0.015)), ylim=c(0, rangey[2]+(rangey[2]*0.2)), col=NA)
polygon(density_list_multisite[[i]], col=add.alpha(parameters_palette_multi[i],0.8))
polygon(priors_list_d[[i]], col=add.alpha("lightgrey", 0.6), lty=2)

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


#ex. coefficient
range<-range(c(density_list_multisite[[9]]$x, priors_list_d[[9]]$x))
rangey<-range(c(density_list_multisite[[9]]$y, priors_list_d[[9]]$y))
plot(density_list_multisite[[9]], main=parameter_names[9],xlim=c(range[1]-(range[1]*0.015), range[2]+(range[2]*0.015)), ylim=c(0, rangey[2]+(rangey[2]*0.2)))
polygon(density_list_multisite[[9]], col=add.alpha(parameters_palette_multi[i],0.8))
polygon(priors_list_d[[9]], col=add.alpha("lightgrey", 0.6), lty=2)

#stubbles ratio
range<-range(c(density_list_multisite[[10]]$x, priors_list_d[[10]]$x))
mean<-mean(c(density_list_multisite[[i]]$x, priors_list_d[[10]]$x))
rangey<-range(c(density_list_multisite[[10]]$y, priors_list_d[[10]]$y))
plot(density_list_multisite[[10]], main=parameter_names[10],xlim=c(range[1]-(range[1]*0.015), range[2]+(range[2]*0.015)), ylim=c(0, rangey[2]+(rangey[2]*0.2)))
polygon(density_list_multisite[[10]], col=add.alpha(parameters_palette_multi[i],0.8))
polygon(priors_list_d[[10]], col=add.alpha("lightgrey", 0.1), lty=2)


#stubbles ratio maize
range<-range(c(density_list_multisite[[11]]$x, priors_list_d[[11]]$x))
mean<-mean(c(density_list_multisite[[i]]$x, priors_list_d[[11]]$x))
rangey<-range(c(density_list_multisite[[11]]$y, priors_list_d[[11]]$y))
plot(density_list_multisite[[11]], main=parameter_names[11],xlim=c(range[1]-(range[1]*0.015), range[2]+(range[2]*0.015)), ylim=c(0, rangey[2]+(rangey[2]*0.2)))
polygon(density_list_multisite[[11]], col=add.alpha(parameters_palette_multi[i],0.8))
polygon(priors_list_d[[11]], col=add.alpha("lightgrey", 0.1), lty=2)




for(i in 12:15){
  range<-range(c(density_list_multisite[[i]]$x, priors_list_d[[i]]$x))
  mean<-mean(c(density_list_multisite[[i]]$x, priors_list_d[[i]]$x))
  rangey<-range(c(density_list_multisite[[i]]$y, priors_list_d[[i]]$y))
  meany<-mean(c(density_list_multisite[[i]]$y, priors_list_d[[i]]$y))
  plot(density_list_multisite[[i]], main=parameter_names[i],xlim=c(range[1]-(range[1]*0.015), range[2]+(range[2]*0.015)), ylim=c(0, rangey[2]+(rangey[2]*0.2)))
  polygon(density_list_multisite[[i]], col=add.alpha(parameters_palette_multi[i],0.8))
  polygon(priors_list_d[[i]], col=add.alpha("lightgrey", 0.6), lty=2)
}



dev.off()



#residuals
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
for(i in 1:15){
plot(residuals[i,], axes=FALSE, ylab="mean residuals", xlab="", main=paste("Treatment", LETTERS[i]), pch=16, ylim=range(residuals, na.rm=T))
axis(1, at=seq(1, length(residuals_type)), labels= residuals_type[i,], las=2, cex.axis=0.5)
axis(2) #default way
abline(h=0, lty=2)
box()
}
dev.off()


png("residuals_freescale.png", width = 5000, height=4000, res=300)
par(mfrow=c(5,3))
for(i in 1:15){
  plot(residuals[i,], axes=FALSE, ylab="mean residuals", xlab="", main=paste("Treatment", LETTERS[i]), pch=16)
  axis(1, at=seq(1, length(residuals_type)), labels= residuals_type[i,], las=2, cex.axis=0.5)
  axis(2) #default way
  abline(h=0, lty=2)
  box()
}
dev.off()

png("residuals_mean.png", width = 2000, height=2000, res=300)
barplot(rowMeans(residuals, na.rm=T))
dev.off()

png("residuals_scatterplots.png", width = 4000, height=2000, res=300)
par(mfrow=c(1,2))
yields<-(calib_data_Ultuna$Yields_cereals_Ultuna+calib_data_Ultuna$Yields_root_crops_Ultuna+calib_data_Ultuna$Yields_oilseeds_Ultuna+calib_data_Ultuna$Yields_maize_Ultuna)
plot(unlist(yields),unlist(residuals), main="Residuals vs yields", col=Ultuna_crop_palette[as.numeric(as.factor(unlist(residuals_type)))], pch=as.numeric(as.factor(unlist(residuals_type))), xlab="Yields", ylab="Residuals")
legend("bottomright", levels(as.factor(unlist(residuals_type))), pch=seq(1,4), col = Ultuna_crop_palette, bty="n")

plot(unlist(calib_data_Ultuna$re_Ultuna),unlist(residuals), main="Residuals vs reclim", col=Ultuna_crop_palette[as.numeric(as.factor(unlist(residuals_type)))], pch=as.numeric(as.factor(unlist(residuals_type))), xlab="reclim", ylab="Residuals")
legend("bottomright", levels(as.factor(unlist(residuals_type))), pch=seq(1,4), col = Ultuna_crop_palette, bty="n")
dev.off()


