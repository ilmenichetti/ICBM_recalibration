library(RColorBrewer)
library(viridis)
library(hydroGOF)

library(zoo)
library(RColorBrewer)
library(tsoutliers)
library(plyr)
library(imputeTS)
library(lubridate)


library(rjags)
library(coda)
library(modeest)

library(readxl)
library(writexl)
library(readODS)


Lanna_treat_palette<-rev(inferno(10)[2:10])
Ultuna_treat_palette<-rev(viridis(16)[2:16])

SR=c(11, 8, 2, 1/0.0339)
names(SR)<-c("Cereals", "Oilseeds", "Forage", "Roots")

#root plant (like potatoes) values from:
#Bolinder, M. A., K?tterer, T., Poeplau, C., B?rjesson, G., & Parent, L. E. (2015). Net primary productivity and below-ground crop residue inputs for root crops: Potato ( Solanum tuberosum L.) and sugar beet ( Beta vulgaris L.). Canadian Journal of Soil Science, 95(2), 87-93. https://doi.org/10.4141/cjss-2014-091
# the ratio is in fact yields to roots

#h values for amendments from
#Peltre, C., Christensen, B. T., Dragon, S., Icard, C., K?tterer, T., & Houot, S. (2012). RothC simulation of carbon accumulation in soil after repeated application of widely different organic amendments. Soil Biology and Biochemistry, 52, 49-60. https://doi.org/10.1016/j.soilbio.2012.03.023
HF<-c(27.4, 11,48.3,50.3,61.2,65.2,50)/100
HF_sd<-(c(2,4.8,3.7,7.2,15.6,2.9,7.5)/100)*3
names(HF)<-c("h_STR", "h_GM", "h_SAW", "h_FYM", "h_SLU", "h_PEA", "Rescale")
names(HF_sd)<-names(HF)
HF_max<-HF+HF_sd
HF_min<-HF-HF_sd
HF_max[HF_max>1]=1
HF_min[HF_min<0]=0


#list for choosing which factor to apply
amend_treat_list_Ultuna<-c("h_a",   "h_a",   "h_a",  "h_a",   "h_a",
                    "h_STR", "h_STR", "h_GM", "h_PEA", "h_FYM",
                    "h_FYM", "h_SAW", "h_PEA","h_SAW", "h_SLU")

### FUNCTIONS LOADING
#the fitness functions come from this package
library(hydroGOF)

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
Lanna_treat_palette_alpha<-add.alpha(Lanna_treat_palette,0.5)



####################################
####### DATA LOADING ULTUNA ########
####################################

yields_Ultuna<-read_excel("./input_data_0_36.xlsx", sheet=1)
yields_Ultuna_thomas<-read_excel("./Thomas_yields_Ultuna_no1959.xlsx", sheet=1)
re_Ultuna<-read.csv("../re_clim/Ultuna_FRAME56_calculated_re.csv")[,2:16]
years_diff_Ultuna<-tail(unique(yields_Ultuna$year),1)-unique(yields_Ultuna$year)[1]

yields_Ultuna<-cbind(yields_Ultuna, yields_Ultuna$grain_dm_kg_ha+yields_Ultuna$straw_dm_kg_ha)
colnames(yields_Ultuna)[15]<-"tot_yields_thomas"

for(i in 1:length(unique(yields_Ultuna_thomas$year))){
          for(j in 1:15){
            yields_Ultuna[yields_Ultuna$year==unique(yields_Ultuna$year)[i] & yields_Ultuna$treat == j, ]$tot_yields_thomas<-
            rep(yields_Ultuna_thomas[yields_Ultuna_thomas$year==unique(yields_Ultuna$year)[i] & yields_Ultuna_thomas$treat == letters[j], ]$total_yield,4)
          }
        }

mean(yields_Ultuna$tot_yields_thomas, na.rm=T)
Ultuna_years_first<-unique(yields_Ultuna$year)[1]
Ultuna_years_last<-tail(unique(yields_Ultuna$year),1)


# prepare the bulk density vectors by linear interpolation
Ultuna_BD_start=1.44
Ultuna_BD_end_list=c(1.43,1.40,1.28,1.21,1.25,1.38,1.21,1.34,1.12,1.24,1.20,1.28,1.05,1.023,1.02)
Ultuna_BD_timeseries<-(mat.or.vec(15, (years_diff_Ultuna+1)))
Ultuna_BD_timeseries[,1]<-Ultuna_BD_start
Ultuna_BD_timeseries[,(years_diff_Ultuna+1)]<-Ultuna_BD_end_list
Ultuna_BD_timeseries[Ultuna_BD_timeseries==0]<-NA
for(i in 1:15){
  Ultuna_BD_timeseries[i,]<-na.approx(Ultuna_BD_timeseries[i,])
}
rownames(Ultuna_BD_timeseries)<-LETTERS[1:15]

Ultuna_BD_timeseries_long<-(mat.or.vec(60, (years_diff_Ultuna+1)))
rownames(Ultuna_BD_timeseries_long)<-LETTERS[yields_Ultuna[yields_Ultuna$year==1956,]$treat]
for(i in 1:60){
  treat<-rownames(Ultuna_BD_timeseries_long)[i]
  Ultuna_BD_timeseries_long[i,]<-Ultuna_BD_timeseries[rownames(Ultuna_BD_timeseries)==treat,]
}



#fill in some missing data with a loop, year 1976 is missing entrely
yields_Ultuna<-as.data.frame(yields_Ultuna,useNA=T)
which_to_fill<-yields_Ultuna$year==1976
for(i in 1:60){
  yields_Ultuna[yields_Ultuna$year==1976 & yields_Ultuna$plot==i,7:10]=as.numeric(colMeans(yields_Ultuna[yields_Ultuna$plot==i,7:10], na.rm=T))
}



#rework all the time series in a table
Ultuna_average_CN<-ddply(yields_Ultuna, .(year,plot),summarise,
                    N    = length(yields_Ultuna),
                    C_conc = mean(C_conc, na.rm=T),
                    C_conc_SD = sd(C_conc, na.rm=T),
                    N_conc = mean(N_conc, na.rm=T),
                    N_conc_SD = sd(N_conc, na.rm=T))

Ultuna_C_timeseries_long<-(mat.or.vec(60, (years_diff_Ultuna+1)))
rownames(Ultuna_C_timeseries_long)<-yields_Ultuna[yields_Ultuna$year==1956,]$plot
treat_vec<-c()
for(i in 1:60){
  plot<-rownames(Ultuna_C_timeseries_long)[i]
  Ultuna_C_timeseries_long[i,]<-yields_Ultuna[yields_Ultuna$plot==plot,]$C_conc
  treat_vec[i]<-LETTERS[unique(yields_Ultuna[yields_Ultuna$plot==plot,]$treat)]
}

Ultuna_N_timeseries_long<-(mat.or.vec(60, (years_diff_Ultuna+1)))
rownames(Ultuna_N_timeseries_long)<-yields_Ultuna[yields_Ultuna$year==1956,]$plot
for(i in 1:60){
  plot<-rownames(Ultuna_N_timeseries_long)[i]
  Ultuna_N_timeseries_long[i,]<-yields_Ultuna[yields_Ultuna$plot==plot,]$N_conc
}

Ultuna_yields_timeseries_long<-(mat.or.vec(60, (years_diff_Ultuna+1)))
rownames(Ultuna_yields_timeseries_long)<-yields_Ultuna[yields_Ultuna$year==1956,]$plot
for(i in 1:60){
  plot<-rownames(Ultuna_yields_timeseries_long)[i]
  Ultuna_yields_timeseries_long[i,]<-yields_Ultuna[yields_Ultuna$plot==plot,]$grain_dm_kg_ha+yields_Ultuna[yields_Ultuna$plot==plot,]$straw_dm_kg_ha
}

Ultuna_yields_timeseries_long_Thomas<-(mat.or.vec(60, (years_diff_Ultuna+1)))
rownames(Ultuna_yields_timeseries_long_Thomas)<-yields_Ultuna[yields_Ultuna$year==1956,]$treat
for(i in 1:60){
  treat<-rownames(Ultuna_yields_timeseries_long_Thomas)[i]
  ave<-mean(yields_Ultuna_thomas[yields_Ultuna_thomas$treat==letters[as.numeric(treat)],]$total_yield)
  Ultuna_yields_timeseries_long_Thomas[i,]<-c(yields_Ultuna_thomas[yields_Ultuna_thomas$treat==letters[as.numeric(treat)],]$total_yield, ave,ave)
}

plot(yields_Ultuna_thomas[yields_Ultuna_thomas$treat=="b",]$total_yield, type="l")
plot(Ultuna_yields_timeseries_long_Thomas[5,], type="l")



### CREATE THE TIME VECTORS for reference
Ultuna_date_vector_days<-seq.Date(from = as.Date(ISOdate(head(Ultuna_years_first,1),1,1)),
                           to=as.Date(ISOdate(head(Ultuna_years_last,1),12,31)), by="days")
Ultuna_date_vector_year<-seq.Date(from = as.Date(ISOdate(head(Ultuna_years_first,1),1,1)),
                           to=as.Date(ISOdate(head(Ultuna_years_last,1),12,31)), by="year")


# treatment  of the outliers and interpolation
Ultuna_C_timeseries_long_interp<-Ultuna_C_timeseries_long
for(i in 1:60){
  Ultuna_C_timeseries_long_interp[i,]<-na_interpolation(Ultuna_C_timeseries_long_interp[i,], option="stine")
}

library(viridis)
Ultuna_treat_palette<-rev(viridis(15))

treat_num<-as.numeric(as.factor(treat_vec))
pch_list<-seq(1:15)

png("Ultuna_C_procent_data.png", height=1800, width=2500, res=380)
plot(year(Ultuna_date_vector_year),Ultuna_C_timeseries_long_interp[1,],
     col=Ultuna_treat_palette[as.numeric(as.factor(treat_vec[1]))], type="l", ylim=c(0,5.7), xlab="year", ylab="C (%)")
for(i in 1:60){
  lines(year(Ultuna_date_vector_year),Ultuna_C_timeseries_long_interp[i,], col=Ultuna_treat_palette[as.numeric(as.factor(treat_vec))[i]], type="l")
  points(year(Ultuna_date_vector_year),Ultuna_C_timeseries_long[i,],
         col=Ultuna_treat_palette[as.numeric(as.factor(treat_vec))[i]], pch=pch_list[as.numeric(as.factor(treat_vec))[i]])
}
legend("topleft", LETTERS[1:15], col=Ultuna_treat_palette, bty="n", lty=1, pch=seq(1:15),cex=0.7)
dev.off()


#depth of soil interested
depth=20#depth of soil considered in cm

#Ctot in tons per ha
Ultuna_SOC_timeseries_long<-Ultuna_C_timeseries_long_interp*depth*Ultuna_BD_timeseries_long


png("Ultuna_SOC_data.png", height=1800, width=2500, res=380)
plot(year(Ultuna_date_vector_year),Ultuna_SOC_timeseries_long[1,],
     col=Ultuna_treat_palette[as.numeric(as.factor(treat_vec[1]))], type="l", ylim=c(0,100), xlab="year", ylab=expression('SOC (Mg ha' ^ -1~')'))
for(i in 1:60){
  lines(year(Ultuna_date_vector_year),Ultuna_SOC_timeseries_long[i,], col=Ultuna_treat_palette[as.numeric(as.factor(treat_vec))[i]], type="l")
  points(year(Ultuna_date_vector_year),Ultuna_SOC_timeseries_long[i,],
         col=Ultuna_treat_palette[as.numeric(as.factor(treat_vec))[i]], pch=pch_list[as.numeric(as.factor(treat_vec))[i]])
}
legend("topleft", LETTERS[1:15], col=Ultuna_treat_palette, bty="n", lty=1, pch=seq(1:15),cex=0.7)
dev.off()


#re matrix
 re_Ultuna_long<-as.data.frame(mat.or.vec(60, (years_diff_Ultuna+1)))
 for(i in 0:14){
   re_crop<-re_Ultuna[i+2,]

   re_Ultuna_long[1+(i*4),]<-re_crop
   re_Ultuna_long[2+(i*4),]<-re_crop
   re_Ultuna_long[3+(i*4),]<-re_crop
   re_Ultuna_long[4+(i*4),]<-re_crop
 }
 colnames(re_Ultuna_long)<-re_Ultuna[1,]+1955
#re_Ultuna_long<-t(re_Ultuna[1:(years_diff_Ultuna+1),-1])
colnames(re_Ultuna_long)<- year(Ultuna_date_vector_year)




#table for treatments
treat_Ultuna<-data.frame(yields_Ultuna$treat ,yields_Ultuna$plot, yields_Ultuna$amended_bool)
colnames(treat_Ultuna)<-c("treat","plot","manure")
treat_Ultuna_long_bool<-(mat.or.vec(60, (years_diff_Ultuna+1)))
rownames(treat_Ultuna_long_bool)<-yields_Ultuna[yields_Ultuna$year==1956,]$plot
for(i in 1:60){
  treat<-treat_Ultuna[i,]$treat
  treat_Ultuna_long_bool[i,]<- (amend_treat_list_Ultuna)[treat]
}


for(i in 1:dim(treat_Ultuna_long_bool)[2]){
  if(i %% 2 == 0){treat_Ultuna_long_bool[,i]<-treat_Ultuna_long_bool[,i]}else{
    treat_Ultuna_long_bool[,i]<-"h_a"
  }
}


#table for manure
h_FYM_Ultuna_long<- (treat_Ultuna_long_bool=="h_FYM")*1900*2
h_GM_Ultuna_long<- (treat_Ultuna_long_bool=="h_GM")*1760*2
h_PEA_Ultuna_long<- (treat_Ultuna_long_bool=="h_PEA")*1970*2
h_SAW_Ultuna_long<- (treat_Ultuna_long_bool=="h_SAW")*1840*2
h_SLU_Ultuna_long<- (treat_Ultuna_long_bool=="h_SLU")*1840*2
h_STR_Ultuna_long<- (treat_Ultuna_long_bool=="h_STR")*1770*2



# create one yields vector for each class of Martin's S:R
cereals_vector<-yields_Ultuna$crop=="fodder" | yields_Ultuna$crop=="winter_small_grains" | yields_Ultuna$crop=="spring_small_grains"
roots_vector<-yields_Ultuna$crop=="root_crop"
oil_vector<-yields_Ultuna$crop=="winter_oil_seeds"
maize_vector<-yields_Ultuna$crop=="maize"


# table of which crop is where
Ultuna_crop_table<-cereals_vector
Ultuna_crop_table[cereals_vector]<-"cereals"
Ultuna_crop_table[roots_vector]<-"root_crop"
Ultuna_crop_table[oil_vector]<-"winter_oil_seeds"
Ultuna_crop_table[maize_vector]<-"maize"



#table for cereals
#yields_Ultuna_cereals<-data.frame(yields_Ultuna$treat ,yields_Ultuna$plot, yields_Ultuna$grain_dm_kg_ha+yields_Ultuna$straw_dm_kg_ha)
yields_Ultuna_cereals<-data.frame(yields_Ultuna$treat ,yields_Ultuna$plot, yields_Ultuna$`yields_Ultuna$grain_dm_kg_ha + yields_Ultuna$straw_dm_kg_ha`)
colnames(yields_Ultuna_cereals)<-c("treat","plot","yields")
yields_Ultuna_cereals[yields_Ultuna_cereals$plot==5,]$yields


yields_Ultuna_cereals$yields[!cereals_vector]<-0
Ultuna_yields_timeseries_long_cereals<-(mat.or.vec(60, (years_diff_Ultuna+1)))
rownames(Ultuna_yields_timeseries_long_cereals)<-yields_Ultuna[yields_Ultuna$year==1956,]$plot
for(i in 1:60){
  plot<-rownames(Ultuna_yields_timeseries_long)[i]
  Ultuna_yields_timeseries_long_cereals[i,]<- na.approx(yields_Ultuna_cereals[yields_Ultuna_cereals$plot==plot,]$yields)
}



#table for roots
#yields_Ultuna_roots<-data.frame(yields_Ultuna$treat ,yields_Ultuna$plot, yields_Ultuna$grain_dm_kg_ha+yields_Ultuna$straw_dm_kg_ha)
yields_Ultuna_roots<-data.frame(yields_Ultuna$treat ,yields_Ultuna$plot, yields_Ultuna$`yields_Ultuna$grain_dm_kg_ha + yields_Ultuna$straw_dm_kg_ha`)
colnames(yields_Ultuna_roots)<-c("treat","plot","yields")
yields_Ultuna_roots$yields[!roots_vector]<-0
Ultuna_yields_timeseries_long_roots<-(mat.or.vec(60, (years_diff_Ultuna+1)))
rownames(Ultuna_yields_timeseries_long_roots)<-yields_Ultuna[yields_Ultuna$year==1956,]$plot
for(i in 1:60){
  plot<-rownames(Ultuna_yields_timeseries_long)[i]
  Ultuna_yields_timeseries_long_roots[i,]<- na.approx(yields_Ultuna_roots[yields_Ultuna_roots$plot==plot,]$yields)
}

#table for oil
#yields_Ultuna_oil<-data.frame(yields_Ultuna$treat ,yields_Ultuna$plot, yields_Ultuna$grain_dm_kg_ha+yields_Ultuna$straw_dm_kg_ha)
yields_Ultuna_oil<-data.frame(yields_Ultuna$treat ,yields_Ultuna$plot, yields_Ultuna$`yields_Ultuna$grain_dm_kg_ha + yields_Ultuna$straw_dm_kg_ha`)
colnames(yields_Ultuna_oil)<-c("treat","plot","yields")
yields_Ultuna_oil$yields[!oil_vector]<-0
Ultuna_yields_timeseries_long_oil<-(mat.or.vec(60, (years_diff_Ultuna+1)))
rownames(Ultuna_yields_timeseries_long_oil)<-yields_Ultuna[yields_Ultuna$year==1956,]$plot
for(i in 1:60){
  plot<-rownames(Ultuna_yields_timeseries_long)[i]
  Ultuna_yields_timeseries_long_oil[i,]<- na.approx(yields_Ultuna_oil[yields_Ultuna_oil$plot==plot,]$yields)
}

#table for maize
#yields_Ultuna_maize<-data.frame(yields_Ultuna$treat ,yields_Ultuna$plot, yields_Ultuna$grain_dm_kg_ha+yields_Ultuna$straw_dm_kg_ha)
yields_Ultuna_maize<-data.frame(yields_Ultuna$treat ,yields_Ultuna$plot, yields_Ultuna$`yields_Ultuna$grain_dm_kg_ha + yields_Ultuna$straw_dm_kg_ha`)
colnames(yields_Ultuna_maize)<-c("treat","plot","yields")

yields_Ultuna_maize$yields[!maize_vector]<-0
Ultuna_yields_timeseries_long_maize<-(mat.or.vec(60, (years_diff_Ultuna+1)))
rownames(Ultuna_yields_timeseries_long_maize)<-yields_Ultuna[yields_Ultuna$year==1956,]$plot
for(i in 1:60){
  plot<-rownames(Ultuna_yields_timeseries_long)[i]
  Ultuna_yields_timeseries_long_maize[i,]<- na.approx(yields_Ultuna_maize[yields_Ultuna_maize$plot==plot,]$yields)
}





#function to extend the inputs by N years
append_NA<-function(mat, N)
{
  NA_vec<-as.numeric(rep(NA, dim(mat)[1]))
  mat_new<-cbind(mat,matrix(rep(NA_vec,each=N), ncol=N, byrow=TRUE))
  return(mat_new)
}

append_mean<-function(mat, N)
{
  mean_vec<-rowMeans(mat)
  mat_new<-cbind(mat,matrix(rep(mean_vec,each=N), ncol=N, byrow=TRUE))
  return(mat_new)
}

#estimate the error in the bare fallows
#this term will be added to the calibration additively
Ultuna_treat<-LETTERS[yields_Ultuna[yields_Ultuna$year==1956,]$treat]
error_BF_Ultuna<-mean(Ultuna_yields_timeseries_long[Ultuna_treat=="B",])*0.05



newyears=10

#find the RMSE for Ultuna, to be used as variance of SOC time series
Ultuna_ls_list<-list()
Ultuna_RMSE_byplot<-c()
Ultuna_max_err_byplot<-c()
for(i in 1:dim(Ultuna_SOC_timeseries_long)[1]){
  Ultuna_ls_list[[i]]<-lm(Ultuna_SOC_timeseries_long[i,] ~ seq(1:length(Ultuna_SOC_timeseries_long[i,])))
  summary(Ultuna_ls_list[[i]])
  RSS <- c(crossprod(Ultuna_ls_list[[i]]$residuals))
  MSE <- RSS / length(Ultuna_ls_list[[i]]$residuals)
  Ultuna_RMSE_byplot[i] <- sqrt(MSE)
  Ultuna_max_err_byplot[i]<- max(Ultuna_ls_list[[i]]$residuals)
  
}

################### Read eqiovalent soil depth calculations #########################
Ultuna_SOC_timeseries_long2<-read_excel("./Ultuna_SOC_equiv_depth_TK.xlsx", sheet="Ultuna_SOC_stocks")[,-1]
Ultuna_SOC_timeseries_long_pre<-Ultuna_SOC_timeseries_long
Ultuna_SOC_timeseries_long<-as.matrix(Ultuna_SOC_timeseries_long2)




png("Thomas_SOC_data.png", height=1800, width=2500, res=380)
plot(colnames(Ultuna_SOC_timeseries_long),Ultuna_SOC_timeseries_long[1,],
     col=Ultuna_treat_palette[as.numeric(as.factor(treat_vec[1]))], type="l", ylim=c(0,120), xlab="year", ylab=expression('SOC (Mg ha' ^ -1~')'))
for(i in 1:60){
  lines(colnames(Ultuna_SOC_timeseries_long),na.approx(Ultuna_SOC_timeseries_long[i,]), col=Ultuna_treat_palette[as.numeric(as.factor(Ultuna_treat))[i]], type="l")
  points(colnames(Ultuna_SOC_timeseries_long),(Ultuna_SOC_timeseries_long[i,]),
         col=Ultuna_treat_palette[as.numeric(as.factor(Ultuna_treat))[i]], pch=pch_list[as.numeric(as.factor(Ultuna_treat))[i]])
}
legend("topleft", LETTERS[1:15], col=Ultuna_treat_palette, bty="n", lty=1, pch=seq(1:15),cex=0.7)
dev.off()



##add back NAs where the original measurement is missing (no interpolation)
Ultuna_SOC_timeseries_long[is.na(Ultuna_C_timeseries_long)]<-NA
aggregate(Ultuna_SOC_timeseries_long, by=list(Ultuna_treat), FUN=mean)[,-1]


calib_data_Ultuna = list(##Ultuna loop
  'SOC_Ultuna' = append_NA(aggregate(Ultuna_SOC_timeseries_long, by=list(Ultuna_treat), FUN=mean)[,-1], N=newyears),
  'SOC_init_Ultuna' = aggregate(Ultuna_SOC_timeseries_long[,1], by=list(Ultuna_treat), FUN=mean)[,-1],
  'error_SOC_Ultuna'= aggregate(Ultuna_max_err_byplot, by=list(Ultuna_treat), FUN=mean)[,-1],
  'error_SOC_multiplier_Ultuna' = c(1, rep(1,14)),
  'Yields_cereals_Ultuna'= append_mean(aggregate(Ultuna_yields_timeseries_long_cereals*10^-3, by=list(Ultuna_treat), FUN=mean)[,-1], N=newyears),
  'Yields_root_crops_Ultuna'= append_mean(aggregate(Ultuna_yields_timeseries_long_roots*10^-3, by=list(Ultuna_treat), FUN=mean)[,-1], N=newyears),
  'Yields_oilseeds_Ultuna'= append_mean(aggregate(Ultuna_yields_timeseries_long_oil*10^-3, by=list(Ultuna_treat), FUN=mean)[,-1], N=newyears),
  'Yields_maize_Ultuna'= append_mean(aggregate(Ultuna_yields_timeseries_long_maize*10^-3, by=list(Ultuna_treat), FUN=mean)[,-1], N=newyears),
  'I_FYM_Ultuna'= append_mean(aggregate(h_FYM_Ultuna_long*10^-3, by=list(Ultuna_treat), FUN=mean)[,-1], N=newyears),
  'I_GM_Ultuna'= append_mean(aggregate(h_GM_Ultuna_long*10^-3, by=list(Ultuna_treat), FUN=mean)[,-1], N=newyears),
  'I_PEA_Ultuna'= append_mean(aggregate(h_PEA_Ultuna_long*10^-3, by=list(Ultuna_treat), FUN=mean)[,-1], N=newyears),
  'I_SAW_Ultuna'= append_mean(aggregate(h_SAW_Ultuna_long*10^-3, by=list(Ultuna_treat), FUN=mean)[,-1], N=newyears),
  'I_SLU_Ultuna'= append_mean(aggregate(h_SLU_Ultuna_long*10^-3, by=list(Ultuna_treat), FUN=mean)[,-1], N=newyears),
  'I_STR_Ultuna'= append_mean(aggregate(h_STR_Ultuna_long*10^-3, by=list(Ultuna_treat), FUN=mean)[,-1], N=newyears),
  'N_Ultuna' = dim(append_NA(aggregate(Ultuna_SOC_timeseries_long, by=list(Ultuna_treat), FUN=mean)[,-1], N=newyears))[2],
  're_Ultuna'= append_mean(as.data.frame(t(re_Ultuna)), N=newyears),
  'J_Ultuna'=dim(append_NA(aggregate(Ultuna_SOC_timeseries_long, by=list(Ultuna_treat), FUN=mean)[,-1], N=newyears))[1])

Ultuna_date_vector_year_long<-seq.Date(from = as.Date(ISOdate(head(Ultuna_years_first,1),1,1)),
                                       to=as.Date(ISOdate(head(Ultuna_years_last+newyears,1),12,31)), by="year")
names(calib_data_Ultuna[[1]])<-year(Ultuna_date_vector_year_long)
names(calib_data_Ultuna[[5]])<-year(Ultuna_date_vector_year_long)
names(calib_data_Ultuna[[6]])<-year(Ultuna_date_vector_year_long)
names(calib_data_Ultuna[[7]])<-year(Ultuna_date_vector_year_long)
names(calib_data_Ultuna[[8]])<-year(Ultuna_date_vector_year_long)
names(calib_data_Ultuna[[9]])<-year(Ultuna_date_vector_year_long)
names(calib_data_Ultuna[[10]])<-year(Ultuna_date_vector_year_long)
names(calib_data_Ultuna[[11]])<-year(Ultuna_date_vector_year_long)
names(calib_data_Ultuna[[12]])<-year(Ultuna_date_vector_year_long)
names(calib_data_Ultuna[[13]])<-year(Ultuna_date_vector_year_long)
names(calib_data_Ultuna[[14]])<-year(Ultuna_date_vector_year_long)
names(calib_data_Ultuna[[16]])<-year(Ultuna_date_vector_year_long)

rownames(calib_data_Ultuna[[1]])<-LETTERS[1:15]
names(calib_data_Ultuna[[2]])<-LETTERS[1:15]
names(calib_data_Ultuna[[3]])<-LETTERS[1:15]
names(calib_data_Ultuna[[4]])<-LETTERS[1:15]
rownames(calib_data_Ultuna[[5]])<-LETTERS[1:15]
rownames(calib_data_Ultuna[[6]])<-LETTERS[1:15]
rownames(calib_data_Ultuna[[7]])<-LETTERS[1:15]
rownames(calib_data_Ultuna[[8]])<-LETTERS[1:15]
rownames(calib_data_Ultuna[[9]])<-LETTERS[1:15]
rownames(calib_data_Ultuna[[10]])<-LETTERS[1:15]
rownames(calib_data_Ultuna[[11]])<-LETTERS[1:15]
rownames(calib_data_Ultuna[[12]])<-LETTERS[1:15]
rownames(calib_data_Ultuna[[13]])<-LETTERS[1:15]
rownames(calib_data_Ultuna[[14]])<-LETTERS[1:15]
rownames(calib_data_Ultuna[[16]])<-LETTERS[1:15]
names(calib_data_Ultuna[[15]])<-"Simulation length"
names(calib_data_Ultuna[[17]])<-"number of treatments"


plot(colnames(Ultuna_SOC_timeseries_long),Ultuna_SOC_timeseries_long[1,],
     col=Ultuna_treat_palette[as.numeric(as.factor(treat_vec[1]))], type="l", ylim=c(0,120), xlab="year", ylab=expression('SOC (Mg ha' ^ -1~')'))
for(i in 1:60){
  lines(colnames(Ultuna_SOC_timeseries_long),na.approx(Ultuna_SOC_timeseries_long[i,]), col=Ultuna_treat_palette[as.numeric(as.factor(Ultuna_treat))[i]], type="l")
  points(colnames(Ultuna_SOC_timeseries_long),(Ultuna_SOC_timeseries_long[i,]),
         col=Ultuna_treat_palette[as.numeric(as.factor(Ultuna_treat))[i]], pch=pch_list[as.numeric(as.factor(Ultuna_treat))[i]])
}
legend("topleft", LETTERS[1:15], col=Ultuna_treat_palette, bty="n", lty=1, pch=seq(1:15),cex=0.7)




png("Ultuna_SOC_data_used.png", height=1800, width=2500, res=380)

plot(year(Ultuna_date_vector_year_long),calib_data_Ultuna[[1]][1,],
     col=Ultuna_treat_palette[as.numeric(as.factor(treat_vec[1]))], type="l", ylim=c(0,120), xlab="year", ylab=expression('SOC (Mg ha' ^ -1~')'))
for(i in 1:15){
  lines(year(Ultuna_date_vector_year),na.approx(as.numeric(calib_data_Ultuna[[1]][i,])), col=Ultuna_treat_palette[[i]], type="l")
  points(year(Ultuna_date_vector_year_long),calib_data_Ultuna[[1]][i,],
         col=Ultuna_treat_palette[[i]], pch=pch_list[[i]])
}
legend("topleft", LETTERS[1:15], col=Ultuna_treat_palette, bty="n", lty=1, pch=seq(1:15),cex=0.7)

dev.off()

data_descr_Ultuna<-data.frame(names(calib_data_Ultuna), 
                       c("SOC, 15 treatments x lenght of simulation years",
                       "SOC initialization values, 15 treatments",
                       "SOC mean errors, 15 treatments", 
                       "SOC error weight",
                       "Cereals yields, 15 treatments x lenght of simulation years",
                       "Root crops yields, 15 treatments x lenght of simulation years",
                       "Oil seeds yields, 15 treatments x lenght of simulation years",
                       "Maize yields, 15 treatments x lenght of simulation years",
                       "Inputs farmyard manure, 15 treatments x lenght of simulation years",
                       "Inputs green manure, 15 treatments x lenght of simulation years",
                       "Inputs peat, 15 treatments x lenght of simulation years",
                       "Inputs sawdust, 15 treatments x lenght of simulation years",
                       "Inputs sludges, 15 treatments x lenght of simulation years",
                       "Inputs straw, 15 treatments x lenght of simulation years",
                       "Lenght of the simulation",
                       "re_clim, 15 treatments x lenght of simulation years",
                       "Number of treatments"))
colnames(data_descr_Ultuna)<-c("Sheet name", "Description")


write_ods(as.data.frame(calib_data_Ultuna[[1]]), "Input_Ultuna.ods", sheet=names(calib_data_Ultuna)[1], col_names = T, row_names = T)
write_ods(as.data.frame(calib_data_Ultuna[[2]]), "Input_Ultuna.ods", sheet=names(calib_data_Ultuna)[2], append = T, col_names = T, row_names = T)
write_ods(as.data.frame(calib_data_Ultuna[[3]]), "Input_Ultuna.ods", sheet=names(calib_data_Ultuna)[3], append = T, col_names = T, row_names = T)
write_ods(as.data.frame(calib_data_Ultuna[[4]]), "Input_Ultuna.ods", sheet=names(calib_data_Ultuna)[4], append = T, col_names = T, row_names = T)
write_ods(as.data.frame(calib_data_Ultuna[[5]]), "Input_Ultuna.ods", sheet=names(calib_data_Ultuna)[5], append = T, col_names = T, row_names = T)
write_ods(as.data.frame(calib_data_Ultuna[[6]]), "Input_Ultuna.ods", sheet=names(calib_data_Ultuna)[6], append = T, col_names = T, row_names = T)
write_ods(as.data.frame(calib_data_Ultuna[[7]]), "Input_Ultuna.ods", sheet=names(calib_data_Ultuna)[7], append = T, col_names = T, row_names = T)
write_ods(as.data.frame(calib_data_Ultuna[[8]]), "Input_Ultuna.ods", sheet=names(calib_data_Ultuna)[8], append = T, col_names = T, row_names = T)
write_ods(as.data.frame(calib_data_Ultuna[[9]]), "Input_Ultuna.ods", sheet=names(calib_data_Ultuna)[9], append = T, col_names = T, row_names = T)
write_ods(as.data.frame(calib_data_Ultuna[[10]]), "Input_Ultuna.ods", sheet=names(calib_data_Ultuna)[10], append = T, col_names = T, row_names = T)
write_ods(as.data.frame(calib_data_Ultuna[[11]]), "Input_Ultuna.ods", sheet=names(calib_data_Ultuna)[11], append = T, col_names = T, row_names = T)
write_ods(as.data.frame(calib_data_Ultuna[[12]]), "Input_Ultuna.ods", sheet=names(calib_data_Ultuna)[12], append = T, col_names = T, row_names = T)
write_ods(as.data.frame(calib_data_Ultuna[[13]]), "Input_Ultuna.ods", sheet=names(calib_data_Ultuna)[13], append = T, col_names = T, row_names = T)
write_ods(as.data.frame(calib_data_Ultuna[[14]]), "Input_Ultuna.ods", sheet=names(calib_data_Ultuna)[14], append = T, col_names = T, row_names = T)
write_ods(as.data.frame(calib_data_Ultuna[[15]]), "Input_Ultuna.ods", sheet=names(calib_data_Ultuna)[15], append = T, col_names = T, row_names = T)
write_ods(as.data.frame(calib_data_Ultuna[[16]]), "Input_Ultuna.ods", sheet=names(calib_data_Ultuna)[16], append = T, col_names = T, row_names = T)
write_ods(as.data.frame(calib_data_Ultuna[[17]]), "Input_Ultuna.ods", sheet=names(calib_data_Ultuna)[17], append = T, col_names = T, row_names = T)















####################################
####### DATA LOADING LANNA ########
####################################

yields_Lanna<-read_excel("./input_data_0_36.xlsx", sheet=3)
re_Lanna<-read.csv("../re_clim/Lanna_calculated_re.csv")[,2:10]
years_diff_Lanna<-tail(unique(yields_Lanna$year),1)-unique(yields_Lanna$year)[1]

mean(unlist(re_Lanna))

Lanna_years_first<-unique(yields_Lanna$year)[1]
Lanna_years_last<-tail(unique(yields_Lanna$year),1)

#substituting the treatment names with Ultuna equivalents
#yields_Lanna$treat<-yields_Lanna$Treatment_ult

# prepare the bulk density vectors
Lanna_BD_start=1.33
Lanna_BD_end_list=c(1.38,1.36,1.33,1.33,1.32,1.30,1.27,1.33,1.38)
Lanna_BD_timeseries<-(mat.or.vec(9, (years_diff_Lanna+1)))
Lanna_BD_timeseries[,1]<-Lanna_BD_start
Lanna_BD_timeseries[,(years_diff_Lanna+1)]<-Lanna_BD_end_list
Lanna_BD_timeseries[Lanna_BD_timeseries==0]<-NA
for(i in 1:9){
  Lanna_BD_timeseries[i,]<-na.approx(Lanna_BD_timeseries[i,])
}
rownames(Lanna_BD_timeseries)<-LETTERS[1:9]

Lanna_BD_timeseries_long<-(mat.or.vec(36, (years_diff_Lanna+1)))
rownames(Lanna_BD_timeseries_long)<-LETTERS[yields_Lanna[yields_Lanna$year==1996,]$treat]
for(i in 1:36){
  treat<-rownames(Lanna_BD_timeseries_long)[i]
  Lanna_BD_timeseries_long[i,]<-Lanna_BD_timeseries[rownames(Lanna_BD_timeseries)==treat,]
}
#rownames(Lanna_BD_timeseries_long)<-LETTERS[yields_Lanna[yields_Lanna$year==1956,]$plot]

#rework all the time series in a table
Lanna_average_CN<-ddply(yields_Lanna, .(year,plot),summarise,
                         N    = length(yields_Lanna),
                         C_conc = mean(C_conc, na.rm=T),
                         C_conc_SD = sd(C_conc, na.rm=T),
                         N_conc = mean(N_conc, na.rm=T),
                         N_conc_SD = sd(N_conc, na.rm=T))

Lanna_average_CN_bytreat<-ddply(yields_Lanna, .(year,treat),summarise,
                        N    = length(yields_Lanna),
                        C_conc = mean(C_conc, na.rm=T),
                        C_conc_SD = sd(C_conc, na.rm=T),
                        N_conc = mean(N_conc, na.rm=T),
                        N_conc_SD = sd(N_conc, na.rm=T))

Lanna_C_timeseries_long<-(mat.or.vec(36, (years_diff_Lanna+1)))
rownames(Lanna_C_timeseries_long)<-yields_Lanna[yields_Lanna$year==1996,]$plot
treat_vec_Lanna<-c()
for(i in 1:36){
  plot<-rownames(Lanna_C_timeseries_long)[i]
  Lanna_C_timeseries_long[i,]<-yields_Lanna[yields_Lanna$plot==plot,]$C_conc
  treat_vec_Lanna[i]<-LETTERS[yields_Lanna$Treatment_ult[i]]
}

Lanna_N_timeseries_long<-(mat.or.vec(36, (years_diff_Lanna+1)))
rownames(Lanna_N_timeseries_long)<-yields_Lanna[yields_Lanna$year==1996,]$plot
for(i in 1:36){
  plot<-rownames(Lanna_N_timeseries_long)[i]
  Lanna_N_timeseries_long[i,]<-yields_Lanna[yields_Lanna$plot==plot,]$N_conc
}

Lanna_yields_timeseries_long<-(mat.or.vec(36, (years_diff_Lanna+1)))
rownames(Lanna_yields_timeseries_long)<-yields_Lanna[yields_Lanna$year==1996,]$plot
for(i in 1:36){
  plot<-rownames(Lanna_yields_timeseries_long)[i]
  Lanna_yields_timeseries_long[i,]<-yields_Lanna[yields_Lanna$plot==plot,]$grain_dm_kg_ha+yields_Lanna[yields_Lanna$plot==plot,]$straw_dm_kg_ha
}


### CREATE THE TIME VECTORS for reference
Lanna_date_vector_days<-seq.Date(from = as.Date(ISOdate(head(Lanna_years_first,1),1,1)),
                                  to=as.Date(ISOdate(head(Lanna_years_last,1),12,31)), by="days")
Lanna_date_vector_year<-seq.Date(from = as.Date(ISOdate(head(Lanna_years_first,1),10,1)),
                                  to=as.Date(ISOdate(head(Lanna_years_last,1),12,31)), by="year")

year(Lanna_date_vector_year)==unique(yields_Lanna$year)

# treatment  of the outliers and interpolation
Lanna_C_timeseries_long_interp<-Lanna_C_timeseries_long
for(i in 1:36){
  Lanna_C_timeseries_long_interp[i,]<-na_interpolation(Lanna_C_timeseries_long_interp[i,], option="stine")
}


treat_num_Lanna<-as.numeric(as.factor(treat_vec_Lanna))
pch_list<-seq(1:15)

png("Lanna_C_procent_data.png", height=1800, width=2500, res=380)
plot(year(Lanna_date_vector_year),Lanna_C_timeseries_long_interp[1,],
     col=Lanna_treat_palette[as.numeric(as.factor(treat_vec_Lanna[1]))], type="l", ylim=c(0,5.7), xlab="year", ylab="C (%)")
for(i in 1:36){
  lines(year(Lanna_date_vector_year),Lanna_C_timeseries_long_interp[i,], col=Lanna_treat_palette[as.numeric(as.factor(treat_vec_Lanna))[i]], type="l")
  points(year(Lanna_date_vector_year),Lanna_C_timeseries_long[i,],
         col=Lanna_treat_palette[as.numeric(as.factor(treat_vec))[i]], pch=pch_list[as.numeric(as.factor(treat_vec_Lanna))[i]])
}
legend("topleft", unique(treat_vec_Lanna), col=Lanna_treat_palette[unique(as.numeric(as.factor(treat_vec_Lanna)))], bty="n", lty=1, pch=seq(1:15),cex=0.7)
dev.off()


#depth of soil interested
depth=20#depth of soil considered in cm

#Ctot in tons per ha
Lanna_SOC_timeseries_long<-Lanna_C_timeseries_long_interp*depth*Lanna_BD_timeseries_long

 re_Lanna_long<-as.data.frame(mat.or.vec(36, (years_diff_Lanna+1)))
 for(i in 0:8){
   re_crop<-re_Lanna[i+2,]

   re_Lanna_long[1+(i*4),]<-re_crop
   re_Lanna_long[2+(i*4),]<-re_crop
   re_Lanna_long[3+(i*4),]<-re_crop
   re_Lanna_long[4+(i*4),]<-re_crop
 }
 colnames(re_Lanna_long)<-re_Lanna[1,]+1995

#re_Lanna_long<-t(re_Lanna[1:(years_diff_Lanna+1),-1])
colnames(re_Lanna_long)<- year(Lanna_date_vector_year)


#table for manure
manure_Lanna<-data.frame(yields_Lanna$treat ,yields_Lanna$plot, yields_Lanna$amendment_C_kg_ha)
colnames(manure_Lanna)<-c("treat","plot","manure")
manure_Lanna_long<-(mat.or.vec(36, (years_diff_Lanna+1)))
for(i in 1:36){
  plot<-rownames(Lanna_C_timeseries_long_interp)[i]
  manure_Lanna_long[i,]<- manure_Lanna[manure_Lanna$plot==plot,]$manure
}

amend_treat_list_Lanna<-c("h_a","h_a","h_a","h_STR","h_FYM","h_SLU","h_SLU","h_CO","h_a")

#table for treatments
treat_Lanna<-data.frame(yields_Lanna$treat ,yields_Lanna$plot, yields_Lanna$amended_bool)
colnames(treat_Lanna)<-c("treat","plot","manure")
treat_Lanna_long_bool<-(mat.or.vec(36, (years_diff_Lanna+1)))
rownames(treat_Lanna_long_bool)<-yields_Lanna[yields_Lanna$year==1996,]$plot
for(i in 1:36){
  treat<-treat_Lanna[i,]$treat
  treat_Lanna_long_bool[i,]<- (amend_treat_list_Lanna)[treat]
}


#tables for specific amendments
h_FYM_Lanna_long<-mat.or.vec(36,22)
h_GM_Lanna_long<- mat.or.vec(36,22)
h_PEA_Lanna_long<- mat.or.vec(36,22)
h_SAW_Lanna_long<- mat.or.vec(36,22)
h_SLU_Lanna_long<- mat.or.vec(36,22)
h_STR_Lanna_long<- mat.or.vec(36,22)
h_CO_Lanna_long<- mat.or.vec(36,22)
for(i in 1:36){h_FYM_Lanna_long[i,]<-manure_Lanna_long[i,]*as.numeric(treat_Lanna_long_bool[i,]=="h_FYM")*1000}
for(i in 1:36){h_GM_Lanna_long[i,]<-manure_Lanna_long[i,]*as.numeric(treat_Lanna_long_bool[i,]=="h_GM")*1000}
for(i in 1:36){h_PEA_Lanna_long[i,]<-manure_Lanna_long[i,]*as.numeric(treat_Lanna_long_bool[i,]=="h_PEA")*1000}
for(i in 1:36){h_SAW_Lanna_long[i,]<-manure_Lanna_long[i,]*as.numeric(treat_Lanna_long_bool[i,]=="h_SAW")*1000}
for(i in 1:36){h_SLU_Lanna_long[i,]<-manure_Lanna_long[i,]*as.numeric(treat_Lanna_long_bool[i,]=="h_SLU")*1000}
for(i in 1:36){h_STR_Lanna_long[i,]<-manure_Lanna_long[i,]*as.numeric(treat_Lanna_long_bool[i,]=="h_STR")*1000}
for(i in 1:36){h_CO_Lanna_long[i,]<-manure_Lanna_long[i,]*as.numeric(treat_Lanna_long_bool[i,]=="h_CO")*1000}



# create one yields vector for each class of Martin's S:R
cereals_vector<-yields_Lanna$crop_id=="fodder" | yields_Lanna$crop_id=="winter_small_grains" | yields_Lanna$crop_id=="spring_small_grains"
roots_vector<-yields_Lanna$crop_id=="root_crop"
oil_vector<-yields_Lanna$crop_id=="winter_oil_seeds"

#table for cereals
yields_Lanna_cereals<-data.frame(yields_Lanna$treat ,yields_Lanna$plot, yields_Lanna$grain_dm_kg_ha+yields_Lanna$straw_dm_kg_ha, yields_Lanna$straw_dm_kg_ha)
colnames(yields_Lanna_cereals)<-c("treat","plot","yields", "straw")
yields_Lanna_cereals$yields[!cereals_vector]<-0
Lanna_yields_timeseries_long_cereals<-(mat.or.vec(36, (years_diff_Lanna+1)))
rownames(Lanna_yields_timeseries_long_cereals)<-yields_Lanna[yields_Lanna$year==1996,]$plot
Lanna_yields_timeseries_long_cereals_straw<-(mat.or.vec(36, (years_diff_Lanna+1)))
rownames(Lanna_yields_timeseries_long_cereals_straw)<-yields_Lanna[yields_Lanna$year==1996,]$plot
for(i in 1:36){
  plot<-rownames(Lanna_yields_timeseries_long)[i]
  Lanna_yields_timeseries_long_cereals[i,]<- yields_Lanna_cereals[yields_Lanna_cereals$plot==plot,]$yields
  Lanna_yields_timeseries_long_cereals_straw[i,]<- yields_Lanna_cereals[yields_Lanna_cereals$plot==plot,]$straw
}


write.csv(Lanna_BD_timeseries_long, file="Lanna_BD_timeseries.csv")

#################### Read eqiovalent soil depth calculations #########################

Lanna_SOC_timeseries_long2<-read_excel("./Lanna_SOC_equiv_depth.xlsx", sheet="Lanna_SOC_stocks")[,-1]
Lanna_SOC_timeseries_long<-as.matrix(Lanna_SOC_timeseries_long2)

##add back NAs where the original measurement is missing (no interpolation)
Lanna_SOC_timeseries_long[is.na(Lanna_C_timeseries_long)]<-NA


#find the RMSE for Lanna, to be used as variance of SOC time series
Lanna_ls_list<-list()
Lanna_RMSE_byplot<-c()
Lanna_max_err_byplot<-c()
for(i in 1:dim(Lanna_SOC_timeseries_long)[1]){
  Lanna_ls_list[[i]]<-lm(Lanna_SOC_timeseries_long[i,] ~ seq(1:length(Lanna_SOC_timeseries_long[i,])))
  summary(Lanna_ls_list[[i]])
  RSS <- c(crossprod(Lanna_ls_list[[i]]$residuals))
  MSE <- RSS / length(Lanna_ls_list[[i]]$residuals)
  Lanna_RMSE_byplot[i] <- sqrt(MSE)
  Lanna_max_err_byplot[i]<- max(Lanna_ls_list[[i]]$residuals)
}


#take out the growing trend from the Lanna BF
Lanna_SOC_timeseries_long_backup<-Lanna_SOC_timeseries_long

min_Lanna_point<-min(which.min(Lanna_SOC_timeseries_long[1,]),
                    which.min(Lanna_SOC_timeseries_long[2,]),
                    which.min(Lanna_SOC_timeseries_long[3,]),
                    which.min(Lanna_SOC_timeseries_long[4,]))

Lanna_SOC_timeseries_long[1:12,min_Lanna_point:dim(Lanna_SOC_timeseries_long)[2]]<-NA
Lanna_SOC_timeseries_long[33:36,min_Lanna_point:dim(Lanna_SOC_timeseries_long)[2]]<-NA


Lanna_treat=c(rep("A",4), #1
              rep("C",4), #2
              rep("D",4), #3
              rep("G",4), #4 Maybe F
              rep("J",4), #5
              rep("O",4), #6
              rep("O",4), #7
              rep("X",4), #8 Compost, must be treated separately
              rep("B",4)) #9





png("Lanna_re_used.png", height=1800, width=2500, res=350)
re_Lanna_plot<- aggregate(re_Lanna_long, by=list(Lanna_treat), FUN=mean)[,-1]
#re_Lanna_plot<- re_Lanna_long
lanna_tr_nr<-unique(as.numeric(as.factor(Lanna_treat)))
plot(as.numeric(names(re_Lanna_plot[1,])), as.numeric(re_Lanna_plot[1,]), type="l", ylim=c(0.3,2.5), col=Lanna_treat_palette[1], ylab="re", xlab="years")
for (i in 1:9){
  lines(as.numeric(names(re_Lanna_plot[i,])), as.numeric(re_Lanna_plot[i,]),  col=Lanna_treat_palette[lanna_tr_nr[i]])
  points(as.numeric(names(re_Lanna_plot[i,])),as.numeric(re_Lanna_plot[i,]), cex=0.5,
         col=Lanna_treat_palette[lanna_tr_nr[i]], pch=pch_list[i])

}
legend("bottomright",unique(Lanna_treat), col=Lanna_treat_palette[lanna_tr_nr], bty="n", lty=1, pch=seq(1:9),cex=0.6)
dev.off()


png("Ultuna_re_used.png", height=1800, width=2500, res=350)
re_Ultuna_plot<- calib_data_Ultuna$re_Ultuna
#re_Ultuna_plot<- re_Ultuna_long
plot(as.numeric(names(re_Ultuna_plot[1,])), (re_Ultuna_plot[1,]), type="l", ylim=c(0.3,1.4), col=Ultuna_treat_palette[1], ylab="re", xlab="years")
for (i in 1:15){
  lines(as.numeric(names(re_Ultuna_plot[i,])), as.numeric(re_Ultuna_plot[i,]),  col=Ultuna_treat_palette[i])
  points(as.numeric(names(re_Ultuna_plot[i,])),as.numeric(re_Ultuna_plot[i,]), cex=0.5,
         col=Ultuna_treat_palette[i], pch=seq(1:15)[i])

}
legend("bottomright", LETTERS[1:15], col=Ultuna_treat_palette, bty="n", lty=1, pch=seq(1:15),cex=0.6)
dev.off()


png("Inputs_used_barplot.png", height=2800, width=2500, res=350)
par(mfrow=c(2,1))
barplot(Ultuna_yields_timeseries_long[,], ylim=c(0, max(Ultuna_yields_timeseries_long)), ylab="Ultuna inputs", beside=T)
barplot(Lanna_yields_timeseries_long[,], ylim=c(0, max(Lanna_yields_timeseries_long)), ylab="Lanna inputs", beside=T)
dev.off()





lanna_order<-as.numeric(as.factor(unique(Lanna_treat))) #to rearrange Lanna according to the list of treatments (otherwise the averages are messed up)

lanna_which<-c(1,2,3,4,5,6,7,8) #to esclude some treatments in Lanna


Yields_Lanna_test<-aggregate(Lanna_yields_timeseries_long_cereals*10^-3, by=list(Lanna_treat), FUN=mean)[lanna_order,]


##LANNA loop
calib_data_Lanna = list('SOC_Lanna' = append_NA(aggregate(Lanna_SOC_timeseries_long, by=list(Lanna_treat), FUN=mean)[lanna_order,][,-1], N=newyears)[lanna_which,],
                        'SOC_init_Lanna' = c(aggregate(Lanna_SOC_timeseries_long[,1], by=list(Lanna_treat), FUN=mean, na.rm=T)[lanna_order,][,-1])[lanna_which],
                        'error_SOC_Lanna'= c(aggregate(Lanna_max_err_byplot, by=list(Lanna_treat), FUN=max)[lanna_order,][,-1], rep(NA,newyears)),
                        'error_SOC_multiplier_Lanna' = c(1, rep(1,8))[lanna_which],
                        'Yields_cereals_Lanna'= append_mean(aggregate(Lanna_yields_timeseries_long_cereals*10^-3, by=list(Lanna_treat), FUN=mean)[lanna_order,][,-1], N=newyears)[lanna_which,],
                        'Straw_cereals_Lanna'= append_mean(aggregate(Lanna_yields_timeseries_long_cereals_straw*10^-3, by=list(Lanna_treat), FUN=mean)[lanna_order,][,-1], N=newyears)[lanna_which,],
                        'I_FYM_Lanna'= append_mean(aggregate(h_FYM_Lanna_long*10^-3, by=list(Lanna_treat), FUN=mean)[lanna_order,][,-1], N=newyears)[lanna_which,],
                        'I_CO_Lanna'= append_mean(aggregate(h_CO_Lanna_long*10^-3, by=list(Lanna_treat), FUN=mean)[lanna_order,][,-1], N=newyears)[lanna_which,],
                        'I_SLU_Lanna'= append_mean(aggregate(h_SLU_Lanna_long*10^-3, by=list(Lanna_treat), FUN=mean)[lanna_order,][,-1], N=newyears)[lanna_which,],
                        'I_STR_Lanna'= append_mean(aggregate(h_STR_Lanna_long*10^-3, by=list(Lanna_treat), FUN=mean)[lanna_order,][,-1], N=newyears)[lanna_which,],
                        'N_Lanna' = dim(append_NA(aggregate(Lanna_SOC_timeseries_long, by=list(Lanna_treat), FUN=mean)[lanna_order,][,-1], N=newyears))[2],
                        're_Lanna'= append_mean((aggregate(re_Lanna_long, by=list(Lanna_treat), FUN=mean)[lanna_order,][,-1]), N=newyears)[lanna_which,],
                        'J_Lanna'= dim(append_NA(aggregate(Lanna_SOC_timeseries_long, by=list(Lanna_treat), FUN=mean)[lanna_order,][,-1], N=newyears)[lanna_which,])[1])




Lanna_date_vector_year_long<-seq.Date(from = as.Date(ISOdate(head(Lanna_years_first,1),1,1)),
                                       to=as.Date(ISOdate(head(Lanna_years_last+newyears,1),12,31)), by="year")
names(calib_data_Lanna[[1]])<-year(Lanna_date_vector_year_long)
names(calib_data_Lanna[[5]])<-year(Lanna_date_vector_year_long)
names(calib_data_Lanna[[6]])<-year(Lanna_date_vector_year_long)
names(calib_data_Lanna[[7]])<-year(Lanna_date_vector_year_long)
names(calib_data_Lanna[[8]])<-year(Lanna_date_vector_year_long)
names(calib_data_Lanna[[9]])<-year(Lanna_date_vector_year_long)
names(calib_data_Lanna[[10]])<-year(Lanna_date_vector_year_long)
names(calib_data_Lanna[[12]])<-year(Lanna_date_vector_year_long)

rownames(calib_data_Lanna[[1]])<-LETTERS[1:8]
names(calib_data_Lanna[[2]])<-LETTERS[1:8]
names(calib_data_Lanna[[3]])<-LETTERS[1:8]
names(calib_data_Lanna[[4]])<-LETTERS[1:8]
rownames(calib_data_Lanna[[5]])<-LETTERS[1:8]
rownames(calib_data_Lanna[[6]])<-LETTERS[1:8]
rownames(calib_data_Lanna[[7]])<-LETTERS[1:8]
rownames(calib_data_Lanna[[8]])<-LETTERS[1:8]
rownames(calib_data_Lanna[[9]])<-LETTERS[1:8]
rownames(calib_data_Lanna[[10]])<-LETTERS[1:8]
rownames(calib_data_Lanna[[12]])<-LETTERS[1:8]

data_descr_Lanna<-data.frame(names(calib_data_Lanna), 
                       c("SOC, 15 treatments x lenght of simulation years",
                       "SOC initialization values, 15 treatments",
                       "SOC mean errors, 15 treatments", 
                       "SOC error weight",
                       "Cereals yields, 8 treatments x lenght of simulation years",
                       "Cereals straw, 8 treatments x lenght of simulation years",
                       "Inputs farmyard manure, 15 treatments x lenght of simulation years",
                       "Inputs compost, 8 treatments x lenght of simulation years",
                       "Inputs sludges, 8 treatments x lenght of simulation years",
                       "Inputs straw, 8 treatments x lenght of simulation years",
                       "Lenght of the simulation",
                       "re_clim, 8 treatments x lenght of simulation years",
                       "number of treatments"))

length(calib_data_Lanna)

write_ods(as.data.frame(calib_data_Lanna[[1]]), "Input_Lanna.ods", sheet=names(calib_data_Lanna)[1], col_names = T, row_names = T)
write_ods(as.data.frame(calib_data_Lanna[[2]]), "Input_Lanna.ods", sheet=names(calib_data_Lanna)[2], append = T, col_names = T, row_names = T)
write_ods(as.data.frame(calib_data_Lanna[[3]]), "Input_Lanna.ods", sheet=names(calib_data_Lanna)[3], append = T, col_names = T, row_names = T)
write_ods(as.data.frame(calib_data_Lanna[[4]]), "Input_Lanna.ods", sheet=names(calib_data_Lanna)[4], append = T, col_names = T, row_names = T)
write_ods(as.data.frame(calib_data_Lanna[[5]]), "Input_Lanna.ods", sheet=names(calib_data_Lanna)[5], append = T, col_names = T, row_names = T)
write_ods(as.data.frame(calib_data_Lanna[[6]]), "Input_Lanna.ods", sheet=names(calib_data_Lanna)[6], append = T, col_names = T, row_names = T)
write_ods(as.data.frame(calib_data_Lanna[[7]]), "Input_Lanna.ods", sheet=names(calib_data_Lanna)[7], append = T, col_names = T, row_names = T)
write_ods(as.data.frame(calib_data_Lanna[[8]]), "Input_Lanna.ods", sheet=names(calib_data_Lanna)[8], append = T, col_names = T, row_names = T)
write_ods(as.data.frame(calib_data_Lanna[[9]]), "Input_Lanna.ods", sheet=names(calib_data_Lanna)[9], append = T, col_names = T, row_names = T)
write_ods(as.data.frame(calib_data_Lanna[[10]]), "Input_Lanna.ods", sheet=names(calib_data_Lanna)[10], append = T, col_names = T, row_names = T)
write_ods(as.data.frame(calib_data_Lanna[[11]]), "Input_Lanna.ods", sheet=names(calib_data_Lanna)[11], append = T, col_names = T, row_names = T)
write_ods(as.data.frame(calib_data_Lanna[[12]]), "Input_Lanna.ods", sheet=names(calib_data_Lanna)[12], append = T, col_names = T, row_names = T)
write_ods(as.data.frame(calib_data_Lanna[[13]]), "Input_Lanna.ods", sheet=names(calib_data_Lanna)[13], append = T, col_names = T, row_names = T)

colnames(data_descr_Lanna)<-c("Sheet name", "Description")













######################################### Checking the data


png("SOC_Yields_relationship.png", height = 3500, width = 2000, res=300)
par(mfrow=c(2,1))
Lanna_names<- unique(Lanna_treat)
Lanna_names[7]<-"X"
plot(rowMeans(calib_data_Lanna[["SOC_Lanna"]],  na.rm = T), rowMeans(calib_data_Lanna[["Yields_cereals_Lanna"]]),  col=Lanna_treat_palette[lanna_order], pch=16, ylab="Yields (cereals), total (Mg ha-1)", xlab="SOC (Mg ha-1)", main="Lanna", ylim=c(-0.5,9))
text(rowMeans(calib_data_Lanna[["SOC_Lanna"]],  na.rm = T), rowMeans(calib_data_Lanna[["Yields_cereals_Lanna"]])-.35, Lanna_names,  col=Lanna_treat_palette[lanna_order])
plot(rowMeans(calib_data_Ultuna[["SOC_Ultuna"]],  na.rm = T), rowMeans(calib_data_Ultuna[["Yields_cereals_Ultuna"]]),  col=Ultuna_treat_palette[seq(1:15)], pch=16, ylab="Yields (cereals), total (Mg ha-1)", xlab="SOC (Mg ha-1)", main="Ultuna", ylim=c(-0.5,5))
text(rowMeans(calib_data_Ultuna[["SOC_Ultuna"]],  na.rm = T), rowMeans(calib_data_Ultuna[["Yields_cereals_Ultuna"]])-.35, LETTERS[seq(1:15)],  col=Ultuna_treat_palette[seq(1:15)])
dev.off()


png("SOC_procent_Yields_relationship.png", height = 3500, width = 2000, res=300)
par(mfrow=c(2,1))
Lanna_treat2<-Lanna_treat
Lanna_treat2[25:28]<-"O2"
Lanna_treat2[29:32]<-"X"
Lanna_ave_grains<-aggregate(yields_Lanna$grain_dm_kg_ha, by=list(yields_Lanna$treat), FUN=mean)
Lanna_ave_SOC<-aggregate(yields_Lanna$C_conc, by=list(yields_Lanna$treat), FUN=mean,na.rm = TRUE)
plot(Lanna_ave_SOC[,2], Lanna_ave_grains[,2]*10^-3,  col=Lanna_treat_palette[Lanna_ave_SOC[,1]], pch=16, ylab="Yields (cereals), grains (Mg ha-1)", xlab="SOC (%)", main="Lanna", ylim=c(-0.5,5))
text(Lanna_ave_SOC[,2], Lanna_ave_grains[,2]*10^-3-.35, unique(Lanna_treat2),  col=Lanna_treat_palette[Lanna_ave_SOC[,1]])
Ultuna_ave_grains<-aggregate(yields_Ultuna$grain_dm_kg_ha, by=list(yields_Ultuna$treat), FUN=mean)
Ultuna_ave_SOC<-aggregate(yields_Ultuna$C_conc, by=list(yields_Ultuna$treat), FUN=mean,na.rm = TRUE)
plot(Ultuna_ave_SOC[,2], Ultuna_ave_grains[,2]*10^-3,  col=Ultuna_treat_palette[Ultuna_ave_SOC[,1]], pch=16, ylab="Yields (cereals), grains (Mg ha-1)", xlab="SOC (%)", main="Ultuna", ylim=c(-0.5,5))
text(Ultuna_ave_SOC[,2], Ultuna_ave_grains[,2]*10^-3-.35, LETTERS[seq(1:15)],  col=Ultuna_treat_palette[Ultuna_ave_SOC[,1]])
dev.off()





save.image(file="Pretreatment.Rdata")
