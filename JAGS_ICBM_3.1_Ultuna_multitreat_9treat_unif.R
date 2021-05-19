model{
  
  ####Loop for treatments
  for(j in 1:J_Ultuna)
  {
    
    Y_R_Ultuna[j,1]<-(SOC_init_Ultuna[j]*(1-Init_ratio_Ultuna))*0.5
    Y_S_Ultuna[j,1]<-(SOC_init_Ultuna[j]*(1-Init_ratio_Ultuna))*0.5
    Y_FYM_Ultuna[j,1]<-0
    Y_GM_Ultuna[j,1]<-0
    Y_PEA_Ultuna[j,1]<-0
    Y_SAW_Ultuna[j,1]<-0
    Y_SLU_Ultuna[j,1]<-0
    Y_STR_Ultuna[j,1]<-0
    O_Ultuna[j,1]  <-SOC_init_Ultuna[j]*Init_ratio_Ultuna
 
    # loop for years
    for (i in 1:(N_Ultuna)){
      
      alpha_1[j,i]=ifelse(j==1, 0, alpha) #ifelse to have the intercept zero in case of bare fallow
      alpha_maize_1[j,i]=ifelse(j==1, 0, alpha_maize) #ifelse to have the intercept zero in case of bare fallow
      

      #Inputs R (roots), with different allometric functions for crops
      #depth coefficients from: Fan, J., McConkey, B., Wang, H., & Janzen, H. (2016). Root distribution by depth for 
      #temperate agricultural crops. 
      #Field Crops Research, 189, 68–74. https://doi.org/10.1016/j.fcr.2016.02.013
      I_R_cereals_Ultuna[j,i]      <-  (1+exudates_coeff)*e_depth_cer*0.5560437*C_percent*
                                                  ((alpha_1[j,i]+Yields_cereals_Ultuna[j,i])*(1/(SR_cereals_ult)))
      I_R_root_crops_Ultuna[j,i]   <-  (1+exudates_coeff)*e_depth_root*0.5902446*0.32*C_percent*
                                                  ((alpha_1[j,i]+Yields_root_crops_Ultuna[j,i])*(1/(SR_root_crops_ult)))#SR_root_crops_ult))
      I_R_oilseeds_Ultuna[j,i]     <-  (1+exudates_coeff)*e_depth_oil*0.6367459*C_percent*
                                                  ((alpha_1[j,i]+Yields_oilseeds_Ultuna[j,i])*(1/(SR_oilseeds_ult)))#SR_oilseeds_ult))
      I_R_maize_Ultuna[j,i]        <-  (1+exudates_coeff)*e_depth_maize*0.5981638*C_percent*
                                                  ((alpha_maize_1[j,i]+Yields_maize_Ultuna[j,i])*(1/(SR_maize_ult)))#SR_maize_ult))
      I_R_Ultuna[j,i]              <-  I_R_cereals_Ultuna[j,i]+I_R_root_crops_Ultuna[j,i]+I_R_oilseeds_Ultuna[j,i]+I_R_maize_Ultuna[j,i]
      #Inputs S
      I_S_Ultuna[j,i]              <- ((Yields_cereals_Ultuna[j,i]+Yields_root_crops_Ultuna[j,i]+Yields_oilseeds_Ultuna[j,i])*
                                                  stubbles_ratio_Ultuna+Yields_maize_Ultuna[j,i]*stubbles_ratio_Ultuna_maize)*C_percent
      
      #Young R
      Y_R_Ultuna[j,i+1] 		<-  (I_R_Ultuna[j,i]+Y_R_Ultuna[j,i])*exp(-k1_ult*re_Ultuna[j,i])
      Y_S_Ultuna[j,i+1] 		<-  (I_S_Ultuna[j,i]+Y_S_Ultuna[j,i])*exp(-k1_ult*re_Ultuna[j,i])
      
      Y_FYM_Ultuna[j,i+1] 		<-  (I_FYM_Ultuna[j,i]+Y_FYM_Ultuna[j,i])*exp(-k1_ult*re_Ultuna[j,i])
      Y_GM_Ultuna[j,i+1] 		  <-  (I_GM_Ultuna[j,i]+Y_GM_Ultuna[j,i])*exp(-k1_ult*re_Ultuna[j,i])
      Y_PEA_Ultuna[j,i+1] 		<-  (I_PEA_Ultuna[j,i]+Y_PEA_Ultuna[j,i])*exp(-k1_ult*re_Ultuna[j,i])
      Y_SAW_Ultuna[j,i+1] 		<-  (I_SAW_Ultuna[j,i]+Y_SAW_Ultuna[j,i])*exp(-k1_ult*re_Ultuna[j,i])
      Y_SLU_Ultuna[j,i+1] 		<-  (I_SLU_Ultuna[j,i]+Y_SLU_Ultuna[j,i])*exp(-k1_ult*re_Ultuna[j,i])
      Y_STR_Ultuna[j,i+1] 		<-  (I_STR_Ultuna[j,i]+Y_STR_Ultuna[j,i])*exp(-k1_ult*re_Ultuna[j,i])
      
      #Old
      
      #old flux
      fluxR_Ultuna[j,i] 		  <-  h_R_ult*((k1_ult*(Y_R_Ultuna[j,i]+I_R_Ultuna[j,i]))/(k2_ult-k1_ult))
      fluxS_Ultuna[j,i] 		  <-  h_S_ult*((k1_ult*(Y_S_Ultuna[j,i]+I_S_Ultuna[j,i]))/(k2_ult-k1_ult))

      flux_FYM_Ultuna[j,i] 		  <-  h_FYM_ult *((k1_ult*(Y_FYM_Ultuna[j,i]+I_FYM_Ultuna[j,i]))/(k2_ult-k1_ult))
      flux_GM_Ultuna[j,i] 		  <-  h_S_ult   *((k1_ult*(Y_GM_Ultuna[j,i] +I_GM_Ultuna[j,i])) /(k2_ult-k1_ult))
      flux_PEA_Ultuna[j,i] 		  <-  h_PEA_ult *((k1_ult*(Y_PEA_Ultuna[j,i]+I_PEA_Ultuna[j,i]))/(k2_ult-k1_ult))
      flux_SAW_Ultuna[j,i] 		  <-  h_SAW_ult *((k1_ult*(Y_SAW_Ultuna[j,i]+I_SAW_Ultuna[j,i]))/(k2_ult-k1_ult))
      flux_SLU_Ultuna[j,i] 		  <-  h_SLU_ult *((k1_ult*(Y_SLU_Ultuna[j,i]+I_SLU_Ultuna[j,i]))/(k2_ult-k1_ult))
      flux_STR_Ultuna[j,i] 		  <-  h_S_ult   *((k1_ult*(Y_STR_Ultuna[j,i]+I_STR_Ultuna[j,i]))/(k2_ult-k1_ult))
      
      flux_sum_Ultuna[j,i]<-(fluxR_Ultuna[j,i]+
                               fluxS_Ultuna[j,i]+
                               flux_FYM_Ultuna[j,i]+
                               flux_GM_Ultuna[j,i]+
                               flux_PEA_Ultuna[j,i]+
                               flux_SAW_Ultuna[j,i]+
                               flux_SLU_Ultuna[j,i]+
                               flux_STR_Ultuna[j,i])
      
      
      O_Ultuna[j,i+1]   	<-  (O_Ultuna[j,i]-flux_sum_Ultuna[j,i])*exp(-k2_ult*re_Ultuna[j,i]) + 
                                    flux_sum_Ultuna[j,i]*exp(-k1_ult*re_Ultuna[j,i])
      
      
      #Total C
      Y_tot_Ultuna[j,i] <- Y_R_Ultuna[j,i] +
                            Y_S_Ultuna[j,i] +
                            Y_FYM_Ultuna[j,i] +
                            Y_GM_Ultuna[j,i] +
                            Y_PEA_Ultuna[j,i] +
                            Y_SAW_Ultuna[j,i] +
                            Y_SLU_Ultuna[j,i] +
                            Y_STR_Ultuna[j,i]
      
      Tot_Ultuna[j,i] <- Y_tot_Ultuna[j,i] + O_Ultuna[j,i]
      
      #Error of the measurement (assumed proportional to the measurement)
      SOC_Ultuna[j,i]  ~ dnorm(Tot_Ultuna[j,i],1/(error_SOC_Ultuna[j]*error_SOC[j,i]))
      #SOC_Ultuna[j,i]  ~ dnorm(Tot_Ultuna[j,i],1/(error_SOC_Ultuna[j]*1.5))
      #SOC_Ultuna[j,i]  ~ dunif(Tot_Ultuna[j,i]-Tot_Ultuna[j,i]*error_SOC_Ultuna[j],Tot_Ultuna[j,i]+Tot_Ultuna[j,i]*error_SOC_Ultuna[j])
      
      error_SOC[j,i]~ dunif(0.75,1.25)
      
    }
    
  }
  
  
  
  ##Parameters Ultuna
  # Xie, Yajun. 2020. “A Meta-Analysis of Critique of Litterbag Method Used in Examining Decomposition of Leaf Litters.
  #” Journal of Soils and Sediments 20 (4): 1881–86. https://doi.org/10.1007/s11368-020-02572-9.
  
  k1_ult    ~ dunif(0.01, 1)
  k2_ult    ~ dunif(0,0.03)
  
   
  h_S_ult     ~ dunif(0.125-0.125*limits_h,0.125+0.125*limits_h)
  h_R_ult     ~ dunif(0.35-0.35*limits_h, 0.35+0.35*limits_h)
  h_FYM_ult   ~ dunif(0.27-0.27*limits_h,0.27+0.27*limits_h)
  h_PEA_ult   ~ dunif(0.59-0.59*limits_h, 0.59+0.59*limits_h)
  h_SAW_ult   ~ dunif(0.25-0.25*limits_h,0.25+0.25*limits_h)
  h_SLU_ult   ~ dunif(0.41-0.41*limits_h,0.41+0.41*limits_h)
  
  #root/shoot ratios priors
  SR_cereals_ult    ~ dunif(3.6,27.9) #range from martin data
  #SR_cereals_ult    ~ dnorm(11,1/2)
  SR_root_crops_ult ~ dunif(29.49853-29.49853*limit_SR,29.49853+29.49853*limit_SR)
  SR_oilseeds_ult   ~ dunif(8-8*limit_SR,8+8*limit_SR)
  SR_maize_ult      ~ dunif(6.25-6.25*limit_SR,6.25+6.25*limit_SR)
  
  #alpha ~dunif(0,1)
  #alpha_maize ~dunif(0,1)
  alpha ~dunif(0,0.001)
  alpha_maize ~dunif(0,0.001)
  
  
  exudates_coeff ~ dnorm(1.65, 1/0.4125) #10% error
  e_depth_cer ~ dnorm(1, 1/0.25) #10% error
  e_depth_root ~ dnorm(1, 1/0.25) #10% error
  e_depth_oil ~ dnorm(1, 1/0.25) #10% error
  e_depth_maize ~ dnorm(1, 1/0.25) #10% error
  
  #Init_ratio_Ultuna ~ dnorm(0.9291667,1/(0.9291667*0.2)) T(0.8,0.98)
  Init_ratio_Ultuna ~ dunif(0.8,0.98)
  
  stubbles_ratio_Ultuna_maize ~ dunif(0.01,0.08)
  stubbles_ratio_Ultuna ~ dunif(0.01,0.08)
  
  C_percent ~ dunif(0.40, 0.51)
  
  limits_h<-1
  limits_k<-1
  limit_SR<-1
  
}
