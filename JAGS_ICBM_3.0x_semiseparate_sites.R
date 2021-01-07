model{


  ####Loop for Ultuna
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


    for (i in 1:(N_Ultuna)){

      #Inputs R (roots), with different allometric functions for crops
      I_R_cereals_Ultuna[j,i]      <-  (1+exudates_coeff)*((Yields_cereals_Ultuna[j,i])*0.7*C_percent*(1/SR_cereals_ult))
      I_R_root_crops_Ultuna[j,i]   <-  (1+exudates_coeff)*((Yields_root_crops_Ultuna[j,i])*0.7*0.32*C_percent*(1/SR_root_crops_ult))
      I_R_oilseeds_Ultuna[j,i]     <-  (1+exudates_coeff)*((Yields_oilseeds_Ultuna[j,i])*0.7*C_percent*(1/SR_oilseeds_ult))
      I_R_maize_Ultuna[j,i]        <-  (1+exudates_coeff)*((Yields_maize_Ultuna[j,i])*0.7*C_percent*(1/SR_maize_ult))
      I_R_Ultuna[j,i]              <-  I_R_cereals_Ultuna[j,i]+I_R_root_crops_Ultuna[j,i]+I_R_oilseeds_Ultuna[j,i]+I_R_maize_Ultuna[j,i]
      #Inputs S
      I_S_Ultuna[j,i]              <- (Yields_cereals_Ultuna[j,i]+Yields_root_crops_Ultuna[j,i]+Yields_oilseeds_Ultuna[j,i]+Yields_maize_Ultuna[j,i])*stubbles_ratio_Ultuna*C_percent

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
      fluxR_Ultuna[j,i] 		  <-  h_R_ult*((k1_ult*(Y_R_Ultuna[j,i]+I_R_Ultuna[j,i]))/(k2_ult-k1_ult))
      fluxS_Ultuna[j,i] 		  <-  h_S_ult*((k1_ult*(Y_S_Ultuna[j,i]+I_S_Ultuna[j,i]))/(k2_ult-k1_ult))

      #old flux manure
      flux_FYM_Ultuna[j,i] 		  <-  h_FYM_ult*((k1_ult*(Y_FYM_Ultuna[j,i]+I_FYM_Ultuna[j,i]))/(k2_ult-k1_ult))
      flux_GM_Ultuna[j,i] 		  <-  h_S_ult*((k1_ult*(Y_GM_Ultuna[j,i]+I_GM_Ultuna[j,i]))/(k2_ult-k1_ult))
      flux_PEA_Ultuna[j,i] 		  <-  h_PEA_ult*((k1_ult*(Y_PEA_Ultuna[j,i]+I_PEA_Ultuna[j,i]))/(k2_ult-k1_ult))
      flux_SAW_Ultuna[j,i] 		  <-  h_SAW_ult*((k1_ult*(Y_SAW_Ultuna[j,i]+I_SAW_Ultuna[j,i]))/(k2_ult-k1_ult))
      flux_SLU_Ultuna[j,i] 		  <-  h_SLU_ult*((k1_ult*(Y_SLU_Ultuna[j,i]+I_SLU_Ultuna[j,i]))/(k2_ult-k1_ult))
      flux_STR_Ultuna[j,i] 		  <-  h_S_ult*((k1_ult*(Y_STR_Ultuna[j,i]+I_STR_Ultuna[j,i]))/(k2_ult-k1_ult))

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
      SOC_Ultuna[j,i]  ~ dnorm(Tot_Ultuna[j,i],1/(error_SOC_Ultuna[j]*error_SOC_multiplier_Ultuna[j]))


    }

  }


  ####Loop for Lanna
  for(k in 1:J_Lanna)
  {

    Y_R_Lanna[k,1]<-(SOC_init_Lanna[k]*(1-Init_ratio_Lanna))*0.5
    Y_S_Lanna[k,1]<-(SOC_init_Lanna[k]*(1-Init_ratio_Lanna))*0.5
    Y_FYM_Lanna[k,1]<-0
    Y_GM_Lanna[k,1]<-0
    Y_CO_Lanna[k,1]<-0
    Y_SLU_Lanna[k,1]<-0
    Y_STR_Lanna[k,1]<-0
    O_Lanna[k,1]  <-SOC_init_Lanna[k]*Init_ratio_Lanna

    for (w in 1:(N_Lanna-1)){

      #Inputs R (roots), with different allometric functions for crops
      I_R_cereals_Lanna[k,w]      <-  (1+exudates_coeff)*((Yields_cereals_Lanna[k,w])*0.7*C_percent*(1/SR_cereals_lan))
      I_R_Lanna[k,w]              <-  I_R_cereals_Lanna[k,w]
      #Inputs S
      I_S_Lanna[k,w]              <-   (Straw_cereals_Lanna[k,w])*stubbles_ratio_Lanna*C_percent #in Lanna the stubbles ratio (33%) refers to the straw

      #Young R
      Y_R_Lanna[k,w+1] 		<-  (I_R_Lanna[k,w]+Y_R_Lanna[k,w])*exp(-k1_ult*re_Lanna[k,w])
      Y_S_Lanna[k,w+1] 		<-  (I_S_Lanna[k,w]+Y_S_Lanna[k,w])*exp(-k1_ult*re_Lanna[k,w])

      Y_FYM_Lanna[k,w+1] 		<-  (I_FYM_Lanna[k,w]+Y_FYM_Lanna[k,w])*exp(-k1_ult*re_Lanna[k,w])
      Y_CO_Lanna[k,w+1] 		<-  (I_CO_Lanna[k,w]+Y_CO_Lanna[k,w])*exp(-k1_ult*re_Lanna[k,w])
      Y_SLU_Lanna[k,w+1] 		<-  (I_SLU_Lanna[k,w]+Y_SLU_Lanna[k,w])*exp(-k1_ult*re_Lanna[k,w])
      Y_STR_Lanna[k,w+1] 		<-  (I_STR_Lanna[k,w]+Y_STR_Lanna[k,w])*exp(-k1_ult*re_Lanna[k,w])

      #Old
      fluxR_Lanna[k,w] 		  <-  h_R_lan*((k1_ult*(Y_R_Lanna[k,w]+I_R_Lanna[k,w]))/(k2_ult-k1_ult))
      fluxS_Lanna[k,w] 		  <-  h_S_lan*((k1_ult*(Y_S_Lanna[k,w]+I_S_Lanna[k,w]))/(k2_ult-k1_ult))

      #old flux manure
      flux_FYM_Lanna[k,w] 		  <-  h_FYM_lan*((k1_ult*(Y_FYM_Lanna[k,w]+I_FYM_Lanna[k,w]))/(k2_ult-k1_ult))
      flux_CO_Lanna[k,w] 		    <-  h_CO_lan*((k1_ult*(Y_CO_Lanna[k,w]+I_CO_Lanna[k,w]))/(k2_ult-k1_ult))
      flux_SLU_Lanna[k,w] 		  <-  h_SLU_lan*((k1_ult*(Y_SLU_Lanna[k,w]+I_SLU_Lanna[k,w]))/(k2_ult-k1_ult))
      flux_STR_Lanna[k,w] 		  <-  h_S_lan*((k1_ult*(Y_STR_Lanna[k,w]+I_STR_Lanna[k,w]))/(k2_ult-k1_ult))

      flux_sum_Lanna[k,w]<-(fluxR_Lanna[k,w]+
                   fluxS_Lanna[k,w]+
                   flux_FYM_Lanna[k,w]+
                   flux_CO_Lanna[k,w]+
                   flux_SLU_Lanna[k,w]+
                   flux_STR_Lanna[k,w])

      O_Lanna[k, w+1]   	<-  (O_Lanna[k,w]-flux_sum_Lanna[k,w])*exp(-k2_ult*re_Lanna[k,w]) +
                                       flux_sum_Lanna[k,w]*exp(-k1_ult*re_Lanna[k,w])


      #Total C
      Y_tot_Lanna[k,w] <- Y_R_Lanna[k,w] +
                            Y_S_Lanna[k,w] +
                            Y_FYM_Lanna[k,w] +
                            Y_CO_Lanna[k,w] +
                            Y_SLU_Lanna[k,w] +
                            Y_STR_Lanna[k,w]

      Tot_Lanna[k,w] <- Y_tot_Lanna[k,w]+ O_Lanna[k,w]

      #Error of the measurement (assumed proportional to the measurement)
      SOC_Lanna[k,w]  ~ dnorm(Tot_Lanna[k,w],1/(error_SOC_Lanna[k]*error_SOC_multiplier_Lanna[k]))

    }

  }



  ##Parameters Ultuna
  # Xie, Yajun. 2020. “A Meta-Analysis of Critique of Litterbag Method Used in Examining Decomposition of Leaf Litters.” Journal of Soils and Sediments 20 (4): 1881–86. https://doi.org/10.1007/s11368-020-02572-9.
  k1_ult    ~ dunif(0.78, 1)
  k2_ult    ~ dnorm(0.00605, 1/(0.00605*error_h)) T(0.00605-0.00605*0.5,0.00605+0.00605*0.5)

  h_S_ult     ~ dnorm(0.15,1/(0.15*error_h)) T(0.15-limits_h,0.15+limits_h)
  h_R_ult     ~  dnorm(0.35,1/(0.35*error_h)) T(0.35-limits_h, 1)
  h_FYM_ult   ~ dnorm(0.27,1/(0.27*error_h)) T(0.27-limits_h,0.27+limits_h)
  h_PEA_ult   ~ dnorm(0.59,1/(0.59*error_h)) T(0.59-limits_h, 0.59+limits_h)
  h_SAW_ult   ~ dnorm(0.25,1/(0.25*error_h)) T(0.25-limits_h,0.25+limits_h)
  h_SLU_ult   ~ dnorm(0.41,1/(0.41*error_h)) T(0.41-limits_h,0.41+limits_h)

  #root/shoot ratios priors
  SR_cereals_ult    ~ dnorm(11, 1/(11*error_SR)) T(11-11*error_SR,11+11*error_SR)
  SR_cereals_lan    ~ dnorm(11, 1/(11*error_SR)) T(11-11*error_SR,11+11*error_SR)
  SR_root_crops_ult ~ dnorm(29.49853, 1/(29.49853*error_SR)) T(29.49853-29.49853*error_SR,29.49853+29.49853*error_SR)
  SR_oilseeds_ult   ~ dnorm(8, 1/(8*error_SR)) T(8-8*error_SR,8+8*error_SR)
  SR_maize_ult      ~ dnorm(6.25, 1/(6.25*error_SR)) T(6.25-6.25*error_SR,6.25+6.25*error_SR)

  ##Parameters Lanna
  h_S_lan     ~ dnorm(0.15,1/(0.15*error_h)) T(0.15-limits_h,0.15+limits_h)
  h_R_lan     ~ dnorm(0.35,1/(0.35*error_h)) T(0.35-limits_h,1)
  h_FYM_lan   ~ dnorm(0.27,1/(0.27*error_h)) T(0.27-limits_h,0.27+limits_h)
  h_SLU_lan   ~ dnorm(0.41,1/(0.41*error_h)) T(0.41-limits_h,0.41+limits_h)
  h_CO_lan    ~ dnorm(0.80,1/(0.80*error_h)) T(0.15,0.95)



  exudates_coeff ~ dnorm(1.65,1/(1.65*0.1)) T(1.65*0.95,1.65*1.05)

  Init_ratio_Ultuna ~ dnorm(0.9291667,1/(0.9291667*0.2)) T(0.8,0.98)
  Init_ratio_Lanna ~ dnorm(0.9291667,1/(0.9291667*0.2)) T(0.8,0.98)

  stubbles_ratio_Ultuna ~ dnorm(0.04,1/0.01) T(0.02,0.06)
  stubbles_ratio_Lanna ~ dnorm(0.33,1/0.05) T(0.1,0.45)

  C_percent ~ dunif(0.40, 0.51)

  error_h<-0.2
  limits_h<-0.3
  error_SR<-0.25


}
