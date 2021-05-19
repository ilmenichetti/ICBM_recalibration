ICBM_single<-function(BF,
                      k1,
                      k2,
                      h_r,
                      h_s,
                      h_fym,
                      h_pea,
                      h_saw,
                      h_slu,
                      stubbles_others,
                      stubbles_maize,
                      S_R_cereals,
                      S_R_roots,
                      S_R_oilseeds,
                      S_R_maize,
                      OY_ratio,
                      alpha_others,
                      alpha_maize,
                      C_percent,   
                      Yields_cereals,
                      Yields_root_crops,
                      Yields_oilseeds,
                      Yields_maize,
                      I_fym,
                      I_gm,
                      I_pea,
                      I_saw,
                      I_slu,
                      I_str, 
                      SOC_init, 
                      re){

exudates_coeff <- 1.65
  
#Inputs R (roots), with differe[i]nt allometric functions for crops
if(BF==T){
alpha_others=0
alpha_maize=0
}

I_R_cereals      <-  (1+exudates_coeff)*0.7*C_percent*((alpha_others+Yields_cereals)*(1/(S_R_cereals)))
I_R_root_crops   <-  (1+exudates_coeff)*0.7*0.32*C_percent*((alpha_others+Yields_root_crops)*(1/(S_R_roots)))#SR_root_crops))
I_R_oilseeds     <-  (1+exudates_coeff)*0.7*C_percent*((alpha_others+Yields_oilseeds)*(1/(S_R_oilseeds)))#SR_oilseeds))
I_R_maize        <-  (1+exudates_coeff)*0.7*C_percent*((alpha_maize+Yields_maize)*(1/(S_R_maize)))#SR_maize))
I_R              <-  I_R_cereals+I_R_root_crops+I_R_oilseeds+I_R_maize

#Inputs S (shoots)
I_S              <- ((Yields_cereals+Yields_root_crops+Yields_oilseeds)*stubbles_others+Yields_maize*stubbles_maize)*C_percent


Y_R<-c()
Y_S<-c()
Y_FYM<-c()
Y_GM<-c()
Y_PEA<-c()
Y_SAW<-c()
Y_SLU<-c()
Y_STR<-c()
Y_tot<-c()
O<-c()
Tot<-c()

Y_R[1]<-(SOC_init*(1-OY_ratio))*0.5
Y_S[1]<-(SOC_init*(1-OY_ratio))*0.5
Y_FYM[1]<-0
Y_GM[1]<-0
Y_PEA[1]<-0
Y_SAW[1]<-0
Y_SLU[1]<-0
Y_STR[1]<-0
O[1]<-(SOC_init*OY_ratio)

#loop running for each yer for Y and O
for(i in 1:length(I_R)){

###Young pools
Y_R[i+1] 		<-  (I_R[i]+Y_R[i])*exp(-k1*re[i])
Y_S[i+1] 		<-  (I_S[i]+Y_S[i])*exp(-k1*re[i])

Y_FYM[i+1] 		<-  (I_fym[i]+Y_FYM[i])*exp(-k1*re[i])
Y_GM[i+1] 		<-  (I_gm[i]+Y_GM[i])*exp(-k1*re[i])
Y_PEA[i+1] 		<-  (I_pea[i]+Y_PEA[i])*exp(-k1*re[i])
Y_SAW[i+1] 		<-  (I_saw[i]+Y_SAW[i])*exp(-k1*re[i])
Y_SLU[i+1] 		<-  (I_slu[i]+Y_SLU[i])*exp(-k1*re[i])
Y_STR[i+1] 		<-  (I_str[i]+Y_STR[i])*exp(-k1*re[i])

###Old pool

#fluxes
fluxR    		  <-  h_r*((k1*(Y_R[i]+I_R[i]))/(k2-k1))
fluxS 	  	  <-  h_s*((k1*(Y_S[i]+I_S[i]))/(k2-k1))
flux_FYM 		  <-  h_fym*((k1*(Y_FYM[i]+I_fym[i]))/(k2-k1))
flux_GM 		  <-  h_s*((k1*(Y_GM[i]+I_gm[i]))/(k2-k1))
flux_PEA 		  <-  h_pea*((k1*(Y_PEA[i]+I_pea[i]))/(k2-k1))
flux_SAW 		  <-  h_saw*((k1*(Y_SAW[i]+I_saw[i]))/(k2-k1))
flux_SLU 		  <-  h_slu*((k1*(Y_SLU[i]+I_slu[i]))/(k2-k1))
flux_STR 		  <-  h_s*((k1*(Y_STR[i]+I_str[i]))/(k2-k1))

flux_sum<-(fluxR+
           fluxS+
           flux_FYM+
           flux_GM+
           flux_PEA+
           flux_SAW+
           flux_SLU+
           flux_STR)
#O 
O[i+1]   	<-  (O[i]-flux_sum)*exp(-k2*re[i]) + flux_sum*exp(-k1*re[i])


#Total C
Y_tot[i] <- Y_R[i] +
            Y_S[i] +
            Y_FYM[i] +
            Y_GM[i] +
            Y_PEA[i] +
            Y_SAW[i] +
            Y_SLU[i] +
            Y_STR[i]

Tot[i] <- Y_tot[i] + O[i]

}

return(Tot)
}
