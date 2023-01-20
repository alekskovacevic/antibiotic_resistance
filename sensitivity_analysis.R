library(tidyverse)
library(wesanderson)

### NOTE: Must load this filepath for accompanying figure file 
filepath_project_group = "~/Desktop/COVID_AMR_GitHub/"

source(paste0(filepath_project_group, "S0_no_cov.R"))
source(paste0(filepath_project_group, "S0_Snc.R"))
source(paste0(filepath_project_group, "S1_SAi.R"))
source(paste0(filepath_project_group, "S2_Sad.R"))
source(paste0(filepath_project_group, "S3_SL.R"))
source(paste0(filepath_project_group, "S4_SAiL.R"))
source(paste0(filepath_project_group, "S5_SadL.R"))

vec_R_init = seq(0,1,0.1)
vec_R0 = seq(0,10,1)

outputSA = data.frame()
rate_IPD = parameters["p"]

for(R_init_i in vec_R_init){
  for(R0_i in vec_R0){
    print(paste0("R_init is ", R_init_i," and R0_i is ",R0_i))
    
    # translate R0 into betaC
    
    betaC_i = R0_i*parameters['gammaC']
    
    # update recovered population at simulation outset
    initial_state_values_i = initial_state_values
    initial_state_values_i['S'] = initial_state_values['S']*(1-R_init_i)
    initial_state_values_i['Ss'] = initial_state_values['Ss']*(1-R_init_i)
    initial_state_values_i['Sr'] = initial_state_values['Sr']*(1-R_init_i)
    initial_state_values_i['Ssr'] = initial_state_values['Ssr']*(1-R_init_i)
    initial_state_values_i['R'] = initial_state_values['S']*(R_init_i)
    initial_state_values_i['Rs'] = initial_state_values['Ss']*(R_init_i)
    initial_state_values_i['Rr'] = initial_state_values['Sr']*(R_init_i)
    initial_state_values_i['Rsr'] = initial_state_values['Ssr']*(R_init_i)
    
    # update transmission rate
    parameters_i = parameters
    parameters_i['betaC'] = betaC_i
    
    outputSA0_i = as.data.frame(ode(y = initial_state_values_i, 
                                    times = times_model, 
                                    func = model_0,
                                    parms = parameters_i,
                                    method = "euler"))%>%
      mutate(IPD_S = (Ss + Ssr*0.5 + Sas + Sasr*0.5 + Es + Esr*0.5 + Eas + Easr*0.5 + Is + Isr*0.5 + Ias + Iasr*0.5 + Rs + Rsr*0.5 + Ras + Rasr*0.5)*rate_IPD,
             IPD_R = (Sr + Ssr*0.5 + Sar + Sasr*0.5 + Er + Esr*0.5 + Ear + Easr*0.5 + Ir + Isr*0.5 + Iar + Iasr*0.5 + Rr + Rsr*0.5 + Rar + Rasr*0.5)*rate_IPD,
             IPD = IPD_S + IPD_R,
             R_init = R_init_i,
             R0 = R0_i,
             model = "model_0")%>%
      dplyr::select(-c(atb_0, beta_0))
    
    time_start = 4000#4120#4000
    time_end_12mo = 4365#4485#4365
    
    outputSA_i = rbind(outputSA0_i)%>%
      dplyr::select(time, R_init, R0, model, IPD_S, IPD_R, IPD)%>%
      pivot_longer(-c(time, R_init, R0, model),
                   names_to = "strain",
                   values_to = "incidence")%>%
      filter(time > time_start,
             time <= time_end_12mo)%>%
      group_by(R_init, R0, model, strain)%>%
      summarise(IPD_total_12mo = sum(incidence))%>%
      as.data.frame()
    
    outputSA = rbind(outputSA, outputSA_i)
  }
}

save(outputSA, file = paste0(filepath_project_group, "outputSA_2120_2485.Rdata"))
