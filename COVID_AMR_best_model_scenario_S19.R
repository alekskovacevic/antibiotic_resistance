# Author: Aleksandra Kovacevic 12/07/2023
# Co-circulation model of SARS-CoV-2 transmission and infection, and pneumococcal colonization
# Table 1, scenario S19

# load the packages:
library(deSolve)
library(reshape2)
library(ggplot2)
library(scales)
library(readxl)
library(shiny)
library(xlsx)
#library(ggh4x)
library(dplyr)
library(cowplot)
library(ggpubr)
library(ggpattern)

#devtools::install_github("zeehio/facetscales")
#library(g)
library(facetscales)


# fixed parameters for pneumococcal carriage
q=1/2        # relative infectiousness with each strain for dually colonized
c=1/2        # fraction of dually colonized returning to single-colonized upon reinfection with S (R) 
k=1/2        # probability of acquiring secondary bacterial carriage
ps=1/2       # probability of transmitting S strain
psingle=1/2  # probability of a single infection

# Model inputs: 
# Total population size
N <- 100000

# pneumoccocal carriage assumed to be 20% where resistance represents 4/20 (20% AR)
initial_state <- c(S = N - 0.20*N, Ss = 0.16*N, Sr = 0.04*N, Ssr = 0, Sss = 0, Srr = 0,
                   Sa = 0, Ssa = 0, Sra = 0, Ssra = 0, Sssa = 0, Srra = 0,
                   E = 0, Es = 0, Er = 0, Esr = 0, Ess = 0, Err = 0,
                   Ea = 0, Esa = 0, Era = 0, Esra = 0, Essa = 0, Erra = 0,
                   I = 0, Is = 0, Ir = 0, Isr = 0, Iss = 0, Irr = 0,
                   Ia = 0, Isa = 0, Ira = 0, Isra = 0, Issa = 0, Irra = 0,
                   Iaz = 0, Isaz = 0, Iraz = 0, Israz = 0, Issaz = 0, Irraz = 0, # infected + prophylaxis
                   R = 0, Rs = 0, Rr = 0, Rsr = 0, Rss = 0, Rrr = 0,
                   Ra = 0, Rsa = 0, Rra = 0, Rsra = 0, Rssa = 0, Rrra = 0,
                   Raz = 0, Rsaz = 0, Rraz = 0, Rsraz = 0, Rssaz = 0, Rrraz = 0, # recovered + prophylaxis
                   IPDs = 0, IPDr = 0, IPD = 0, ATB = 0) 

parameters_10 <- c(betaS = 0.042,   # pneumococcal transmission rate
                  gammaS = 1/30,    # pneumococcal recovery rate (sensitive strain) 
                  gammaR = 1/30,    # pneumococcal recovery rate (resistant strain) - MECHANISM 5
                  betaC = 0.46,     # covid transmission rate
                  alpha = 1/5,      # covid incubation rate 5 days, latent period
                  gammaC = 1/7,     # covid recovery rate 7 days
                  fit= 0.949,       # fitness cost to be at pre-pandemic equilibrium
                  w = 1/3,          # rate of antibiotic action; takes 3 days for antibiotics to clear susceptible carriage
                  tau = 0.0014,     # antibiotic prescription rate (per capita), daily rate, calculated from Bara et al 2022 
                  r = 1/7,          # baseline antibiotic treatment duration
                  r_az = 1/11.5    # 11.5 (if 3-day treatment) or 13.5 if 5 day treatment
)

times <- seq(0, 3650, by = 1)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parameters to change according to the scenario:
p_az_10 = 0.10                     # MECHANISM 4 - proportion of infected individuals getting azithromycin treatment (0, 0.5, 0.1, 0.15, and 0.20)
psummer_10 = 1e-06                 # summer invasion rate
pwinter_10 = 2.5e-06               # winter invasion rate 
a_lockdown_10 = 0.77               # reduction in the nb of antibiotic prescriptions  
covid_transm_lockdown_10 = 0.23    # reduction of the covid transmission during lockdown
IPD_risk_lockdown_10 = 0.2         # risk of IPD during lockdown (0.2)
IPD_risk_rest_10 = 0.4             # risk of IPD during the rest of the year (0.4)
covid.present_10 = 1               # presence (1)/absence(0) of covid (pre-pandemic and pandemic period)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Interpolating functions of mechanisms implemented:

# MECHANISM 1:
# Reduced community antibiotic prescribing:
atb_10 <- as.data.frame(list(times = times, a_10 = rep(0, length(times))))   
atb_10$a_10 <- ifelse(atb_10$times > 1074 & atb_10$times < 1365, a_lockdown_10, 1)
input_10 <- approxfun(atb_10)
input_10(seq(from = 0, to = 3650, by = 1)) 

# MECHANISM 2:
# Lockdown effect on reducing transmission of S. pneumoniae:
beta_10p <- as.data.frame(list(times = times, thetab_10p = rep(0, length(times))))
beta_10p$thetab_10p <- ifelse((beta_10p$times > 1074 & beta_10p$times < 1136|beta_10p$times > 1304 & beta_10p$times < 1366), 1, 1) 
inputbeta_10p <- approxfun(beta_10p)
inputbeta_10p(seq(from = 0, to = 3650, by = 1))

# MECHANISM 3:
# Reduced IPD risk, IPD disease risk reduction during NPI implementation before, during, and after the lockdown (1, 0.2, and 0.4) 
risk_10 <- as.data.frame(list(times = times, IPDrisk_10 = rep(0, length(times))))
risk_10$IPDrisk_10 <- ifelse(risk_10$times < 1075, 1,
                             ifelse(risk_10$times > 1074 & risk_10$times < 1136, IPD_risk_lockdown_10, IPD_risk_rest_10)) 
inputrisk_10 <- approxfun(risk_10)
inputrisk_10(seq(from = 0, to = 3650, by = 1))

sim_time <- seq(from = 0, to = 3650, by = 1)  # model simulation time



# Pneumococcal invasion rate (summer vs winter invasion rate) (p):
p_10 <- as.data.frame(list(times = times, p_inv_10 = rep(0, length(times))))
p_10$p_inv_10 <- ifelse((p_10$times > 1089 & p_10$times < 1271), (psummer_10), (pwinter_10)) # winter invasion rate > summer invasion rate
invasion_10 <- approxfun(p_10)   #Create the interpolating function
invasion_10(seq(from = 0, to = 3650, by = 1))

# SARS-CoV-2 transmission risk:
beta_10 <- as.data.frame(list(times = times, thetabeta_10 = rep(0, length(times))))
beta_10$thetabeta_10 <- ifelse(beta_10$times < 1075, 1,
                               ifelse((beta_10$times > 1074 & beta_10$times < 1136), 
                                      covid_transm_lockdown_10, 
                                      ifelse(beta_10$times > 1135 & beta_10$times < 1231, 0.45, 0.85))) 
inputbeta_10 <- approxfun(beta_10)
inputbeta_10(seq(from = 0, to = 3650, by = 1))




################################### MODEL FUNCTION: ######################################

#MODEL FUNCTION: 

model_10 <- function(time, state, parameters) {  
  with(as.list(c(state, parameters)), { 
    covid = if(floor(time) == 1000) { covid = covid.present_10 } else {covid = 0}  # COVID
    
    inv_rate_10 <- invasion_10(time)         # invasion rate varies between summer and winter
    a_10 <- input_10(time)                   # antibiotic use in the community
    thetabeta_10 <- inputbeta_10(time)       # covid transmission risk 
    thetab_10 <- inputbeta_10p(time)         # carriage transmission risk
    IPDrisk_10 <- inputrisk_10(time)         # IPD risk (varies with the present or absence of ILIs)
    
    N <- S+Ss+Sr+Ssr+Sss+Srr+Sa+Ssa+Sra+Ssra+Sssa+Srra+
         E+Es+Er+Esr+Ess+Err+Ea+Esa+Era+Esra+Essa+Erra+ 
         I+Is+Ir+Isr+Iss+Irr+Ia+Isa+Ira+Isra+Issa+Irra+Iaz+Isaz+Iraz+Israz+Issaz+Irraz+
         R+Rs+Rr+Rsr+Rss+Rrr+Ra+Rsa+Rra+Rsra+Rssa+Rrra+Raz+Rsaz+Rraz+Rsraz+Rssaz+Rrraz
    
    #print(sum(N))
    
    # Forces of infection:
    lambdaS <- betaS*thetab_10*(Ss + 2*q*ps*Ssr + 2*q*Sss + Ssa + 2*q*ps*Ssra + 2*q*Sssa + 
                                Es + 2*q*ps*Esr + 2*q*Ess + Esa + 2*q*ps*Esra + 2*q*Essa + 
                                Is + 2*q*ps*Isr + 2*q*Iss + Isa + 2*q*ps*Isra + 2*q*Issa +
                                Isaz + 2*q*ps*Israz + 2*q*Issaz +
                                Rs + 2*q*ps*Rsr + 2*q*Rss + Rsa + 2*q*ps*Rsra + 2*q*Rssa +
                                Rsaz + 2*q*ps*Rsraz + 2*q*Rssaz)/N
    
    
    lambdaR <- betaS*fit*thetab_10*(Sr + 2*q*(1-ps)*Ssr + 2*q*Srr + Sra + 2*q*(1-ps)*Ssra + 2*q*Srra +
                                     Er + 2*q*(1-ps)*Esr + 2*q*Err + Era + 2*q*(1-ps)*Esra + 2*q*Erra +
                                     Ir + 2*q*(1-ps)*Isr + 2*q*Irr + Ira + 2*q*(1-ps)*Isra + 2*q*Irra +
                                     Iraz + 2*q*(1-ps)*Israz + 2*q*Irraz +
                                     Rr + 2*q*(1-ps)*Rsr + 2*q*Rrr + Rra + 2*q*(1-ps)*Rsra + 2*q*Rrra + 
                                     Rraz + 2*q*(1-ps)*Rsraz + 2*q*Rrraz)/N
    
    
    lambdaXS <- betaS*thetab_10*(Ss + psingle*(2*q*ps*Ssr + 2*q*Sss) + Ssa + psingle*(2*q*ps*Ssra + 2*q*Sssa) + 
                                  Es + psingle*(2*q*ps*Esr + 2*q*Ess) + Esa + psingle*(2*q*ps*Esra + 2*q*Essa) + 
                                  Is + psingle*(2*q*ps*Isr + 2*q*Iss) + Isa + psingle*(2*q*ps*Isra + 2*q*Issa) + 
                                  Isaz + psingle*(2*q*ps*Israz + 2*q*Issaz) +
                                  Rs + psingle*(2*q*ps*Rsr + 2*q*Rss) + Rsa + psingle*(2*q*ps*Rsra + 2*q*Rssa) + 
                                  Rsaz + psingle*(2*q*ps*Rsraz + 2*q*Rssaz))/N
    
    
    lambdaXR <- betaS*fit*thetab_10*(Sr + psingle*(2*q*(1-ps)*Ssr + 2*q*Srr) + Sra + psingle*(2*q*(1-ps)*Ssra + 2*q*Srra) + 
                                      Er + psingle*(2*q*(1-ps)*Esr + 2*q*Err) + Era + psingle*(2*q*(1-ps)*Esra + 2*q*Erra) + 
                                      Ir + psingle*(2*q*(1-ps)*Isr + 2*q*Irr) + Ira + psingle*(2*q*(1-ps)*Isra + 2*q*Irra) + 
                                      Iraz + psingle*(2*q*(1-ps)*Israz + 2*q*Irraz) + 
                                      Rr + psingle*(2*q*(1-ps)*Rsr + 2*q*Rrr) + Rra + psingle*(2*q*(1-ps)*Rsra + 2*q*Rrra) + 
                                      Rraz + psingle*(2*q*(1-ps)*Rsraz + 2*q*Rrraz))/N
    
    
    lambdaXSR <- ((1-psingle)*(betaS*thetab_10*ps + betaS*thetab_10*fit*(1-ps))*2*q*(Ssr+Ssra) + 
                    (1-psingle)*(betaS*thetab_10*ps + betaS*thetab_10*fit*(1-ps))*2*q*(Esr+Esra) + 
                    (1-psingle)*(betaS*thetab_10*ps + betaS*thetab_10*fit*(1-ps))*2*q*(Isr+Isra+Israz) +
                    (1-psingle)*(betaS*thetab_10*ps + betaS*thetab_10*fit*(1-ps))*2*q*(Rsr+Rsra+Rsraz))/N
    
    lambdaXSS <- ((1-psingle)*betaS*thetab_10*2*q*(Sss+Sssa) + 
                    (1-psingle)*betaS*thetab_10*2*q*(Ess+Essa) + 
                    (1-psingle)*betaS*thetab_10*2*q*(Iss+Issa+Issaz) + 
                    (1-psingle)*betaS*thetab_10*2*q*(Rss+Rssa+Rssaz))/N
    
    lambdaXRR <- ((1-psingle)*betaS*thetab_10*fit*2*q*(Srr+Srra) + 
                    (1-psingle)*betaS*thetab_10*fit*2*q*(Err+Erra) + 
                    (1-psingle)*betaS*thetab_10*fit*2*q*(Irr+Irra+Irraz) + 
                    (1-psingle)*betaS*thetab_10*fit*2*q*(Rrr+Rrra+Rrraz))/N
    
    # Force of infection (SARS-CoV-2)
    lambdaC <- betaC*thetabeta_10*((I+Is+Ir+Isr+Iss+Irr+Ia+Isa+Ira+Isra+Issa+Irra+
                                     Iaz+Isaz+Iraz+Israz+Issaz+Irraz)/N)
    
    ########################################## Model equations: #############################################
    
    # Susceptible to SARS-CoV-2: 
    dS <- -lambdaXS*S - lambdaXR*S - lambdaXSR*S - lambdaXSS*S - lambdaXRR*S + 
      gammaS*Ss + gammaS*Sss + gammaR*Sr + gammaR*Ssr + gammaR*Srr + 
      gammaR*Sra + gammaR*Srra + (gammaS+w)*Ssa + (gammaS+w)*Sssa - 
      tau*a_10*S + r*Sa - lambdaC*S - covid  
    
    dSs <- lambdaXS*S - gammaS*Ss - k*lambdaS*Ss - k*lambdaR*Ss - tau*a_10*Ss - lambdaC*Ss  
    dSr <- lambdaXR*S - gammaR*Sr - k*lambdaS*Sr - k*lambdaR*Sr - tau*a_10*Sr + r*Sra + (gammaS+w)*Ssra - lambdaC*Sr 
    
    dSsr <- lambdaXSR*S - gammaR*Ssr + k*lambdaR*Ss + k*lambdaS*Sr - 
      k*c*lambdaS*Ssr - k*c*lambdaR*Ssr + 2*k*c*lambdaR*Sss + 2*k*c*lambdaS*Srr - 
      tau*a_10*Ssr - lambdaC*Ssr
    
    dSss <- lambdaXSS*S - gammaS*Sss + k*lambdaS*Ss + k*c*lambdaS*Ssr - 2*k*c*lambdaR*Sss - 
      tau*a_10*Sss - lambdaC*Sss
    
    dSrr <- lambdaXRR*S - gammaR*Srr + k*lambdaR*Sr + k*c*lambdaR*Ssr - 2*k*c*lambdaS*Srr - 
      tau*a_10*Srr + r*Srra - lambdaC*Srr
    
    # Antibiotic exposed compartments:
    dSa <- -lambdaXR*Sa - lambdaXRR*Sa + tau*a_10*S - r*Sa - lambdaC*Sa
    dSsa <- -(gammaS+w)*Ssa + tau*a_10*Ss - k*lambdaR*Ssa - lambdaC*Ssa
    dSra <- - gammaR*Sra + lambdaXR*Sa + tau*a_10*Sr - r*Sra - lambdaC*Sra
    dSsra <- -(gammaS+w)*Ssra + k*lambdaR*Ssa + 2*k*c*lambdaR*Sssa + tau*a_10*Ssr - lambdaC*Ssra
    dSssa <- -(gammaS+w)*Sssa - 2*k*c*lambdaR*Sssa + tau*a_10*Sss - lambdaC*Sssa
    dSrra <-  lambdaXRR*Sa + tau*a_10*Srr - r*Srra - gammaR*Srra - lambdaC*Srra
    
    
    # Exposed to SARS-CoV-2: 
    dE <- -lambdaXS*E - lambdaXR*E - lambdaXSR*E - lambdaXSS*E - lambdaXRR*E + 
      gammaS*Es + gammaS*Ess + gammaR*Er + gammaR*Esr + gammaR*Err +
      (gammaS+w)*Esa + (gammaS+w)*Essa + gammaR*Era + gammaR*Erra - 
      tau*a_10*E + r*Ea + lambdaC*S - alpha*E   
    
    dEs <- lambdaXS*E - gammaS*Es - k*lambdaS*Es - k*lambdaR*Es - tau*a_10*Es + lambdaC*Ss - alpha*Es
    dEr <- lambdaXR*E - gammaR*Er - k*lambdaS*Er - k*lambdaR*Er - tau*a_10*Er + r*Era + lambdaC*Sr - alpha*Er + (gammaS+w)*Esra
    
    dEsr <- lambdaXSR*E - gammaR*Esr + k*lambdaR*Es + k*lambdaS*Er - 
      k*c*lambdaS*Esr - k*c*lambdaR*Esr + 2*k*c*lambdaR*Ess + 2*k*c*lambdaS*Err - 
      tau*a_10*Esr + lambdaC*Ssr - alpha*Esr 
    
    dEss <- lambdaXSS*E - gammaS*Ess + k*lambdaS*Es + k*c*lambdaS*Esr - 2*k*c*lambdaR*Ess - 
      tau*a_10*Ess + lambdaC*Sss - alpha*Ess
    
    dErr <- lambdaXRR*E - gammaR*Err + k*lambdaR*Er + k*c*lambdaR*Esr - 2*k*c*lambdaS*Err - 
      tau*a_10*Err + r*Erra + lambdaC*Srr - alpha*Err
    
    # Antibiotic exposed compartments:
    dEa <- -lambdaXR*Ea - lambdaXRR*Ea + tau*a_10*E - r*Ea + lambdaC*Sa - alpha*Ea
    dEsa <- -(gammaS+w)*Esa - k*lambdaR*Esa + tau*a_10*Es + lambdaC*Ssa - alpha*Esa
    dEra <- lambdaXR*Ea - gammaR*Era + tau*a_10*Er - r*Era + lambdaC*Sra - alpha*Era
    dEsra <- -(gammaS+w)*Esra + k*lambdaR*Esa + 2*k*c*lambdaR*Essa + tau*a_10*Esr + lambdaC*Ssra - alpha*Esra
    dEssa <- -(gammaS+w)*Essa - 2*k*c*lambdaR*Essa + tau*a_10*Ess + lambdaC*Sssa - alpha*Essa
    dErra <- lambdaXRR*Ea + tau*a_10*Err - r*Erra - gammaR*Erra + lambdaC*Srra - alpha*Erra
    
    
    # Infected with SARS-CoV-2: 
    dI <- -lambdaXS*I - lambdaXR*I - lambdaXSR*I - lambdaXSS*I - lambdaXRR*I + 
      gammaS*Is + gammaS*Iss + gammaR*Ir + gammaR*Isr + gammaR*Irr + 
      (gammaS+w)*Isa + (gammaS+w)*Issa + gammaR*Ira + gammaR*Irra - 
      (1-p_az_10)*tau*a_10*I + r*Ia + alpha*E + covid - gammaC*I - p_az_10*I + 
      gammaR*Iraz + gammaR*Irraz
    
    dIs <- lambdaXS*I - gammaS*Is - k*lambdaS*Is - k*lambdaR*Is - (1-p_az_10)*tau*a_10*Is + alpha*Es - gammaC*Is - p_az_10*Is
    dIr <- lambdaXR*I - gammaR*Ir - k*lambdaS*Ir - k*lambdaR*Ir + (gammaS+w)*Isra - (1-p_az_10)*tau*a_10*Ir + r*Ira + alpha*Er - gammaC*Ir - p_az_10*Ir 
    
    dIsr <- lambdaXSR*I - gammaR*Isr + k*lambdaR*Is + k*lambdaS*Ir - 
      k*c*lambdaS*Isr - k*c*lambdaR*Isr + 2*k*c*lambdaR*Iss + 2*k*c*lambdaS*Irr - 
      (1-p_az_10)*tau*a_10*Isr + alpha*Esr - gammaC*Isr - p_az_10*Isr 
    
    dIss <- lambdaXSS*I - gammaS*Iss + k*lambdaS*Is + k*c*lambdaS*Isr - 2*k*c*lambdaR*Iss - 
      (1-p_az_10)*tau*a_10*Iss + alpha*Ess - gammaC*Iss - p_az_10*Iss 
    
    dIrr <- lambdaXRR*I - gammaR*Irr + k*lambdaR*Ir + k*c*lambdaR*Isr - 2*k*c*lambdaS*Irr - 
      (1-p_az_10)*tau*a_10*Irr + r*Irra + alpha*Err - gammaC*Irr - p_az_10*Irr 
    
    # Antibiotic exposed compartments:
    dIa <- -lambdaXR*Ia - lambdaXRR*Ia + (1-p_az_10)*tau*a_10*I - r*Ia + alpha*Ea - gammaC*Ia 
    dIsa <- -(gammaS+w)*Isa - k*lambdaR*Isa + (1-p_az_10)*tau*a_10*Is + alpha*Esa - gammaC*Isa
    dIra <- lambdaXR*Ia - gammaR*Ira + (1-p_az_10)*tau*a_10*Ir - r*Ira + alpha*Era - gammaC*Ira
    dIsra <- -(gammaS+w)*Isra + k*lambdaR*Isa + 2*k*c*lambdaR*Issa + (1-p_az_10)*tau*a_10*Isr + alpha*Esra - gammaC*Isra 
    dIssa <- -(gammaS+w)*Issa - 2*k*c*lambdaR*Issa + (1-p_az_10)*tau*a_10*Iss + alpha*Essa - gammaC*Issa  
    dIrra <-  lambdaXRR*Ia + (1-p_az_10)*tau*a_10*Irr - r*Irra - gammaR*Irra + alpha*Erra - gammaC*Irra 
    
    # Additional azithromycin (az) exposed compartments:
    dIaz <- p_az_10*I - gammaC*Iaz - lambdaXR*Iaz - lambdaXRR*Iaz + (gammaS+w)*Isaz + (gammaS+w)*Issaz 
    dIsaz <- p_az_10*Is - gammaC*Isaz - (gammaS+w)*Isaz - k*lambdaR*Isaz   
    dIraz <- p_az_10*Ir - gammaC*Iraz + lambdaXR*Iaz - gammaR*Iraz + (gammaS+w)*Israz 
    dIsraz <- p_az_10*Isr - gammaC*Israz - (gammaS+w)*Israz + k*lambdaR*Isaz + 2*k*c*lambdaR*Issaz  
    dIssaz <- p_az_10*Iss - gammaC*Issaz - (gammaS+w)*Issaz - 2*k*c*lambdaR*Issaz   
    dIrraz <- p_az_10*Irr - gammaC*Irraz + lambdaXRR*Iaz - gammaR*Irraz 
    
    
    # Recovered from SARS-CoV-2: 
    dR <- -lambdaXS*R - lambdaXR*R - lambdaXSR*R - lambdaXSS*R - lambdaXRR*R + 
      gammaS*Rs + gammaS*Rss + gammaR*Rr + gammaR*Rsr + gammaR*Rrr + 
      gammaR*Rra + gammaR*Rrra + (gammaS+w)*Rsa + (gammaS+w)*Rssa - 
      tau*a_10*R + r*Ra + gammaC*I + r_az*Raz + 
      gammaR*Rraz + gammaR*Rrraz
    
    dRs <- lambdaXS*R - gammaS*Rs - k*lambdaS*Rs - k*lambdaR*Rs - tau*a_10*Rs + gammaC*Is + r_az*Rsaz
    dRr <- lambdaXR*R - gammaR*Rr - k*lambdaS*Rr - k*lambdaR*Rr + (gammaS+w)*Rsra - tau*a_10*Rr + r*Rra + gammaC*Ir + r_az*Rraz
    
    dRsr <- lambdaXSR*R - gammaR*Rsr + k*lambdaR*Rs + k*lambdaS*Rr - 
      k*c*lambdaS*Rsr - k*c*lambdaR*Rsr + 2*k*c*lambdaR*Rss + 2*k*c*lambdaS*Rrr - 
      tau*a_10*Rsr + gammaC*Isr + r_az*Rsraz
    
    dRss <- lambdaXSS*R - gammaS*Rss + k*lambdaS*Rs + k*c*lambdaS*Rsr - 2*k*c*lambdaR*Rss - 
      tau*a_10*Rss + gammaC*Iss + r_az*Rssaz
    
    dRrr <- lambdaXRR*R - gammaR*Rrr + k*lambdaR*Rr + k*c*lambdaR*Rsr - 2*k*c*lambdaS*Rrr - 
      tau*a_10*Rrr + r*Rrra + gammaC*Irr + r_az*Rrraz
    
    # Antibiotic exposed compartments:
    dRa <- gammaC*Ia + tau*a_10*R - r*Ra - lambdaXR*Ra - lambdaXRR*Ra 
    dRsa <- gammaC*Isa + tau*a_10*Rs - (gammaS+w)*Rsa - k*lambdaR*Rsa
    dRra <- gammaC*Ira + tau*a_10*Rr - r*Rra + lambdaXR*Ra - gammaR*Rra
    dRsra <- gammaC*Isra + tau*a_10*Rsr - (gammaS+w)*Rsra + k*lambdaR*Rsa + 2*k*c*lambdaR*Rssa 
    dRssa <- gammaC*Issa + tau*a_10*Rss - (gammaS+w)*Rssa - 2*k*c*lambdaR*Rssa
    dRrra <- gammaC*Irra + tau*a_10*Rrr - r*Rrra + lambdaXRR*Ra - gammaR*Rrra 
    
    # Additional antibiotic (az) exposed compartments:
    dRaz <- gammaC*Iaz - r_az*Raz - lambdaXR*Raz - lambdaXRR*Raz + (gammaS+w)*Rsaz + (gammaS+w)*Rssaz
    dRsaz <- gammaC*Isaz - r_az*Rsaz - (gammaS+w)*Rsaz - k*lambdaR*Rsaz
    dRraz <- gammaC*Iraz - r_az*Rraz + lambdaXR*Raz - gammaR*Rraz + (gammaS+w)*Rsraz
    dRsraz <- gammaC*Israz - r_az*Rsraz - (gammaS+w)*Rsraz + k*lambdaR*Rsaz + 2*k*c*lambdaR*Rssaz 
    dRssaz <- gammaC*Issaz - r_az*Rssaz - (gammaS+w)*Rssaz - 2*k*c*lambdaR*Rssaz
    dRrraz <- gammaC*Irraz - r_az*Rrraz + lambdaXRR*Raz - gammaR*Rrraz  
    
    # Cumulative incidence of invasive pneumococcal disease:
    
    # IPD caused by S cannot be caused by individuals receiving antibiotic treatment:
    dIPDs <- inv_rate_10*IPDrisk_10*(Ss+Sss+0.5*Ssr+Es+Ess+0.5*Esr+Is+Iss+0.5*Isr+Rs+Rss+0.5*Rsr)
    
    # IPD caused by R can be caused by all individuals carrying R:
    dIPDr <- inv_rate_10*IPDrisk_10*(Sr+Srr+0.5*Ssr+Er+Err+0.5*Esr+Ir+Irr+0.5*Isr+Rr+Rrr+0.5*Rsr+
                                     Sra+Srra+Ssra+Era+Erra+Esra+Ira+Irra+Isra+Rra+Rrra+Rsra+
                                     Iraz+Irraz+Israz+Rraz+Rrraz+Rsraz)
    
    dIPD <- inv_rate_10*IPDrisk_10*(Ss+Sss+0.5*Ssr+Es+Ess+0.5*Esr+Is+Iss+0.5*Isr+Rs+Rss+0.5*Rsr)+
      inv_rate_10*IPDrisk_10*(Sr+Srr+0.5*Ssr+Er+Err+0.5*Esr+Ir+Irr+0.5*Isr+Rr+Rrr+0.5*Rsr+
                              Sra+Srra+Ssra+Era+Erra+Esra+Ira+Irra+Isra+Rra+Rrra+Rsra+
                              Iraz+Irraz+Israz+Rraz+Rrraz+Rsraz)
    
    # Cumulative number of antibiotic exposed individuals:
    dATB <- tau*a_10*(S+Ss+Ssr+Ssr+Sss+Srr+E+Es+Er+Esr+Ess+Err) +
      (1-p_az_10)*tau*a_10*(I+Is+Ir+Isr+Iss+Irr) + 
      p_az_10*tau*a_10*(I+Is+Ir+Isr+Iss+Irr) + 
      tau*a_10*(R+Rs+Rr+Rsr+Rss+Rrr) 
    
    return(list(c(dS, dSs, dSr, dSsr, dSss, dSrr, dSa, dSsa, dSra, dSsra, dSssa, dSrra,  
                  dE, dEs, dEr, dEsr, dEss, dErr, dEa, dEsa, dEra, dEsra, dEssa, dErra,
                  dI, dIs, dIr, dIsr, dIss, dIrr, dIa, dIsa, dIra, dIsra, dIssa, dIrra,
                  dIaz, dIsaz, dIraz, dIsraz, dIssaz, dIrraz,
                  dR, dRs, dRr, dRsr, dRss, dRrr, dRa, dRsa, dRra, dRsra, dRssa, dRrra, 
                  dRaz, dRsaz, dRraz, dRsraz, dRssaz, dRrraz, 
                  dIPDs, dIPDr, dIPD, dATB))) 
  })
  
}

output_10 <- as.data.frame(ode(y = initial_state, 
                              times = sim_time, 
                              func = model_10, 
                              parms = parameters_10, 
                              method = "euler"))
#####################################################################################################################



##################### Pneumococcal carriage prevalence over time:
# Total pneumococcal carriage
out10 <- dplyr::mutate(output_10, PC_10 = Ss + Sr + Ssr + Sss + Srr + Ssa + Sra + Ssra + Sssa + Srra + 
                        Es + Er + Esr + Ess + Err + Esa + Era + Esra + Essa + Erra + 
                        Is + Ir + Isr + Iss + Irr + Isa + Ira + Isra + Issa + Irra + 
                        Isaz + Iraz + Israz + Issaz + Irraz +
                        Rs + Rr + Rsr + Rss + Rrr + Rsa + Rra + Rsra + Rssa + Rrra +
                        Rsaz + Rraz + Rsraz + Rssaz + Rrraz,
                      R_10 = Sr + 0.5*Ssr + Srr + Sra + 0.5*Ssra + Srra + 
                        Er + 0.5*Esr + Err + Era + 0.5*Esra + Erra + 
                        Ir + 0.5*Isr + Irr + Ira + 0.5*Isra + Irra + 
                        Iraz + 0.5*Israz + Irraz +
                        Rr + 0.5*Rsr + Rrr + Rra + 0.5*Rsra + Rrra + 
                        Rraz + 0.5*Rsraz + Rrraz,
                      S_10 = Ss + 0.5*Ssr + Sss + Ssa + 0.5*Ssra + Sssa + 
                        Es + 0.5*Esr + Ess + Esa + 0.5*Esra + Essa + 
                        Is + 0.5*Isr + Iss + Isa + 0.5*Isra + Issa + 
                        Isaz + 0.5*Israz + Issaz + 
                        Rs + 0.5*Rsr + Rss + Rsa + 0.5*Rsra + Rssa + 
                        Rsaz + 0.5*Rsraz + Rssaz)



# Cumulative IPD incidence short and long term:------------------------------------
# Cumulative ipd incidence (S, R, SR):
# short term and long term
total_short_10 <- output_10$IPD[1135] - output_10$IPD[1075] # per 100,000 people cumulative ipd incidence during the 60-day period
total_long_10 <- output_10$IPD[1365] - output_10$IPD[1000] # per 100,000 people cumulative ipd incidence during one year since the outbreak

# IPD incidence caused by susceptible strain only (S + 0.5*SR):
# short term
S_short_10 <- output_10$IPDs[1135] - output_10$IPDs[1075] # per 100,000 people cumulative ipd incidence during the 60-day period
S_long_10 <- output_10$IPDs[1365] - output_10$IPDs[1000] # per 100,000 people cumulative ipd incidence during one year since the outbreak

# IPD incidence caused by resistant strain only (R + 0.5*SR):
# short term
R_short_10 <- output_10$IPDr[1135] - output_10$IPDr[1075] # per 100,000 people cumulative ipd incidence during the 60-day period
R_long_10 <- output_10$IPDr[1365] - output_10$IPDr[1000] # per 100,000 people cumulative ipd incidence during one year since the outbreak



#:::::::::::::::::::: Pneumococcal trends in 2020 using scenario S19: ::::::::::::
##################################################################################
######### IPD inc. ################
round(total_long_10, digits = 1)                  # Annual total IPD incidence in 2020
R_long_10                                         # Annual resistant IPD incidence in 2020
######### AR (%) #################
AR_a10 <- (R_long_10/total_long_10)*100  
round(AR_a10, digits = 1)                         # ANTIBIOTIC RESISTANCE (%) in 2020
######### Sp. (%) #################
# change in the prevalence of TOTAL bacterial carriage at the end of the first 60-day lockdown period
relative_percent_change_10 <- ((out10$PC_10[1135] - out10$PC_10[1075])/out10$PC_10[1075])*100
print('% change in the prevalence of TOTAL pneumococcal carriage at the end of the 60-day period')
round(relative_percent_change_10, digits =1)
##################################################################################

# The output above corresponds to Table 1: scenario S19; all other scenarios can be run by activating/deactivating mechanisms (1-5)
# labelled throughout