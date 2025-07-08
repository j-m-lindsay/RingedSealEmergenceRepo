# Header---------------------------
# Code associated with the paper "Ringed seal (Pusa hispida) haul-out behavior and emergence timing in the Bering, Chukchi, and Beaufort seas"
# Authors: Jessica M. Lindsay, Paul B. Conn, Peter L. Boveng, Justin A. Crawford, Lori T. Quakenbush, Andrew L. Von Duyke, & Kristin L. Laidre
# Contact: jessica.lindsay@noaa.gov
# Date Created: 2025-06-24
# Description: code to fit Hidden Markov Models (HMMs) to haul-out data from subadult ringed seals, estimate emergence dates for individual 
# seals, and generate the corresponding figures and summary statistics reported in the manuscript. The code for the 
# analysis is identical to the code within the adult script. The main differences between this script and the one for adults
# are a) slightly different figures, since the covariates from the best subadult model are different from the those in the best adult model,
# and b) no regression performed on the relationship between emergence date and latitude due to the smaller subadult sample size. 
#
# Notes for use --------------------
# The code below uses the here() package and assumes you have a "data" folder as indicated in the GitHub repo structure
# Code was written in R version 4.3.0 "Already Tomorrow" and RStudio version 2023.12.1+402 "Ocean Storm"
#


#Load packages, functions, data-----------------
rm(list = ls()) #clear environment

library(pacman)
p_load(tidyverse, momentuHMM, here, ggplot2, RColorBrewer, circular, gridExtra, VGAM, lubridate)

#functions for converting back and forth from hours to radians
hrs_to_rad <- function(hrs){
  return((pi/12)*(hrs-12))
}

rad_to_hrs <- function(rad){
  12*rad/pi+12
}


#function for converting DOY to "day month" format for easier interpretation
doytodate <- function(x){
  paste(day(as.Date(round(x,0), origin = "2010-12-31"))," ",
        lubridate::month(as.Date(round(x,0), origin = "2010-12-31"),label=TRUE))
}

#load the cleaned data
tags <- read.csv(here("data","ringed_ho_data_clean.csv")) %>%
  mutate_at(c("speno_yr","age_class","sex","data_source"), factor)

#HMMs---------------
#___format data for momentuHMM & set HMM basics-----------
# Heavily based on the manual and associated vignettes for the momentuHMM package 
# https://cran.r-project.org/web/packages/momentuHMM/index.html and https://github.com/bmcclintock/momentuHMM 

#format tag data for momentuHMM
dat <- tags %>%
 filter(age_class=="Subadult")%>%
 mutate(daily_prop_ho=daily_prop_ho/100)%>% #convert percentage to proportion in order to be compatible with beta distribution
 select(c(speno_yr, doy, iceconc, CM_index, EM_index, daily_prop_ho, peak_ho_hr, age_class, sex, year, lat, air2m, lon, daylength, rolling_temp, atdd ))%>%
 rename(ID=speno_yr) #momentuHMM automatically recognizes "ID" as the unique animal identifier

dat <- prepData(data=dat, covNames = c("doy","iceconc","CM_index","EM_index","age_class",
                                       "sex","year","lat","air2m","lon","daylength","rolling_temp","atdd"), coordNames=NULL)
#need to specifiy covNames because momentuHMM automatically interprets any columns that aren't "ID" or in the covNames list as observed variables
#the covNames list above leaves our observed variables as daily_prop_ho and peak_ho_hr

#specify the probability distribution types for each of the observed variables
dist <- list(daily_prop_ho="beta",peak_ho_hr="vm") #beta & von mises

#label state 1 "lair" and state 2 "emerged"
stateNames <- c("lair","emerged")

#set starting parameter values for the state-dependent probability distributions 
Par0_m1 <- list(daily_prop_ho=c(1,8, #lair shape1, emerged shape1
                                 8,4, #lair shape2, emerged shape2
                                 0.75,0.01, #lair zero-mass, emerged zero-mass
                                 0.01,0.1), #lair one-mass, emerged one-mass
                peak_ho_hr=c(pi,0, #lair mu, emerged mu
                             0.5,0.2)) #lair sd, emerged sd

#label all data from February (doy <60) as belonging to the lair state (i.e., state 1)
knownStates <- rep(NA,nrow(dat))
knownStates[which(dat$doy < 60)] <- 1

#prevent seals from switching from the emerged state to the lair state by setting the 2-->1 transition probability to a large negative number (-100)
#(see momentuHMM vignette)
fixedbeta <- matrix(c(NA,-100),nrow=1) #for null model with no covariates on transition probability
fixedbeta1 <- matrix(c(NA,-100,
                       NA, 0),nrow=2, byrow=TRUE) #for models with one covariate on transition probability
fixedbeta2 <- matrix(c(NA,-100,
                       NA, 0,
                       NA, 0),nrow=3, byrow=TRUE) #for models with two covariates on transition probability
fixedbeta3 <- matrix(c(NA,-100,
                       NA, 0,
                       NA, 0,
                       NA, 0),nrow=4, byrow=TRUE) #for models with three covariates on transition probability


#___fit candidate model set-----------
#fit null model (no covariates)
m_null <- fitHMM(data = dat, nbStates = 2, dist = dist, Par0 = Par0_m1, estAngleMean = list(peak_ho_hr=TRUE), 
             stateNames = stateNames, beta0=NULL, knownStates = knownStates, fixPar=list(beta=fixedbeta)) 

#now can add covariates on emergence probability

#doy as covariate
formula_doy <- ~ doy #formula for transition probabilities, now influenced by covariate
Par0_doy <- getPar0(model=m_null, formula=formula_doy) #obtains the parameter values from the null model, reformatted to now have doy as a covariate on transition probability
m_doy <- fitHMM(data = dat, nbStates = 2, dist = dist, Par0 = Par0_doy$Par,
             estAngleMean = list(peak_ho_hr=TRUE), stateNames = stateNames, beta0=Par0_doy$beta, delta0=Par0_doy$delta, formula=formula_doy, knownStates = knownStates,
             fixPar=list(beta=fixedbeta1))

#CM as covariate
formula_CM <- ~ CM_index  
Par0_CM <- getPar0(model=m_null, formula=formula_CM)
m_CM <- fitHMM(data = dat, nbStates = 2, dist = dist, Par0 = Par0_CM$Par,
             estAngleMean = list(peak_ho_hr=TRUE), stateNames = stateNames, beta0=Par0_CM$beta, delta0=Par0_CM$delta, formula=formula_CM, knownStates = knownStates,fixPar=list(beta=fixedbeta1))

#EM as covariate
formula_EM <- ~ EM_index  
Par0_EM <- getPar0(model=m_null, formula=formula_EM)
m_EM <- fitHMM(data = dat, nbStates = 2, dist = dist, Par0 = Par0_EM$Par,
             estAngleMean = list(peak_ho_hr=TRUE), stateNames = stateNames, beta0=Par0_EM$beta, delta0=Par0_EM$delta, formula=formula_EM, knownStates = knownStates,fixPar=list(beta=fixedbeta1))

#iceconc as covariate
formula_ice <- ~ iceconc  
Par0_ice <- getPar0(model=m_null, formula=formula_ice)
m_iceconc <- fitHMM(data = dat, nbStates = 2, dist = dist, Par0 = Par0_ice$Par,
             estAngleMean = list(peak_ho_hr=TRUE), stateNames = stateNames, beta0=Par0_ice$beta, delta0=Par0_ice$delta, formula=formula_ice, knownStates = knownStates,fixPar=list(beta=fixedbeta1))

#air2m as covariate
formula_air <- ~ air2m  
Par0_air <- getPar0(model=m_null, formula=formula_air)
m_air2m <- fitHMM(data = dat, nbStates = 2, dist = dist, Par0 = Par0_air$Par,
             estAngleMean = list(peak_ho_hr=TRUE), stateNames = stateNames, beta0=Par0_air$beta, delta0=Par0_air$delta, formula=formula_air, knownStates = knownStates,fixPar=list(beta=fixedbeta1))

#daylength as covariate
formula_DL <- ~ daylength  
Par0_DL <- getPar0(model=m_null, formula=formula_DL)
m_daylength <- fitHMM(data = dat, nbStates = 2, dist = dist, Par0 = Par0_DL$Par,
                  estAngleMean = list(peak_ho_hr=TRUE), stateNames = stateNames, beta0=Par0_DL$beta, delta0=Par0_DL$delta, formula=formula_DL, knownStates = knownStates,fixPar=list(beta=fixedbeta1))

#rolling_temp as covariate
formula_rollT <- ~ rolling_temp  
Par0_rollT <- getPar0(model=m_null, formula=formula_rollT)
m_rolling_temp <- fitHMM(data = dat, nbStates = 2, dist = dist, Par0 = Par0_rollT$Par,
                         estAngleMean = list(peak_ho_hr=TRUE), stateNames = stateNames, beta0=Par0_rollT$beta, delta0=Par0_rollT$delta, formula=formula_rollT, knownStates = knownStates,fixPar=list(beta=fixedbeta1))
#ATDD as covariate
formula_atdd <- ~ atdd  
Par0_atdd <- getPar0(model=m_null, formula=formula_atdd)
m_atdd <- fitHMM(data = dat, nbStates = 2, dist = dist, Par0 = Par0_atdd$Par,
                   estAngleMean = list(peak_ho_hr=TRUE), stateNames = stateNames, beta0=Par0_atdd$beta, delta0=Par0_atdd$delta, formula=formula_atdd, knownStates = knownStates,fixPar=list(beta=fixedbeta1))

#daylength+air2m as covariates
formula_DL_air <- ~ air2m+daylength  
Par0_DL_air <- getPar0(model=m_daylength, formula=formula_DL_air)
m_daylength_air2m <- fitHMM(data = dat, nbStates = 2, dist = dist, Par0 = Par0_DL_air$Par,
                            estAngleMean = list(peak_ho_hr=TRUE), stateNames = stateNames, beta0=Par0_DL_air$beta, formula=formula_DL_air, knownStates = knownStates,fixPar=list(beta=fixedbeta2))

#daylength+roling_temp as covariate
formula_DL_rollT <- ~ daylength+rolling_temp  
Par0_DL_rollT <- getPar0(model=m_daylength, formula=formula_DL_rollT)
m_daylength_temp <- fitHMM(data = dat, nbStates = 2, dist = dist, Par0 = Par0_DL_rollT$Par,
                           estAngleMean = list(peak_ho_hr=TRUE), stateNames = stateNames, beta0=Par0_DL_rollT$beta, formula=formula_DL_rollT, knownStates = knownStates,fixPar=list(beta=fixedbeta2))

#doy*sex  as covariate
formula_doy_sex <- ~ doy+sex  
Par0_doy_sex <- getPar0(model=m_doy, formula=formula_doy_sex)
m_doy_sex <- fitHMM(data = dat, nbStates = 2, dist = dist, Par0 = Par0_doy_sex$Par,
                    estAngleMean = list(peak_ho_hr=TRUE), stateNames = stateNames, beta0=Par0_doy_sex$beta, delta0=Par0_doy_sex$delta,formula=formula_doy_sex, knownStates = knownStates,fixPar=list(beta=fixedbeta2))
formula_doyxsex <- ~ doy*sex  
Par0_doyxsex <- getPar0(model=m_doy_sex, formula=formula_doyxsex)
m_doyxsex <- fitHMM(data = dat, nbStates = 2, dist = dist, Par0 = Par0_doyxsex$Par,
                 estAngleMean = list(peak_ho_hr=TRUE), stateNames = stateNames, beta0=Par0_doyxsex$beta, delta0=Par0_doyxsex$delta,formula=formula_doyxsex, knownStates = knownStates,fixPar=list(beta=fixedbeta3))

#doy*lat as covariate
formula_doy_lat <- ~ doy+lat  
Par0_doy_lat <- getPar0(model=m_doy, formula=formula_doy_lat)
m_doy_lat <- fitHMM(data = dat, nbStates = 2, dist = dist, Par0 = Par0_doy_lat$Par,
                    estAngleMean = list(peak_ho_hr=TRUE), stateNames = stateNames, beta0=Par0_doy_lat$beta, delta0=Par0_doy_lat$delta, formula=formula_doy_lat, knownStates = knownStates,fi_Par=list(beta=fi_edbeta2))
formula_doyxlat <- ~ doy*lat  
Par0_doyxlat <- getPar0(model=m_doy_lat, formula=formula_doyxlat)
m_doyxlat <- fitHMM(data = dat, nbStates = 2, dist = dist, Par0 = Par0_doyxlat$Par,
                estAngleMean = list(peak_ho_hr=TRUE), stateNames = stateNames, beta0=Par0_doyxlat$beta, delta0=Par0_doyxlat$delta, formula=formula_doyxlat, knownStates = knownStates,fixPar=list(beta=fixedbeta3))

#daylength*lat as covariate
formula_DL_lat <- ~ daylength+lat  
Par0_DL_lat <- getPar0(model=m_daylength, formula=formula_DL_lat)
m_daylength_lat <- fitHMM(data = dat, nbStates = 2, dist = dist, Par0 = Par0_DL_lat$Par,
                          estAngleMean = list(peak_ho_hr=TRUE), stateNames = stateNames, beta0=Par0_DL_lat$beta, delta0=Par0_DL_lat$delta, formula=formula_DL_lat, knownStates = knownStates,fixPar=list(beta=fixedbeta2))
formula_DLxlat <- ~ daylength*lat  
Par0_DLxlat <- getPar0(model=m_daylength_lat, formula=formula_DLxlat)
m_daylengthxlat <- fitHMM(data = dat, nbStates = 2, dist = dist, Par0 = Par0_DLxlat$Par,
                     estAngleMean = list(peak_ho_hr=TRUE), stateNames = stateNames, beta0=Par0_DLxlat$beta, delta0=Par0_DLxlat$delta, formula=formula_DLxlat, knownStates = knownStates,fixPar=list(beta=fixedbeta3))


#doy*year as covariate
formula_doy_yr <- ~ doy+year  
Par0_doy_yr <- getPar0(model=m_doy, formula=formula_doy_yr)
m_doy_year <- fitHMM(data = dat, nbStates = 2, dist = dist, Par0 = Par0_doy_yr$Par,
                     estAngleMean = list(peak_ho_hr=TRUE), stateNames = stateNames, beta0=Par0_doy_yr$beta, delta0=Par0_doy_yr$delta,formula=formula_doy_yr, knownStates = knownStates,fixPar=list(beta=fixedbeta2))
formula_doyxyr <- ~ doy*year  
Par0_doyxyr <- getPar0(model=m_doy_year, formula=formula_doyxyr)
m_doyxyear <- fitHMM(data = dat, nbStates = 2, dist = dist, Par0 = Par0_doyxyr$Par,
                    estAngleMean = list(peak_ho_hr=TRUE), stateNames = stateNames, beta0=Par0_doyxyr$beta, delta0=Par0_doyxyr$delta,formula=formula_doyxyr, knownStates = knownStates,fixPar=list(beta=fixedbeta3))


#___AIC---------------
AICtab <- AIC(m_null,m_doy, m_daylength, m_CM, m_EM, m_atdd, m_air2m, m_rolling_temp, m_iceconc,m_doyxlat,m_doyxyear,m_doyxsex,
           m_daylengthxlat,m_daylength_air2m,m_daylength_temp) %>%
  mutate(dAIC=round(AIC-min(AIC),2))%>%
  arrange(dAIC)
AICtab$npar <- NA
for (i in 1:nrow(AICtab)){
  AICtab$npar[i] <- length(get((AICtab$Model[i]))$mod$wpar) #extract number of parameters estimated by each candidate HMM to include this as a column
}

best_mod <- get(AICtab$Model[which(AICtab$dAIC==0)])
dat$states <- viterbi(best_mod) #use Viterbi algorithm to estimate state sequence for each seal's haul-out record


#Stats to report from best model-------------
#___Table S3: state-dependent probability distributions---------
best_mod$CIreal$daily_prop_ho$est #estimated parameters
best_mod$CIreal$daily_prop_ho$se #se for those estimates
best_mod$CIreal$peak_ho_hr$est #estimated parameters 
best_mod$CIreal$peak_ho_hr$se #se for those estimates

#___Table S4: estimated transition probabilities & regression coefficients-----
best_mod$CIreal$gamma$est #state transition probabilities when covariates are at their mean values
best_mod$CIreal$gamma$se #se for those estimates
best_mod$CIbeta$beta$est #estimated regression coefficients (intercept & covariates)
best_mod$CIbeta$beta$se #se for those estimates

#___Other stats to report in the main text-----------
lairstate <- dat %>% filter(states==1)
emergedstate <- dat %>% filter(states==2)
mean(lairstate$daily_prop_ho,na.rm=TRUE) #mean daily_prop_ho for seals in the lair state
sd(lairstate$daily_prop_ho,na.rm=TRUE)

mean(emergedstate$daily_prop_ho,na.rm=TRUE) #mean daily_prop_ho for seals in the emerged state
sd(emergedstate$daily_prop_ho,na.rm=TRUE)

#convert the mean peak haul-out hour for each state from radians to hours for easier interpretation
#mean for the lair state
paste0(floor(rad_to_hrs(best_mod$CIreal$peak_ho_hr$est[1,"lair"])),":",
       round(60*(rad_to_hrs(best_mod$CIreal$peak_ho_hr$est[1,1])-floor(rad_to_hrs(best_mod$CIreal$peak_ho_hr$est[1,1]))),0)) 
#mean for the emerged state
paste0(floor(rad_to_hrs(best_mod$CIreal$peak_ho_hr$est[1,"emerged"])),":",
       round(60*(rad_to_hrs(best_mod$CIreal$peak_ho_hr$est[1,1])-floor(rad_to_hrs(best_mod$CIreal$peak_ho_hr$est[1,1]))),0)) 


#HMM figures (Fig 6) -------------
#momentuHMM has nice built-in plots, but we preferred custom formatting that replicates those built-in plots with some tweaks 
#(e.g., different colors, labels, shaded confidence intervals instead of bars, rug to show distribution of covariate values, etc.)
#to view the built-in plots, you can use "plot(best_mod,plotCI=TRUE)"

theme_set(theme_bw()+theme(panel.background = element_blank(),  panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),panel.border=element_blank(),axis.line = element_line(colour = "black"),
                           axis.text.y = element_text(angle=360, vjust=0.5, hjust=1),plot.tag=element_text(size=21),
                           text=element_text(size=21), plot.margin = unit(c(20,20,20,20), "pt")))

#___Fig6a: Subadults; proportion of day hauled out ---------------

#extract parameter estimates for state-dependent probability distributions from the best model
#lair state 
shape1_lair <- best_mod$CIreal$daily_prop_ho$est[1] #shape1 estimate for the lair state
shape2_lair <- best_mod$CIreal$daily_prop_ho$est[2] #shape2 estimate for the lair state
zm_lair <- best_mod$CIreal$daily_prop_ho$est[3] #zero-mass estimate for the lair state
om_lair <- best_mod$CIreal$daily_prop_ho$est[4] #one-mass estimate for the lair state
#emerged state
shape1_emerged <- best_mod$CIreal$daily_prop_ho$est[5] #shape1 estimate for the emerged state
shape2_emerged <- best_mod$CIreal$daily_prop_ho$est[6] #shape2 estimate for the emerged state
zm_emerged <- best_mod$CIreal$daily_prop_ho$est[7] #extract zero-mass estimate for the emerged state
om_emerged <- best_mod$CIreal$daily_prop_ho$est[8] #extract one-mass estimate for the emerged state

#use those parameter estimates to generate density curves
#"_a" suffix indicates that objects are for panel a of Fig4
xseq_a <- seq(0,1,by=0.001)
dcurv_propho_lair <- data.frame(x=xseq_a,y=dzoabeta(x=xseq_a,shape1=shape1_lair,shape2=shape2_lair, pobs0=zm_lair,pobs1=om_lair))
dcurv_propho_emerged <- data.frame(x=xseq_a,y=dzoabeta(x=xseq_a,shape1=shape1_emerged,shape2=shape2_emerged, pobs0=zm_emerged,pobs1=om_emerged))

#convert the density curves into curves of expected frequency
nobs_a <- sum(complete.cases(dat$daily_prop_ho)) #how many observations of daily_prop_ho are actually included in the panel a histogram (histogram only plots non-NA values of daily_prop_ho)  
frac_lair_a <- length(dat$states[which(dat$states==1 & is.na(dat$daily_prop_ho) == FALSE)])/nobs_a #what fraction of those non-NA oaservations were estimated to have belonged to the lair state 
frac_emerged_a <- length(dat$states[which(dat$states==2 & is.na(dat$daily_prop_ho) == FALSE)])/nobs_a #what fraction of those non-NA oaservations were estimated to have belonged to the emerged state
nbins_a <- 10 #number of bins to use in panel a histogram
ht_a <- 1/nbins_a*nobs_a #height modifier for panel a curves. the y-scale for the curves is density, whereas the y-scale of our histogram 
#is frequency. To convert density to frequency so that the curves can appropriately overlay the histogram, we need to multiply the y-values
#from the density curves above by the area of the histogram. this is calculated as the binwidth of the histogram (1/nbins_a) multiplied by 
#the number of observations in the histogram. We then have to multiply this by the fraction of those observations that the HMM assigned to each
#behavioral state (frac_lair_a or frac_emerged_a)
fcurv_propho_lair <- dcurv_propho_lair %>% mutate(y=y*ht_a*frac_lair_a)
fcurv_propho_emerged <- dcurv_propho_emerged %>% mutate(y=y*ht_a*frac_emerged_a)
fcurv_propho_sum <- bind_cols(x=xseq_a,y=fcurv_propho_lair$y+fcurv_propho_emerged$y) #sum of the two frequency curves

#the density & frequency curves generated above do not reflect the zero-mass or one-mass, so we need to draw that separately if we want it in our plot
nzero_lair <- nobs_a*frac_lair_a*zm_lair #number of exact zeros predicted in lair state 
none_lair <- nobs_a*frac_lair_a*om_lair #number of exact ones predicted in lair state 
nzero_emerged <- nobs_a*frac_emerged_a*zm_emerged #number of exact zeros predicted in emerged state 
none_emerged <- nobs_a*frac_emerged_a*om_emerged #number of exact ones predicted in emerged state  

fig6a <- dat %>%
  filter(complete.cases(daily_prop_ho))%>%
  ggplot()+
  labs(tag="a")+
  theme(plot.tag.position = c(0.2,1))+
  geom_histogram(aes(x=daily_prop_ho),boundary=0,binwidth=0.1,closed="right",fill="gray88",col="gray80") +
  scale_x_continuous(limits=c(0,1),breaks=seq(0,1,by=0.2),expand=c(0.02,0))+
  scale_y_continuous(limits=c(0,650), expand=c(0,0)) +
  geom_line(data=fcurv_propho_emerged, aes(x=x,y=y),col="orange1",lwd=1)+
  geom_line(data=fcurv_propho_sum, aes(x=x,y=y),col="black",lwd=1)+
  geom_line(data=fcurv_propho_lair, aes(x=x,y=y),col="deepskyblue1",lwd=1)+ 
  ylab("Frequency")+
  xlab("Proportion of day hauled out")+
  annotate("segment",x=0,xend=0,y=fcurv_propho_sum$y[1], yend=nzero_lair+nzero_emerged, col="black",lwd=1.5)+ #segment for expected exact 0s, both states
  annotate("segment",x=0+0.004,xend=0+0.004,y=fcurv_propho_lair$y[1],yend=nzero_lair,col="deepskyblue1",lwd=1)+ #segment for expected exact 0s, lair state. jittering the x-values a little bit so it isn't covered up by the black line
  annotate("segment",x=0,xend=0,y=fcurv_propho_emerged$y[1],yend=nzero_emerged,col="orange1",lwd=1)+ #segment for expected exact 0s, emerged state
  annotate("segment",x=1-0.004,xend=1-0.004,y=fcurv_propho_sum$y[length(xseq_a)], yend=none_lair+none_emerged, col="black",lwd=1.25)+ #segment for expected exact 1s, both states. jittering the x-values a little bit so it isn't covered up by the orange line
  annotate("segment",x=1,xend=1,y=fcurv_propho_emerged$y[length(xseq_a)],yend=none_emerged,col="orange1",lwd=1)+ #segment for expected exact 1s, emerged state
  annotate("segment",x=1,xend=1,y=fcurv_propho_lair$y[length(xseq_a)],yend=none_lair,col="deepskyblue1",lwd=1)+ #segment for expected exact 1s, lair state
  coord_cartesian(clip = 'off',expand=TRUE)

#___Fig6b: Subadults; peak haul-out hour ------------

#extract parameter estimates for state-dependent probability distributions from the best model
mu_lair <- best_mod$CIreal$peak_ho_hr$est[1] #mu estimate for the lair state
conc_lair <- best_mod$CIreal$peak_ho_hr$est[2] #concentration estimate for the lair state
mu_emerged <- best_mod$CIreal$peak_ho_hr$est[3] #mu estimate for the emerged state
conc_emerged <- best_mod$CIreal$peak_ho_hr$est[4] #concentration estimate for the emerged state

#use those parameter estimates to generate density curves
#"_b" suffix indicates that objects are for panel b of Fig4
xseq_b <- seq(-pi,pi,by=0.001)
dcurv_peakhr_lair <- data.frame(x=xseq_b,y=dvonmises(x=xseq_b,mu=mu_lair,kappa=conc_lair)) #will warn you that it is using default values for some arguments, but the defaults are appropriate so OK to disregard
dcurv_peakhr_emerged <- data.frame(x=xseq_b,y=dvonmises(x=xseq_b,mu=mu_emerged,kappa=conc_emerged)) #will warn you that it is using default values for some arguments, but the defaults are appropriate so OK to disregard 

#convert the density curves into curves of expected frequency
nobs_b <- sum(complete.cases(dat$peak_ho_hr)) #how many observations of peak_ho_hr are actually included in the panel b histogram (histogram only plots non-NA values of peak_ho_hr)
frac_lair_b <- length(dat$states[which(dat$states==1 & is.na(dat$peak_ho_hr) == FALSE)])/nobs_b #what fraction of those non-NA observations were estimated to have belonged to the lair state
frac_emerged_b <- length(dat$states[which(dat$states==2 & is.na(dat$peak_ho_hr) == FALSE)])/nobs_b #what fraction of those non-NA observations were estimated to have belonged to the emerged state
nbins_b <- 13 #number of bins to use in panel b histogram
ht_b <- (2*pi)/nbins_b*nobs_b #height modifier for panel b curves. the y-scale for the curves is density, whereas the y-scale of our histogram 
#is frequency. To convert density to frequency so that the curves can appropriately overlay the histogram, we need to multiply the y-values
#from the density curves above by the area of the histogram. this is calculated as the binwidth of the histogram (2*pi/nbins_b) multiplied by 
#the number of observations in the histogram. We then have to multiply this by the fraction of those observations that the HMM assigned to each
#behavioral state (frac_lair_b or frac_emerged_b)
fcurv_peakhr_lair <- dcurv_peakhr_lair %>% mutate(y=y*ht_b*frac_lair_b)
fcurv_peakhr_emerged <- dcurv_peakhr_emerged %>% mutate(y=y*ht_b*frac_emerged_b)
fcurv_peakhr_sum <- bind_cols(x=xseq_b,y=fcurv_peakhr_lair$y+fcurv_peakhr_emerged$y) #sum of the two frequency curves

fig6b <- dat %>%
  filter(complete.cases(peak_ho_hr))%>%
  ggplot()+
  labs(tag="b")+
  theme(plot.tag.position = c(0.175,1))+
  geom_histogram(aes(x=peak_ho_hr),binwidth = (2*pi/nbins_b),boundary=-pi,closed="right",fill="gray88",col="gray80")+
  scale_x_continuous(limits=c(-pi,pi),expand=c(0,0),breaks=hrs_to_rad(seq(0,24,by=6)),labels=c("12am","6am","12pm","6pm","12am"))+ 
  scale_y_continuous(expand=c(0,0))+
  geom_line(data=fcurv_peakhr_lair, aes(x=x,y=y),col="deepskyblue1",lwd=1)+
  geom_line(data=fcurv_peakhr_emerged, aes(x=x,y=y),col="orange1",lwd=1)+
  geom_line(data=fcurv_peakhr_sum, aes(x=x,y=y),col="black",lwd=1)+
  ylab("Frequency")+
  xlab("Peak haul-out hour")




#___Fig6c: Subadults; effects of covariates on emergence---------------------
DLSeq_sub <- (data.frame(daylength=seq(min(dat$daylength),max(dat$daylength),length=101))) #make a sequence of daylengths
estCI_DL_sub <- vector('list',3)
names(estCI_DL_sub) <- c("est","lower","upper")
estCI_DL_sub[names(estCI_DL_sub)] <- list(list())
for(i in 1:nrow(DLSeq_sub)){ #at each daylength value, estimate emergence probability using the best HMM for subadults
  out <- CIreal(best_mod, covs=data.frame(daylength=DLSeq_sub[i,]))$gamma
  estCI_DL_sub$est[[i]] <- out$est
  estCI_DL_sub$lower[[i]] <- out$lower
  estCI_DL_sub$upper[[i]] <- out$upper
}
par(mfrow=c(1,1))
DL_emerg_sub <- data.frame(DLSeq_sub$daylength,unlist(lapply(estCI_DL_sub$est,function(x) x[1,2])),unlist(lapply(estCI_DL_sub$upper,function(x) x[1,2])),unlist(lapply(estCI_DL_sub$lower,function(x) x[1,2])))
colnames(DL_emerg_sub) <- c("daylength","est","CIupper","CIlower")

fig6c <- DL_emerg_sub %>%
  ggplot()+
  labs(tag="c")+
  theme(plot.tag.position = c(0.2,1),axis.text.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = -20)))+
  geom_ribbon(aes(ymin=CIlower,ymax=CIupper,x=daylength), fill="gray", alpha=0.5) +
  geom_line(aes(x=daylength,y=est),lwd=1)+
  xlab("Daylength (hrs)")+
  ylab("Emergence probability\n")+
  geom_rug(data=dat, aes(x=daylength), col="maroon4", alpha=0.1, size=1.0, sides="b") + 
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0), limits=c(-0.05,1))

grid.arrange(fig6a, fig6b, fig6c, ncol=2) #view
fig6 <- arrangeGrob(fig6a, fig6b, fig6c, ncol=2)
ggsave(fig6, file=here("figures","fig6.png"),width=12,height=8,units="in") #save

#Emergence dates from best model--------------------------
#we have already run "dat$states <- Viterbi(best_mod)" above, which estimates the state sequence for each individual seal
#using that, here we identify when seals switched from the lair state to the emerged state
dat$IDnum <- as.numeric(dat$ID) #make a numeric version of seal ID
splitlist <- split(dat, dat$ID)
emergence <- data.frame(matrix(ncol=7,nrow=length(splitlist)))
colnames(emergence) <- c("IDnum","speno_yr","year","emerg","emerglat","emergtemp","emergDL")
for (i in 1:length(splitlist)) {
  emergence$IDnum[i] <- i
  emergence$year[i] <- mean(dat$year[which(dat$IDnum==i)])
  emergence$speno_yr[i] <- as.character(first(dat$ID[which(dat$IDnum==i)]))
  if (any(dat$doy[which((dat$states > dplyr::lag(dat$states)) & (dat$IDnum == dplyr::lag(dat$IDnum)) & (dat$IDnum == i))]) == TRUE) {
    emergence$emerg[i] <- dat$doy[which((dat$states > dplyr::lag(dat$states)) & (dat$IDnum == dplyr::lag(dat$IDnum)) & (dat$IDnum == i))]
    #if there is a day where a seal is estimated to be in state 2 (emerged) when it was in state 1 (lair) on the previous day, then this the seal's emergence date
    }   else { 
    emergence$emerg[i] <- NA #if not, then the seal doesn't have an emergence date (set emergence date to NA)
  }
   if(is.na(emergence$emerg[i])==FALSE){ #if the seal has an emergence date, then extract...
    emergence$emerglat[i] <- dat$lat[which(dat$doy == emergence$emerg[i] & dat$IDnum == emergence$IDnum[i])] #its latitude on the day that it emerged
    emergence$emergtemp[i] <- dat$air2m[which(dat$doy == emergence$emerg[i] & dat$IDnum == emergence$IDnum[i])] #the air temp on the day that it emerged
    emergence$emergDL[i] <- dat$daylength[which(dat$doy == emergence$emerg[i] & dat$IDnum == emergence$IDnum[i])] #the daylength on the day that it emerged
   }
}
emergence <- emergence %>% filter(complete.cases(emerg)) #remove seals without estimated emergence dates

#summary statistics for emergence date
mean(emergence$emerg) #DOY format
median(emergence$emerg) #DOY format
sd(emergence$emerg) #days
lapply(c(mean(emergence$emerg),median(emergence$emerg)), doytodate) #day month format

#emergence date & daylength
mean(emergence$emergDL,na.rm=TRUE) #on average, what was the daylength when seals emerged
sd(emergence$emergDL,na.rm=TRUE)

#___FigS3------------------
#supplementary figs showing haul-out data & emergence dates for individual seals

dat2 <- dat %>% left_join((emergence[ , c("emerg","speno_yr")] %>% mutate(ID=speno_yr)), by="ID")%>%
  mutate(ID = factor(ID, levels=unique(ID)))
#for drawing purposes, extract the first & last day that the seal was in each state. we'll use these to shade the plot blue for the lair state
#and orange for the emerged state
splitlist <- split(dat2, dat2$ID, drop = TRUE) #sep by speno
for (i in 1:length(splitlist)) {
  splitlist[[i]] <- as.data.frame(splitlist[[i]]) %>%
    mutate(lairmin=min(doy[which(states==1)]),
           lairmax=max(doy[which(states==1)]),
           emergmin=min(doy[which(states==2)]),
           emergmax=max(doy[which(states==2)]))
}
dat2 <- do.call("rbind",splitlist)

dat2$lairmin[which(dat2$lairmin == Inf)] <- NA
dat2$emergmin[which(dat2$emergmin == Inf)] <- NA
dat2$lairmax[which(dat2$lairmax == -Inf)] <- NA
dat2$emergmax[which(dat2$emergmax == -Inf)] <- NA

splitlist <- split(dat2, dat2$ID, drop = TRUE) #separate by unique seal ID

#then, for each of those seals, plot the following...
plots <- lapply(splitlist, function(x) {
  #labels for seal ID, sex, age, & emergence date (if available)
  label_plot <- ggplot(x, aes(x=doy, y=states))+
    ggtitle(paste(x$age_class[1],x$sex[1],x$ID[1],sep=" "),subtitle=paste("Emergence DOY:",x$emerg[1],sep=" "))+ #
    theme_minimal()+
    theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
          axis.text.y=element_blank(),axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),title = element_text(vjust = -2))
  
  #daily proportion hauled out
  prop_ho_plot <- ggplot()+
    annotate("rect", xmin=(as.numeric(paste(x$lairmin[1]))), xmax=(as.numeric(paste(x$lairmax[1]))), ymin=0, ymax=1, alpha=0.2, fill="skyblue1")+
    annotate("rect", xmin=(as.numeric(paste(x$emergmin[1]))), xmax=(as.numeric(paste(x$emergmax[1]))), ymin=0, ymax=1, alpha=0.2, fill="orange")+
    geom_point(data=x, aes(x=doy, y=daily_prop_ho))+
    scale_x_continuous(limits=c(30,170),expand=c(0,0))+
    scale_y_continuous(expand=c(0,0),limits=c(0,1))+
    theme_bw()+
    ylab("Daily proportion hauled out")+
    xlab("DOY")+
    theme(panel.background = element_blank(),  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border=element_blank(),axis.line = element_line(colour = "black"),strip.background =element_blank(),legend.title=element_blank())
  
  #peak haul-out hour
  peak_ho_hr_plot <- ggplot()+
    annotate("rect", xmin=(as.numeric(paste(x$lairmin[1]))), xmax=(as.numeric(paste(x$lairmax[1]))), ymin=-pi, ymax=pi, alpha=0.2, fill="skyblue1")+
    annotate("rect", xmin=(as.numeric(paste(x$emergmin[1]))), xmax=(as.numeric(paste(x$emergmax[1]))), ymin=-pi, ymax=pi, alpha=0.2, fill="orange")+
    geom_point(data=x, aes(x=doy, y=peak_ho_hr))+
    scale_y_continuous(limits=c(-pi,pi),breaks=hrs_to_rad(seq(0,24,by=6)),labels=c("12 am","6 am","12 pm","6 pm","12 am"),expand=c(0,0))+ 
    xlab("DOY")+
    scale_x_continuous(limits=c(30,166),expand=c(0,0))+
    theme_bw()+
    ylab("Peak haul-out hour")+
    theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border=element_blank(),axis.line = element_line(colour = "black"),strip.background =element_blank(),legend.title=element_blank())
  
  arrangeGrob(label_plot,prop_ho_plot,peak_ho_hr_plot,ncol=3)
})

FigS3 <- arrangeGrob(grobs=plots,nrow=length(plots),heights=rep(0.25,length(plots)))
ggsave(here::here("figures","figS3.png"),plot=FigS3,width=12,height=length(plots)*3,limitsize=FALSE )

