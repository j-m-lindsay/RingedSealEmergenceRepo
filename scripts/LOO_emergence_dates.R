# Header---------------------------
# Code associated with the paper "Ringed seal (Pusa hispida) haul-out behavior and emergence timing in the Bering, Chukchi, and Beaufort seas"
# Authors: Jessica M. Lindsay, Paul B. Conn, Peter L. Boveng, Justin A. Crawford, Lori T. Quakenbush, Andrew L. Von Duyke, & Kristin L. Laidre
# Contact: jessica.lindsay@noaa.gov
# Date Created: 2025-06-27 
# Description: code to run leave-one-out sensitivy analysis on ringed seal emergence dates estimated from HMMs and generate figs S4 & S5
#
#
# Notes for use --------------------
# The code below uses the here() package and assumes you have a "data" folder as indicated in the GitHub repo structure 
# Code was written in R version 4.3.0 "Already Tomorrow" and RStudio version 2023.12.1+402 "Ocean Storm"
#

#Load packages, functions, data-----------------
rm(list = ls()) #clear environment

library(pacman)
p_load(tidyverse, momentuHMM, here, ggplot2, RColorBrewer, circular, gridExtra, VGAM)

#functions for converting back and forth from hours to radians
hrs_to_rad <- function(hrs){
  return((pi/12)*(hrs-12))
}

rad_to_hrs <- function(rad){
  12*rad/pi+12
}

tags <- read.csv(here("data","ringed_ho_data_clean.csv")) %>%
  mutate_at(c("speno_yr","age_class","sex","data_source"), factor)

tags_adults <- tags %>% filter(age_class=="Adult")
tags_subs <- tags %>% filter(age_class=="Subadult")

#Leave-one-out sensitivity analysis ------------------------------------
#___Adults-----------
seals_adults <- unique(tags_adults$speno_yr) #list of unique seal IDs
dat_adults <- tags_adults %>%
  mutate(daily_prop_ho=daily_prop_ho/100)%>% #convert percentage to proportion in order to be compatible with beta distribution
  select(c(speno_yr, doy, iceconc, CM_index, EM_index, daily_prop_ho, peak_ho_hr, age_class, sex, year, lat, air2m, lon, daylength, rolling_temp, atdd ))%>%
  rename(ID=speno_yr) #momentuHMM automatically recognizes "ID" as the unique animal identifier

#make an empty data frame to store LOO results in
NA_vec_adults <- rep(NA,length(seals_adults))
LOO_output_adults <- data.frame(seal_removed=NA_vec_adults, emerg_dates=NA_vec_adults,emerg_mean=NA_vec_adults,emerg_med=NA_vec_adults,emerg_sd=NA_vec_adults)

#sequentially remove each seal from dat_LOO, re-fit HMM, re-estimate emergence dates for remaining seals, & store in LOO_output
for (j in 1:length(seals_adults)){
  seal_removed <- seals_adults[j]
  LOO_output_adults$seal_removed[j] <- seal_removed
  dat_LOO <- dat_adults %>%
    filter(ID!=seal_removed)
  
  dat_LOO <- prepData(data=dat_LOO, covNames = c("doy","iceconc","CM_index","EM_index","age_class",
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
  knownStates <- rep(NA,nrow(dat_LOO))
  knownStates[which(dat_LOO$doy < 60)] <- 1
  
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
  

  #no covars
  m_null <- fitHMM(data = dat_LOO, nbStates = 2, dist = dist, Par0 = Par0_m1,
                   estAngleMean = list(peak_ho_hr=TRUE), beta0=NULL, stateNames = stateNames, knownStates = knownStates, fixPar=list(beta=fixedbeta))
  
  
  #daylength+air2m 
  formula_bestmod <- ~ daylength+air2m #formula for transition probabilities, now influenced by covariate
  Par0_bestmod <- getPar0(model=m_null, formula=formula_bestmod)
  best_mod_adults <- fitHMM(data = dat_LOO, nbStates = 2, dist = dist, Par0 = Par0_bestmod$Par,
                              estAngleMean = list(peak_ho_hr=TRUE), stateNames = stateNames, beta0=Par0_bestmod$beta, delta0=Par0_bestmod$delta, formula=formula_bestmod, knownStates = knownStates,fixPar=list(beta=fixedbeta2))
  dat_LOO$states <- viterbi(best_mod_adults)
  
  
  #what are the emergence dates?
  splitlist <- split(dat_LOO, dat_LOO$ID)
  emergence <- data.frame(matrix(ncol=4,nrow=length(splitlist)))
  colnames(emergence) <- c("IDnum","speno_yr","emerg","year")
  dat_LOO$IDnum <- as.numeric(dat_LOO$ID)
  dat_LOO$speno_yr <- dat_LOO$ID
  for (i in 1:length(splitlist)) {
    emergence$IDnum[i] <- i
    emergence$year[i] <- mean(dat_LOO$year[which(dat_LOO$IDnum==i)])
    emergence$speno_yr[i] <- as.character(first(dat_LOO$ID[which(dat_LOO$IDnum==i)]))
    if (any(dat_LOO$doy[which((dat_LOO$states > dplyr::lag(dat_LOO$states)) & (dat_LOO$IDnum == dplyr::lag(dat_LOO$IDnum)) & (dat_LOO$IDnum == i))]) == TRUE) {
      emergence$emerg[i] <- dat_LOO$doy[which((dat_LOO$states > dplyr::lag(dat_LOO$states)) & (dat_LOO$IDnum == dplyr::lag(dat_LOO$IDnum)) & (dat_LOO$IDnum == i))]
    }   else { 
      emergence$emerg[i] <- NA
    }
  }
  
  LOO_output_adults$emerg_dates[j] <- list(emergence)
  LOO_output_adults$emerg_mean[j] <- mean(emergence$emerg,na.rm=TRUE)
  LOO_output_adults$emerg_med[j] <- median(emergence$emerg,na.rm=TRUE)
  LOO_output_adults$emerg_sd[j] <- sd(emergence$emerg,na.rm=TRUE)
 
  print(paste0(j,"/",length(seals_adults)," done")) #print status update to keep track of progress
}


#make a dataframe of all the emergence dates for all the individual seals
all_emergs_adults <- data.frame()
for (i in 1:length(seals_adults)){
  all_emergs_adults <- rbind(all_emergs_adults,data.frame(LOO_output_adults$emerg_dates[[i]]))
}

#take a look at how consistent (low sd_e, small range_e) or inconsistent (high sd_e, large range_e) the emergence dates for individual seals were across LOO iterations
all_emergs_summarized_adults <- all_emergs_adults %>%
  filter(complete.cases(emerg))%>%
  mutate(speno_yr=as.factor(speno_yr))%>%
  group_by(speno_yr)%>%
  summarize(mean_e=mean(emerg,na.rm=TRUE),med_e=median(emerg,na.rm=TRUE),min_e=min(emerg,na.rm=TRUE),max_e=max(emerg,na.rm=TRUE),
            range_e=max_e-min_e,sd_e=sd(emerg,na.rm=TRUE),count_e=sum(!is.na(emerg)),year=mean(year))


#____FigS4 ------------
figS4 <- all_emergs_adults %>%
  arrange(year)%>%
  filter(complete.cases(emerg))%>%
  ggplot()+
  geom_histogram(aes(x=emerg))+
  scale_x_continuous(name="Emergence date (DOY)")+
  facet_wrap(vars(speno_yr),ncol=4)+
  theme_gray()
figS4 #view
ggsave(here::here("figures","figS4.png"),figS4,width=8,height=10,units="in",limitsize=FALSE )

#___Subadults-----------
seals_subs <- unique(tags_subs$speno_yr) #list of unique seal IDs
dat_subs <- tags_subs %>%
  mutate(daily_prop_ho=daily_prop_ho/100)%>% #convert percentage to proportion in order to be compatible with beta distribution
  select(c(speno_yr, doy, iceconc, CM_index, EM_index, daily_prop_ho, peak_ho_hr, age_class, sex, year, lat, air2m, lon, daylength, rolling_temp, atdd ))%>%
  rename(ID=speno_yr) #momentuHMM automatically recognizes "ID" as the unique animal identifier

#make an empty data frame to store LOO results in
NA_vec_subs <- rep(NA,length(seals_subs))
LOO_output_subs <- data.frame(seal_removed=NA_vec_subs, emerg_dates=NA_vec_subs,emerg_mean=NA_vec_subs,emerg_med=NA_vec_subs,emerg_sd=NA_vec_subs)

#sequentially remove each seal from dat_LOO, re-fit HMM, re-estimate emergence dates for remaining seals, & store in LOO_output
for (j in 1:length(seals_subs)){
  seal_removed <- seals_subs[j]
  LOO_output_subs$seal_removed[j] <- seal_removed
  dat_LOO <- dat_subs %>%
    filter(ID!=seal_removed)
  
  dat_LOO <- prepData(data=dat_LOO, covNames = c("doy","iceconc","CM_index","EM_index","age_class",
                                                 "sex","year","lat","air2m","lon","daylength","rolling_temp","atdd"), coordNames=NULL)
  #need to specifiy covNames because momentuHMM automatically interprets any columns that aren't "ID" or in the covNames list as observed variables
  #the covNames list above leaves our observed variables as daily_prop_ho and peak_ho_hr
  
  #specify the probability distribution types for each of the observed variables
  dist <- list(daily_prop_ho="beta",peak_ho_hr="vm") #beta & von mises
  
  #label state 1 "lair" and state 2 "emerged"
  stateNames <- c("lair","emerged")
  
  #for subsubs, PH2007_6021_2007 is the only seal with exact 1s for daily prop HO. when it is this seal's turn 
  #to be removed and there are no exact 1s remaining, HMM gets confused about why there is a one-mass shape param and no
  #1s. So for this one case, need to fit the HMM with only shape1, shape2, and zero-mass params. 
  if (seal_removed=="PH2007_6021_2007"){
    Par0_m1 <- list(daily_prop_ho=c(2,4, #lair shape1, emerged shape1
                                    8,3, #lair shape2, emerged shape2
                                    0.25,0.01), #lair zero-mass, emerged zero-mass
                    peak_ho_hr=c(pi,0, #lair mu, emerged mu
                                 0.5,0.2)) #lair sd, emerged sd
    
  } else if (seal_removed!="PH2007_6021_2007") { 
    Par0_m1 <- list(daily_prop_ho=c(1,8, #lair shape1, emerged shape1
                                    8,4, #lair shape2, emerged shape2
                                    0.75,0.01, #lair zero-mass, emerged zero-mass
                                    0.01,0.1), #lair one-mass, emerged one-mass
                    peak_ho_hr=c(pi,0, #lair mu, emerged mu
                                 0.5,0.2)) #lair sd, emerged sd
  }
  
  #label all data from February (doy <60) as belonging to the lair state (i.e., state 1)
  knownStates <- rep(NA,nrow(dat_LOO))
  knownStates[which(dat_LOO$doy < 60)] <- 1
  
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
  
  
  #no covars
  m_null <- fitHMM(data = dat_LOO, nbStates = 2, dist = dist, Par0 = Par0_m1,
                   estAngleMean = list(peak_ho_hr=TRUE), beta0=NULL, stateNames = stateNames, knownStates = knownStates, fixPar=list(beta=fixedbeta))
  
  #daylength
  formula_DL <- ~ daylength #formula for transition probabilities, now influenced by covariate
  Par0_DL <- getPar0(model=m_null, formula=formula_DL)
  best_mod_subs <- fitHMM(data = dat_LOO, nbStates = 2, dist = dist, Par0 = Par0_DL$Par,
                        estAngleMean = list(peak_ho_hr=TRUE), stateNames = stateNames, beta0=Par0_DL$beta, formula=formula_DL, knownStates = knownStates,fixPar=list(beta=fixedbeta1))
  dat_LOO$states <- viterbi(best_mod_subs)
  
  
  #what are the emergence dates?
  splitlist <- split(dat_LOO, dat_LOO$ID)
  emergence <- data.frame(matrix(ncol=4,nrow=length(splitlist)))
  colnames(emergence) <- c("IDnum","speno_yr","emerg","year")
  dat_LOO$IDnum <- as.numeric(dat_LOO$ID)
  dat_LOO$speno_yr <- dat_LOO$ID
  for (i in 1:length(splitlist)) {
    emergence$IDnum[i] <- i
    emergence$year[i] <- mean(dat_LOO$year[which(dat_LOO$IDnum==i)])
    emergence$speno_yr[i] <- as.character(first(dat_LOO$ID[which(dat_LOO$IDnum==i)]))
    if (any(dat_LOO$doy[which((dat_LOO$states > dplyr::lag(dat_LOO$states)) & (dat_LOO$IDnum == dplyr::lag(dat_LOO$IDnum)) & (dat_LOO$IDnum == i))]) == TRUE) {
      emergence$emerg[i] <- dat_LOO$doy[which((dat_LOO$states > dplyr::lag(dat_LOO$states)) & (dat_LOO$IDnum == dplyr::lag(dat_LOO$IDnum)) & (dat_LOO$IDnum == i))]
    }   else { 
      emergence$emerg[i] <- NA
    }
  }
  
  LOO_output_subs$emerg_dates[j] <- list(emergence)
  LOO_output_subs$emerg_mean[j] <- mean(emergence$emerg,na.rm=TRUE)
  LOO_output_subs$emerg_med[j] <- median(emergence$emerg,na.rm=TRUE)
  LOO_output_subs$emerg_sd[j] <- sd(emergence$emerg,na.rm=TRUE)
  
  print(paste0(j,"/",length(seals_subs)," done")) #print status update to keep track of progress
}


#make a dataframe of all the emergence dates for all the individual seals
all_emergs_subs <- data.frame()
for (i in 1:length(seals_subs)){
  all_emergs_subs <- rbind(all_emergs_subs,data.frame(LOO_output_subs$emerg_dates[[i]]))
}

#take a look at how consistent (low sd_e, small range_e) or inconsistent (high sd_e, large range_e) the emergence dates for individual seals were across LOO iterations
all_emergs_summarized_subs <- all_emergs_subs %>%
  filter(complete.cases(emerg))%>%
  mutate(speno_yr=as.factor(speno_yr))%>%
  group_by(speno_yr)%>%
  summarize(mean_e=mean(emerg,na.rm=TRUE),med_e=median(emerg,na.rm=TRUE),min_e=min(emerg,na.rm=TRUE),max_e=max(emerg,na.rm=TRUE),
            range_e=max_e-min_e,sd_e=sd(emerg,na.rm=TRUE),count_e=sum(!is.na(emerg)),year=mean(year))


#____FigS5 ------------
figS5 <- all_emergs_subs %>%
  arrange(year)%>%
  filter(complete.cases(emerg))%>%
  ggplot()+
  geom_histogram(aes(x=emerg))+
  scale_x_continuous(name="Emergence date (DOY)")+
  facet_wrap(vars(speno_yr),ncol=4)+
  theme_gray()

figS5 #view
ggsave(here("figures","figS5.png"),figS5, width=8, height=3, units="in",limitsize=FALSE ) #save
