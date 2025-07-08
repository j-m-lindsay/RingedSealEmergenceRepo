# Header---------------------------
# Code associated with the paper "Ringed seal (Pusa hispida) haul-out behavior and emergence timing in the Bering, Chukchi, and Beaufort seas"
# Authors: Jessica M. Lindsay, Paul B. Conn, Peter L. Boveng, Justin A. Crawford, Lori T. Quakenbush, Andrew L. Von Duyke, & Kristin L. Laidre
# Contact: jessica.lindsay@noaa.gov
# Date Created: 2025-06-27
# Description: code to compare seasonal patterns in haul-out behavior by sex & age
#
#
# Notes for use --------------------
# The code below uses the here() package and assumes you have a "data" folder folder as indicated in the GitHub repo structure  
# Code was written in R version 4.3.0 "Already Tomorrow" and RStudio version 2023.12.1+402 "Ocean Storm"
#
#Load packages, functions, data-----------------
rm(list = ls()) #clear environment

library(pacman)
p_load(tidyverse, here, ggplot2, circular, gridExtra, ggpubr)

#functions for converting back and forth from hours to radians
hrs_to_rad <- function(hrs){
  return((pi/12)*(hrs-12))
}

rad_to_hrs <- function(rad){
  12*rad/pi+12
}

tags <- read.csv(here("data","ringed_ho_data_clean.csv")) %>%
  mutate_at(c("speno_yr","age_class","sex","data_source"), factor)

#Demographic comparison figures (Figs 3, S1) -----------
theme_set(theme_bw()+theme(panel.background = element_blank(),  panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),panel.border=element_blank(),axis.line = element_line(colour = "black"),
                           axis.text.y = element_text(angle=360, vjust=0.5, hjust=1),text=element_text(size=21), 
                           plot.margin = unit(c(20,20,20,20), "pt")))


#___Fig3: seasonal patterns by age class-------------
#Temporal patterns in haul-out behavior by day of year (DOY) for adult (n = 55) and subadult (n = 16) ringed seals during February–June 
#panel a: daily % hauled out by age
fig3a <- tags %>%
  group_by(doy,age_class,sex)%>%
  summarize(daily_prop_ho_doy=mean(daily_prop_ho,na.rm=TRUE),sd_daily=sd(daily_prop_ho,na.rm=TRUE))%>%
  ggplot()+
  geom_smooth(aes(x=doy,y=daily_prop_ho_doy,group=age_class,fill=age_class,col=age_class),method="loess",alpha=0.5,
              show.legend = FALSE,na.rm=TRUE)+
  scale_color_manual(values=c("darkblue","turquoise"))+
  scale_fill_manual(values=c("darkblue","turquoise"))+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  xlab(label="DOY")+
  ylab(label="Daily % hauled out")

#panels b,c: peak haul-out hour by age
tags2 <- tags %>%
  mutate(age_class_doy=paste0(age_class,doy))%>%
  group_by(age_class_doy) %>%
  summarise(peak_ho_hrs=list(peak_ho_hr),doy=mean(doy),age_class=first(age_class))
tags2$peak_ho_hr <- NA
splitlist <- split(tags2, tags2$age_class_doy, drop = TRUE)
for (i in 1:length(splitlist)){
  splitlist[[i]] <- as.data.frame(splitlist[[i]]) %>%
    arrange(age_class_doy,doy) 
  tags2$peak_ho_hr[[i]] <- mean.circular(unlist(splitlist[[i]]$peak_ho_hrs),na.rm=TRUE)
} #mean.circular will warn you that it is using default values for some arguments, but the defaults are appropriate so OK to disregard
fig3bc <- tags2 %>%
  ggplot()+
  geom_point(aes(x=doy,y=peak_ho_hr,col=age_class),show.legend = FALSE,na.rm=TRUE)+
  scale_color_manual(values=c("darkblue","turquoise"))+
  facet_wrap(vars(age_class))+
  geom_hline(yintercept=0,lty="dashed")+ #dashed line for solar noon
  scale_y_continuous(breaks=hrs_to_rad(seq(0,24,by=6)),labels=c("12 am","6 am","12 pm","6 pm","12 am"))+ 
  xlab("DOY")+
  ylab("Peak hour of day hauled out")+
  theme(panel.spacing = unit(1,"cm"),strip.background = element_blank())

#rm(splitlist,tags2) #removing tags2 as it is not needed outside of creating this figure

grid.arrange(fig3a,fig3bc,widths=c(1,2),    #view
             top = text_grob("b                                         c", y=-1.35,x=0.45,size=21,just="left"),
             bottom = text_grob("a",y=12.8,x=.085,size=21,just="left"))
fig3 <- arrangeGrob(fig3a,fig3bc,widths=c(1,2),      top = text_grob("b                                         c", y=-1.35,x=0.45,size=21,just="left"),
                    bottom = text_grob("a",y=13,x=.085,size=21,just="left"))
ggsave(fig3, file=here("figures","fig3.png"),width=13,height=4.75,units="in") #save

#___FigS1: seasonal patterns by sex----------
#Exploratory plots of seasonal patterns in haul-out behavior by sex for combined age classes (a-c; n = 27 females, n = 44 males) 
#and for adults only (d-f; n = 22 adult females, n = 33 males).

#panel a: daily % hauled out by sex (combined age classes)
figS1a <- tags %>%
  group_by(doy,sex)%>%
  summarize(daily_prop_ho_doy=mean(daily_prop_ho,na.rm=TRUE),sd_daily=sd(daily_prop_ho,na.rm=TRUE))%>%
  ggplot()+
  labs(tag="a")+
  theme(plot.tag.position= c(0.15,1), plot.tag=element_text(size=21))+
  geom_smooth(aes(x=doy,y=daily_prop_ho_doy,group=sex,fill=sex,col=sex),
              method="loess",alpha=0.5,na.rm=TRUE)+
  scale_color_manual(values=c("magenta4","springgreen3"))+
  scale_fill_manual(values=c("magenta4","springgreen3"))+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  xlab(label="DOY")+
  ylab(label="Daily % hauled out")

#panel b,c: peak haul-out hour by sex (combined age classes)
tags3<- tags %>%
  mutate(sex_doy=paste0(sex,doy))%>%
  group_by(sex_doy) %>%
  summarise(peak_ho_hrs=list(peak_ho_hr),doy=mean(doy),sex=first(sex))
tags3$peak_ho_hr <- NA
splitlist <- split(tags3, tags3$sex_doy, drop = TRUE)
for (i in 1:length(splitlist)){
  splitlist[[i]] <- as.data.frame(splitlist[[i]]) %>%
    arrange(sex_doy,doy) 
  tags3$peak_ho_hr[[i]] <- mean.circular(unlist(splitlist[[i]]$peak_ho_hrs),na.rm=TRUE)
}
figS1bc <- tags3 %>%
  ggplot()+
  geom_point(aes(x=doy,y=peak_ho_hr,col=sex),show.legend = FALSE,alpha=0.5,stroke=NA,cex=2,na.rm=TRUE)+
  scale_color_manual(values=c("magenta4","springgreen3"))+
  facet_wrap(vars(sex))+
  geom_hline(yintercept=0,lty="dashed")+ #dashed line for solar noon
  scale_y_continuous(breaks=hrs_to_rad(seq(0,24,by=6)),labels=c("12 am","6 am","12 pm","6 pm","12 am"),limits=c(-pi,pi))+
  xlab("DOY")+
  ylab("Peak hour of day hauled out")+
  theme(panel.spacing = unit(1,"cm"),strip.background = element_blank(),legend.title=element_blank())
rm(splitlist, tags3)

#panel d: daily % hauled out by sex (adults only)
figS1d <- tags %>%
  filter(age_class=="Adult")%>%
  group_by(doy,sex)%>%
  summarize(daily_prop_ho_doy=mean(daily_prop_ho,na.rm=TRUE),sd_daily=sd(daily_prop_ho,na.rm=TRUE))%>%
  ggplot()+
  labs(tag="d")+
  theme(plot.tag.position=c(0.15,1),plot.tag=element_text(size=21))+
  geom_smooth(aes(x=doy,y=daily_prop_ho_doy,group=sex,fill=sex,col=sex),
              method="loess",alpha=0.5,na.rm=TRUE)+
  scale_color_manual(values=c("magenta4","springgreen3"))+
  scale_fill_manual(values=c("magenta4","springgreen3"))+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  xlab(label="DOY")+
  ylab(label="Daily % hauled out")

#panel e,f: peak haul-out hour by sex (adults only)
tags4<- tags %>%
  filter(age_class=="Adult")%>%
  mutate(sex_doy=paste0(sex,doy))%>%
  group_by(sex_doy) %>%
  summarise(peak_ho_hrs=list(peak_ho_hr),doy=mean(doy),sex=first(sex))
tags4$peak_ho_hr <- NA
splitlist <- split(tags4, tags4$sex_doy, drop = TRUE)
for (i in 1:length(splitlist)){
  splitlist[[i]] <- as.data.frame(splitlist[[i]]) %>%
    arrange(sex_doy,doy) 
  tags4$peak_ho_hr[[i]] <- mean.circular(unlist(splitlist[[i]]$peak_ho_hrs),na.rm=TRUE)
}
figS1ef <- tags4 %>%
  mutate(sex=fct_recode(sex, "Adult Female" = "Female")) %>%
  mutate(sex=fct_recode(sex, "Adult Male" = "Male")) %>%
  ggplot()+
  geom_point(aes(x=doy,y=peak_ho_hr,col=sex),show.legend = FALSE,alpha=0.5,stroke=NA,cex=2,na.rm=TRUE)+
  scale_color_manual(values=c("magenta4","springgreen3"))+
  facet_wrap(vars(sex))+
  geom_hline(yintercept=0,lty="dashed")+ #dashed line for solar noon
  scale_y_continuous(breaks=hrs_to_rad(seq(0,24,by=6)),labels=c("12 am","6 am","12 pm","6 pm","12 am"),limits=c(-pi,pi))+
  xlab("DOY")+
  ylab("Peak hour of day hauled out")+
  theme(panel.spacing = unit(1,"cm"),strip.background = element_blank(),legend.title=element_blank())
rm(splitlist, tags4)

grid.arrange(figS1a, figS1bc, figS1d, figS1ef, nrow=2, widths=c(1.5,2),     #view
             top = text_grob("b                                         c", y=-1,x=0.53,size=21,just="left"),
             bottom = text_grob("e                                         f", y=14.5,x=0.53,size=21,just="left")) #view  
figS1 <- arrangeGrob(figS1a, figS1bc, figS1d, figS1ef, nrow=2, widths=c(1.5,2),   
                     top = text_grob("b                                         c", y=-1,x=0.53,size=21,just="left"),
                     bottom = text_grob("e                                         f", y=14.5,x=0.53,size=21,just="left")) 
ggsave(figS1, file=here("figures","figS1.png"),width=14,height=9.5,units="in") #save
#legend for panel d in the manuscript version was altered outside of R to ensure the plot width was the same as panel a
 