# rm(list=ls())

# #### Load libraries
library(MCMCpack)
library(R2jags)
library(data.table)
library(dplyr)
library(ggmcmc)
library(gridExtra)
library(ggthemes)
library(coda)
library(tidyverse)
library(lubridate)

library(lmtest) # granger causality
library(bruceR) # granger causality
library(vars) # granger causality
library(NlinTS) # nonlinear granger causality
library(tseries) # ADF stationarity test

dat <- read_csv(file.path("G:",
                          "Shared drives",
                          "Hansen Lab",
                          "RESEARCH PROJECTS",
                          "MAISRC Mercury Pathways 53.2",
                          "Data",
                          "growth_model_wae.csv"))
dim(dat)

dat <- dat %>% 
  rename('Lat' = `latitude_lake_centroid`,
         'Long' = `longitude_lake_centroid`,
         'ModelReach' = `lake_id`,
         'YEAR' = `year`,
         'length' = length_1) %>% 
mutate(ModelReach = factor(ModelReach))
       
str(dat)

head(dat)
# Number of unique reaches to model (n = 10)
length(unique(dat$ModelReach))

# Number of observations per site per year
n_by_year <- dat %>% 
  count(ModelReach, YEAR, sort = TRUE)
summary(n_by_year$n)

# Plot n by year and site
ggplot(data = n_by_year, mapping = aes(x=YEAR, y=n)) + 
  facet_wrap(~ModelReach) +
  geom_point() + 
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 13)) +
  theme(plot.margin = margin(0, .5, .5, .5, "cm")) +
  labs(title="", y="Number of observations", 
       x="Year") 
# ggsave("figures/abundance_vs_temp.png", height=4, width=8, units="in")
dim(dat)
# Remove unreasonable some length at ages

age_length<-dat%>%count(est_age, length, sort=TRUE)
min(age_length$n)
#there are some samples that may be too large at age 1 to be realistic, might need to remove.
#Need to look into walleye lengths at ages, using the same filter used for Smallmouth Bass
dat <- dat %>% 
  filter(!(est_age==0 & length > 350)) %>% 
  filter(!(est_age==1 & length > 400))

#Checking number of samples for each year
n_by_year <- dat %>% 
  count(ModelReach, YEAR, sort = TRUE)
summary(n_by_year$n) #still more than 20 for each 

# Length at age across all years by site
ggplot(data = dat, mapping = aes(x=est_age, y=length)) + 
  facet_wrap(~ModelReach) +
  geom_point() + 
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 13)) +
  theme(plot.margin = margin(0, .5, .5, .5, "cm")) +
  labs(title="", y="Length (mm)", 
       x="Age (yr)") 
# ggsave("figures/abundance_vs_temp.png", height=4, width=8, units="in")



# Length at age
# US1 <- dat %>% 
#   filter(ModelReach == "US" & Age==1)
# head(US1)
# hist(US1$Length, breaks = 60)
# summary(US1$Length)
# 
# US2 <- dat %>% 
#   filter(ModelReach == "US" & Age==2)
# head(US2)
# hist(US2$Length, breaks = 60)
# summary(US2$Length)
##### !!!!!!!!!!! Rename to match Li et al 2018 model code 
 
# number rivers so we can look back at nhd estimates after running model
# dat$River <- as.numeric(as.factor(as.numeric(dat$ModelReach)))
#write out this data so we have nhd matched to river number for each growth estimate
#write.csv(dat, file = "10yearsRiverNo.csv")

# Rename coumns to match model code
setnames(dat, old=c("ModelReach","YEAR", "est_age", "length"), new=c("River","YR", "SiteAge" ,"SiteFishLength"))
head(dat)

# # Select reaches of interest
# unique(dat$River)
# sites <- c('NB', 'LS', 'US')
# dat <- dat %>% 
#   filter(River %in% sites) %>% 
#   droplevels()

unique(dat$River)

# Name dataset to match code
dd.age <- dat

#only keep these columns
dd.age=dd.age[,c('River','YR','SiteAge','SiteFishLength')]
head(dd.age)


#jags need river to be a numeric
dd.age$River <- as.numeric(as.factor(as.numeric(dd.age$River)))
unique(dd.age$River) #these need to be consecutive order 

head(dd.age)

dim(dd.age)
summary(dd.age)
# dd.age=na.omit(dd.age)
dim(dd.age)

#order years
dd.age <- dd.age[order(dd.age$YR),]

###fixed life history and popn parameters
sim.yr=min(dd.age$YR):max(dd.age$YR)
n.river=length(unique(dd.age$River))
n.age=nrow(dd.age)
n.simyr=length(sim.yr)
n.join=2 #num of para modeled jointly as mutivaiable dist
ww=matrix(c(1,-0.5,-0.5,0.5),nrow=2,ncol=2,byrow=T)


yrlab=dd.age$YR #add a column to label yr in length-age data for model coding
for(y in 1:length(sim.yr)){
  yrlab[which(dd.age$YR==((min(dd.age$YR)-1)+y))]=y
}


Linf.name=rep(NA,n.river*n.simyr)
for(r in 1:n.river){
  for(y in 1:n.simyr){
    Linf.name[(r-1)*n.simyr+y]=paste('LinfYr[',y,',',r,']',sep='')
  }
}


k.name=rep(NA,n.river*n.simyr)
for(r in 1:n.river){
  for(y in 1:n.simyr){
    k.name[(r-1)*n.simyr+y]=paste('kYr[',y,',',r,']',sep='')
  }
}

t0.name=rep(NA,n.river)
for(r in 1:n.river){
  t0.name[r]=paste('t0Yr[',r,']',sep='')
}

sigmaJoin.name=rep(NA,n.join)
for(i in 1:n.join){
  sigmaJoin.name[i]=paste('sigmaJoin[',i,']',sep='')
}

corr.name=rep(NA,sum(1:(n.join-1)))
for(i in 1:(n.join-1)){
  if(i==1){
    aa=0}else{
      aa=1
    }
  for(j in (i+1):n.join){
    corr.name[aa*sum((n.join-i+1):(n.join-1))+(j-i)]=paste('corr[',i,',',j,']',sep='')
  }
}

sigmaJoinYr.name=rep(NA,n.join)
for(i in 1:n.join){
  sigmaJoinYr.name[i]=paste('sigmaJoinYr[',i,']',sep='')
}

corrYr.name=rep(NA,sum(1:(n.join-1)))
for(i in 1:(n.join-1)){
  if(i==1){
    aa=0}else{
      aa=1
    }
  for(j in (i+1):n.join){
    corrYr.name[aa*sum((n.join-i+1):(n.join-1))+(j-i)]=paste('corrYr[',i,',',j,']',sep='')
  }
}

linfk.name=rep(NA,n.river*n.simyr)
for(r in 1:n.river){
  for(y in 1:n.simyr){
    linfk.name[(r-1)*n.simyr+y]=paste('linfk[',y,',',r,']',sep='')
  }
}


######jags model
model=function(){
  #growth	
  for(i in 1:n.age){
    dd.age[i,4]~dlnorm(mu.age[i],tau.age)
    mu.age[i]=log(Linf.ind[i]*(1-exp(-k.ind[i]*(dd.age[i,3]-t0.ind[i]))))
  }	
  tau.age=pow(sigmaAge,-2)
  
  for(i in 1:n.age){
    Linf.ind[i]=LinfYr[yrlab[i],dd.age[i,1]] #yr,river label
    k.ind[i]=kYr[yrlab[i],dd.age[i,1]]
    t0.ind[i]=t0Yr[dd.age[i,1]]
  }
  
  for(r in 1:n.river){
    for(y in 1:n.simyr){
      LinfYr[y,r]=exp(join[y,(r-1)*n.join+1]) #first 2 colns are river 1, 2nd 2 are river 2
      kYr[y,r]=exp(join[y,(r-1)*n.join+2])
      linfk[y,r]=LinfYr[y,r]*kYr[y,r] #track Linf*K
    }
  }
  
  for(r in 1:n.river){
    for (y in 1:(n.simyr-1)){
      join[y+1,((r-1)*n.join+1):((r-1)*n.join+n.join)]~dmnorm(join[y,((r-1)*n.join+1):((r-1)*n.join+n.join)],tau.joinYr[,])
    }
  }
  
  for(r in 1:n.river){
    join[1,((r-1)*n.join+1):((r-1)*n.join+n.join)]~dmnorm(grandmu.join[1:n.join],tau.join[,])
  }
  grandmu.join[1]=log(Linfbar)
  grandmu.join[2]=log(kbar)
  
  for(r in 1:n.river){
    t0Yr[r]~dnorm(t0bar,tau.t0) %_% T(-2,0)
  }
  tau.t0=pow(sigmat0,-2)
  
  tau.join[1:n.join,1:n.join]~dwish(ww[,],df) #first year spatial variation
  df=n.join+1
  varcov[1:n.join,1:n.join]=inverse(tau.join[,])
  for(i in 1:n.join){
    for(j in 1:n.join){
      corr[i,j]=varcov[i,j]/sqrt(varcov[i,i]*varcov[j,j])
    }
    sigmaJoin[i]=sqrt(varcov[i,i])
  }	
  
  tau.joinYr[1:n.join,1:n.join]~dwish(ww[,],df) #each year witin year variation
  varcovYr[1:n.join,1:n.join]=inverse(tau.joinYr[,])
  for(i in 1:n.join){
    for(j in 1:n.join){
      corrYr[i,j]=varcovYr[i,j]/sqrt(varcovYr[i,i]*varcovYr[j,j])
    }
    sigmaJoinYr[i]=sqrt(varcovYr[i,i])
  }	
  
  
  #priors
  Linfbar~dunif(300,875) # grand-mean Linf
  kbar~dunif(0.01,0.5) # Grand-mean k
  t0bar~dunif(-2,0)
  sigmaAge~dunif(0.001,10)
  sigmat0~dunif(0.001,1)
}

datafit=list('dd.age'=dd.age,'n.age'=n.age,'n.river'=n.river,'n.simyr'=n.simyr,'yrlab'=yrlab,'n.join'=n.join,'ww'=ww)
para=c('Linfbar','kbar','t0bar','sigmaAge','sigmat0',Linf.name,k.name,t0.name,sigmaJoin.name,corr.name,sigmaJoinYr.name,corrYr.name,linfk.name)

# Starting values for MCMC
initial=list(list('Linfbar'=800,'kbar'=0.2,'t0bar'=-1.4),
             list('Linfbar'=825,'kbar'=0.2,'t0bar'=-1.5),
             list('Linfbar'=815,'kbar'=0.15,'t0bar'=-1.6))



fit=jags(data=datafit,inits=initial,parameters.to.save=para,model.file=model,
         n.chains=3,n.iter=100000,n.burnin=50000,n.thin=10,DIC=T)


#save fit as rds file-
 saveRDS(fit, file='model.rds')

# Read back into R for further processing
fit <- readRDS(file='model.rds') #test 2 used inital values

# ###############assess convergence
# # Rhat > 1.1
 which(fit$BUGSoutput$summary[, c("Rhat")] > 1.1) 
# # 
# #
# #
# # ####traceplots
 out.mcmc <- as.mcmc(fit)
# # # # Creates full report
S <- ggmcmc::ggs(out.mcmc)
ggmcmc(S) #run this line if you want to see traceplots but takes awhile 
ggmcmc::ggmcmc(S)

# 
# 
# 
# ################################## below is all code to write out txt files of model output, dont need to run this, 
# ####################################just go to Temporal plot K script to plot output
# ###################################extract para values after burn-in and thin
Linf.name.2=rep(NA,n.river*n.simyr)
for(r in 1:n.river){
  for(y in 1:n.simyr){
    Linf.name.2[(r-1)*n.simyr+y]=paste('Linf',r,'_',sim.yr[y],sep='')
  }
}

k.name.2=rep(NA,n.river*n.simyr)
for(r in 1:n.river){
  for(y in 1:n.simyr){
    k.name.2[(r-1)*n.simyr+y]=paste('K',r,'_',sim.yr[y],sep='')
  }
}

t0.name.2=rep(NA,n.river)
for(hh in 1:n.river){
  t0.name.2[hh]=paste('t0',hh,sep='')
}

sigmaJoin.name.2=c('sigmaLinf','sigmaK')
corr.name.2=c('corr')
sigmaJoinYr.name.2=c('sigmaLinfYr','sigmaKYr')
corrYr.name.2=c('corrYr')

linfk.name.2=rep(NA,n.river*n.simyr)
for(r in 1:n.river){
  for(y in 1:n.simyr){
    linfk.name.2[(r-1)*n.simyr+y]=paste('linfk',r,'_',sim.yr[y],sep='')
  }
}

output.array=fit$BUGSoutput$sims.array #[,3,342]
output=rbind(output.array[,1,],output.array[,2,],output.array[,3,]) #convert to matrix; the matrix directly from fit output (fit$BUGSoutput$sims.matrix) does not sort by chain

output=output[,c(para,'deviance')] #re-organize the parameter order in the matrix
colnames(output)=c(para[1:5],Linf.name.2,k.name.2,t0.name.2,sigmaJoin.name.2,corr.name.2,sigmaJoinYr.name.2,corrYr.name.2,linfk.name.2,'deviance')
write.table(output,'estimate_par_VB_walkG_R1.txt',row.names=F,sep='\t',quote=F)


##################summary HPD values for para
para.est=matrix(NA,nrow=ncol(output)+2,ncol=4)
rownames(para.est)=c(colnames(output),'pD','DIC')
colnames(para.est)=c('HPD','lower','upper',"River")

for(i in 1:ncol(output)){
  den=density(output[,i])
  para.est[i,'HPD']=round(den$x[which.max(den$y)],3)
  para.est[i,'upper']=round(quantile(output[,i],probs=0.95),3)
  para.est[i,'lower']=round(quantile(output[,i],probs=0.05),3)

}

para.est['pD',1]=round(fit$BUGSoutput$pD,3)
para.est['DIC',1]=round(fit$BUGSoutput$DIC,3)


write.table(para.est,'summary_par_VB_walkG_R1.txt',row.names=T,sep='\t',quote=F)

########################################
####### Calculate Omega estimates
########################################
dat.out <- fread('estimate_par_VB_walkG_R1.txt')
str(dat.out)
head(dat.out)
dim(dat.out)

# mean(as.matrix(dat.out[,1]))

# Number of rivers * number of years for subsetting growth parameters
total_ests <- n.river * length(sim.yr)
# Number of initial parameters not to use for this
n_skip <- 5

Linfs <- as.matrix(dat.out[,(n_skip+1):(total_ests + n_skip)]) # shift by 6 due to tracking pop average parameters
Ks <- as.matrix(dat.out[,(total_ests + n_skip + 1):(total_ests + n_skip + total_ests)])
# colnames(Linfs)
dim(Linfs)
dim(Ks)

# omega.test1 <- as.numeric(Linfs[,6] * Ks[,6])
# test.den <- density(omega.test1)
# round(test.den$x[which.max(test.den$y)],3)

omegas <- Ks * Linfs
# dim(omegas)
# test.den2 <- density(omegas[,6])
# round(test.den2$x[which.max(test.den2$y)],3)


# New code
# HPD for omegas
##################summary HPD values for para
para.est <- matrix(NA,nrow=ncol(omegas),ncol=4)
colnames(para.est)=c('Omega_HPD','lower','upper','sd')

for(i in 1:ncol(omegas)){
  den=density(omegas[,i])
  para.est[i,'Omega_HPD']=round(den$x[which.max(den$y)],3)
  para.est[i,'upper']=round(quantile(omegas[,i],probs=0.95),3)
  para.est[i,'lower']=round(quantile(omegas[,i],probs=0.05),3)
  para.est[i,'sd']=sd(omegas[,i])
  
}
dim(para.est)
head(para.est)


#Old code#
# HPD for omegas
##################summary HPD values for para
#para.est <- matrix(NA,nrow=ncol(omegas),ncol=3)
#colnames(para.est)=c('Omega_HPD','lower','upper')

#for(i in 1:ncol(omegas)){
 # den=density(omegas[,i])
  #para.est[i,'Omega_HPD']=round(den$x[which.max(den$y)],3)
#  para.est[i,'upper']=round(quantile(omegas[,i],probs=0.95),3)
 # para.est[i,'lower']=round(quantile(omegas[,i],probs=0.05),3)
  
#}
#dim(para.est)
#head(para.est)
# lake <- rep(1:67, 36)

## Ty new code - test, related to site order to match jags output
dat2 <- dat[order(dat$YR),] # Sort by year like in dd.age above
sites2 <- unique(dat2$River) # Compare this site list to the one above using non-sorted data (i.e., sites)
as.numeric(as.factor(unique(dat2$River)))
as.numeric(as.factor(sites2))
#

lake <- vector()
for(i in 1:n.river){
  temp <- rep(i, length=length(sim.yr))
  lake <- append(lake, temp)
}

length(lake)

year <- rep(min(sim.yr):max(sim.yr), n.river)
length(year)

para.est2 <- data.frame(para.est, lake, year)
head(para.est2)
dim(para.est2)

unique(dat$River) 
sites <- sort(unique(dat$River))
as.numeric(as.factor(unique(dat$River)))
as.numeric(as.factor(sites))

#denver adding a lake look up/cross walk
lake_lookup <- data.frame(
  lake = seq_along(sites),
  River = sites
)

para.est2 <- para.est2 %>% 
  left_join(lake_lookup) %>% 
  select(-lake) %>% 
  rename(lake = River)

# para.est2 <- para.est2 %>% 
#   mutate(lake = factor(ifelse(lake==1, "01006200", ifelse(lake==2, '01013700', 
#                         ifelse(lake==3, '01015900', ifelse(lake==4, '03057600',
#                         ifelse(lake==5, '04003000', ifelse(lake==6, '04003501',
#                         ifelse(lake==7, '06000200', ifelse(lake==8, '06015200',
#                         ifelse(lake==9, '10005900', ifelse(lake==10, '11005900', 
#                         ifelse(lake==11,'11014700', ifelse(lake==12,'11016700', 
#                         ifelse(lake==13,'11020300', ifelse(lake==14,'11023400',
#                         ifelse(lake==15,'11030500', ifelse(lake==16,'16038400', 
#                         ifelse(lake==17,'18030800', ifelse(lake==18,'18031000', 
#                         ifelse(lake==19,'18037200', ifelse(lake==20,'18037300',
#                         ifelse(lake==21,'25000100', ifelse(lake==22,'27013300',
#                         ifelse(lake==23,'27017600', ifelse(lake==24,'31082600',
#                         ifelse(lake==25,'32006900', ifelse(lake==26,'33002800',
#                         ifelse(lake==27,'34007900', ifelse(lake==28,'37004600',
#                         ifelse(lake==29,'39000200', ifelse(lake==30,'41011000',
#                         ifelse(lake==31,'47004600', ifelse(lake==32,'47004901',
#                         ifelse(lake==33,'51006300', ifelse(lake==34,'53002800',
#                         ifelse(lake==35,'69037300', ifelse(lake==36,'69037800',
#                         ifelse(lake==37,'69049100', ifelse(lake==38,'69061700',
#                         ifelse(lake==39,'69069300', ifelse(lake==40,'69084500',
#                         ifelse(lake==41,'69129100', ifelse(lake==42,'78002500',
#                         ifelse(lake==43,'83003600', ifelse(lake==44,'86029300','87018000'
#                         ))))))))))))))))))))))))))))))))))))))))))))))

head(para.est2)
write.csv(para.est2, 'omega_estimates.csv', row.names = F)
para.est2 <- read_csv('GrowthModel/omega_estimates.csv')

ggplot() +
  geom_pointrange(data=para.est2,mapping=aes(x=year, y=Omega_HPD, ymin=lower,ymax=upper),
                  position="identity", size=0.1) +
  facet_wrap(~lake) +
  # ylab(expression(paste('Estimated ', omega, ' (mm ',year^-1,')'))) +
  ylab(expression(paste('Early life growth (mm ',year^-1,')'))) +
  xlab('Year') +
  theme(panel.grid.major = element_blank()) +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 11),
        legend.title = element_text(size=11), 
        legend.text=element_text(size=11),strip.text.y = element_text(size=11)) 
ggsave("Growth_trends_1.pdf", height = 8, width = 10, units="in")

#### Adding in lake data and invasion data ####
dim(para.est2)
str(para.est2)


dat #Want lat, long, year_infested, invaded, lake_status added to para.est2 

dat3<-dat %>% 
  select(River,
         YR,
         Lat,
         Long,
         zm_lag,
         exposure)

dat3<-dat3%>%distinct(River,YR, .keep_all=TRUE) #this has the years that are actually in the data and not the estimated years that the omegas have
dat3<-dat3%>%rename(lake=River)


TEST<-left_join(para.est2,dat3, by=join_by(lake==lake, year==YR))
head(TEST)
GM_Omega_Data_FINAL<-TEST
write_csv(GM_Omega_Data_FINAL,"GM_Omega_Data_FINAL.csv")

#checking the dataframes and making sure they line up 
GM_Omega_Data_FINAL%>%filter(lake=='3028700')
dat3%>%filter(lake=='3028700')

GM_Omega_Data_FINAL%>%filter(lake=='47004600')
dat3%>%filter(lake=='47004600')


GM_Omega_Data_FINAL%>%filter(lake=='18037300')
dat3%>%filter(lake=='18037300')

GM_Omega_Data_FINAL%>%filter(lake=='3057600')
dat3%>%filter(lake=='3057600')



#Appear to be lining up in the years that we have data and years that we don't having NAs


#######################################################
#Denver looking at omegas quick
GM_Omega_Data_FINAL %>% 
  filter(!is.na(exposure)) %>% 
  ggplot(aes(x = year, y = Omega_HPD)) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.2
  ) +
  geom_point(aes(color = exposure)) +
  geom_vline(
    data = GM_Omega_Data_FINAL %>% 
      distinct(lake, zm_lag) %>%
      filter(!is.na(zm_lag)),
    aes(xintercept = zm_lag),
    linetype = "dashed"
  ) +
  theme(legend.position = "bottom") +
  facet_wrap(~lake)
ggsave("omega_exposure.jpg", height = 7, width = 11)

GM_Omega_Data_FINAL %>% 
  filter(!is.na(exposure)) %>% 
  ggplot(aes(y = year, x = Omega_HPD)) +
  geom_errorbar(
    aes(xmin = lower, xmax = upper),
    width = 0.2
  ) +
  geom_point(aes(color = exposure)) +
  geom_hline(
    data = GM_Omega_Data_FINAL %>% 
      distinct(lake, zm_lag) %>%
      filter(!is.na(zm_lag)),
    aes(yintercept = zm_lag),
    linetype = "dashed"
  ) +
  facet_wrap(~lake)

#lets explore that poor convergence
bad_rhat <- fit$BUGSoutput$summary[
  fit$BUGSoutput$summary[, "Rhat"] > 1.1,
  ,
  drop = FALSE
]

bad_rhat <- data.frame(
  parameter = rownames(bad_rhat),
  bad_rhat,
  row.names = NULL
)

bad_rhat <- bad_rhat %>%
  extract(
    parameter,
    into = c("parameter", "year_index", "river_index"),
    regex = "^([A-Za-z]+)\\[(\\d+),(\\d+)\\]$",
    convert = TRUE
  )

bad_rhat <- bad_rhat %>%
  mutate(
    YR = sim.yr[year_index],
    River = sites[river_index]
  ) %>% 
  distinct(YR, River, .keep_all = T)

bad_rhat %>% 
  left_join(dat, by = c("River", "YR")) %>% 
  filter(!is.na(River)) %>% 
  ggplot() +
  geom_point(aes(SiteAge, SiteFishLength)) +
  facet_wrap(YR~River, scales = "free")

bad_rhat %>% 
  left_join(dat, by = c("River", "YR")) %>% 
  filter(!is.na(River)) %>% 
  group_by(YR, River) %>% 
  summarise(n = n()) %>% 
  print(n = nrow(.))
ggsave("bad_rhats.jpg", height = 7, width = 11)

dat %>% 
  group_by(YR, River) %>% 
  count() %>% 
  ungroup() %>% 
  summarise(min = min(n),
            median = median(n),
            mean = mean(n),
            max = max(n))

#######################################################
# River temp data
# sus_flow <- read_csv('data/flow/AllsitesSUS.dv.csv')
# head(sus_flow)
# str(sus_flow)
# # Remove empty column
# sus_flow <- sus_flow[,-1]
# 
# sus_flow <- sus_flow %>% 
#   mutate(year = year(Date),
#          month = month(Date),
#          day = day(Date)) %>% 
#   filter(year > 1985 & year < 2023) %>% 
#   select(site_no, Date, Wtemp, Flow, year, month, day) %>% 
#   mutate(site_no = as.numeric(site_no))
# summary(sus_flow$year)
# dim(sus_flow)
# 
# site_info <- read_csv('data/flow/Sites.Coords.csv')
# 
# sus_flow <- sus_flow %>% 
#   left_join(site_info, by="site_no")
# dim(sus_flow)
# str(sus_flow)
# 
# # Grab summer months
# sus_flow <- sus_flow %>% 
#   filter(month > 6 & month < 10)
# summary(sus_flow$month)
# 
# # Grab site for test
# lower_sus <- sus_flow %>% 
#   filter(site_no == 1570500 | 1578310 | 1579550)
# dim(lower_sus)
# 
# # Annual temps
# temps <- lower_sus %>% 
#   group_by(year) %>% 
#   summarize(mean_temp = mean(Wtemp, na.rm = TRUE))
# summary(temps)

# # Read in WQ Portal data
# wqTemps <- read_csv('data/flow/FlowTempSUS.WQP.csv')
# str(wqTemps)
# 
# wqTemps <- wqTemps %>% 
#   filter(CharacteristicName=="Temperature, water") %>% 
#   mutate(year = year(ActivityStartDateTime),
#          month = month(ActivityStartDateTime),
#          day = day(ActivityStartDateTime)) %>% 
#   select(ResultMeasureValue, year, month, day) %>% 
#   filter(month > 6 & month < 10) %>% 
#   filter(year > 1985 & year < 2023) 
# head(wqTemps)
# dim(wqTemps)
# 
# temps2 <- wqTemps %>% 
#   group_by(year) %>% 
#   summarize(mean_temp = mean(ResultMeasureValue, na.rm = TRUE),
#             max_temp = max(ResultMeasureValue, na.rm = TRUE))
# summary(temps2)
# 
# ggplot(temps2, aes(x=year,y=mean_temp)) +
#   geom_point()
# 
# ggplot(temps2, aes(x=year,y=max_temp)) +
#   geom_point()
# 
# # Merge temps and growth
# ls_growth <- para.est2 %>% 
#   filter(river=='LS')
# str(ls_growth)
# 
# ls_growth <- ls_growth %>%
#   left_join(temps2, by="year")
# 
# ggplot(ls_growth, aes(x=mean_temp,y=Omega_HPD)) +
#   geom_point() +
#   geom_smooth()
# 
# ggplot(ls_growth, aes(x=max_temp,y=Omega_HPD)) +
#   geom_point() +
#   geom_smooth()
# 
# # Stdz function
# zscore <- function(x){
#   (x - mean(x))/sd(x)
# }
# 
# # Create standardized version of growth and temperature to visualize
# ls_growth <- ls_growth %>% 
#   mutate(Omega_HPD_z = zscore(Omega_HPD),
#          max_temp_z = zscore(max_temp))
# 
# # Plot standardized values
# graph_dat <- ls_growth %>% 
#   select(Omega_HPD_z, year, max_temp_z) %>% 
#   pivot_longer(cols = c('Omega_HPD_z','max_temp_z'),
#                names_to = "variable", values_to = "value")
# ggplot(graph_dat, aes(x = year, y = value, group=variable)) +
#   geom_line(aes(color = variable), linewidth = 0.7) +
#   scale_color_manual(values = c("#00AFBB", "#E7B800")) +
#   theme_minimal()+
#   labs(title = "")+
#   theme(text = element_text(size = 15))
# 
# 
# #######################################################
# # Causal analysis
# # Granger: The term causality is sometimes missleading, and could rather be thought as "predicatibility".
# # Granger-causality tests were performed in order to understand
# # if the time series x is predictive of the future values of the time series y.
# # granger_test(Omega_HPD ~ mean_temp, data = ls_growth)
# # p1 = ccf_plot(Omega_HPD ~ mean_temp, data = ls_growth)
# # put time-series into ts object
# # https://github.com/nicolarighetti/Granger-causality-test-with-R/blob/main/Granger-test-with-R.md
# # https://lost-stats.github.io/Time_Series/Granger_Causality.html
# 
# # ADF test
# # We cannot reject the null hypothesis because the p-value is not smaller than 0.05.
# # This indicates that the time series is non-stationary. To put it another way, 
# # it has some time-dependent structure and does not exhibit constant variance over time.
# # adf.test(ls_growth$Omega_HPD) # non-stationary
# # adf.test(ls_growth$max_temp) # stationary
# # newY <- diff(log(ls_growth$Omega_HPD), lag=2)
# # adf.test(newY)
# 
# growth <- ts(ls_growth$Omega_HPD)
# temp <- ts(ls_growth$max_temp)
# 
# # Differenced growth data
# # growth <- ts(newY)
# # temp <- ts(ls_growth$max_temp[3:37])
# 
# tsDat <- ts.union(growth, temp) 
# 
# # plot(tsDat)
# # Using var package
# tsVAR <- vars::VAR(tsDat, p = 1)
# c1 <- vars::causality(tsVAR, cause = "temp")$Granger
# c1
# summary(tsVAR)
# 
# # Using lmtest package
# lmtest::grangertest(x=growth, y=temp, order=1)
# lmtest::grangertest(x=temp, y=growth, order=1)
# # We see that the effect of lags of number of X is significant, 
# # and conclude that X predicts the future of Y. The null hypothesis is not
# # rejected for the converse relationship. Thus, we conclude that X Granger causes Y
# 
# # Patrick code
# tsVARselect <- vars::VAR(tsDat, lag.max = 5, ic = "SC")
# summary(tsVARselect)$varresult
# as.numeric(causality(tsVARselect, cause = c("temp"))$Granger[[1]])
# as.numeric(causality(tsVARselect, cause = c("temp"))$Granger[[3]])
# tsVARselect$p # lag chosen based on AIC
# 
# 
# # Nonlinear causality test
# # The null hypothesis of this test is that the second time series does not cause the first one.
# # model = nlin_causality.test (growth, temp, lag=1,LayersUniv = 2, LayersBiv = 1, 5000, 0.01, "sgd", 30, TRUE, 0)
# # model$summary()
# 
# # Impulse response analysis
# # https://www.r-econometrics.com/timeseries/irf/
# # Save the model as a vector autoregression (VAR)
# # model = VAR(tsDat, type = "const", p=1)
# # # Create the impulse response function
# # model_IRF = irf(model, n.ahead = 10, ortho = TRUE, runs = 1000)
# # # Plot 
# # plot(model_IRF)
# 
# ############ CCM
# # # Detrend data
# # ccmdat <- ls_growth %>% 
# #   mutate(Omega_HPD_z_dt = astsa::detrend(Omega_HPD_z) )
# # 
# # # biotic variables to examine 
# # spp <- c('Omega_HPD_z_dt')
# # # physical variables in the study
# # 
# # # physical ariables in the study
# # phys.vars <- c("max_temp_z")
# # 
# # # FUNCTION - Calculate E 
# # simplex_extra_fun<-function(x){
# #   output <- simplex(time_series = x$Omega_HPD_z_dt,
# #                     E=1:10)
# #   output$E[which.max(output$rho)]
# # }
# # 
# # 
# # # FUNCTION = CALCULATE THETA - NONLINEARITY
# # theta_fun <- function(x){
# #   output <- s_map(time_series = x$Omega_HPD_z_dt,
# #                   E = min(x$E))
# #   output$theta[which.max(output$rho)]
# # }
# # 
# # # FUNCTION = CALCULATE NONLINEARITY FROM THETA's MAE
# # # randomization procedure to categorize nonlinearity
# # nonlin_fun <- function(x){
# #   output <- s_map(time_series = x$Omega_HPD_z_dt,
# #                   E = min(x$E))
# #   as.numeric(output$mae)[which(as.numeric(output$theta)==0)] - min(as.numeric(output$mae))
# # }
# # 
# # # CALCULATE E and Theta;
# # dat3 <- ccmdat %>%
# #   mutate(E=simplex_extra_fun(.data)) %>%
# #   mutate(theta=theta_fun(.data)) 
# # 
# # 
# # # original function - to calculate original mae's
# # # Time series were then classified as
# # # nonlinear if the change in mean absolute error (MAE) from a linear to nonlinear
# # # model (ΔMAE = MAEθ=0 − MAEmin) was positive and significant at P ≤ 0.05
# # MAEs <- dat3 %>%
# #   summarize(mae=nonlin_fun(.data))
# # 
# # ###########################################
# # # Calculate null distribution of ΔMAE to compare our original ΔMAE against
# # ###########################################
# # # generate phase-randomized surrogates for assessing significance of nonlinear vs. linear model
# # # container for random datasets
# # n.random.datasets <- 1000
# # ran_series.list <- list()
# #   sur1 <- ccmdat %>% 
# #     select(matches(spp)) %>% 
# #     as.data.frame()
# #   temp.series <- surrogates(sur1,nsim=n.random.datasets)
# #   ran_series.list <- as.data.frame(temp.series)
# # 
# # 
# # # Container for ΔMAE
# # null_maes <- vector(length=n.random.datasets)
# #   for(j in 1:n.random.datasets){
# #     # Simplex to get E
# #     t_simp <- simplex(time_series = ran_series.list[,j],
# #                       E=1:10)
# #     E <- t_simp$E[which.max(t_simp$rho)]
# #     # s_map to get ΔMAE
# #     t_s_map <- s_map(time_series = ran_series.list[,j],
# #                      E = E)
# #     null_maes[j] <- as.numeric(t_s_map$mae)[which(as.numeric(t_s_map$theta)==0)] - min(as.numeric(t_s_map$mae))
# #     
# #   }
# # 
# # # Get 95th %-tile of null distribution of ΔMAEs
# # quant <- quantile(null_maes, 0.95)
# # # Merge with estimated MAEs
# # MAEsSig <- data.frame(MAEs, quant)
# # # Categorize lakes as linear or nonlinear based on significance of ΔMAE compared to null dist'n of ΔMAEs
# # linear_nonlinear <- MAEsSig %>%
# #   mutate(nonlin=ifelse(
# #     mae > quant,
# #     "nonlinear",
# #     "linear"))
# # 
# # 
# # ###################################
# # # Calculate at predictability 
# # ###################################
# # # FUNCTION: CALCULATE Tp - FORECAST SKILL
# # 
# # forecast_fun <- function(x) {
# #   output <- s_map(time_series = x$Omega_HPD_z_dt,
# #                   E = min(x$E),
# #                   theta = min(x$theta),
# #                   tp = 1)
# #   
# #   as.numeric(output$rho)[which(output$tp==1)]
# # }
# # 
# # 
# # #FUNCTION: Calculate p values for forecast skill
# # pred_fun <- function(x) {
# #   output <- s_map(time_series = x$Omega_HPD_z_dt,
# #                   E = min(x$E),
# #                   theta = min(x$theta),
# #                   tp = 1)
# #   
# #   output$p_val[which(output$tp==1)]
# # }
# # 
# # 
# # # Calculate predictability and p-value
# # dat4 <- dat3 %>%
# #   mutate(rho=forecast_fun(.data)) %>%
# #   mutate(p=pred_fun(.data))
# # 
# # # Get one record per lake for merging with significance of nonlinearity
# # dat5 <- dat4 %>%
# #   filter(row_number()==1)
# # 
# # dat5 <- cbind(dat5, linear_nonlinear)
# # 
# # dat5 <- dat5 %>% 
# #   mutate(nonlin = factor(nonlin)) %>% 
# #   mutate(significant = factor(ifelse(p < 0.05, "sig", "not sig")) )
