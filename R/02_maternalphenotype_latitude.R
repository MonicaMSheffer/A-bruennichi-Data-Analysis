#this script uses our data on maternal and phenotypes (body size, pigmentation, fecundity, hatching success, offspring body size) to test for relationships between the measured traits and latitude

#load libraries
library(ggplot2)
library(glmmTMB)
library(bbmle)
library(DHARMa)
library(sjPlot)
library(ggeffects)
library(performance)
library(MuMIn)
set_theme(theme_bw())

#our starting data should have been read in in the 01_import_data.R script
#we will use the following data frames:
#View(ddrad_phenotypesMothers)
#View(ddrad_phenotypesOffspring)


#prepare pigmentation data for binary data analysis
temp_pig=cbind(dark=ddrad_phenotypesMothers$dark,light=ddrad_phenotypesMothers$light)
ddrad_phenotypesMothers$Pig_binary=temp_pig
ddrad_phenotypesMothers$allPixels=ddrad_phenotypesMothers$light+ddrad_phenotypesMothers$dark
ddrad_phenotypesMothers$dark.light_ratio=ddrad_phenotypesMothers$dark/ddrad_phenotypesMothers$allPixels

#get sample sizes (may be different, some photos had to be thrown out due to poor positioning or other issues; all should have leg lengths measured)
length(na.omit(ddrad_phenotypesMothers$dark.light_ratio))
length(na.omit(ddrad_phenotypesMothers$LegLengthM))
length(na.omit(ddrad_phenotypesMothers$PC1))

plot(ddrad_phenotypesMothers$PC1~ddrad_phenotypesMothers$GPS_N,xlab="GPS °N",ylab="PC1",pch=19)
text(PC1 ~GPS_N, labels=population_name,data=ddrad_phenotypesMothers, cex=0.6,pos=4)

#fix dates, calculate oviposition days
ddrad_phenotypesMothers$OvipoistionDate=as.Date(ddrad_phenotypesMothers$OvipoistionDate,format="%m/%d/%Y",origin="12/30/1899")
ddrad_phenotypesOffspring$OvipoistionDate=as.Date(ddrad_phenotypesOffspring$OvipoistionDate,format="%m/%d/%Y",origin="12/30/1899")
ddrad_phenotypesMothers$CollectingDate=as.Date(ddrad_phenotypesMothers$CollectingDate,format="%m/%d/%Y",origin="12/30/1899")
ddrad_phenotypesOffspring$CollectingDate=as.Date(ddrad_phenotypesOffspring$CollectingDate,format="%m/%d/%Y",origin="12/30/1899")

ddrad_phenotypesMothers$OviDays=ddrad_phenotypesMothers$OvipoistionDate-ddrad_phenotypesMothers$CollectingDate
ddrad_phenotypesOffspring$OviDays=ddrad_phenotypesOffspring$OvipoistionDate-ddrad_phenotypesOffspring$CollectingDate

#get rid of duplicate columns, only keep necessary ones for plotting
mother_toPlot=ddrad_phenotypesMothers[,c(4:6,10:65,67:80)]
for(i in 1:length(mother_toPlot$GPS_E)){
  if(mother_toPlot$GPS_E[i]<11.2) {
    mother_toPlot$Group[i]="southwestern"
  }  else {
    mother_toPlot$Group[i]="northeastern"
  } 
}
levels(as.factor(mother_toPlot$population_name))

#visualize data to check for any weirdness
ggplot(data=mother_toPlot,aes(x=GPS_N,y=PC1,color=Group))+
  geom_point(aes(shape=Country),size=2,alpha=.5)+
  scale_color_manual(values=c("#0082c8","#f58231","grey36"))+
  scale_shape_manual(values=c(17,19,15))+
  theme_bw()+
  xlab("Latitude (°N)")

ggplot(data=mother_toPlot,aes(x=GPS_N,y=Total,color=Group))+
  geom_point(aes(shape=Country),size=2,alpha=.5)+
  scale_color_manual(values=c("#0082c8","#f58231","grey36"))+
  scale_shape_manual(values=c(17,19,15))+
  theme_bw()+
  xlab("Latitude (°N)")

ggplot(data=mother_toPlot,aes(x=GPS_N,y=LegLengthM,color=Group))+
  geom_point(aes(shape=Country),size=2,alpha=.5)+
  scale_color_manual(values=c("#0082c8","#f58231","grey36"))+
  scale_shape_manual(values=c(17,19,15))+
  theme_bw()+
  xlab("Latitude (°N)")


##########lmm for adult female leg length##########
hist(mother_toPlot$LegLengthM)
m0_motherLL_lat=glmmTMB(LegLengthM~(1|population_name),data=mother_toPlot,family=gaussian())
m1_motherLL_lat=glmmTMB(LegLengthM~GPS_N+(1|population_name),data=mother_toPlot,family=gaussian())
AICctab(m0_motherLL_lat,m1_motherLL_lat)

#check residuals using DHARMa
hist(resid(m1_motherLL_lat))
M1gauss_simulationOutput <- simulateResiduals(m1_motherLL_lat, n = 500)
testOutliers(m1_motherLL_lat,type="bootstrap",nBoot=100)
testResiduals(M1gauss_simulationOutput)
testDispersion(M1gauss_simulationOutput)
plot(M1gauss_simulationOutput) 
testQuantiles(M1gauss_simulationOutput) 
#All are acceptable

#get model summary stats
performance::r2(m1_motherLL_lat)
length(na.omit(mother_toPlot$LegLengthM))
confint(m1_motherLL_lat)
summary(m1_motherLL_lat)

#initial plots of model estimates
plot_model(m1_motherLL_lat,type="re")
plot_model(m1_motherLL_lat,type="est",sort.est = F,show.values=T,p.adjust="fdr")
plot_model(m1_motherLL_lat,type="pred",terms="GPS_N")

pred_motherLL_lat=ggpredict(m1_motherLL_lat,terms="GPS_N")
plot(pred_motherLL_lat,add.data=T)

#create plot for figure in manuscript
mother_LL_plot_lat=ggplot(dat=pred_motherLL_lat)+
  geom_point(dat=mother_toPlot,aes(x=GPS_N,y=LegLengthM/1000,color=Group,shape=Country),alpha=0.6)+
  scale_color_manual(values=c("#0082c8","#f58231"))+
  scale_shape_manual(values=c(17,16,15))+
  geom_ribbon(aes(x=x,ymin=conf.low/1000,ymax=conf.high/1000),alpha=0.25)+
  geom_line(aes(x=x,y=predicted/1000),linewidth=.75)+
  xlab("Latitude (°N)")+
  ylab("Adult female leg length (mm)")+
  theme_bw()+
  xlim(43,58)+
  theme(axis.text.x=element_text(size=6),axis.text.y=element_text(size=6))+
  theme(axis.title=element_text(size=7))+
  theme(axis.title.y = element_text(margin = margin(r = 4)),axis.title.x = element_text(margin = margin(t = 4)))+ 
  theme(legend.position="none")
mother_LL_plot_lat
pdf("output/Figures/motherLegLength_final.pdf",width=2.5,height=2.5)
mother_LL_plot_lat
dev.off()
#####################glmm for adult female pigmentation########################
hist(mother_toPlot$dark.light_ratio)
m1_motherPig_betabin=glmmTMB(Pig_binary~GPS_N+(1|population_name),data=mother_toPlot,family=betabinomial())
m1_motherPig_binLogit=glmmTMB(Pig_binary~GPS_N+(1|population_name),data=mother_toPlot,family=binomial(link="logit"))
m1_motherPig_binProbit=glmmTMB(Pig_binary~GPS_N+(1|population_name),data=mother_toPlot,family=binomial(link="probit"))
AICctab(m1_motherPig_betabin,m1_motherPig_binLogit,m1_motherPig_binProbit)#betabinomial is the best by far

m0_motherPig_lat=glmmTMB(Pig_binary~(1|population_name),data=mother_toPlot,family=betabinomial())
m1_motherPig_lat=glmmTMB(Pig_binary~GPS_N+(1|population_name),data=mother_toPlot,family=betabinomial()) 
AICctab(m0_motherPig_lat,m1_motherPig_lat)#null model fits slightly better, almost no difference between the two

hist(resid(m1_motherPig_lat))
m1_pig_simulationOutput <- simulateResiduals(m1_motherPig_lat, n = 500)
testOutliers(m1_pig_simulationOutput,type="bootstrap",nBoot=100)
testResiduals(m1_pig_simulationOutput)
testDispersion(m1_pig_simulationOutput)
plot(m1_pig_simulationOutput) 
testQuantiles(m1_pig_simulationOutput)#all acceptable

performance::r2_nakagawa(m1_motherPig_lat,tolerance=1e-5)#can't reliably calculate for betabinomial models; use rough calculation

pig_response=apply(get.response(model.frame(m1_motherPig_betabin)),1,function(x)x[1]/sum(x))
cor(pig_response,predict(m1_motherPig_betabin,type="response"))^2


length(na.omit(mother_toPlot$Pig_binary))/2
confint(m1_motherPig_lat)
summary(m1_motherPig_lat)

plot_model(m1_motherPig_lat,type="re")
plot_model(m1_motherPig_lat,type="est",sort.est = F,show.values=T,p.adjust="fdr")
plot_model(m1_motherPig_lat,type="pred",terms="GPS_N")

pred_motherpig_lat=ggpredict(m1_motherPig_lat,terms="GPS_N")

mother_pig_plot_lat=ggplot(dat=pred_motherpig_lat)+
  geom_point(dat=mother_toPlot,aes(x=GPS_N,y=dark.light_ratio,color=Group,shape=Country),alpha=0.6)+
  scale_color_manual(values=c("#0082c8","#f58231"))+
  scale_shape_manual(values=c(17,16,15))+
  geom_ribbon(aes(x=x,ymin=conf.low,ymax=conf.high),alpha=0.25)+
  geom_line(aes(x=x,y=predicted),linewidth=.75)+
  xlab("Latitude (°N)")+
  ylab("Proportion of dark pigmentation")+
  theme_bw()+
  xlim(43,58)+
  theme(axis.text.x=element_text(size=6),axis.text.y=element_text(size=6))+
  theme(axis.title=element_text(size=7))+
  theme(axis.title.y=element_text(margin=margin(r=4)),axis.title.x = element_text(margin = margin(t = 4)))+ 
  theme(legend.position="none")
mother_pig_plot_lat

pdf("output/Figures/motherPigmentation_final.pdf",width=2.5,height=2.5)
mother_pig_plot_lat
dev.off()
###########glmm for adult female fecundity##########
hist(mother_toPlot$Total)

#center and scale predictor variables
mother_toPlot$CS_lat=scale(mother_toPlot$GPS_N,center=T,scale=T)
mother_toPlot$CS_LegLengthM=scale(mother_toPlot$LegLengthM,center=T,scale=T)
mother_toPlot$CS_OviDays=scale(mother_toPlot$OviDays,center=T,scale=T)

m0_motherFecund_lat=glmmTMB(Total~(1|population_name),data=mother_toPlot,family=poisson())
m1_motherFecund_lat=glmmTMB(Total~CS_lat+(1|population_name),data=mother_toPlot,family=poisson())
m2_motherFecund_lat=glmmTMB(Total~CS_LegLengthM+(1|population_name),data=mother_toPlot,family=poisson())
m3_motherFecund_lat=glmmTMB(Total~CS_OviDays+(1|population_name),data=mother_toPlot,family=poisson())
m4_motherFecund_lat=glmmTMB(Total~CS_lat+CS_LegLengthM+(1|population_name),data=mother_toPlot,family=poisson())
m5_motherFecund_lat=glmmTMB(Total~CS_lat+CS_OviDays+(1|population_name),data=mother_toPlot,family=poisson())
m6_motherFecund_lat=glmmTMB(Total~CS_LegLengthM+CS_OviDays+(1|population_name),data=mother_toPlot,family=poisson())
m7_motherFecund_lat=glmmTMB(Total~CS_lat+CS_LegLengthM+CS_OviDays+(1|population_name),data=mother_toPlot,family=poisson(link="log"))
m8_motherFecund_lat=glmmTMB(Total~CS_lat+CS_LegLengthM+CS_OviDays+(1|population_name),data=mother_toPlot,family=nbinom1())
m9_motherFecund_lat=glmmTMB(Total~CS_lat+CS_LegLengthM+CS_OviDays+(1|population_name),data=mother_toPlot,family=nbinom2())
AICctab(m0_motherFecund_lat,m1_motherFecund_lat,m2_motherFecund_lat,m3_motherFecund_lat,m4_motherFecund_lat,m5_motherFecund_lat,m6_motherFecund_lat,m7_motherFecund_lat,m8_motherFecund_lat,m9_motherFecund_lat)#negative binomial model with latitude, leg length, and oviposition latency has lowest AICc

hist(resid(m9_motherFecund_lat))
M9_fecund_simulationOutput <- simulateResiduals(m9_motherFecund_lat, n = 500)
testOutliers(M9_fecund_simulationOutput,type="bootstrap",nBoot=100)
testResiduals(M9_fecund_simulationOutput)
testDispersion(M9_fecund_simulationOutput)
plot(M9_fecund_simulationOutput) 
testQuantiles(M9_fecund_simulationOutput)

performance::r2(m9_motherFecund_lat)
length(na.omit(ddrad_phenotypesMothers$Total))
confint(m9_motherFecund_lat)
summary(m9_motherFecund_lat)

plot_model(m9_motherFecund_lat,type="re")
plot_model(m9_motherFecund_lat,type="pred",terms=c("CS_lat"))
plot_model(m9_motherFecund_lat,type="pred",terms=c("CS_OviDays"))
plot_model(m9_motherFecund_lat,type="pred",terms=c("CS_LegLengthM"))

pred_fecund_lat=ggpredict(m9_motherFecund_lat,terms="CS_lat")
pred_fecund_lat$x.orig=pred_fecund_lat$x*attr(mother_toPlot$CS_lat, 'scaled:scale')+attr(mother_toPlot$CS_lat,'scaled:center')

mother_fecund_plot_lat=ggplot(data=pred_fecund_lat,aes(x=x.orig))+
  geom_point(data=mother_toPlot,aes(x=GPS_N,y=Total,color=Group,shape=Country),alpha=0.6)+
  scale_color_manual(values=c("#0082c8","#f58231"))+
  scale_shape_manual(values=c(17,16,15))+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.25)+
  geom_line(aes(y=predicted),linewidth=.75)+
  xlab("Latitude (°N)")+
  ylab("Clutch size")+
  theme_bw()+
  xlim(43,58)+
  theme(axis.text.x=element_text(size=6),axis.text.y=element_text(size=6))+
  theme(axis.title=element_text(size=7))+
  theme(axis.title.y = element_text(margin = margin(r = 4)),axis.title.x = element_text(margin = margin(t = 4)))+ 
  theme(legend.position="none")
mother_fecund_plot_lat

pdf("output/Figures/motherFecundity_final.pdf",width=2.5,height=2.5)
mother_fecund_plot_lat
dev.off()
##################glmm for hatching success #############
temp_HS=cbind(hatched=mother_toPlot$Dead_Alive,unhatched=mother_toPlot$Eggs)
mother_toPlot$HatchingSuccess=temp_HS
plot(((Dead+Alive)/Total)~GPS_N,data=mother_toPlot)

hist(mother_toPlot$HatchingSuccess) #likely zero-inflated
mother_toPlot$HatchProportion<-(mother_toPlot$Dead+mother_toPlot$Alive)/mother_toPlot$Total
hist(mother_toPlot$HatchProportion)
mean(mother_toPlot$HatchProportion)
sd(mother_toPlot$HatchProportion)

m1_hatch_binLog=glmmTMB(HatchingSuccess~CS_lat+CS_LegLengthM+CS_OviDays+(1|population_name),data=mother_toPlot,family=binomial(link="logit"))
m1_hatch_binProb=glmmTMB(HatchingSuccess~CS_lat+CS_LegLengthM+CS_OviDays+(1|population_name),data=mother_toPlot,family=binomial(link="probit"))
m1_hatch_betabin=glmmTMB(HatchingSuccess~CS_lat+CS_LegLengthM+CS_OviDays+(1|population_name),data=mother_toPlot,family=betabinomial())
AICctab(m1_hatch_binLog,m1_hatch_binProb,m1_hatch_betabin)

m0_hatch_lat=glmmTMB(HatchingSuccess~(1|population_name),data=mother_toPlot,family=betabinomial())
m1_hatch_lat=glmmTMB(HatchingSuccess~CS_lat+(1|population_name),data=mother_toPlot,family=betabinomial())
m2_hatch_lat=glmmTMB(HatchingSuccess~CS_LegLengthM+(1|population_name),data=mother_toPlot,family=betabinomial())
m3_hatch_lat=glmmTMB(HatchingSuccess~CS_OviDays+(1|population_name),data=mother_toPlot,family=betabinomial())
m4_hatch_lat=glmmTMB(HatchingSuccess~CS_lat+CS_LegLengthM+(1|population_name),data=mother_toPlot,family=betabinomial())
m5_hatch_lat=glmmTMB(HatchingSuccess~CS_lat+CS_OviDays+(1|population_name),data=mother_toPlot,family=betabinomial())
m6_hatch_lat=glmmTMB(HatchingSuccess~CS_LegLengthM+CS_OviDays+(1|population_name),data=mother_toPlot,family=betabinomial())
m7_hatch_lat=glmmTMB(HatchingSuccess~CS_lat+CS_LegLengthM+CS_OviDays+(1|population_name),data=mother_toPlot,family=betabinomial())
AICctab(m0_hatch_lat,m1_hatch_lat,m2_hatch_lat,m3_hatch_lat,m4_hatch_lat,m5_hatch_lat,m6_hatch_lat,m7_hatch_lat)

hist(resid(m7_hatch_lat))
m7_hatch_simulationOutput <- simulateResiduals(m7_hatch_lat, n = 500)
testOutliers(m7_hatch_simulationOutput,type="bootstrap",nBoot=100)
testResiduals(m7_hatch_simulationOutput)
testDispersion(m7_hatch_simulationOutput)
plot(m7_hatch_simulationOutput) 
testQuantiles(m7_hatch_simulationOutput)
testZeroInflation(m7_hatch_simulationOutput)
check_singularity(m7_hatch_lat)

plot_model(m7_hatch_lat,type="re")
plot_model(m7_hatch_lat,type="pred",terms=c("CS_lat"))
plot_model(m7_hatch_lat,type="pred",terms=c("CS_OviDays"))
plot_model(m7_hatch_lat,type="pred",terms=c("CS_LegLengthM"))

performance::r2(m7_hatch_lat)

HS_response=apply(get.response(model.frame(m7_hatch_lat)),1,function(x)x[1]/sum(x))
cor(HS_response,predict(m7_hatch_lat,type="response"))^2
length(na.omit(mother_toPlot$HatchingSuccess))/2
confint(m7_hatch_lat)
summary(m7_hatch_lat)

pred_hatch_lat=ggpredict(m7_hatch_lat,terms="CS_lat")
pred_hatch_lat$x.orig=pred_hatch_lat$x*attr(mother_toPlot$CS_lat, 'scaled:scale')+attr(mother_toPlot$CS_lat,'scaled:center')

mother_hatch_plot_lat=ggplot(data=pred_hatch_lat,aes(x=x.orig))+
  geom_point(data=mother_toPlot,aes(x=GPS_N,y=((Dead+Alive)/Total),color=Group,shape=Country),alpha=0.6)+
  scale_color_manual(values=c("#0082c8","#f58231"))+
  scale_shape_manual(values=c(17,16,15))+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.25)+
  geom_line(aes(y=predicted),linewidth=.75)+
  xlab("Latitude (°N)")+
  ylab("Hatching proportion")+
  theme_bw()+
  xlim(43,58)+
  theme(axis.text.x=element_text(size=6),axis.text.y=element_text(size=6))+
  theme(axis.title=element_text(size=7))+
  theme(axis.title.y=element_text(margin=margin(r=4)),axis.title.x=element_text(margin=margin(t = 4)))+ 
  theme(legend.position="none")
mother_hatch_plot_lat
pdf("output/Figures/motherHatchingSuccess_final.pdf",width=2.5,height=2.5)
mother_hatch_plot_lat
dev.off()

##########lmm for offspring leg length##########
#get rid of duplicate columns, only keep necessary ones for modeling and plotting
offspring_toPlot=ddrad_phenotypesOffspring[,c(1,3:6,10:13,18,21,30,77)]
for(i in 1:length(offspring_toPlot$GPS_E)){
  if(offspring_toPlot$GPS_E[i]<11.2) {
    offspring_toPlot$Group[i]="southwestern"
  }  else {
    offspring_toPlot$Group[i]="northeastern"
  } 
}
hist(offspring_toPlot$LegLengthS)
mean(offspring_toPlot$LegLengthS)
sd(offspring_toPlot$LegLengthS)
#center and scale predictor variables
offspring_toPlot$CS_lat=scale(offspring_toPlot$GPS_N,center=T,scale=T)
offspring_toPlot$CS_LegLengthM=scale(offspring_toPlot$LegLengthM,center=T,scale=T)
offspring_toPlot$CS_OviDays=scale(offspring_toPlot$OviDays,center=T,scale=T)
offspring_toPlot$CS_famFecundity=scale(offspring_toPlot$Total,center=T,scale=T)

m0_offspringLL_lat=glmmTMB(LegLengthS~(1|population_name/Mother),data=offspring_toPlot,family=gaussian())
m7_offspringLL_lat=glmmTMB(LegLengthS~CS_lat+CS_LegLengthM+CS_OviDays+CS_famFecundity+(1|population_name/Mother),data=offspring_toPlot,family=gaussian())
AICctab(m0_offspringLL_lat,m7_offspringLL_lat)#null model actually better

hist(resid(m7_offspringLL_lat))
m7_offspringLL_simulationOutput <- simulateResiduals(m7_offspringLL_lat, n = 500)
testOutliers(m7_offspringLL_simulationOutput,type="bootstrap",nBoot=100)
testResiduals(m7_offspringLL_simulationOutput)
testDispersion(m7_offspringLL_simulationOutput)
plot(m7_offspringLL_simulationOutput) 
testQuantiles(m7_offspringLL_simulationOutput)
check_singularity(m7_offspringLL_lat)

performance::r2(m7_offspringLL_lat)
length(na.omit(offspring_toPlot$LegLengthS))
confint(m7_offspringLL_lat)
summary(m7_offspringLL_lat)

plot_model(m7_offspringLL_lat,type="re")
plot_model(m7_offspringLL_lat,type="pred",terms=c("CS_lat"))
plot_model(m7_offspringLL_lat,type="pred",terms=c("CS_OviDays"))
plot_model(m7_offspringLL_lat,type="pred",terms=c("CS_LegLengthM"))
plot_model(m7_offspringLL_lat,type="pred",terms=c("CS_famFecundity"))

pred_offspringLL_lat=ggpredict(m7_offspringLL_lat,terms="CS_lat")
pred_offspringLL_lat$x.orig=pred_offspringLL_lat$x*attr(offspring_toPlot$CS_lat, 'scaled:scale')+attr(offspring_toPlot$CS_lat,'scaled:center')

offspringLL_plot_lat=ggplot(data=pred_offspringLL_lat,aes(x=x.orig))+
  geom_point(data=offspring_toPlot,aes(x=GPS_N,y=LegLengthS,color=Group,shape=Country),alpha=0.3)+
  scale_color_manual(values=c("#0082c8","#f58231"))+
  scale_shape_manual(values=c(17,16,15))+
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.25)+
  geom_line(aes(y=predicted),linewidth=1.2)+
  xlab("Latitude (°N)")+
  ylab("Offspring leg length (µm)")+
  theme_bw()+
  xlim(43,58)+
  theme(axis.text.x=element_text(size=6),axis.text.y=element_text(size=6))+
  theme(axis.title=element_text(size=7))+
  theme(axis.title.y=element_text(margin=margin(r=4)),axis.title.x=element_text(margin=margin(t = 4)))+ 
  theme(legend.position="none")
offspringLL_plot_lat
pdf("output/Figures/offspringBodySize_final.pdf",width=2.5,height=2.5)
offspringLL_plot_lat
dev.off()
##################glmm for oviposition latency #############
mother_toPlot$OviDays=as.numeric(mother_toPlot$OviDays)
hist(mother_toPlot$OviDays)
range(mother_toPlot$OviDays)
plot(OviDays~GPS_N,data=mother_toPlot)

#center and scale predictor variables
mother_toPlot$CS_lat=scale(mother_toPlot$GPS_N,center=T,scale=T)
mother_toPlot$CS_LegLengthM=scale(mother_toPlot$LegLengthM,center=T,scale=T)
mother_toPlot$OviDays=as.numeric(mother_toPlot$OviDays)

m1_oviDays_gaus=glmmTMB(OviDays~CS_lat+CS_LegLengthM+(1|population_name),data=mother_toPlot,family=gaussian())
m1_oviDays_pois=glmmTMB(OviDays~CS_lat+CS_LegLengthM+(1|population_name),data=mother_toPlot,family=poisson())
m1_oviDays_nb1=glmmTMB(OviDays~CS_lat+CS_LegLengthM+(1|population_name),data=mother_toPlot,family=nbinom1())
m1_oviDays_nb2=glmmTMB(OviDays~CS_lat+CS_LegLengthM+(1|population_name),data=mother_toPlot,family=nbinom2())

AICctab(m1_oviDays_gaus,m1_oviDays_pois,m1_oviDays_nb1,m1_oviDays_nb2)#negative binomial best, 1 slightly better than 2

m0_oviDays_nb1=glmmTMB(OviDays~(1|population_name),data=mother_toPlot,family=nbinom1())
AICctab(m1_oviDays_nb1,m0_oviDays_nb1)


hist(resid(m1_oviDays_nb1))
m1_oviDays_simulationOutput <- simulateResiduals(m1_oviDays_nb1, n = 500)
testOutliers(m1_oviDays_simulationOutput,type="bootstrap",nBoot=100)
testResiduals(m1_oviDays_simulationOutput)
testDispersion(m1_oviDays_simulationOutput)
plot(m1_oviDays_simulationOutput) 
testQuantiles(m1_oviDays_simulationOutput)
testZeroInflation(m1_oviDays_simulationOutput)#significant "zero inflation" because there are no zeroes in the dataset
check_singularity(m1_oviDays_simulationOutput)

plot_model(m1_oviDays_nb1,type="re")
plot_model(m1_oviDays_nb1,type="est",show.values=T,p.adjust="fdr")
plot_model(m1_oviDays_nb1,type="pred",terms=c("CS_lat"))
plot_model(m1_oviDays_nb1,type="pred",terms=c("CS_LegLengthM"))

performance::r2(m1_oviDays_nb1)

cor(model.response(model.frame(m1_oviDays_nb1)),predict(m1_oviDays_nb1,type="response"))^2
length(na.omit(mother_toPlot$OviDays))
confint(m1_oviDays_nb1)
summary(m1_oviDays_nb1)

pred_ovi_lat=ggpredict(m1_oviDays_nb1,terms="CS_lat")
pred_ovi_lat$x.orig=pred_ovi_lat$x*attr(mother_toPlot$CS_lat, 'scaled:scale')+attr(mother_toPlot$CS_lat,'scaled:center')

mother_ovi_plot_lat=ggplot(data=pred_ovi_lat,aes(x=x.orig))+
  geom_point(data=mother_toPlot,aes(x=GPS_N,y=OviDays,color=Group,shape=Country),alpha=0.6)+
  scale_color_manual(values=c("#0082c8","#f58231"))+
  scale_shape_manual(values=c(17,16,15))+
  geom_ribbon(aes(x=x.orig,ymin=conf.low,ymax=conf.high),alpha=0.25)+
  geom_line(aes(y=predicted),linewidth=1.2)+
  xlab("Latitude (°N)")+
  ylab("Oviposition latency (days)")+
  theme_bw()+
  xlim(43,58)+
  theme(axis.text.x=element_text(size=6),axis.text.y=element_text(size=6))+
  theme(axis.title=element_text(size=7))+
  theme(axis.title.y=element_text(margin=margin(r=4)),axis.title.x=element_text(margin=margin(t = 4)))+ 
  theme(legend.position="none")
mother_ovi_plot_lat
pdf("output/Figures/ovipositionLatency_final.pdf",width=2.5,height=2.5)
mother_ovi_plot_lat
dev.off()

