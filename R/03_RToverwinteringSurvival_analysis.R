### This script analyzes the overwintering survival of Argiope bruennichi offspring reared in a reciprocal common garden experiment
#data should be imported and organized in the 01_import_data.R script, so it is ready to start with here

#load libraries
library(glmmTMB)
library(DHARMa)
library(performance)
library(bbmle) #for AICtab
library(sjPlot)
library(ggplot2)
library(MuMIn)
library(effects)
library(ggeffects)

set_theme(theme_bw())

str(count)
count=count[which(count$Total>0),] #simplify data frame name and remove zeros, which can't exist in reality, only by accidents (i.e. eggs dropped and not wrapped into egg sac)

### only looking at 1st egg sacs, no 2nd egg sacs
count$Cocoon_number=as.character(count$Cocoon_number)
count=count[count$Cocoon_number=="1",]


#calculate variables needed to then make matrix (success/failure) response variables (hatching success, survival) from count data
count$Hatched=count$Alive.Number+count$Dead.Number

hatch_temp=cbind(hatched=count$Hatched,unhatched=count$Egg.Number)
count$HatchingSuccess=hatch_temp
count$Hatching_proportion=count$Hatched/count$Total

survival_temp=cbind(alive=count$Alive.Number,dead=count$Dead.Number)
count$Survival=survival_temp
count$Survival_proportion=count$Alive.Number/count$Hatched
count$Survival_proportion[is.nan(count$Survival_proportion)]<-0

#add individual-level ID for individual-level random effect
count$Mother_ES=paste0(count$MotherID,sep="_",count$Cocoon_number)
count$Counted_ES=paste0(count$Counted.by,sep="_",count$Mother_ES)
count$Origin_ES=paste0(count$Origin,"_",count$Cocoon_number)

#calculate a column that shows the differences in oviposition dates as numeric values (days), so that I can use them in a model as a fixed effect
count$Oviposition_days=as.numeric(as.Date(count$Oviposition_date,format="%d/%m/%Y")-as.Date(count$CollectionDate,format="%d/%m/%Y"),units="days")

#scale predictor variables
count$CS_Mother_LegLength_um=scale(count$Mother_LegLength_um)
count$CS_Oviposition_days=scale(count$Oviposition_days)
count$CS_squared_OviDays=scale(I(count$Oviposition_days)^2)
count$CS_divided_OviDays=scale(I(1/count$Oviposition_days))
count$CS_famFecundity=scale(count$Total)
#check if continuous variables are correlated
cor.test(count$CS_Mother_LegLength_um,count$CS_Oviposition_days,method="pearson")

count$Origin=factor(count$Origin, levels=c("core","edge"))
count$Winter.treatment=factor(count$Winter.treatment, levels=c("warm","cold"))

####glm for offspring survival####
hist(count$Survival_proportion)
length(count$Survival_proportion)
aggregate(Survival_proportion~Winter.treatment*Origin,data=count,FUN=length)

surv_m1_bin_logit=glmmTMB(Survival~Winter.treatment*Origin+CS_Mother_LegLength_um+CS_Oviposition_days+CS_famFecundity, data=count,family=binomial(link = "logit"))
surv_m1_betabin=glmmTMB(Survival~Winter.treatment*Origin+CS_Mother_LegLength_um+CS_Oviposition_days+CS_famFecundity, data=count,family=betabinomial(link = "logit"))
surv_m1_bin_probit=glmmTMB(Survival~Winter.treatment*Origin+CS_Mother_LegLength_um+CS_Oviposition_days+CS_famFecundity, data=count,family=binomial(link = probit))
surv_m1_bin_cloglog=glmmTMB(Survival~Winter.treatment*Origin+CS_Mother_LegLength_um+CS_Oviposition_days+CS_famFecundity, data=count,family=binomial(link = cloglog))

AICtab(surv_m1_bin_logit,surv_m1_betabin,surv_m1_bin_probit,surv_m1_bin_cloglog)##betabinomial is the best by far

surv_m0_betabin_logit=glmmTMB(Survival~(1), data=count,family=betabinomial(link = "logit"))
surv_m1_betabin_logit=glmmTMB(Survival~Winter.treatment*Origin+CS_Mother_LegLength_um+CS_Oviposition_days+CS_famFecundity, data=count,family=betabinomial(link = "logit"))
surv_m2_betabin_logit=glmmTMB(Survival~Winter.treatment+Origin+CS_Mother_LegLength_um+CS_Oviposition_days+CS_famFecundity, data=count,family=betabinomial(link = "logit"))


AICtab(surv_m0_betabin_logit,surv_m1_betabin_logit,surv_m2_betabin_logit)#model with interaction is within delta2 AICc so including that since the interaction is part of the experiment

plot_model(surv_m1_betabin_logit,type="pred",terms=c("Winter.treatment","Origin"))

hist(resid(surv_m1_betabin_logit))
check_collinearity(surv_m1_betabin_logit)
m1_betabin_surv_simulationOutput <- simulateResiduals(surv_m1_betabin_logit, n = 500)
plot(m1_betabin_surv_simulationOutput) 
testResiduals(m1_betabin_surv_simulationOutput)
testDispersion(m1_betabin_surv_simulationOutput) 
testZeroInflation(m1_betabin_surv_simulationOutput) #definite problem with zero inflation
testQuantiles(m1_betabin_surv_simulationOutput) #quantiles look fine

#try zero inflation formula with oviposition days since survival tends to be worse the later they were laid, to deal with zero inflation
surv_m1_betabin_ZIoviDays=update(surv_m1_betabin_logit,.~.,ziformula=~CS_Oviposition_days)


AICtab(surv_m1_betabin_ZIoviDays,surv_m1_betabin_logit) #zero inflation factor improves AICc quite a bit

hist(resid(surv_m1_betabin_ZIoviDays))
check_collinearity(surv_m1_betabin_ZIoviDays)#
m1_betabin_surv_ZI <- simulateResiduals(surv_m1_betabin_ZIoviDays, n = 500)
plot(m1_betabin_surv_ZI) 
testResiduals(m1_betabin_surv_ZI)
testDispersion(m1_betabin_surv_ZI) 
testZeroInflation(m1_betabin_surv_ZI) #fixed zero inflation!
testQuantiles(m1_betabin_surv_ZI) #quantiles look even better now

summary(surv_m1_betabin_ZIoviDays)
performance::r2(surv_m1_betabin_ZIoviDays)# doesn't work for models with zero inflation factor
surv_response=apply(get.response(model.frame(surv_m1_betabin_ZIoviDays)),1,function(x)x[1]/sum(x))
surv_response[which(surv_response=="NaN")]=0
cor(surv_response,predict(surv_m1_betabin_ZIoviDays,type="response"))^2

length(na.omit(count$Survival))/2
confint(surv_m1_betabin_ZIoviDays)

plot_model(surv_m1_betabin_ZIoviDays,type="re")
plot_model(surv_m1_betabin_ZIoviDays,type="est",p.adjust="fdr")
plot_model(surv_m1_betabin_ZIoviDays,type="pred",terms=c("Winter.treatment","Origin"))
plot_model(surv_m1_betabin_ZIoviDays,type="pred",terms="Origin")

pred_surv_treatOrigin=ggpredict(surv_m1_betabin_ZIoviDays,terms=c("Winter.treatment","Origin"))

plot(pred_surv_treatOrigin,add.data=T)+
  ylim(0,1) 

pd1=position_dodge(0.6)

finalSurvPlot=ggplot(data=pred_surv_treatOrigin,aes(x=x,y=predicted)) +
  geom_violin(data=count,aes(x=Winter.treatment,y=Survival_proportion,color=Origin,fill=Origin),alpha=0.1,position=pd1,linewidth=0.2)+
  geom_point(data=count,aes(x=Winter.treatment,y=Survival_proportion,color=Origin,fill=Origin,shape=Origin),alpha=0.5,position=position_jitterdodge(jitter.width = 0.2,dodge.width = 0.6),size=0.5)+
  geom_point(data=pred_surv_treatOrigin,aes(x=x,y=predicted,shape=group),position=pd1,size=1.3)+
  geom_errorbar(data=pred_surv_treatOrigin,aes(ymin=conf.low,ymax=conf.high,group=group),width=0.1,position=pd1,linewidth=0.35)+
  scale_color_manual(values=c("#f58231", "#0082c8"))+
  scale_fill_manual(values=c("#f58231", "#0082c8"))+
  ylab(label="Survival proportion")+
  xlab(label="Winter treatment")+
  theme_bw()+
  theme(axis.text.x=element_text(size=5),axis.text.y=element_text(size=5))+
  theme(axis.title=element_text(size=6))+
  theme(axis.title.y=element_text(margin=margin(r=4)),axis.title.x=element_text(margin=margin(t=4)))+ 
  theme(legend.position="none")
finalSurvPlot

pdf("output/Figures/overwinteringsurvival_final.pdf",width=1.8,height=2)
finalSurvPlot
dev.off()
