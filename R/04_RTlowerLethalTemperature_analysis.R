### This script analyzes the low temperature survival of Argiope bruennichi offspring reared in a reciprocal common garden experiment
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

str(llt)
#make survival data model-able with binomial model 
temp_survival=cbind(alive=llt$test_live,dead=llt$test_dead)
llt$Survival=temp_survival
#View(llt$Survival)
llt$survival_proportion=llt$test_live/llt$test_total
hist(llt$survival_proportion)

#remove NAs from survival column
llt1=llt[which(!is.na(llt$survival_proportion)),]
#remove second egg sacs
llt=llt1[llt1$Cocoon_number=="1",]
#get sample sizes
aggregate(survival_proportion~Winter.treatment*Origin,data=llt,FUN=length)
#add individual-level ID for individual-level random effect if needed

llt$Mother_ES=paste0(llt$MotherID,sep="_",llt$Cocoon_number)
llt$counted_ES=paste0(llt$Counted.by,sep="_",llt$Mother_ES)
llt$Origin_ES=paste0(llt$Origin,"_",llt$Cocoon_number)

#calculate a column that shows the differences in oviposition dates as numeric values, so that I can use them in a model as a fixed effect
llt$Oviposition_days=as.numeric(llt$Oviposition_date-as.Date("02/08/2018",format="%d/%m/%Y"),units="days")

#check for mortality in control groups
sum(llt$control_dead)
sum(llt$control_total)
#mortality was so low, I'm not concerned about accounting for it in the formal analysis
#scale predictor variables
llt$CS_Mother_LegLength_um=scale(llt$Mother_LegLength_um)
llt$CS_Oviposition_days=scale(llt$Oviposition_days)

m1_llt_simple_betabin=glmmTMB(Survival~LT50.Temp*Origin*Winter.treatment+CS_Mother_LegLength_um+CS_Oviposition_days,data=llt,family=betabinomial())
m1_llt_simple_bin=glmmTMB(Survival~LT50.Temp*Origin*Winter.treatment+CS_Mother_LegLength_um+CS_Oviposition_days,data=llt,family=binomial())
AICctab(m1_llt_simple_betabin,m1_llt_simple_bin)#betabinomial best by far

m0_llt_simple=glmmTMB(Survival~1,data=llt,family=betabinomial())
m1_llt_simple=glmmTMB(Survival~LT50.Temp*Origin*Winter.treatment+CS_Mother_LegLength_um+CS_Oviposition_days,data=llt,family=betabinomial())
m2_llt_simple=glmmTMB(Survival~LT50.Temp+Origin*Winter.treatment+CS_Mother_LegLength_um+CS_Oviposition_days,data=llt,family=betabinomial())
m3_llt_simple=glmmTMB(Survival~LT50.Temp*Origin+Winter.treatment+CS_Mother_LegLength_um+CS_Oviposition_days,data=llt,family=betabinomial())
m4_llt_simple=glmmTMB(Survival~LT50.Temp+Origin+Winter.treatment+CS_Mother_LegLength_um+CS_Oviposition_days,data=llt,family=betabinomial())

AICtab(m0_llt_simple,m1_llt_simple,m2_llt_simple,m3_llt_simple,m4_llt_simple)
AICctab(m1_llt_simple,m2_llt_simple,m3_llt_simple,m4_llt_simple)#model with all interactions best

hist(resid(m1_llt_simple))
M1betabin_simulationOutput <- simulateResiduals(m1_llt_simple, n = 500)
testResiduals(M1betabin_simulationOutput)
testOutliers(M1betabin_simulationOutput,type="bootstrap")
testDispersion(M1betabin_simulationOutput)
plot(M1betabin_simulationOutput)
testZeroInflation(M1betabin_simulationOutput)
#everything looks good!

summary(m1_llt_simple)
performance::r2(m1_llt_simple,tolerance=1e-10)

llt_response=apply(get.response(model.frame(m1_llt_simple)),1,function(x)x[1]/sum(x))
cor(llt_response,predict(m1_llt_simple,type="response"))^2

confint(m1_llt_simple)
length(na.omit(llt$survival_proportion))
length(na.omit(llt$Survival))/2

plot_model(m1_llt_simple,type="est",sort.est=T,show.values=T,p.adjust="fdr")
plot_model(m1_llt_simple,type="pred",terms=c("LT50.Temp","Origin","Winter.treatment"))
plot_model(m1_llt_simple,type="pred",terms=c("LT50.Temp","Origin"))
plot_model(m1_llt_simple,type="pred",terms=c("Origin"))
plot_model(m1_llt_simple,type="pred",terms=c("Winter.treatment"))
plot_model(m1_llt_simple,type="pred",terms=c("LT50.Temp","Origin"))
plot_model(m1_llt_simple,type="pred",terms=c("LT50.Temp"))

pred_llt_temoOrigin=ggpredict(m1_llt_simple,terms=c("LT50.Temp","Origin"))
pred_llt_TempTreatOrigin=ggpredict(m1_llt_simple,terms=c("LT50.Temp","Origin","Winter.treatment"))
colnames(pred_llt_TempTreatOrigin)<-c("x","predicted","std.error","conf.low","conf.high","Origin","Winter.treatment") 

llt_finalPlot=ggplot()+
  geom_jitter(data=llt,aes(x=LT50.Temp,y=survival_proportion,color=Origin,fill=Origin,shape=Origin),alpha=0.5,width=0.5,show.legend=T,size=1)+
  scale_shape_manual(values=c(16,17))+
  facet_grid(cols=vars(factor(Winter.treatment,levels=c("warm","cold"))))+
  geom_ribbon(data=pred_llt_TempTreatOrigin,mapping=aes(x=x,ymin=conf.low,ymax=conf.high,fill=factor(Origin,levels=c("core","edge"))),alpha=0.15,show.legend=T)+
  geom_line(data=pred_llt_TempTreatOrigin,aes(x=x,y=predicted,color=factor(Origin,levels=c("core","edge")),linetype=factor(Origin,levels=c("core","edge"))),show.legend=T,linewidth=.5)+ 
  theme(legend.position="bottom")+
  scale_color_manual(values=c("#f58231","#0082c8"))+
  scale_fill_manual(values=c("#f58231","#0082c8"))+
  scale_linetype_manual(values=c("dotdash","solid"))+
  ylab(label="LLT (proportion surviving)")+
  xlab(label="Exposure Temperature (°C)")+
  theme_bw()+
  theme(axis.text.x=element_text(size=5),axis.text.y=element_text(size=5))+
  theme(axis.title=element_text(size=6))+
  theme(axis.title.y=element_text(margin=margin(r=4)),axis.title.x=element_text(margin=margin(t=4)))+
  theme(strip.text.x=element_text(size=6))+
  theme(legend.position="none")
llt_finalPlot

pdf("output/Figures/lowerlethaltemperature_final.pdf",width=3.5,height=2)
llt_finalPlot
dev.off()
