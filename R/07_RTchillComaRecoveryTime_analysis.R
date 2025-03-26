### This script analyzes the chill coma recovery time of Argiope bruennichi offspring reared in a reciprocal common garden experiment
#data should be imported and organized in the 01_import_data.R script, so it is ready to start with here. Some variables related to number of offspring within a cocoon are calculated in the 03_RToverwinteringSurvival_analysis.R script, so that should be run before this one


library(beanplot)
#install.packages("glmmTMB", repos="https://glmmTMB.github.io/glmmTMB/repos")
library(glmmTMB)
library(DHARMa)
library(performance)
library(bbmle) #for AICtab
library(sjPlot)
library(ggplot2)
library(MuMIn)
library(effects)
library(ggeffects)

set_theme(base=theme_bw())

str(ccr)
levels(as.factor(ccr$Notes))
hist(ccr$CCR.time.seconds)

#remove NAs from log_CCR_seconds column
ccr=ccr[which(!is.na(ccr$CCR.time.seconds)),]

#add individual-level ID for individual-level random effect if needed
ccr$Mother_ES=paste0(ccr$MotherID,sep="_",ccr$Cocoon_number)
ccr$counted_ES=paste0(ccr$ccred.by,sep="_",ccr$Mother_ES)
ccr$Origin_ES=paste0(ccr$Origin,"_",ccr$Cocoon_number)

#add fecundity of family as a column
index_mother_fecundity=match(ccr$Mother_ES,count$Mother_ES)
ccr=cbind(ccr,count$Total[index_mother_fecundity])
colnames(ccr)=c(colnames(ccr[1:30]),"fam_fecundity")

#calculate a column that shows the differences in oviposition dates as numeric values, so that I can use them in a model as a fixed effect
ccr$Oviposition_days=as.numeric(ccr$Oviposition_date-as.Date("02/08/2018",format="%d/%m/%Y"),units="days")

#include only first egg sacs:
ccr=ccr[ccr$Cocoon_number=="1",]

hist(ccr$CCR.time.seconds)
hist(log(ccr$CCR.time.seconds))
beanplot(ccr$CCR.time.seconds~ccr$Origin)
plot(CCR.time.seconds~Oviposition_days,data=ccr)
plot(CCR.time.seconds~Mother_LegLength_um,data=ccr)
plot(CCR.time.seconds~fam_fecundity,data=ccr)
beanplot(CCR.time.seconds~Origin*Winter.treatment,data=ccr)

boxplot(ccr$CCR.time.seconds)
boxplot(log(ccr$CCR.time.seconds))
#scale predictor variables
ccr$CS_Mother_LegLength_um=scale(ccr$Mother_LegLength_um)
ccr$CS_Oviposition_days=scale(ccr$Oviposition_days)
ccr$CS_famFecundity=scale(ccr$fam_fecundity)

ccr$ccr_minutes=ccr$CCR.time.seconds/60

#transform CCR 
ccr$log_CCR_minutes=log(ccr$ccr_minutes)

#scale predictor variables
ccr$CS_Mother_LegLength_um=scale(ccr$Mother_LegLength_um)
ccr$CS_Oviposition_days=scale(ccr$Oviposition_days)
ccr$CS_famFecundity=scale(ccr$fam_fecundity)

#check sample sizes
aggregate(log_CCR_minutes~Winter.treatment*Origin,data=ccr,FUN=length)

m1_ccr_simple_gauslog=glmmTMB(log_CCR_minutes~Winter.treatment*Origin+CS_Mother_LegLength_um+CS_Oviposition_days+CS_famFecundity+(1|MotherID),data=ccr,family=gaussian(link="log"))
m1_ccr_simple_gausident=glmmTMB(log_CCR_minutes~Winter.treatment*Origin+CS_Mother_LegLength_um+CS_Oviposition_days+CS_famFecundity+(1|MotherID),data=ccr,family=gaussian(link="identity"))

bbmle::AICctab(m1_ccr_simple_gausident,m1_ccr_simple_gauslog)

m1_ccr_simple=glmmTMB(log_CCR_minutes~Winter.treatment*Origin+CS_Mother_LegLength_um+CS_Oviposition_days+CS_famFecundity+(1|MotherID),data=ccr,family=gaussian(link="log"))
m2_ccr_simple=glmmTMB(log_CCR_minutes~Winter.treatment+Origin+CS_Mother_LegLength_um+CS_Oviposition_days+CS_famFecundity+(1|MotherID),data=ccr,family=gaussian(link="log"))
m0_ccr_simple=glmmTMB(log_CCR_minutes~(1|MotherID),data=ccr,family=gaussian(link="log"))
AICtab(m1_ccr_simple,m0_ccr_simple,m2_ccr_simple)
m3_ccr_simple=glmmTMB(log_CCR_minutes~Winter.treatment*Origin+CS_Mother_LegLength_um+CS_Oviposition_days+CS_famFecundity,data=ccr,family=gaussian(link="log"))

AICctab(m1_ccr_simple,m3_ccr_simple)

hist(resid(m1_ccr_simple))
check_collinearity(m1_ccr_simple)
M1gauss_simulationOutput <- simulateResiduals(m1_ccr_simple, n = 500)
testResiduals(M1gauss_simulationOutput)
testDispersion(M1gauss_simulationOutput)
plot(M1gauss_simulationOutput) 
testQuantiles(M1gauss_simulationOutput)#significant deviations but visually tolerable

summary(m1_ccr_simple)
length(na.omit(ccr$ccr_minutes))
performance::r2_nakagawa(m1_ccr_simple,tolerance=1e-11)

cor(model.response(model.frame(m1_ccr_simple)),predict(m1_ccr_simple,type="response"))^2 #conditional R2 (seemingly poor estimate)
cor(model.response(model.frame(m3_ccr_simple)),predict(m3_ccr_simple,type="response"))^2 #marginal R2 (seemingly poor estimate)
r.squaredGLMM(m1_ccr_simple)

confint(m1_ccr_simple)

plot_model(m1_ccr_simple,type="est",sort.est = T,show.values=T,p.adjust="fdr")

pred_ccr_TreatOrigin=ggpredict(m1_ccr_simple,terms=c("Winter.treatment","Origin"))
plot(pred_ccr_TreatOrigin)

transformed_pred_ccr_TreatOrigin=pred_ccr_TreatOrigin
transformed_pred_ccr_TreatOrigin$predicted=exp(transformed_pred_ccr_TreatOrigin$predicted)
transformed_pred_ccr_TreatOrigin$conf.low=exp(transformed_pred_ccr_TreatOrigin$conf.low)
transformed_pred_ccr_TreatOrigin$conf.high=exp(transformed_pred_ccr_TreatOrigin$conf.high)
transformed_pred_ccr_TreatOrigin

pred_ccr_Treat=ggpredict(m1_ccr_simple,terms=c("Winter.treatment"))
plot(pred_ccr_Treat)

pd1=position_dodge(0.6)
finalCCRplot=ggplot(data=pred_ccr_TreatOrigin,aes(x=x,y=predicted)) +
  geom_violin(data=ccr,aes(x=factor(Winter.treatment,levels=c('warm','cold')),y=exp(log_CCR_minutes),color=Origin,fill=Origin),alpha=0.15,position=pd1,linewidth=0.2)+
  geom_point(data=ccr,aes(x=factor(Winter.treatment,levels=c('warm','cold')),y=exp(log_CCR_minutes),color=Origin,fill=Origin,shape=Origin),alpha=0.5,position=position_jitterdodge(jitter.width = 0.2,dodge.width = 0.6),size=0.5)+
  geom_point(data=pred_ccr_TreatOrigin,aes(x=factor(x,levels=c("warm","cold")),y=exp(predicted),shape=factor(group,levels=c("core","edge"))),position=pd1,size=1.3)+
  geom_errorbar(data=pred_ccr_TreatOrigin,aes(ymin=exp(conf.low),ymax=exp(conf.high),group=factor(group,levels=c("core","edge"))),width=0.1,position=pd1,linewidth=.35)+
  scale_color_manual(values=c("#f58231", "#0082c8"))+
  scale_fill_manual(values=c("#f58231", "#0082c8"))+
  ylab(label="CCR Time (minutes)")+
  xlab(label="Winter treatment")+
  theme_bw()+
  theme(axis.text.x=element_text(size=5),axis.text.y=element_text(size=5))+
  theme(axis.title=element_text(size=6))+
  theme(axis.title.y=element_text(margin=margin(r=6)),axis.title.x=element_text(margin=margin(t=4)))+
  theme(legend.position="none")
finalCCRplot

pdf("output/Figures/chillcomarecovery_final.pdf",width=1.8,height=2)
finalCCRplot
dev.off()
