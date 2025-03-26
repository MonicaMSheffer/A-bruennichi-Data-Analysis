### This script analyzes the body size and condition of Argiope bruennichi offspring reared in a reciprocal common garden experiment
#data should be imported and organized in the 01_import_data.R script, so it is ready to start with here. Some variables related to number of offspring within a cocoon are calculated in the 03_RToverwinteringSurvival_analysis.R script, so that should be run before this one

#load libraries
library(beanplot)
library(glmmTMB)
library(DHARMa)
library(performance)
library(bbmle) #for AICtab
library(sjPlot)
library(ggplot2)
library(MuMIn)
library(effects)
library(ggeffects)
library(ggpubr)
library(lmtest)

str(size)

#add column for observation-level random effect, if necessary
size$Mother_ES=paste0(size$MotherID,sep="_",size$Cocoon_number)
size$Cocoon_number=as.character(size$Cocoon_number)

#add fecundity of family as a column
index_mother_fecundity=match(size$Mother_ES,count$Mother_ES)
size=cbind(size,count$Total[index_mother_fecundity])
colnames(size)=c(colnames(size[1:29]),"fam_fecundity")

#add number of live spiderlings as a column
index_mother_liveTotal=match(size$Mother_ES,count$Mother_ES)
size=cbind(size,count$Alive.Number[index_mother_liveTotal])
colnames(size)=c(colnames(size[1:30]),"fam_liveTotal")

#calculate oviposition days again
size$Oviposition_days=as.numeric(size$Oviposition_date-as.Date("02/08/2018",format="%d/%m/%Y"),units="days")

beanplot(Mass..ug.~Origin,data=size)

#create a column with yes/no categories for deformation
size$deformed_yesno=size$Notes
levels(as.factor(size$deformed_yesno))
size$deformed_yesno[which(size$deformed_yesno=="deformed"|size$deformed_yesno=="deformed legs"|size$deformed_yesno=="deformed legs, photo of left leg"|size$deformed_yesno=="deformed, left"|size$deformed_yesno=="deformed; photo of left leg"|size$deformed_yesno=="freezer confusion, photo of left leg"|size$deformed_yesno=="mini spider")]="yes"
size$deformed_yesno[which(size$deformed_yesno!="yes" | is.na(size$deformed_yesno))]="no"
size$deformed_yesno=as.factor(size$deformed_yesno)

#subset to only include post-winter 1st egg sacs
size=size[size$Cocoon_number=="1",]

#center and scale continuous predictor variables
size$CS_Oviposition_days=scale(size$Oviposition_days)
size$CS_Mother_LegLength_um=scale(size$Mother_LegLength_um)
size$CS_famFecundity=scale(size$fam_fecundity)
size$CS_famLiveTotal=scale(size$fam_liveTotal)
size$CS_famLiveTotal=scale(size$fam_liveTotal)

boxplot(size$Mass..ug.~size$Origin)
boxplot(size$Leg.length_corrected~size$Origin)
boxplot(size$Leg.length_corrected~size$Winter.treatment*size$Origin)
colnames(size)[16]="Mass_ug"

beanplot(Mass_ug~Winter.treatment*Origin,data=size)

###only modeling mass changes
m1_mass_gausIdent=glmmTMB(Mass_ug~Winter.treatment*Origin+CS_Mother_LegLength_um+CS_Oviposition_days+CS_famFecundity+(1|MotherID),data=size,family=gaussian(link="identity"))
m1_mass_gausLog=glmmTMB(Mass_ug~Winter.treatment*Origin+CS_Mother_LegLength_um+CS_Oviposition_days+CS_famFecundity+(1|MotherID),data=size,family=gaussian(link="log"))
AICctab(m1_mass_gausIdent,m1_mass_gausLog)

m0_mass_gausLog=glmmTMB(Mass_ug~(1|MotherID),data=size,family=gaussian(link="log"))
m2_mass_gausLog=glmmTMB(Mass_ug~Winter.treatment+Origin+CS_Mother_LegLength_um+CS_Oviposition_days+CS_famFecundity+(1|MotherID),data=size,family=gaussian(link="log"))
m3_mass_gausLog=glmmTMB(Mass_ug~Winter.treatment*Origin+CS_Mother_LegLength_um+CS_Oviposition_days+CS_famFecundity,data=size,family=gaussian(link="log"))#model without random effects so we can get a rough estimate of marginal vs conditional pseudoR2

AICtab(m0_mass_gausLog,m1_mass_gausLog,m2_mass_gausLog,m3_mass_gausLog)

hist(resid(m1_mass_gausLog))
check_collinearity(m1_mass_gausLog)
M1_mass_gauss_simulationOutput <- simulateResiduals(m1_mass_gausLog, n = 500)
testResiduals(M1_mass_gauss_simulationOutput)
testDispersion(M1_mass_gauss_simulationOutput)
plot(M1_mass_gauss_simulationOutput) 
testQuantiles(M1_mass_gauss_simulationOutput) #significant deviations - need to figure out what that's about, but maybe okay to ignore given that the "qqplot" style thing looks alright to me

summary(m1_mass_gausLog)
length(na.omit(size$Mass_ug))
cor(model.response(model.frame(m1_mass_gausLog)),predict(m1_mass_gausLog,type="response"))^2 #conditional R2
cor(model.response(model.frame(m3_mass_gausLog)),predict(m3_mass_gausLog,type="response"))^2 #marginal R2

confint(m1_mass_gausLog)

plot_model(m1_mass_gausLog,type="re")
plot_model(m1_mass_gausLog,type="est",p.adjust="fdr",show.values=T)

pred_mass_TreatOrigin=ggpredict(m1_mass_gausLog,terms=c("Winter.treatment","Origin"))
plot(pred_mass_TreatOrigin)
size$Winter.treatment=factor(size$Winter.treatment,levels=c("warm","cold"))
pd1=position_dodge(0.6)
final_mass_plot=ggplot(data=pred_mass_TreatOrigin,aes(x=x,y=predicted)) +
  geom_violin(data=size,aes(x=Winter.treatment,y=Mass_ug,color=Origin,fill=Origin),alpha=0.15,position=pd1,linewidth=0.2)+
  geom_point(data=size,aes(x=Winter.treatment,y=Mass_ug,color=Origin,fill=Origin,shape=Origin),alpha=0.5,position=position_jitterdodge(jitter.width = 0.2,dodge.width = 0.6),size=.5)+
  geom_point(data=pred_mass_TreatOrigin,aes(x=factor(x,levels=c("warm","cold")),y=predicted,shape=factor(group,levels=c("core","edge"))),position=pd1,size=1.3)+
  geom_errorbar(data=pred_mass_TreatOrigin,aes(ymin=conf.low,ymax=conf.high,group=factor(group,levels=c("core","edge"))),width=0.1,position=pd1,linewidth=0.35)+
  scale_color_manual(values=c("#f58231", "#0082c8"))+
  scale_fill_manual(values=c("#f58231", "#0082c8"))+
  ylab(label="Mass (µg)")+
  xlab(label="Winter treatment")+
  theme(axis.text.x=element_text(size=5),axis.text.y=element_text(size=5))+
  theme(axis.title=element_text(size=6))+
  theme(axis.title.y=element_text(margin=margin(r=6)),axis.title.x=element_text(margin=margin(t=4)))+
  theme(legend.position="none")
final_mass_plot

pdf("output/Figures/offspringmass_final.pdf",width=1.8,height=2)
final_mass_plot
dev.off()


###look at body condition index (ends up in supplement)
hist(size$Leg.length_corrected)
hist(size$Mass_ug)
plot(Mass_ug~Leg.length_corrected,data=size)

#remove individuals with no leg length measurement
size_noNA=size[!is.na(size$Leg.length_corrected),]
size_noNA=size_noNA[!is.na(size_noNA$Mass_ug),]

#linear regression of mass against leg length
lm_bodyCond <- lm(Mass_ug ~ Leg.length_corrected, data = size_noNA)
#check for patterns in the residuals
size_noNA$resi <- lm_bodyCond$residuals
#library(lmtest)
ggplot(data = size_noNA, aes(y = resi, x = Leg.length_corrected)) + geom_point(col = 'blue') + geom_abline(slope = 0)
bptest(lm_bodyCond)#no heteroskedasticity 

#use the residuals of the linear regression as the response variable
hist(size_noNA$resi)
m1_bci_gausIdent=glmmTMB(resi~Winter.treatment*Origin+CS_Mother_LegLength_um+CS_Oviposition_days+CS_famFecundity+(1|MotherID),data=size_noNA,family=gaussian(link="identity"))
#identity link is the only possible one since values can be negative
m0_bci_gausIdent=glmmTMB(resi~(1|MotherID),data=size_noNA,family=gaussian(link="identity"))
m2_bci_gausIdent=glmmTMB(resi~Winter.treatment+Origin+CS_Mother_LegLength_um+CS_Oviposition_days+CS_famFecundity+(1|MotherID),data=size_noNA,family=gaussian(link="identity"))
AICtab(m1_bci_gausIdent,m0_bci_gausIdent,m2_bci_gausIdent)

hist(resid(m1_bci_gausIdent))
m1_BCI_gauss_simulationOutput <- simulateResiduals(m1_bci_gausIdent, n = 500)
testResiduals(m1_BCI_gauss_simulationOutput)
testDispersion(m1_BCI_gauss_simulationOutput)
plot(m1_BCI_gauss_simulationOutput) 
testQuantiles(m1_BCI_gauss_simulationOutput) #significant deviations - need to figure out what that's about, but maybe okay to ignore given that nothing is too terribly severe

summary(m1_bci_gausIdent)
length(na.omit(size_noNA$resi))
performance::r2_nakagawa(m1_bci_gausIdent)

confint(m1_bci_gausIdent)

plot_model(m1_bci_gausIdent,type="re")
plot_model(m1_bci_gausIdent,type="est",p.adjust="fdr",show.values=T)

pred_bci_TreatOrigin=ggpredict(m1_bci_gausIdent,terms=c("Winter.treatment","Origin"))
plot(pred_bci_TreatOrigin)

size_noNA$Winter.treatment=factor(size_noNA$Winter.treatment,levels=c("warm","cold"))

pd1=position_dodge(0.6)
final_bci_plot=ggplot(data=pred_bci_TreatOrigin,aes(x=x,y=predicted)) +
  geom_violin(data=size_noNA,aes(x=Winter.treatment,y=resi,color=Origin,fill=Origin),alpha=0.1,position=pd1,linewidth=0.1)+
  geom_point(data=size_noNA,aes(x=Winter.treatment,y=resi,color=Origin,fill=Origin,shape=Origin),alpha=0.5,position=position_jitterdodge(jitter.width = 0.2,dodge.width = 0.6),size=0.5)+
  geom_point(data=pred_bci_TreatOrigin,aes(x=x,y=predicted,shape=factor(group,levels=c("core","edge"))),position=pd1,size=1.3)+
  geom_errorbar(data=pred_bci_TreatOrigin,aes(ymin=conf.low,ymax=conf.high,group=factor(group,levels=c("core","edge"))),width=0.1,position=pd1,linewidth=0.35)+
  scale_color_manual(values=c("#f58231", "#0082c8"))+
  scale_fill_manual(values=c("#f58231", "#0082c8"))+
  ylab(label="Body condition index")+
  xlab(label="Winter treatment")+
  theme(axis.text.x=element_text(size=5),axis.text.y=element_text(size=5))+
  theme(axis.title=element_text(size=6))+
  theme(axis.title.y=element_text(margin=margin(r=6)),axis.title.x=element_text(margin=margin(t=4)))+
  theme(legend.text=element_text(size=5),legend.title=element_text(size=6),legend.position="none")
final_bci_plot
pdf("output/Figures/offspringcondition_supplement_final.pdf",width=2.5,height=2.5)
final_bci_plot
dev.off()