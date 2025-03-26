### This script analyzes the supercooling points of Argiope bruennichi offspring reared in a reciprocal common garden experiment
#data should be imported and organized in the 01_import_data.R script, so it is ready to start with here. Some variables related to number of offspring within a cocoon are calculated in the 03_RToverwinteringSurvival_analysis.R script, so that should be run before this one

#load libraries
library(beanplot)
library(glmmTMB)
library(DHARMa)
library(performance)
library(bbmle) #for AICctab
library(sjPlot)
library(ggplot2)
library(MuMIn)
library(effects)
library(ggeffects)
library(ggpubr)

set_theme(theme_bw())

str(scp)
scp=scp[which(scp$Notes!="crushed spider" & scp$Notes!="crushed with thermocouple"& scp$Notes!="Datalogger reading off "& scp$Notes!="thermocouple problem"& scp$Notes!="thermocouple problem!"& scp$Notes!="thermocouple problem?"& scp$Notes!="too crushed to measure"& scp$Notes!="two spiders in vial"& scp$Notes!="unmeasurable" | is.na(scp$Notes)),] #get rid of problematic measurements, i.e. when the spiderling was damaged when we put it in the tube
scp=scp[which(!is.na(scp$SCP)),]

scp$Cocoon_number=as.character(scp$Cocoon_number)

#add individual-level ID for individual-level random effect if needed
scp$Mother_ES=paste0(scp$MotherID,sep="_",scp$Cocoon_number)
scp$counted_ES=paste0(scp$scped.by,sep="_",scp$Mother_ES)
scp$Origin_ES=paste0(scp$Origin,"_",scp$Cocoon_number)

#add fecundity of family as a column
index_mother_fecundity=match(scp$Mother_ES,count$Mother_ES)
scp=cbind(scp,count$Total[index_mother_fecundity])
colnames(scp)=c(colnames(scp[1:30]),"fam_fecundity")

##include only 1st egg sacs
scp=scp[scp$Cocoon_number=="1",]

#look at distribution of raw data
beanplot(SCP~Winter.treatment*Origin,data=scp)#multimodal supercooling points
median(scp$SCP)
scp$Oviposition_days=as.numeric(scp$Oviposition_date-as.Date("02/08/2018",format="%d/%m/%Y"),units="days")

#center and scale continuous predictor variables
scp$CS_Mother_LegLength_um=scale(scp$Mother_LegLength_um)
scp$CS_Oviposition_days=scale(scp$Oviposition_days)
scp$CS_famFecundity=scale(scp$fam_fecundity)

#try modeling with all the data regardless of multimodality
m1_scp_gausIdent=glmmTMB(SCP~Winter.treatment*Origin+CS_Mother_LegLength_um+CS_Oviposition_days+CS_famFecundity+(1|MotherID),data=scp,family=gaussian(link="identity"))
m1_scp_gausLog=glmmTMB(-SCP~Winter.treatment*Origin+CS_Mother_LegLength_um+CS_Oviposition_days+CS_famFecundity+(1|MotherID),data=scp,family=gaussian(link="log"))
m1_scp_gausInverse=glmmTMB(-SCP~Winter.treatment*Origin+CS_Mother_LegLength_um+CS_Oviposition_days+CS_famFecundity+(1|MotherID),data=scp,family=gaussian(link="inverse"))#doesn't converge

AICctab(m1_scp_gausIdent,m1_scp_gausLog)

m0_scp_final=glmmTMB(SCP~(1|MotherID),data=scp,family=gaussian(link="identity"))
m1_scp_final=glmmTMB(SCP~Winter.treatment*Origin+CS_Mother_LegLength_um+CS_Oviposition_days+CS_famFecundity+(1|MotherID),data=scp,family=gaussian(link="identity"))
m2_scp_final=glmmTMB(SCP~Winter.treatment+Origin+CS_Mother_LegLength_um+CS_Oviposition_days+CS_famFecundity+(1|MotherID),data=scp,family=gaussian(link="identity"))
AICctab(m0_scp_final,m1_scp_final,m2_scp_final)
AICctab(m1_scp_final,m2_scp_final)

hist(resid(m1_scp_final))
check_collinearity(m1_scp_final)
M1gauss_simulationOutput <- simulateResiduals(m1_scp_final, n = 500)
testResiduals(M1gauss_simulationOutput)
testDispersion(M1gauss_simulationOutput)
plot(M1gauss_simulationOutput) 
testQuantiles(M1gauss_simulationOutput)
###normality violations are way too bad - need to go with splitting into groups

####Establish multiple groups of data to deal with multimodal distribution of SCPs
beanplot(scp$SCP~scp$Origin,ll=0.01,beanlines="median",overallline = "median",main="All data with cutoff lines")
abline(h=-17.5,lty="dashed",col="red")
abline(h=-25.5,lty="dashed",col="blue")

scp_HG=scp[which(scp$SCP>(-18.5)),]
scp_HMG=scp[which(scp$SCP>(-25.5)),]
scp_MG=scp[which(scp$SCP>(-25.5)&scp$SCP<(-18.5)),]
scp_LG=scp[which(scp$SCP<(-25.5)),]

scp_HG$group="HG"
scp_MG$group="MG"
scp_LG$group="LG"
scp_HMG$group="HMG"

beanplot(scp_HG$SCP~scp_HG$Origin*scp_HG$Winter.treatment,ll=0.05,beanlines="median",overallline = "median",main="High SCP Group >-18.5°C")
beanplot(scp_MG$SCP~scp_MG$Origin*scp_MG$Winter.treatment,ll=0.05,beanlines="median",overallline = "median",main="-25.5°C < Middle SCP Group > -18.5°C")
beanplot(scp_LG$SCP~scp_LG$Origin*scp_LG$Winter.treatment,ll=0.05,beanlines="median",overallline = "median",main="Low SCP Group < -25.5°C")


beanplot(SCP~Origin*Winter.treatment,data=scp,ll=0.05,beanlines="median",overallline = "median",ylab="SCP (°C)")
abline(h=-25.5,lty="dashed",col="red",lwd=2)

aggregate(SCP~Winter.treatment*Origin,data=scp_LG,FUN=length)

median(scp_LG$SCP)

#model low group only
m1_LGscp_gausIdent=glmmTMB(SCP~Winter.treatment*Origin+CS_Mother_LegLength_um+CS_Oviposition_days+CS_famFecundity+(1|MotherID),data=scp_LG,family=gaussian(link="identity"))
m1_LGscp_gausLog=glmmTMB(-SCP~Winter.treatment*Origin+CS_Mother_LegLength_um+CS_Oviposition_days+CS_famFecundity+(1|MotherID),data=scp_LG,family=gaussian(link="log"))
m1_LGscp_gausInverse=glmmTMB(-SCP~Winter.treatment*Origin+CS_Mother_LegLength_um+CS_Oviposition_days+CS_famFecundity+(1|MotherID),data=scp_LG,family=gaussian(link="inverse"))

AICctab(m1_LGscp_gausIdent,m1_LGscp_gausLog)

m0_LGscp_final=glmmTMB(SCP~(1|MotherID),data=scp_LG,family=gaussian(link="identity"))
m1_LGscp_final=glmmTMB(SCP~Winter.treatment*Origin+CS_Mother_LegLength_um+CS_Oviposition_days+CS_famFecundity+(1|MotherID),data=scp_LG,family=gaussian(link="identity"))
m2_LGscp_final=glmmTMB(SCP~Winter.treatment+Origin+CS_Mother_LegLength_um+CS_Oviposition_days+CS_famFecundity+(1|MotherID),data=scp_LG,family=gaussian(link="identity"))
AICtab(m0_LGscp_final,m1_LGscp_final,m2_LGscp_final)
AICctab(m1_LGscp_final,m2_LGscp_final)

hist(resid(m1_LGscp_final))
check_collinearity(m1_LGscp_final)
M1LGgauss_simulationOutput <- simulateResiduals(m1_LGscp_final, n = 500)
testResiduals(M1LGgauss_simulationOutput)
testDispersion(M1LGgauss_simulationOutput)
plot(M1LGgauss_simulationOutput) 
testQuantiles(M1LGgauss_simulationOutput)#looks great!

summary(m1_LGscp_final)

length(na.omit(scp_LG$SCP))
performance::r2_nakagawa(m1_LGscp_final,tolerance=1e-11)
confint(m1_LGscp_final)

plot_model(m1_LGscp_final,type="est",sort.est = T,show.values=T,p.adjust="fdr")
plot_model(m1_LGscp_final,type="re")
pred_LGscp_TreatOrigin=ggpredict(m1_LGscp_final,terms=c("Winter.treatment","Origin"))
pred_LGscp_TreatOrigin
plot(pred_LGscp_TreatOrigin)

pred_LGscp_Treat=ggpredict(m1_LGscp_final,terms=c("Winter.treatment"))
plot(pred_LGscp_Treat)
pred_LGscp_origin=ggpredict(m1_LGscp_final,terms=c("Origin"))
plot(pred_LGscp_origin)


pd1=position_dodge(0.6)
scp_LG_plot=ggplot(data=pred_LGscp_TreatOrigin,aes(x=x,y=predicted)) +
  geom_violin(data=scp_LG,aes(x=factor(Winter.treatment,levels=c('warm','cold')),y=SCP,color=Origin,fill=Origin),alpha=0.15,position=pd1,linewidth=0.2)+
  geom_point(data=scp_LG,aes(x=factor(Winter.treatment,levels=c('warm','cold')),y=SCP,color=Origin,fill=Origin,shape=Origin),alpha=0.5,position=position_jitterdodge(jitter.width = 0.2,dodge.width = 0.6),size=0.5)+
  geom_point(data=pred_LGscp_TreatOrigin,aes(x=factor(x,levels=c("warm","cold")),y=predicted,shape=factor(group,levels=c("core","edge"))),position=pd1,size=1.3)+
  geom_errorbar(data=pred_LGscp_TreatOrigin,aes(ymin=conf.low,ymax=conf.high,group=factor(group,levels=c("core","edge"))),width=0.1,position=pd1,linewidth=.35)+
  scale_color_manual(values=c("#f58231", "#0082c8"))+
  scale_fill_manual(values=c("#f58231", "#0082c8"))+
  ylab(label="LG SCP (°C)")+
  xlab(label="Winter treatment")+
  theme_bw()+
  theme(axis.text.x=element_text(size=5),axis.text.y=element_text(size=5))+
  theme(axis.title=element_text(size=5))+
  theme(axis.title.y=element_text(margin=margin(r=4)),axis.title.x=element_text(margin=margin(t=4)))+
  theme(legend.position="none")
scp_LG_plot

pdf("output/Figures/LGSCP_final.pdf",width=1.8,height=2)
scp_LG_plot
dev.off()


###calculate ratio of MG:LG 
scp_LGHMG=rbind(scp_LG,scp_HMG)

motherES_ratio=aggregate(group~Mother_ES,data=scp_LGHMG,function(x){
  table(x)
})
columns_ratio=do.call(rbind,lapply(motherES_ratio$group,function(x){
  if(length(unlist(x))==1){
    newname=c("HMG","LG")[!c("HMG","LG")%in%names(x)]
    x=c(x,0)
    names(x)[2]=newname
  }
  return(x[order(names(x))])
}))
motherES_ratio_df=as.data.frame(cbind(motherES_ratio[,1],columns_ratio))
colnames(motherES_ratio_df)=c("Mother_ES","HMG","LG")
motherES_ratio_df$LG=as.numeric(motherES_ratio_df$LG)
motherES_ratio_df$HMG=as.numeric(motherES_ratio_df$HMG)
motherES_ratio_idx=match(motherES_ratio_df$Mother_ES,scp_LGHMG$Mother_ES)
scp_LGHMG_ratio=cbind(scp_LGHMG[motherES_ratio_idx,],motherES_ratio_df)
ratio_temp=cbind(LG=as.numeric(scp_LGHMG_ratio$LG),HMG=as.numeric(scp_LGHMG_ratio$HMG))
scp_LGHMG_ratio$LG_HMG_ratioResponse=ratio_temp
#View(scp_LGHMG_ratio$LG_HMG_ratioResponse)
scp_LGHMG_ratio$nTested=scp_LGHMG_ratio$LG+scp_LGHMG_ratio$HMG
scp_LGHMG_ratio$proportion_LGHMG=scp_LGHMG_ratio$LG/scp_LGHMG_ratio$nTested
boxplot(scp_LGHMG_ratio$proportion_LGHMG~scp_LGHMG_ratio$Origin*scp_LGHMG_ratio$Winter.treatment)

###now model with the newly created ratio response variable
aggregate(LG_HMG_ratioResponse~Winter.treatment*Origin,data=scp_LGHMG_ratio,FUN=length)

m1_scpRatio_bin_logit=glmmTMB(LG_HMG_ratioResponse~Winter.treatment*Origin+CS_Mother_LegLength_um+CS_Oviposition_days+CS_famFecundity+(1|MotherID), data=scp_LGHMG_ratio,family=binomial(link=logit))
m1_scpRatio_bin_probit=glmmTMB(LG_HMG_ratioResponse~Winter.treatment*Origin+CS_Mother_LegLength_um+CS_Oviposition_days+CS_famFecundity+(1|MotherID), data=scp_LGHMG_ratio,family=binomial(link=probit))
m1_scpRatio_bin_cloglog=glmmTMB(LG_HMG_ratioResponse~Winter.treatment*Origin+CS_Mother_LegLength_um+CS_Oviposition_days+CS_famFecundity+(1|MotherID), data=scp_LGHMG_ratio,family=binomial(link=cloglog))
m1_scpRatio_betabin_logit=glmmTMB(LG_HMG_ratioResponse~Winter.treatment*Origin+CS_Mother_LegLength_um+CS_Oviposition_days+CS_famFecundity+(1|MotherID), data=scp_LGHMG_ratio,family=betabinomial(link="logit"), control=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS")))

AICctab(m1_scpRatio_bin_logit,m1_scpRatio_bin_probit,m1_scpRatio_bin_cloglog,m1_scpRatio_betabin_logit)#binomial with logit link is best

#compare with and without fixed effects and with decreasing interactions
m0_scpRatio_bin_logit=glmmTMB(LG_HMG_ratioResponse~(1|MotherID), data=scp_LGHMG_ratio,family=binomial(link="logit"))
m1_scpRatio_bin_logit=glmmTMB(LG_HMG_ratioResponse~Winter.treatment*Origin+CS_Mother_LegLength_um+CS_Oviposition_days+CS_famFecundity+(1|MotherID), data=scp_LGHMG_ratio,family=binomial(link=logit))
m2_scpRatio_bin_logit=glmmTMB(LG_HMG_ratioResponse~Winter.treatment+Origin+CS_Mother_LegLength_um+CS_Oviposition_days+CS_famFecundity+(1|MotherID), data=scp_LGHMG_ratio,family=binomial(link=logit))

AICtab(m0_scpRatio_bin_logit,m1_scpRatio_bin_logit,m2_scpRatio_bin_logit)
AICctab(m1_scpRatio_bin_logit,m2_scpRatio_bin_logit)

hist(resid(m1_scpRatio_bin_logit))
check_collinearity(m1_scpRatio_bin_logit)
M1_scpRatio_gauss_simulationOutput <- simulateResiduals(m1_scpRatio_bin_logit, n = 500)
testResiduals(M1_scpRatio_gauss_simulationOutput)
testDispersion(M1_scpRatio_gauss_simulationOutput)
plot(M1_scpRatio_gauss_simulationOutput) 
testQuantiles(M1_scpRatio_gauss_simulationOutput)#acceptable

summary(m1_scpRatio_bin_logit)
length(na.omit(scp_LGHMG_ratio$LG_HMG_ratioResponse))/2
performance::r2(m1_scpRatio_bin_logit)
scp_response=apply(get.response(model.frame(m1_scpRatio_bin_logit)),1,function(x)x[1]/sum(x))
cor(scp_response,predict(m1_scpRatio_bin_logit,type="response"))^2
confint(m1_scpRatio_bin_logit)


plot_model(m1_scpRatio_bin_logit,type="re")
plot_model(m1_scpRatio_bin_logit,type="est",p.adjust="fdr")

pred_Ratioscp_TreatOrigin=ggpredict(m1_scpRatio_bin_logit,terms=c("Winter.treatment","Origin"))
pred_Ratioscp_TreatOrigin
plot(pred_Ratioscp_TreatOrigin)

pd1=position_dodge(0.6)
ratioSCP_plottable=scp_LGHMG_ratio[c("Winter.treatment","Origin","proportion_LGHMG")]
ratioSCP_plot=ggplot(data=pred_Ratioscp_TreatOrigin,aes(x=x,y=predicted)) +
  geom_violin(data=ratioSCP_plottable,aes(x=factor(Winter.treatment,levels=c("warm","cold")),y=proportion_LGHMG,color=Origin,fill=Origin),alpha=0.15,position=pd1,linewidth=0.2)+
  geom_point(data=ratioSCP_plottable,aes(x=factor(Winter.treatment,levels=c("warm","cold")),y=proportion_LGHMG,color=Origin,fill=Origin,shape=Origin),alpha=0.5,position=position_jitterdodge(jitter.width = 0.2,dodge.width = 0.6),size=0.5)+
  geom_point(data=pred_Ratioscp_TreatOrigin,aes(x=factor(x,levels=c("warm","cold")),y=predicted,shape=factor(group,levels=c("core","edge"))),position=pd1,size=1.3)+
  geom_errorbar(data=pred_Ratioscp_TreatOrigin,aes(ymin=conf.low,ymax=conf.high,group=factor(group,levels=c("core","edge"))),width=0.1,position=pd1,linewidth=.35)+
  scale_color_manual(values=c("#f58231", "#0082c8"))+
  scale_fill_manual(values=c("#f58231", "#0082c8"))+
  ylab(label="Proportion in LG SCP")+
  xlab(label="Winter treatment")+
  theme_bw()+
  theme(axis.text.x=element_text(size=5),axis.text.y=element_text(size=5))+
  theme(axis.title=element_text(size=6))+
  theme(axis.title.y=element_text(margin=margin(r=4)),axis.title.x=element_text(margin=margin(t=4)))+
  theme(legend.position="none")
ratioSCP_plot

pdf("output/Figures/SCPratio_final.pdf",width=1.8,height=2)
ratioSCP_plot
dev.off()
