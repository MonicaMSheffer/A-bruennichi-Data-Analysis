###This script analyzes the metabolomic profiles of Argiope bruennichi spiderlings following overwintering in a reciprocal common garden experiment

#load libraries
library(ggplot2)
library(tidyr)
library(ggeffects)
library(sjPlot)
library(bbmle)
library(ggpubr)

theme_set(theme_bw())

dat<-read.csv("data/metabolomics/abundance_2.csv",header=TRUE,row.names = 1)
meta<-read.csv("data/metabolomics/sample_meta.csv",header=T,row.names = 1)
compound<-read.csv("data/metabolomics/compound.csv",header=T)

#View(dat)
#View(meta)
#View(compound)

colnames(dat)<-compound$english_name

meta$Winter.treatment[meta$Winter.treatment=="Estonia"]="cold"
meta$Winter.treatment[meta$Winter.treatment=="France"]="warm"

meta$Origin[meta$Origin=="Estonia"]="edge"
meta$Origin[meta$Origin=="France"]="core"

meta_compound=cbind(meta,dat)
#only including post-winter egg sacs
post_only=meta_compound[meta_compound$Opening.round=="post",]

#from prior explorations, sample 191111_CA_35 seems to be an outlier - exclude and see if that helps things
post_only=post_only[post_only$TubeID!="CA_35",]

#get sample sizes
aggregate(Spider~Winter.treatment*Origin,data=post_only,NROW)
#get distribution of input biomass
hist(post_only$Mass,breaks=10)

#make a PDF of boxplots to look at all compounds
pdf("output/metabolomics_overview_allCompounds_boxplots.pdf",width=7,height=5)
for (i in 22:length(colnames(post_only))) {
  print(post_only%>%
          ggplot(aes(x=.[[1]],y=.[[i]],fill=Origin))+
          geom_boxplot(aes(fill=Origin))+
          geom_point(position=position_dodge(width=0.75),aes(group=Origin))+
          scale_fill_manual(values=c("#f58231","#0082c8"))+
          xlab(label="Winter treatment")+
          ylab(label=colnames(post_only[i])))
}
dev.off()


#write a loop to run the model (compound ~ Winter.treatment*Origin+Extracted.by) for each compound, and make a table of the results
compounds<-colnames(post_only[22:length(colnames(post_only))])
explanatoryvariables=c("Winter.treatment","Origin","Extracted.by","Winter.treatment:Origin")
stats= array(data=NA, dim=c(length(compounds),4,5),
             dimnames=list(compounds, c("Winter.treatment","Origin","Extracted.by","Winter.treatment:Origin"), c("coef","SE","Sum Sq","F value","Pr(>F)")) ) 

pdf("output/metabolomics_model_plots_allCompounds.pdf",width=7,height=5)
for(compound.i in 22:length(colnames(post_only))){
  mod.lm <- lm(post_only[,compound.i]~Winter.treatment*Origin+Extracted.by,data=post_only)
  print(plot_model(mod.lm,title="",type="pred",terms=c("Winter.treatment","Origin","Extracted.by"))+
          scale_color_manual(values=c("#f58231","#0082c8"))+
          ylab(label=colnames(post_only[compound.i])))
  stats[compound.i-21,,]=cbind(coef(summary(mod.lm))[2:nrow(coef(summary(mod.lm))),1:2],
                               as.matrix(anova(mod.lm))[1:4,c("Sum Sq","F value","Pr(>F)")])
}
dev.off()

#write a matrix of extracted statistics for each one
stat.mat= rbind(stats[1,,],stats[2,,],stats[3,,],stats[4,,],stats[5,,],stats[6,,],stats[7,,],stats[8,,],stats[9,,],stats[10,,],stats[11,,],stats[12,,],stats[13,,],stats[14,,],stats[15,,],stats[16,,],stats[17,,],stats[18,,],stats[19,,],stats[20,,],stats[21,,],stats[22,,],stats[23,,],stats[24,,],stats[25,,],stats[26,,],stats[27,,],stats[28,,],stats[29,,],stats[30,,],stats[31,,],stats[32,,],stats[33,,],stats[34,,],stats[35,,],stats[36,,],stats[37,,],stats[38,,],stats[39,,],stats[40,,],stats[41,,],stats[42,,],stats[43,,],stats[44,,],stats[45,,],stats[46,,],stats[47,,],stats[48,,],stats[49,,] )
stat.mat= as.data.frame(stat.mat)
stat.mat=cbind(rep(compounds, each=4), stat.mat)
stat.mat$sig=""
stat.mat$sig[stat.mat$`Pr(>F)`<0.01]="."
stat.mat$sig[stat.mat$`Pr(>F)`<0.05]="*"
stat.mat$sig[stat.mat$`Pr(>F)`<0.01]="**"
stat.mat$sig[stat.mat$`Pr(>F)`<0.001]="***"

stat.mat$var=rep(explanatoryvariables,times=49)

write.csv(stat.mat,"output/metabolomics_all_model_table.csv")


#from above, make a list of compounds that showed a significant effect of either winter treatment, origin, or their interaction
int_comp=c("Valine","Leucine","Isoleucine","Proline","Lysine","Citrate","Itaconate","Aconitate","3-Phosphoglycerate","myo-Inositol1")
necessary_cols=c("Winter.treatment","Origin","Extracted.by","Valine","Leucine","Isoleucine","Proline","Lysine","Citrate","Itaconate","Aconitate","3-Phosphoglycerate","myo-Inositol1")
post_int_comp<-subset(post_only,select=necessary_cols)
post_int_comp$Winter.treatment=factor(post_int_comp$Winter.treatment,levels=c("warm","cold"))

###plot the raw data and model estimates for all of the 'interesting' compounds
pd1=position_dodge(0.6)

pdf("output/metabolomics_pretty_modelPlots_interestingCompounds.pdf",width=7,height=5)
for (i in 4:length(colnames(post_int_comp))) {
  mod.lm <- lm(post_int_comp[,i]~Winter.treatment*Origin+Extracted.by,data=post_int_comp)
  pred.model<-ggpredict(mod.lm,terms=c("Winter.treatment","Origin"))

  print(ggplot(data=pred.model,aes(x=x,y=predicted)) +
          geom_violin(data=post_int_comp,aes(x=Winter.treatment,y=.data[[colnames(post_int_comp)[i]]],color=Origin,fill=Origin),alpha=0.1,position=pd1)+
          geom_point(data=post_int_comp,aes(x=Winter.treatment,y=.data[[colnames(post_int_comp)[i]]],color=Origin,fill=Origin,shape=Origin),alpha=0.5,position=position_jitterdodge(jitter.width = 0.2,dodge.width = 0.6))+
          geom_point(data=pred.model,aes(x=factor(x,levels=c("warm","cold")),y=predicted,shape=factor(group,levels=c("core","edge"))),position=pd1,size=3)+
          geom_errorbar(data=pred.model,aes(ymin=conf.low,ymax=conf.high,group=factor(group,levels=c("core","edge"))),width=0.1,position=pd1,linewidth=0.5)+
          scale_color_manual(values=c("#f58231", "#0082c8"))+
          scale_fill_manual(values=c("#f58231", "#0082c8"))+
          ylab(label=paste("Relative",colnames(post_int_comp[i]),"Concentration",sep=" "))+
          xlab(label="Winter treatment")+
          theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=10))+
          theme(axis.title=element_text(size=13,face="bold"))+
          theme(axis.title.y=element_text(margin=margin(r=4)),axis.title.x=element_text(margin=margin(t=4)))+
          theme(legend.text=element_text(size=12),legend.title=element_text(size=13),legend.position="bottom"))

  }
dev.off()

###single compound plots for talks/paper####
pd1=position_dodge(0.6)

mod.lm <- lm(`myo-Inositol1`~Winter.treatment*Origin+Extracted.by,data=post_int_comp)
pred.model<-ggpredict(mod.lm,terms=c("Winter.treatment","Origin"))

print(ggplot(data=pred.model,aes(x=x,y=predicted)) +
          geom_point(data=pred.model,aes(x=factor(x,levels=c("warm","cold")),y=predicted,shape=factor(group,levels=c("core","edge")),color=factor(group,levels=c("core","edge"))),position=pd1,size=3)+
          geom_errorbar(data=pred.model,aes(ymin=conf.low,ymax=conf.high,group=factor(group,levels=c("core","edge")),color=factor(group,levels=c("core","edge"))),width=0.1,position=pd1,linewidth=0.5)+
        geom_point(data=post_int_comp,aes(x=Winter.treatment,y=`myo-Inositol1`,color=Origin,shape=Origin),position=position_jitterdodge(jitter.width = 0.05,dodge.width = 0.6))+
        geom_boxplot(data=post_only,aes(x=Winter.treatment,y=`myo-Inositol1`,fill=Origin),alpha=0.4,width=0.6)+
          scale_color_manual(values=c("#f58231", "#0082c8"))+
          scale_fill_manual(values=c("#f58231", "#0082c8"))+
          ylab(label="Relative myo-inositol concentration")+
          xlab(label="Winter treatment")+
          theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=10))+
          theme(axis.title=element_text(size=13,face="bold"))+
          theme(axis.title.y=element_text(margin=margin(r=4)),axis.title.x=element_text(margin=margin(t=4)))+
          theme(legend.text=element_text(size=12),legend.title=element_text(size=13),legend.position="bottom"))


mod.lm <- lm(Proline~Winter.treatment*Origin+Extracted.by,data=post_int_comp)
pred.model<-ggpredict(mod.lm,terms=c("Winter.treatment","Origin"))

print(ggplot(data=pred.model,aes(x=x,y=predicted)) +
        geom_point(data=pred.model,aes(x=factor(x,levels=c("warm","cold")),y=predicted,shape=factor(group,levels=c("core","edge")),color=factor(group,levels=c("core","edge"))),position=pd1,size=3)+
        geom_errorbar(data=pred.model,aes(ymin=conf.low,ymax=conf.high,group=factor(group,levels=c("core","edge")),color=factor(group,levels=c("core","edge"))),width=0.1,position=pd1,linewidth=0.5)+
        geom_point(data=post_int_comp,aes(x=Winter.treatment,y=Proline,color=Origin,shape=Origin),position=position_jitterdodge(jitter.width = 0.05,dodge.width = 0.6))+
        scale_color_manual(values=c("#f58231", "#0082c8"))+
        scale_fill_manual(values=c("#f58231", "#0082c8"))+
        ylab(label="Relative proline concentration")+
        xlab(label="Winter treatment")+
        theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=10))+
        theme(axis.title=element_text(size=13,face="bold"))+
        theme(axis.title.y=element_text(margin=margin(r=4)),axis.title.x=element_text(margin=margin(t=4)))+
        theme(legend.text=element_text(size=12),legend.title=element_text(size=13),legend.position="bottom"))


pd1=position_dodge(0.6)

mod.lm <- lm(`myo-Inositol1`~Winter.treatment*Origin+Extracted.by,data=post_int_comp)
pred.model<-ggpredict(mod.lm,terms=c("Winter.treatment","Origin"))
myo_box<-ggplot(data=post_only,aes(y=`myo-Inositol1`,x=factor(Winter.treatment,levels=c("warm","cold")),fill=Origin))+
          geom_boxplot(data=post_only,aes(fill=Origin),alpha=0.15,width=0.3,linewidth=0.2)+
  geom_jitter(data=post_only,position=position_jitterdodge(dodge.width=0.3,jitter.width=0.1),aes(group=Origin,color=Origin,shape=Origin),size=1,alpha=0.5)+
          scale_fill_manual(values=c("#f58231","#0082c8"))+
          scale_color_manual(values=c("#f58231","#0082c8"))+
          xlab(label="Winter treatment")+
          ylab(label="Relative myo-inositol concentration")+ylim(0,1.75)+
  theme_bw()+
  theme(axis.text.x=element_text(size=5),axis.text.y=element_text(size=5))+
  theme(axis.title=element_text(size=6))+
  theme(axis.title.y=element_text(margin=margin(r=4)),axis.title.x=element_text(margin=margin(t=4)))+ 
  theme(legend.position="none")
myo_box

pdf("output/Figures/myo-inositol_final.pdf",width=1.8,height=2)
myo_box
dev.off()

mod.lm <- lm(Proline~Winter.treatment*Origin+Extracted.by,data=post_int_comp)
pred.model<-ggpredict(mod.lm,terms=c("Winter.treatment","Origin"))
pro_box<-ggplot(data=post_only,aes(y=Proline,x=factor(Winter.treatment,levels=c("warm","cold")),fill=Origin))+
  geom_boxplot(data=post_only,aes(fill=Origin),alpha=0.15,width=0.3,linewidth=0.2,outliers=FALSE)+
  geom_jitter(data=post_only,position=position_jitterdodge(dodge.width=0.3,jitter.width=0.1),aes(group=Origin,color=Origin,shape=Origin),size=1,alpha=0.5)+
  scale_fill_manual(values=c("#f58231","#0082c8"))+
  scale_color_manual(values=c("#f58231","#0082c8"))+
  xlab(label="Winter treatment")+
  ylab(label="Relative proline concentration")+ylim(0,1.75)+
  theme_bw()+
  theme(axis.text.x=element_text(size=5),axis.text.y=element_text(size=5))+
  theme(axis.title=element_text(size=6))+
  theme(axis.title.y=element_text(margin=margin(r=4)),axis.title.x=element_text(margin=margin(t=4)))+ 
  theme(legend.position="none")
pro_box

pdf("output/Figures/proline_final.pdf",width=1.8,height=2)
pro_box
dev.off()


###get confidence intervals for the model estimates for the 10 main compounds
#myo-inositol
mod.lm<-lm(`myo-Inositol1`~Winter.treatment*Origin+Extracted.by,data=post_int_comp)
confint(mod.lm)
summary(mod.lm)
stats_myoInositol1<-cbind(confint(mod.lm),coef(summary(mod.lm))[,1])
#copy rows 2-5 and all columns to the clipboard to paste into a table, without column or row names
install.packages("clipr")
library(clipr)
write_clip(stats_myoInositol1[2:5,],col.names=FALSE,row.names=FALSE)

#Valine
mod.lm<-lm(`Valine`~Winter.treatment*Origin+Extracted.by,data=post_int_comp)
stats_Valine<-cbind(confint(mod.lm),coef(summary(mod.lm))[,1])
#copy rows 2-5 and all columns to the clipboard to paste into a table, without column or row names
write_clip(stats_Valine[2:5,],col.names=FALSE,row.names=FALSE)

#Leucine
mod.lm<-lm(`Leucine`~Winter.treatment*Origin+Extracted.by,data=post_int_comp)
stats_Leucine<-cbind(confint(mod.lm),coef(summary(mod.lm))[,1])
#copy rows 2-5 and all columns to the clipboard to paste into a table, without column or row names
write_clip(stats_Leucine[2:5,],col.names=FALSE,row.names=FALSE)

#Isoleucine
mod.lm<-lm(`Isoleucine`~Winter.treatment*Origin+Extracted.by,data=post_int_comp)
stats_Isoleucine<-cbind(confint(mod.lm),coef(summary(mod.lm))[,1])
#copy rows 2-5 and all columns to the clipboard to paste into a table, without column or row names
write_clip(stats_Isoleucine[2:5,],col.names=FALSE,row.names=FALSE)

#Proline
mod.lm<-lm(`Proline`~Winter.treatment*Origin+Extracted.by,data=post_int_comp)
stats_Proline<-cbind(confint(mod.lm),coef(summary(mod.lm))[,1])
#copy rows 2-5 and all columns to the clipboard to paste into a table, without column or row names
write_clip(stats_Proline[2:5,],col.names=FALSE,row.names=FALSE)

#Lysine
mod.lm<-lm(`Lysine`~Winter.treatment*Origin+Extracted.by,data=post_int_comp)
stats_Lysine<-cbind(confint(mod.lm),coef(summary(mod.lm))[,1])
#copy rows 2-5 and all columns to the clipboard to paste into a table, without column or row names
write_clip(stats_Lysine[2:5,],col.names=FALSE,row.names=FALSE)

#Citrate
mod.lm<-lm(`Citrate`~Winter.treatment*Origin+Extracted.by,data=post_int_comp)
stats_Citrate<-cbind(confint(mod.lm),coef(summary(mod.lm))[,1])
#copy rows 2-5 and all columns to the clipboard to paste into a table, without column or row names
write_clip(stats_Citrate[2:5,],col.names=FALSE,row.names=FALSE)

#Itaconate
mod.lm<-lm(`Itaconate`~Winter.treatment*Origin+Extracted.by,data=post_int_comp)
stats_Itaconate<-cbind(confint(mod.lm),coef(summary(mod.lm))[,1])
#copy rows 2-5 and all columns to the clipboard to paste into a table, without column or row names
write_clip(stats_Itaconate[2:5,],col.names=FALSE,row.names=FALSE)

#Aconitate
mod.lm<-lm(`Aconitate`~Winter.treatment*Origin+Extracted.by,data=post_int_comp)
stats_Aconitate<-cbind(confint(mod.lm),coef(summary(mod.lm))[,1])
#copy rows 2-5 and all columns to the clipboard to paste into a table, without column or row names
write_clip(stats_Aconitate[2:5,],col.names=FALSE,row.names=FALSE)

#3-Phosphoglycerate
mod.lm<-lm(`3-Phosphoglycerate`~Winter.treatment*Origin+Extracted.by,data=post_int_comp)
stats_3Phosphoglycerate<-cbind(confint(mod.lm),coef(summary(mod.lm))[,1])
#copy rows 2-5 and all columns to the clipboard to paste into a table, without column or row names
write_clip(stats_3Phosphoglycerate[2:5,],col.names=FALSE,row.names=FALSE)
