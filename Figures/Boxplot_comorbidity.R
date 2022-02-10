## final plot of allergy prediction paper

rm(list=ls())
## ===============================================
## 1. boxplots of different groups in PIAMA
## ===============================================
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsci)

## load data
setwd("~/Documents/Projects/PIAMA/MIcompany/allergy prediction/")
## nasal data
load("./data/beta_nasal_sig.Rdata")
pheno<-read.csv("./data/allergy_newid.csv")
## select samples used in prediction
load("./data/sample_id_used_prediction.Rdata")
## load phenotype
df.pheno<-pheno[match(sample.id,pheno$newid),]
## load cpg id
cpgs<-c("cg20372759","cg01870976","cg24224501")

df<-t(beta_nasal_sig[match(cpgs,rownames(beta_nasal_sig)),sample.id])
df.plot<-as.data.frame(df)
M_matrix<-t(df)

PHENO<-pheno[match(rownames(df.plot),pheno$newid),]
PHENO$newid<-as.character(PHENO$newid)
rownames(PHENO)<-PHENO$newid

## ======================================
##  a. boxplot allergy non-allergy

data1<-cbind(PHENO$Allergy,df.plot)
colnames(data1)[1]<-"Allergy"
data1$Allergy<-as.factor(data1$Allergy)

## 3 CpG
df1<-gather(data1,"cg24224501","cg20372759","cg01870976",key = "cpg",value = "beta_value")
df1$cpg<-as.factor(df1$cpg)
levels(df1$Allergy)<-c("no","yes")

## save as width = 5, height=4
pdf("./replication+plots/PIAMA_3cpg.pdf",width = 5,height = 5)
ggplot(df1, aes(x=cpg,y=beta_value,fill=Allergy))+
  geom_boxplot()+
  xlab("")+ylab("DNA methylation (beta value)")+ 
  ylim(0.1,1) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
  theme(text = element_text(size=15,face="bold"),axis.text.y = element_text(size=15,face="bold"),
        axis.text.x = element_text(size=15,face="bold",angle = 30,vjust = 0.5,hjust = 0.5))+
  scale_color_nejm()+
  ggtitle("PIAMA")
dev.off()

## =====================================================
## b. comorbidity

PHENO$N_disease<-PHENO$Eczema + PHENO$Rhinitis + PHENO$Asthma
PHENO$allergy1<-0
PHENO$allergy1[which(PHENO$N_disease==1 & PHENO$Allergy==1)]<-1
PHENO$allergy1[which(PHENO$N_disease>1 & PHENO$Allergy==1)]<-2

##  strafity no allergy by no diseaes only, no IgE only, no disease and no IgE
PHENO$allergy1[which(PHENO$N_disease==0 & PHENO$sentization==1)]<-3
PHENO$allergy1[which(PHENO$N_disease>0 & PHENO$sentization==0)]<-4
PHENO$allergy1[which(PHENO$N_disease==0 & PHENO$sentization==0)]<-5

data1<-cbind(PHENO$allergy1,df.plot)
colnames(data1)[1]<-"Allergy"
data1$Allergy<-as.factor(data1$Allergy)
levels(data1$Allergy)<-c("IgE+one_symptom","IgE+2or3_symptom","IgE+symptom-","IgE-symptom+","IgE-symptpm-")

df1<-gather(data1,"cg24224501","cg20372759","cg01870976",key = "cpg",value = "beta_value")
df1$cpg<-as.factor(df1$cpg)

## save
pdf("./replication+plots/PIAMA_comorbidity.pdf",width = 12,height = 4)
ggplot(df1,aes(x=Allergy,y=beta_value,fill=Allergy))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width=0.3,alpha=0.4,colour="black",shape=21,aes(colour=Allergy))+
  stat_compare_means(comparisons = list(c("IgE+one_symptom","IgE+2or3_symptom"),
                                        c("IgE+one_symptom","IgE+symptom-"),
                                        c("IgE+2or3_symptom","IgE+symptom-")),
                     label = "p.signif",method = "t.test")+
  facet_wrap(~cpg)+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
  theme(text = element_text(size=15,face="bold"),axis.text.y = element_text(size=15,face="bold"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  scale_y_continuous(breaks=c(0.4,0.6,0.8,1.0))+
  ylab("DNA methylation (beta value)")+ 
  scale_fill_manual(values = c("#F39B7FB2","#E64B35B2","#4DBBD5B2","#00A087B2","#3C5488B2"))+
  scale_color_manual(values = c("#F39B7FB2","#E64B35B2","#4DBBD5B2","#00A087B2","#3C5488B2"))
dev.off()


## =====================================================
## c. sym-asym

PHENO$allergy2<-NA
PHENO$allergy2[which(PHENO$sentization==1 & PHENO$N_disease>0)]<-1
PHENO$allergy2[which(PHENO$sentization==1 & PHENO$N_disease==0)]<-2
PHENO$allergy2[which(PHENO$sentization==0 & PHENO$N_disease>0)]<-3
PHENO$allergy2[which(PHENO$sentization==0 & PHENO$N_disease==0)]<-4

data1<-cbind(PHENO$allergy2,df.plot)
colnames(data1)[1]<-"Allergy"
data1$Allergy<-as.factor(data1$Allergy)
levels(data1$Allergy)<-c("IgE+symptom+","IgE+symptom-","IgE-symptom+","IgE-symptom-")

df1<-gather(data1,"cg24224501","cg20372759","cg01870976",key = "cpg",value = "beta_value")
df1$cpg<-as.factor(df1$cpg)

## save plot
pdf("./replication+plots/PIAMA_asymptomatic.pdf",width = 11,height = 4)
ggplot(df1,aes(x=Allergy,y=beta_value,fill=Allergy))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width=0.3,alpha=0.4,colour="black",shape=21,aes(colour=Allergy))+
  stat_compare_means(comparisons = list(c("IgE+symptom+","IgE+symptom-"),
                                        c("IgE-symptom+","IgE-symptom-"),
                                        c("IgE+symptom+","IgE-symptom+"),
                                        c("IgE+symptom-","IgE-symptom-")),
                     label = "p.signif",method = "t.test")+
  facet_wrap(~cpg)+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
  theme(text = element_text(size=15,face="bold"),axis.text.y = element_text(size=15,face="bold"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  scale_y_continuous(breaks=c(0.4,0.6,0.8,1.0))+
  ylab("DNA methylation (beta value)")+ 
  scale_fill_manual(values = c("#DC0000B2","#4DBBD5B2","#00A087B2","#3C5488B2"))+
  scale_color_manual(values = c("#DC0000B2","#4DBBD5B2","#00A087B2","#3C5488B2"))
dev.off()


