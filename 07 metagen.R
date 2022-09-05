
library(openxlsx)
library(ggplot2)
library(ggpubr)
setwd('E:/study2/microbiome/study_new/total/9Data Review')
metadata = read.xlsx("Data/data_total.xlsx", sheet=1,rowNames = F)

metadata2 = read.table("Data/metagen/ID.txt", header = T,sep = "\t")
#metadata2<-subset(metadata2,dyslipimedia_24=="Control")

ID<-as.character(metadata2$number)
metadata_new<-subset(metadata,number%in%ID)


#一、 species -----------------------------------------------------------------


dat_species<-read.csv("Data/metagen/Species.csv")

rownames(dat_species)<-dat_species$Species
dat_species<-dat_species[,-c(1:7)]

abundance = 0.01
idx = rowMeans(dat_species) > abundance
filtered_data = dat_species[idx,]

detectionRate=10

detect<-function(x){
  length(which(x > 0))*100/length(x)
}

DR<-apply(filtered_data, 1, detect)
idx3= apply(filtered_data, 1, detect) >detectionRate

filtered_data2<-filtered_data[idx3,]

species_name<-rownames(filtered_data2)
write.csv(filtered_data2,file = "Data/metagen/Species_filtered.csv",row.names = T)

filtered_data2<-as.data.frame(t(filtered_data2))
filtered_data2$ID<-rownames(filtered_data2)
#rownames(metadata2)<-metadata2$SampleID
dat_species2<-merge(metadata2,filtered_data2,by.x = "SampleID",by.y = "ID")
dat_species2<-dat_species2[,-c(1,3,4)]
dat_species2<-merge(metadata_new,dat_species2,by="number")
write.csv(dat_species2,file = "Data/metagen/dat_species.csv",row.names = T)



#1. linear regression -------------------------------------------------------
dir.create("Result/metagen")
library(fdrtool)
x_var<-c("TG_24","CHOL_24","TG_32","CHOL_32","delta_TG","delta_CHOL")
y_var<-species_name

cov<-c("pre_BMI","Parity","Age","gesweek_24")
dir.create("metagen")
dir.create("Result/metagen")

# linear ------------------------------------------------------------------
data_combine<-dat_species2
data_combine[species_name]<-log(data_combine[species_name]+0.001)

# unadjusted ----------------------------------------------------------------


###unadjusted
result<-data.frame()
for (i in x_var){
  for(j in y_var){
    # xnam_i1<-paste(i,j,sep="*")
    # xnam_i<-paste(j,"pre_BMI","Parity","Age","gesweek_24",sep="+")
    fit_i<-lm(as.formula(paste(i,"~",j)),data=data_combine)
    p_i<-tryCatch({summary(fit_i)$coefficients[2,4]},error=function(e){"NA"})
    b_i<-tryCatch({coef(fit_i)[2]},error=function(e){"NA"})
    blow_i<-tryCatch({confint(fit_i)[2,1]},error=function(e){"NA"})
    bhigh_i<-tryCatch({confint(fit_i)[2,2]},error=function(e){"NA"})
    result_i<-c(i,j,b_i,blow_i,bhigh_i,p_i)  
    result_i<-as.data.frame(result_i) 
    result_i<-as.data.frame(t(result_i)) 
    result<-rbind(result,result_i) 
  }
}

names(result)<-c("outcome","species","beta","lower","upper","pvalue")   

result_new<-data.frame()
for (i in x_var) {
  result2<-subset(result,outcome==i)
  p<-as.numeric(as.character(result2$pvalue))
  fdr=fdrtool(p,statistic="pvalue",plot = F)
  qval<-data.frame(fdr$qval)
  p2<-p.adjust(p,method = "fdr",length(p))
  p2<-data.frame(p2)
  result2<-cbind(result2,qval,p2)
  result_new<-rbind(result_new,result2)
  
}

group<-rep("Unadjusted",times=nrow(result_new))

result1<-data.frame(result_new,group)



# adjusted ----------------------------------------------------------------



###adjusted
result<-data.frame()
for (i in x_var){
  for(j in y_var){
    # xnam_i1<-paste(i,j,sep="*")
    xnam_i<-paste(j,"pre_BMI","Parity","Age","gesweek_24",sep="+")
    fit_i<-lm(as.formula(paste(i,"~",xnam_i)),data=data_combine)
    p_i<-tryCatch({summary(fit_i)$coefficients[2,4]},error=function(e){"NA"})
    b_i<-tryCatch({coef(fit_i)[2]},error=function(e){"NA"})
    blow_i<-tryCatch({confint(fit_i)[2,1]},error=function(e){"NA"})
    bhigh_i<-tryCatch({confint(fit_i)[2,2]},error=function(e){"NA"})
    result_i<-c(i,j,b_i,blow_i,bhigh_i,p_i)  
    result_i<-as.data.frame(result_i) 
    result_i<-as.data.frame(t(result_i)) 
    result<-rbind(result,result_i) 
  }
}


names(result)<-c("outcome","species","beta","lower","upper","pvalue")   

result_new<-data.frame()
for (i in x_var) {
  result2<-subset(result,outcome==i)
  p<-as.numeric(as.character(result2$pvalue))
  fdr=fdrtool(p,statistic="pvalue",plot = F)
  qval<-data.frame(fdr$qval)
  p2<-p.adjust(p,method = "fdr",length(p))
  p2<-data.frame(p2)
  result2<-cbind(result2,qval,p2)
  result_new<-rbind(result_new,result2)
  
}


group<-rep("Adjusted",times=nrow(result_new))

result2<-data.frame(result_new,group)

result_full<-rbind(result1,result2)
write.csv(result_full,file="Result/metagen/linear_all.csv",row.names = F) 


#2. binary ------------------------------------------------------------------


# unadjusted -------------------------------------------------------------


###binary
x_var<-c("Hyperlipid_24","TG_24_High","CHOL_24_High",
         "Hyperlipid_32","TG_32_High","CHOL_32_High")


result<-data.frame()
for (i in x_var){
  for(j in y_var){
    # xnam_i1<-paste(i,j,sep="*")
    #xnam_i<-paste(j,"pre_BMI","Parity","Age","gesweek_24",sep="+")
    fit_i<-glm(as.formula(paste(i,"~",j)),family = "binomial",data=data_combine)
    p_i<-tryCatch({summary(fit_i)$coefficients[2,4]},error=function(e){"NA"})
    b_i<-tryCatch({exp(coef(fit_i)[2])},error=function(e){"NA"})
    blow_i<-tryCatch({exp(confint(fit_i)[2,1])},error=function(e){"NA"})
    bhigh_i<-tryCatch({exp(confint(fit_i)[2,2])},error=function(e){"NA"})
    result_i<-c(i,j,b_i,blow_i,bhigh_i,p_i)
    result_i<-as.data.frame(result_i)
    result_i<-as.data.frame(t(result_i))
    result<-rbind(result,result_i)
  }
}

names(result)<-c("outcome","species","OR","low","upper","pvalueֵ")

result_new<-data.frame()
for (i in x_var) {
  result2<-subset(result,outcome==i)
  p<-as.numeric(as.character(result2$pvalue))
  fdr=fdrtool(p,statistic="pvalue",plot = F)
  qval<-data.frame(fdr$qval)
  p2<-p.adjust(p,method = "fdr",length(p))
  p2<-data.frame(p2)
  result2<-cbind(result2,qval,p2)
  result_new<-rbind(result_new,result2)
  
}
group<-rep("Unadjusted",times=nrow(result_new))

result1<-data.frame(result_new,group)


# adjusted ----------------------------------------------------------------


result<-data.frame()
for (i in x_var){
  for(j in y_var){
    # xnam_i1<-paste(i,j,sep="*")
    xnam_i<-paste(j,"pre_BMI","Parity","Age","gesweek_24",sep="+")
    fit_i<-glm(as.formula(paste(i,"~",xnam_i)),family = "binomial",data=data_combine)
    p_i<-tryCatch({summary(fit_i)$coefficients[2,4]},error=function(e){"NA"})
    b_i<-tryCatch({exp(coef(fit_i)[2])},error=function(e){"NA"})
    blow_i<-tryCatch({exp(confint(fit_i)[2,1])},error=function(e){"NA"})
    bhigh_i<-tryCatch({exp(confint(fit_i)[2,2])},error=function(e){"NA"})
    result_i<-c(i,j,b_i,blow_i,bhigh_i,p_i)
    result_i<-as.data.frame(result_i)
    result_i<-as.data.frame(t(result_i))
    result<-rbind(result,result_i)
  }
}

names(result)<-c("outcome","species","OR","low","upper","pvalueֵ")


result_new<-data.frame()
for (i in x_var) {
  result2<-subset(result,outcome==i)
  p<-as.numeric(as.character(result2$pvalue))
  fdr=fdrtool(p,statistic="pvalue",plot = F)
  qval<-data.frame(fdr$qval)
  p2<-p.adjust(p,method = "fdr",length(p))
  p2<-data.frame(p2)
  result2<-cbind(result2,qval,p2)
  result_new<-rbind(result_new,result2)
  
}


group<-rep("Adjusted",times=nrow(result_new))

result2<-data.frame(result_new,group)

result_full<-rbind(result1,result2)

write.csv(result_full,file="Result/metagen/linear_all_group.csv",row.names = F)




#3. multiple category ----------------------------------------------------------

library(DALEX)
library(iBreakDown)
library(nnet)
library(questionr)
library(dplyr)
outcome<-"dyslipimedia_32"
#data<-na.omit(data_combine[,c(feature,outcome)])
#dir.create(outcome)

library(fdrtool)
coef<-data.frame()
p<-data.frame()
result<-data.frame()

data_combine$dyslipimedia_32<-ordered(data_combine$dyslipimedia_32,levels=c("Control", "Hypercholesterolemia",
                                                                    "Hyperlipidemia", "Hypertriglyceridemia"))
data_combine$delta_TG_group1<-ordered(data_combine$delta_TG_group1,levels=c("Group1","Group2","Group3","Group4"))
data_combine$delta_CHOL_group1<-ordered(data_combine$delta_TG_group1,levels=c("Group1","Group2","Group3","Group4"))
data_combine$delta_TG_group2<-ordered(data_combine$delta_TG_group1,levels=c("Group1","Group2","Group3"))
data_combine$delta_CHOL_group2<-ordered(data_combine$delta_TG_group1,levels=c("Group1","Group2","Group3"))



#3.1 four category ------------------------------------------------------------


# unadjusted --------------------------------------------------------------



x_var<-c("dyslipimedia_32","delta_TG_group1","delta_CHOL_group1")
result<-data.frame()
for (i in x_var){
  for(j in y_var){
    # xnam_i1<-paste(i,j,sep="*")
    #xnam_i<-paste(j,"pre_BMI","Parity","Age","gesweek_24",sep="+")
    test<-multinom(as.formula(paste(i,"~",j,sep = "")),data=data_combine)
    r<-coef(test)
    se<-summary(test)$standard.errors
    z <- summary(test)$coefficients/summary(test)$standard.errors
    p_i <- (1 - pnorm(abs(z), 0, 1)) * 2
    result_i<-c(r[,2],se[,2])
    result_i<-as.data.frame(result_i)
    result_i<-as.data.frame(t(result_i))
    #name_i<-names(result_i)[-c(1:2)]
    names(result_i)<-c("coef_2","coef_3","coef_4",
                       "se_2","se_3","se_4")
    for (k in 1:ncol(result_i)) {
      result_i[,k]<-as.numeric(result_i[,k])
    }
    
    ci_2<-paste0("(",round(exp(result_i$coef_2-1.96*result_i$se_2),3),", ",
                 round(exp(result_i$coef_2+1.96*result_i$se_2),3),")")
    ci_3<-paste0("(",round(exp(result_i$coef_3-1.96*result_i$se_3),3),", ",
                 round(exp(result_i$coef_3+1.96*result_i$se_3),3),")")
    ci_4<-paste0("(",round(exp(result_i$coef_4-1.96*result_i$se_4),3),", ",
                 round(exp(result_i$coef_4+1.96*result_i$se_4),3),")")
    
    
    ci<-c(ci_2,ci_3,ci_4)
    
    result_1i<-c(i,j,exp(r[,2]),ci,p_i[,2])
    result_1i<-as.data.frame(result_1i)
    result_1i<-as.data.frame(t(result_1i))
    result<-rbind(result,result_1i)
    
  }
}

names(result)<-c("outcome","genus","OR_2","OR_3","OR_4",
                 "CI_2","CI_3","CI_4",
                 "p_2","p_3","p_4")


result_new<-data.frame()
for (i in x_var) {
  resulT3<-subset(result,outcome==i)
  
  p<-as.numeric(as.character(resulT3$p_2))
  # fdr=fdrtool(p,statistic="pvalue",plot = F)
  # qval<-data.frame(fdr$qval)
  p2.adj<-p.adjust(p,method = "fdr",length(p))
  p2.adj<-data.frame(p2.adj)
  
  p<-as.numeric(as.character(resulT3$p_3))
  # fdr=fdrtool(p,statistic="pvalue",plot = F)
  # qval<-data.frame(fdr$qval)
  p3.adj<-p.adjust(p,method = "fdr",length(p))
  p3.adj<-data.frame(p3.adj)
  
  p<-as.numeric(as.character(resulT3$p_4))
  # fdr=fdrtool(p,statistic="pvalue",plot = F)
  # qval<-data.frame(fdr$qval)
  p4.adj<-p.adjust(p,method = "fdr",length(p))
  p4.adj<-data.frame(p4.adj)
  
  resulT3<-cbind(resulT3,p2.adj,p3.adj,p4.adj)
  result_new<-rbind(result_new,resulT3)
  
}
group<-rep("Unadjusted",times=nrow(result_new))

result1<-data.frame(result_new,group)


# adjusted ----------------------------------------------------------------



x_var<-c("dyslipimedia_32","delta_TG_group1","delta_CHOL_group1")
result<-data.frame()
for (i in x_var){
  for(j in y_var){
    # xnam_i1<-paste(i,j,sep="*")
    xnam_i<-paste(j,"pre_BMI","Parity","Age","gesweek_24",sep="+")
    test<-multinom(as.formula(paste(i,"~",xnam_i,sep = "")),data=data_combine)
    r<-coef(test)
    se<-summary(test)$standard.errors
    z <- summary(test)$coefficients/summary(test)$standard.errors
    p_i <- (1 - pnorm(abs(z), 0, 1)) * 2
    result_i<-c(r[,2],se[,2])
    result_i<-as.data.frame(result_i)
    result_i<-as.data.frame(t(result_i))
    #name_i<-names(result_i)[-c(1:2)]
    names(result_i)<-c("coef_2","coef_3","coef_4",
                       "se_2","se_3","se_4")
    for (k in 1:ncol(result_i)) {
      result_i[,k]<-as.numeric(result_i[,k])
    }
    
    ci_2<-paste0("(",round(exp(result_i$coef_2-1.96*result_i$se_2),3),", ",
                 round(exp(result_i$coef_2+1.96*result_i$se_2),3),")")
    ci_3<-paste0("(",round(exp(result_i$coef_3-1.96*result_i$se_3),3),", ",
                 round(exp(result_i$coef_3+1.96*result_i$se_3),3),")")
    ci_4<-paste0("(",round(exp(result_i$coef_4-1.96*result_i$se_4),3),", ",
                 round(exp(result_i$coef_4+1.96*result_i$se_4),3),")")
    
    
    ci<-c(ci_2,ci_3,ci_4)
    
    result_1i<-c(i,j,exp(r[,2]),ci,p_i[,2])
    result_1i<-as.data.frame(result_1i)
    result_1i<-as.data.frame(t(result_1i))
    result<-rbind(result,result_1i)
    
  }
}

names(result)<-c("outcome","genus","OR_2","OR_3","OR_4",
                 "CI_2","CI_3","CI_4",
                 "p_2","p_3","p_4")


result_new<-data.frame()
for (i in x_var) {
  resulT3<-subset(result,outcome==i)
  
  p<-as.numeric(as.character(resulT3$p_2))
  # fdr=fdrtool(p,statistic="pvalue",plot = F)
  # qval<-data.frame(fdr$qval)
  p2.adj<-p.adjust(p,method = "fdr",length(p))
  p2.adj<-data.frame(p2.adj)
  
  p<-as.numeric(as.character(resulT3$p_3))
  # fdr=fdrtool(p,statistic="pvalue",plot = F)
  # qval<-data.frame(fdr$qval)
  p3.adj<-p.adjust(p,method = "fdr",length(p))
  p3.adj<-data.frame(p3.adj)
  
  p<-as.numeric(as.character(resulT3$p_4))
  # fdr=fdrtool(p,statistic="pvalue",plot = F)
  # qval<-data.frame(fdr$qval)
  p4.adj<-p.adjust(p,method = "fdr",length(p))
  p4.adj<-data.frame(p4.adj)
  
  resulT3<-cbind(resulT3,p2.adj,p3.adj,p4.adj)
  result_new<-rbind(result_new,resulT3)
  
}
group<-rep("Adjusted",times=nrow(result_new))

result2<-data.frame(result_new,group)

result_full<-rbind(result1,result2)

write.csv(result_full,file="Result/metagen/linear_T2micro_T3_multiple_group1.csv",row.names = F)



#3.2 three category ------------------------------------------------------------


# unadjusted --------------------------------------------------------------



x_var<-c("delta_TG_group2","delta_CHOL_group2")
result<-data.frame()
for (i in x_var){
  for(j in y_var){
    # xnam_i1<-paste(i,j,sep="*")
    #xnam_i<-paste(j,"pre_BMI","Parity","Age","gesweek_24",sep="+")
    test<-multinom(as.formula(paste(i,"~",j,sep = "")),data=data_combine)
    r<-coef(test)
    se<-summary(test)$standard.errors
    z <- summary(test)$coefficients/summary(test)$standard.errors
    p_i <- (1 - pnorm(abs(z), 0, 1)) * 2
    result_i<-c(r[,2],se[,2])
    result_i<-as.data.frame(result_i)
    result_i<-as.data.frame(t(result_i))
    #name_i<-names(result_i)[-c(1:2)]
    names(result_i)<-c("coef_2","coef_3",
                       "se_2","se_3")
    for (k in 1:ncol(result_i)) {
      result_i[,k]<-as.numeric(result_i[,k])
    }
    
    ci_2<-paste0("(",round(exp(result_i$coef_2-1.96*result_i$se_2),3),", ",
                 round(exp(result_i$coef_2+1.96*result_i$se_2),3),")")
    ci_3<-paste0("(",round(exp(result_i$coef_3-1.96*result_i$se_3),3),", ",
                 round(exp(result_i$coef_3+1.96*result_i$se_3),3),")")
    
    ci<-c(ci_2,ci_3)
    
    result_1i<-c(i,j,exp(r[,2]),ci,p_i[,2])
    result_1i<-as.data.frame(result_1i)
    result_1i<-as.data.frame(t(result_1i))
    result<-rbind(result,result_1i)
    
  }
}

names(result)<-c("outcome","genus","OR_2","OR_3",
                 "CI_2","CI_3",
                 "p_2","p_3")


result_new<-data.frame()
for (i in x_var) {
  resulT3<-subset(result,outcome==i)
  
  p<-as.numeric(as.character(resulT3$p_2))
  # fdr=fdrtool(p,statistic="pvalue",plot = F)
  # qval<-data.frame(fdr$qval)
  p2.adj<-p.adjust(p,method = "fdr",length(p))
  p2.adj<-data.frame(p2.adj)
  
  p<-as.numeric(as.character(resulT3$p_3))
  # fdr=fdrtool(p,statistic="pvalue",plot = F)
  # qval<-data.frame(fdr$qval)
  p3.adj<-p.adjust(p,method = "fdr",length(p))
  p3.adj<-data.frame(p3.adj)
  
  
  resulT3<-cbind(resulT3,p2.adj,p3.adj)
  result_new<-rbind(result_new,resulT3)
  
}
group<-rep("Unadjusted",times=nrow(result_new))

result1<-data.frame(result_new,group)


# adjusted ----------------------------------------------------------------



x_var<-c("delta_TG_group2","delta_CHOL_group2")
result<-data.frame()
for (i in x_var){
  for(j in y_var){
    # xnam_i1<-paste(i,j,sep="*")
    xnam_i<-paste(j,"pre_BMI","Parity","Age","gesweek_24",sep="+")
    test<-multinom(as.formula(paste(i,"~",xnam_i,sep = "")),data=data_combine)
    r<-coef(test)
    se<-summary(test)$standard.errors
    z <- summary(test)$coefficients/summary(test)$standard.errors
    p_i <- (1 - pnorm(abs(z), 0, 1)) * 2
    result_i<-c(r[,2],se[,2])
    result_i<-as.data.frame(result_i)
    result_i<-as.data.frame(t(result_i))
    #name_i<-names(result_i)[-c(1:2)]
    names(result_i)<-c("coef_2","coef_3",
                       "se_2","se_3")
    for (k in 1:ncol(result_i)) {
      result_i[,k]<-as.numeric(result_i[,k])
    }
    
    ci_2<-paste0("(",round(exp(result_i$coef_2-1.96*result_i$se_2),3),", ",
                 round(exp(result_i$coef_2+1.96*result_i$se_2),3),")")
    ci_3<-paste0("(",round(exp(result_i$coef_3-1.96*result_i$se_3),3),", ",
                 round(exp(result_i$coef_3+1.96*result_i$se_3),3),")")
    
    ci<-c(ci_2,ci_3)
    
    result_1i<-c(i,j,exp(r[,2]),ci,p_i[,2])
    result_1i<-as.data.frame(result_1i)
    result_1i<-as.data.frame(t(result_1i))
    result<-rbind(result,result_1i)
    
  }
}

names(result)<-c("outcome","genus","OR_2","OR_3",
                 "CI_2","CI_3",
                 "p_2","p_3")


result_new<-data.frame()
for (i in x_var) {
  resulT3<-subset(result,outcome==i)
  
  p<-as.numeric(as.character(resulT3$p_2))
  # fdr=fdrtool(p,statistic="pvalue",plot = F)
  # qval<-data.frame(fdr$qval)
  p2.adj<-p.adjust(p,method = "fdr",length(p))
  p2.adj<-data.frame(p2.adj)
  
  p<-as.numeric(as.character(resulT3$p_3))
  # fdr=fdrtool(p,statistic="pvalue",plot = F)
  # qval<-data.frame(fdr$qval)
  p3.adj<-p.adjust(p,method = "fdr",length(p))
  p3.adj<-data.frame(p3.adj)
  
  
  resulT3<-cbind(resulT3,p2.adj,p3.adj)
  result_new<-rbind(result_new,resulT3)
  
}
group<-rep("Adjusted",times=nrow(result_new))

result2<-data.frame(result_new,group)

result_full<-rbind(result1,result2)

write.csv(result_full,file="Result/metagen/linear_T2micro_T3_multiple_group2.csv",row.names = F)



#2. RF_prediction -----------------------------------------------------------

#2.1 RF_cate -------------------------------------------------------------


#data_rf_third<-subset(data_rf_third,number %in% data_rf_second$number)

Hypercholesterolemia_third<-data_combine$CHOL_32_High
Hypertriglyceridemia_third<-data_combine$TG_32_High
hyperlipid_third<-data_combine$Hyperlipid_32

data_rf_second<-data.frame(data_combine,Hypercholesterolemia_third,hyperlipid_third,Hypertriglyceridemia_third)
dir.create("Result/metagen/ML")
setwd("E:/study2/microbiome/study_new/total/9Data Review/Result/metagen/ML")
#随机森林分类
library(randomForest)
library(caret)
library(ggplot2)
library(bootstrap)
library(ROCR)
library(pROC)
set.seed(315)

data_rf_second$Hypercholesterolemia_third<-as.factor(data_rf_second$Hypercholesterolemia_third)
data_rf_second$Hypertriglyceridemia_third<-as.factor(data_rf_second$Hypertriglyceridemia_third)

data_rf_second$hyperlipid_third<-as.factor(data_rf_second$hyperlipid_third)

#test_name<-colnames(data_rf_second[,72:91])
test_name<-c("BUN", "CREA", "UA", "TP", "ALB", "CHOL", "TG", "ALT", "AST", 
             "ALP", "LDH", "GGT", "TBIL", "DBIL", "TBA", "FMN", "FBG", "RBP", 
             "cystain", "HCY")
test_name<-paste0(test_name,"_24")
#impute missing value
library(dplyr)
# data_rf_second<-data_rf_second %>% 
#   mutate_all(function(x){x[is.na(x)] <- mean(x)
#   x})

for (i in test_name) {
  
  data_rf_second[,i][is.na(data_rf_second[,i])==TRUE]<-mean(data_rf_second[,i],na.rm = T)
}


# 
# 2.1.1. species alone ------------------------------------------------------------
rf_classify<-function(i,j){
  # i<-"Hypercholesterolemia_third"
  dir.create(i)
  dir.create(paste0(i,"/species_alone"))
  varia<-c(i,species_name)
  m<-100
  result<-data.frame()
  plot_data<-data.frame()
  pred<-data.frame()
  set.seed(123)
  accura<-data.frame()
  result_m<-data.frame()
  imp<-data.frame()
  b=1
  ratio=j
  for (b in 1:m){
    
    Train <- createDataPartition(data_rf_second[,i], p=ratio, list=FALSE)
    training <- data_rf_second[ Train, varia]
    testing <- data_rf_second[ -Train, varia]
    
    rf = randomForest(as.formula(paste0(i,"~.") ), data=training,importance=TRUE, proximity=TRUE, ntree = 1000)
    
    pred_b=predict(rf,testing)
    df_b = table(observed = testing[,i], predict = pred_b)
    
    df_b<-as.data.frame(df_b)
    sensi_b<-df_b[4,3]/(df_b[2,3]+df_b[4,3]) 
    speci_b<-df_b[1,3]/(df_b[1,3]+df_b[3,3])
    
    pred_rf_b<-predict(rf,testing,type = "prob")
    pred_2b_b<-prediction(pred_rf_b[,2],testing[,i])
    perf_b<-performance(pred_2b_b,"tpr","fpr")
    auc_b<-performance(pred_2b_b,"auc")@y.values
    auc_b<-as.numeric(auc_b)
    
    predictions2 <- predict(rf, testing, type = "class")
    result_matrix<-confusionMatrix(predictions2, testing[,i])
    accura_b<-result_matrix$overall
    
    accura_b<-c(accura_b[1],accura_b[3],accura_b[4],auc_b,b)
    names(accura_b)<-c("Accuracy","Lower","Upper","AUC","Times")
    accura<-rbind(accura,accura_b)
    
    sensi<-result_matrix$byClass
    sensi<-as.data.frame(sensi)
    sensi<-as.data.frame(t(sensi))
    times<-rep(b,times=dim(sensi)[1])
    result_b<-data.frame(sensi,times)
    result_m<-rbind(result_m,result_b)
    
    pred_2b<-data.frame(pred_b,pred_rf_b)
    times<-rep(b,times=dim(pred_2b)[1])
    
    groundtruth<-testing[,i]
    pred_2b<-data.frame(pred_2b,groundtruth,times)
    pred<-rbind(pred,pred_2b)
    # rownames(pred)<-rownames(testing)
    #ROC曲线
    
    result_b<-c(b,sensi_b,speci_b,auc_b)
    result_b<-as.data.frame(result_b) 
    result_b<-as.data.frame(t(result_b)) 
    result<-rbind(result,result_b) 
    
    sensi_value_b<-perf_b@y.values[[1]]
    speci_value_b<-perf_b@x.values[[1]]
    plot_data_b<-data.frame(sensi_value_b,speci_value_b)
    group<-rep(b,times=dim(plot_data_b)[1])
    plot_data_b<-data.frame(plot_data_b,group)
    plot_data<-rbind(plot_data,plot_data_b)
    
    imp_i<-importance(rf)
    imp_i<-as.data.frame(imp_i)
    imp_i$variable<-rownames(imp_i)
    times2<-rep(b,times=nrow(imp_i))
    imp_i<-data.frame(imp_i,times2)
    imp<-rbind(imp,imp_i)
  } 
  names(result)<-c("times","sensitivity","specifity","auc")
  names(plot_data)<-c("sensitivity","specificity","Group")
  names(accura)<-c("Accuracy","Lower","higher","AUC","Time")
  
  write.csv(accura,file=paste0(i,"/species_alone/accura_",ratio,".csv"),row.names = T)
  write.csv(imp,file=paste0(i,"/species_alone/imp_",ratio,".csv"),row.names = T)
  
  write.csv(result,file=paste0(i,"/species_alone/result_",ratio,".csv"),row.names = T)
  write.csv(plot_data,file=paste0(i,"/species_alone/plot_data_",ratio,".csv"),row.names = T)
  write.csv(pred,file=paste0(i,"/species_alone/pred_",ratio,".csv"),row.names = T)
}



for (i in c("Hypercholesterolemia_third","Hypertriglyceridemia_third","hyperlipid_third")) {
  for (j in c(0.8)) {
    rf_classify(i,j)
  }
  
}


# 2.1.2.species+test -------------------------------------------------------------

rf_classify<-function(i,j){
  # i<-"Hypercholesterolemia_third"
  dir.create(i)
  dir.create(paste0(i,"/species_test"))
  varia<-c(i,species_name,test_name)
  m<-100
  result<-data.frame()
  plot_data<-data.frame()
  pred<-data.frame()
  set.seed(123)
  accura<-data.frame()
  result_m<-data.frame()
  imp<-data.frame()
  b=1
  ratio=j
  for (b in 1:m){
    
    Train <- createDataPartition(data_rf_second[,i], p=ratio, list=FALSE)
    training <- data_rf_second[ Train, varia]
    testing <- data_rf_second[ -Train, varia]
    
    rf = randomForest(as.formula(paste0(i,"~.") ), data=training,importance=TRUE, proximity=TRUE, ntree = 1000)
    
    pred_b=predict(rf,testing)
    df_b = table(observed = testing[,i], predict = pred_b)
    
    df_b<-as.data.frame(df_b)
    sensi_b<-df_b[4,3]/(df_b[2,3]+df_b[4,3]) 
    speci_b<-df_b[1,3]/(df_b[1,3]+df_b[3,3])
    
    pred_rf_b<-predict(rf,testing,type = "prob")
    pred_2b_b<-prediction(pred_rf_b[,2],testing[,i])
    perf_b<-performance(pred_2b_b,"tpr","fpr")
    auc_b<-performance(pred_2b_b,"auc")@y.values
    auc_b<-as.numeric(auc_b)
    
    predictions2 <- predict(rf, testing, type = "class")
    result_matrix<-confusionMatrix(predictions2, testing[,i])
    accura_b<-result_matrix$overall
    
    accura_b<-c(accura_b[1],accura_b[3],accura_b[4],auc_b,b)
    names(accura_b)<-c("Accuracy","Lower","Upper","AUC","Times")
    accura<-rbind(accura,accura_b)
    
    sensi<-result_matrix$byClass
    sensi<-as.data.frame(sensi)
    sensi<-as.data.frame(t(sensi))
    times<-rep(b,times=dim(sensi)[1])
    result_b<-data.frame(sensi,times)
    result_m<-rbind(result_m,result_b)
    
    pred_2b<-data.frame(pred_b,pred_rf_b)
    times<-rep(b,times=dim(pred_2b)[1])
    
    groundtruth<-testing[,i]
    pred_2b<-data.frame(pred_2b,groundtruth,times)
    pred<-rbind(pred,pred_2b)
    # rownames(pred)<-rownames(testing)
    #ROC曲线
    
    result_b<-c(b,sensi_b,speci_b,auc_b)
    result_b<-as.data.frame(result_b) 
    result_b<-as.data.frame(t(result_b)) 
    result<-rbind(result,result_b) 
    
    sensi_value_b<-perf_b@y.values[[1]]
    speci_value_b<-perf_b@x.values[[1]]
    plot_data_b<-data.frame(sensi_value_b,speci_value_b)
    group<-rep(b,times=dim(plot_data_b)[1])
    plot_data_b<-data.frame(plot_data_b,group)
    plot_data<-rbind(plot_data,plot_data_b)
    
    imp_i<-importance(rf)
    imp_i<-as.data.frame(imp_i)
    imp_i$variable<-rownames(imp_i)
    times2<-rep(b,times=nrow(imp_i))
    imp_i<-data.frame(imp_i,times2)
    imp<-rbind(imp,imp_i)
  } 
  names(result)<-c("times","sensitivity","specifity","auc")
  names(plot_data)<-c("sensitivity","specificity","Group")
  names(accura)<-c("Accuracy","Lower","higher","AUC","Time")
  
  write.csv(accura,file=paste0(i,"/species_test/accura_",ratio,".csv"),row.names = T)
  write.csv(imp,file=paste0(i,"/species_test/imp_",ratio,".csv"),row.names = T)
  
  write.csv(result,file=paste0(i,"/species_test/result_",ratio,".csv"),row.names = T)
  write.csv(plot_data,file=paste0(i,"/species_test/plot_data_",ratio,".csv"),row.names = T)
  write.csv(pred,file=paste0(i,"/species_test/pred_",ratio,".csv"),row.names = T)
}



for (i in c("Hypercholesterolemia_third","Hypertriglyceridemia_third","hyperlipid_third")) {
  for (j in c(0.8)) {
    rf_classify(i,j)
  }
  
}

#2.1.3. ROC_95%CI ---------------------------------------------------------------

dir_root<-"E:/study2/microbiome/study_new/total/9Data Review/Result/metagen/ML"
setwd(dir_root)
#2.1.3.1. ROC curve----------------------------------------------------

plot_result<-function(i,j,k){
  dir1<-paste0(dir_root,"/",i,"/",j)
  
  result<-read.csv(paste0(dir1,"/result_",k,".csv"),header = T,sep = ",")
  plot_data<-read.csv(paste0(dir1,"/plot_data_",k,".csv"),header = T,sep = ",")
  pred<-read.csv(paste0(dir1,"/pred_",k,".csv"),header = T,sep = ",")
  
  result1<-result[order(result$auc,decreasing = T),]
  imp<-read.csv(paste0(dir1,"/imp_",k,".csv"),header = T,sep = ",")
  
  library(pROC)
  time<-result1[50,2]
  set.seed(315)
  pred1<-pred[pred$times==time,]
  proc_obj_signature<-roc(pred1$groundtruth,pred1$X1,auc = T,ci= T,print.auc=T)
  sensi.ci<-ci.se(proc_obj_signature,specificities = seq(0, 1, .01))
  
  auc<-auc(proc_obj_signature)[1]
  auc_low<-ci(proc_obj_signature,of="auc")[1]
  auc_high<-ci(proc_obj_signature,of="auc")[3]
  auc_signature<-paste(j," AUC:",round(auc,digits = 3),"(",
                       round(auc_low,digits = 3),",",round(auc_high,digits = 3),")",sep = "")
  
  data_ci<-sensi.ci[1:101,1:3]
  data_ci<-as.data.frame(data_ci)
  y1<-data_ci[,1]
  y2<-data_ci[,3]
  x<-seq(0,1,0.01)
  data_ci_signature<-data.frame(y1,y2,x)
  sensitivity<-proc_obj_signature$sensitivities
  specificity<-proc_obj_signature$specificities
  Group<-rep(j,times=length(sensitivity))
  plot_signature<-data.frame(sensitivity,specificity,Group)
  
  #pdf(file="Hypercholesterolemia/family_test/roc_ci_2.pdf",width = 5,height = 5)
  ggroc(proc_obj_signature,color="red",size=1)+theme_bw()+
    geom_segment(aes(x = 1, y = 0, xend = 0, yend = 1), 
                 colour='grey', linetype = 'dotdash') +
    geom_ribbon(data = data_ci,aes(x=x,ymin=`2.5%`,ymax=`97.5%`), fill = 'lightblue',alpha=0.5)+
    theme(plot.title = element_text(hjust = 0.5), 
          legend.justification=c(1, 0), legend.position=c(.95, .05),
          #legend.title=title, 
          legend.background = element_rect(fill=NULL, size=0.5, 
                                           linetype="solid", colour ="black"))+
    labs(x="Specificity",y="Sensitivity",title = auc_signature)
  ggsave(paste0(dir_root,"/",i,"/",j,"/",i,"_",j,"_",k,"_roc_ci_median.pdf"),width=5,height=5)
  write.csv(data_ci,file = paste0(dir_root,"/",i,"/",j,"/",i,"_",j,"_",k,"_roc_ci_median.csv"))
  
  #dev.off()
  ###top value
  library(pROC)
  time<-result1[1,2]
  set.seed(315)
  pred1<-pred[pred$times==time,]
  proc_obj_signature<-roc(pred1$groundtruth,pred1$X1,auc = T,ci= T,print.auc=T)
  sensi.ci<-ci.se(proc_obj_signature,specificities = seq(0, 1, .01))
  
  auc<-auc(proc_obj_signature)[1]
  auc_low<-ci(proc_obj_signature,of="auc")[1]
  auc_high<-ci(proc_obj_signature,of="auc")[3]
  auc_signature<-paste(j," AUC:",round(auc,digits = 3),"(",
                       round(auc_low,digits = 3),",",round(auc_high,digits = 3),")",sep = "")
  
  data_ci<-sensi.ci[1:101,1:3]
  data_ci<-as.data.frame(data_ci)
  y1<-data_ci[,1]
  y2<-data_ci[,3]
  x<-seq(0,1,0.01)
  data_ci_signature<-data.frame(y1,y2,x)
  sensitivity<-proc_obj_signature$sensitivities
  specificity<-proc_obj_signature$specificities
  Group<-rep(j,times=length(sensitivity))
  plot_signature<-data.frame(sensitivity,specificity,Group)
  
  #pdf(file="Hypercholesterolemia/family_test/roc_ci_2.pdf",width = 5,height = 5)
  ggroc(proc_obj_signature,color="red",size=1)+theme_bw()+
    geom_segment(aes(x = 1, y = 0, xend = 0, yend = 1), 
                 colour='grey', linetype = 'dotdash') +
    geom_ribbon(data = data_ci,aes(x=x,ymin=`2.5%`,ymax=`97.5%`), fill = 'lightblue',alpha=0.5)+
    theme(plot.title = element_text(hjust = 0.5), 
          legend.justification=c(1, 0), legend.position=c(.95, .05),
          #legend.title=title, 
          legend.background = element_rect(fill=NULL, size=0.5, 
                                           linetype="solid", colour ="black"))+
    labs(x="Specificity",y="Sensitivity",title = auc_signature)
  ggsave(paste0(dir_root,"/",i,"/",j,"/",i,"_",j,"_",k,"_roc_ci_top.pdf"),width=5,height=5)
  write.csv(data_ci,file = paste0(dir_root,"/",i,"/",j,"/",i,"_",j,"_",k,"_roc_ci_top.csv"))
  
  ##importance
  ggplot(imp,aes(x=variable,y=MeanDecreaseGini))+
    geom_boxplot()+
    theme_bw()+
    # theme(axis.text.x=element_text(angle=45,size=8,vjust = 1,hjust = 1))+
    coord_flip()
  ggsave(paste0(dir_root,"/",i,"/",j,"/",i,"_",j,"_",k,"_varimp.pdf"),width=10,height=10)
  
  
}

i<-"Hypercholesterolemia_third"
j<-"species_test"

outome<-c("Hypercholesterolemia_third","Hypertriglyceridemia_third","hyperlipid_third")
#var2<-c("Family_alone","Family_test","species_test","species_alone","CAG_test","CAG_alone","test_alone")
var2<-c("species_test","species_alone")

ratio<-c(0.8)
for (i in outome) {
  for (j in var2) {
    for (k in ratio) {
      plot_result(i,j,k)
    }
    
  }
  
}

#2.1.3.2. varImportance -----------------------------------------------------------

plot_imp<-function(i,j,k){
  
  dir1<-paste0(dir_root,"/",i,"/",j)
  
  #result<-read.csv(paste0(dir1,"/result_",k,".csv"),header = T,sep = ",")
  #plot_data<-read.csv(paste0(dir1,"/plot_data_",k,".csv"),header = T,sep = ",")
  #pred<-read.csv(paste0(dir1,"/pred_",k,".csv"),header = T,sep = ",")
  
  imp<-read.csv(paste0(dir1,"/imp_",k,".csv"),header = T,sep = ",")
  imp$variable<-as.character(imp$variable)
  t<-tapply(imp$MeanDecreaseGini,imp$variable,mean)##get mean value of imp
  t<-t[order(t,decreasing = T)]
  var_name<-names(t)[1:30]
  ggplot(imp,aes(x=variable,y=MeanDecreaseGini))+
    geom_boxplot()+
    theme_bw()+
    scale_x_discrete(limits=rev(var_name))+
    # theme(axis.text.x=element_text(angle=45,size=8,vjust = 1,hjust = 1))+
    coord_flip()
  ggsave(paste0(dir_root,"/",i,"/",j,"/",i,"_",j,"_",k,"_varimp_30.pdf"),width=10,height=10)
  
}

outome<-c("Hypercholesterolemia_third","Hypertriglyceridemia_third","hyperlipid_third")
#var2<-c("Family_alone","Family_test","species_test","species_alone")
ratio<-c(0.8)
for (i in outome){
  for (j in var2) {
    for (k in ratio) {
      plot_imp(i,j,k)
    }
  }
}




# 2.1.3.3.combine ROC curve & varimportance plot-----------------------------------------------------

i<-"Hypercholesterolemia_third"
#j<-"Family_test"

outome<-c("Hypercholesterolemia_third","Hypertriglyceridemia_third","hyperlipid_third")
#var2<-c("Family_alone","Family_test","species_test","species_alone","CAG_test","CAG_alone","test_alone")

imp_total<-data.frame()
plot_signature_top<-data.frame()
plot_signature_median<-data.frame()
auc_median<-data.frame()
auc_top<-data.frame()
data_ci_median<-data.frame()
data_ci_top<-data.frame()
for (i in outome) {
  for (j in var2) {
    for (k in c(0.8)) {
      dir1<-paste0(dir_root,"/",i,"/",j)
      
      result<-read.csv(paste0(dir1,"/result_",k,".csv"),header = T,sep = ",")
      plot_data<-read.csv(paste0(dir1,"/plot_data_",k,".csv"),header = T,sep = ",")
      pred<-read.csv(paste0(dir1,"/pred_",k,".csv"),header = T,sep = ",")
      #data_ci<-read.csv(paste0(dir1,"/",i,"_",j,"_",k,"_roc_ci_top.csv"),header = T,sep = ",")
      ####imp
      imp<-read.csv(paste0(dir1,"/imp_",k,".csv"),header = T,sep = ",")
      Group<-rep(i,times=nrow(imp))
      variab<-rep(j,times=nrow(imp))
      rati<-rep(k,times=nrow(imp))
      imp_i<-data.frame(imp,Group,variab,rati)
      imp_total<-rbind(imp_total,imp_i)
      
      result1<-result[order(result$auc,decreasing = T),]
      
      library(pROC)
      time<-result1[50,2]
      set.seed(315)
      pred1<-pred[pred$times==time,]
      proc_obj_signature<-roc(pred1$groundtruth,pred1$X1,auc = T,ci= T,print.auc=T)
      sensi.ci<-ci.se(proc_obj_signature,specificities = seq(0, 1, .01))
      
      auc<-auc(proc_obj_signature)[1]
      auc_low<-ci(proc_obj_signature,of="auc")[1]
      auc_high<-ci(proc_obj_signature,of="auc")[3]
      auc_signature<-paste(j," AUC:",round(auc,digits = 3),"(",
                           round(auc_low,digits = 3),",",round(auc_high,digits = 3),")",sep = "")
      
      auc_median1<-c(auc,auc_low,auc_high,i,j,k)
      auc_median1<-as.data.frame(auc_median1)
      auc_median1<-as.data.frame(t(auc_median1))
      auc_median<-rbind(auc_median,auc_median1)
      
      data_ci<-sensi.ci[1:101,1:3]
      data_ci<-as.data.frame(data_ci)
      y1<-data_ci[,1]
      y2<-data_ci[,3]
      x<-seq(0,1,0.01)
      data_ci_signature<-data.frame(y1,y2,x)
      Group<-rep(i,times=nrow(data_ci))
      variab<-rep(j,times=nrow(data_ci))
      rati<-rep(k,times=nrow(data_ci))
      data_ci<-data.frame(data_ci,Group,variab,rati)
      data_ci_median<-rbind(data_ci_median,data_ci)
      
      
      sensitivity<-proc_obj_signature$sensitivities
      specificity<-proc_obj_signature$specificities
      Group<-rep(i,times=length(sensitivity))
      variab<-rep(j,times=length(sensitivity))
      rati<-rep(k,times=length(sensitivity))
      plot_signature_i<-data.frame(sensitivity,specificity,Group,variab,rati)
      plot_signature_median<-rbind(plot_signature_median,plot_signature_i)
      #dev.off()
      ###top value
      library(pROC)
      time<-result1[1,2]
      set.seed(315)
      pred1<-pred[pred$times==time,]
      proc_obj_signature<-roc(pred1$groundtruth,pred1$X1,auc = T,ci= T,print.auc=T)
      sensi.ci<-ci.se(proc_obj_signature,specificities = seq(0, 1, .01))
      
      auc<-auc(proc_obj_signature)[1]
      auc_low<-ci(proc_obj_signature,of="auc")[1]
      auc_high<-ci(proc_obj_signature,of="auc")[3]
      auc_signature<-paste(j," AUC:",round(auc,digits = 3),"(",
                           round(auc_low,digits = 3),",",round(auc_high,digits = 3),")",sep = "")
      
      auc_high1<-c(auc,auc_low,auc_high,i,j,k)
      auc_high1<-as.data.frame(auc_high1)
      auc_high1<-as.data.frame(t(auc_high1))
      auc_top<-rbind(auc_top,auc_high1)
      
      data_ci<-sensi.ci[1:101,1:3]
      data_ci<-as.data.frame(data_ci)
      y1<-data_ci[,1]
      y2<-data_ci[,3]
      x<-seq(0,1,0.01)
      data_ci_signature<-data.frame(y1,y2,x)
      Group<-rep(i,times=nrow(data_ci))
      variab<-rep(j,times=nrow(data_ci))
      rati<-rep(k,times=nrow(data_ci))
      data_ci<-data.frame(data_ci,Group,variab,rati)
      data_ci_top<-rbind(data_ci_top,data_ci)
      
      
      sensitivity<-proc_obj_signature$sensitivities
      specificity<-proc_obj_signature$specificities
      Group<-rep(i,times=length(sensitivity))
      variab<-rep(j,times=length(sensitivity))
      rati<-rep(k,times=length(sensitivity))
      plot_signature_1i<-data.frame(sensitivity,specificity,Group,variab,rati)
      plot_signature_top<-rbind(plot_signature_top,plot_signature_1i)
      
      
    }
  }
}


write.csv(plot_signature_top,file = "plot_signature_top.csv")
write.csv(plot_signature_median,file = "plot_signature_median.csv")
write.csv(imp_total,file = "imp_total.csv")
write.csv(data_ci_median,file = "data_ci_median.csv")
write.csv(data_ci_top,file = "data_ci_top.csv")
write.csv(auc_top,file = "auc_top.csv")
write.csv(auc_median,file = "auc_median.csv")

imp_total<-read.csv("imp_total.csv",header = T,sep = ",")
plot_signature_median<-read.csv("plot_signature_median.csv",header = T,sep = ",")
plot_signature_top<-read.csv("plot_signature_top.csv",header = T,sep = ",")

# data_ci_median<-read.csv("data_ci_median.csv",header = T,sep = ",")

#pdf(file="Hypercholesterolemia/family_test/roc_ci_2.pdf",width = 5,height = 5)

require(ggplot2)
library(ggsci)
for (i in outome) {
  for (j in c(0.8)) {
    plot_data<-subset(plot_signature_top,Group==i&rati==j)
    
    plot_data<-subset(plot_signature_top,Group==i&rati==j)
    p<-ggplot(plot_data, aes(x = 1-specificity, y=sensitivity)) +
      geom_path(aes(color = variab), size=1) +
      geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
                   colour='grey', linetype = 'dotdash') +
      theme_bw() + 
      theme(plot.title = element_text(hjust = 0.5), 
            legend.justification=c(1, 0), legend.position=c(.95, .05),
            legend.title=element_blank(), 
            legend.background = element_rect(fill=NULL, size=0.5, 
                                             linetype="solid", colour ="black"))+
      labs(title = i)+
      scale_color_npg()
    ggsave(p,filename = paste0(i,"_",j,".pdf"),width = 8,height = 8)
    
  }
  
}


# 
# ggplot(imp,aes(x=variable,y=MeanDecreaseGini))+
#   geom_boxplot()+
#   theme_bw()+
#   scale_x_discrete(limits=rev(var_name))+
#   # theme(axis.text.x=element_text(angle=45,size=8,vjust = 1,hjust = 1))+
#   coord_flip()
# ggsave(paste0(paste0(dir_root,"/",i,"/",j,"/",i,"_",j,"_varimp.pdf")),width=10,height=10)
# 





#2.1.3.4. auc --------------------------------------------------------------------

dir_root<-"E:/study2/microbiome/study_new/total/9Data Review/Result/metagen/ML"
setwd(dir_root)
outome<-c("Hypercholesterolemia_third","Hypertriglyceridemia_third","hyperlipid_third")
#var2<-c("Family_alone","Family_test","species_test","species_alone","CAG_test","CAG_alone","test_alone")
ratio<-c(0.8)

i<-outome[1]
j<-var2[1]


data_full<-data.frame()
data_full_total<-data.frame()
for (i in outome) {
  for (j in var2) {
    dir1<-paste0(dir_root,"/",i,"/",j)
    files2<-list.files(dir1,pattern = "accura*.")
    setwd(dir1)
    # data_full<-data.frame()
    for (k in files2){
      dat<-read.csv(k,header = T,sep = ",")
      group<-rep(j,times=dim(dat)[1])
      group2<-rep(i,times=dim(dat)[1])
      group3<-rep(k,times=dim(dat)[1])
      dat<-data.frame(dat,group,group2,group3)
      data_full<-rbind(data_full,dat)
    }
    # category<-rep(i,times=nrow(data_full))
    # data_full2<-data.frame(data_full,category)
    data_full_total<-rbind(data_full_total,data_full)
    #write.csv(data_full,file = "data_full.csv",row.names = F)
  }
  
}

data_full_2<-subset(data_full_total,group3=="accura_0.8.csv")


setwd("E:/study2/microbiome/study_new/total/9Data Review/Result/metagen/ML/combine")
for (i in var2) {
  plot_dat<-subset(data_full_2,group==i)
  p1<-ggboxplot(plot_dat,x="group2",y="AUC",#add="iitter",
                color = "group2")+
    # stat_compare_means(method = "kruskal.test")+
    theme_classic()+
    theme(legend.position = "none")+
    theme(panel.border = element_rect(color = "black",fill = NA,size = 1))+
    
    scale_color_npg()+
    labs(x="",y="AUC")
  #  facet_wrap(~group_2,nrow = 2)
  p2<-ggboxplot(plot_dat,x="group2",y="Accuracy",#add="iitter",
                color = "group2")+
    # stat_compare_means(method = "kruskal.test")+
    theme_classic()+
    theme(legend.position = "none")+
    theme(panel.border = element_rect(color = "black",fill = NA,size = 1))+
    
    scale_color_npg()+
    # scale_fill_manual(values = cbPalette)+
    
    labs(x="",y="Accuracy")
  ggsave(p1,filename =paste(i,"_AUC.pdf"),width =6,height = 4 )
  ggsave(p2,filename =paste(i,"_Accuracy.pdf"),width = 6,height = 4 )
  
}




#2.2 multi_cate --------------------------------------------------------------

Dyslipidemia<-data_combine$dyslipimedia_32

data_rf_second<-data.frame(data_combine,Dyslipidemia)
for (i in test_name) {
  
  data_rf_second[,i][is.na(data_rf_second[,i])==TRUE]<-mean(data_rf_second[,i],na.rm = T)
}

dir_root<-"E:/study2/microbiome/study_new/total/9Data Review/Result/metagen/ML"
setwd(dir_root)
library(openxlsx)
library(randomForest)
library(ggplot2)
library(nnet)
library(car)

library(DALEX)
library(iBreakDown)
library(nnet)
library(questionr)
library(dplyr)
#2.2.1.species_alone -------------------------------------------------------
data_rf_second$Dyslipidemia<-as.factor(data_rf_second$Dyslipidemia)

library(randomForest)
library(bootstrap)
library(ROCR)
library(pROC)
library(dplyr)
library(caret)
library(multiROC)

rf_multiclassify<-function(i,j){
  m<-100
  result<-data.frame()
  plot_data<-data.frame()
  pred<-data.frame()
  accura<-data.frame()
  set.seed(123)
  
  auc_cate<-data.frame()
  #variable
  i<-"Dyslipidemia"
  ratio=j
  dir.create(i)
  dir.create(paste0(i,"/species_alone"))
  #varia<-c(i,index)
  varia<-c(i,species_name)
  b=1
  imp<-data.frame()
  data<-data_rf_second
  for (b in 1:m){
    Train <- createDataPartition(data[,i], p=ratio, list=FALSE)
    training <- data[ Train, varia]
    testing <- data[ -Train, varia]
    
    rf=randomForest(as.formula(paste(i,"~.",sep = "")),data=training,importance=TRUE, proximity=TRUE, ntree = 500)
    predictions2 <- predict(rf, testing, type = "class")
    result_matrix<-confusionMatrix(predictions2, testing[,i])
    accura_b<-result_matrix$overall
    
    sensi<-result_matrix$byClass
    sensi<-as.data.frame(sensi)
    times<-rep(b,times=dim(sensi)[1])
    result_b<-data.frame(sensi,times)
    result<-rbind(result,result_b)
    
    groundtruth<-testing[,i]
    times<-rep(b,times=length(groundtruth))
    pred_2b<-data.frame(predictions2,groundtruth,times)
    pred<-rbind(pred,pred_2b)
    # rownames(pred)<-rownames(testing)
    
    #AUC
    true_label<-dummies::dummy(testing[,i],sep=".")
    true_label <- data.frame(true_label)
    colnames(true_label) <- gsub(".*?\\.", "", colnames(true_label))
    colnames(true_label) <- paste(colnames(true_label), "_true")
    
    pred_value <- predict(rf, testing, type = "prob")
    label<-colnames(pred_value)
    colnames(pred_value)<-paste(label,"_pred_value")
    #合并数据集
    final_value<-cbind(true_label,pred_value)
    #multiROC and multiPR
    roc_res_value <- multi_roc(final_value, force_diag=T)
    pr_res_value <- multi_pr(final_value, force_diag=T)
    
    # auc_inflamed<-roc_res_value$AUC$value$`1 `
    # auc_Exclude<-roc_res_value$AUC$value$`2`
    # auc_Desert<-roc_res_value$AUC$value$`3`
    # auc_average<-roc_res_value$AUC$value$macro
    auc_value<-roc_res_value$AUC$value
    auc_name<-names(unlist(auc_value))
    auc_num<-unlist(auc_value)
    auc_cate_i<-c(auc_num, b)
    auc_cate_i<-as.data.frame(auc_cate_i)
    auc_cate_i<-as.data.frame(t(auc_cate_i))
    auc_cate<-rbind(auc_cate,auc_cate_i)
    
    #Plot
    plot_roc_value <- plot_roc_data(roc_res_value)
    plot_pr_value <- plot_pr_data(pr_res_value)
    
    times2<-rep(b,times=nrow(plot_roc_value))
    plot_roc_value2<-data.frame(plot_roc_value,times2)
    plot_data<-rbind(plot_data,plot_roc_value2)
    
    plot_roc_value<-plot_roc_value[plot_roc_value$Group=="Macro",]
    auc_b<-plot_roc_value[1,4]
    auc_b<-round(auc_b,digits = 3)
    
    accura_b<-c(accura_b[1],accura_b[3],accura_b[4],auc_b,b)
    names(accura_b)<-c("Accuracy","Lower","Upper","AUC","Times")
    accura<-rbind(accura,accura_b)
    
    imp_i<-importance(rf)
    imp_i<-as.data.frame(imp_i)
    imp_i$variable<-rownames(imp_i)
    times2<-rep(b,times=nrow(imp_i))
    imp_i<-data.frame(imp_i,times2)
    imp<-rbind(imp,imp_i)
    
  } 
  names(accura)<-c("Accuracy","Lower","higher","AUC","Time")
  #names(auc_cate)<-c("1","2","3","4","5","macro","micro","times")
  
  names(accura)<-c("Accuracy","Lower","higher","AUC","Time")
  #names(result)<-c("times","sensitivity","specifity","auc")
  #names(plot_data)<-c("sensitivity","specificity","Group")
  
  write.csv(accura,file=paste0(i,"/species_alone/accura_",ratio,"_full.csv"),row.names = T)
  write.csv(result,file=paste0(i,"/species_alone/result_",ratio,"_full.csv"),row.names = T)
  write.csv(auc_cate,file=paste0(i,"/species_alone/auc_cate",ratio,"_full.csv"),row.names = T)
  
  #write.csv(result_m,file=paste0(i,"/species_alone/result_m_",ratio,"_full.csv"),row.names = T)
  #write.csv(plot_data,file="Result2/plot_data_bootstrap.csv",row.names = T)
  write.csv(pred,file=paste0(i,"/species_alone/pred_",ratio,"_full.csv"),row.names = T)
  write.csv(plot_data,file=paste0(i,"/species_alone/plot_data_",ratio,"_full.csv"),row.names = T)
  write.csv(imp,file=paste0(i,"/species_alone/imp",ratio,"_full.csv"),row.names = T)
  
}


for (i in "Dyslipidemia") {
  for (j in c(0.8)) {
    rf_multiclassify(i,j)
  }
  
}

#2.2.2.species_test -------------------------------------------------------

library(randomForest)
library(bootstrap)
library(ROCR)
library(pROC)
library(dplyr)
library(caret)
library(multiROC)

rf_multiclassify<-function(i,j){
  m<-100
  result<-data.frame()
  plot_data<-data.frame()
  pred<-data.frame()
  accura<-data.frame()
  set.seed(123)
  
  auc_cate<-data.frame()
  #variable
  i<-"Dyslipidemia"
  ratio=j
  dir.create(i)
  dir.create(paste0(i,"/species_test"))
  #varia<-c(i,index)
  varia<-c(i,species_name,test_name)
  b=1
  imp<-data.frame()
  data<-data_rf_second
  for (b in 1:m){
    Train <- createDataPartition(data[,i], p=ratio, list=FALSE)
    training <- data[ Train, varia]
    testing <- data[ -Train, varia]
    
    rf=randomForest(as.formula(paste(i,"~.",sep = "")),data=training,importance=TRUE, proximity=TRUE, ntree = 500)
    predictions2 <- predict(rf, testing, type = "class")
    result_matrix<-confusionMatrix(predictions2, testing[,i])
    accura_b<-result_matrix$overall
    
    sensi<-result_matrix$byClass
    sensi<-as.data.frame(sensi)
    times<-rep(b,times=dim(sensi)[1])
    result_b<-data.frame(sensi,times)
    result<-rbind(result,result_b)
    
    groundtruth<-testing[,i]
    times<-rep(b,times=length(groundtruth))
    pred_2b<-data.frame(predictions2,groundtruth,times)
    pred<-rbind(pred,pred_2b)
    # rownames(pred)<-rownames(testing)
    
    #AUC
    true_label<-dummies::dummy(testing[,i],sep=".")
    true_label <- data.frame(true_label)
    colnames(true_label) <- gsub(".*?\\.", "", colnames(true_label))
    colnames(true_label) <- paste(colnames(true_label), "_true")
    
    pred_value <- predict(rf, testing, type = "prob")
    label<-colnames(pred_value)
    colnames(pred_value)<-paste(label,"_pred_value")
    #合并数据集
    final_value<-cbind(true_label,pred_value)
    #multiROC and multiPR
    roc_res_value <- multi_roc(final_value, force_diag=T)
    pr_res_value <- multi_pr(final_value, force_diag=T)
    
    # auc_inflamed<-roc_res_value$AUC$value$`1 `
    # auc_Exclude<-roc_res_value$AUC$value$`2`
    # auc_Desert<-roc_res_value$AUC$value$`3`
    # auc_average<-roc_res_value$AUC$value$macro
    auc_value<-roc_res_value$AUC$value
    auc_name<-names(unlist(auc_value))
    auc_num<-unlist(auc_value)
    auc_cate_i<-c(auc_num, b)
    auc_cate_i<-as.data.frame(auc_cate_i)
    auc_cate_i<-as.data.frame(t(auc_cate_i))
    auc_cate<-rbind(auc_cate,auc_cate_i)
    
    #Plot
    plot_roc_value <- plot_roc_data(roc_res_value)
    plot_pr_value <- plot_pr_data(pr_res_value)
    
    times2<-rep(b,times=nrow(plot_roc_value))
    plot_roc_value2<-data.frame(plot_roc_value,times2)
    plot_data<-rbind(plot_data,plot_roc_value2)
    
    plot_roc_value<-plot_roc_value[plot_roc_value$Group=="Macro",]
    auc_b<-plot_roc_value[1,4]
    auc_b<-round(auc_b,digits = 3)
    
    accura_b<-c(accura_b[1],accura_b[3],accura_b[4],auc_b,b)
    names(accura_b)<-c("Accuracy","Lower","Upper","AUC","Times")
    accura<-rbind(accura,accura_b)
    
    imp_i<-importance(rf)
    imp_i<-as.data.frame(imp_i)
    imp_i$variable<-rownames(imp_i)
    times2<-rep(b,times=nrow(imp_i))
    imp_i<-data.frame(imp_i,times2)
    imp<-rbind(imp,imp_i)
    
  } 
  names(accura)<-c("Accuracy","Lower","higher","AUC","Time")
  #names(auc_cate)<-c("1","2","3","4","5","macro","micro","times")
  
  names(accura)<-c("Accuracy","Lower","higher","AUC","Time")
  #names(result)<-c("times","sensitivity","specifity","auc")
  #names(plot_data)<-c("sensitivity","specificity","Group")
  
  write.csv(accura,file=paste0(i,"/species_test/accura_",ratio,"_full.csv"),row.names = T)
  write.csv(result,file=paste0(i,"/species_test/result_",ratio,"_full.csv"),row.names = T)
  write.csv(auc_cate,file=paste0(i,"/species_test/auc_cate",ratio,"_full.csv"),row.names = T)
  
  #write.csv(result_m,file=paste0(i,"/species_test/result_m_",ratio,"_full.csv"),row.names = T)
  #write.csv(plot_data,file="Result2/plot_data_bootstrap.csv",row.names = T)
  write.csv(pred,file=paste0(i,"/species_test/pred_",ratio,"_full.csv"),row.names = T)
  write.csv(plot_data,file=paste0(i,"/species_test/plot_data_",ratio,"_full.csv"),row.names = T)
  write.csv(imp,file=paste0(i,"/species_test/imp",ratio,"_full.csv"),row.names = T)
  
}


for (i in "Dyslipidemia") {
  for (j in c(0.8)) {
    rf_multiclassify(i,j)
  }
  
}

#2.3. combine -----------------------------------------------------------------
#2.3.1. auc --------------------------------------------------------------------
library(ggplot2)
library(ggpubr)
library(ggsci)
var2<-c("species_test","species_alone")
data_full_total<-data.frame()

for (i in var2) {
  dir_root<-paste("E:/study2/microbiome/study_new/total/9Data Review/Result/metagen/ML/Dyslipidemia/",i,sep = "")
  files2<-list.files(dir_root,pattern = "accura*.")
  setwd(dir_root)
  data_full<-data.frame()
  for (j in files2){
    dat<-read.csv(j,header = T,sep = ",")
    group<-rep(j,times=dim(dat)[1])
    dat<-data.frame(dat,group)
    data_full<-rbind(data_full,dat)
  }
  category<-rep(i,times=nrow(data_full))
  data_full2<-data.frame(data_full,category)
  data_full_total<-rbind(data_full_total,data_full2)
  #write.csv(data_full,file = "data_full.csv",row.names = F)
}


data_full<-splitstackshape::cSplit(data_full_total, splitCols="group", sep="_")
data_full$group_3<-gsub("[_*.csv]", "", data_full$group_3)

data_full$group_2<-as.factor(data_full$group_2)

p<-ggboxplot(data_full,x="category",y="AUC",#add="iitter",
             fill = "category",notch = TRUE)+
  stat_compare_means(method = "kruskal.test")+
  theme(legend.position = "none")+
  scale_fill_npg()+
  # scale_fill_manual(values = cbPalette)+
  theme_bw()+
  labs(x="",y="AUC")+
  facet_wrap(~group_2,nrow = 2)
setwd("E:/study2/microbiome/study_new/total/9Data Review/Result/metagen/ML/Dyslipidemia")
# p<-ggplot() +
#   geom_boxplot( data = data_full,aes(x=group_2,y=AUC,fill=factor(group_2)),outlier.colour = NA) +
#   theme_bw()+
#  # scale_x_discrete(limits=c("full","index","glm"))+
#   labs(x="")
ggsave(p,filename = "boxplot_auc.pdf",width = 10,height = 8) 

p<-ggboxplot(data_full,x="category",y="Accuracy",#add="iitter",
             fill = "category",notch = TRUE)+
  stat_compare_means(method = "kruskal.test")+
  theme(legend.position = "none")+
  scale_fill_npg()+
  # scale_fill_manual(values = cbPalette)+
  theme_bw()+
  labs(x="",y="Accuracy")+
  facet_wrap(~group_2,nrow = 2)
setwd("E:/study2/microbiome/study_new/total/9Data Review/Result/metagen/ML/Dyslipidemia")
# p<-ggplot() +
#   geom_boxplot( data = data_full,aes(x=group_2,y=AUC,fill=factor(group_2)),outlier.colour = NA) +
#   theme_bw()+
#  # scale_x_discrete(limits=c("full","index","glm"))+
#   labs(x="")
ggsave(p,filename = "boxplot_Accuracy.pdf",width = 10,height = 8) 


auc<-with(data_full,do.call(rbind, tapply(AUC, category, function(x) c(Q = quantile(x,probs = c(0.5,0.025,0.975))))))
accuracy<-with(data_full, do.call(rbind, tapply(Accuracy, category, function(x) c(Q = quantile(x,probs = c(0.5,0.025,0.975))))))

# auc_0.6<-with(subset(data_full,group_2==0.6), do.call(rbind, tapply(AUC, category, function(x) c(Q = quantile(x,probs = c(0.5,0.025,0.975))))))
# auc_0.8<-with(subset(data_full,group_2==0.8), do.call(rbind, tapply(AUC, category, function(x) c(Q = quantile(x,probs = c(0.5,0.025,0.975))))))
# 
# accuracy_0.6<-with(subset(data_full,group_2==0.6), do.call(rbind, tapply(Accuracy, category, function(x) c(Q = quantile(x,probs = c(0.5,0.025,0.975))))))
# accuracy_0.8<-with(subset(data_full,group_2==0.8), do.call(rbind, tapply(Accuracy, category, function(x) c(Q = quantile(x,probs = c(0.5,0.025,0.975))))))

result_combine1<-data.frame(auc,accuracy)
names(result_combine1)<-c("auc","auc_l","auc_u",
                          "accuracy","accuracy_l","accuracy_u")
write.csv(result_combine1,file = "result_combine_auc.csv")



#2.3.2. auc_cate ----------------------------------------------------------------

library(ggplot2)
library(ggpubr)
library(ggsci)
#var2<-c("Family_alone","Family_test","species_test","species_alone","CAG_test","CAG_alone","test_alone")
data_full_total<-data.frame()

for (i in var2) {
  dir_root<-paste("E:/study2/microbiome/study_new/total/9Data Review/Result/metagen/ML/Dyslipidemia/",i,sep = "")
  files2<-list.files(dir_root,pattern = "auc_cate*.")
  setwd(dir_root)
  data_full<-data.frame()
  for (j in files2){
    dat<-read.csv(j,header = T,sep = ",")
    group<-rep(j,times=dim(dat)[1])
    dat<-data.frame(dat,group)
    data_full<-rbind(data_full,dat)
  }
  category<-rep(i,times=nrow(data_full))
  data_full2<-data.frame(data_full,category)
  data_full_total<-rbind(data_full_total,data_full2)
  #write.csv(data_full,file = "data_full.csv",row.names = F)
}
setwd("E:/study2/microbiome/study_new/total/9Data Review/Result/metagen/ML/Dyslipidemia")

write.csv(data_full_total,file = "data_full_auc_cate.csv",row.names = F)

data_full<-splitstackshape::cSplit(data_full_total, splitCols="group", sep="_")
data_full$group_3<-gsub("[_*.csv]", "", data_full$group_3)

dir.create("Auc_figure")
library(psych)
library(reshape2)
dir.create("Auc_figure")
library(psych)
library(reshape2)

def_plot<-function(i){
  plot_data1<-subset(data_full,category==i)
  plot_data1<-melt(plot_data1[,c(2:5,7)])
  p<-ggplot(data=plot_data1,aes(x=variable,y=value,color=variable))+
    geom_boxplot(alpha=1, outlier.size=0, size=0.9, width=0.6)+
    labs(x="", y="AUC")+theme_bw()+
    theme(panel.grid.major=element_blank(),legend.position = "none",panel.grid.minor=element_blank())+
    #  scale_colour_manual(values=cbPalette)
    scale_color_npg()
  ggsave(p,filename = paste("AUC_figure//",i,"_boxplot.pdf",sep = ""),width = 10,height = 5) 
}

for (i in var2) {
  def_plot(i)
  
}

d2<-melt(data_full[,c(2:7,9)])
auc_combine<-data.frame()
for (i in var2) {
  auc_i<-with(subset(d2,category==i), do.call(rbind, tapply(value, variable, function(x) c(Q = quantile(x,probs = c(0.5,0.025,0.975))))))
  auc_i<-as.data.frame(auc_i)
  auc_i$cate<-rownames(auc_i)
  group<-rep(i,times=nrow(auc_i))
  auc_i<-data.frame(auc_i,group)
  auc_combine<-rbind(auc_combine,auc_i)
}

write.csv(auc_combine,file = "auc_seperate.csv")



#2.3.3. confusion matrix --------------------------------------------------------

# 2.3.3.1 according to all type ---------------------------------------------------
group_version<-c("Control","Hypercholesterolemia","Hyperlipidemia","Hypertriglyceridemia")
#var2<-c("Family_alone","Family_test","species_test","species_alone","CAG_test","CAG_alone","test_alone")

library(reshape2)
for (i in var2) {
  dir_root<-paste("E:/study2/microbiome/study_new/total/9Data Review/Result/metagen/ML/Dyslipidemia/",i,sep = "")
  files2<-list.files(dir_root,pattern = "result_0.*")
  setwd(dir_root)
  data_full<-data.frame()
  for (t in files2){
    dat<-read.csv(t,header = T,sep = ",")
    group<-rep(t,times=dim(dat)[1])
    dat<-data.frame(dat,group)
    data_full<-rbind(data_full,dat)
  }
  
  write.csv(data_full,file = "result_matrix_full.csv",row.names = F)
  
  data_full<-splitstackshape::cSplit(data_full, splitCols="group", sep="_")
  data_full$group_3<-gsub("[_*.csv]", "", data_full$group_3)
  
  category<-rep(group_version,times=dim(data_full)[1]/4)
  data_full<-data.frame(category,data_full)
  
  
  variabs<-c("Sensitivity","Specificity","Precision","Recall","F1")
  
  sensi<-with(subset(data_full), do.call(rbind, tapply(Sensitivity, category, function(x) c(Q = quantile(x,probs = c(0.5,0.025,0.975))))))
  speci<-with(subset(data_full), do.call(rbind, tapply(Specificity, category, function(x) c(Q = quantile(x,probs = c(0.5,0.025,0.975))))))
  preci<-with(subset(data_full), do.call(rbind, tapply(Precision, category, function(x) c(Q = quantile(x,probs = c(0.5,0.025,0.975),na.rm = T)))))
  reca<-with(subset(data_full), do.call(rbind, tapply(Recall, category, function(x) c(Q = quantile(x,probs = c(0.5,0.025,0.975))))))
  f1<-with(subset(data_full), do.call(rbind, tapply(F1, category, function(x) c(Q = quantile(x,probs = c(0.5,0.025,0.975),na.rm = T)))))
  
  result_combine1<-data.frame(sensi,speci,preci,reca,f1)
  names(result_combine1)<-c("sensi_m","sensi_l","sensi_h","speci_m","speci_l","speci_h","preci_m","preci_l","preci_h",
                            "recall_m","recall_l","recall_h","F1_m","F1_l","F1_h")
  write.csv(result_combine1,file = "result_sensis.csv")
  
  
  def_plot<-function(s,u){
    variab<-c("category",s)
    plot_data<-subset(data_full,group_2==u)
    plot_data1<-plot_data[variab]
    p<-ggplot() +
      geom_boxplot( data = plot_data1,aes(x=category,y=plot_data1[,s],
                                          color=factor(category)),
                    outlier.colour = NA) +
      theme_bw()+
      # scale_x_discrete(limits=group_version)+
      theme(axis.text.x=element_text(angle=0,size=8,vjust = 1,hjust = 1))+
      labs(x="",y=s)+
      scale_color_npg()
    ggsave(p,filename = paste(s,"_",u,"_boxplot.pdf",sep = ""),width = 10,height = 6) 
  }
  
  for (j in variabs){
    for (k in c(0.8)) {
      def_plot(j,k)
    }
    
  }
}



#2.3.4. roc curve ---------------------------------------------------------------


# full --------------------------------------------------------------------


for (i in var2) {
  dir_root<-paste("E:/study2/microbiome/study_new/total/9Data Review/Result/metagen/ML/Dyslipidemia/",i,sep = "")
  
  setwd(dir_root)
  result<-read.csv("accura_0.8_full.csv",header = T)
  plot_data<-read.csv("plot_data_0.8_full.csv",header = T)
  plot_data$Group<-as.character(plot_data$Group)
  #plot_data<-subset(plot_data,Group!="Micro")
  # plot_data$Group2[plot_data$Group==1]<-"Inflamed"
  # plot_data$Group2[plot_data$Group==2]<-"Exclude"
  # plot_data$Group2[plot_data$Group==3]<-"Desert"
  # plot_data$Group2[plot_data$Group=="Macro"]<-"Average"
  
  
  pred<-read.csv("pred_0.8_full.csv",header = T)
  
  result1<-result[order(result$AUC,decreasing = F),]
  
  library(pROC)
  time<-result1[51,6]
  set.seed(315)
  pred1<-pred[pred$times==time,]
  colnames(plot_data)<-c("X","Sensitivity","Specificity","Group",
                         "cate","value","times2")
  plot_data2<-subset(plot_data,times2==time)
  
  #pdf(file="ROC_0.8_full.pdf")
  p<-ggplot(plot_data2, aes(x = 1-Specificity, y=Sensitivity)) +
    geom_path(aes(color = Group), size=1) +
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
                 colour='grey', linetype = 'dotdash') +
    theme_bw() + 
    labs(main="")+
    theme(plot.title = element_text(hjust = 0.5), 
          legend.justification=c(1, 0), legend.position=c(.95, .05),
          legend.title=element_blank(), 
          legend.background = element_rect(fill=NULL, size=0.5, 
                                           linetype="solid", colour ="black"))+
    # scale_colour_manual(values=cbbPalette)
    scale_color_npg()
  #dev.off()
  ggsave(file="ROC_0.8_full.pdf",width = 8,height = 8)
}


# imp ---------------------------------------------------------------------
dir_root<-"E:/study2/microbiome/study_new/total/9Data Review/Result/metagen/ML"

plot_imp<-function(i,j,k){
  
  dir1<-paste0(dir_root,"/",i,"/",j)
  
  #result<-read.csv(paste0(dir1,"/result_",k,".csv"),header = T,sep = ",")
  #plot_data<-read.csv(paste0(dir1,"/plot_data_",k,".csv"),header = T,sep = ",")
  #pred<-read.csv(paste0(dir1,"/pred_",k,".csv"),header = T,sep = ",")
  
  imp<-read.csv(paste0(dir1,"/imp",k,"_full.csv"),header = T,sep = ",")
  imp$variable<-as.character(imp$variable)
  t<-tapply(imp$MeanDecreaseGini,imp$variable,mean)##get mean value of imp
  t<-t[order(t,decreasing = T)]
  var_name<-names(t)[1:30]
  ggplot(imp,aes(x=variable,y=MeanDecreaseGini))+
    geom_boxplot()+
    theme_bw()+
    scale_x_discrete(limits=rev(var_name))+
    # theme(axis.text.x=element_text(angle=45,size=8,vjust = 1,hjust = 1))+
    coord_flip()
  ggsave(paste0(dir_root,"/",i,"/",j,"/",i,"_",j,"_",k,"_varimp_30.pdf"),width=10,height=10)
  
}

outome<-c("Dyslipidemia")
var2<-c("species_test","species_alone")
ratio<-c(0.8)

for (i in outome){
  for (j in var2) {
    for (k in ratio) {
      plot_imp(i,j,k)
    }
  }
}



dir_root<-"E:/study2/microbiome/study_new/total/9Data Review/Result/metagen/ML"
setwd(dir_root)



plot_result<-function(i,j,k){
  dir1<-paste0(dir_root,"/",i,"/",j)
  
  result<-read.csv(paste0(dir1,"/result_",k,".csv"),header = T,sep = ",")
  plot_data<-read.csv(paste0(dir1,"/plot_data_",k,".csv"),header = T,sep = ",")
  pred<-read.csv(paste0(dir1,"/pred_",k,".csv"),header = T,sep = ",")
  
  result1<-result[order(result$auc,decreasing = T),]
  imp<-read.csv(paste0(dir1,"/imp_",k,".csv"),header = T,sep = ",")
  
  library(pROC)
  time<-result1[50,2]
  set.seed(315)
  pred1<-pred[pred$times==time,]
  proc_obj_signature<-roc(pred1$groundtruth,pred1$X1,auc = T,ci= T,print.auc=T)
  sensi.ci<-ci.se(proc_obj_signature,specificities = seq(0, 1, .01))
  
  auc<-auc(proc_obj_signature)[1]
  auc_low<-ci(proc_obj_signature,of="auc")[1]
  auc_high<-ci(proc_obj_signature,of="auc")[3]
  auc_signature<-paste(j," AUC:",round(auc,digits = 3),"(",
                       round(auc_low,digits = 3),",",round(auc_high,digits = 3),")",sep = "")
  
  data_ci<-sensi.ci[1:101,1:3]
  data_ci<-as.data.frame(data_ci)
  y1<-data_ci[,1]
  y2<-data_ci[,3]
  x<-seq(0,1,0.01)
  data_ci_signature<-data.frame(y1,y2,x)
  sensitivity<-proc_obj_signature$sensitivities
  specificity<-proc_obj_signature$specificities
  Group<-rep(j,times=length(sensitivity))
  plot_signature<-data.frame(sensitivity,specificity,Group)
  
  #pdf(file="Hypercholesterolemia/family_test/roc_ci_2.pdf",width = 5,height = 5)
  ggroc(proc_obj_signature,color="red",size=1)+theme_bw()+
    geom_segment(aes(x = 1, y = 0, xend = 0, yend = 1), 
                 colour='grey', linetype = 'dotdash') +
    geom_ribbon(data = data_ci,aes(x=x,ymin=`2.5%`,ymax=`97.5%`), fill = 'lightblue',alpha=0.5)+
    theme(plot.title = element_text(hjust = 0.5), 
          legend.justification=c(1, 0), legend.position=c(.95, .05),
          #legend.title=title, 
          legend.background = element_rect(fill=NULL, size=0.5, 
                                           linetype="solid", colour ="black"))+
    labs(x="Specificity",y="Sensitivity",title = auc_signature)
  ggsave(paste0(dir_root,"/",i,"/",j,"/",i,"_",j,"_",k,"_roc_ci_median.pdf"),width=5,height=5)
  write.csv(data_ci,file = paste0(dir_root,"/",i,"/",j,"/",i,"_",j,"_",k,"_roc_ci_median.csv"))
  
  #dev.off()
  ###top value
  library(pROC)
  time<-result1[1,2]
  set.seed(315)
  pred1<-pred[pred$times==time,]
  proc_obj_signature<-roc(pred1$groundtruth,pred1$X1,auc = T,ci= T,print.auc=T)
  sensi.ci<-ci.se(proc_obj_signature,specificities = seq(0, 1, .01))
  
  auc<-auc(proc_obj_signature)[1]
  auc_low<-ci(proc_obj_signature,of="auc")[1]
  auc_high<-ci(proc_obj_signature,of="auc")[3]
  auc_signature<-paste(j," AUC:",round(auc,digits = 3),"(",
                       round(auc_low,digits = 3),",",round(auc_high,digits = 3),")",sep = "")
  
  data_ci<-sensi.ci[1:101,1:3]
  data_ci<-as.data.frame(data_ci)
  y1<-data_ci[,1]
  y2<-data_ci[,3]
  x<-seq(0,1,0.01)
  data_ci_signature<-data.frame(y1,y2,x)
  sensitivity<-proc_obj_signature$sensitivities
  specificity<-proc_obj_signature$specificities
  Group<-rep(j,times=length(sensitivity))
  plot_signature<-data.frame(sensitivity,specificity,Group)
  
  #pdf(file="Hypercholesterolemia/family_test/roc_ci_2.pdf",width = 5,height = 5)
  ggroc(proc_obj_signature,color="red",size=1)+theme_bw()+
    geom_segment(aes(x = 1, y = 0, xend = 0, yend = 1), 
                 colour='grey', linetype = 'dotdash') +
    geom_ribbon(data = data_ci,aes(x=x,ymin=`2.5%`,ymax=`97.5%`), fill = 'lightblue',alpha=0.5)+
    theme(plot.title = element_text(hjust = 0.5), 
          legend.justification=c(1, 0), legend.position=c(.95, .05),
          #legend.title=title, 
          legend.background = element_rect(fill=NULL, size=0.5, 
                                           linetype="solid", colour ="black"))+
    labs(x="Specificity",y="Sensitivity",title = auc_signature)
  ggsave(paste0(dir_root,"/",i,"/",j,"/",i,"_",j,"_",k,"_roc_ci_top.pdf"),width=5,height=5)
  write.csv(data_ci,file = paste0(dir_root,"/",i,"/",j,"/",i,"_",j,"_",k,"_roc_ci_top.csv"))
  
  ##importance
  ggplot(imp,aes(x=variable,y=MeanDecreaseGini))+
    geom_boxplot()+
    theme_bw()+
    # theme(axis.text.x=element_text(angle=45,size=8,vjust = 1,hjust = 1))+
    coord_flip()
  ggsave(paste0(dir_root,"/",i,"/",j,"/",i,"_",j,"_",k,"_varimp.pdf"),width=10,height=10)
  
  
}

i<-"Hypercholesterolemia_third"
#j<-"Family_test"

outome<-c("Hypercholesterolemia_third","Hypertriglyceridemia_third","hyperlipid_third")
#var2<-c("Family_alone","Family_test","Genus_test","Genus_alone","CAG_test","CAG_alone","test_alone")

var2<-c("species_test","species_alone")
ratio<-c(0.8)
for (i in outome) {
  for (j in var2) {
    for (k in ratio) {
      plot_result(i,j,k)
    }
    
  }
  
}


# 2.4、ROC combine -----------------------------------------------------------




# top ---------------------------------------------------------------------


i<-"Hypertriglyceridemia_third"
j<-"species_test"
k<-0.8

dir1<-paste0(dir_root,"/",i,"/",j)

result<-read.csv(paste0(dir1,"/result_",k,".csv"),header = T,sep = ",")
plot_data<-read.csv(paste0(dir1,"/plot_data_",k,".csv"),header = T,sep = ",")
pred<-read.csv(paste0(dir1,"/pred_",k,".csv"),header = T,sep = ",")

result1<-result[order(result$auc,decreasing = T),]
imp<-read.csv(paste0(dir1,"/imp_",k,".csv"),header = T,sep = ",")

library(pROC)
time<-result1[1,2]
set.seed(315)
pred1<-pred[pred$times==time,]
Hypertriglyceridemia<-roc(pred1$groundtruth,pred1$X1,auc = T,ci= T,print.auc=T)



i<-"Hypercholesterolemia_third"
j<-"species_test"
k<-0.8

dir1<-paste0(dir_root,"/",i,"/",j)

result<-read.csv(paste0(dir1,"/result_",k,".csv"),header = T,sep = ",")
plot_data<-read.csv(paste0(dir1,"/plot_data_",k,".csv"),header = T,sep = ",")
pred<-read.csv(paste0(dir1,"/pred_",k,".csv"),header = T,sep = ",")

result1<-result[order(result$auc,decreasing = T),]
imp<-read.csv(paste0(dir1,"/imp_",k,".csv"),header = T,sep = ",")

library(pROC)
time<-result1[1,2]
set.seed(315)
pred1<-pred[pred$times==time,]
Hypercholesterolemia<-roc(pred1$groundtruth,pred1$X1,auc = T,ci= T,print.auc=T)





i<-"hyperlipid_third"
j<-"species_test"
k<-0.8

dir1<-paste0(dir_root,"/",i,"/",j)

result<-read.csv(paste0(dir1,"/result_",k,".csv"),header = T,sep = ",")
plot_data<-read.csv(paste0(dir1,"/plot_data_",k,".csv"),header = T,sep = ",")
pred<-read.csv(paste0(dir1,"/pred_",k,".csv"),header = T,sep = ",")

result1<-result[order(result$auc,decreasing = T),]
imp<-read.csv(paste0(dir1,"/imp_",k,".csv"),header = T,sep = ",")

library(pROC)
library(plotROC)
time<-result1[1,2]
set.seed(315)
pred1<-pred[pred$times==time,]
hyperlipid<-roc(pred1$groundtruth,pred1$X1,auc = T,ci= T,print.auc=T)


rocLs<-c(Hypertriglyceridemia,Hypercholesterolemia,hyperlipid)

rocLs = list('Hypertriglyceridemia'=Hypertriglyceridemia, 
             'Hypercholesterolemia' = Hypercholesterolemia,"Hyperlipid"=hyperlipid)

drawROC<-function(rocLs,conf.level=0.95){
  rocDat=purrr::map2(rocLs,names(rocLs),function(x,GrpName){
    dat=data.frame(
      'predictor'=x$predictor,
      'response'=x$response,
      'Grp'=GrpName
    )
    return(dat)
  })
  rocDat=rlist::list.rbind(rocDat)
  rocciDat=purrr::map2(rocLs,names(rocLs),function(x,GrpName){
    dat=ci.se(x,x$specificities,conf.level=conf.level,progress="none")
    sps=rownames(dat)%>%as.numeric()
    dat=data.frame(dat)
    dat$Grp=GrpName
    dat$Specificity=sps
    return(dat)
  })
  rocciDat=rlist::list.rbind(rocciDat)
  rocciDat$Grp<-ordered(rocciDat$Grp,levels=c("Hypertriglyceridemia","Hypercholesterolemia",
                                              "Hyperlipid"))
  p=ggplot(data=rocDat)+
    geom_roc(aes(m=predictor,d=response,color=Grp),n.cuts=0)+
    geom_ribbon(data=rocciDat,aes(x=1-Specificity,ymin=X2.5.,ymax=X97.5.,fill=Grp),alpha=0.1)+
    theme_bw()+
    # geom_segment(aes(x = 1, y = 0, xend = 0, yend = 1), 
    #              colour='grey', linetype = 'dotdash') +
    
    theme(panel.border = element_rect(color = "black",fill = NA,size = 1))+
    guides(color=guide_legend(title = ''),fill=guide_legend(title = ''))+
    scale_color_npg()+
    scale_fill_npg()
  return(p)
}

Rocpics=drawROC(rocLs)
Rocpics

setwd(dir_root)
dir.create("combine")
ggsave(Rocpics,filename = "combine/roc_combine_species_top.pdf",width = 10,height = 8)

# median ---------------------------------------------------------------------


i<-"Hypertriglyceridemia_third"
j<-"species_test"
k<-0.8

dir1<-paste0(dir_root,"/",i,"/",j)

result<-read.csv(paste0(dir1,"/result_",k,".csv"),header = T,sep = ",")
plot_data<-read.csv(paste0(dir1,"/plot_data_",k,".csv"),header = T,sep = ",")
pred<-read.csv(paste0(dir1,"/pred_",k,".csv"),header = T,sep = ",")

result1<-result[order(result$auc,decreasing = T),]
imp<-read.csv(paste0(dir1,"/imp_",k,".csv"),header = T,sep = ",")

library(pROC)
time<-result1[51,2]
set.seed(315)
pred1<-pred[pred$times==time,]
Hypertriglyceridemia<-roc(pred1$groundtruth,pred1$X1,auc = T,ci= T,print.auc=T)



i<-"Hypercholesterolemia_third"
j<-"species_test"
k<-0.8

dir1<-paste0(dir_root,"/",i,"/",j)

result<-read.csv(paste0(dir1,"/result_",k,".csv"),header = T,sep = ",")
plot_data<-read.csv(paste0(dir1,"/plot_data_",k,".csv"),header = T,sep = ",")
pred<-read.csv(paste0(dir1,"/pred_",k,".csv"),header = T,sep = ",")

result1<-result[order(result$auc,decreasing = T),]
imp<-read.csv(paste0(dir1,"/imp_",k,".csv"),header = T,sep = ",")

library(pROC)
time<-result1[51,2]
set.seed(315)
pred1<-pred[pred$times==time,]
Hypercholesterolemia<-roc(pred1$groundtruth,pred1$X1,auc = T,ci= T,print.auc=T)





i<-"hyperlipid_third"
j<-"species_test"
k<-0.8

dir1<-paste0(dir_root,"/",i,"/",j)

result<-read.csv(paste0(dir1,"/result_",k,".csv"),header = T,sep = ",")
plot_data<-read.csv(paste0(dir1,"/plot_data_",k,".csv"),header = T,sep = ",")
pred<-read.csv(paste0(dir1,"/pred_",k,".csv"),header = T,sep = ",")

result1<-result[order(result$auc,decreasing = T),]
imp<-read.csv(paste0(dir1,"/imp_",k,".csv"),header = T,sep = ",")

library(pROC)
library(plotROC)
time<-result1[51,2]
set.seed(315)
pred1<-pred[pred$times==time,]
hyperlipid<-roc(pred1$groundtruth,pred1$X1,auc = T,ci= T,print.auc=T)


rocLs<-c(Hypertriglyceridemia,Hypercholesterolemia,hyperlipid)

rocLs = list('Hypertriglyceridemia'=Hypertriglyceridemia, 
             'Hypercholesterolemia' = Hypercholesterolemia,"Hyperlipid"=hyperlipid)

drawROC<-function(rocLs,conf.level=0.95){
  rocDat=purrr::map2(rocLs,names(rocLs),function(x,GrpName){
    dat=data.frame(
      'predictor'=x$predictor,
      'response'=x$response,
      'Grp'=GrpName
    )
    return(dat)
  })
  rocDat=rlist::list.rbind(rocDat)
  rocciDat=purrr::map2(rocLs,names(rocLs),function(x,GrpName){
    dat=ci.se(x,x$specificities,conf.level=conf.level,progress="none")
    sps=rownames(dat)%>%as.numeric()
    dat=data.frame(dat)
    dat$Grp=GrpName
    dat$Specificity=sps
    return(dat)
  })
  rocciDat=rlist::list.rbind(rocciDat)
  rocciDat$Grp<-ordered(rocciDat$Grp,levels=c("Hypertriglyceridemia","Hypercholesterolemia",
                                              "Hyperlipid"))
  p=ggplot(data=rocDat)+
    geom_roc(aes(m=predictor,d=response,color=Grp),n.cuts=0)+
    geom_ribbon(data=rocciDat,aes(x=1-Specificity,ymin=X2.5.,ymax=X97.5.,fill=Grp),alpha=0.1)+
    theme_bw()+
    # geom_segment(aes(x = 1, y = 0, xend = 0, yend = 1), 
    #              colour='grey', linetype = 'dotdash') +
    
    theme(panel.border = element_rect(color = "black",fill = NA,size = 1))+
    guides(color=guide_legend(title = ''),fill=guide_legend(title = ''))+
    scale_color_npg()+
    scale_fill_npg()
  return(p)
}

Rocpics=drawROC(rocLs)
Rocpics

ggsave(Rocpics,filename = "combine/roc_combine_species_median.pdf",width = 10,height = 8)





# plot_reg ----------------------------------------------------------------

# plot  -------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(ggsci)
library(openxlsx)
setwd("E:/study2/microbiome/study_new/total/9Data Review")
data_plot<-read.xlsx("Result/metagen/linear_all.xlsx",sheet = 1)
#data_plot<-splitstackshape::cSplit(data_plot, splitCols="species", sep="_")

data_plot_new<-subset(data_plot,pvalue<0.05)
species_selected<-as.character(unique(data_plot_new$species))
data_plot2<-subset(data_plot,species%in%species_selected)


data_plot_unadjusted<-subset(data_plot2,group=="Unadjusted")
data_plot_adjusted<-subset(data_plot2,group=="Adjusted")

beta2<-data_plot_unadjusted$beta
data_plot2<-data.frame(data_plot_adjusted,beta2)

data_plot_1<-subset(data_plot2,outcome=="TG_24")

data_plot_1<-arrange(data_plot_1, desc(beta))

name<-as.character(unique(data_plot_1$species))

ggplot(data_plot_1,aes(x=species,y=beta))+
  geom_errorbar(aes(ymin=lower,ymax=upper),position=position_dodge(width=0.8),width=0.1)+
  geom_point(size=3.5,position=position_dodge(width=0.8),color="#00afbb")+
  geom_point(aes(x=species,y=beta2),size=3.5,shape=17,color="#e7b800")+
  # geom_point(size=3.5,shape=18,position=position_dodge(width=0.8))+
  geom_hline(yintercept=0,color="red",linetype="dashed")+
  labs(x="",y="β (95%CI)")+
  theme_bw()+
  scale_color_npg()+
  theme(axis.text.y = element_text(face = "italic"))+
  scale_x_discrete(limits=name)+
  # facet_wrap(~variable1,nrow = 1,scales = "free")+
  coord_flip()
ggsave(filename = "Result/metagen/plot_LME_linear_TG_24.pdf",width = 8,height = 8)

data_plot_1<-subset(data_plot2,outcome=="TG_32")

data_plot_1<-arrange(data_plot_1, desc(beta))

name<-as.character(unique(data_plot_1$species))

ggplot(data_plot_1,aes(x=species,y=beta))+
  geom_errorbar(aes(ymin=lower,ymax=upper),position=position_dodge(width=0.8),width=0.1)+
  geom_point(size=3.5,position=position_dodge(width=0.8),color="#00afbb")+
  geom_point(aes(x=species,y=beta2),size=3.5,shape=17,color="#e7b800")+
  # geom_point(size=3.5,shape=18,position=position_dodge(width=0.8))+
  geom_hline(yintercept=0,color="red",linetype="dashed")+
  labs(x="",y="β (95%CI)")+
  theme_bw()+
  scale_color_npg()+
  theme(axis.text.y = element_text(face = "italic"))+
  scale_x_discrete(limits=name)+
  # facet_wrap(~variable1,nrow = 1,scales = "free")+
  coord_flip()
ggsave(filename = "Result/metagen/plot_LME_linear_TG_32.pdf",width = 8,height = 8)



data_plot_1<-subset(data_plot2,outcome=="CHOL_24")

data_plot_1<-arrange(data_plot_1, desc(beta))

name<-as.character(unique(data_plot_1$species))

ggplot(data_plot_1,aes(x=species,y=beta))+
  geom_errorbar(aes(ymin=lower,ymax=upper),position=position_dodge(width=0.8),width=0.1)+
  geom_point(size=3.5,position=position_dodge(width=0.8),color="#00afbb")+
  geom_point(aes(x=species,y=beta2),size=3.5,shape=17,color="#e7b800")+
  # geom_point(size=3.5,shape=18,position=position_dodge(width=0.8))+
  geom_hline(yintercept=0,color="red",linetype="dashed")+
  labs(x="",y="β (95%CI)")+
  theme_bw()+
  scale_color_npg()+
  scale_x_discrete(limits=name)+
  theme(axis.text.y = element_text(face = "italic"))+
  # facet_wrap(~variable1,nrow = 1,scales = "free")+
  coord_flip()
ggsave(filename = "Result/metagen/plot_LME_linear_CHOL_24.pdf",width = 8,height = 8)

data_plot_1<-subset(data_plot2,outcome=="CHOL_32")

data_plot_1<-arrange(data_plot_1, desc(beta))

name<-as.character(unique(data_plot_1$species))

ggplot(data_plot_1,aes(x=species,y=beta))+
  geom_errorbar(aes(ymin=lower,ymax=upper),position=position_dodge(width=0.8),width=0.1)+
  geom_point(size=3.5,position=position_dodge(width=0.8),color="#00afbb")+
  geom_point(aes(x=species,y=beta2),size=3.5,shape=17,color="#e7b800")+
  # geom_point(size=3.5,shape=18,position=position_dodge(width=0.8))+
  geom_hline(yintercept=0,color="red",linetype="dashed")+
  labs(x="",y="β (95%CI)")+
  theme_bw()+
  scale_color_npg()+
  scale_x_discrete(limits=name)+
  theme(axis.text.y = element_text(face = "italic"))+
  # facet_wrap(~variable1,nrow = 1,scales = "free")+
  coord_flip()
ggsave(filename = "Result/metagen/plot_LME_linear_CHOL_32.pdf",width = 8,height = 8)



# 二、KO --------------------------------------------------------------------
setwd('E:/study2/microbiome/study_new/total/9Data Review')
library(openxlsx)
metadata = read.xlsx("Data/data_total.xlsx", sheet=1,rowNames = F)

metadata2 = read.table("Data/metagen/ID.txt", header = T,sep = "\t")
#metadata2<-subset(metadata2,dyslipimedia_24=="Control")

#ID<-as.character(metadata2$number)
metadata_new<-inner_join(metadata2,metadata,by="number")


# 一、KO ---------------------------------------------------------------

dat_KO<-read.table("Data/metagen/ko_relab_unstratified.tsv",
                   header = T,row.names = 1,sep = "\t")
KO_name<-rownames(dat_KO)
dat_KO_full<-as.data.frame(t(dat_KO))
dat_KO_full$SampleID<-rownames(dat_KO_full)
#rownames(metadata2)<-metadata2$SampleID

dat_KO_full<-merge(metadata_new,dat_KO_full,by= "SampleID")
write.csv(dat_KO_full,file = "Data/metagen/dat_KO_full.csv",row.names = T)

# 
# transform_RA<-function(x) x*100 / sum(x)
# dat_KO<-dat_KO
# dat_KO2<-apply(dat_KO, 2, transform_RA)


abundance = 1e-04

idx = rowMeans(dat_KO) > abundance
filtered_data = dat_KO[idx,]

# detectionRate=10
# 
# detect<-function(x){
#   length(which(x > 0))*100/length(x)
# }
# 
# DR<-apply(filtered_data, 1, detect)
# idx3= apply(filtered_data, 1, detect) >detectionRate
# 
# filtered_data2<-filtered_data[idx3,]
KO_name<-rownames(filtered_data)

write.csv(filtered_data,file = "Data/metagen/KO_filtered.csv",row.names = T)

filtered_data<-as.data.frame(t(filtered_data))
filtered_data$ID<-rownames(filtered_data)
#rownames(metadata2)<-metadata2$SampleID
dat_KO2<-merge(metadata2,filtered_data,by.x = "SampleID",by.y = "ID")
dat_KO2<-dat_KO2[,-c(1,3,4)]
dat_KO2<-merge(metadata_new,dat_KO2,by="number")
write.csv(dat_KO2,file = "Data/metagen/dat_KO.csv",row.names = T)



#1. linear regression -------------------------------------------------------

library(fdrtool)
x_var<-c("TG_24","CHOL_24","TG_32","CHOL_32")
y_var<-paste0("path",seq(1,length(KO_name),by=1))

cov<-c("pre_BMI","Parity","Age","gesweek_24")
dir.create("Result/metagen")
dir.create("Result/metagen/association")

# linear ------------------------------------------------------------------
data_combine<-dat_KO2
data_combine[KO_name]<-log(data_combine[KO_name]+1e-05)
#rename column name
colnames(data_combine)[colnames(data_combine) %in% KO_name] <- y_var
# unadjusted ----------------------------------------------------------------


###unadjusted
result<-data.frame()
for (i in x_var){
  for(j in y_var){
    # xnam_i1<-paste(i,j,sep="*")
    # xnam_i<-paste(j,"pre_BMI","Parity","Age",sep="+")
    fit_i<-lm(as.formula(paste(i,"~",j)),data=data_combine)
    p_i<-tryCatch({summary(fit_i)$coefficients[2,4]},error=function(e){"NA"})
    b_i<-tryCatch({coef(fit_i)[2]},error=function(e){"NA"})
    blow_i<-tryCatch({confint(fit_i)[2,1]},error=function(e){"NA"})
    bhigh_i<-tryCatch({confint(fit_i)[2,2]},error=function(e){"NA"})
    result_i<-c(i,j,b_i,blow_i,bhigh_i,p_i)  
    result_i<-as.data.frame(result_i) 
    result_i<-as.data.frame(t(result_i)) 
    result<-rbind(result,result_i) 
  }
}
result<-data.frame(KO_name,result)
names(result)<-c("KO_name","outcome","KO","beta","lower","upper","pvalue")   

result_new<-data.frame()
for (i in x_var) {
  result2<-subset(result,outcome==i)
  p<-as.numeric(as.character(result2$pvalue))
  fdr=fdrtool(p,statistic="pvalue",plot = F)
  qval<-data.frame(fdr$qval)
  p2<-p.adjust(p,method = "fdr",length(p))
  p2<-data.frame(p2)
  result2<-cbind(result2,qval,p2)
  result_new<-rbind(result_new,result2)
  
}

group<-rep("Unadjusted",times=nrow(result_new))

result1<-data.frame(result_new,group)



# adjusted ----------------------------------------------------------------



###adjusted
result<-data.frame()
for (i in x_var){
  for(j in y_var){
    # xnam_i1<-paste(i,j,sep="*")
    xnam_i<-paste(j,"pre_BMI","Parity","Age","gesweek_24",sep="+")
    fit_i<-lm(as.formula(paste(i,"~",xnam_i)),data=data_combine)
    p_i<-tryCatch({summary(fit_i)$coefficients[2,4]},error=function(e){"NA"})
    b_i<-tryCatch({coef(fit_i)[2]},error=function(e){"NA"})
    blow_i<-tryCatch({confint(fit_i)[2,1]},error=function(e){"NA"})
    bhigh_i<-tryCatch({confint(fit_i)[2,2]},error=function(e){"NA"})
    result_i<-c(i,j,b_i,blow_i,bhigh_i,p_i)  
    result_i<-as.data.frame(result_i) 
    result_i<-as.data.frame(t(result_i)) 
    result<-rbind(result,result_i) 
  }
}


result<-data.frame(KO_name,result)
names(result)<-c("KO_name","outcome","KO","beta","lower","upper","pvalue")   

result_new<-data.frame()
for (i in x_var) {
  result2<-subset(result,outcome==i)
  p<-as.numeric(as.character(result2$pvalue))
  fdr=fdrtool(p,statistic="pvalue",plot = F)
  qval<-data.frame(fdr$qval)
  p2<-p.adjust(p,method = "fdr",length(p))
  p2<-data.frame(p2)
  result2<-cbind(result2,qval,p2)
  result_new<-rbind(result_new,result2)
  
}


group<-rep("Adjusted",times=nrow(result_new))

result2<-data.frame(result_new,group)

result_full<-rbind(result1,result2)
write.csv(result_full,file="Result/metagen/association/linear_all_KO.csv",row.names = F) 


# binary ------------------------------------------------------------------


# uunadjusted -------------------------------------------------------------


###binary
x_var<-c("Hyperlipid_24","TG_24_High","CHOL_24_High",
         "Hyperlipid_32","TG_32_High","CHOL_32_High")


result<-data.frame()
for (i in x_var){
  for(j in y_var){
    # xnam_i1<-paste(i,j,sep="*")
    #xnam_i<-paste(j,"pre_BMI","Parity","Age",sep="+")
    fit_i<-glm(as.formula(paste(i,"~",j)),family = "binomial",data=data_combine)
    p_i<-tryCatch({summary(fit_i)$coefficients[2,4]},error=function(e){"NA"})
    b_i<-tryCatch({exp(coef(fit_i)[2])},error=function(e){"NA"})
    blow_i<-tryCatch({exp(confint(fit_i)[2,1])},error=function(e){"NA"})
    bhigh_i<-tryCatch({exp(confint(fit_i)[2,2])},error=function(e){"NA"})
    result_i<-c(i,j,b_i,blow_i,bhigh_i,p_i)
    result_i<-as.data.frame(result_i)
    result_i<-as.data.frame(t(result_i))
    result<-rbind(result,result_i)
  }
}
result<-data.frame(KO_name,result)

names(result)<-c("KO_name","outcome","KO","OR","low","upper","pvalueֵ")

result_new<-data.frame()
for (i in x_var) {
  result2<-subset(result,outcome==i)
  p<-as.numeric(as.character(result2$pvalue))
  fdr=fdrtool(p,statistic="pvalue",plot = F)
  qval<-data.frame(fdr$qval)
  p2<-p.adjust(p,method = "fdr",length(p))
  p2<-data.frame(p2)
  result2<-cbind(result2,qval,p2)
  result_new<-rbind(result_new,result2)
  
}
group<-rep("Unadjusted",times=nrow(result_new))

result1<-data.frame(result_new,group)


# adjusted ----------------------------------------------------------------


result<-data.frame()
for (i in x_var){
  for(j in y_var){
    # xnam_i1<-paste(i,j,sep="*")
    xnam_i<-paste(j,"pre_BMI","Parity","Age","gesweek_24",sep="+")
    fit_i<-glm(as.formula(paste(i,"~",xnam_i)),family = "binomial",data=data_combine)
    p_i<-tryCatch({summary(fit_i)$coefficients[2,4]},error=function(e){"NA"})
    b_i<-tryCatch({exp(coef(fit_i)[2])},error=function(e){"NA"})
    blow_i<-tryCatch({exp(confint(fit_i)[2,1])},error=function(e){"NA"})
    bhigh_i<-tryCatch({exp(confint(fit_i)[2,2])},error=function(e){"NA"})
    result_i<-c(i,j,b_i,blow_i,bhigh_i,p_i)
    result_i<-as.data.frame(result_i)
    result_i<-as.data.frame(t(result_i))
    result<-rbind(result,result_i)
  }
}
result<-data.frame(KO_name,result)

names(result)<-c("KO_name","outcome","KO","OR","low","upper","pvalueֵ")


result_new<-data.frame()
for (i in x_var) {
  result2<-subset(result,outcome==i)
  p<-as.numeric(as.character(result2$pvalue))
  fdr=fdrtool(p,statistic="pvalue",plot = F)
  qval<-data.frame(fdr$qval)
  p2<-p.adjust(p,method = "fdr",length(p))
  p2<-data.frame(p2)
  result2<-cbind(result2,qval,p2)
  result_new<-rbind(result_new,result2)
  
}


group<-rep("Adjusted",times=nrow(result_new))

result2<-data.frame(result_new,group)

result_full<-rbind(result1,result2)

write.csv(result_full,file="Result/metagen/association/linear_all_KO_group.csv",row.names = F)



# plot --------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(ggsci)
library(openxlsx)
library(dplyr)
setwd("E:/study2/microbiome/study_new/total/9Data Review")
data_plot<-read.csv("Result/metagen/association/linear_all_KO.csv", header = T,sep = ",")
data_plot<-subset(data_plot,group=="Adjusted")
#data_plot<-splitstackshape::cSplit(data_plot, splitCols="species", sep="_")
ko_description<-read.xlsx("Data/metagen/KEGG/ko_description.xlsx",sheet = 1)
#ko_description<-splitstackshape::cSplit(ko_description, splitCols="Description", sep=";")

data_plot_new<-inner_join(data_plot,ko_description,by=c("KO_name"="KO"))
write.csv(data_plot_new,file = "Result/metagen/ko_description_combine.csv",row.names = F)
data_plot_new1<-subset(data_plot_new,pvalue<0.05)
KO_selected<-as.character(unique(data_plot_new1$KO))
data_plot2<-subset(data_plot_new,KO%in%KO_selected)

data_plot2$name<-paste(data_plot2$KO_name,data_plot2$Description1,sep = ":")
data_plot2$logp<- -log10(data_plot2$pvalue)
# dotplot <- 
#   ggplot(data_plot2,aes(x=outcome, y = name, color = beta, size = logp)) + 
#   
#   geom_point()+
#   geom_point(shape=1,color="black") + 
#   cowplot::theme_cowplot() + 
#   theme(axis.line  = element_blank()) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#   xlab('')+
#   ylab('') +
#   theme(axis.ticks = element_blank()) +
#   scale_color_gradient2(low="blue", high = "red", name = 'regression coefficient')+
#   #  scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,10), oob = scales::squish, name = '-log (P value)') +
#   scale_y_discrete(position = "right")


bl <- colorRampPalette(c("navy","royalblue","lightskyblue"))(200)                      
re <- colorRampPalette(c("mistyrose", "red2","darkred"))(200)
data_plot2$outcome<-ordered(data_plot2$outcome,levels=c("TG_24","CHOL_24","TG_32","CHOL_32"))
data_plot2$label<-ifelse(data_plot2$pvalue<0.05,round(data_plot2$pvalue,3),"")


getSig <- function(dc) {
  sc <- ''
  if (dc < 0.001) sc <- '***'
  else if (dc < 0.01) sc <- '**'
  else if (dc < 0.05) sc <- '*'
  sc
}

data_plot2$label2<-sapply(data_plot2$pvalue,getSig)
library("scales")
p<-ggplot(data_plot2, aes(outcome, KO_name)) +
  geom_tile(aes(fill = beta)) +
#  scale_fill_gradient2(low="blue", high = "red", name = 'regression coefficient')+
  scale_fill_gradientn(colours=c("#4DBBD5","white", "#E64B35"), 
                       values = rescale(c(-.8,0,.5)),
                       guide = "colorbar", limits=c(-.8,.5),
                       na.value = "grey98")+
  geom_text(aes(label=label2))+
  labs(x="",y="")
ggsave(p,filename = "Result/metagen/plot_KO.pdf",width = 5,height = 8)


p2<-ggplot(data_plot2, aes(name,outcome)) +
  geom_tile(aes(fill = beta)) +
  #  scale_fill_gradient2(low="blue", high = "red", name = 'regression coefficient')+
  scale_fill_gradientn(colours=c("#4DBBD5","white", "#E64B35"), 
                       values = rescale(c(-.8,0,.5)),
                       guide = "colorbar", limits=c(-.8,.5),
                       na.value = "grey98")+
  geom_text(aes(label=label2))+
  labs(x="",y="")
ggsave(p2,filename = "Result/metagen/plot_KO_name.pdf",width = 8,height = 8)
#    limits = c(-1, 1))+

# beta2<-data_plot_unadjusted$beta
# data_plot2<-data.frame(data_plot_adjusted,beta2)
# 
# data_plot_1<-subset(data_plot2,outcome=="TG_24")
# 
# data_plot_1<-arrange(data_plot_1, desc(beta))
# 
# name<-as.character(unique(data_plot_1$KO))
# 
# ggplot(data_plot_1,aes(x=KO,y=beta))+
#   geom_errorbar(aes(ymin=lower,ymax=upper),position=position_dodge(width=0.8),width=0.1)+
#   geom_point(size=3.5,position=position_dodge(width=0.8),color="#00afbb")+
#   geom_point(aes(x=KO,y=beta2),size=3.5,shape=17,color="#e7b800")+
#   # geom_point(size=3.5,shape=18,position=position_dodge(width=0.8))+
#   geom_hline(yintercept=0,color="red",linetype="dashed")+
#   labs(x="",y="β (95%CI)")+
#   theme_bw()+
#   scale_color_npg()+
#   theme(axis.text.y = element_text(face = "italic"))+
#   scale_x_discrete(limits=name)+
#   # facet_wrap(~variable1,nrow = 1,scales = "free")+
#   coord_flip()
# ggsave(filename = "Result/metagen/plot_LME_linear_TG_24.pdf",width = 8,height = 8)
# 
# 
# plot-heatmap2 --------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(ggsci)
library(openxlsx)
library(dplyr)
setwd("E:/study2/microbiome/study_new/total/9Data Review")
data_plot<-read.csv("Result/metagen/association/linear_all_KO.csv", header = T,sep = ",")
data_plot<-subset(data_plot,group=="Adjusted")
#data_plot<-splitstackshape::cSplit(data_plot, splitCols="species", sep="_")
#ko_description1<-read.xlsx("Data/metagen/KEGG/KO1-4.xlsx",sheet = 1)
ko_description<-read.xlsx("Data/metagen/KEGG/ko_description.xlsx",sheet = 1)
#ko_description<-splitstackshape::cSplit(ko_description, splitCols="Description", sep=";")

data_plot_new<-inner_join(data_plot,ko_description,by=c("KO_name"="KO"))
write.csv(data_plot_new,file = "Result/metagen/ko_description_combine.csv",row.names = F)
data_plot_new1<-subset(data_plot_new,pvalue<0.05)
KO_selected<-as.character(unique(data_plot_new1$KO))
data_plot2<-subset(data_plot_new,KO%in%KO_selected)

#data_plot2$name<-paste(data_plot2$KO_name,data_plot2$Description1,sep = ":")
data_plot2$logp<- -log10(data_plot2$pvalue)

dat_heat_beta<-dcast(data_plot2[,c("KO_name","outcome","beta")],KO_name~outcome,value.var = "beta")
rownames(dat_heat_beta)<-dat_heat_beta$KO_name
dat_heat_beta<-dat_heat_beta[,-1]
dat_heat_beta<-as.matrix(dat_heat_beta)

dat_heat_p<-dcast(data_plot2[,c("KO_name","outcome","pvalue")],KO_name~outcome,value.var = "pvalue")
rownames(dat_heat_p)<-dat_heat_p$KO_name
dat_heat_p<-dat_heat_p[,-1]
dat_heat_p<-as.matrix(dat_heat_p)

row_anno<-subset(data_plot2,outcome=="TG_24")[,c("KO_name","PathwayL1","PathwayL2","Description1")]
rownames(row_anno)<-row_anno$KO_name
row_anno<-row_anno[,c("PathwayL1","PathwayL2")]
#row_anno<-as.data.frame(row_anno)
# ann_colors = list(
#   # Season=c(Summer = "#8491B4", Winter = "#91D1C2"),
#   PathwayL1 = c(`Brite Hierarchies` = "#E64B35", 
#                 `Genetic Information Processing` = "#4DBBD5", 
#                 Metabolism="#91D1C2",
#                 `Not Included in Pathway or Brite`="F39B7F",
#                 `Organismal Systems`="#8491B4"),
#   PathwayL2=c(`HNSC_HPV-` ="#F39B7F",`HNSC_HPV+` ="#3C5488"))

getSig <- function(dc) {
  sc <- ''
  if (dc < 0.001) sc <- '***'
  else if (dc < 0.01) sc <- '**'
  else if (dc < 0.05) sc <- '*'
  sc
}
sig.mat <- matrix(sapply(dat_heat_p, getSig), nrow=nrow(dat_heat_p))
str(sig.mat)
pheatmap::pheatmap(dat_heat_beta, scale="none",
                   cluster_rows=T,
                 #  breaks = breaksList,
                 #  annotation_col=row_anno,
                   color=colorRampPalette(c("#4DBBD5","white","#E64B35"))(6),
                 #  annotation_colors = ann_colors[1],
                   display_numbers = sig.mat,border=FALSE,
                   angle_col="270",
                   filename="Result/metagen/pheatmap_anno_filter.pdf",width = 18,height = 8)

library(pheatmap)
p<-pheatmap(dat_heat_beta, scale="none",cluster_row = T,cluster_cols = T,
         display_numbers = sig.mat,
         #  annotation_col = 
         color=colorRampPalette(c("#4DBBD5","white","#E64B35"))(6),
         border=FALSE)
ggsave(p,filename = "Result/metagen/plot_KO_heatmap.pdf",width = 4,height = 8)



hmcol <- function(n){
  colfun1 <- colorRampPalette(rev(brewer.pal(9, "Blues")))
  colfun2 <- colorRampPalette(brewer.pal(9, "Reds"))
  p <- (0 - min(mat)) / (max(mat) - min(mat))
  gap <- n * p
  c(colfun1(floor(gap)), colfun2(n - gap))
}
mat<-dat_heat_beta
p<-pheatmap(t(mat), scale="none",cluster_row = T,cluster_cols = T,
            display_numbers = t(sig.mat),
            fontface="italic",
            #  annotation_col = 
            #  color=colorRampPalette(c("#4DBBD5","white","#E64B35"))(6),
            color=hmcol(50),
            border=FALSE)
ggsave(p,filename = "Result/metagen/plot_KO_heatmap2.pdf",width = 8,height = 3)

# KO_species --------------------------------------------------------------
library(data.table)
library(dplyr)
dir.create("Result/metagen/abundance")
dat<-read.table("Data/metagen/ko_relab_stratified.tsv",
                header = T,sep = "\t")


dat_1<-splitstackshape::cSplit(dat, splitCols="Gene.Family", sep="|")
sampleID<-as.character(dat_KO_full$SampleID)
vari<-c(sampleID,"Gene.Family_1","Gene.Family_2")
dat_2<-dat_1[,..vari]

ko_remain<-data_plot2$KO_name

dat_combine<-data.frame()
for (i in ko_remain) {
 # i<-ko_remain[1]
  
  dat_plot_i<-subset(dat_2,Gene.Family_1==i)
  dat_plot_i<-as.data.frame(dat_plot_i)
  
  `%notin%` <- Negate(`%in%`)
  dat_plot_1i<-subset(dat_plot_i,Gene.Family_2%notin%"unclassified")
  
  #t<-apply(dat_plot_i[,1:(ncol(dat_plot_i)-2)], 1, sum)
  dat_plot_1i$sum_<-rowSums(dat_plot_1i[,1:(ncol(dat_plot_1i)-2)])
#  dat_plot_1i<-dat_plot_1i[order(dat_plot_1i$sum_,decreasing = T),]
  
  dat_plot_1i<-arrange(dat_plot_1i, desc(sum_))
  
  top10<-as.character(dat_plot_1i$Gene.Family_2)[1:10]
  getSig <- function(dc) {
    sc <- ''
    if (dc %in% top10) sc <- dc
    else 
      sc <- 'Others'
    # else if (dc < 0.05) sc <- '*'
    sc
  }
  dat_plot_1i$Gene.Family_2<-as.character(dat_plot_1i$Gene.Family_2)
  dat_plot_1i$species<-sapply(dat_plot_1i$Gene.Family_2,getSig)
  
  Data_sum = aggregate(dat_plot_1i[,sampleID], by=dat_plot_1i["species"], FUN=sum)
  
  Data_sum<-splitstackshape::cSplit(Data_sum, splitCols="species", sep=".")
  Data_sum<-splitstackshape::cSplit(Data_sum, splitCols="species_1", sep="__")
  Data_sum<-splitstackshape::cSplit(Data_sum, splitCols="species_2", sep="__")
  
  species_new<-paste(Data_sum$species_1_2,Data_sum$species_2_2,sep = " ")  
  species_new<-c(species_new[-11],"Others")
  
  Data_sum<-as.data.frame(Data_sum)
  data_sum_1<-Data_sum[sampleID]
  rownames(data_sum_1)<-species_new
  data_mean_2<-rowMeans(data_sum_1)
  data_sum_2<-rowSums(data_sum_1)
  dat_com<-data.frame(data_mean_2,data_sum_2)
  dat_com$species<-rownames(dat_com)
  dat_com$KO<-rep(i,times=nrow(dat_com))
  dat_combine<-rbind(dat_combine,dat_com)
  
  dat_com<-arrange(dat_com, desc(data_mean_2))
  level_new<-dat_com$species
  dat_com$species<-ordered(dat_com$species,levels=level_new)
  
  p<-ggplot(dat_com,aes(x=KO,y=data_mean_2,fill=species))+geom_bar(stat="identity")+
    #scale_fill_npg()+
    labs(y="Relative abundance (%)")+
    #geom_text(aes(label=level_new))+
    labs(x="")
  ggsave(p,filename = paste0("Result/metagen/abundance/",i,"_mean.pdf"),width = 6,height = 8)
  
  p1<-ggplot(dat_com,aes(x=KO,y=data_sum_2,fill=species))+geom_bar(stat="identity")+
    #scale_fill_npg()+
    labs(y="Relative abundance (%)")+
    #geom_text(aes(label=level_new))+
    labs(x="")
  ggsave(p1,filename = paste0("Result/metagen/abundance/",i,"_sum.pdf"),width = 6,height = 8)
  
}


write.csv(dat_combine,file = "Result/metagen/abundance/abundance_combine.csv")
# data_2=data_sum_1[rev(order(data_sum_1$`Akkermansia Akkermansia_muciniphila`)),]
# data_2$ID<-rownames(data_2)
# plot_data=melt(data_2,id.vars="ID")
# 
# library(ggsci)
# library(ggplot2)
# 
# feature_idx<-data_2$ID
# plot_data$ID<-factor(plot_data$ID,levels=rev(feature_idx))
# plot_data$variable<-factor(plot_data$variable,levels = species_new)
# pdf(file="taxa-plot\\taxplot_24W.pdf")
# ggplot(plot_data,aes(x=ID,y=value,fill=variable))+geom_bar(stat="identity")+
#   scale_fill_npg()+labs(x="",y="Relative abundance (%)")+
#   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))#+
#  # scale_x_discrete(limits=level)
# dev.off()
# #data_sum_1<-as.data.frame(t(data_sum_1))
# 
# data_plot<-melt(data_sum_1,id.vars = c("species_new"))
# data_plot$species_new<-factor(data_plot$species_new,
#                          levels=rev(unique(data_plot$species_new)))
# 
# pdf(file="taxa-plot\\tax-plot-trimster.pdf")
# ggplot(data_plot,aes(x=variable,y=value,fill=species_new))+
#   geom_bar(stat="identity")+#scale_fill_npg()+
#   labs(y="Relative abundance (%)")
# dev.off()

# 二、pathway ---------------------------------------------------------------

dat_pathway<-read.table("Data/metagen/ko.pathway.raw.txt",
                        header = T,row.names = 1,sep = "\t")

colSums(dat_pathway)


dat_pathway_full<-as.data.frame(t(dat_pathway))
dat_pathway_full$SampleID<-rownames(dat_pathway_full)
#rownames(metadata2)<-metadata2$SampleID
dat_pathway_full<-merge(metadata_new,dat_pathway_full,by= "SampleID")
write.csv(dat_pathway_full,file = "Data/metagen/dat_pathway_full.csv",row.names = T)


transform_RA<-function(x) x*100 / sum(x)
dat_pathway<-dat_pathway
dat_pathway2<-apply(dat_pathway, 2, transform_RA)


abundance = 0.01

idx = rowMeans(dat_pathway2) > abundance
filtered_data = dat_pathway2[idx,]

detectionRate=10

detect<-function(x){
  length(which(x > 0))*100/length(x)
}

DR<-apply(filtered_data, 1, detect)
idx3= apply(filtered_data, 1, detect) >detectionRate

filtered_data2<-filtered_data[idx3,]
pathway_name<-rownames(filtered_data2)

write.csv(filtered_data2,file = "Data/metagen/pathway_filtered.csv",row.names = T)

filtered_data2<-as.data.frame(t(filtered_data2))
filtered_data2$ID<-rownames(filtered_data2)
#rownames(metadata2)<-metadata2$SampleID
dat_pathway2<-merge(metadata2,filtered_data2,by.x = "SampleID",by.y = "ID")
dat_pathway2<-dat_pathway2[,-c(1,3,4)]
dat_pathway2<-merge(metadata_new,dat_pathway2,by="number")
write.csv(dat_pathway2,file = "Data/metagen/dat_pathway.csv",row.names = T)



#1. linear regression -------------------------------------------------------

library(fdrtool)
x_var<-c("TG_24","CHOL_24","TG_32","CHOL_32")
y_var<-paste0("path",seq(1,length(pathway_name),by=1))

cov<-c("pre_BMI","Parity","Age","gesweek_24")
#dir.create("metagen")
dir.create("Result/metagen")

# linear ------------------------------------------------------------------
data_combine<-dat_pathway2
data_combine[pathway_name]<-log(data_combine[pathway_name]+0.001)
#rename column name
colnames(data_combine)[colnames(data_combine) %in% pathway_name] <- y_var
# unadjusted ----------------------------------------------------------------


###unadjusted
result<-data.frame()
for (i in x_var){
  for(j in y_var){
    # xnam_i1<-paste(i,j,sep="*")
    # xnam_i<-paste(j,"pre_BMI","Parity","Age","gesweek_24",sep="+")
    fit_i<-lm(as.formula(paste(i,"~",j)),data=data_combine)
    p_i<-tryCatch({summary(fit_i)$coefficients[2,4]},error=function(e){"NA"})
    b_i<-tryCatch({coef(fit_i)[2]},error=function(e){"NA"})
    blow_i<-tryCatch({confint(fit_i)[2,1]},error=function(e){"NA"})
    bhigh_i<-tryCatch({confint(fit_i)[2,2]},error=function(e){"NA"})
    result_i<-c(i,j,b_i,blow_i,bhigh_i,p_i)  
    result_i<-as.data.frame(result_i) 
    result_i<-as.data.frame(t(result_i)) 
    result<-rbind(result,result_i) 
  }
}
result<-data.frame(pathway_name,result)
names(result)<-c("pathway_name","outcome","pathway","beta","lower","upper","pvalue")   

result_new<-data.frame()
for (i in x_var) {
  result2<-subset(result,outcome==i)
  p<-as.numeric(as.character(result2$pvalue))
  fdr=fdrtool(p,statistic="pvalue",plot = F)
  qval<-data.frame(fdr$qval)
  p2<-p.adjust(p,method = "fdr",length(p))
  p2<-data.frame(p2)
  result2<-cbind(result2,qval,p2)
  result_new<-rbind(result_new,result2)
  
}

group<-rep("Unadjusted",times=nrow(result_new))

result1<-data.frame(result_new,group)



# adjusted ----------------------------------------------------------------



###adjusted
result<-data.frame()
for (i in x_var){
  for(j in y_var){
    # xnam_i1<-paste(i,j,sep="*")
    xnam_i<-paste(j,"pre_BMI","Parity","Age","gesweek_24",sep="+")
    fit_i<-lm(as.formula(paste(i,"~",xnam_i)),data=data_combine)
    p_i<-tryCatch({summary(fit_i)$coefficients[2,4]},error=function(e){"NA"})
    b_i<-tryCatch({coef(fit_i)[2]},error=function(e){"NA"})
    blow_i<-tryCatch({confint(fit_i)[2,1]},error=function(e){"NA"})
    bhigh_i<-tryCatch({confint(fit_i)[2,2]},error=function(e){"NA"})
    result_i<-c(i,j,b_i,blow_i,bhigh_i,p_i)  
    result_i<-as.data.frame(result_i) 
    result_i<-as.data.frame(t(result_i)) 
    result<-rbind(result,result_i) 
  }
}


result<-data.frame(pathway_name,result)
names(result)<-c("pathway_name","outcome","pathway","beta","lower","upper","pvalue")   

result_new<-data.frame()
for (i in x_var) {
  result2<-subset(result,outcome==i)
  p<-as.numeric(as.character(result2$pvalue))
  fdr=fdrtool(p,statistic="pvalue",plot = F)
  qval<-data.frame(fdr$qval)
  p2<-p.adjust(p,method = "fdr",length(p))
  p2<-data.frame(p2)
  result2<-cbind(result2,qval,p2)
  result_new<-rbind(result_new,result2)
  
}


group<-rep("Adjusted",times=nrow(result_new))

result2<-data.frame(result_new,group)

result_full<-rbind(result1,result2)
write.csv(result_full,file="Result/metagen/linear_all_pathway.csv",row.names = F) 


# binary ------------------------------------------------------------------


# uunadjusted -------------------------------------------------------------


###binary
x_var<-c("Hyperlipid_24","TG_24_High","CHOL_24_High",
         "Hyperlipid_32","TG_32_High","CHOL_32_High")


result<-data.frame()
for (i in x_var){
  for(j in y_var){
    # xnam_i1<-paste(i,j,sep="*")
    #xnam_i<-paste(j,"pre_BMI","Parity","Age","gesweek_24",sep="+")
    fit_i<-glm(as.formula(paste(i,"~",j)),family = "binomial",data=data_combine)
    p_i<-tryCatch({summary(fit_i)$coefficients[2,4]},error=function(e){"NA"})
    b_i<-tryCatch({exp(coef(fit_i)[2])},error=function(e){"NA"})
    blow_i<-tryCatch({exp(confint(fit_i)[2,1])},error=function(e){"NA"})
    bhigh_i<-tryCatch({exp(confint(fit_i)[2,2])},error=function(e){"NA"})
    result_i<-c(i,j,b_i,blow_i,bhigh_i,p_i)
    result_i<-as.data.frame(result_i)
    result_i<-as.data.frame(t(result_i))
    result<-rbind(result,result_i)
  }
}
result<-data.frame(pathway_name,result)

names(result)<-c("pathway_name","outcome","pathway","OR","low","upper","pvalueֵ")

result_new<-data.frame()
for (i in x_var) {
  result2<-subset(result,outcome==i)
  p<-as.numeric(as.character(result2$pvalue))
  fdr=fdrtool(p,statistic="pvalue",plot = F)
  qval<-data.frame(fdr$qval)
  p2<-p.adjust(p,method = "fdr",length(p))
  p2<-data.frame(p2)
  result2<-cbind(result2,qval,p2)
  result_new<-rbind(result_new,result2)
  
}
group<-rep("Unadjusted",times=nrow(result_new))

result1<-data.frame(result_new,group)


# adjusted ----------------------------------------------------------------


result<-data.frame()
for (i in x_var){
  for(j in y_var){
    # xnam_i1<-paste(i,j,sep="*")
    xnam_i<-paste(j,"pre_BMI","Parity","Age","gesweek_24",sep="+")
    fit_i<-glm(as.formula(paste(i,"~",xnam_i)),family = "binomial",data=data_combine)
    p_i<-tryCatch({summary(fit_i)$coefficients[2,4]},error=function(e){"NA"})
    b_i<-tryCatch({exp(coef(fit_i)[2])},error=function(e){"NA"})
    blow_i<-tryCatch({exp(confint(fit_i)[2,1])},error=function(e){"NA"})
    bhigh_i<-tryCatch({exp(confint(fit_i)[2,2])},error=function(e){"NA"})
    result_i<-c(i,j,b_i,blow_i,bhigh_i,p_i)
    result_i<-as.data.frame(result_i)
    result_i<-as.data.frame(t(result_i))
    result<-rbind(result,result_i)
  }
}
result<-data.frame(pathway_name,result)

names(result)<-c("pathway_name","outcome","pathway","OR","low","upper","pvalueֵ")


result_new<-data.frame()
for (i in x_var) {
  result2<-subset(result,outcome==i)
  p<-as.numeric(as.character(result2$pvalue))
  fdr=fdrtool(p,statistic="pvalue",plot = F)
  qval<-data.frame(fdr$qval)
  p2<-p.adjust(p,method = "fdr",length(p))
  p2<-data.frame(p2)
  result2<-cbind(result2,qval,p2)
  result_new<-rbind(result_new,result2)
  
}


group<-rep("Adjusted",times=nrow(result_new))

result2<-data.frame(result_new,group)

result_full<-rbind(result1,result2)

write.csv(result_full,file="Result/metagen/linear_all_pathway_group.csv",row.names = F)



# plot --------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(ggsci)
library(openxlsx)
library(dplyr)
setwd("E:/study2/microbiome/study_new/total/9Data Review")
data_plot<-read.xlsx("Result/metagen/linear_all_pathway.xlsx", sheet = 1)
data_plot<-subset(data_plot,group=="Adjusted")
data_plot_new1<-subset(data_plot,pvalue<0.05)
pathway_selected<-as.character(unique(data_plot_new1$pathway_name))
data_plot2<-subset(data_plot,pathway_name%in%pathway_selected)

#data_plot2$name<-paste(data_plot2$pathway_name,data_plot2$Description1,sep = ":")
data_plot2$logp<- -log10(data_plot2$pvalue)
# dotplot <- 
#   ggplot(data_plot2,aes(x=outcome, y = name, color = beta, size = logp)) + 
#   
#   geom_point()+
#   geom_point(shape=1,color="black") + 
#   cowplot::theme_cowplot() + 
#   theme(axis.line  = element_blank()) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#   xlab('')+
#   ylab('') +
#   theme(axis.ticks = element_blank()) +
#   scale_color_gradient2(low="blue", high = "red", name = 'regression coefficient')+
#   #  scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,10), oob = scales::squish, name = '-log (P value)') +
#   scale_y_discrete(position = "right")


# bl <- colorRampPalette(c("navy","royalblue","lightskyblue"))(200)                      
# re <- colorRampPalette(c("mistyrose", "red2","darkred"))(200)
data_plot2$outcome<-ordered(data_plot2$outcome,levels=c("TG_24","CHOL_24","TG_32","CHOL_32"))
data_plot2$label<-ifelse(data_plot2$pvalue<0.05,round(data_plot2$pvalue,3),"")


getSig <- function(dc) {
  sc <- ''
  if (dc < 0.001) sc <- '***'
  else if (dc < 0.01) sc <- '**'
  else if (dc < 0.05) sc <- '*'
  sc
}

data_plot2$label2<-sapply(data_plot2$pvalue,getSig)
library("scales")
p<-ggplot(data_plot2, aes( pathway_name.1,outcome)) +
  geom_tile(aes(fill = beta)) +
  #  scale_fill_gradient2(low="blue", high = "red", name = 'regression coefficient')+
  scale_fill_gradientn(colours=c("#4DBBD5","white", "#E64B35"), 
                       values = rescale(c(-.8,0,.5)),
                       guide = "colorbar", limits=c(-.8,.5),
                       na.value = "grey98")+
  geom_text(aes(label=label2))+
  labs(x="",y="")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(p,filename = "Result/metagen/plot_pathway.pdf",width = 9,height = 5)


###heatmap
dat_heat_beta<-dcast(data_plot2[,c("pathway_name.1","outcome","beta")],pathway_name.1~outcome,value.var = "beta")
rownames(dat_heat_beta)<-dat_heat_beta$pathway_name.1
dat_heat_beta<-dat_heat_beta[,-1]
dat_heat_beta<-as.matrix(dat_heat_beta)

dat_heat_p<-dcast(data_plot2[,c("pathway_name.1","outcome","pvalue")],pathway_name.1~outcome,value.var = "pvalue")
rownames(dat_heat_p)<-dat_heat_p$pathway_name.1
dat_heat_p<-dat_heat_p[,-1]
dat_heat_p<-as.matrix(dat_heat_p)

row_anno<-subset(data_plot2,outcome=="TG_24")[,c("pathway_name.1")]
# rownames(row_anno)<-row_anno$KO_name
# row_anno<-row_anno[,c("PathwayL1","PathwayL2")]
#row_anno<-as.data.frame(row_anno)
# ann_colors = list(
#   # Season=c(Summer = "#8491B4", Winter = "#91D1C2"),
#   PathwayL1 = c(`Brite Hierarchies` = "#E64B35", 
#                 `Genetic Information Processing` = "#4DBBD5", 
#                 Metabolism="#91D1C2",
#                 `Not Included in Pathway or Brite`="F39B7F",
#                 `Organismal Systems`="#8491B4"),
#   PathwayL2=c(`HNSC_HPV-` ="#F39B7F",`HNSC_HPV+` ="#3C5488"))

getSig <- function(dc) {
  sc <- ''
  if (dc < 0.001) sc <- '***'
  else if (dc < 0.01) sc <- '**'
  else if (dc < 0.05) sc <- '*'
  sc
}
sig.mat <- matrix(sapply(dat_heat_p, getSig), nrow=nrow(dat_heat_p))
str(sig.mat)
# pheatmap::pheatmap(dat_heat_beta, scale="none",
#                    cluster_rows=T,
#                    #  breaks = breaksList,
#                    #  annotation_col=row_anno,
#                    color=colorRampPalette(c("#4DBBD5","white","#E64B35"))(6),
#                    #  annotation_colors = ann_colors[1],
#                    display_numbers = sig.mat,border=FALSE,
#                    angle_col="270",
#                    filename="Result/metagen/pheatmap_anno_filter.pdf",width = 18,height = 8)

library(pheatmap)
library(gplots)
library(RColorBrewer)
hmcol <- function(n){
  colfun1 <- colorRampPalette(rev(brewer.pal(9, "Blues")))
  colfun2 <- colorRampPalette(brewer.pal(9, "Reds"))
  p <- (0 - min(mat)) / (max(mat) - min(mat))
  gap <- n * p
  c(colfun1(floor(gap)), colfun2(n - gap))
}
mat<-dat_heat_beta
p<-pheatmap(mat, scale="none",cluster_row = T,cluster_cols = T,
            display_numbers = sig.mat,
            fontface="italic",
            #  annotation_col = 
          #  color=colorRampPalette(c("#4DBBD5","white","#E64B35"))(6),
            color=hmcol(50),
            border=FALSE)
ggsave(p,filename = "Result/metagen/plot_path_heatmap2.pdf",width = 8,height = 8)
