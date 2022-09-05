setwd('E:/study2/microbiome/study_new/total/9Data Review')
library(Maaslin2)
library(ggplot2)
library(lme4)
library(nlme)
library(fdrtool)
library(openxlsx)
genus<-read.xlsx("Data/16S/genus_filter.xlsx",sheet = 2)
rownames(genus)<-genus$Genus
genus_name<-genus$Genus
genus<-genus[,-c(1:6)]
data<-read.xlsx("Data/data_total.xlsx",sheet =2)
# combine data ------------------------------------------------------------
genus_filter<-as.data.frame(t(genus))

genus_filter<-genus_filter+1
genus_filter<-genus_filter/(32007+98)*100

genus_filter_log<-log10(genus_filter)
genus_filter_log$ID<-rownames(genus_filter_log)
data_all<-merge(genus_filter_log,data,by.x = "ID",by.y = "ID")

x_var<-c("TG","CHOL")
y_var<-genus_name

cov<-c("pre_BMI","Parity","Age","gesweek_seq")

# linear mixed model ------------------------------------------------------

dir.create("Result/LME")

# Unadjsuted --------------------------------------------------------------


data_all$Time<-as.factor(data_all$Time)
result<-data.frame()
for (i in x_var){
  for(j in y_var){
    # xnam_i1<-paste(i,j,sep="*")
  #  xnam_i<-paste(j,"pre_BMI","Parity","Age","Time",sep="+")
    fit_i<-lme(as.formula(paste(i,"~",j)),random = ~ 1|number,data=data_all)
    p_i<-tryCatch({summary(fit_i)$tTable[2,5]},error=function(e){"NA"})
    b_i<-tryCatch({summary(fit_i)$tTable[2,1]},error=function(e){"NA"})
    
    se_i<-tryCatch({summary(fit_i)$tTable[2,2]},error=function(e){"NA"})
    blow_i<-round(b_i-1.96*se_i,3)
    bhigh_i<-round(b_i+1.96*se_i,3)
    b_i<-round(b_i,3)
    
    b<-paste0(b_i," (",blow_i,",",bhigh_i,")")
    result_i<-c(i,j,b_i,blow_i,bhigh_i,b,p_i)  
    result_i<-as.data.frame(result_i) 
    result_i<-as.data.frame(t(result_i)) 
    result<-rbind(result,result_i) 
  }
}

names(result)<-c("variable1","variable2","beta","lower","upper","CI","pvalue")   
result_TG<-subset(result,variable1=="TG")
result_CHOL<-subset(result,variable1=="CHOL")
p<-as.numeric(as.character(result_TG$pvalue))

fdr=fdrtool(p,statistic="pvalue",plot = F)
qval<-data.frame(fdr$qval)
p2<-p.adjust(p,method = "fdr",length(p))
p2<-data.frame(p2)
result_TG<-cbind(result_TG,qval,p2)

p<-as.numeric(as.character(result_CHOL$pvalue))

fdr=fdrtool(p,statistic="pvalue",plot = F)
qval<-data.frame(fdr$qval)
p2<-p.adjust(p,method = "fdr",length(p))
p2<-data.frame(p2)
result_CHOL<-cbind(result_CHOL,qval,p2)

result<-rbind(result_TG,result_CHOL)
group<-rep("Unadjusted",times=nrow(result))
result1<-data.frame(result,group)


# adjusted ----------------------------------------------------------------

data_all$Time<-as.factor(data_all$Time)
result<-data.frame()
for (i in x_var){
  for(j in y_var){
    # xnam_i1<-paste(i,j,sep="*")
    xnam_i<-paste(j,"pre_BMI","Parity","Age","Time",sep="+")
    fit_i<-lme(as.formula(paste(i,"~",xnam_i)),random = ~ 1|number,data=data_all)
    p_i<-tryCatch({summary(fit_i)$tTable[2,5]},error=function(e){"NA"})
    b_i<-tryCatch({summary(fit_i)$tTable[2,1]},error=function(e){"NA"})
    
    se_i<-tryCatch({summary(fit_i)$tTable[2,2]},error=function(e){"NA"})
    blow_i<-round(b_i-1.96*se_i,3)
    bhigh_i<-round(b_i+1.96*se_i,3)
    b_i<-round(b_i,3)
    
    b<-paste0(b_i," (",blow_i,",",bhigh_i,")")
    result_i<-c(i,j,b_i,blow_i,bhigh_i,b,p_i)  
    result_i<-as.data.frame(result_i) 
    result_i<-as.data.frame(t(result_i)) 
    result<-rbind(result,result_i) 
  }
}

names(result)<-c("variable1","variable2","beta","lower","upper","CI","pvalue")   
result_TG<-subset(result,variable1=="TG")
result_CHOL<-subset(result,variable1=="CHOL")
p<-as.numeric(as.character(result_TG$pvalue))

fdr=fdrtool(p,statistic="pvalue",plot = F)
qval<-data.frame(fdr$qval)
p2<-p.adjust(p,method = "fdr",length(p))
p2<-data.frame(p2)
result_TG<-cbind(result_TG,qval,p2)

p<-as.numeric(as.character(result_CHOL$pvalue))

fdr=fdrtool(p,statistic="pvalue",plot = F)
qval<-data.frame(fdr$qval)
p2<-p.adjust(p,method = "fdr",length(p))
p2<-data.frame(p2)
result_CHOL<-cbind(result_CHOL,qval,p2)

result<-rbind(result_TG,result_CHOL)
group<-rep("Adjusted",times=nrow(result))
result2<-data.frame(result,group)



result_full<-rbind(result1,result2)
write.csv(result_full,file="Result/LME/lmer_lipid_log10_new.csv",row.names = F) 





# draw regression plot ----------------------------------------------------


dir.create("Result/LME/regplot")

plot_reg<-function(i,j){
  
  p<-ggplot(data=data_all,aes(y =data_all[,i],x =data_all[,j],color=Time))+
    geom_point(size=1.5,alpha=0.3)+
  #  geom_line(method="loess")+
    geom_smooth(method="lm")+#geom_jitter()+
    labs(y =i, x =j)+
    theme_bw()+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
    #scale_y_continuous(limits=c(0.5, 0.9))+
    #scale_color_npg()+
    theme(axis.text.x = element_text(face = "italic"))+
    scale_color_manual(values = c("#00AFBB", "#E7B800"))
  ggsave(p,filename = paste0("Result/LME/regplot/",i,"_",j,".pdf"),width = 6,height = 5)
  
}

for (i in c("TG","CHOL")) {
  for (j in genus_name) {
    plot_reg(i,j)
  }
  
}



#sig_result
dir.create("Result/LME/regplot2")

result_new<-subset(result,fdr.qval<0.2)
genus_new<-as.character(result_new$variable2)

outco<-as.character(unique(result_new$variable1))

plot_reg<-function(i,j){
  
  p<-ggplot(data=data_all,aes(y =data_all[,i],x =data_all[,j],color=Time))+
    geom_point(size=1.5,alpha=0.3)+
    #  geom_line(method="loess")+
    geom_smooth(method="lm")+#geom_jitter()+
    labs(y =i, x =j)+
    theme_bw()+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
    #scale_y_continuous(limits=c(0.5, 0.9))+
    #scale_color_npg()+
    theme(axis.text.x = element_text(face = "italic"))+
    scale_color_manual(values = c("#00AFBB", "#E7B800"))
  ggsave(p,filename = paste0("Result/LME/regplot2/",i,"_",j,".pdf"),width = 6,height = 5)
  
}

for (i in outco) {
  for (j in genus_new) {
    plot_reg(i,j)
  }
  
}




# draw regression plot2 ----------------------------------------------------


dir.create("Result/LME/regplot_all")

plot_reg<-function(i,j){
  
  p<-ggplot(data=data_all,aes(y =data_all[,i],x =data_all[,j]))+
    geom_point(size=2,alpha=0.3,color='#2980B9')+
    #  geom_line(method="loess")+
    geom_smooth(method="lm",color='#2980B9')+#geom_jitter()+
    labs(y =i, x =j)+
    theme_bw()+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
    theme(axis.text.x = element_text(face = "italic"))
    #scale_y_continuous(limits=c(0.5, 0.9))+
    #scale_color_npg()+
   # scale_color_manual(values = c("#00AFBB", "#E7B800"))
  ggsave(p,filename = paste0("Result/LME/regplot_all/",i,"_",j,".pdf"),width = 5,height = 5)
  
}

for (i in c("TG","CHOL")) {
  for (j in genus_name) {
    plot_reg(i,j)
  }
  
}



#sig_result
dir.create("Result/LME/regplot_sig")

result_new<-subset(result,fdr.qval<0.2)
genus_new<-as.character(result_new$variable2)

outco<-as.character(unique(result_new$variable1))

plot_reg<-function(i,j){
  
  p<-ggplot(data=data_all,aes(y =data_all[,i],x =data_all[,j]))+
    geom_point(size=2,alpha=0.3,color='#2980B9')+
    #  geom_line(method="loess")+
    geom_smooth(method="lm",color='#2980B9')+#geom_jitter()+
    labs(y =i, x =j)+
    theme_bw()+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
    theme(axis.text.x = element_text(face = "italic"))
    #scale_y_continuous(limits=c(0.5, 0.9))+
    #scale_color_npg()+
    scale_color_manual(values = c("#00AFBB", "#E7B800"))
  ggsave(p,filename = paste0("Result/LME/regplot_sig/",i,"_",j,".pdf"),width = 5,height = 5)
  
}

for (i in outco) {
  for (j in genus_new) {
    plot_reg(i,j)
  }
  
}



# draw plot 3 -------------------------------------------------------------

dir.create("Result/LME/regplot_sig2")

require(ggplot2)
require(reshape2)
library(dplyr)
library(openxlsx)
genus_select<-read.xlsx("Result/LME/result_combine.xlsx",sheet = 8)
genus_select<-as.character(na.omit(genus_select$intersect))

###TG
mtcars2 = data_all[,c("TG",genus_select)]
mtcars2<-melt(mtcars2,id.vars = "TG")
iris_groups<- group_by(mtcars2, variable)


df_summarise<- summarise(iris_groups, mean(TG), mean(value))

lab = mean(mtcars2$TG)

#----我们使用的都是同样的y做相关,所以呢都是相同的
ggplot2::ggplot(mtcars2,aes(value,TG)) +
  ggplot2::geom_point(color='#2980B9') +
  ggpubr::stat_cor(label.y=lab*1.1)+
  ggpubr::stat_regline_equation(label.y=lab) +
  facet_wrap(~variable, scales="free_x",nrow = 1) +
  geom_smooth(aes(value,TG), color='#2980B9',method=lm, se=T)+
  theme_grey()+
  theme_bw()


ggsave(filename = "Result/LME/regplot_sig2/TG_combine.pdf",width = 15,height = 4)
  
###CHOL
mtcars2 = data_all[,c("CHOL",genus_select)]
mtcars2<-melt(mtcars2,id.vars = "CHOL")
iris_groups<- group_by(mtcars2, variable)


df_summarise<- summarise(iris_groups, mean(CHOL), mean(value))

lab = mean(mtcars2$CHOL)

#----我们使用的都是同样的y做相关,所以呢都是相同的
ggplot2::ggplot(mtcars2,aes(value,CHOL)) +
  ggplot2::geom_point(color='#2980B9') +
  ggpubr::stat_cor(label.y=lab*1.1)+
  ggpubr::stat_regline_equation(label.y=lab) +
  facet_wrap(~variable, scales="free_x",nrow = 1) +
  geom_smooth(aes(value,CHOL), color='#2980B9',method=lm, se=T)+
  theme_grey()+
  theme_bw()


ggsave(filename = "Result/LME/regplot_sig2/CHOL_combine.pdf",width = 15,height = 4)


# maAslin ------------------------------------------------------------------
rownames(data_all)<-data_all$ID
dir.create("Result/LME/maAslin_log")
data_all$CHOL_High<-as.factor(data_all$CHOL_High)
data_all$TG_High<-as.factor(data_all$TG_High)

#dir.create("maAslin")
maaslin_out<-function(i){
  Maaslin2(
    genus_filter_log, data_all[c("number",i,cov)], 
    # min_prevalence = 0.1,
    # min_abundance = 0.01,
    paste("Result/LME/maAslin_log/",i,sep = ""),
    #  'demo_output_maaslin_GDM', #transform = "AST",
    fixed_effects = c("Age","pre_BMI","gesweek_seq",
                      "Parity",i),
    random_effects = c('number'),
    normalization = 'NONE',
    analysis_method = "LM",
    max_significance=0.2,
    transform= "NONE",
    correction = "BH",
    standardize = FALSE)
}

for (i in c("TG","CHOL","TG_High","CHOL_High","Hyperlipid","dyslipimedia")){
  maaslin_out(i)
}


rownames(data_all)<-data_all$ID

data_all$CHOL_High<-as.factor(data_all$CHOL_High)
data_all$TG_High<-as.factor(data_all$TG_High)

dir.create("Result/LME/maAslin")
maaslin_out<-function(i){
  Maaslin2(
    genus_filter, data_all[c("number","TG","TG_High","CHOL","CHOL_High","Hyperlipid","dyslipimedia",cov)], 
    # min_prevalence = 0.1,
    # min_abundance = 0.01,
    paste("Result/LME/maAslin/",i,sep = ""),
    #  'demo_output_maaslin_GDM', #transform = "AST",
    fixed_effects = c("Age","pre_BMI","gesweek_seq",
                      "Parity",i),
    random_effects = c('number'),
    normalization = 'NONE',
    analysis_method = "LM",
    max_significance=0.2,
    transform= "NONE",
    correction = "BH",
    standardize = FALSE)
}

for (i in c("TG","CHOL","TG_High","CHOL_High","Hyperlipid","dyslipimedia")){
  maaslin_out(i)
}


# maaslin_new -------------------------------------------------------------

# 
# dir.create("maAslin_new")
# maaslin_out<-function(i){
#   Maaslin2(
#     data_all[c("TG","CHOL")], data_all[c("number",cov,genus_name)], 
#     # min_prevalence = 0.1,
#     # min_abundance = 0.01,
#     paste("maAslin_new/",i,sep = ""),
#     #  'demo_output_maaslin_GDM', #transform = "AST",
#     fixed_effects = c("Age","pre_BMI","gesweek_seq",
#                       "Parity",i),
#     random_effects = c('number'),
#     normalization = 'NONE',
#     analysis_method = "LM",
#     max_significance=0.2,
#     transform= "NONE",
#     correction = "BH",
#     standardize = FALSE)
# }
# 
# for (i in genus_name){
#   maaslin_out(i)
# }



#Part 2. arcsine transformation -----------------------------------------------------------

# arcsine transformation in r
#> asin(sqrt(0.5))


genus_filter<-as.data.frame(t(genus))
#colnames(genus_filter)<-as.character(genus$Genus)
genus_filter<-genus_filter+1
genus_filter<-genus_filter/(32007+98)

# for (i in 1:ncol(genus_filter)){
#   genus_filter[,i][genus_filter[,i]==0] <-1/32007*100
# }

arc<-function(x){
  asin(sqrt(x))
}
genus_filter_arc<-apply(genus_filter,2,arc)
genus_filter_arc<-as.data.frame(genus_filter_arc)

#genus_filter_log<-log10(genus_filter)
genus_filter_arc$ID<-rownames(genus_filter_arc)
data_all<-merge(genus_filter_arc,data,by.x = "ID",by.y = "ID")

x_var<-c("TG","CHOL")
y_var<-genus_name

cov<-c("pre_BMI","Parity","Age","gesweek_seq")

# linear mixed model ------------------------------------------------------

data_all$Time<-as.factor(data_all$Time)
result<-data.frame()
for (i in x_var){
  for(j in y_var){
    # xnam_i1<-paste(i,j,sep="*")
    xnam_i<-paste(j,"pre_BMI","Parity","Age","Time",sep="+")
    fit_i<-lme(as.formula(paste(i,"~",xnam_i)),random = ~ 1|number,data=data_all)
    p_i<-tryCatch({summary(fit_i)$tTable[2,5]},error=function(e){"NA"})
    b_i<-tryCatch({summary(fit_i)$tTable[2,1]},error=function(e){"NA"})
    
    se_i<-tryCatch({summary(fit_i)$tTable[2,2]},error=function(e){"NA"})
    blow_i<-round(b_i-1.96*se_i,3)
    bhigh_i<-round(b_i+1.96*se_i,3)
    b_i<-round(b_i,3)
    
    b<-paste0(b_i," (",blow_i,",",bhigh_i,")")
    result_i<-c(i,j,b_i,blow_i,bhigh_i,b,p_i)  
    result_i<-as.data.frame(result_i) 
    result_i<-as.data.frame(t(result_i)) 
    result<-rbind(result,result_i) 
  }
}

#names(result)<-c("variable1","variable2","beta","lower","upper","CI","pvalue")    
names(result)<-c("variable1","variable2","beta","lower","upper","CI","pvalue")   
result_TG<-subset(result,variable1=="TG")
result_CHOL<-subset(result,variable1=="CHOL")
p<-as.numeric(as.character(result_TG$pvalue))

fdr=fdrtool(p,statistic="pvalue",plot = F)
qval<-data.frame(fdr$qval)
p2<-p.adjust(p,method = "fdr",length(p))
p2<-data.frame(p2)
result_TG<-cbind(result_TG,qval,p2)

p<-as.numeric(as.character(result_CHOL$pvalue))

fdr=fdrtool(p,statistic="pvalue",plot = F)
qval<-data.frame(fdr$qval)
p2<-p.adjust(p,method = "fdr",length(p))
p2<-data.frame(p2)
result_CHOL<-cbind(result_CHOL,qval,p2)

result<-rbind(result_TG,result_CHOL)

write.csv(result,file="Result/LME/lmer_lipid_arc.csv",row.names = F) 

# maAslin ------------------------------------------------------------------
rownames(data_all)<-data_all$ID
dir.create("Result/LME/maAslin_arc")
data_all$CHOL_High<-as.factor(data_all$CHOL_High)
data_all$TG_High<-as.factor(data_all$TG_High)

maaslin_out<-function(i){
  Maaslin2(
    genus_filter_log, data_all[c("number",i,cov)], 
    # min_prevalence = 0.1,
    # min_abundance = 0.01,
    paste("Result/LME/maAslin_arc/",i,sep = ""),
    #  'demo_output_maaslin_GDM', #transform = "AST",
    fixed_effects = c("Age","pre_BMI","gesweek_seq",
                      "Parity",i),
    random_effects = c('number'),
    normalization = 'NONE',
    analysis_method = "LM",
    max_significance=0.2,
    transform= "NONE",
    correction = "BH",
    standardize = FALSE)
}

for (i in c("TG","CHOL","TG_High","CHOL_High","Hyperlipid","dyslipimedia")){
  maaslin_out(i)
}


rownames(data_all)<-data_all$ID

data_all$CHOL_High<-as.factor(data_all$CHOL_High)
data_all$TG_High<-as.factor(data_all$TG_High)



# maaslin_new_arc -------------------------------------------------------------

# 
# dir.create("maAslin_new_arc")
# maaslin_out<-function(i){
#   Maaslin2(
#     data_all[c("TG","CHOL")], data_all[c("number",cov,genus_name)], 
#     # min_prevalence = 0.1,
#     # min_abundance = 0.01,
#     paste("maAslin_new_arc/",i,sep = ""),
#     #  'demo_output_maaslin_GDM', #transform = "AST",
#     fixed_effects = c("Age","pre_BMI","gesweek_seq",
#                       "Parity",i),
#     random_effects = c('number'),
#     normalization = 'NONE',
#     analysis_method = "LM",
#     max_significance=0.2,
#     transform= "NONE",
#     correction = "BH",
#     standardize = FALSE)
# }
# 
# for (i in genus_name){
#   maaslin_out(i)
# }


# part 3 combine result------------------------------------------------------------------

dir<-'E:/study2/microbiome/study_new/total/9Data Review/Result/LME'
# folder<-c("maAslin","maAslin_log","maAslin_new",
#           "maAslin_arc","maAslin_new_arc"
# )

folder<-c("maAslin","maAslin_log",
          "maAslin_arc"
)

sig_result<-data.frame()
for (i in folder) {
  dir_i<-paste(dir,i,sep = "/")
  setwd(dir_i)
  file_2i<-list.files(dir_i)
  for (j in file_2i){
    dir_j<-paste(dir_i,j,sep = "/")
    setwd(dir_j)
    sigresult_j<-read.table("significant_results.tsv",sep = '\t', header = TRUE)
    sigresult_j<-subset(sigresult_j,value==j)
    group<-rep(i,times=nrow(sigresult_j))
    sigresult_j<-data.frame(sigresult_j,group)
    sig_result<-rbind(sig_result,sigresult_j)
  }
  
}

setwd(dir)
dir.create("maAslin_combine3")

write.csv(sig_result,file = "maAslin_combine3/sigresult_full.csv")





#dir<-"E:/study2/microbiome/study_new/total/Linear mixed model"
# folder<-c("maAslin","maAslin_log","maAslin_new",
#           "maAslin_arc","maAslin_new_arc"
# )

folder<-c("maAslin","maAslin_log",
          "maAslin_arc"
)

sig_result<-data.frame()
for (i in folder) {
  dir_i<-paste(dir,i,sep = "/")
  setwd(dir_i)
  file_2i<-list.files(dir_i)
  for (j in file_2i){
    dir_j<-paste(dir_i,j,sep = "/")
    setwd(dir_j)
    sigresult_j<-read.table("all_results.tsv",sep = '\t', header = TRUE)
    sigresult_j<-subset(sigresult_j,value==j)
    group<-rep(i,times=nrow(sigresult_j))
    sigresult_j<-data.frame(sigresult_j,group)
    sig_result<-rbind(sig_result,sigresult_j)
  }
  
}

setwd(dir)
dir.create("maAslin_combine3")

write.csv(sig_result,file = "maAslin_combine3/allresult_full.csv")


#Part 4. linear regression -------------------------------------------------------

source_dir<-"E:/study2/microbiome/study_new/total/9Data Review"
setwd(source_dir)
#setwd("E:/study2/microbiome/study_new/total")
library(Maaslin2)
library(ggplot2)
library(lme4)
library(nlme)
library(fdrtool)
library(openxlsx)
genus<-read.xlsx("Data/16S/genus_filter.xlsx",sheet = 2)
rownames(genus)<-genus$Genus
genus_name<-genus$Genus
genus<-genus[,-c(1:6)]

# 1.total -------------------------------------------------------------------
data<-read.xlsx("Data/data_total.xlsx",sheet =2)

genus_filter<-as.data.frame(t(genus))
#colnames(genus_filter)<-as.character(genus$Genus)
#genus_filter<-genus_filter/32007*100

#write.csv(genus_filter,file = "Data/16S/genus_filter_RA.csv")

# for (i in 1:ncol(genus_filter)){
#   genus_filter[,i][genus_filter[,i]==0] <-1/32007*100
# }
genus_filter<-genus_filter+1
genus_filter<-genus_filter/(32007+98)*100

genus_filter_log<-log10(genus_filter)
genus_filter_log$ID<-rownames(genus_filter_log)
# 
# ID_24<-data$ID_24
# genus_filter_log<-genus_filter_log[ID_24,]

data_all<-merge(genus_filter_log,data,by.x = "ID",by.y = "ID")

x_var<-c("TG","CHOL")
y_var<-genus_name

test_name<-c("BUN", "CREA", "UA", "TP", "ALB", "CHOL", "TG", "ALT", "AST", 
             "ALP", "LDH", "GGT", "TBIL", "DBIL", "TBA", "FMN", "FBG", "RBP", 
             "cystain", "HCY")

#impute missing value
library(dplyr)
# data_rf_second<-data_rf_second %>% 
#   mutate_all(function(x){x[is.na(x)] <- mean(x)
#   x})

for (i in test_name) {
  
  data_all[,i][is.na(data_all[,i])==TRUE]<-mean(data_all[,i],na.rm = T)
}
cov<-c("pre_BMI","Parity","Age","gesweek_seq")
dir.create("Result/LME/association")


#1.1 linear ------------------------------------------------------------------


# unadjusted ----------------------------------------------------------------


###unadjusted
result<-data.frame()
for (i in x_var){
  for(j in y_var){
    # xnam_i1<-paste(i,j,sep="*")
    # xnam_i<-paste(j,"pre_BMI","Parity","Age","gesweek_seq",sep="+")
    fit_i<-lm(as.formula(paste(i,"~",j)),data=data_all)
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

names(result)<-c("outcome","genus","beta","lower","upper","pvalue")   

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
    xnam_i<-paste(j,"pre_BMI","Parity","Age","gesweek_seq",sep="+")
    fit_i<-lm(as.formula(paste(i,"~",xnam_i)),data=data_all)
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


names(result)<-c("outcome","genus","beta","lower","upper","pvalue")   

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
write.csv(result_full,file="Result/LME/association/linear_all.csv",row.names = F) 


#1.2 binary ------------------------------------------------------------------


# uunadjusted -------------------------------------------------------------


###binary
x_var<-c("Hyperlipid","TG_High","CHOL_High")


result<-data.frame()
for (i in x_var){
  for(j in y_var){
    # xnam_i1<-paste(i,j,sep="*")
    #xnam_i<-paste(j,"pre_BMI","Parity","Age","gesweek_seq",sep="+")
    fit_i<-glm(as.formula(paste(i,"~",j)),family = "binomial",data=data_all)
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

names(result)<-c("outcome","genus","OR","low","upper","pvalueֵ")

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
    xnam_i<-paste(j,"pre_BMI","Parity","Age","gesweek_seq",sep="+")
    fit_i<-glm(as.formula(paste(i,"~",xnam_i)),family = "binomial",data=data_all)
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

names(result)<-c("outcome","genus","OR","low","upper","pvalueֵ")


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

write.csv(result_full,file="Result/LME/association/linear_all_group.csv",row.names = F)





#2. T2 and T2 ---------------------------------------------------------------

data<-read.xlsx("Data/data_total.xlsx",sheet =1)

genus_filter<-as.data.frame(t(genus))
#colnames(genus_filter)<-as.character(genus$Genus)
#genus_filter<-genus_filter/32007*100

#write.csv(genus_filter,file = "Data/16S/genus_filter_RA.csv")

# for (i in 1:ncol(genus_filter)){
#   genus_filter[,i][genus_filter[,i]==0] <-1/32007*100
# }
genus_filter<-genus_filter+1
genus_filter<-genus_filter/(32007+98)*100

genus_filter_log<-log10(genus_filter)
genus_filter_log$ID<-rownames(genus_filter_log)

ID_24<-data$ID_24
genus_filter_log<-genus_filter_log[ID_24,]

data_all<-merge(genus_filter_log,data,by.x = "ID",by.y = "ID_24")

x_var<-c("TG_24","CHOL_24")
y_var<-genus_name

cov<-c("pre_BMI","Parity","Age","gesweek_seq")
dir.create("Result/LME/association")

#2.1 linear ------------------------------------------------------------------


# unadjusted --------------------------------------------------------------


###unadjusted
result<-data.frame()
for (i in x_var){
  for(j in y_var){
    # xnam_i1<-paste(i,j,sep="*")
    xnam_i<-paste(j,"pre_BMI","Parity","Age","gesweek_24",sep="+")
    fit_i<-lm(as.formula(paste(i,"~",j)),data=data_all)
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


names(result)<-c("outcome","genus","beta","lower","upper","pvalue")   

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



###Adjusted
result<-data.frame()
for (i in x_var){
  for(j in y_var){
    # xnam_i1<-paste(i,j,sep="*")
    xnam_i<-paste(j,"pre_BMI","Parity","Age","gesweek_24",sep="+")
    fit_i<-lm(as.formula(paste(i,"~",xnam_i)),data=data_all)
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


names(result)<-c("outcome","genus","beta","lower","upper","pvalue")   

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
write.csv(result_full,file="Result/LME/association/linear_T2.csv",row.names = F) 



#2.2 binary ------------------------------------------------------------------


# unadjusted --------------------------------------------------------------


x_var<-c("Hyperlipid_24","TG_24_High","CHOL_24_High")
result<-data.frame()
for (i in x_var){
  for(j in y_var){
    # xnam_i1<-paste(i,j,sep="*")
    xnam_i<-paste(j,"pre_BMI","Parity","Age","gesweek_24",sep="+")
    fit_i<-glm(as.formula(paste(i,"~",j)),family = "binomial",data=data_all)
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

names(result)<-c("outcome","genus","OR","low","upper","pvalueֵ")


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

x_var<-c("Hyperlipid_24","TG_24_High","CHOL_24_High")
result<-data.frame()
for (i in x_var){
  for(j in y_var){
    # xnam_i1<-paste(i,j,sep="*")
    xnam_i<-paste(j,"pre_BMI","Parity","Age","gesweek_24",sep="+")
    fit_i<-glm(as.formula(paste(i,"~",xnam_i)),family = "binomial",data=data_all)
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

names(result)<-c("outcome","genus","OR","low","upper","pvalueֵ")


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

write.csv(result_full,file="Result/LME/association/linear_T2_group.csv",row.names = F)





#3. T3 and T3 ---------------------------------------------------------------

data<-read.xlsx("Data/data_total.xlsx",sheet =1)

genus_filter<-as.data.frame(t(genus))
#colnames(genus_filter)<-as.character(genus$Genus)
#genus_filter<-genus_filter/32007*100

#write.csv(genus_filter,file = "Data/16S/genus_filter_RA.csv")

# for (i in 1:ncol(genus_filter)){
#   genus_filter[,i][genus_filter[,i]==0] <-1/32007*100
# }
genus_filter<-genus_filter+1
genus_filter<-genus_filter/(32007+98)*100

genus_filter_log<-log10(genus_filter)
genus_filter_log$ID<-rownames(genus_filter_log)

ID_32<-data$ID_32
genus_filter_log<-genus_filter_log[ID_32,]

data_all<-merge(genus_filter_log,data,by.x = "ID",by.y = "ID_32")

x_var<-c("TG_32","CHOL_32")
y_var<-genus_name

cov<-c("pre_BMI","Parity","Age","gesweek_seq")
dir.create("Result/LME/association")

#3.1 linear ------------------------------------------------------------------


# unadjusted --------------------------------------------------------------


###unadjusted
result<-data.frame()
for (i in x_var){
  for(j in y_var){
    # xnam_i1<-paste(i,j,sep="*")
    xnam_i<-paste(j,"pre_BMI","Parity","Age","gesweek_32",sep="+")
    fit_i<-lm(as.formula(paste(i,"~",j)),data=data_all)
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


names(result)<-c("outcome","genus","beta","lower","upper","pvalue")   

result_new<-data.frame()
for (i in x_var) {
  resulT3<-subset(result,outcome==i)
  p<-as.numeric(as.character(resulT3$pvalue))
  fdr=fdrtool(p,statistic="pvalue",plot = F)
  qval<-data.frame(fdr$qval)
  p2<-p.adjust(p,method = "fdr",length(p))
  p2<-data.frame(p2)
  resulT3<-cbind(resulT3,qval,p2)
  result_new<-rbind(result_new,resulT3)
  
}
group<-rep("Unadjusted",times=nrow(result_new))

result1<-data.frame(result_new,group)


# adjusted ----------------------------------------------------------------



###Adjusted
result<-data.frame()
for (i in x_var){
  for(j in y_var){
    # xnam_i1<-paste(i,j,sep="*")
    xnam_i<-paste(j,"pre_BMI","Parity","Age","gesweek_32",sep="+")
    fit_i<-lm(as.formula(paste(i,"~",xnam_i)),data=data_all)
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


names(result)<-c("outcome","genus","beta","lower","upper","pvalue")   

result_new<-data.frame()
for (i in x_var) {
  resulT3<-subset(result,outcome==i)
  p<-as.numeric(as.character(resulT3$pvalue))
  fdr=fdrtool(p,statistic="pvalue",plot = F)
  qval<-data.frame(fdr$qval)
  p2<-p.adjust(p,method = "fdr",length(p))
  p2<-data.frame(p2)
  resulT3<-cbind(resulT3,qval,p2)
  result_new<-rbind(result_new,resulT3)
  
}

group<-rep("Adjusted",times=nrow(result_new))

resulT3<-data.frame(result_new,group)

result_full<-rbind(result1,resulT3)
write.csv(result_full,file="Result/LME/association/linear_T3.csv",row.names = F) 



#3.2 binary ------------------------------------------------------------------


# unadjusted --------------------------------------------------------------


x_var<-c("Hyperlipid_32","TG_32_High","CHOL_32_High")
result<-data.frame()
for (i in x_var){
  for(j in y_var){
    # xnam_i1<-paste(i,j,sep="*")
    xnam_i<-paste(j,"pre_BMI","Parity","Age","gesweek_32",sep="+")
    fit_i<-glm(as.formula(paste(i,"~",j)),family = "binomial",data=data_all)
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

names(result)<-c("outcome","genus","OR","low","upper","pvalueֵ")


result_new<-data.frame()
for (i in x_var) {
  resulT3<-subset(result,outcome==i)
  p<-as.numeric(as.character(resulT3$pvalue))
  fdr=fdrtool(p,statistic="pvalue",plot = F)
  qval<-data.frame(fdr$qval)
  p2<-p.adjust(p,method = "fdr",length(p))
  p2<-data.frame(p2)
  resulT3<-cbind(resulT3,qval,p2)
  result_new<-rbind(result_new,resulT3)
  
}
group<-rep("Unadjusted",times=nrow(result_new))

result1<-data.frame(result_new,group)



# adjusted ----------------------------------------------------------------

x_var<-c("Hyperlipid_32","TG_32_High","CHOL_32_High")
result<-data.frame()
for (i in x_var){
  for(j in y_var){
    # xnam_i1<-paste(i,j,sep="*")
    xnam_i<-paste(j,"pre_BMI","Parity","Age","gesweek_32",sep="+")
    fit_i<-glm(as.formula(paste(i,"~",xnam_i)),family = "binomial",data=data_all)
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

names(result)<-c("outcome","genus","OR","low","upper","pvalueֵ")


result_new<-data.frame()
for (i in x_var) {
  resulT3<-subset(result,outcome==i)
  p<-as.numeric(as.character(resulT3$pvalue))
  fdr=fdrtool(p,statistic="pvalue",plot = F)
  qval<-data.frame(fdr$qval)
  p2<-p.adjust(p,method = "fdr",length(p))
  p2<-data.frame(p2)
  resulT3<-cbind(resulT3,qval,p2)
  result_new<-rbind(result_new,resulT3)
  
}

group<-rep("Adjusted",times=nrow(result_new))

resulT3<-data.frame(result_new,group)

result_full<-rbind(result1,resulT3)

write.csv(result_full,file="Result/LME/association/linear_T3_group.csv",row.names = F)







#4. T2 and TG_T3 ------------------------------------------------------------
data<-read.xlsx("Data/data_total.xlsx",sheet =1)

genus_filter<-as.data.frame(t(genus))
#colnames(genus_filter)<-as.character(genus$Genus)
#genus_filter<-genus_filter/32007*100

#write.csv(genus_filter,file = "Data/16S/genus_filter_RA.csv")

# for (i in 1:ncol(genus_filter)){
#   genus_filter[,i][genus_filter[,i]==0] <-1/32007*100
# }
genus_filter<-genus_filter+1
genus_filter<-genus_filter/(32007+98)*100

genus_filter_log<-log10(genus_filter)
genus_filter_log$ID<-rownames(genus_filter_log)

ID_24<-data$ID_24
genus_filter_log<-genus_filter_log[ID_24,]

data_all<-merge(genus_filter_log,data,by.x = "ID",by.y = "ID_24")

x_var<-c("TG_32","CHOL_32","delta_TG","delta_CHOL")
y_var<-genus_name

cov<-c("pre_BMI","Parity","Age","gesweek_seq")
dir.create("Result/LME/association")

#4.1 linear ------------------------------------------------------------------


# unadjusted --------------------------------------------------------------


###unadjusted
result<-data.frame()
for (i in x_var){
  for(j in y_var){
    # xnam_i1<-paste(i,j,sep="*")
    xnam_i<-paste(j,"pre_BMI","Parity","Age","gesweek_32",sep="+")
    fit_i<-lm(as.formula(paste(i,"~",j)),data=data_all)
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


names(result)<-c("outcome","genus","beta","lower","upper","pvalue")   

result_new<-data.frame()
for (i in x_var) {
  resulT3<-subset(result,outcome==i)
  p<-as.numeric(as.character(resulT3$pvalue))
  fdr=fdrtool(p,statistic="pvalue",plot = F)
  qval<-data.frame(fdr$qval)
  p2<-p.adjust(p,method = "fdr",length(p))
  p2<-data.frame(p2)
  resulT3<-cbind(resulT3,qval,p2)
  result_new<-rbind(result_new,resulT3)
  
}
group<-rep("Unadjusted",times=nrow(result_new))

result1<-data.frame(result_new,group)


# adjusted ----------------------------------------------------------------



###Adjusted
result<-data.frame()
for (i in x_var){
  for(j in y_var){
    # xnam_i1<-paste(i,j,sep="*")
    xnam_i<-paste(j,"pre_BMI","Parity","Age","gesweek_32",sep="+")
    fit_i<-lm(as.formula(paste(i,"~",xnam_i)),data=data_all)
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


names(result)<-c("outcome","genus","beta","lower","upper","pvalue")   

result_new<-data.frame()
for (i in x_var) {
  resulT3<-subset(result,outcome==i)
  p<-as.numeric(as.character(resulT3$pvalue))
  fdr=fdrtool(p,statistic="pvalue",plot = F)
  qval<-data.frame(fdr$qval)
  p2<-p.adjust(p,method = "fdr",length(p))
  p2<-data.frame(p2)
  resulT3<-cbind(resulT3,qval,p2)
  result_new<-rbind(result_new,resulT3)
  
}

group<-rep("Adjusted",times=nrow(result_new))

resulT3<-data.frame(result_new,group)

result_full<-rbind(result1,resulT3)
write.csv(result_full,file="Result/LME/association/linear_T2micro_T3.csv",row.names = F) 



#4.2 binary ------------------------------------------------------------------


# unadjusted --------------------------------------------------------------


x_var<-c("Hyperlipid_32","TG_32_High","CHOL_32_High")
result<-data.frame()
for (i in x_var){
  for(j in y_var){
    # xnam_i1<-paste(i,j,sep="*")
    xnam_i<-paste(j,"pre_BMI","Parity","Age","gesweek_32",sep="+")
    fit_i<-glm(as.formula(paste(i,"~",j)),family = "binomial",data=data_all)
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

names(result)<-c("outcome","genus","OR","low","upper","pvalueֵ")


result_new<-data.frame()
for (i in x_var) {
  resulT3<-subset(result,outcome==i)
  p<-as.numeric(as.character(resulT3$pvalue))
  fdr=fdrtool(p,statistic="pvalue",plot = F)
  qval<-data.frame(fdr$qval)
  p2<-p.adjust(p,method = "fdr",length(p))
  p2<-data.frame(p2)
  resulT3<-cbind(resulT3,qval,p2)
  result_new<-rbind(result_new,resulT3)
  
}
group<-rep("Unadjusted",times=nrow(result_new))

result1<-data.frame(result_new,group)



# adjusted ----------------------------------------------------------------

x_var<-c("Hyperlipid_32","TG_32_High","CHOL_32_High")
result<-data.frame()
for (i in x_var){
  for(j in y_var){
    # xnam_i1<-paste(i,j,sep="*")
    xnam_i<-paste(j,"pre_BMI","Parity","Age","gesweek_32",sep="+")
    fit_i<-glm(as.formula(paste(i,"~",xnam_i)),family = "binomial",data=data_all)
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

names(result)<-c("outcome","genus","OR","low","upper","pvalueֵ")


result_new<-data.frame()
for (i in x_var) {
  resulT3<-subset(result,outcome==i)
  p<-as.numeric(as.character(resulT3$pvalue))
  fdr=fdrtool(p,statistic="pvalue",plot = F)
  qval<-data.frame(fdr$qval)
  p2<-p.adjust(p,method = "fdr",length(p))
  p2<-data.frame(p2)
  resulT3<-cbind(resulT3,qval,p2)
  result_new<-rbind(result_new,resulT3)
  
}

group<-rep("Adjusted",times=nrow(result_new))

resulT3<-data.frame(result_new,group)

result_full<-rbind(result1,resulT3)

write.csv(result_full,file="Result/LME/association/linear_T2micro_T3_group.csv",row.names = F)




#4.3 multiple category -------------------------------------------------------

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

data_all$dyslipimedia_32<-ordered(data_all$dyslipimedia_32,levels=c("Control", "Hypercholesterolemia",
                                                                    "Hyperlipidemia", "Hypertriglyceridemia"))
data_all$delta_TG_group1<-ordered(data_all$delta_TG_group1,levels=c("Group1","Group2","Group3","Group4"))
data_all$delta_CHOL_group1<-ordered(data_all$delta_TG_group1,levels=c("Group1","Group2","Group3","Group4"))
data_all$delta_TG_group2<-ordered(data_all$delta_TG_group1,levels=c("Group1","Group2","Group3"))
data_all$delta_CHOL_group2<-ordered(data_all$delta_TG_group1,levels=c("Group1","Group2","Group3"))



#4.3.1 four category ------------------------------------------------------------


# unadjusted --------------------------------------------------------------



x_var<-c("dyslipimedia_32","delta_TG_group1","delta_CHOL_group1")
result<-data.frame()
for (i in x_var){
  for(j in y_var){
    # xnam_i1<-paste(i,j,sep="*")
    #xnam_i<-paste(j,"pre_BMI","Parity","Age","gesweek_32",sep="+")
    test<-multinom(as.formula(paste(i,"~",j,sep = "")),data=data_all)
    r<-coef(test)
    z <- summary(test)$coefficients/summary(test)$standard.errors
    p_i <- (1 - pnorm(abs(z), 0, 1)) * 2
    result_i<-c(r[,2],z[,2])
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
    xnam_i<-paste(j,"pre_BMI","Parity","Age","gesweek_32",sep="+")
    test<-multinom(as.formula(paste(i,"~",xnam_i,sep = "")),data=data_all)
    r<-coef(test)
    z <- summary(test)$coefficients/summary(test)$standard.errors
    p_i <- (1 - pnorm(abs(z), 0, 1)) * 2
    result_i<-c(r[,2],z[,2])
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

write.csv(result_full,file="Result/LME/association/linear_T2micro_T3_multiple_group1.csv",row.names = F)



#4.3.2 three category ------------------------------------------------------------


# unadjusted --------------------------------------------------------------



x_var<-c("delta_TG_group2","delta_CHOL_group2")
result<-data.frame()
for (i in x_var){
  for(j in y_var){
    # xnam_i1<-paste(i,j,sep="*")
    #xnam_i<-paste(j,"pre_BMI","Parity","Age","gesweek_32",sep="+")
    test<-multinom(as.formula(paste(i,"~",j,sep = "")),data=data_all)
    r<-coef(test)
    z <- summary(test)$coefficients/summary(test)$standard.errors
    p_i <- (1 - pnorm(abs(z), 0, 1)) * 2
    result_i<-c(r[,2],z[,2])
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
    xnam_i<-paste(j,"pre_BMI","Parity","Age","gesweek_32",sep="+")
    test<-multinom(as.formula(paste(i,"~",xnam_i,sep = "")),data=data_all)
    r<-coef(test)
    z <- summary(test)$coefficients/summary(test)$standard.errors
    p_i <- (1 - pnorm(abs(z), 0, 1)) * 2
    result_i<-c(r[,2],z[,2])
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

write.csv(result_full,file="Result/LME/association/linear_T2micro_T3_multiple_group2.csv",row.names = F)



# Part 5 plot_reg ---------------------------------------------------------


setwd(source_dir)
# plot  -------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(ggsci)
library(openxlsx)
dir.create("Result/LME/errorbar")

genus_select<-read.xlsx('Result/LME/result_combine.xlsx',sheet = 3)
genus_select<-genus_select$unique_genus


# LME ---------------------------------------------------------------------

data_plot<-read.csv("Result/LME/lmer_lipid_log10_new.csv",header = T,sep = ",")

data_plot2<-subset(data_plot,variable2%in%genus_select)

data_plot_unadjusted<-subset(data_plot2,group=="Unadjusted")
data_plot_adjusted<-subset(data_plot2,group=="Adjusted")

beta2<-data_plot_unadjusted$beta
data_plot2<-data.frame(data_plot_adjusted,beta2)

data_plot_1<-subset(data_plot2,variable1=="TG")

data_plot_1<-arrange(data_plot_1, desc(beta))

name<-as.character(unique(data_plot_1$variable2))

ggplot(data_plot_1,aes(x=variable2,y=beta))+
  geom_errorbar(aes(ymin=lower,ymax=upper),position=position_dodge(width=0.8),width=0.1)+
  geom_point(size=3.5,position=position_dodge(width=0.8),color="#00afbb")+
  geom_point(aes(x=variable2,y=beta2),size=3.5,shape=17,color="#e7b800")+
  # geom_point(size=3.5,shape=18,position=position_dodge(width=0.8))+
  geom_hline(yintercept=0,color="red",linetype="dashed")+
  labs(x="",y="β (95%CI)")+
  theme_bw()+
  theme(axis.text.y = element_text(face = "italic"))+
  scale_color_npg()+
  scale_x_discrete(limits=name)+
  # facet_wrap(~variable1,nrow = 1,scales = "free")+
  coord_flip()
ggsave(filename = "Result/LME/errorbar/plot_LME_linear_TG.pdf",width = 8,height = 8)





data_plot_1<-subset(data_plot2,variable1=="CHOL")
data_plot_1<-arrange(data_plot_1, desc(beta))

name2<-as.character(unique(data_plot_1$variable2))


ggplot(data_plot_1,aes(x=variable2,y=beta))+
  geom_errorbar(aes(ymin=lower,ymax=upper),position=position_dodge(width=0.8),width=0.1)+
  geom_point(size=3.5,position=position_dodge(width=0.8),color="#00afbb")+
  geom_point(aes(y=beta2),size=3.5,shape=17,color="#e7b800")+
  geom_hline(yintercept=0,color="red",linetype="dashed")+
  labs(x="",y="β (95%CI)")+
  theme_bw()+
  theme(axis.text.y = element_text(face = "italic"))+
  scale_color_npg()+
  scale_x_discrete(limits=name2)+
  # facet_wrap(~variable1,nrow = 1,scales = "free")+
  coord_flip()
ggsave(filename = "Result/LME/errorbar/plot_LME_linear_CHOL.pdf",width = 8,height = 8)



# linear ------------------------------------------------------------------


file_name1<-c('linear_T2')
file_name2<-c('linear_T3','linear_T2micro_T3')



for (i in file_name1) {
  data_plot<-read.csv(paste0("Result/LME/association/",i,".csv"),header = T,sep = ",")
  data_plot2<-subset(data_plot,genus%in%genus_select)
  
  data_plot_unadjusted<-subset(data_plot2,group=="Unadjusted")
  data_plot_adjusted<-subset(data_plot2,group=="Adjusted")
  
  beta2<-data_plot_unadjusted$beta
  data_plot2<-data.frame(data_plot_adjusted,beta2)
  
  data_plot_1<-subset(data_plot2,outcome=="TG_24")
  # data_plot_1<-arrange(data_plot_1, desc(beta))
  # name<-as.character(unique(data_plot_1$genus))
  
  ggplot(data_plot_1,aes(x=genus,y=beta))+
    geom_errorbar(aes(ymin=lower,ymax=upper),position=position_dodge(width=0.8),width=0.1)+
    geom_point(size=3.5,position=position_dodge(width=0.8),color="#00afbb")+
    geom_point(aes(y=beta2),size=3.5,shape=17,color="#e7b800")+
    # geom_point(size=3.5,shape=18,position=position_dodge(width=0.8))+
    geom_hline(yintercept=0,color="red",linetype="dashed")+
    labs(x="",y="β (95%CI)")+
    theme_bw()+
    scale_color_npg()+
    scale_x_discrete(limits=name)+
    theme(axis.text.y = element_text(face = "italic"))+
    # facet_wrap(~outcome,nrow = 1,scales = "free")+
    coord_flip()
  ggsave(filename = paste0("Result/LME/errorbar/",i,"_TG.pdf"),width = 8,height = 8)
  
  data_plot_1<-subset(data_plot2,outcome=="CHOL_24")
  # data_plot_1<-arrange(data_plot_1, desc(beta))
  # name2<-as.character(unique(data_plot_1$genus))
  
  ggplot(data_plot_1,aes(x=genus,y=beta))+
    geom_errorbar(aes(ymin=lower,ymax=upper),position=position_dodge(width=0.8),width=0.1)+
    geom_point(size=3.5,position=position_dodge(width=0.8),color="#00afbb")+
    geom_point(aes(y=beta2),size=3.5,shape=17,color="#e7b800")+
    # geom_point(size=3.5,shape=18,position=position_dodge(width=0.8))+
    geom_hline(yintercept=0,color="red",linetype="dashed")+
    labs(x="",y="β (95%CI)")+
    theme_bw()+
    scale_color_npg()+
    scale_x_discrete(limits=name2)+
    theme(axis.text.y = element_text(face = "italic"))+
    # facet_wrap(~outcome,nrow = 1,scales = "free")+
    coord_flip()
  ggsave(filename = paste0("Result/LME/errorbar/",i,"_CHOL.pdf"),width = 8,height = 8)
  
}


for (i in file_name2) {
  data_plot<-read.csv(paste0("Result/LME/association/",i,".csv"),header = T,sep = ",")
  data_plot2<-subset(data_plot,genus%in%genus_select)
  
  data_plot_unadjusted<-subset(data_plot2,group=="Unadjusted")
  data_plot_adjusted<-subset(data_plot2,group=="Adjusted")
  
  beta2<-data_plot_unadjusted$beta
  data_plot2<-data.frame(data_plot_adjusted,beta2)
  
  data_plot_1<-subset(data_plot2,outcome=="TG_32")
  # data_plot_1<-arrange(data_plot_1, desc(beta))
  # name<-as.character(unique(data_plot_1$genus))
  
  ggplot(data_plot_1,aes(x=genus,y=beta))+
    geom_errorbar(aes(ymin=lower,ymax=upper),position=position_dodge(width=0.8),width=0.1)+
    geom_point(size=3.5,position=position_dodge(width=0.8),color="#00afbb")+
    geom_point(aes(y=beta2),size=3.5,shape=17,color="#e7b800")+
    # geom_point(size=3.5,shape=18,position=position_dodge(width=0.8))+
    geom_hline(yintercept=0,color="red",linetype="dashed")+
    labs(x="",y="β (95%CI)")+
    theme_bw()+
    scale_color_npg()+
    scale_x_discrete(limits=name)+
    theme(axis.text.y = element_text(face = "italic"))+
    # facet_wrap(~outcome,nrow = 1,scales = "free")+
    coord_flip()
  ggsave(filename = paste0("Result/LME/errorbar/",i,"_TG.pdf"),width = 8,height = 8)
  
  data_plot_1<-subset(data_plot2,outcome=="CHOL_32")
  # data_plot_1<-arrange(data_plot_1, desc(beta))
  # name2<-as.character(unique(data_plot_1$genus))
  
  ggplot(data_plot_1,aes(x=genus,y=beta))+
    geom_errorbar(aes(ymin=lower,ymax=upper),position=position_dodge(width=0.8),width=0.1)+
    geom_point(size=3.5,position=position_dodge(width=0.8),color="#00afbb")+
    geom_point(aes(y=beta2),size=3.5,shape=17,color="#e7b800")+
    # geom_point(size=3.5,shape=18,position=position_dodge(width=0.8))+
    geom_hline(yintercept=0,color="red",linetype="dashed")+
    labs(x="",y="β (95%CI)")+
    theme_bw()+
    scale_color_npg()+
    scale_x_discrete(limits=name2)+
    theme(axis.text.y = element_text(face = "italic"))+
    # facet_wrap(~outcome,nrow = 1,scales = "free")+
    coord_flip()
  ggsave(filename = paste0("Result/LME/errorbar/",i,"_CHOL.pdf"),width = 8,height = 8)
  
}

# binomial ------------------------------------------------------------------


file_name1<-c('linear_T2_group')

file_name2<-c('linear_T3_group','linear_T2micro_T3_group')


for (i in file_name1) {
  data_plot<-read.csv(paste0("Result/LME/association/",i,".csv"),header = T,sep = ",")
  data_plot2<-subset(data_plot,genus%in%genus_select)
  
  data_plot_unadjusted<-subset(data_plot2,group=="Unadjusted")
  data_plot_adjusted<-subset(data_plot2,group=="Adjusted")
  
  OR2<-data_plot_unadjusted$OR
  data_plot2<-data.frame(data_plot_adjusted,OR2)
  
  data_plot_1<-subset(data_plot2,outcome=="CHOL_24_High")
  # data_plot_1<-arrange(data_plot_1, desc(OR))
  # name2<-as.character(unique(data_plot_1$genus))
  
  ggplot(data_plot_1,aes(x=genus,y=OR))+
    geom_errorbar(aes(ymin=low,ymax=upper),position=position_dodge(width=0.8),width=0.1)+
    geom_point(size=3.5,position=position_dodge(width=0.8),color="#00afbb")+
    geom_point(aes(y=OR2),size=3.5,shape=17,color="#e7b800")+
    # geom_point(size=3.5,shape=18,position=position_dodge(width=0.8))+
    geom_hline(yintercept=1,color="red",linetype="dashed")+
    labs(x="",y="OR (95%CI)")+
    theme_bw()+
    scale_color_npg()+
    scale_x_discrete(limits=name2)+
    theme(axis.text.y = element_text(face = "italic"))+
    # facet_wrap(~outcome,nrow = 1,scales = "free")+
    coord_flip()
  ggsave(filename = paste0("Result/LME/errorbar/",i,"_CHOL_High.pdf"),width = 8,height = 8)
  
  
  data_plot_1<-subset(data_plot2,outcome=="TG_24_High")
  # data_plot_1<-arrange(data_plot_1, desc(OR))
  # name2<-as.character(unique(data_plot_1$genus))
  
  ggplot(data_plot_1,aes(x=genus,y=OR))+
    geom_errorbar(aes(ymin=low,ymax=upper),position=position_dodge(width=0.8),width=0.1)+
    geom_point(size=3.5,position=position_dodge(width=0.8),color="#00afbb")+
    geom_point(aes(y=OR2),size=3.5,shape=17,color="#e7b800")+
    # geom_point(size=3.5,shape=18,position=position_dodge(width=0.8))+
    geom_hline(yintercept=1,color="red",linetype="dashed")+
    labs(x="",y="OR (95%CI)")+
    theme_bw()+
    scale_color_npg()+
    scale_x_discrete(limits=name)+
    theme(axis.text.y = element_text(face = "italic"))+
    # facet_wrap(~outcome,nrow = 1,scales = "free")+
    coord_flip()
  ggsave(filename = paste0("Result/LME/errorbar/",i,"_TG_High.pdf"),width = 8,height = 8)
  
  
  data_plot_1<-subset(data_plot2,outcome=="Hyperlipid_24")
  # data_plot_1<-arrange(data_plot_1, desc(OR))
  # name2<-as.character(unique(data_plot_1$genus))
  
  ggplot(data_plot_1,aes(x=genus,y=OR))+
    geom_errorbar(aes(ymin=low,ymax=upper),position=position_dodge(width=0.8),width=0.1)+
    geom_point(size=3.5,position=position_dodge(width=0.8),color="#00afbb")+
    geom_point(aes(y=OR2),size=3.5,shape=17,color="#e7b800")+
    # geom_point(size=3.5,shape=18,position=position_dodge(width=0.8))+
    geom_hline(yintercept=1,color="red",linetype="dashed")+
    labs(x="",y="OR (95%CI)")+
    theme_bw()+
    scale_color_npg()+
    scale_x_discrete(limits=name2)+
    theme(axis.text.y = element_text(face = "italic"))+
    # facet_wrap(~outcome,nrow = 1,scales = "free")+
    coord_flip()
  ggsave(filename = paste0("Result/LME/errorbar/",i,"Hyperlipid_24.pdf"),width = 8,height = 8)
}



for (i in file_name2) {
  data_plot<-read.csv(paste0("Result/LME/association/",i,".csv"),header = T,sep = ",")
  data_plot2<-subset(data_plot,genus%in%genus_select)
  
  data_plot_unadjusted<-subset(data_plot2,group=="Unadjusted")
  data_plot_adjusted<-subset(data_plot2,group=="Adjusted")
  
  OR2<-data_plot_unadjusted$OR
  data_plot2<-data.frame(data_plot_adjusted,OR2)
  
  data_plot_1<-subset(data_plot2,outcome=="CHOL_32_High")
  # data_plot_1<-arrange(data_plot_1, desc(OR))
  # name2<-as.character(unique(data_plot_1$genus))
  
  ggplot(data_plot_1,aes(x=genus,y=OR))+
    geom_errorbar(aes(ymin=low,ymax=upper),position=position_dodge(width=0.8),width=0.1)+
    geom_point(size=3.5,position=position_dodge(width=0.8),color="#00afbb")+
    geom_point(aes(y=OR2),size=3.5,shape=17,color="#e7b800")+
    # geom_point(size=3.5,shape=18,position=position_dodge(width=0.8))+
    geom_hline(yintercept=1,color="red",linetype="dashed")+
    labs(x="",y="OR (95%CI)")+
    theme_bw()+
    scale_color_npg()+
    scale_x_discrete(limits=name2)+
    theme(axis.text.y = element_text(face = "italic"))+
    # facet_wrap(~outcome,nrow = 1,scales = "free")+
    coord_flip()
  ggsave(filename = paste0("Result/LME/errorbar/",i,"_CHOL_High.pdf"),width = 8,height = 8)
  
  
  data_plot_1<-subset(data_plot2,outcome=="TG_32_High")
  # data_plot_1<-arrange(data_plot_1, desc(OR))
  # name2<-as.character(unique(data_plot_1$genus))
  
  ggplot(data_plot_1,aes(x=genus,y=OR))+
    geom_errorbar(aes(ymin=low,ymax=upper),position=position_dodge(width=0.8),width=0.1)+
    geom_point(size=3.5,position=position_dodge(width=0.8),color="#00afbb")+
    geom_point(aes(y=OR2),size=3.5,shape=17,color="#e7b800")+
    # geom_point(size=3.5,shape=18,position=position_dodge(width=0.8))+
    geom_hline(yintercept=1,color="red",linetype="dashed")+
    labs(x="",y="OR (95%CI)")+
    theme_bw()+
    scale_color_npg()+
    scale_x_discrete(limits=name)+
    theme(axis.text.y = element_text(face = "italic"))+
    # facet_wrap(~outcome,nrow = 1,scales = "free")+
    coord_flip()
  ggsave(filename = paste0("Result/LME/errorbar/",i,"_TG_High.pdf"),width = 8,height = 8)
  
  
  data_plot_1<-subset(data_plot2,outcome=="Hyperlipid_32")
  # data_plot_1<-arrange(data_plot_1, desc(OR))
  # name2<-as.character(unique(data_plot_1$genus))
  
  ggplot(data_plot_1,aes(x=genus,y=OR))+
    geom_errorbar(aes(ymin=low,ymax=upper),position=position_dodge(width=0.8),width=0.1)+
    geom_point(size=3.5,position=position_dodge(width=0.8),color="#00afbb")+
    geom_point(aes(y=OR2),size=3.5,shape=17,color="#e7b800")+
    # geom_point(size=3.5,shape=18,position=position_dodge(width=0.8))+
    geom_hline(yintercept=1,color="red",linetype="dashed")+
    labs(x="",y="OR (95%CI)")+
    theme_bw()+
    scale_color_npg()+
    scale_x_discrete(limits=name2)+
    theme(axis.text.y = element_text(face = "italic"))+
    # facet_wrap(~outcome,nrow = 1,scales = "free")+
    coord_flip()
  ggsave(filename = paste0("Result/LME/errorbar/",i,"Hyperlipid_32.pdf"),width = 8,height = 8)
}




# boxplot-1 -----------------------------------------------------------------

#setwd("E:/study2/microbiome/study_new/total")
library(Maaslin2)
library(ggplot2)
library(lme4)
library(nlme)
library(fdrtool)
library(openxlsx)
library(ggpubr)
library(ggsci)
genus<-read.xlsx("Data/16S/genus_filter.xlsx",sheet = 2)
rownames(genus)<-genus$Genus
genus_name<-genus$Genus
genus<-genus[,-c(1:6)]
data<-read.xlsx("Data/data_total.xlsx",sheet =2)
genus_filter<-as.data.frame(t(genus))

genus_filter<-genus_filter+1
genus_filter<-genus_filter/(32007+98)*100

genus_filter_log<-log10(genus_filter)
genus_filter_log$ID<-rownames(genus_filter_log)
data_all<-merge(genus_filter_log,data,by.x = "ID",by.y = "ID")

x_var<-c("TG","CHOL")
y_var<-genus_name

cov<-c("pre_BMI","Parity","Age","gesweek_seq")

dir.create("Result/LME/boxplot")

genus_select<-read.xlsx('Result/LME/result_combine.xlsx',sheet = 3)
genus_select<-genus_select$unique_genus




# dyslipidemia ------------------------------------------------------------
my_comparisons <- list( c("Control", "Hypercholesterolemia"), c("Control", "Hypertriglyceridemia"),
                        c("Control", "Hyperlipidemia") )

data_all$dyslipimedia<-ordered(data_all$dyslipimedia,levels=c("Control", "Hypercholesterolemia","Hypertriglyceridemia","Hyperlipidemia"))
dir.create("Result/LME/boxplot/dyslipidemia")
def_boxplot_plot<-function(i){
  ggboxplot(data_all,x="dyslipimedia",y=i,#add="jitter",
            fill = "dyslipimedia",notch = FALSE)+
    stat_compare_means(comparisons = my_comparisons)+
    stat_compare_means()+
    # stat_compare_means(label = "p.signif", method = "wilcox.test",
    #                    ref.group = "Control") +
    #scale_y_continuous(limits=c(25, 180))+
    theme(legend.position = "none")+
    scale_fill_npg()+
    theme_bw()+
    theme(axis.text.y = element_text(face = "italic"))+
    theme(axis.text.x=element_text(angle=45,size=8,vjust = 1,hjust = 1))+
    labs(x="",y=i)#+
  # facet_wrap(~Time,nrow = 1)
  ggsave(filename = paste0("Result/LME/boxplot/dyslipidemia/",i,".pdf"),width=6,height=5)
  
}

for (i in genus_select) {
  def_boxplot_plot(i)
}



def_boxplot_plot<-function(i){
  ggboxplot(data_all,x="dyslipimedia",y=i,#add="jitter",
            fill = "dyslipimedia",notch = FALSE)+
    stat_compare_means(comparisons = my_comparisons)+
    stat_compare_means()+
    # stat_compare_means(label = "p.signif", method = "wilcox.test",
    #                    ref.group = "Control") +
    #scale_y_continuous(limits=c(25, 180))+
    theme(legend.position = "none")+
    scale_fill_npg()+
    theme_bw()+
    theme(axis.text.x=element_text(angle=45,size=8,vjust = 1,hjust = 1))+
    labs(x="",y=i)+
    theme(axis.text.y = element_text(face = "italic"))+
    facet_wrap(~Time,nrow = 1)
  ggsave(filename = paste0("Result/LME/boxplot/dyslipidemia/",i,"_by_time.pdf"),width=12,height=5)
  
}

for (i in genus_select) {
  def_boxplot_plot(i)
}

# TG_high ------------------------------------------------------------
# my_comparisons <- list( c("Control", "Hypercholesterolemia"), c("Control", "Hypertriglyceridemia"),
#                         c("Control", "Hyperlipidemia") )

data_all$TG_High<-ordered(data_all$TG_High,levels=c(0,1),labels=c("Control", "Hypertriglyceridemia"))
dir.create("Result/LME/boxplot/TG_High")
def_boxplot_plot<-function(i){
  ggboxplot(data_all,x="TG_High",y=i,#add="jitter",
            fill = "TG_High",notch = FALSE)+
    #  stat_compare_means(comparisons = my_comparisons)+
    stat_compare_means()+
    # stat_compare_means(label = "p.signif", method = "wilcox.test",
    #                    ref.group = "Control") +
    #scale_y_continuous(limits=c(25, 180))+
    theme(legend.position = "none")+
    # scale_fill_npg()+
    scale_fill_manual(values = c("#4DBBD5FF", "#E64B35FF"))+
    theme_bw()+
    theme(axis.text.y = element_text(face = "italic"))+
    theme(axis.text.x=element_text(angle=45,size=8,vjust = 1,hjust = 1))+
    labs(x="",y=i)#+
  # facet_wrap(~Time,nrow = 1)
  ggsave(filename = paste0("Result/LME/boxplot/TG_High/",i,".pdf"),width=6,height=5)
  
}

for (i in genus_select) {
  def_boxplot_plot(i)
}



def_boxplot_plot<-function(i){
  ggboxplot(data_all,x="TG_High",y=i,#add="jitter",
            fill = "TG_High",notch = FALSE)+
    #  stat_compare_means(comparisons = my_comparisons)+
    stat_compare_means()+
    # stat_compare_means(label = "p.signif", method = "wilcox.test",
    #                    ref.group = "Control") +
    #scale_y_continuous(limits=c(25, 180))+
    theme(legend.position = "none")+
    #scale_fill_npg()+
    scale_fill_manual(values = c("#4DBBD5FF", "#E64B35FF"))+
    theme_bw()+
    theme(axis.text.y = element_text(face = "italic"))+
    theme(axis.text.x=element_text(angle=45,size=8,vjust = 1,hjust = 1))+
    labs(x="",y=i)+
    facet_wrap(~Time,nrow = 1)
  ggsave(filename = paste0("Result/LME/boxplot/TG_High/",i,"_by_Time.pdf"),width=8,height=5)
  
}

for (i in genus_select) {
  def_boxplot_plot(i)
}


# CHOL_High ------------------------------------------------------------
# my_comparisons <- list( c("Control", "Hypercholesterolemia"), c("Control", "Hypertriglyceridemia"),
#                         c("Control", "Hyperlipidemia") )

data_all$CHOL_High<-ordered(data_all$CHOL_High,levels=c(0,1),labels=c("Control", "Hypercholesterolemia"))
dir.create("Result/LME/boxplot/CHOL_High")
def_boxplot_plot<-function(i){
  ggboxplot(data_all,x="CHOL_High",y=i,#add="jitter",
            fill = "CHOL_High",notch = FALSE)+
    #  stat_compare_means(comparisons = my_comparisons)+
    stat_compare_means()+
    # stat_compare_means(label = "p.signif", method = "wilcox.test",
    #                    ref.group = "Control") +
    #scale_y_continuous(limits=c(25, 180))+
    theme(legend.position = "none")+
    #  scale_fill_npg()+
    scale_fill_manual(values = c("#4DBBD5FF", "#E64B35FF"))+
    theme_bw()+
    theme(axis.text.y = element_text(face = "italic"))+
    theme(axis.text.x=element_text(angle=45,size=8,vjust = 1,hjust = 1))+
    labs(x="",y=i)#+
  # facet_wrap(~Time,nrow = 1)
  ggsave(filename = paste0("Result/LME/boxplot/CHOL_High/",i,".pdf"),width=6,height=5)
  
}

for (i in genus_select) {
  def_boxplot_plot(i)
}



def_boxplot_plot<-function(i){
  ggboxplot(data_all,x="CHOL_High",y=i,#add="jitter",
            fill = "CHOL_High",notch = FALSE)+
    #  stat_compare_means(comparisons = my_comparisons)+
    stat_compare_means()+
    # stat_compare_means(label = "p.signif", method = "wilcox.test",
    #                    ref.group = "Control") +
    #scale_y_continuous(limits=c(25, 180))+
    theme(legend.position = "none")+
    #scale_fill_npg()+
    scale_fill_manual(values = c("#4DBBD5FF", "#E64B35FF"))+
    # scale_color_manual(values = c("#4DBBD5FF", "#E64B35FF"))+
    theme_bw()+
    theme(axis.text.y = element_text(face = "italic"))+
    theme(axis.text.x=element_text(angle=45,size=8,vjust = 1,hjust = 1))+
    labs(x="",y=i)+
    facet_wrap(~Time,nrow = 1)
  ggsave(filename = paste0("Result/LME/boxplot/CHOL_High/",i,"_by_Time.pdf"),width=8,height=5)
  
}

for (i in genus_select) {
  def_boxplot_plot(i)
}


# boxplot-2 -----------------------------------------------------------------

#setwd("E:/study2/microbiome/study_new/total")
library(Maaslin2)
library(ggplot2)
library(lme4)
library(nlme)
library(fdrtool)
library(openxlsx)
library(ggpubr)
library(ggsci)
genus<-read.xlsx("Data/16S/genus_filter.xlsx",sheet = 2)
rownames(genus)<-genus$Genus
genus_name<-genus$Genus
genus<-genus[,-c(1:6)]
data<-read.xlsx("Data/data_total.xlsx",sheet =2)
genus_filter<-as.data.frame(t(genus))

genus_filter<-genus_filter+1
genus_filter<-genus_filter/(32007+98)*100

genus_filter_log<-log10(genus_filter)
genus_filter_log$ID<-rownames(genus_filter_log)
data_all<-merge(genus_filter_log,data,by.x = "ID",by.y = "ID")


data_24<-subset(data_all,Time=="24week")
data_32<-subset(data_all,Time=="32week")


dat_T2_T3<-data.frame(data_24[genus_name],data_32[c("dyslipimedia","TG_High","CHOL_High")])
x_var<-c("TG","CHOL")
y_var<-genus_name

cov<-c("pre_BMI","Parity","Age","gesweek_seq")

dir.create("Result/LME/boxplot2")

genus_select<-read.xlsx('Result/LME/result_combine.xlsx',sheet = 3)
genus_select<-genus_select$unique_genus




# dyslipidemia ------------------------------------------------------------
my_comparisons <- list( c("Control", "Hypercholesterolemia"), c("Control", "Hypertriglyceridemia"),
                        c("Control", "Hyperlipidemia") )

dat_T2_T3$dyslipimedia<-ordered(dat_T2_T3$dyslipimedia,levels=c("Control", "Hypercholesterolemia","Hypertriglyceridemia","Hyperlipidemia"))
dir.create("Result/LME/boxplot2/dyslipidemia")
def_boxplot_plot<-function(i){
  ggboxplot(dat_T2_T3,x="dyslipimedia",y=i,#add="jitter",
            fill = "dyslipimedia",notch = FALSE)+
    stat_compare_means(comparisons = my_comparisons)+
    stat_compare_means()+
    # stat_compare_means(label = "p.signif", method = "wilcox.test",
    #                    ref.group = "Control") +
    #scale_y_continuous(limits=c(25, 180))+
    theme(legend.position = "none")+
    scale_fill_npg()+
    theme_bw()+
    theme(axis.text.y = element_text(face = "italic"))+
    theme(axis.text.x=element_text(angle=45,size=8,vjust = 1,hjust = 1))+
    labs(x="",y=i)#+
  # facet_wrap(~Time,nrow = 1)
  ggsave(filename = paste0("Result/LME/boxplot2/dyslipidemia/",i,".pdf"),width=6,height=5)
  
}

for (i in genus_name) {
  def_boxplot_plot(i)
}


# 
# def_boxplot_plot<-function(i){
#   ggboxplot(dat_T2_T3,x="dyslipimedia",y=i,#add="jitter",
#             fill = "dyslipimedia",notch = FALSE)+
#     stat_compare_means(comparisons = my_comparisons)+
#     stat_compare_means()+
#     # stat_compare_means(label = "p.signif", method = "wilcox.test",
#     #                    ref.group = "Control") +
#     #scale_y_continuous(limits=c(25, 180))+
#     theme(legend.position = "none")+
#     scale_fill_npg()+
#     theme_bw()+
#     theme(axis.text.x=element_text(angle=45,size=8,vjust = 1,hjust = 1))+
#     labs(x="",y=i)+
#     facet_wrap(~Time,nrow = 1)
#   ggsave(filename = paste0("Result/LME/boxplot2/dyslipidemia/",i,"_by_time.pdf"),width=12,height=5)
#   
# }
# 
# for (i in genus_name) {
#   def_boxplot_plot(i)
# }

# TG_high ------------------------------------------------------------
# my_comparisons <- list( c("Control", "Hypercholesterolemia"), c("Control", "Hypertriglyceridemia"),
#                         c("Control", "Hyperlipidemia") )

dat_T2_T3$TG_High<-ordered(dat_T2_T3$TG_High,levels=c(0,1),labels=c("Control", "Hypertriglyceridemia"))
dir.create("Result/LME/boxplot2/TG_High")
def_boxplot_plot<-function(i){
  ggboxplot(dat_T2_T3,x="TG_High",y=i,#add="jitter",
            fill = "TG_High",notch = FALSE)+
    #  stat_compare_means(comparisons = my_comparisons)+
    stat_compare_means()+
    # stat_compare_means(label = "p.signif", method = "wilcox.test",
    #                    ref.group = "Control") +
    #scale_y_continuous(limits=c(25, 180))+
    theme(legend.position = "none")+
    # scale_fill_npg()+
    scale_fill_manual(values = c("#4DBBD5FF", "#E64B35FF"))+
    theme_bw()+
    theme(axis.text.y = element_text(face = "italic"))+
    theme(axis.text.x=element_text(angle=45,size=8,vjust = 1,hjust = 1))+
    labs(x="",y=i)#+
  # facet_wrap(~Time,nrow = 1)
  ggsave(filename = paste0("Result/LME/boxplot2/TG_High/",i,".pdf"),width=6,height=5)
  
}

for (i in genus_name) {
  def_boxplot_plot(i)
}

# 
# 
# def_boxplot_plot<-function(i){
#   ggboxplot(dat_T2_T3,x="TG_High",y=i,#add="jitter",
#             fill = "TG_High",notch = FALSE)+
#     #  stat_compare_means(comparisons = my_comparisons)+
#     stat_compare_means()+
#     # stat_compare_means(label = "p.signif", method = "wilcox.test",
#     #                    ref.group = "Control") +
#     #scale_y_continuous(limits=c(25, 180))+
#     theme(legend.position = "none")+
#     #scale_fill_npg()+
#     scale_fill_manual(values = c("#4DBBD5FF", "#E64B35FF"))+
#     theme_bw()+
#     theme(axis.text.x=element_text(angle=45,size=8,vjust = 1,hjust = 1))+
#     labs(x="",y=i)+
#     facet_wrap(~Time,nrow = 1)
#   ggsave(filename = paste0("Result/LME/boxplot2/TG_High/",i,"_by_Time.pdf"),width=8,height=5)
#   
# }
# 
# for (i in genus_name) {
#   def_boxplot_plot(i)
# }
# 

# CHOL_High ------------------------------------------------------------
# my_comparisons <- list( c("Control", "Hypercholesterolemia"), c("Control", "Hypertriglyceridemia"),
#                         c("Control", "Hyperlipidemia") )

dat_T2_T3$CHOL_High<-ordered(dat_T2_T3$CHOL_High,levels=c(0,1),
                             labels=c("Control", "Hypercholesterolemia"))
dir.create("Result/LME/boxplot2/CHOL_High")
def_boxplot_plot<-function(i){
  ggboxplot(dat_T2_T3,x="CHOL_High",y=i,#add="jitter",
            fill = "CHOL_High",notch = FALSE)+
    #  stat_compare_means(comparisons = my_comparisons)+
    stat_compare_means()+
    # stat_compare_means(label = "p.signif", method = "wilcox.test",
    #                    ref.group = "Control") +
    #scale_y_continuous(limits=c(25, 180))+
    theme(legend.position = "none")+
    #  scale_fill_npg()+
    scale_fill_manual(values = c("#4DBBD5FF", "#E64B35FF"))+
    theme_bw()+
    theme(axis.text.y = element_text(face = "italic"))+
    theme(axis.text.x=element_text(angle=45,size=8,vjust = 1,hjust = 1))+
    labs(x="",y=i)#+
  # facet_wrap(~Time,nrow = 1)
  ggsave(filename = paste0("Result/LME/boxplot2/CHOL_High/",i,".pdf"),width=6,height=5)
  
}

for (i in genus_name) {
  def_boxplot_plot(i)
}


# 
# def_boxplot_plot<-function(i){
#   ggboxplot(dat_T2_T3,x="CHOL_High",y=i,#add="jitter",
#             fill = "CHOL_High",notch = FALSE)+
#     #  stat_compare_means(comparisons = my_comparisons)+
#     stat_compare_means()+
#     # stat_compare_means(label = "p.signif", method = "wilcox.test",
#     #                    ref.group = "Control") +
#     #scale_y_continuous(limits=c(25, 180))+
#     theme(legend.position = "none")+
#     #scale_fill_npg()+
#     scale_fill_manual(values = c("#4DBBD5FF", "#E64B35FF"))+
#     # scale_color_manual(values = c("#4DBBD5FF", "#E64B35FF"))+
#     theme_bw()+
#     theme(axis.text.x=element_text(angle=45,size=8,vjust = 1,hjust = 1))+
#     labs(x="",y=i)+
#     facet_wrap(~Time,nrow = 1)
#   ggsave(filename = paste0("Result/LME/boxplot2/CHOL_High/",i,"_by_Time.pdf"),width=8,height=5)
#   
# }
# 
# for (i in genus_name) {
#   def_boxplot_plot(i)
# }

