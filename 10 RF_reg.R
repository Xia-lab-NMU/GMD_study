#setwd("D:/BaiduNetdiskDownload/amplicon/33MachineLearning/RF_classification")
library(openxlsx)
library(ggplot2)

setwd("E:/study2/microbiome/study_new/total/9Data Review")
metadata = read.xlsx("Data/data_total.xlsx", sheet=2,rowNames = F)
test_name<-c("TG","CHOL","ALP", "RBP")

#impute missing value
library(dplyr)
# data_rf_second<-data_rf_second %>% 
#   mutate_all(function(x){x[is.na(x)] <- mean(x)
#   x})

for (i in test_name) {
  
  metadata[,i][is.na(metadata[,i])==TRUE]<-mean(metadata[,i],na.rm = T)
}

#二、 genus_level -------------------------------------------------------------
source_dir<-"E:/study2/microbiome/study_new/total/9Data Review"
setwd(source_dir)
genus_filter<-read.csv("Data/16S/genus_filter_RA.csv",header = T,row.names = 1)
genus_name<-colnames(genus_filter)
# for (i in 1:ncol(genus_filter)){
#   genus_filter[,i][genus_filter[,i]==0] <-1/20993*100
# }
# 
# genus_filter_log<-log10(genus_filter)
genus_filter$ID<-rownames(genus_filter)
data_genus_RA<-merge(genus_filter,metadata,by.x = "ID",by.y = "ID")

#rownames(data_genus_RA)<-data_genus_RA$ID
data_rf_second<-subset(data_genus_RA,Time=="24week")
#data_rf_second<-subset(data_rf_second,CHOL_High==0 & TG_High==0)

data_rf_third<-subset(data_genus_RA,Time=="32week")
#data_rf_third<-subset(data_rf_third,number %in% data_rf_second$number)

CHOL_third<-data_rf_third$CHOL
TG_third<-data_rf_third$TG

data_rf_second<-data.frame(data_rf_second,CHOL_third,TG_third)
setwd("E:/study2/microbiome/study_new/total/9Data Review/Result/ML")
#随机森林分类
library(randomForest)
library(caret)
library(ggplot2)
library(bootstrap)
library(ROCR)
library(pROC)
set.seed(315)
# 
# data_rf_second$CHOL_third<-as.factor(data_rf_second$CHOL_third)
# data_rf_second$TG_third<-as.factor(data_rf_second$TG_third)
# 
# data_rf_second$hyperlipid_third<-as.factor(data_rf_second$hyperlipid_third)
# 
# #test_name<-colnames(data_rf_second[,72:91])
# genus_selected2<-c("Alistipes", "Bacteroides","Paraprevotella",
#                    "Christensenellaceae_R_7_group", "UCG_002","Clostridia_UCG_014")
# 
# test_name<-c("BUN", "CREA", "UA", "TP", "ALB", "CHOL", "TG", "ALT", "AST", 
#              "ALP", "LDH", "GGT", "TBIL", "DBIL", "TBA", "FMN", "FBG", "RBP", 
#              "cystain", "HCY")
# #impute missing value
# library(dplyr)
# # data_rf_second<-data_rf_second %>% 
# #   mutate_all(function(x){x[is.na(x)] <- mean(x)
# #   x})
# 
# for (i in test_name) {
#   
#   data_rf_second[,i][is.na(data_rf_second[,i])==TRUE]<-mean(data_rf_second[,i],na.rm = T)
# }


# 
# #1. genus alone ------------------------------------------------------------
genus_selected2<-c("Alistipes", "Bacteroides","Paraprevotella",
                   "Christensenellaceae_R_7_group", "UCG_002","Clostridia_UCG_014")



rf_classify<-function(i){
  # i<-"CHOL_third"
  dir.create(i)
  dir.create(paste0(i,"/genus_selected2"))
  varia<-c(i,genus_selected2)
  m<-100
  result<-data.frame()
  plot_data<-data.frame()
  pred<-data.frame()
  set.seed(123)
  accura<-data.frame()
  result_m<-data.frame()
  imp<-data.frame()
  result<-data.frame()
  b=1
  df<-data.frame()
  for (b in 1:m){
    
    Train <- createDataPartition(data_rf_second[,i], p=0.8, list=FALSE)
    training <- data_rf_second[ Train, varia]
    testing <- data_rf_second[ -Train, varia]
    
    rf = randomForest(as.formula(paste0(i,"~.") ), data=training,importance=TRUE, proximity=TRUE, ntree = 1000)
    
    pred_b=predict(rf,testing,type="response")
    df_b = data.frame(observed = testing[,i], predict = pred_b)
    times2<-rep(b,times=nrow(df_b))
    
    
    cor = cor.test(df_b[,1], df_b[,2], method = "spearman")
    df2 = df_b
    colnames(df2) = c("x", "y")
    m = lm(y ~ x, df2)
    rho<-round(cor$estimate, digits = 3)
    pvalue<-round(cor$p.value, digits = 3)
    R2<-round(summary(m)$r.squared, digits = 3)
    result_i<-c(rho,pvalue,R2)
    result<-rbind(result,result_i)
    
    df_b<-data.frame(df_b,times2)
    df<-rbind(df,df_b)
    # p = ggplot(df, aes(predict, observed)) +
    #   geom_point() +
    #   geom_smooth(method = "lm") +
    #   labs(title = paste("rho = " , round(cor$estimate, digits = 3), ", P = " , signif(cor$p.value, digits = 3), ", R2 = ", round(summary(m)$r.squared, digits = 3) , sep = "")) +
    #   theme_bw()
    # p
    # ggsave(paste("test_predict.pdf", sep=""), p, width = 4, height = 2.5)
    # 
    # imp= as.data.frame(rf$importance)
    # imp = imp[order(imp[,1],decreasing = T),]
    
    imp_i<-importance(rf)
    imp_i<-as.data.frame(imp_i)
    imp_i$variable<-rownames(imp_i)
    times2<-rep(b,times=nrow(imp_i))
    imp_i<-data.frame(imp_i,times2)
    imp<-rbind(imp,imp_i)
  } 
  names(result)<-c("rho","pvalue","R2")
  # names(result)<-c("times","sensitivity","specifity","auc")
  # names(plot_data)<-c("sensitivity","specificity","Group")
  # names(accura)<-c("Accuracy","Lower","higher","AUC","Time")
  # 
  # write.csv(accura,file=paste0(i,"/genus_selected2/accura.csv"),row.names = T)
  write.csv(imp,file=paste0(i,"/genus_selected2/imp.csv"),row.names = T)
  write.csv(df,file=paste0(i,"/genus_selected2/df.csv"),row.names = T)
  
  write.csv(result,file=paste0(i,"/genus_selected2/result.csv"),row.names = T)
  # write.csv(plot_data,file=paste0(i,"/genus_selected2/plot_data.csv"),row.names = T)
  # write.csv(pred,file=paste0(i,"/genus_selected2/pred.csv"),row.names = T)
}


for (i in c("CHOL_third","TG_third")) {
  rf_classify(i)
}


# 2.genus+test -------------------------------------------------------------
test_name<-c("TG","CHOL","ALP", "RBP")
rf_classify<-function(i){
  # i<-"CHOL_third"
  dir.create(i)
  dir.create(paste0(i,"/genus_selected2_test"))
  varia<-c(i,genus_selected2,test_name)
  m<-100
  result<-data.frame()
  plot_data<-data.frame()
  pred<-data.frame()
  set.seed(123)
  accura<-data.frame()
  result_m<-data.frame()
  imp<-data.frame()
  b=1
  df<-data.frame()
  for (b in 1:m){
    
    Train <- createDataPartition(data_rf_second[,i], p=0.8, list=FALSE)
    training <- data_rf_second[ Train, varia]
    testing <- data_rf_second[ -Train, varia]
    
    rf = randomForest(as.formula(paste0(i,"~.") ), data=training,importance=TRUE, proximity=TRUE, ntree = 1000)
    
    pred_b=predict(rf,testing,type="response")
    df_b = data.frame(observed = testing[,i], predict = pred_b)
    times2<-rep(b,times=nrow(df_b))
    
    
    cor = cor.test(df_b[,1], df_b[,2], method = "spearman")
    df2 = df_b
    colnames(df2) = c("x", "y")
    m = lm(y ~ x, df2)
    rho<-round(cor$estimate, digits = 3)
    pvalue<-round(cor$p.value, digits = 3)
    R2<-round(summary(m)$r.squared, digits = 3)
    result_i<-c(rho,pvalue,R2)
    result<-rbind(result,result_i)
    
    df_b<-data.frame(df_b,times2)
    df<-rbind(df,df_b)
    # p = ggplot(df, aes(predict, observed)) +
    #   geom_point() +
    #   geom_smooth(method = "lm") +
    #   labs(title = paste("rho = " , round(cor$estimate, digits = 3), ", P = " , signif(cor$p.value, digits = 3), ", R2 = ", round(summary(m)$r.squared, digits = 3) , sep = "")) +
    #   theme_bw()
    # p
    # ggsave(paste("test_predict.pdf", sep=""), p, width = 4, height = 2.5)
    # 
    # imp= as.data.frame(rf$importance)
    # imp = imp[order(imp[,1],decreasing = T),]
    
    imp_i<-importance(rf)
    imp_i<-as.data.frame(imp_i)
    imp_i$variable<-rownames(imp_i)
    times2<-rep(b,times=nrow(imp_i))
    imp_i<-data.frame(imp_i,times2)
    imp<-rbind(imp,imp_i)
  } 
  names(result)<-c("rho","pvalue","R2")
  # names(result)<-c("times","sensitivity","specifity","auc")
  # names(plot_data)<-c("sensitivity","specificity","Group")
  # names(accura)<-c("Accuracy","Lower","higher","AUC","Time")
  # 
  # write.csv(accura,file=paste0(i,"/genus_selected2/accura.csv"),row.names = T)
  write.csv(imp,file=paste0(i,"/genus_selected2_test/imp.csv"),row.names = T)
  write.csv(df,file=paste0(i,"/genus_selected2_test/df.csv"),row.names = T)
  
  write.csv(result,file=paste0(i,"/genus_selected2_test/result.csv"),row.names = T)
  # write.csv(plot_data,file=paste0(i,"/genus_selected2_test/plot_data.csv"),row.names = T)
  # write.csv(pred,file=paste0(i,"/genus_selected2_test/pred.csv"),row.names = T)
}


for (i in c("CHOL_third","TG_third")) {
  rf_classify(i)
}

#三、 test alone --------------------------------------------------------------
# 
# test_name<-c("BUN", "CREA", "UA", "TP", "ALB", "CHOL", "TG", "ALT", "AST", 
#              "ALP", "LDH", "GGT", "TBIL", "DBIL", "TBA", "FMN", "FBG", "RBP", 
#              "cystain", "HCY")
#impute missing value
library(dplyr)
# data_rf_second<-data_rf_second %>% 
#   mutate_all(function(x){x[is.na(x)] <- mean(x)
#   x})

for (i in test_name) {
  
  data_rf_second[,i][is.na(data_rf_second[,i])==TRUE]<-mean(data_rf_second[,i],na.rm = T)
}

rf_classify<-function(i){
  # i<-"CHOL_third"
  dir.create(i)
  dir.create(paste0(i,"/test_alone"))
  varia<-c(i,test_name)
  m<-100
  result<-data.frame()
  plot_data<-data.frame()
  pred<-data.frame()
  set.seed(123)
  accura<-data.frame()
  result_m<-data.frame()
  imp<-data.frame()
  b=1
  df<-data.frame()
  for (b in 1:m){
    
    Train <- createDataPartition(data_rf_second[,i], p=0.8, list=FALSE)
    training <- data_rf_second[ Train, varia]
    testing <- data_rf_second[ -Train, varia]
    
    rf = randomForest(as.formula(paste0(i,"~.") ), data=training,importance=TRUE, proximity=TRUE, ntree = 1000)
    
    pred_b=predict(rf,testing,type="response")
    df_b = data.frame(observed = testing[,i], predict = pred_b)
    times2<-rep(b,times=nrow(df_b))
    
    
    cor = cor.test(df_b[,1], df_b[,2], method = "spearman")
    df2 = df_b
    colnames(df2) = c("x", "y")
    m = lm(y ~ x, df2)
    rho<-round(cor$estimate, digits = 3)
    pvalue<-round(cor$p.value, digits = 3)
    R2<-round(summary(m)$r.squared, digits = 3)
    result_i<-c(rho,pvalue,R2)
    result<-rbind(result,result_i)
    
    df_b<-data.frame(df_b,times2)
    df<-rbind(df,df_b)
    # p = ggplot(df, aes(predict, observed)) +
    #   geom_point() +
    #   geom_smooth(method = "lm") +
    #   labs(title = paste("rho = " , round(cor$estimate, digits = 3), ", P = " , signif(cor$p.value, digits = 3), ", R2 = ", round(summary(m)$r.squared, digits = 3) , sep = "")) +
    #   theme_bw()
    # p
    # ggsave(paste("test_predict.pdf", sep=""), p, width = 4, height = 2.5)
    # 
    # imp= as.data.frame(rf$importance)
    # imp = imp[order(imp[,1],decreasing = T),]
    
    imp_i<-importance(rf)
    imp_i<-as.data.frame(imp_i)
    imp_i$variable<-rownames(imp_i)
    times2<-rep(b,times=nrow(imp_i))
    imp_i<-data.frame(imp_i,times2)
    imp<-rbind(imp,imp_i)
  } 
  names(result)<-c("rho","pvalue","R2")
  # names(result)<-c("times","sensitivity","specifity","auc")
  # names(plot_data)<-c("sensitivity","specificity","Group")
  # names(accura)<-c("Accuracy","Lower","higher","AUC","Time")
  # 
  # write.csv(accura,file=paste0(i,"/genus_selected2/accura.csv"),row.names = T)
  write.csv(imp,file=paste0(i,"/test_alone/imp.csv"),row.names = T)
  write.csv(df,file=paste0(i,"/test_alone/df.csv"),row.names = T)
  
  write.csv(result,file=paste0(i,"/test_alone/result.csv"),row.names = T)
  # write.csv(plot_data,file=paste0(i,"/test_alone/plot_data.csv"),row.names = T)
  # write.csv(pred,file=paste0(i,"/test_alone/pred.csv"),row.names = T)
}


for (i in c("CHOL_third","TG_third")) {
  rf_classify(i)
}


#四、 CAG level ---------------------------------------------------------------
setwd("E:/study2/microbiome/study_new/total/9Data Review")
#setwd('E:/study2/microbiome/study_new/new2/10WGCNA')
data_CAG_combine<-read.csv("Data/16S/CAG/data_CAG_combine.csv",header = T,sep = ",")
rownames(data_CAG_combine)<-data_CAG_combine$X
data_CAG_combine$ID<-data_CAG_combine$X
# metadata<-read.xlsx("metadata.xlsx",sheet = 4)
dat_CAG<-merge(metadata,data_CAG_combine,by.x = "ID",by.y = "X")
#write.csv(dat_CAG,file = "dat_CAG.csv")
dat_CAG_24<-subset(dat_CAG,Time=="24week")
dat_CAG_32<-subset(dat_CAG,Time=="32week")

CHOL_third<-dat_CAG_32$CHOL
TG_third<-dat_CAG_32$TG
#hyperlipid_third<-dat_CAG_32$Hyperlipid
dat_CAG_24<-data.frame(dat_CAG_24,CHOL_third,TG_third)
CAG_name<-paste("CAG",seq(1,5),sep = "")

# dat_CAG_24$CHOL_third<-as.factor(dat_CAG_24$CHOL_third)
# dat_CAG_24$TG_third<-as.factor(dat_CAG_24$TG_third)
# dat_CAG_24$hyperlipid_third<-as.factor(dat_CAG_24$hyperlipid_third)
# 
# test_name<-c("BUN", "CREA", "UA", "TP", "ALB", "CHOL", "TG", "ALT", "AST", 
#              "ALP", "LDH", "GGT", "TBIL", "DBIL", "TBA", "FMN", "FBG", "RBP", 
#              "cystain", "HCY")
#impute missing value
library(dplyr)
# data_rf_second<-data_rf_second %>% 
#   mutate_all(function(x){x[is.na(x)] <- mean(x)
#   x})

for (i in test_name) {
  
  dat_CAG_24[,i][is.na(dat_CAG_24[,i])==TRUE]<-mean(dat_CAG_24[,i],na.rm = T)
}
dir_root<-"E:/study2/microbiome/study_new/total/9Data Review/Result/ML"
setwd(dir_root)
# #1. CAG alone ------------------------------------------------------------
rf_classify<-function(i){
  # i<-"CHOL_third"
  dir.create(i)
  dir.create(paste0(i,"/CAG_alone"))
  varia<-c(i,CAG_name)
  m<-100
  result<-data.frame()
  plot_data<-data.frame()
  pred<-data.frame()
  set.seed(123)
  accura<-data.frame()
  result_m<-data.frame()
  imp<-data.frame()
  b=1
  df<-data.frame()
  for (b in 1:m){
    
    Train <- createDataPartition(dat_CAG_24[,i], p=0.8, list=FALSE)
    training <- dat_CAG_24[ Train, varia]
    testing <- dat_CAG_24[ -Train, varia]
    
    rf = randomForest(as.formula(paste0(i,"~.") ), data=training,importance=TRUE, proximity=TRUE, ntree = 1000)
    
    pred_b=predict(rf,testing,type="response")
    df_b = data.frame(observed = testing[,i], predict = pred_b)
    times2<-rep(b,times=nrow(df_b))
    
    
    cor = cor.test(df_b[,1], df_b[,2], method = "spearman")
    df2 = df_b
    colnames(df2) = c("x", "y")
    m = lm(y ~ x, df2)
    rho<-round(cor$estimate, digits = 3)
    pvalue<-round(cor$p.value, digits = 3)
    R2<-round(summary(m)$r.squared, digits = 3)
    result_i<-c(rho,pvalue,R2)
    result<-rbind(result,result_i)
    
    df_b<-data.frame(df_b,times2)
    df<-rbind(df,df_b)
    # p = ggplot(df, aes(predict, observed)) +
    #   geom_point() +
    #   geom_smooth(method = "lm") +
    #   labs(title = paste("rho = " , round(cor$estimate, digits = 3), ", P = " , signif(cor$p.value, digits = 3), ", R2 = ", round(summary(m)$r.squared, digits = 3) , sep = "")) +
    #   theme_bw()
    # p
    # ggsave(paste("test_predict.pdf", sep=""), p, width = 4, height = 2.5)
    # 
    # imp= as.data.frame(rf$importance)
    # imp = imp[order(imp[,1],decreasing = T),]
    
    imp_i<-importance(rf)
    imp_i<-as.data.frame(imp_i)
    imp_i$variable<-rownames(imp_i)
    times2<-rep(b,times=nrow(imp_i))
    imp_i<-data.frame(imp_i,times2)
    imp<-rbind(imp,imp_i)
  } 
  names(result)<-c("rho","pvalue","R2")
  # names(result)<-c("times","sensitivity","specifity","auc")
  # names(plot_data)<-c("sensitivity","specificity","Group")
  # names(accura)<-c("Accuracy","Lower","higher","AUC","Time")
  # 
  # write.csv(accura,file=paste0(i,"/genus_selected2/accura.csv"),row.names = T)
  write.csv(imp,file=paste0(i,"/CAG_alone/imp.csv"),row.names = T)
  write.csv(df,file=paste0(i,"/CAG_alone/df.csv"),row.names = T)
  
  write.csv(result,file=paste0(i,"/CAG_alone/result.csv"),row.names = T)
  # write.csv(plot_data,file=paste0(i,"/CAG_alone/plot_data.csv"),row.names = T)
  # write.csv(pred,file=paste0(i,"/CAG_alone/pred.csv"),row.names = T)
}


for (i in c("CHOL_third","TG_third")) {
  rf_classify(i)
}


# 2.CAG+test -------------------------------------------------------------

rf_classify<-function(i){
  # i<-"CHOL_third"
  dir.create(i)
  dir.create(paste0(i,"/CAG_test"))
  varia<-c(i,CAG_name,test_name)
  m<-100
  result<-data.frame()
  plot_data<-data.frame()
  pred<-data.frame()
  set.seed(123)
  accura<-data.frame()
  result_m<-data.frame()
  imp<-data.frame()
  b=1
  df<-data.frame()
  for (b in 1:m){
    
    Train <- createDataPartition(dat_CAG_24[,i], p=0.8, list=FALSE)
    training <- dat_CAG_24[ Train, varia]
    testing <- dat_CAG_24[ -Train, varia]
    
    rf = randomForest(as.formula(paste0(i,"~.") ), data=training,importance=TRUE, proximity=TRUE, ntree = 1000)
    
    pred_b=predict(rf,testing,type="response")
    df_b = data.frame(observed = testing[,i], predict = pred_b)
    times2<-rep(b,times=nrow(df_b))
    
    
    cor = cor.test(df_b[,1], df_b[,2], method = "spearman")
    df2 = df_b
    colnames(df2) = c("x", "y")
    m = lm(y ~ x, df2)
    rho<-round(cor$estimate, digits = 3)
    pvalue<-round(cor$p.value, digits = 3)
    R2<-round(summary(m)$r.squared, digits = 3)
    result_i<-c(rho,pvalue,R2)
    result<-rbind(result,result_i)
    
    df_b<-data.frame(df_b,times2)
    df<-rbind(df,df_b)
    # p = ggplot(df, aes(predict, observed)) +
    #   geom_point() +
    #   geom_smooth(method = "lm") +
    #   labs(title = paste("rho = " , round(cor$estimate, digits = 3), ", P = " , signif(cor$p.value, digits = 3), ", R2 = ", round(summary(m)$r.squared, digits = 3) , sep = "")) +
    #   theme_bw()
    # p
    # ggsave(paste("test_predict.pdf", sep=""), p, width = 4, height = 2.5)
    # 
    # imp= as.data.frame(rf$importance)
    # imp = imp[order(imp[,1],decreasing = T),]
    
    imp_i<-importance(rf)
    imp_i<-as.data.frame(imp_i)
    imp_i$variable<-rownames(imp_i)
    times2<-rep(b,times=nrow(imp_i))
    imp_i<-data.frame(imp_i,times2)
    imp<-rbind(imp,imp_i)
  } 
  names(result)<-c("rho","pvalue","R2")
  # names(plot_data)<-c("sensitivity","specificity","Group")
  # names(accura)<-c("Accuracy","Lower","higher","AUC","Time")
  # 
  # write.csv(accura,file=paste0(i,"/genus_selected2/accura.csv"),row.names = T)
  write.csv(imp,file=paste0(i,"/CAG_test/imp.csv"),row.names = T)
  write.csv(df,file=paste0(i,"/CAG_test/df.csv"),row.names = T)
  
  write.csv(result,file=paste0(i,"/CAG_test/result.csv"),row.names = T)
  # write.csv(plot_data,file=paste0(i,"/CAG_test/plot_data.csv"),row.names = T)
  # write.csv(pred,file=paste0(i,"/CAG_test/pred.csv"),row.names = T)
}


for (i in c("CHOL_third","TG_third")) {
  rf_classify(i)
}


#五、 ROC_95%CI ---------------------------------------------------------------
library(ggplot2)
library(pROC)
dir_root<-"E:/study2/microbiome/study_new/total/9Data Review/Result/ML"
setwd(dir_root)
dir.create("combine")
#1. ROC curve----------------------------------------------------

plot_result<-function(i,j){
  dir1<-paste0(dir_root,"/",i,"/",j)
  
  result<-read.csv(paste0(dir1,"/result.csv"),header = T,sep = ",")
 # plot_data<-read.csv(paste0(dir1,"/plot_data_",k,".csv"),header = T,sep = ",")
  pred<-read.csv(paste0(dir1,"/df.csv"),header = T,sep = ",")
  
  result1<-result[order(result$rho,decreasing = T),]
  imp<-read.csv(paste0(dir1,"/imp.csv"),header = T,sep = ",")

  time<-result1[50,1]
  set.seed(315)
  pred1<-pred[pred$times2==time,]
  cor = cor.test(pred1$predict, pred1$observed, method = "spearman")
  df2 = pred1[,c("predict","observed")]
  colnames(df2) = c("x", "y")
  m = lm(y ~ x, df2)
  p = ggplot(pred1, aes(predict, observed)) +
    geom_point() +
    geom_smooth(method = "lm",color="darkblue") +
    labs(title = paste("rho = " , round(cor$estimate, digits = 3), ", P = " , signif(cor$p.value, digits = 3), ", R2 = ", round(summary(m)$r.squared, digits = 3) , sep = "")) +
    theme_bw()
  p
  ggsave(paste0(dir_root,"/",i,"/",j,"/",i,"_",j,"_reg_median.pdf"),width=5,height=5)
  
  
 ###top value
  time<-result1[1,1]
  set.seed(315)
  pred1<-pred[pred$times2==time,]
  cor = cor.test(pred1$predict, pred1$observed, method = "spearman")
  df2 = pred1[,c("predict","observed")]
  colnames(df2) = c("x", "y")
  m = lm(y ~ x, df2)
  p = ggplot(pred1, aes(predict, observed)) +
    geom_point() +
    geom_smooth(method = "lm",color="darkblue") +
    labs(title = paste("rho = " , round(cor$estimate, digits = 3), ", P = " , signif(cor$p.value, digits = 3), ", R2 = ", round(summary(m)$r.squared, digits = 3) , sep = "")) +
    theme_bw()
  p
  ggsave(paste0(dir_root,"/",i,"/",j,"/",i,"_",j,"_reg_top.pdf"),width=5,height=5)
  
  ##importance
  
 # imp<-read.csv(paste0(dir1,"/imp.csv"),header = T,sep = ",")
  t<-tapply(imp$X.IncMSE,imp$variable,mean)##get mean value of imp
  t<-t[order(t,decreasing = T)]
  var_name<-names(t)
  ggplot(imp,aes(x=variable,y=X.IncMSE))+
    geom_boxplot()+
    theme_bw()+
    scale_x_discrete(limits=rev(var_name))+
    # theme(axis.text.x=element_text(angle=45,size=8,vjust = 1,hjust = 1))+
    coord_flip()
  ggsave(paste0(dir_root,"/",i,"/",j,"/",i,"_",j,"_varimp.pdf"),width=10,height=10)
  
  
}

i<-"TG_third"
#j<-"Family_alone"

outome<-c("TG_third","CHOL_third")
var2<-c("genus_selected2_test","genus_selected2","CAG_test","CAG_alone","test_alone")

for (i in outome) {
  for (j in var2) {
      plot_result(i,j)
  }
}

# 2.combine reg linear plot & varimportance plot-----------------------------------------------------

i<-"TG_third"
#j<-"Family_alone"

outome<-c("TG_third","CHOL_third")
var2<-c("genus_selected2_test","genus_selected2","CAG_test","CAG_alone","test_alone")

imp_total<-data.frame()
pred_top<-data.frame()
pred_median<-data.frame()

for (i in outome) {
  for (j in var2) {
    dir1<-paste0(dir_root,"/",i,"/",j)
    
    result<-read.csv(paste0(dir1,"/result.csv"),header = T,sep = ",")
    # plot_data<-read.csv(paste0(dir1,"/plot_data_",k,".csv"),header = T,sep = ",")
    pred<-read.csv(paste0(dir1,"/df.csv"),header = T,sep = ",")
    
    result1<-result[order(result$rho,decreasing = T),]
    ####imp
    imp<-read.csv(paste0(dir1,"/imp.csv"),header = T,sep = ",")
    Group<-rep(i,times=nrow(imp))
    variab<-rep(j,times=nrow(imp))
    imp_i<-data.frame(imp,Group,variab)
    imp_total<-rbind(imp_total,imp_i)
    time<-result1[50,1]
    set.seed(315)
    pred1<-pred[pred$times2==time,]
    
    Group<-rep(i,times=nrow(pred1))
    variab<-rep(j,times=nrow(pred1))
    imp_i<-data.frame(pred1,Group,variab)
    
    pred_median<-rbind(pred_median,pred1)
    
    time<-result1[1,1]
    set.seed(315)
    pred1<-pred[pred$times2==time,]
    Group<-rep(i,times=nrow(pred1))
    variab<-rep(j,times=nrow(pred1))
    imp_i<-data.frame(pred1,Group,variab)
    pred_top<-rbind(pred_top,pred1)
    
    
    
  }
}


write.csv(pred_top,file = "combine//reg_pred_top.csv")
write.csv(pred_median,file = "combine//reg_pred_median.csv")
write.csv(imp_total,file = "combine//reg_imp_total.csv")

library(ggsci)

for (j in var2) {
  imp_plot<-subset(imp_total,variab==j)
  ggplot(imp_plot,aes(x=variable,y=X.IncMSE,fill=Group))+
    geom_boxplot()+
    theme_bw()+
    scale_fill_nejm()+
   # facet_wrap(~variab,ncol = 4)+
    #  scale_x_discrete(limits=rev(var_name))+
    # theme(axis.text.x=element_text(angle=45,size=8,vjust = 1,hjust = 1))+
    coord_flip()
  ggsave(paste0("combine/varimp_reg_",j,".pdf"),width=10,height=12)
  
}


