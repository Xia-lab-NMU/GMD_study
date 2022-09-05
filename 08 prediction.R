library(openxlsx)
library(randomForest)
library(ggplot2)
library(ggsci)
library(openxlsx)
library(ggplot2)

setwd("E:/study2/microbiome/study_new/total/9Data Review")
dir.create("Result/ML")
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

# 一、genus_level_selected1 ------------------------------------------------------------
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

Hypercholesterolemia_third<-data_rf_third$CHOL_High
Hypertriglyceridemia_third<-data_rf_third$TG_High
hyperlipid_third<-data_rf_third$Hyperlipid
#Hypercholesterolemia_third<-data_rf_third$CHOL_High


data_rf_second<-data.frame(data_rf_second,Hypercholesterolemia_third,Hypertriglyceridemia_third,
                           hyperlipid_third )
data_rf_second$Hypercholesterolemia_third<-as.factor(data_rf_second$Hypercholesterolemia_third)
data_rf_second$Hypertriglyceridemia_third<-as.factor(data_rf_second$Hypertriglyceridemia_third)
data_rf_second$hyperlipid_third<-as.factor(data_rf_second$hyperlipid_third)

dir_root<-"E:/study2/microbiome/study_new/total/9Data Review/Result/ML"
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
library(caret)
library(ROCR)
library(pROC)
#1.Genus_selected2 -------------------------------------------------------

# genus_selected<-c("Bacteroides", "Lactobacillus", "Monoglobus", "Oscillospiraceae.uncultured", 
#                   "Agathobacter", "Alistipes", "Allisonella", "Butyricimonas", 
#                   "Christensenellaceae_R_7_group", "Collinsella", "Enterococcus", 
#                   "Erysipelatoclostridium", "Eubacterium_ruminantium_group", "Prevotella", 
#                   "TM7x", "UCG_010", "Butyricicoccus", "Klebsiella", "Megasphaera")

genus_selected2<-c("Alistipes", "Bacteroides","Paraprevotella",
                   "Christensenellaceae_R_7_group", "UCG_002","Clostridia_UCG_014")


# 
# #1. genus alone ------------------------------------------------------------
rf_classify<-function(i,j){
  # i<-"Hypercholesterolemia_third"
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
    err_b<-rf$err.rate[500,1]
    result_b<-c(b,sensi_b,speci_b,auc_b,err_b)
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
  names(result)<-c("times","sensitivity","specifity","auc","error")
  names(plot_data)<-c("sensitivity","specificity","Group")
  names(accura)<-c("Accuracy","Lower","higher","AUC","Time")
  
  write.csv(accura,file=paste0(i,"/genus_selected2/accura_",ratio,".csv"),row.names = T)
  write.csv(imp,file=paste0(i,"/genus_selected2/imp_",ratio,".csv"),row.names = T)
  
  write.csv(result,file=paste0(i,"/genus_selected2/result_",ratio,".csv"),row.names = T)
  write.csv(plot_data,file=paste0(i,"/genus_selected2/plot_data_",ratio,".csv"),row.names = T)
  write.csv(pred,file=paste0(i,"/genus_selected2/pred_",ratio,".csv"),row.names = T)
}



for (i in c("Hypercholesterolemia_third","Hypertriglyceridemia_third","hyperlipid_third")) {
  for (j in c(0.8)) {
    rf_classify(i,j)
  }
  
}


# 2.genus+test -------------------------------------------------------------

rf_classify<-function(i,j){
  # i<-"Hypercholesterolemia_third"
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
    err_b<-rf$err.rate[500,1]
    result_b<-c(b,sensi_b,speci_b,auc_b,err_b)
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
  names(result)<-c("times","sensitivity","specifity","auc","error")
  names(plot_data)<-c("sensitivity","specificity","Group")
  names(accura)<-c("Accuracy","Lower","higher","AUC","Time")
  
  write.csv(accura,file=paste0(i,"/genus_selected2_test/accura_",ratio,".csv"),row.names = T)
  write.csv(imp,file=paste0(i,"/genus_selected2_test/imp_",ratio,".csv"),row.names = T)
  
  write.csv(result,file=paste0(i,"/genus_selected2_test/result_",ratio,".csv"),row.names = T)
  write.csv(plot_data,file=paste0(i,"/genus_selected2_test/plot_data_",ratio,".csv"),row.names = T)
  write.csv(pred,file=paste0(i,"/genus_selected2_test/pred_",ratio,".csv"),row.names = T)
}



for (i in c("Hypercholesterolemia_third","Hypertriglyceridemia_third","hyperlipid_third")) {
  for (j in c(0.8)) {
    rf_classify(i,j)
  }
  
}

#三、 test alone --------------------------------------------------------------

# test_name<-c("BUN", "CREA", "UA", "TP", "ALB", "CHOL", "TG", "ALT", "AST", 
#              "ALP", "LDH", "GGT", "TBIL", "DBIL", "TBA", "FMN", "FBG", "RBP", 
#              "cystain", "HCY")
#impute missing value
# library(dplyr)
# # data_rf_second<-data_rf_second %>% 
# #   mutate_all(function(x){x[is.na(x)] <- mean(x)
# #   x})
# 
# for (i in test_name) {
#   
#   data_rf_second[,i][is.na(data_rf_second[,i])==TRUE]<-mean(data_rf_second[,i],na.rm = T)
# }
rf_classify<-function(i,j){
  # i<-"Hypercholesterolemia_third"
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
    err_b<-rf$err.rate[500,1]
    result_b<-c(b,sensi_b,speci_b,auc_b,err_b)
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
  names(result)<-c("times","sensitivity","specifity","auc","error")
  names(plot_data)<-c("sensitivity","specificity","Group")
  names(accura)<-c("Accuracy","Lower","higher","AUC","Time")
  
  write.csv(accura,file=paste0(i,"/test_alone/accura_",ratio,".csv"),row.names = T)
  write.csv(imp,file=paste0(i,"/test_alone/imp_",ratio,".csv"),row.names = T)
  
  write.csv(result,file=paste0(i,"/test_alone/result_",ratio,".csv"),row.names = T)
  write.csv(plot_data,file=paste0(i,"/test_alone/plot_data_",ratio,".csv"),row.names = T)
  write.csv(pred,file=paste0(i,"/test_alone/pred_",ratio,".csv"),row.names = T)
}



for (i in c("Hypercholesterolemia_third","Hypertriglyceridemia_third","hyperlipid_third")) {
  for (j in c(0.8)) {
    rf_classify(i,j)
  }
  
}


#四、 CAG level ---------------------------------------------------------------

setwd(source_dir)
data_CAG_combine<-read.csv("Data/16S/CAG/data_CAG_combine.csv",header = T,sep = ",")
rownames(data_CAG_combine)<-data_CAG_combine$X
data_CAG_combine$ID<-data_CAG_combine$X
# metadata<-read.xlsx("metadata.xlsx",sheet = 4)
dat_CAG<-merge(metadata,data_CAG_combine,by.x = "ID",by.y = "X")
#write.csv(dat_CAG,file = "dat_CAG.csv")
dat_CAG_24<-subset(dat_CAG,Time=="24week")
dat_CAG_32<-subset(dat_CAG,Time=="32week")

Hypercholesterolemia_third<-dat_CAG_32$CHOL_High
Hypertriglyceridemia_third<-dat_CAG_32$TG_High
hyperlipid_third<-dat_CAG_32$Hyperlipid
dat_CAG_24<-data.frame(dat_CAG_24,Hypercholesterolemia_third,Hypertriglyceridemia_third,hyperlipid_third)
CAG_name<-paste("CAG",seq(1,5),sep = "")

dat_CAG_24$Hypercholesterolemia_third<-as.factor(dat_CAG_24$Hypercholesterolemia_third)
dat_CAG_24$Hypertriglyceridemia_third<-as.factor(dat_CAG_24$Hypertriglyceridemia_third)
dat_CAG_24$hyperlipid_third<-as.factor(dat_CAG_24$hyperlipid_third)

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


setwd(dir_root)

# #1. CAG alone ------------------------------------------------------------
rf_classify<-function(i,j){
  # i<-"Hypercholesterolemia_third"
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
  ratio=j
  for (b in 1:m){
    
    Train <- createDataPartition(dat_CAG_24[,i], p=ratio, list=FALSE)
    training <- dat_CAG_24[ Train, varia]
    testing <- dat_CAG_24[ -Train, varia]
    
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
    err_b<-rf$err.rate[500,1]
    result_b<-c(b,sensi_b,speci_b,auc_b,err_b)
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
  names(result)<-c("times","sensitivity","specifity","auc","error")
  names(plot_data)<-c("sensitivity","specificity","Group")
  names(accura)<-c("Accuracy","Lower","higher","AUC","Time")
  
  write.csv(accura,file=paste0(i,"/CAG_alone/accura_",ratio,".csv"),row.names = T)
  write.csv(imp,file=paste0(i,"/CAG_alone/imp_",ratio,".csv"),row.names = T)
  
  write.csv(result,file=paste0(i,"/CAG_alone/result_",ratio,".csv"),row.names = T)
  write.csv(plot_data,file=paste0(i,"/CAG_alone/plot_data_",ratio,".csv"),row.names = T)
  write.csv(pred,file=paste0(i,"/CAG_alone/pred_",ratio,".csv"),row.names = T)
}



for (i in c("Hypercholesterolemia_third","Hypertriglyceridemia_third","hyperlipid_third")) {
  for (j in c(0.8)) {
    rf_classify(i,j)
  }
  
}

# 2.CAG+test -------------------------------------------------------------

rf_classify<-function(i,j){
  # i<-"Hypercholesterolemia_third"
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
  ratio=j
  for (b in 1:m){
    
    Train <- createDataPartition(dat_CAG_24[,i], p=ratio, list=FALSE)
    training <- dat_CAG_24[ Train, varia]
    testing <- dat_CAG_24[ -Train, varia]
    
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
    err_b<-rf$err.rate[500,1]
    result_b<-c(b,sensi_b,speci_b,auc_b,err_b)
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
  names(result)<-c("times","sensitivity","specifity","auc","error")
  names(plot_data)<-c("sensitivity","specificity","Group")
  names(accura)<-c("Accuracy","Lower","higher","AUC","Time")
  
  write.csv(accura,file=paste0(i,"/CAG_test/accura_",ratio,".csv"),row.names = T)
  write.csv(imp,file=paste0(i,"/CAG_test/imp_",ratio,".csv"),row.names = T)
  
  write.csv(result,file=paste0(i,"/CAG_test/result_",ratio,".csv"),row.names = T)
  write.csv(plot_data,file=paste0(i,"/CAG_test/plot_data_",ratio,".csv"),row.names = T)
  write.csv(pred,file=paste0(i,"/CAG_test/pred_",ratio,".csv"),row.names = T)
}



for (i in c("Hypercholesterolemia_third","Hypertriglyceridemia_third","hyperlipid_third")) {
  for (j in c(0.8)) {
    rf_classify(i,j)
  }
  
}

#五、 ROC_95%CI ---------------------------------------------------------------

dir_root<-"E:/study2/microbiome/study_new/total/9Data Review/Result/ML"
setwd(dir_root)
#1. ROC curve----------------------------------------------------

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
j<-"Family_test"

outome<-c("Hypercholesterolemia_third","Hypertriglyceridemia_third","hyperlipid_third")
var2<-c("genus_selected2_test","genus_selected2","CAG_test","CAG_alone","test_alone")
ratio<-c(0.8)
for (i in outome) {
  for (j in var2) {
    for (k in ratio) {
      plot_result(i,j,k)
    }
    
  }
  
}

#2. varImportance -----------------------------------------------------------

plot_imp<-function(i,j,k){
  
  dir1<-paste0(dir_root,"/",i,"/",j)
  
  #result<-read.csv(paste0(dir1,"/result_",k,".csv"),header = T,sep = ",")
  #plot_data<-read.csv(paste0(dir1,"/plot_data_",k,".csv"),header = T,sep = ",")
  #pred<-read.csv(paste0(dir1,"/pred_",k,".csv"),header = T,sep = ",")
  
  imp<-read.csv(paste0(dir1,"/imp_",k,".csv"),header = T,sep = ",")
  imp$variable<-as.character(imp$variable)
  t<-tapply(imp$MeanDecreaseGini,imp$variable,mean)##get mean value of imp
  t<-t[order(t,decreasing = T)]
  var_name<-names(t)
  ggplot(imp,aes(x=variable,y=MeanDecreaseGini))+
    geom_boxplot()+
    theme_bw()+
    scale_x_discrete(limits=rev(var_name))+
    # theme(axis.text.x=element_text(angle=45,size=8,vjust = 1,hjust = 1))+
    coord_flip()
  ggsave(paste0(dir_root,"/",i,"/",j,"/",i,"_",j,"_",k,"_varimp.pdf"),width=10,height=10)
  
}

outome<-c("Hypercholesterolemia_third","Hypertriglyceridemia_third","hyperlipid_third")
var2<-c("genus_selected2_test","genus_selected2")
ratio<-c(0.8)
for (i in outome){
  for (j in var2) {
    for (k in ratio) {
      plot_imp(i,j,k)
    }
  }
}




# 3.combine ROC curve & varimportance plot-----------------------------------------------------

i<-"Hypercholesterolemia_third"
j<-"Family_test"

outome<-c("Hypercholesterolemia_third","Hypertriglyceridemia_third","hyperlipid_third")
var2<-c("genus_selected2_test","genus_selected2","CAG_test","CAG_alone","test_alone")

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





#4. auc --------------------------------------------------------------------

dir_root<-"E:/study2/microbiome/study_new/total/9Data Review/Result/ML"
setwd(dir_root)
outome<-c("Hypercholesterolemia_third","Hypertriglyceridemia_third","hyperlipid_third")
var2<-c("genus_selected2_test","genus_selected2","CAG_test","CAG_alone","test_alone")
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
library(ggpubr)
setwd("E:/study2/microbiome/study_new/total/9Data Review/Result/ML")
dir.create("combine")
setwd("E:/study2/microbiome/study_new/total/9Data Review/Result/ML/combine")
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



