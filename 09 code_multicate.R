
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

genus_filter<-read.csv("Data/genus_filter_RA.csv",header = T,row.names = 1)
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

Dyslipidemia<-data_rf_third$dyslipimedia

data_rf_second<-data.frame(data_rf_second,Dyslipidemia)

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

#1.Genus_selected2 -------------------------------------------------------

# genus_selected<-c("Bacteroides", "Lactobacillus", "Monoglobus", "Oscillospiraceae.uncultured", 
#                   "Agathobacter", "Alistipes", "Allisonella", "Butyricimonas", 
#                   "Christensenellaceae_R_7_group", "Collinsella", "Enterococcus", 
#                   "Erysipelatoclostridium", "Eubacterium_ruminantium_group", "Prevotella", 
#                   "TM7x", "UCG_010", "Butyricicoccus", "Klebsiella", "Megasphaera")

genus_selected2<-c("Alistipes", "Bacteroides","Paraprevotella",
                   "Christensenellaceae_R_7_group", "UCG_002","Clostridia_UCG_014")

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
  dir.create(paste0(i,"/genus_selected2"))
  #varia<-c(i,index)
  varia<-c(i,genus_selected2)
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
  
  write.csv(accura,file=paste0(i,"/genus_selected2/accura_",ratio,"_full.csv"),row.names = T)
  write.csv(result,file=paste0(i,"/genus_selected2/result_",ratio,"_full.csv"),row.names = T)
  write.csv(auc_cate,file=paste0(i,"/genus_selected2/auc_cate",ratio,"_full.csv"),row.names = T)
  
  #write.csv(result_m,file=paste0(i,"/genus_selected2/result_m_",ratio,"_full.csv"),row.names = T)
  #write.csv(plot_data,file="Result2/plot_data_bootstrap.csv",row.names = T)
  write.csv(pred,file=paste0(i,"/genus_selected2/pred_",ratio,"_full.csv"),row.names = T)
  write.csv(plot_data,file=paste0(i,"/genus_selected2/plot_data_",ratio,"_full.csv"),row.names = T)
  write.csv(imp,file=paste0(i,"/genus_selected2/imp",ratio,"_full.csv"),row.names = T)
  
}


for (i in c("Dyslipidemia")) {
  for (j in c(0.8)) {
    rf_multiclassify(i,j)
  }
  
}

#2.genus_selected2_test -------------------------------------------------------

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
  dir.create(paste0(i,"/genus_selected2_test"))
  #varia<-c(i,index)
  varia<-c(i,genus_selected2,test_name)
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
  
  write.csv(accura,file=paste0(i,"/genus_selected2_test/accura_",ratio,"_full.csv"),row.names = T)
  write.csv(result,file=paste0(i,"/genus_selected2_test/result_",ratio,"_full.csv"),row.names = T)
  write.csv(auc_cate,file=paste0(i,"/genus_selected2_test/auc_cate",ratio,"_full.csv"),row.names = T)
  
  #write.csv(result_m,file=paste0(i,"/genus_selected2_test/result_m_",ratio,"_full.csv"),row.names = T)
  #write.csv(plot_data,file="Result2/plot_data_bootstrap.csv",row.names = T)
  write.csv(pred,file=paste0(i,"/genus_selected2_test/pred_",ratio,"_full.csv"),row.names = T)
  write.csv(plot_data,file=paste0(i,"/genus_selected2_test/plot_data_",ratio,"_full.csv"),row.names = T)
  write.csv(imp,file=paste0(i,"/genus_selected2_test/imp",ratio,"_full.csv"),row.names = T)
  
}


for (i in "Dyslipidemia") {
  for (j in c(0.8)) {
    rf_multiclassify(i,j)
  }
  
}



#3. test_alone --------------------------------------------------------------

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
  dir.create(paste0(i,"/test_alone"))
  #varia<-c(i,index)
  varia<-c(i,test_name)
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
  
  write.csv(accura,file=paste0(i,"/test_alone/accura_",ratio,"_full.csv"),row.names = T)
  write.csv(result,file=paste0(i,"/test_alone/result_",ratio,"_full.csv"),row.names = T)
  write.csv(auc_cate,file=paste0(i,"/test_alone/auc_cate",ratio,"_full.csv"),row.names = T)
  
  #write.csv(result_m,file=paste0(i,"/test_alone/result_m_",ratio,"_full.csv"),row.names = T)
  #write.csv(plot_data,file="Result2/plot_data_bootstrap.csv",row.names = T)
  write.csv(pred,file=paste0(i,"/test_alone/pred_",ratio,"_full.csv"),row.names = T)
  write.csv(plot_data,file=paste0(i,"/test_alone/plot_data_",ratio,"_full.csv"),row.names = T)
  write.csv(imp,file=paste0(i,"/test_alone/imp",ratio,"_full.csv"),row.names = T)
  
}


for (i in "Dyslipidemia") {
  for (j in c(0.8)) {
    rf_multiclassify(i,j)
  }
  
}



# 二、CAG level -------------------------------------------------------------
setwd("E:/study2/microbiome/study_new/total/9Data Review")
data_CAG_combine<-read.csv("Data/16S/CAG/data_CAG_combine.csv",header = T,sep = ",")
rownames(data_CAG_combine)<-data_CAG_combine$X
data_CAG_combine$ID<-data_CAG_combine$X
# metadata<-read.xlsx("metadata.xlsx",sheet = 4)
dat_CAG<-merge(metadata,data_CAG_combine,by.x = "ID",by.y = "X")
#write.csv(dat_CAG,file = "dat_CAG.csv")
dat_CAG_24<-subset(dat_CAG,Time=="24week")
dat_CAG_32<-subset(dat_CAG,Time=="32week")
CAG_name<-paste("CAG",seq(1,5),sep = "")
Dyslipidemia<-dat_CAG_32$dyslipimedia

data_rf_second<-data.frame(dat_CAG_24,Dyslipidemia)

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
#1.CAG_alone -------------------------------------------------------
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
  dir.create(paste0(i,"/CAG_alone"))
  #varia<-c(i,index)
  varia<-c(i,CAG_name)
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
  
  write.csv(accura,file=paste0(i,"/CAG_alone/accura_",ratio,"_full.csv"),row.names = T)
  write.csv(result,file=paste0(i,"/CAG_alone/result_",ratio,"_full.csv"),row.names = T)
  write.csv(auc_cate,file=paste0(i,"/CAG_alone/auc_cate",ratio,"_full.csv"),row.names = T)
  
  #write.csv(result_m,file=paste0(i,"/CAG_alone/result_m_",ratio,"_full.csv"),row.names = T)
  #write.csv(plot_data,file="Result2/plot_data_bootstrap.csv",row.names = T)
  write.csv(pred,file=paste0(i,"/CAG_alone/pred_",ratio,"_full.csv"),row.names = T)
  write.csv(plot_data,file=paste0(i,"/CAG_alone/plot_data_",ratio,"_full.csv"),row.names = T)
  write.csv(imp,file=paste0(i,"/CAG_alone/imp",ratio,"_full.csv"),row.names = T)
  
}


for (i in "Dyslipidemia") {
  for (j in c(0.8)) {
    rf_multiclassify(i,j)
  }
  
}

#2.CAG_test -------------------------------------------------------

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
  dir.create(paste0(i,"/CAG_test"))
  #varia<-c(i,index)
  varia<-c(i,CAG_name,test_name)
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
  
  write.csv(accura,file=paste0(i,"/CAG_test/accura_",ratio,"_full.csv"),row.names = T)
  write.csv(result,file=paste0(i,"/CAG_test/result_",ratio,"_full.csv"),row.names = T)
  write.csv(auc_cate,file=paste0(i,"/CAG_test/auc_cate",ratio,"_full.csv"),row.names = T)
  
  #write.csv(result_m,file=paste0(i,"/CAG_test/result_m_",ratio,"_full.csv"),row.names = T)
  #write.csv(plot_data,file="Result2/plot_data_bootstrap.csv",row.names = T)
  write.csv(pred,file=paste0(i,"/CAG_test/pred_",ratio,"_full.csv"),row.names = T)
  write.csv(plot_data,file=paste0(i,"/CAG_test/plot_data_",ratio,"_full.csv"),row.names = T)
  write.csv(imp,file=paste0(i,"/CAG_test/imp",ratio,"_full.csv"),row.names = T)
  
}


for (i in "Dyslipidemia") {
  for (j in c(0.8)) {
    rf_multiclassify(i,j)
  }
  
}



#4. combine -----------------------------------------------------------------
#1. auc --------------------------------------------------------------------
library(ggplot2)
library(ggpubr)
library(ggsci)
var2<-c("genus_selected2","genus_selected2_test","CAG_test","CAG_alone","test_alone")
data_full_total<-data.frame()

for (i in var2) {
  dir_root<-paste("E:/study2/microbiome/study_new/total/9Data Review/Result/ML/Dyslipidemia/",i,sep = "")
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

#data_full$group_2<-as.factor(data_full$group_2)

p<-ggboxplot(data_full,x="category",y="AUC",#add="iitter",
             fill = "category",notch = TRUE)+
  stat_compare_means(method = "kruskal.test")+
  theme(legend.position = "none")+
  scale_fill_npg()+
  # scale_fill_manual(values = cbPalette)+
  theme_bw()+
  labs(x="",y="AUC")#+
#  facet_wrap(~group_2,nrow = 2)
setwd("E:/study2/microbiome/study_new/total/9Data Review/Result/ML/Dyslipidemia")
# p<-ggplot() +
#   geom_boxplot( data = data_full,aes(x=group_2,y=AUC,fill=factor(group_2)),outlier.colour = NA) +
#   theme_bw()+
#  # scale_x_discrete(limits=c("full","index","glm"))+
#   labs(x="")
ggsave(p,filename = "boxplot_auc.pdf",width = 8,height = 8) 

p<-ggboxplot(data_full,x="category",y="Accuracy",#add="iitter",
             fill = "category",notch = TRUE)+
  stat_compare_means(method = "kruskal.test")+
  theme(legend.position = "none")+
  scale_fill_npg()+
  # scale_fill_manual(values = cbPalette)+
  theme_bw()+
  labs(x="",y="Accuracy")#+
  #facet_wrap(~group_2,nrow = 2)
setwd("E:/study2/microbiome/study_new/total/9Data Review/Result/ML/Dyslipidemia")
# p<-ggplot() +
#   geom_boxplot( data = data_full,aes(x=group_2,y=AUC,fill=factor(group_2)),outlier.colour = NA) +
#   theme_bw()+
#  # scale_x_discrete(limits=c("full","index","glm"))+
#   labs(x="")
ggsave(p,filename = "boxplot_Accuracy.pdf",width = 8,height = 8) 

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


#2. auc_cate ----------------------------------------------------------------

library(ggplot2)
library(ggpubr)
library(ggsci)
var2<-c("genus_selected2","genus_selected2_test","CAG_test","CAG_alone","test_alone")
data_full_total<-data.frame()

for (i in var2) {
  dir_root<-paste("E:/study2/microbiome/study_new/total/9Data Review/Result/ML/Dyslipidemia/",i,sep = "")
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
setwd("E:/study2/microbiome/study_new/total/9Data Review/Result/ML/Dyslipidemia")

write.csv(data_full_total,file = "data_full_auc_cate.csv",row.names = F)

data_full<-splitstackshape::cSplit(data_full_total, splitCols="group", sep="_")
data_full$group_3<-gsub("[_*.csv]", "", data_full$group_3)

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
  ggsave(p,filename = paste("AUC_figure/",i,"_boxplot.pdf",sep = ""),width = 10,height = 5) 
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


#3. confusion matrix --------------------------------------------------------

# 3.1 according to all type ---------------------------------------------------
group_version<-c("Control","Hypercholesterolemia","Hyperlipidemia","Hypertriglyceridemia")

library(reshape2)
for (i in var2) {
  dir_root<-paste("E:/study2/microbiome/study_new/total/9Data Review/Result/ML/Dyslipidemia/",i,sep = "")
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
  
  
  # data_full_06<-subset(data_full,group_2==0.6)
  # data_full_08<-subset(data_full,group_2==0.8)
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
  
  # sensi<-with(subset(data_full_08), do.call(rbind, tapply(Sensitivity, category, function(x) c(Q = quantile(x,probs = c(0.5,0.025,0.975))))))
  # speci<-with(subset(data_full_08), do.call(rbind, tapply(Specificity, category, function(x) c(Q = quantile(x,probs = c(0.5,0.025,0.975))))))
  # preci<-with(subset(data_full_08), do.call(rbind, tapply(Precision, category, function(x) c(Q = quantile(x,probs = c(0.5,0.025,0.975),na.rm = T)))))
  # reca<-with(subset(data_full_08), do.call(rbind, tapply(Recall, category, function(x) c(Q = quantile(x,probs = c(0.5,0.025,0.975))))))
  # f1<-with(subset(data_full_08), do.call(rbind, tapply(F1, category, function(x) c(Q = quantile(x,probs = c(0.5,0.025,0.975),na.rm = T)))))
  # 
  # result_combine1<-data.frame(sensi,speci,preci,reca,f1)
  # names(result_combine1)<-c("sensi_m","sensi_l","sensi_h","speci_m","speci_l","speci_h","preci_m","preci_l","preci_h",
  #                           "recall_m","recall_l","recall_h","F1_m","F1_l","F1_h")
  # write.csv(result_combine1,file = "result_sensis_08.csv")
  # 
  
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


# 3.2 according to specific category --------------------------------------
# 
# 
# def_plot_08<-function(i,j){
#   plot_data1<-subset(data_full_08,group_3==i)
#   variab<-c("category",j)
#   plot_data1<-melt(plot_data1[variab])
#   p<-ggplot(data=plot_data1,aes(x=category,y=value,color=category))+
#     geom_boxplot(alpha=1, outlier.size=0, size=0.9, width=0.6)+
#     labs(x="", y=j)+theme_bw()+
#     theme(panel.grid.major=element_blank(),legend.position = "none",panel.grid.minor=element_blank())+
#     #  scale_colour_manual(values=cbPalette)
#     scale_color_npg()
#   ggsave(p,filename = paste("Auc_figure/",j,"_",i,"_boxplot_4_1.pdf",sep = ""),width = 6,height = 4) 
# }
# 
# for (j in variabs){
#   def_plot_08("full",j)
# }
# 
# for (j in variabs){
#   def_plot_08("index",j)
# }
# 
# for (j in variabs){
#   def_plot_08("glm",j)
# }
# # 

# 
# accura<-read.csv("accura_bootstrap.csv",header = T,sep = ",")
# accuracy_CI<-quantile(accura$Accuracy,probs = c(0.5,0.025,0.975))
# AUC_CI<-quantile(accura$AUC,probs = c(0.5,0.025,0.975))
# 
# accura<-read.csv("accura_bootstrap_index.csv",header = T,sep = ",")
# accuracy_CI_index<-quantile(accura$Accuracy,probs = c(0.5,0.025,0.975))
# AUC_CI_index<-quantile(accura$AUC,probs = c(0.5,0.025,0.975))
# 
# accura<-read.csv("accura_bootstrap_glm.csv",header = T,sep = ",")
# accuracy_CI_glm<-quantile(accura$Accuracy,probs = c(0.5,0.025,0.975))
# AUC_CI_glm<-quantile(accura$AUC,probs = c(0.5,0.025,0.975))
# 
# result1<-data.frame(accuracy_CI,AUC_CI,accuracy_CI_index,AUC_CI_index,accuracy_CI_glm,AUC_CI_glm)
# result1<-as.data.frame(t(result1))
# write.csv(result1,file = "result_auc.csv")
# 
# result2<-read.csv("result_bootstrap.csv",header = T,sep = ",")



#4. roc curve ---------------------------------------------------------------


# full --------------------------------------------------------------------


for (i in var2) {
  dir_root<-paste("E:/study2/microbiome/study_new/total/9Data Review/Result/ML/Dyslipidemia/",i,sep = "")
  
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

dir_root<-"E:/study2/microbiome/study_new/total/9Data Review/Result/ML"
plot_imp<-function(i,j,k){
  
  dir1<-paste0(dir_root,"/",i,"/",j)
  
  #result<-read.csv(paste0(dir1,"/result_",k,".csv"),header = T,sep = ",")
  #plot_data<-read.csv(paste0(dir1,"/plot_data_",k,".csv"),header = T,sep = ",")
  #pred<-read.csv(paste0(dir1,"/pred_",k,".csv"),header = T,sep = ",")
  
  imp<-read.csv(paste0(dir1,"/imp",k,"_full.csv"),header = T,sep = ",")
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

outome<-c("Dyslipidemia")
var2<-c("genus_selected2","genus_selected2_test","CAG_test","CAG_alone","test_alone")

ratio<-c(0.8)

for (i in outome){
  for (j in var2) {
    for (k in ratio) {
      plot_imp(i,j,k)
    }
  }
}


