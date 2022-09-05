

dir_root<-"E:/study2/microbiome/study_new/total/9Data Review/Result/ML"
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


# top ---------------------------------------------------------------------


i<-"Hypertriglyceridemia_third"
j<-"CAG_test"
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
j<-"CAG_test"
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
j<-"CAG_test"
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

ggsave(Rocpics,filename = "combine/roc_combine_CAG_top.pdf",width = 10,height = 8)

#2. median ---------------------------------------------------------------------

# CAG ---------------------------------------------------------------------



i<-"Hypertriglyceridemia_third"
j<-"CAG_test"
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
j<-"CAG_test"
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
j<-"CAG_test"
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

ggsave(Rocpics,filename = "combine/roc_combine_CAG_median.pdf",width = 10,height = 8)

# Genera ---------------------------------------------------------------------


i<-"Hypertriglyceridemia_third"
j<-"genus_selected2_test"
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
j<-"genus_selected2_test"
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
j<-"genus_selected2_test"
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

ggsave(Rocpics,filename = "combine/roc_combine_genera_median.pdf",width = 10,height = 8)

# 3.combine genus,test and genus test  -----------------------------------------------

# 3.1 Hypertriglyceridemia_third ----------------------------------------------

i<-"Hypertriglyceridemia_third"
j<-"genus_selected2"
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
Genera_alone<-roc(pred1$groundtruth,pred1$X1,auc = T,ci= T,print.auc=T)



i<-"Hypertriglyceridemia_third"
j<-"test_alone"
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
Test_alone<-roc(pred1$groundtruth,pred1$X1,auc = T,ci= T,print.auc=T)





i<-"Hypertriglyceridemia_third"
j<-"genus_selected2_test"
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
Combined<-roc(pred1$groundtruth,pred1$X1,auc = T,ci= T,print.auc=T)

rocLs<-c(Genera_alone,Test_alone,Combined)

rocLs = list('Genera_alone'=Genera_alone, 
             'Test_alone' = Test_alone,"Combined"=Combined)

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
  rocciDat$Grp<-ordered(rocciDat$Grp,levels=c("Genera_alone","Test_alone",
                                              "Combined"))
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

ggsave(Rocpics,filename = "combine/roc_combine_Hypertriglyceridemia.pdf",width = 10,height = 8)

# 3.2 Hypercholesterolemia_third ----------------------------------------------

i<-"Hypercholesterolemia_third"
j<-"genus_selected2"
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
Genera_alone<-roc(pred1$groundtruth,pred1$X1,auc = T,ci= T,print.auc=T)



i<-"Hypercholesterolemia_third"
j<-"test_alone"
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
Test_alone<-roc(pred1$groundtruth,pred1$X1,auc = T,ci= T,print.auc=T)





i<-"Hypercholesterolemia_third"
j<-"genus_selected2_test"
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
Combined<-roc(pred1$groundtruth,pred1$X1,auc = T,ci= T,print.auc=T)

rocLs<-c(Genera_alone,Test_alone,Combined)

rocLs = list('Genera_alone'=Genera_alone, 
             'Test_alone' = Test_alone,"Combined"=Combined)

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
  rocciDat$Grp<-ordered(rocciDat$Grp,levels=c("Genera_alone","Test_alone",
                                              "Combined"))
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

ggsave(Rocpics,filename = "combine/roc_combine_Hypercholesterolemia.pdf",width = 10,height = 8)

# 3.3 Hyperlipid_third ----------------------------------------------

i<-"Hyperlipid_third"
j<-"genus_selected2"
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
Genera_alone<-roc(pred1$groundtruth,pred1$X1,auc = T,ci= T,print.auc=T)



i<-"Hyperlipid_third"
j<-"test_alone"
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
Test_alone<-roc(pred1$groundtruth,pred1$X1,auc = T,ci= T,print.auc=T)





i<-"Hyperlipid_third"
j<-"genus_selected2_test"
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
Combined<-roc(pred1$groundtruth,pred1$X1,auc = T,ci= T,print.auc=T)

rocLs<-c(Genera_alone,Test_alone,Combined)

rocLs = list('Genera_alone'=Genera_alone, 
             'Test_alone' = Test_alone,"Combined"=Combined)

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
  rocciDat$Grp<-ordered(rocciDat$Grp,levels=c("Genera_alone","Test_alone",
                                              "Combined"))
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

ggsave(Rocpics,filename = "combine/roc_combine_hyperlipid.pdf",width = 10,height = 8)

