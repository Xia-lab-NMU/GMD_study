
library(openxlsx)
#library(xlsx)

# Figure 1b ---------------------------------------------------------------


Original_data<-read.xlsx("Data/Original_data.xlsx",sheet = 1)

Original_data_filter<-subset(Original_data,GDM_pregnancy==0)#n=1262,Diabetes(n=265)
#Original_data_filter<-subset(Original_data,pre_BMI<28)#n=1262,Diabetes(n=265)
Original_data_filter<-subset(Original_data_filter,Tumor==0)#n=1245,tumor=17
Original_data_filter<-subset(Original_data_filter,Thyroid_dysfunction==0)#n=1116,thyroid=129
Original_data_filter<-subset(Original_data_filter,Eclampsia==0)#n=1097, Ecalmpsia=19
Original_data_filter<-subset(Original_data_filter,HBV==0)#1065, HBV=32
Original_data_filter<-subset(Original_data_filter,ICP==0)#1061, ICP=4

Original_data_filter<-subset(Original_data_filter,Lipid_missing==0)#1033, lipid_missing=28

Original_data_filter<-subset(Original_data_filter,`24wFaeces`==1)#976, n=57
Original_data_filter<-subset(Original_data_filter,`32wFaeces`==1)#542, n=434

Original_data_filter<-subset(Original_data_filter,antibioticUse_3month==0)#528, antibioticUse=14


Original_data_filter<-subset(Original_data_filter,Twins==0)#523, n=5
Original_data_filter<-subset(Original_data_filter,ART==0)#513, n=10

Original_data_filter$CHOL_32<-as.numeric(Original_data_filter$CHOL_32)
Original_data_filter$TG_32<-as.numeric(Original_data_filter$TG_32)


Original_data_filter$CHOL_24_High[Original_data_filter$CHOL_24<6.22]<-0
Original_data_filter$CHOL_24_High[Original_data_filter$CHOL_24>=6.22]<-1

Original_data_filter$TG_24_High[Original_data_filter$TG_24<2.26]<-0
Original_data_filter$TG_24_High[Original_data_filter$TG_24>=2.26]<-1


Original_data_filter$CHOL_32_High[Original_data_filter$CHOL_32<6.22]<-0
Original_data_filter$CHOL_32_High[Original_data_filter$CHOL_32>=6.22]<-1

Original_data_filter$TG_32_High[Original_data_filter$TG_32<2.26]<-0
Original_data_filter$TG_32_High[Original_data_filter$TG_32>=2.26]<-1

Original_data_filter$Hyperlipid_24<-ifelse(Original_data_filter$TG_24_High==1&Original_data_filter$CHOL_24_High==1,1,0)
Original_data_filter$Hyperlipid_32<-ifelse(Original_data_filter$TG_32_High==1&Original_data_filter$CHOL_32_High==1,1,0)

###final sample number: 513

# Table S1 ------------------------------------------------------------------

data2<-read.xlsx("data/data_total.xlsx",sheet = 1)
#data2<-Original_data_filter
data2$Parity<-factor(data2$Parity,levels=c(0,1),labels=c('0','≥1'))
data2$preBMI_group<-factor(data2$preBMI_group,levels=c(1,2,3),
                            labels=c('<18.5','18.5-23.9','≥24'))
data2$Educational_level<-factor(data2$Educational_level,levels=c(1,2,3),
                                labels=c('<High School','High school-Bachelor degree','>Bachelor degree'))
data2$Income_group<-factor(data2$Income_group,levels=c(0,1),labels=c('<200,000','≥200,000'))
data2$SEX<-factor(data2$SEX,levels=c(1,2),labels=c('Boys','Girls'))
data2$PTB[data2$Delivery_week<37]<-1
data2$PTB[data2$Delivery_week>=37]<-0

data2$PTB<-factor(data2$PTB,levels=c(0,1),labels=c('No','Yes'))
data2$LBW<-factor(data2$LBW,levels=c(0,1),labels=c('No','Yes'))
data2$preDrinking<-factor(data2$preDrinking,levels=c(0,1),labels=c('No','Yes'))
data2$passsmokHis_1y<-factor(data2$passsmokHis_1y,levels=c(0,1),labels=c('No','Yes'))
data2$SmokingHis<-factor(data2$SmokingHis,levels=c(0,1),labels=c('No','Yes'))

dir.create("Result/linear")
library(compareGroups)
data2$CHOL_32<-as.numeric(data2$CHOL_32)
data2$TG_32<-as.numeric(data2$TG_32)
bas.t <- descrTable(Hyperlipid_24~Age+pre_BMI+preBMI_group+Parity+
                      Educational_level+Income_group+#income+
                      preDrinking+SmokingHis+passsmokHis_1y+
                      CHOL_24+CHOL_32+TG_24+TG_32+
                      SEX+Delivery_week+PTB+Birthweight+LBW,data=data2,show.p.overall = F,show.p.trend = TRUE)
export2csv(bas.t,file = 'Result/linear/table1_1.csv')
tab1.1 <- read.csv('Result/linear/table1_1.csv')


bas.t <- descrTable(Hyperlipid_32~Age+pre_BMI+preBMI_group+Parity+
                      Educational_level+Income_group+#income+
                      preDrinking+SmokingHis+passsmokHis_1y+
                      CHOL_24+CHOL_32+TG_24+TG_32+
                      SEX+Delivery_week+PTB+Birthweight+LBW,data=data2,show.p.overall = F,show.p.trend = TRUE)
export2csv(bas.t,file = 'Result/linear/table1_2.csv')
tab1.2 <- read.csv('Result/linear/table1_2.csv')


bas.t2 <- descrTable(~Age+pre_BMI+preBMI_group+Parity+
                       Educational_level+Income_group+#income+
                       preDrinking+SmokingHis+passsmokHis_1y+
                       CHOL_24+CHOL_32+TG_24+TG_32+
                       SEX+Delivery_week+PTB+Birthweight+LBW,data=data2,show.p.overall = F,show.p.trend = TRUE)
export2csv(bas.t2,file = 'Result/linear/table1_3.csv')
tab1.3 <- read.csv('Result/linear/table1_3.csv')

# bas.t2 <- descrTable(~Age+pre_BMI+preBMI_group+Parity+
#                        Educational_level+Income_group+#income+
#                        preDrinking+SmokingHis+passsmokHis_1y+
#                        CHOL_24+CHOL_32+TG_24+TG_32+
#                        SEX+Delivery_week+PTB+Birthweight+LBW,data=data2,show.p.overall = F,show.p.trend = TRUE)



tab1 <- cbind(tab1.3[,1:2],tab1.1[,2:4],tab1.2[,2:4])
colnames(tab1) <- c('Characteristic','Overall','Hyperlipid-T2','Hyperlipid+T2','P-value','Hyperlipid-T3','Hyperlipid+T3','P-value')
tab1
write.csv(tab1,file = 'Result/linear/table1.csv',row.names = F)



