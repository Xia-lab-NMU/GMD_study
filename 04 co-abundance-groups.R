setwd('E:/study2/microbiome/study_new/total/9Data Review')
library(openxlsx)
library(ConsensusClusterPlus)
wgcnadata<-read.csv("Data/16S/data_combine_genus2.csv",header = T,sep = ",")
#wgcnadata<-merge(metadata,genus_filter,by.x = "ID",by.y = "ID")
#rownames(wgcnadata)<-wgcnadata$ID
genus_name<-colnames(wgcnadata[,59:ncol(wgcnadata)])

# 
# wgcnadata_24<-wgcnadata[wgcnadata$Time=="24W",]
# wgcnadata_32<-wgcnadata[wgcnadata$Time=="32W",]
# combine -----------------------------------------------------------------

#co-abundance groups
library(made4)
library(vegan)
library(SpiecEasi)
library(bioDist)
microbiome<-wgcnadata[,genus_name]

#delete genus with standard deviation ==0
microbiome<-microbiome[,apply(microbiome, 2, function(x) sd(x)!=0)]

c<-cor(microbiome, method = "kendall")
c_dist<-as.dist(1-cor(t(c), method="spearman"))
clu<-hclust(c_dist, method="ward.D2")
plot(clu)
abline(h = 1.4, col = "red");

tauD = tau.dist(c)
tauD_new<-as.matrix(tauD)

clust = cutree(clu, k=5)
table(clust)
clust<-as.data.frame(clust)
clust$clust<-factor(clust$clust)
#a<-adonis(tauD_new ~clust,data=clust)
a<-adonis(tauD_new ~clust,data=clust)

#define pairwise_adonis function
pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='bonferroni')
{
  library(vegan)
  
  co = combn(unique(as.character(factors)),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  
  for(elem in 1:ncol(co)){
    if(sim.function == 'daisy'){
      library(cluster); x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
    } else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
    
    ad = adonis(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])] );
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.05] <-'.'
  sig[p.adjusted <= 0.01] <-'*'
  sig[p.adjusted <= 0.001] <-'**'
  sig[p.adjusted <= 0.0001] <-'***'
  
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted,sig)
  print("Signif. codes:  0 ¡®***¡¯ 0.001 ¡®**¡¯ 0.01 ¡®*¡¯ 0.05 ¡®.¡¯ 0.1 ¡® ¡¯ 1")
  return(pairw.res)
  
} 
pair_test<-pairwise.adonis(tauD_new,clust$clust)
write.csv(pair_test,file = "CAG/combine/pair_test_combine.csv")

colors = seq(-0.4,0.4,length=100)
my_palette <- colorRampPalette(c("#4DBBD5","white","#E64B35"))(n = 99)

#plot_color = c('#FF9900','#009900','#FC3200','#0066CC','#FFFF00','#FF9999','#9999FF','#E6399B')[clust$clust]
plot_color = c('#FF9900','#009900','#FC3200','#0066CC','#FFFF00','#FF9999')[clust$clust]

#plot correlation heatmap
pdf(file = "Result/linear/FigS2_heatmap_combine.pdf")
heatmap.2(c, col=my_palette, 
          breaks=colors, density.info="none", trace="none", 
          ColSideColors = plot_color,
          RowSideColors = plot_color,
          distfun   = function(x) as.dist(1-cor(t(x), method="spearman")), 
          hclustfun = function(x) hclust(x, method="ward.D2"),
          symm=F,symkey=F,symbreaks=T, scale="none")
dev.off()

data_new<-cbind(t(microbiome),clust)
write.csv(data_new,file = "Data/16S/CAG/data_CAG.csv")
#²»ÓÃ£¬Êµ¼ÊÐèÒªÁ½×é½øÐÐÅä¶Ô¼ìÑé
result<-data.frame()
for (i in 3:15)
{clust = cutree(clu, k=i)
table(clust)
clust<-as.data.frame(clust)
clust$clust<-factor(clust$clust)
a<-adonis(tauD_new ~clust,data=clust)
pairwise.adonis(tauD_new,clust$clust)
r2<-a$aov.tab[1,5]
p<-a$aov.tab[1,6]
result_i<-c(i,r2,p)
result_i<-as.data.frame(result_i)
result_i<-as.data.frame(t(result_i))
result<-rbind(result,result_i)}
write.csv(result,file = "Data/16S/CAG/result_kendall_clk_combine.csv")
#h=1.25£¬get the biggest R2
write.csv(c,file = "Data/16S/CAG/result_kendall_correlation.csv")

#correlation fdr
library(psych)
res2<-corr.test(microbiome,use="complete",method="kendall",adjust="fdr")
r<-res2$r
p<-res2$p
p1<-p
for (i in 1:dim(p1)[1]){
  p1[i,i]<-NA
}

flattenCorrMatrix <- function(cormat, pmat){
  ut <- upper.tri(cormat) 
  data.frame( row = rownames(cormat)[row(cormat)[ut]], column = rownames(cormat)[col(cormat)[ut]], cor =(cormat)[ut], p = pmat[ut] )
}
result<-flattenCorrMatrix(r,p1)

flattenCorrMatrix2 <- function(cormat, pmat){
  ut <- lower.tri(cormat) 
  data.frame( row = rownames(cormat)[row(cormat)[ut]], column = rownames(cormat)[col(cormat)[ut]], cor =(cormat)[ut], p = pmat[ut] )
}
result2<-flattenCorrMatrix2(r,p1)

result$index<-paste(result$row,result$column,sep = "_")
result2$index<-paste(result2$colum,result2$row,sep = "_")
result_combine<-merge(result2,result[,c("index","p")],by.x = "index",by.y = "index")

names(result_combine)<-c("index","row","column","cor","pvalue","qvalue")
write.csv(result_combine,file="Data/16S/CAG/kendall.csv")

write.csv(p1,"Data/16S/CAG/pvalue of correlation.csv")
write.csv(r,"Data/16S/CAG/correlation.csv")

# CAG combine ----------------------------------------------------------
library(pheatmap)
Data_mean = aggregate(data_new[,-dim(data_new)[2]], by=data_new[dim(data_new)[2]], FUN=sum)
# 各组求和 
Data4Pic = as.data.frame(do.call(rbind, Data_mean)[,-dim(data_new)[2]])
names(Data4Pic)<-paste("CAG",seq(1,5),sep = "")
Data4Pic<-Data4Pic[-1,]
write.csv(Data4Pic,file = "Data/16S/CAG/data_CAG_combine.csv")
# data_CAG_combine<-Data4Pic
# #colors = seq(-0.4,0.4,length=100)
# my_palette=colorRampPalette(c("#4DBBD5","white","#E64B35"))(6)
# #correlation
# c<-cor(data_CAG_combine, method = "kendall")
#plot_color = c('#FF9900','#009900','#FC3200','#0066CC','#FFFF00','#FF9999','#9999FF','#E6399B')[clust$clust]
#plot correlation heatmap
# breaksList = seq(-0.15, 0.15, by = 0.05)
# pdf(file = "Data/16S/CAG/heatmap_CAG_combine.pdf")
# # heatmap.2(c, col=my_palette,
# #           breaks=colors, density.info="none", trace="none",
# #           ColSideColors = plot_color,
# #           RowSideColors = plot_color,
# #           distfun   = function(x) as.dist(1-cor(t(x), method="spearman")),
# #           hclustfun = function(x) hclust(x, method="ward.D2"),
# #           symm=F,symkey=F,symbreaks=T, scale="none")
# 
# 
# pheatmap(c, scale="none",cluster_row = T,cluster_cols = T,
#          # display_numbers = sig.mat,
#          color=colorRampPalette(c("#4DBBD5","white","#E64B35"))(6),
#          border=FALSE)
# dev.off()


#exclude CAG1
data_CAG<-read.csv("Data/16S/CAG/data_CAG.csv",header = T,sep = ",")
genus_name<-data_CAG$X
data_CAG_2<-subset(data_CAG,clust!=1)
genus_n<-as.character(data_CAG_2$X)


correlation_genera<-read.csv("Data/16S/CAG/correlation.csv",header = T,sep = ",",row.names = 1)
pvalue_genera<-read.csv("Data/16S/CAG/pvalue of correlation.csv",header = T,sep = ",",row.names = 1)
spearman<-read.csv("Data/16S/CAG/spearman.csv")
correlation_genera<-correlation_genera[genus_n,genus_n]
pvalue_genera<-pvalue_genera[genus_n,genus_n]
spearman<-subset(spearman,row %in%genus_n &column%in%genus_n)

sum_genus<-apply(data_CAG_2[2:(ncol(data_CAG)-1)],1,mean)
sum_genera<-data.frame(genus_n,sum_genus)

mean<-apply(data_CAG[,2:(ncol(data_CAG)-1)],1,mean)
mean_genera<-data.frame(genus_name,mean)

dir.create("Data/16S/CAG/cytoscape")
write.csv(correlation_genera,"Data/16S/CAG/cytoscape/correlation_genera.csv")
write.csv(pvalue_genera,"Data/16S/CAG/cytoscape/pvalue_genera.csv")
write.csv(spearman,"Data/16S/CAG/cytoscape/spearman.csv")
write.csv(sum_genera,"Data/16S/CAG/cytoscape/sum_genera.csv")
write.csv(mean_genera,"Data/16S/CAG/cytoscape/mean_genera_all.csv")




