setwd('E:/study2/microbiome/study_new/total/9Data Review')
library(microbiome)
library(DirichletMultinomial)
library(reshape2)
library(magrittr)
library(dplyr)
library(phyloseq)
library(openxlsx)
library(ConsensusClusterPlus)
library(openxlsx)
library(ggplot2)
library(ggpubr)
#setwd("E:/study2/microbiome/study_new/total")
metadata = read.xlsx("Data/data_total.xlsx", sheet=2,rowNames = F)

metadata2 = read.xlsx("Data/data_total.xlsx", sheet=1,rowNames = F)
metadata2<-subset(metadata2,dyslipimedia_24=="Control")

ID_T2<-metadata2$ID_24
ID_T3<-metadata2$ID_32

#二、 genus_level -------------------------------------------------------------

genus_filter<-read.csv("Data/16S/genus_filter_RA.csv",header = T,row.names = 1)
genus_name<-colnames(genus_filter)

tax<-read.xlsx("Data/16S/genus_filter.xlsx",sheet = 2)
tax<-tax[,1:6]
# for (i in 1:ncol(genus_filter)){
#   genus_filter[,i][genus_filter[,i]==0] <-1/20993*100
# }
# 
# genus_filter_log<-log10(genus_filter)
genus_filter$ID<-rownames(genus_filter)
data_genus_RA<-merge(genus_filter,metadata,by.x = "ID",by.y = "ID")

data_24<-data_genus_RA[data_genus_RA$Time=="24week",]
data_32<-data_genus_RA[data_genus_RA$Time=="32week",]

microb_24<-data_24[genus_name]
rownames(microb_24)<-data_24$ID
microb_24<-as.data.frame(t(microb_24))

microb_32<-data_32[genus_name]
rownames(microb_32)<-data_32$ID
microb_32<-as.data.frame(t(microb_32))


rownames(metadata)<-metadata$ID
metadata_24<-subset(metadata,Time=="24week")
metadata_32<-subset(metadata,Time=="32week")


# total ------------------------------------------------------------------
dir.create("Result/DMM")
tax_total<-tax
rownames(tax_total)<-tax_total$Genus
microb<-data_genus_RA[genus_name]
rownames(microb)<-data_genus_RA$ID
microb<-as.data.frame(t(microb))


microb<-as.matrix(microb)
tax_total<-as.matrix(tax_total)
OTU = otu_table(microb, taxa_are_rows = TRUE)
TAX = tax_table(tax_total)
samples = sample_data(metadata)

carbom <- phyloseq(OTU, TAX, samples)
carbom
head(sample_data(carbom))


# Load example data
# data(dietswap)
# pseq <- dietswap

# To speed up, only consider the core taxa
# that are prevalent at 0.1% relative abundance in 50% of the samples
# (note that this is not strictly correct as information is
# being discarded; one alternative would be to aggregate rare taxa)
# pseq.comp <- microbiome::transform(pseq, "compositional")
# taxa <- core_members(pseq.comp, detection = 0.1/100, prevalence = 50/100)
# pseq <- prune_taxa(taxa, pseq)

# Pick the OTU count matrix
# and convert it into samples x taxa format
dat <- abundances(carbom)
count <- as.matrix(t(dat))
#Fit the DMM model. Let us set the maximum allowed number of community types to 3 to speed up the example.

fit <- lapply(1:10, dmn, count = count, verbose=TRUE)

#Check model fit with different number of mixture components using standard information criteria

lplc <- sapply(fit, laplace) # AIC / BIC / Laplace
aic  <- sapply(fit, AIC) # AIC / BIC / Laplace
bic  <- sapply(fit, BIC) # AIC / BIC / Laplace

p<-plot(lplc, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit")
# lines(aic, type="b", lty = 2)
# lines(bic, type="b", lty = 3)
ggsave(p,filename="Result/DMM/DMM_cluster_lplc.pdf",width=8,height=6)
#Pick the optimal model

best <- fit[[which.min(unlist(lplc))]]
#Mixture parameters pi and theta

mixturewt(best)
# pi     theta
# 1 0.3671777 21.364616
# 2 0.3071238  9.772219
# 3 0.1728329 22.422971
# 4 0.1528656 15.575160

#Sample-component assignments

ass <- apply(mixture(best), 1, which.max)
#Contribution of each taxonomic group to each component

for (k in seq(ncol(fitted(best)))) {
  d <- melt(fitted(best))
  write.csv(d,file=paste0("Result/DMM/Drivers_",k,".csv"))
  colnames(d) <- c("OTU", "cluster", "value")
  d <- subset(d, cluster == k) %>%
    # Arrange OTUs by assignment strength
    arrange(value) %>%
    mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
    # Only show the most important drivers
    filter(abs(value) > quantile(abs(value), 0.8))     
  
  p <- ggplot(d, aes(x = OTU, y = value)) +
    geom_bar(stat = "identity") +
    theme_bw()+
    theme(axis.text.y = element_text(face = "italic"))+
    coord_flip() +
    labs(title = paste("Top drivers: community type", k))
  print(p)
  ggsave(p,filename = paste0("Result/DMM/top_drivers_",k,".pdf"),width = 7,height = 8)
}


###new driver plot

dat<-read.csv("Result/DMM/Drivers_1.csv",header = T,sep = ",")
dat<-dat[,-1]
#write.csv(d,file=paste0("Result/DMM/Drivers_",k,".csv"))
colnames(dat) <- c("OTU", "cluster", "value")

for (k in 1:4) {
  d <- subset(dat, cluster == k) %>%
    # Arrange OTUs by assignment strength
    arrange(value) %>%
    mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
    # Only show the most important drivers
    filter(abs(value) > quantile(abs(value), 0.8))     
  
  p <- ggplot(d, aes(x = OTU, y = value)) +
    geom_bar(stat = "identity") +
    theme_bw()+
    theme(axis.text.y = element_text(face = "italic"))+
    coord_flip() +
    labs(title = paste("Top drivers: community type", k))
  print(p)
  ggsave(p,filename = paste0("Result/DMM/top25_drivers_",k,".pdf"),width = 7,height = 8)
}


clusters=mixture(best,assign = TRUE)
clusters=data.frame(sample=names(clusters),cluster=clusters)
write.table(clusters,file = "Result/DMM/DMM_clusters_group.txt",sep = "\t",quote = F,row.names = F)
library(vegan)
library(ade4)
data.dist<-vegdist(count,method="bray")
obs.pcoa<-dudi.pco(data.dist,scannf=F,nf=3)

pdf(file = "Result/DMM/genus_DMM_PCOA.pdf",width = 6,height = 6)
s.class(obs.pcoa$li,fac = as.factor(clusters$cluster),
        grid = F,sub = "Principal coordinate analysis",
        # col = c(1,2,3,4)
        col = c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF"))

dev.off()


count_T2<-count[metadata_24$ID,]
count_T2<-as.matrix(count_T2)
data.dist<-vegdist(count_T2,method="bray")
obs.pcoa<-dudi.pco(data.dist,scannf=F,nf=3)

clusters_T2<-subset(clusters,sample%in%metadata_24$ID)
pdf(file = "Result/DMM/genus_DMM_PCOA_24week.pdf",width = 6,height = 6)
s.class(obs.pcoa$li,fac = as.factor(clusters_T2$cluster),
        grid = F,sub = "Principal coordinate analysis",
        # col = c(1,2,3,4)
        col = c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF"))

dev.off()



count_T3<-count[metadata_32$ID,]
count_T3<-as.matrix(count_T3)
data.dist<-vegdist(count_T3,method="bray")
obs.pcoa<-dudi.pco(data.dist,scannf=F,nf=3)

clusters_T3<-subset(clusters,sample%in%metadata_32$ID)
pdf(file = "Result/DMM/genus_DMM_PCOA_32week.pdf",width = 6,height = 6)
s.class(obs.pcoa$li,fac = as.factor(clusters_T3$cluster),
        grid = F,sub = "Principal coordinate analysis",
        # col = c(1,2,3,4)
        col = c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF"))

dev.off()

# heatmap -----------------------------------------------------------------

dat_heatmap<-read.xlsx("Result/DMM/change of cluster.xlsx",sheet=2)
  
genera<-data_genus_RA[genus_name]
genera<-as.data.frame(t(genera))
top25 <- head(names(rev(sort(rowSums(genera)))), 25)

colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(n=7, name="RdYlBu")), bias=3)(100)  #use to set colour (name =) and scale (bias =)

anno<-subset(dat_heatmap,Var1%in%top25)
write.csv(anno,file = "Result/DMM//top25_anno.csv")
ggplot(anno, aes(Var2, Var1)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradientn(colours = 
                         colors, 
                       name ="Relative abundance", 
                       na.value = "white") +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.position = "left",
        plot.title = element_blank(),
        axis.title=element_blank(),
        axis.text.y=element_text(size=8,face = "italic"),
        axis.text.x = element_blank())
ggsave(filename = "Result/DMM/heatmap of cluster.pdf",width = 6,height = 8)


ggplot(anno, aes(Var2, Var1)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradientn(colours = 
                         colors, 
                       name ="Relative abundance", 
                       na.value = "white") +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.position = "left",
        plot.title = element_blank(),
        axis.title=element_blank(),
        axis.text.x=element_text(size=8,face = "italic",angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_blank())+
  coord_flip()
ggsave(filename = "Result/DMM/heatmap of cluster_rotate.pdf",width =8,height = 4)

# #anno_new<-anno
# anno_new<-anno$phylum[1:25]
# names(anno_new)<-unique(anno$Var1)[1:25]
# mat <- genera[top25,]
# mat <- t(apply(mat, 1L, scales::rescale))   #uncomment to normalize each taxa (row)
# #anno_new<-anno_new[,-2]
# pheatmap::pheatmap(
#   mat            = mat, 
#   color          = colors,   #uncheck if setting the colour scale manual
# #  annotation_col = anno_new, 
#   show_colnames  = FALSE,
#   cluster_rows   = FALSE,
#   cluster_cols   = FALSE,
# #  gaps_col       = cumsum(unname(table(anno[[CLUSTER_COLUMN]]))),
#   labels_row     = sub("^.*__", "", top25)
# )

# trasition plot ----------------------------------------------------------

library(grid)
library(Gmisc)
library(magrittr)
library(RColorBrewer)
library(networkD3)
dat_DMM_group<-read.xlsx("Result/DMM/change of cluster.xlsx",sheet=1)
dat_DMM_group<-dat_DMM_group[,c("sample","cluster")]
metadata = read.xlsx("data/data_total.xlsx", sheet=2,rowNames = F)

metadata_all<-merge(dat_DMM_group,metadata,by.x="sample",by.y="ID")

kruskal.test(metadata_all$TG~metadata_all$cluster)

metadata_T2<-subset(metadata_all,Time=="24week")
metadata_T3<-subset(metadata_all,Time=="32week")
cluster_T2<-metadata_T2$cluster
cluster_T3<-metadata_T3$cluster

dat_cluster<-data.frame(cluster_T2,cluster_T3)
t<-table(dat_cluster$cluster_T2,dat_cluster$cluster_T3)
t<-as.data.frame(t)
colnames(t) <- c("source", "target", "value")
t$target <- paste(t$target, " ", sep="")

nodes <- data.frame(name=c(as.character(t$source), as.character(t$target)) %>% unique())

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
t$IDsource=match(t$source, nodes$name)-1 
t$IDtarget=match(t$target, nodes$name)-1
# prepare colour scale
ColourScal ='d3.scaleOrdinal() .range(["#FDE725FF","#B4DE2CFF","#6DCD59FF","#35B779FF","#1F9E89FF","#26828EFF","#31688EFF","#3E4A89FF","#482878FF","#440154FF"])'

# Make the Network
p<-sankeyNetwork(Links = t, Nodes = nodes,
                 Source = "IDsource", Target = "IDtarget",
                 Value = "value", NodeID = "name", 
                 sinksRight=FALSE, 
                 # colourScale=ColourScal, 
                 nodeWidth=40, fontSize=13, nodePadding=20)

p2<-htmlwidgets::onRender(p, '
  function(el) { 
    var nodeWidth = this.sankey.nodeWidth();
    var links = this.sankey.links();
        
    links.forEach((d, i) => {
      var startX = d.source.x + nodeWidth;
      var endX = d.target.x;
      
      var startY = d.source.y + d.sy + d.dy / 2;
      var endY = d.target.y + d.ty + d.dy / 2;
      
      d3.select(el).select("svg g")
        .append("text")
        .attr("text-anchor", "middle")
        .attr("alignment-baseline", "middle")
        .attr("x", startX + ((endX - startX) / 2))
        .attr("y", startY + ((endY - startY) / 2))
        .text(d.value);
    })
  }
')
p
ggsave(p,filename = "Result/DMM/sankey.pdf",height = 8,width = 4,device = "pdf")
write.csv(t,file = "Result/DMM/sankey_plot_data.csv",row.names = F)


# diversity ---------------------------------------------------------
alpha<-read.xlsx("Data/alpha_diversity.xlsx",sheet = 1)

metadata_all<-merge(metadata_all,alpha,by.x="sample",by.y="ID")
# 一、alpha-diversity -------------------------------------------------------
library(ggplot2)
library(ggsci)
library(ggsignif)

# 2、boxplot-time  ---------------------------------------------------------

# 配色方案1 -------------------------------------------------------------------
dir.create("Result/DMM/cluster")

my_comparisons <- list( c(1, 4), c(1,3),
                        c(2,4) )

def_boxplot_plot<-function(i){
  ggboxplot(metadata_all,x="cluster",y=i,add="jitter",
            fill = "cluster",notch = TRUE)+
    stat_compare_means(comparisons = my_comparisons)+
    stat_compare_means()+
    # stat_compare_means(label = "p.signif", method = "wilcox.test",
    #                    ref.group = "Control") +
    #scale_y_continuous(limits=c(25, 180))+
    #scale_x_discrete(limits=c(2,4,1,3))+
    theme(legend.position = "none")+
    scale_fill_npg()+
    labs(x="",y=i)+
    facet_wrap(~Time,nrow = 1)
  ggsave(filename = paste0("Result/DMM/cluster/",i,"_alpha.pdf"),width=15,height=10)
  
}
alpha_name<-c("chao1", "simpson", "shannon", "observed_otus", "evenness", 
              "faith_pd")
for (i in alpha_name) {
  def_boxplot_plot(i)
}



# 2、boxplot-TG  ---------------------------------------------------------

# 配色方案1 -------------------------------------------------------------------
dir.create("Result/DMM/cluster")
# 
# my_comparisons <- list( c(1, 4), c(1,3),
#                         c(2,4) )

def_boxplot_plot<-function(i){
  ggboxplot(metadata_all,x="cluster",y=i,add="jitter",
            fill = "cluster",notch = TRUE)+
#    stat_compare_means(comparisons = my_comparisons)+
    stat_compare_means()+
    # stat_compare_means(label = "p.signif", method = "wilcox.test",
    #                    ref.group = "Control") +
    #scale_y_continuous(limits=c(25, 180))+
    #scale_x_discrete(limits=c(2,4,1,3))+
    theme(legend.position = "none")+
    scale_fill_npg()+
    labs(x="",y=i)+
    facet_wrap(~Time,nrow = 1)
  ggsave(filename = paste0("Result/DMM/cluster/",i,".pdf"),width=15,height=10)
  
}
alpha_name<-c("TG", "CHOL")
for (i in alpha_name) {
  def_boxplot_plot(i)
}


# 3、boxplot  ---------------------------------------------------------

# 配色方案1 -------------------------------------------------------------------
dir.create("Result/DMM/cluster")

my_comparisons <- list( c(1, 3), c(1,4),
                        c(2,4) )

def_boxplot_plot<-function(i){
  ggboxplot(metadata_all,x="cluster",y=i,#add="jitter",
            fill = "cluster",notch = TRUE)+
    stat_compare_means(comparisons = my_comparisons)+
    stat_compare_means()+
    # stat_compare_means(label = "p.signif", method = "wilcox.test",
    #                    ref.group = "Control") +
    #scale_y_continuous(limits=c(25, 180))+
    #scale_x_discrete(limits=c(2,4,1,3))+
    theme(legend.position = "none")+
    scale_fill_npg()+
    labs(x="",y=i)#+
   # facet_wrap(~Time,nrow = 1)
  ggsave(filename = paste0("Result/DMM/cluster/",i,"_all.pdf"),width=8,height=8)
  
}
alpha_name<-c("chao1", "simpson", "shannon", "observed_otus", "evenness", 
              "faith_pd")
for (i in alpha_name) {
  def_boxplot_plot(i)
}



# 3、boxplot-TG  ---------------------------------------------------------

# 配色方案1 -------------------------------------------------------------------
dir.create("Result/DMM/cluster")
# 
my_comparisons <- list( c(1, 3), c(1,4),
                        c(2,4) )

def_boxplot_plot<-function(i){
  ggboxplot(metadata_all,x="cluster",y=i,#add="jitter",
            fill = "cluster",notch = TRUE)+
    stat_compare_means(comparisons = my_comparisons)+
    stat_compare_means()+
    # stat_compare_means(label = "p.signif", method = "wilcox.test",
    #                    ref.group = "Control") +
    #scale_y_continuous(limits=c(25, 180))+
    #scale_x_discrete(limits=c(2,4,1,3))+
    theme(legend.position = "none")+
    scale_fill_npg()+
    labs(x="",y=i)#+
  #  facet_wrap(~Time,nrow = 1)
  ggsave(filename = paste0("Result/DMM/cluster/",i,"_all.pdf"),width=8,height=8)
  
}
alpha_name<-c("TG", "CHOL")
for (i in alpha_name) {
  def_boxplot_plot(i)
}




