source_dir<-"E:/study2/microbiome/study_new/total/9Data Review"
setwd(source_dir)
library(openxlsx)
data = read.xlsx("Data/data_total.xlsx", sheet = 2, rowNames =F)
alpha<-read.xlsx("Data/alpha_diversity.xlsx",sheet = 1)

data_all<-merge(alpha,data,by.x = "ID",by.y = "ID")


# Fig2a -------------------------------------------------------------------

wilcox.test(TG~Time,data,paired=T) #V = 228.5, p-value = 1.312e-08
wilcox.test(CHOL~Time,data,paired=T) #V = 227, p-value = 7.584e-09
wilcox.test(TG~Time,data) #V = 228.5, p-value = 1.312e-08
wilcox.test(CHOL~Time,data) #V = 227, p-value = 7.584e-09

pdf(file="Result/linear/Fig2a.pdf",width = 5,height = 5)
ggplot(data=data,aes(x=Time,y=TG))+
  geom_violin(aes(fill=Time),trim = FALSE)+
  geom_boxplot(alpha=1, outlier.size=0, size=0.9, width=0.2,color="black")+
  #  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  geom_line(aes(group=number) ,size=0.8,colour="#9C9C9C",linetype="dotted")+
  labs(x="", y="TG")+theme_bw()+
  scale_y_continuous(limits=c(0.5, 5))+
  #  scale_color_manual(values = c("#00AFBB", "#E7B800"))+
  scale_fill_manual(values = c("#00AFBB", "#E7B800"))
# scale_color_npg()


ggplot(data=data,aes(x=Time,y=CHOL))+
  geom_violin(aes(fill=Time),trim = FALSE)+
  geom_boxplot(alpha=1, outlier.size=0, size=0.9, width=0.2,color="black")+
  #  geom_jitter( position=position_jitter(0.17), size=1, alpha=0.7)+
  geom_line(aes(group=number) ,size=0.8,colour="#9C9C9C",linetype="dotted")+
  labs(x="", y="CHOL")+theme_bw()+
  scale_y_continuous(limits=c(3, 9))+
  #  scale_color_manual(values = c("#00AFBB", "#E7B800"))+
  scale_fill_manual(values = c("#00AFBB", "#E7B800"))

dev.off()


# Fig2b -------------------------------------------------------
library(ggplot2)
library(ggsci)
library(ggsignif)
library(ggpubr)
cbPalette<-c("#00AFBB", "#E7B800")
def_boxplot_plot<-function(i){
 # dir.create(paste0("alpha/",i))
  ggboxplot(data_all,x="Time",y=i,#add="iitter",
            fill = "Time")+#,notch = TRUE)+
    #   stat_compare_means(comparisons = my_comparisons)+
    stat_compare_means(method = "wilcox.test",paired=TRUE)+
    # stat_compare_means(label = "p.signif", method = "wilcox.test",
    #                    ref.group = "Control") +
    #scale_y_continuous(limits=c(25, 180))+
    theme(legend.position = "none")+
    #  scale_fill_npg()+
    scale_fill_manual(values = cbPalette)+
    theme_bw()+
    labs(x="",y=i)#+
    #facet_wrap(~Time,nrow = 1)
  ggsave(filename = paste0("Result/linear/Fig2b_",i,".pdf"),width=4,height=4)
  
}
# alpha_name<-c("chao1", "simpson", "shannon", "observed_otus", "evenness", 
#               "faith_pd")

# 
#   for (i in alpha_name) {
#     def_boxplot_plot(i)
#   }
def_boxplot_plot("shannon")


#Fig2c ----------------------------------------------------------

library(patchwork)
library(ape) 
library(vegan) 
library(ggsci)
library(ggplot2) 
library(openxlsx)
library(ggpubr)
map<-data_all
sample_ID<-data_all$ID

#shape <- c(16,17) 

# cbPalette<-c("#4DBBD5FF","#E64B35FF")
cbPalette<-c("#00AFBB", "#E7B800")
def_beta_plot<-function(i,j){
 # dir.create(paste0("beta/",i))
  
  data<-read.table(paste0("Data/16S/",j,"_distance_matrix.tsv"),header=TRUE,row.names=1)
  data<-data[sample_ID,sample_ID]
  Group <-  factor(map[,i])
  PCOA <- pcoa(data, correction="none", rn=NULL) 
  result <-PCOA$values[,"Relative_eig"]
  pro1 = as.numeric(sprintf("%.3f",result[1]))*100
  pro2 = as.numeric(sprintf("%.3f",result[2]))*100
  x = PCOA$vectors
  sample_names = rownames(x)
  pc = as.data.frame(PCOA$vectors)
  pc$names = sample_names
  legend_title = ""
  group = Group
  pc$group = group
  xlab=paste("PCoA1(",pro1,"%)",sep="") 
  ylab=paste("PCoA2(",pro2,"%)",sep="")
#  cbPalette<-c("#4DBBD5FF","#E64B35FF")
  
  # pdf(file="beta/weighted/Weighted-unifrac-pcoa_boxplot.pdf",width = 4,height=6)
  
  # p2<-ggplot(data=PC,aes(x=Group,y=Axis.1,fill=Group))+
  #   geom_boxplot(alpha=1, outlier.size=0, size=0.9, width=0.6)+
  p2<-  ggboxplot(pc, x="group", y="Axis.1", fill = "group", 
                  palette = "jco")+
    stat_compare_means(aes(group=group),label = "p.format")+
    labs(x=" ", y="Axis.1")+theme_bw()+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position = "none")+
#    scale_y_continuous(limits=c(-0.4, 0.4))+
    theme(axis.ticks.length = unit(0.4,"lines"), 
          axis.ticks = element_line(color='black'),
          axis.line = element_line(colour = "black"), 
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_text(colour='black',size=10,face = "plain"),
          axis.text.x=element_blank(),
          legend.position = "none")+
    scale_fill_manual(values = cbPalette)+coord_flip()
  
  p3<-ggboxplot(pc, x="group", y="Axis.2", fill = "group", 
                palette = "jco")+
    stat_compare_means(aes(group=group),label = "p.format")+
    labs(x=" ", y="Axis.2")+theme_bw()+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position = "none")+
    scale_fill_manual(values = cbPalette)+
    #scale_y_continuous(limits=c(-0.4, 0.4))+
    theme(axis.ticks.length = unit(0.4,"lines"), 
          axis.ticks = element_line(color='black'),
          axis.line = element_line(colour = "black"), 
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_text(colour='black',size=10,angle = 0,
                                   vjust = 1,hjust = 0.5,face = "plain"),
          axis.text.y=element_blank(),
          legend.position = "none")
  
 # cbPalette<-c("#4DBBD5FF","#E64B35FF")
  
  p1<-ggplot(pc,aes(Axis.1,Axis.2))+geom_point(size=2,aes(color=group))+
    labs(x=xlab,y=ylab,color=legend_title,shape=legend_title)+
    geom_hline(yintercept=0,linetype=4,color="grey")+geom_vline(xintercept=0,linetype=4,color="grey")+
    theme_bw()+
    stat_ellipse(aes(Axis.1,Axis.2,color=group))+#theme(legend.title=element_blank())+
  #  scale_x_continuous(limits=c(-0.5, 0.5))+scale_y_continuous(limits=c(-0.3, 0.3))+
    
    theme(panel.background = element_rect(fill='white', colour='black'),
          axis.title.x=element_text(colour='black', size=12),
          axis.title.y=element_text(colour='black', size=12),
          axis.text=element_text(colour='black',size=12),
          legend.title=element_blank(),
          legend.key.height=unit(0.6,"cm"),
          legend.position = c(0.75, 0.95),legend.direction = "horizontal")+
    scale_color_manual(values = cbPalette)
  
  
  otu.adonis=adonis(as.formula(paste0("data~",i)),data=map)
  
  p4 <- ggplot(pc,
               aes(Axis.1, Axis.2))+
    geom_text(aes(x = -0.5,
                  y = 0.6,
                  label = paste("PERMANOVA:\ndf = ",
                                otu.adonis$aov.tab$Df[1],#"\nR2 = ",
                                # round(otu.adonis$aov.tab$R2[1],4),
                                "\np-value = ",
                                otu.adonis$aov.tab$`Pr(>F)`[1],
                                sep = "")),size = 4) +theme_bw() +
    xlab(NULL) + ylab(NULL) +
    theme(panel.grid=element_blank(), 
          axis.title = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank())
  
  
  #pdf(file = paste0("beta/",i,"/dyslipimedia_T2.pdf"),width = 12,height = 12)
  p2+p4+p1+p3 + 
    plot_layout(heights = c(1,4),widths = c(4,1),ncol = 2,nrow = 2)
  # dev.off()
  ggsave(file = paste0("Result/linear/Fig2c_",i,"_",j,".pdf"),width = 8,height = 8)
}

#outcome<-c("unweighted_unifrac","weighted_unifrac","bray_curtis")
outcome<-"unweighted_unifrac"
for (i in c("Time")) {
  for (j in outcome) {
    def_beta_plot(i,j)
  }
}


#Fig2d----------------------------------------------------------------

data_plot<-read.xlsx("Data/16S/Phylum_plot.xlsx",sheet = 2)
library(reshape2)
data_plot<-melt(data_plot,id.vars = c("Phylum"))
data_plot$Phylum<-factor(data_plot$Phylum,
                         levels=rev(unique(data_plot$Phylum)))

pdf(file="Result/linear/Fig2d.pdf")
ggplot(data_plot,aes(x=variable,y=value,fill=Phylum))+geom_bar(stat="identity")+scale_fill_npg()+labs(y="Relative abundance (%)")
dev.off()

