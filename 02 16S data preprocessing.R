library("phyloseq")
library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("dplyr")        # filter and reformat data frames
library(openxlsx)
library(decontam)
source_dir<-"E:/study2/microbiome/study_new/total/9Data Review"
setwd(paste0(source_dir,"/Data/16S/decontam"))
table_with_taxonomy<-read.table("table-with-taxonomy.txt", header = T,sep = "\t")
tax_mat<- read.table("taxonomy.txt", header = T,sep = "\t",row.names = 1)

rownames(table_with_taxonomy)<-rownames(tax_mat)
otu_mat<- table_with_taxonomy[,2:(dim(table_with_taxonomy)[2]-1)]

samples_df <- read.table("metadata_new.txt", header = T,sep = "\t",row.names = 1)
ID<-rownames(samples_df)
otu_mat<-otu_mat[,ID]

#samples_df<-subset(samples_df,Inclusion=="Yes")
#Transform into matrixes otu and tax tables (sample table can be left as data frame)
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

#Transform to phyloseq objects
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)

carbom <- phyloseq(OTU, TAX, samples)
carbom
head(sample_data(carbom))


#Inspect library sizes
df <- as.data.frame(sample_data(carbom))
df$LibrarySize <- sample_sums(carbom)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
pdf(file = "FigureS1a_library_size.pdf")
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Group2)) + geom_point()+
  theme_bw()
dev.off()

#identify contaminants-Prevalence
sample_data(carbom)$is.neg <- sample_data(carbom)$Group2 == "Control"

#1.threshold_defalut=0.1
contamdf.prev <- isContaminant(carbom, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)
# FALSE  TRUE 
# 39575    116  
head(which(contamdf.prev$contaminant))

#2.identify as contaminants all sequences thare are more prevalent in negative controls than in positive samples
# 
# contamdf.prev05 <- isContaminant(carbom, method="prevalence", neg="is.neg", threshold=0.5)
# table(contamdf.prev05$contaminant)
# # FALSE  TRUE 
# # 39514    177 
# #Let’s take a look at the number of times several of these taxa were observed in negative controls and positive samples.
# 
# Make phyloseq object of presence-absence in negative controls and true samples

carbom.pa <- transform_sample_counts(carbom, function(abund) 1*(abund>0))

carbom.pa.neg <- prune_samples(sample_data(carbom.pa)$Group2 == "Control", carbom.pa)
carbom.pa.pos <- prune_samples(sample_data(carbom.pa)$Group2 == "Sample", carbom.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(carbom.pa.pos), pa.neg=taxa_sums(carbom.pa.neg),
                    contaminant=contamdf.prev$contaminant)
pdf(file = "Figure S1b_presence_absence.pdf")
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")+
  theme_bw()
dev.off()


#1. decomtam ----------------------------------------------------------------
#去除污染ASV
conta<-rownames(contamdf.prev[contamdf.prev$contaminant=="TRUE",])
contam<- otu_mat[conta,]
tax_mat_contam<- tax_mat[conta,]

write.table(contam,file="otu_contam.txt",row.names = T,sep="\t")
write.table(tax_mat_contam,file="tax_contam.txt",row.names = T,sep="\t")

samples_df_filter <- samples_df[samples_df$Group2=="Sample",]
true_ASV<-rownames(contamdf.prev[contamdf.prev$contaminant=="FALSE",])
otu_mat_filter<- otu_mat[true_ASV,rownames(samples_df_filter)]
tax_mat_filter<- tax_mat[true_ASV,]

write.table(otu_mat_filter,file="otu_mat_filter.txt",row.names = T,sep="\t")
write.table(tax_mat_filter,file="tax_mat_filter.txt",row.names = T,sep="\t")
write.table(samples_df_filter,file="samples_df_filter.txt",row.names = T,sep="\t")

# 3.1.抽平 ------------------------------------------------------------------

# otu_mat_filter<-read.table("otu_mat_filter.txt",header = T,sep = "\t",row.names = 1)
# tax_mat_filter<-read.table("tax_mat_filter.txt",header = T,sep = "\t",row.names = 1)

package_list <- c("vegan","dplyr")
# 判断R包加载是否成功来决定是否安装后再加载
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos=site)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}
species<-otu_mat_filter
## 2.1 抽样至最少值或指定值Rarefaction
print(paste0("Samples size are:"))
colSums(species)
min = min(colSums(species))
max = max(colSums(species))
print(paste("Rarefaction depth: ", min, sep = ""))
# 抽样数为零则使用最小值，否则使用指定抽样数
#[1] "Rarefaction depth: 32007"
# print(paste0("Rarefaction depth is ", min))
# vegan::rrarefy抽平至最小值或指定值
seed=1234
set.seed(seed)
otu = vegan::rrarefy(t(species), min)
# print(paste0("All sample rarefaction as following"))
# rowSums(otu)

## 2.3 Alpha diversity
# vegan::estimateR计算obs, chao1和ACE指数
estimateR = t(estimateR(otu))[,c(1,2,4)]
colnames(estimateR) = c("richness", "chao1", "ACE")
# vegan::diversity计算多样性指数shannon, simpson和invsimpson
shannon = diversity(otu, index = "shannon")
simpson = diversity(otu, index = "simpson")
invsimpson = diversity(otu, index = "invsimpson")
# 合并6种指数
alpha_div = cbind(estimateR, shannon, simpson, invsimpson)
print(paste0("Calculate six alpha diversities by estimateR and diversity"))
head(alpha_div, n=1)

# 3. 结果输出
# 3.1 保存抽平的物种表
# 保存一个制表符，解决存在行名时，列名无法对齐的问题
#write.table("#OTUID\t", file=paste(opts$normalize,sep=""),append = F, quote = F, eol = "", row.names = F, col.names = F)
# 保存统计结果，有waring正常
#setwd("E:\\study2\\microbiome\\study_new\\new2\\1OTU\\rarefaction")
suppressWarnings(write.table(t(otu), file="otu_rare.txt", append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T))

tax_mat_filter<-as.data.frame(tax_mat_filter)
tax_mat_filter$Sequencing<-as.character(tax_mat_filter$Sequencing)
tax_mat_filter$taxonomy<-as.character(tax_mat_filter$taxonomy)

combine<-cbind(tax_mat_filter[,1],t(otu),tax_mat_filter[,"taxonomy"])
colnames(combine)[1]<-"#OTUID"
tax_id<-dim(combine)[2]
colnames(combine)[tax_id]<-"taxonomy"
dir.create("new")
write.table(combine, file="new\\table-with-taxonomy-rafe.txt", append = F, sep="\t", quote=F, row.names=F, col.names=T)

# 3.2 保存alpha多样性指数
# 保存一个制表符，解决存在行名时，列名无法对齐的问题
#write.table("Alpha_diversity\t", file="otu_rare_0.001.txt",append = F, quote = F, eol = "", row.names = F, col.names = F)
# 保存统计结果，有waring正常
suppressWarnings(write.table(alpha_div, file="new\\Alpha_diversity.txt", append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T))

#print(paste("Name of rarefaction file ", opts$normalize,  sep = ""))
#print(paste("Output alpha diversity filename ", opts$output, sep = ""))



# 3.2.生成各分类级别的表格 --------------------------------------------------------------
thre<-0
norm = otu/colSums(t(otu),na=T)*100
rowSums(norm)
idx = colMeans(norm) > thre
otu_0 = t(otu)[idx,]
tax_0<-tax_mat_filter[rownames(otu_0),]
write.table(otu_0, file=paste("new\\otu_filter_", thre, ".txt", sep = ""), append = F, sep="\t", quote=F, row.names=T, col.names=T)
write.table(tax_0, file=paste("new\\tax_filter_", thre, ".txt", sep = ""), append = F, sep="\t", quote=F, row.names=T, col.names=T)

otu<-otu_0
tax<-tax_0
tax<-as.data.frame(tax)
tax$Phylum=paste(tax$Kingdom,tax$Phylum,sep = "|")
tax$Class=paste(tax$Phylum,tax$Class,sep = "|")
tax$Order=paste(tax$Class,tax$Order,sep = "|")
tax$Family=paste(tax$Order,tax$Family,sep = "|")
tax$Genus=paste(tax$Family,tax$Genus,sep = "|")
tax$Species=paste(tax$Genus,tax$Species,sep = "|")
# head(tax)

# 按Kingdom合并
grp <- tax[rownames(tax), "Kingdom", drop=F]
merge=cbind(otu, grp)
HA_Kingdom = merge %>% group_by(Kingdom) %>% summarise_all(sum)
colnames(HA_Kingdom)[1]="Kingdom"
write.table(HA_Kingdom, file="class2\\otutab_1Kingdom.txt", append = FALSE, sep="\t", quote=F, row.names=F, col.names=T)
colnames(HA_Kingdom)[1]="Class"

# 按Phylum合并
grp <- tax[rownames(tax), "Phylum", drop=F]
merge=cbind(otu, grp)
HA_Phylum = merge %>% group_by(Phylum) %>% summarise_all(sum)
colnames(HA_Phylum)[1]="Phylum"
write.table(HA_Phylum, file="class2\\otutab_2Phylum.txt", append = FALSE, sep="\t", quote=F, row.names=F, col.names=T)
colnames(HA_Phylum)[1]="Class"

# 按Class合并
grp <- tax[rownames(tax), "Class", drop=F]
merge=cbind(otu, grp)
HA_Class = merge %>% group_by(Class) %>% summarise_all(sum)
colnames(HA_Class)[1]="Class"
write.table(HA_Class, file="class2\\otutab_3Class.txt", append = FALSE, sep="\t", quote=F, row.names=F, col.names=T)
colnames(HA_Class)[1]="Class"

# 按Order合并
grp <- tax[rownames(tax), "Order", drop=F]
merge=cbind(otu, grp)
HA_Order = merge %>% group_by(Order) %>% summarise_all(sum)
colnames(HA_Order)[1]="Order"
write.table(HA_Order, file = "class2\\otutab_4Order.txt", append = FALSE, sep="\t", quote=F, row.names=F, col.names=T)
colnames(HA_Order)[1]="Class"

# 按Family合并
grp <- tax[rownames(tax), "Family", drop=F]
merge=cbind(otu, grp)
HA_Family = merge %>% group_by(Family) %>% summarise_all(sum)
colnames(HA_Family)[1]="Family"
write.table(HA_Family, file="class2\\otutab_5Family.txt", append = FALSE, sep="\t", quote=F, row.names=F, col.names=T)
colnames(HA_Family)[1]="Class"

# 按Genus合并
grp <- tax[rownames(tax), "Genus", drop=F]
merge=cbind(otu, grp)
HA_Genus = merge %>% group_by(Genus) %>% summarise_all(sum)
colnames(HA_Genus)[1]="Genus"
write.table(HA_Genus, "class2\\otutab_6Genus.txt", append = FALSE, sep="\t", quote=F, row.names=F, col.names=T)
colnames(HA_Genus)[1]="Class"

# 按Species合并
grp <- tax[rownames(tax), "Species", drop=F]
merge=cbind(otu, grp)
HA_Species = merge %>% group_by(Species) %>% summarise_all(sum)
colnames(HA_Species)[1]="Species"
write.table(HA_Species, file="class2\\otutab_7Species.txt", append = FALSE, sep="\t", quote=F, row.names=F, col.names=T)
colnames(HA_Species)[1]="Class"

# 合并6个分类级
all = rbind(HA_Kingdom, HA_Phylum, HA_Class, HA_Order, HA_Family, HA_Genus) # , HA_Species
write.table(all, file="class2\\otutab_all.txt", append = FALSE, sep="\t", quote=F, row.names=F, col.names=T)



