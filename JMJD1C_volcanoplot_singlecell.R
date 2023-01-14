setwd("D:/Jmjd1c_Treg_Tumor/single_cell_RNAseq/201692_20210126_1_foxp3jmjd1c_mc205_tumour_wt_ko_treg/Seurat2_result/")
diff_gene <- read.table("pbmc.markers_forvolcanoplot.txt",sep = "\t", header=TRUE, row.names=1)
diff_gene=as.data.frame(diff_gene)
gene_list=diff_gene[,c("avg_log2FC","p_val_adj")]
colnames(gene_list)=c("logFC","padj")
gene_list$threshold = as.factor(abs(gene_list$logFC) > 0.25 & gene_list$padj < 0.00005)
colored_point<-gene_list[gene_list$threshold == "TRUE",]
gene_list$threshold<-as.character(gene_list$threshold)

gene_list$threshold[which(gene_list$threshold == "TRUE" & gene_list$logFC >0.25)] <- "UP"
gene_list$threshold[which(gene_list$threshold == "TRUE" & gene_list$logFC < -0.25)] <- "Down"
colnames(gene_list)[3]<-"Significant"
gene_list$Significant[which(gene_list$Significant == "FALSE")]<-"Not Sig"
Mycolors<-c("Black","Gray","Black")
library("ggplot2")
pdf("vocano.pdf")

g = ggplot(data=gene_list, aes(x=logFC, y=-log10(padj),color=Significant)) + geom_point(alpha=0.4, size=1.75)  + xlim(c(-2, 2)) + ylim(c(0, 300)) +xlab("log2 fold change") + ylab("-log10 p-value") + theme_set(theme_bw()) + theme(panel.grid.major=element_line(colour=NA)) + scale_color_manual(values = Mycolors)
print(g)
dev.off()
