library(dplyr)
RNA<-read.table("D:/Jmjd1c_Treg_Tumor/single_cell_RNAseq/201692_20210126_1_foxp3jmjd1c_mc205_tumour_wt_ko_treg/Seurat2_result/pbmc.markers_forvolcanoplot_UP.txt",sep = "\t",header = TRUE,row.names = 1)
UP<-read.table("D:/stat3/20211215/20211215lxh_stat3_cuttag/homer/KO_vs_WT_tag10_F0.00001_P0.05_L0.00001_LP0.05_C5_withannotation_1.txt",sep = "\t",header = TRUE,row.names = 1)
Down<-read.table("D:/stat3/20211215/20211215lxh_stat3_cuttag/homer/WT_vs_KO_tag10_F0.00001_P0.05_L0.00001_LP0.05_C5_1.txt",sep = "\t",header = TRUE,row.names = 1)
colnames(UP)[29:33]<-c("F","P","L","LP","C")
colnames(Down)[10:14]<-c("F","P","L","LP","C")
k<-data.frame(0,0,0,0,0,0,0,0)
colnames(k)<-c("a","b","c","d","f","g","h","i")
for (a in c(1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5)) {
  for (b in c(5e-2,1e-2,5e-3,1e-3,5e-4,1e-4,5e-5,1e-5,5e-6,1e-6,5e-7,1e-7,5e-8,1e-8,5e-9,1e-9,5e-10,1e-10,
              5e-11,1e-11,5e-12,1e-12,5e-13,1e-13,5e-14,1e-14,5e-15,1e-15,5e-16,1e-16,5e-17,1e-17,5e-18,1e-18,
              5e-19,1e-19,5e-20,1e-20,5e-21,1e-21,5e-22,1e-22,5e-23,1e-23,5e-24,1e-24,5e-25,1e-25,5e-26,1e-26,
              5e-27,1e-27,5e-28,1e-28,5e-29,1e-29,5e-30,1e-30,5e-31,1e-31,5e-32,1e-32,5e-33,1e-33,5e-34,1e-34,5e-35,1e-35,
              5e-36,1e-36,5e-37,1e-37,5e-38,1e-38,5e-39,1e-39,5e-40,1e-40)) {
    for (c in c(1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5)) {
      for (d in c(5e-2,1e-2,5e-3,1e-3,5e-4,1e-4,5e-5,1e-5,5e-6,1e-6,5e-7,1e-7,5e-8,1e-8,5e-9,1e-9,5e-10,1e-10,
                  5e-11,1e-11,5e-12,1e-12,5e-13,1e-13,5e-14,1e-14,5e-15,1e-15,5e-16,1e-16,5e-17,1e-17,5e-18,1e-18,
                  5e-19,1e-19,5e-20,1e-20,5e-21,1e-21,5e-22,1e-22,5e-23,1e-23,5e-24,1e-24,5e-25,1e-25,5e-26,1e-26,
                  5e-27,1e-27,5e-28,1e-28,5e-29,1e-29,5e-30,1e-30,5e-31,1e-31,5e-32,1e-32,5e-33,1e-33,5e-34,1e-34,5e-35,1e-35,
                  5e-36,1e-36,5e-37,1e-37,5e-38,1e-38,5e-39,1e-39,5e-40,1e-40)) {
        for (f in c(1.5,2,2.5,3,3.5,4)) {
          UP_new<-UP[which(UP$F > a & UP$P < b & UP$L > c & UP$LP < d & UP$C < f),]
          Down_new<-Down[which(Down$F > a & Down$P < b & Down$L > c & Down$LP < d & Down$C < f),]
          if(length(rownames(UP_new)) > length(rownames(Down_new))){
            for (e in 1:length(rownames(UP_new))) {
              if(! is.na(UP_new$Gene.Name[e])){
                if(UP_new$Gene.Name[e] == "Ifng"){
                  g<-length(rownames(UP_new))
                  h<-length(rownames(Down_new))
                  i<-length(intersect(UP_new$Gene.Name,RNA$gene))
                  j<-data.frame(a,b,c,d,f,g,h,i)
                  colnames(j)<-c("a","b","c","d","f","g","h","i")
                  k<-rbind(k,j)
                  break
                }
              }
            }
          }
        }
      }
    }
  }
}
