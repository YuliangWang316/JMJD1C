b<-read.table("c:/Users/xjmik/Desktop/STAT3_new.txt",sep = "\t",header = TRUE,row.names = 1)

b$n<-b$l / b$m
b_new<-b[which(b$i<50),]
b_new_new<-b_new[which(b_new$m<3000),]
b_new_new_new<-b_new_new[which(b_new_new$l>3000),]

