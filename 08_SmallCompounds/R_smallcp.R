setwd("C:/Users/zliu39/OneDrive - City University of Hong Kong/Desktop/COVID19/09_smallcomp")
remove(list=ls())

library("readxl")
library("tidyverse")
sm_cp<-read_excel("raw_results.xlsx", sheet = 1, na = "NA")
dim(sm_cp)
#[1] 423422     21

sm_cp_fil<-sm_cp[which(sm_cp$norm_cs<0),]
dim(sm_cp_fil)
#[1] 245096     21

sm_cp_TOP2A<-sm_cp_fil[str_which(sm_cp_fil$target_name,"TOP2A"),]
dim(sm_cp_TOP2A)
#[1] 398  21
write.csv(sm_cp_TOP2A, "TOP2A_res.csv",quote = F, row.names = F)

sm_cp_AURKB<-sm_cp_fil[str_which(sm_cp_fil$target_name,"AURKB"),]
dim(sm_cp_AURKB)
#[1] 184  21
write.csv(sm_cp_AURKB, "AURKB_res.csv",quote = F, row.names = F)

sm_cp_KIF11<-sm_cp_fil[str_which(sm_cp_fil$target_name,"KIF11"),]
dim(sm_cp_KIF11)
#[1] 55 21
write.csv(sm_cp_KIF11, "KIF11_res.csv",quote = F, row.names = F)


