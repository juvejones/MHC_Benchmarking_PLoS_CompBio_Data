library(ggplot2)
library(gridExtra)
library(DT)
library(reshape2)
library(plyr)

theme_set(theme_bw(base_size = 16))

DF_alleles <- c("A0201","B0702","B3501","B4403","B5301","B5701")
AB_alleles <- c("A0101","A0201","A0301","A2402","A6802","B3501","B4402")
ST_alleles <- c("A0101","A0201","A0301","B0702","B1501","B2705","B4001")

## Merge data for Dana Farber dataset
## 08/07/2018 add MixMHCpred results
df <- vector(mode="list",length=length(DF_alleles))
for (i in 1:length(DF_alleles)){
    allele = DF_alleles[i]
    data <- read.table(paste("./DFRMLI_netMHC4_pred/",allele,".txt", sep=""),header=T,stringsAsFactors = F)
    df[[i]] <- data.frame(x=as.factor(data$elute), y=as.numeric(data$predict))
    colnames(df[[i]]) <- c("Elute","Prediction")
    # df[[i]]$Prediction <- df[[i]]$Prediction/100.0
    df[[i]]$Class[df[[i]]$Elute == 0] <- "Negative"
    df[[i]]$Class[df[[i]]$Elute == 1] <- "Eluted"
    df[[i]]$Method <- "NetMHC4"
    df[[i]]$Cutoff <- 0.426

    data <- read.table(paste("./DFRMLI_netMHCpan4_pred/output/HLA-",allele,".xls",sep=""),header=T,sep="\t",skip = 1)
    df_tmp <- data.frame(x=as.factor(df[[i]]$Elute), y=as.numeric(data$Ave))
    colnames(df_tmp) <- c("Elute","Prediction")
    # df_tmp$Prediction <- df_tmp$Prediction/100.0
    df_tmp$Class <- df[[i]]$Class
    df_tmp$Method <- "NetMHC4pan"
    df_tmp$Cutoff <- 0.5
    df[[i]] <- rbind(df[[i]], df_tmp)

    if (allele != "B5301"){
      data <- read.table(paste("./DFRMLI_MixMHCpred/",allele,".MixMHCpred",sep=""),stringsAsFactors = F,header=T,comment.char = "#",sep="\t")
      df_tmp <- data.frame(x=as.factor(df[[i]]$Elute), y=as.numeric(data$Max_score))
      # Need to scale the Max score
      colnames(df_tmp) <- c("Elute","Prediction")
      df_tmp$Prediction <- (df_tmp$Prediction - min(df_tmp$Prediction))/(max(df_tmp$Prediction)-min(df_tmp$Prediction))
      df_tmp$Class <- df[[i]]$Class
      df_tmp$Method <- "MixMHCpred"
      df_tmp$Cutoff <- 0.5
      colnames(df_tmp) <- colnames(df[[i]])

      df[[i]] <- rbind(df[[i]],df_tmp)
    }

    df[[i]]$Allele <- allele
}

## Merge data for Sternberg dataset
df <- vector(mode="list",length=length(ST_alleles))
for (i in 1:length(ST_alleles)){
  allele = ST_alleles[i]
  data1 <- read.table(paste("./STB_netMHC4_pred/",allele,"_Sternberg2016_test.out", sep=""),header=T,
                      skip=1,stringsAsFactors = F)
  data2 <- read.table(paste("./STB_netMHC4_pred/",allele,"_Sternberg2016_test", sep=""),header=F,stringsAsFactors = F)
  data2 <- data2[nchar(data2$V1) == 9,]

  df[[i]] <- data.frame(x=as.factor(data2$V2), y=as.numeric(data1$nM))
  colnames(df[[i]]) <- c("Elute","Prediction")
  df[[i]]$Prediction <- 1-log10(df[[i]]$Prediction)/log10(50000.0)
  # df[[i]]$Prediction <- df[[i]]$Prediction/100.0
  df[[i]]$Class[df[[i]]$Elute == 0] <- "Negative"
  df[[i]]$Class[df[[i]]$Elute == 1] <- "Eluted"
  df[[i]]$Method <- "NetMHC4"
  df[[i]]$Cutoff <- 0.426

  data <- read.table(paste("./STB_netMHCpan4_pred/",allele,"_NetMHCpan.xls",sep=""),header=T,
                     sep="\t",skip = 1,stringsAsFactors = F)
  data <- data[nchar(data$Peptide)==9,]
  df_tmp <- data.frame(x=as.factor(df[[i]]$Elute), y=as.numeric(data$Ave))
  colnames(df_tmp) <- c("Elute","Prediction")
  # df_tmp$Prediction <- df_tmp$Prediction/100.0
  df_tmp$Class <- df[[i]]$Class
  df_tmp$Method <- "NetMHC4pan"
  df_tmp$Cutoff <- 0.5
  df[[i]] <- rbind(df[[i]], df_tmp)

  data <- read.table(paste("./STB_MixMHCpred/",allele,".MixMHCpred",sep=""),stringsAsFactors = F,header=T,comment.char = "#",sep="\t")
  data <- data[nchar(data$Peptide)==9,]
  df_tmp <- data.frame(x=as.factor(df[[i]]$Elute), y=as.numeric(data$Max_score))
  # Need to scale the Max score
  colnames(df_tmp) <- c("Elute","Prediction")
  df_tmp$Prediction <- (df_tmp$Prediction - min(df_tmp$Prediction))/(max(df_tmp$Prediction)-min(df_tmp$Prediction))
  df_tmp$Class <- df[[i]]$Class
  df_tmp$Method <- "MixMHCpred"
  df_tmp$Cutoff <- 0.5
  colnames(df_tmp) <- colnames(df[[i]])
  df[[i]] <- rbind(df[[i]], df_tmp)

  df[[i]]$Allele <- allele
}

## Merge data for Abelin dataset
## 08/07/2018 add MixMHCpred results
df <- vector(mode="list",length=length(AB_alleles))
for (i in 1:length(AB_alleles)){
  allele = AB_alleles[i]

  data_true <- read.table(paste("./Abelin_netMHCpan4/",allele,".seq.negative.sub.test", sep=""),header=F,stringsAsFactors = F,sep="\t")
  colnames(data_true) <- c("Peptide","Elute")

  data <- read.table(paste("./Abelin_netMHC4/HLA-",allele,".xls", sep=""),header=T,stringsAsFactors = F,sep="\t",skip=1)

  df[[i]] <- merge(data[,c("Peptide","nM")],data_true,by="Peptide")
  df[[i]]$Prediction <- 1-log10(df[[i]]$nM)/log10(50000.0)
  # df[[i]]$Prediction <- df[[i]]$Rank/100.0
  df[[i]]$Class[df[[i]]$Elute == 0] <- "Negative"
  df[[i]]$Class[df[[i]]$Elute == 1] <- "Eluted"
  df[[i]]$Method <- "NetMHC4"
  df[[i]]$Cutoff <- 0.426

  data <- read.table(paste("./Abelin_netMHCpan4/HLA-",allele,".xls",sep=""),stringsAsFactors = F,header=T,sep="\t",skip = 1)
  df_tmp <- merge(data[,c("Peptide","Ave")],data_true,by="Peptide")
  df_tmp$Prediction <- df_tmp$Ave
  df_tmp$Class[df_tmp$Elute == 0] <- "Negative"
  df_tmp$Class[df_tmp$Elute == 1] <- "Eluted"
  df_tmp$Method <- "NetMHC4pan"
  df_tmp$Cutoff <- 0.5

  colnames(df[[i]]) <- colnames(df_tmp)
  df[[i]] <- rbind(df[[i]], df_tmp)

  data <- read.table(paste("./Abelin_MixMHCpred/",allele,".MixMHCpred",sep=""),stringsAsFactors = F,header=T,comment.char = "#",sep="\t")
  df_tmp <- merge(data[,c("Peptide","Max_score")],data_true,by="Peptide")
  # Need to scale the Max score
  df_tmp$Prediction <- (df_tmp$Max_score - min(df_tmp$Max_score))/(max(df_tmp$Max_score)-min(df_tmp$Max_score))
  df_tmp$Class[df_tmp$Elute == 0] <- "Negative"
  df_tmp$Class[df_tmp$Elute == 1] <- "Eluted"
  df_tmp$Method <- "MixMHCpred"
  df_tmp$Cutoff <- 0.5
  colnames(df_tmp) <- colnames(df[[i]])

  df[[i]] <- rbind(df[[i]],df_tmp)

  df[[i]]$Allele <- allele
}
# df_AB_all <- ldply(df,data.frame)

# df_DF_all <- ldply(df,data.frame)
df_ST_all <- ldply(df, data.frame)

## Calculate FPr and FNr using Prediction Rank
FR_df <- vector("list",length(DF_alleles)+length(AB_alleles)+length(ST_alleles))
for (i in 1:length(DF_alleles)){
  allele = DF_alleles[i]
  data <- read.table(paste("./DFRMLI_netMHC4_pred/",allele,".txt", sep=""),header=T,stringsAsFactors = F)
  elute_label <- data$elute
  df_tmp <- data.frame(x=as.factor(data$elute), y=as.numeric(data$pred_rank))
  colnames(df_tmp) <- c("Elute","Rank")
  # df_tmp$Rank <- df_tmp$Rank/100.0
  fnr1 <- dim(df_tmp[(df_tmp$Elute == 1) & (df_tmp$Rank > 0.02),])[[1]]/dim(df_tmp[(df_tmp$Elute == 1),])[[1]]
  fpr1 <- dim(df_tmp[(df_tmp$Elute == 0) & (df_tmp$Rank < 0.02),])[[1]]/dim(df_tmp[(df_tmp$Rank < 0.02),])[[1]]


  data <- read.table(paste("./DFRMLI_netMHCpan4_pred/output/HLA-",allele,".xls",sep=""),header=T,sep="\t",skip = 1)
  df_tmp <- data.frame(x=as.factor(elute_label), y=as.numeric(data$Rank))
  colnames(df_tmp) <- c("Elute","Rank")
  df_tmp$Rank <- df_tmp$Rank/100.0
  fnr2 <- dim(df_tmp[(df_tmp$Elute == 1) & (df_tmp$Rank > 0.02),])[[1]]/dim(df_tmp[(df_tmp$Elute == 1),])[[1]]
  fpr2 <- dim(df_tmp[(df_tmp$Elute == 0) & (df_tmp$Rank < 0.02),])[[1]]/dim(df_tmp[(df_tmp$Rank < 0.02),])[[1]]

  rateTmp <- c(fnr1,fnr2,fpr1,fpr2)
  methodTmp <- c("NetMHC4","NetMHCpan4","NetMHC4","NetMHCpan4")
  nameTmp <- c("FNr","FNr","FDr","FDr")
  FR_df[[i]]$Rate <- rateTmp
  FR_df[[i]]$Method <- methodTmp
  FR_df[[i]]$Metrics <- nameTmp
  FR_df[[i]]$Allele <- allele
  FR_df[[i]]$Dataset <- "Dana Farber"
}

for (i in 1:length(AB_alleles)){
  allele = AB_alleles[i]
  data_true <- read.table(paste("./Abelin_netMHCpan4/",allele,".seq.negative.sub.test", sep=""),header=F,stringsAsFactors = F,sep="\t")
  colnames(data_true) <- c("Peptide","Elute")
  data <- read.table(paste("./Abelin_netMHC4/HLA-",allele,".xls", sep=""),header=T,stringsAsFactors = F,sep="\t",skip=1)
  df_tmp <- merge(data[,c("Peptide","Rank")],data_true,by="Peptide")
  df_tmp$Rank <- df_tmp$Rank/100.0
  fnr1 <- dim(df_tmp[(df_tmp$Elute == 1) & (df_tmp$Rank > 0.02),])[[1]]/dim(df_tmp[(df_tmp$Elute == 1),])[[1]]
  fpr1 <- dim(df_tmp[(df_tmp$Elute == 0) & (df_tmp$Rank < 0.02),])[[1]]/dim(df_tmp[(df_tmp$Rank < 0.02),])[[1]]

  data <- read.table(paste("./Abelin_netMHCpan4/HLA-",allele,".xls",sep=""),stringsAsFactors = F,header=T,sep="\t",skip = 1)
  df_tmp <- merge(data[,c("Peptide","Rank")],data_true,by="Peptide")
  df_tmp$Rank <- df_tmp$Rank/100.0
  fnr2 <- dim(df_tmp[(df_tmp$Elute == 1) & (df_tmp$Rank > 0.02),])[[1]]/dim(df_tmp[(df_tmp$Elute == 1),])[[1]]
  fpr2 <- dim(df_tmp[(df_tmp$Elute == 0) & (df_tmp$Rank < 0.02),])[[1]]/dim(df_tmp[(df_tmp$Rank < 0.02),])[[1]]

  rateTmp <- c(fnr1,fnr2,fpr1,fpr2)
  methodTmp <- c("NetMHC4","NetMHCpan4","NetMHC4","NetMHCpan4")
  nameTmp <- c("FNr","FNr","FDr","FDr")

  FR_df[[i+length(DF_alleles)]]$Rate <- rateTmp
  FR_df[[i+length(DF_alleles)]]$Method <- methodTmp
  FR_df[[i+length(DF_alleles)]]$Metrics <- nameTmp
  FR_df[[i+length(DF_alleles)]]$Allele <- allele
  FR_df[[i+length(DF_alleles)]]$Dataset <- "Abelin"
}

for (i in 1:length(ST_alleles)){
  allele = ST_alleles[i]
  data1 <- read.table(paste("./STB_netMHC4_pred/",allele,"_Sternberg2016_test.out", sep=""),header=T,
                      skip=1,stringsAsFactors = F)
  data2 <- read.table(paste("./STB_netMHC4_pred/",allele,"_Sternberg2016_test", sep=""),header=F,stringsAsFactors = F)
  data2 <- data2[nchar(data2$V1) == 9,]
  df_tmp <- data.frame(x=as.factor(data2$V2), y=as.numeric(data1$Rank))
  colnames(df_tmp)<-c("Elute","Rank")
  df_tmp$Rank <- df_tmp$Rank/100.0

  # only use a subset of negative
  df_tmp <- df_tmp[1:(length(df_tmp[df_tmp$Elute == 1,])*50),]

  fnr1 <- dim(df_tmp[(df_tmp$Elute == 1) & (df_tmp$Rank > 0.02),])[[1]]/dim(df_tmp[(df_tmp$Elute == 1),])[[1]]
  fpr1 <- dim(df_tmp[(df_tmp$Elute == 0) & (df_tmp$Rank < 0.02),])[[1]]/dim(df_tmp[(df_tmp$Rank < 0.02),])[[1]]

  data <- read.table(paste("./STB_netMHCpan4_pred/",allele,"_NetMHCpan.xls",sep=""),header=T,
                     sep="\t",skip = 1,stringsAsFactors = F)
  data <- data[nchar(data$Peptide)==9,]
  df_tmp <- data.frame(x=as.factor(data2$V2), y=as.numeric(data$Rank))
  colnames(df_tmp)<-c("Elute","Rank")
  df_tmp$Rank <- df_tmp$Rank/100.0

  df_tmp <- df_tmp[1:(length(df_tmp[df_tmp$Elute == 1,])*50),]

  fnr2 <- dim(df_tmp[(df_tmp$Elute == 1) & (df_tmp$Rank > 0.02),])[[1]]/dim(df_tmp[(df_tmp$Elute == 1),])[[1]]
  fpr2 <- dim(df_tmp[(df_tmp$Elute == 0) & (df_tmp$Rank < 0.02),])[[1]]/dim(df_tmp[(df_tmp$Rank < 0.02),])[[1]]

  rateTmp <- c(fnr1,fnr2,fpr1,fpr2)
  methodTmp <- c("NetMHC4","NetMHCpan4","NetMHC4","NetMHCpan4")
  nameTmp <- c("FNr","FNr","FDr","FDr")

  FR_df[[i+length(DF_alleles)+length(AB_alleles)]]$Rate <- rateTmp
  FR_df[[i+length(DF_alleles)+length(AB_alleles)]]$Method <- methodTmp
  FR_df[[i+length(DF_alleles)+length(AB_alleles)]]$Metrics <- nameTmp
  FR_df[[i+length(DF_alleles)+length(AB_alleles)]]$Allele <- allele
  FR_df[[i+length(DF_alleles)+length(AB_alleles)]]$Dataset <- "Sternberg"
}

FR_all <- ldply(FR_df, data.frame)
FR_all$Endpoint <- "Pred.Rank"

# ## Calculate FPr and FNr using Prediction Score
FR_df_score <- vector("list",length(DF_alleles)+length(AB_alleles)+length(ST_alleles))
for (i in 1:length(DF_alleles)){
  allele = DF_alleles[i]
  data <- read.table(paste("./DFRMLI_netMHC4_pred/",allele,".txt", sep=""),header=T,stringsAsFactors = F)
  df_tmp <- data.frame(x=as.factor(data$elute), y=as.numeric(data$predict))
  elute_label <- data$elute
  colnames(df_tmp) <- c("Elute","Prediction")
  fnr1 <- dim(df_tmp[(df_tmp$Elute == 1) & (df_tmp$Prediction < 0.426),])[[1]]/dim(df_tmp[(df_tmp$Elute == 1),])[[1]]
  fpr1 <- dim(df_tmp[(df_tmp$Elute == 0) & (df_tmp$Prediction > 0.426),])[[1]]/dim(df_tmp[(df_tmp$Prediction > 0.426),])[[1]]


  data <- read.table(paste("./DFRMLI_netMHCpan4_pred/output/HLA-",allele,".xls",sep=""),header=T,sep="\t",skip = 1)
  df_tmp <- data.frame(x=as.factor(elute_label), y=as.numeric(data$Ave))
  colnames(df_tmp) <- c("Elute","Prediction")
  fnr2 <- dim(df_tmp[(df_tmp$Elute == 1) & (df_tmp$Prediction < 0.5),])[[1]]/dim(df_tmp[(df_tmp$Elute == 1),])[[1]]
  fpr2 <- dim(df_tmp[(df_tmp$Elute == 0) & (df_tmp$Prediction > 0.5),])[[1]]/dim(df_tmp[(df_tmp$Prediction > 0.5),])[[1]]

  if (allele != "B5301") {
    data <- read.table(paste("./DFRMLI_MixMHCpred/",allele,".MixMHCpred",sep=""),stringsAsFactors = F,header=T,comment.char = "#",sep="\t")
    df_tmp <- data.frame(x=as.factor(elute_label), y=as.numeric(data$Max_score))
    colnames(df_tmp) <- c("Elute","Prediction")
    df_tmp$Prediction <- (df_tmp$Prediction - min(df_tmp$Prediction))/(max(df_tmp$Prediction)-min(df_tmp$Prediction))
    fnr3 <- dim(df_tmp[(df_tmp$Elute == 1) & (df_tmp$Prediction < 0.75),])[[1]]/dim(df_tmp[(df_tmp$Elute == 1),])[[1]]
    fpr3 <- dim(df_tmp[(df_tmp$Elute == 0) & (df_tmp$Prediction > 0.75),])[[1]]/dim(df_tmp[(df_tmp$Prediction > 0.75),])[[1]]
  }
  else {
    fnr3 <- NA
    fpr3 <- NA
  }

  rateTmp <- c(fnr1,fnr2,fnr3,fpr1,fpr2,fpr3)
  methodTmp <- c("NetMHC4","NetMHCpan4","MixMHCpred","NetMHC4","NetMHCpan4","MixMHCpred")
  nameTmp <- c("FNr","FNr","FNr","FDr","FDr","FDr")
  FR_df_score[[i]]$Rate <- rateTmp
  FR_df_score[[i]]$Method <- methodTmp
  FR_df_score[[i]]$Metrics <- nameTmp
  FR_df_score[[i]]$Allele <- allele
  FR_df_score[[i]]$Dataset <- "Dana Farber"
}

for (i in 1:length(AB_alleles)){
  allele = AB_alleles[i]
  data_true <- read.table(paste("./Abelin_netMHCpan4/",allele,".seq.negative.sub.test", sep=""),header=F,stringsAsFactors = F,sep="\t")
  colnames(data_true) <- c("Peptide","Elute")
  data <- read.table(paste("./Abelin_netMHC4/HLA-",allele,".xls", sep=""),header=T,stringsAsFactors = F,sep="\t",skip=1)
  df_tmp <- merge(data[,c("Peptide","nM")],data_true,by="Peptide")
  df_tmp$Prediction <- 1-log10(df_tmp$nM)/log10(50000.0)
  fnr1 <- dim(df_tmp[(df_tmp$Elute == 1) & (df_tmp$Prediction < 0.426),])[[1]]/dim(df_tmp[(df_tmp$Elute == 1),])[[1]]
  fpr1 <- dim(df_tmp[(df_tmp$Elute == 0) & (df_tmp$Prediction > 0.426),])[[1]]/dim(df_tmp[(df_tmp$Prediction > 0.426),])[[1]]

  data <- read.table(paste("./Abelin_netMHCpan4/HLA-",allele,".xls",sep=""),stringsAsFactors = F,header=T,sep="\t",skip = 1)
  df_tmp <- merge(data[,c("Peptide","Ave")],data_true,by="Peptide")
  df_tmp$Prediction <- df_tmp$Ave
  fnr2 <- dim(df_tmp[(df_tmp$Elute == 1) & (df_tmp$Ave < 0.5),])[[1]]/dim(df_tmp[(df_tmp$Elute == 1),])[[1]]
  fpr2 <- dim(df_tmp[(df_tmp$Elute == 0) & (df_tmp$Ave > 0.5),])[[1]]/dim(df_tmp[(df_tmp$Ave > 0.5),])[[1]]

  data <- read.table(paste("./Abelin_MixMHCpred/",allele,".MixMHCpred",sep=""),stringsAsFactors = F,header=T,comment.char = "#",sep="\t")
  df_tmp <- merge(data[,c("Peptide","Max_score")],data_true,by="Peptide")
  df_tmp$Prediction <- df_tmp$Max_score
  df_tmp$Prediction <- (df_tmp$Prediction - min(df_tmp$Prediction))/(max(df_tmp$Prediction)-min(df_tmp$Prediction))
  fnr3 <- dim(df_tmp[(df_tmp$Elute == 1) & (df_tmp$Prediction < 0.75),])[[1]]/dim(df_tmp[(df_tmp$Elute == 1),])[[1]]
  fpr3 <- dim(df_tmp[(df_tmp$Elute == 0) & (df_tmp$Prediction > 0.75),])[[1]]/dim(df_tmp[(df_tmp$Prediction > 0.75),])[[1]]

  rateTmp <- c(fnr1,fnr2,fnr3,fpr1,fpr2,fpr3)
  methodTmp <- c("NetMHC4","NetMHCpan4","MixMHCpred","NetMHC4","NetMHCpan4","MixMHCpred")
  nameTmp <- c("FNr","FNr","FNr","FDr","FDr","FDr")

  FR_df_score[[i+length(DF_alleles)]]$Rate <- rateTmp
  FR_df_score[[i+length(DF_alleles)]]$Method <- methodTmp
  FR_df_score[[i+length(DF_alleles)]]$Metrics <- nameTmp
  FR_df_score[[i+length(DF_alleles)]]$Allele <- allele
  FR_df_score[[i+length(DF_alleles)]]$Dataset <- "Abelin"
}

for (i in 1:length(ST_alleles)){
  allele = ST_alleles[i]
  data1 <- read.table(paste("./STB_netMHC4_pred/",allele,"_Sternberg2016_test.out", sep=""),header=T,
                      skip=1,stringsAsFactors = F)
  data2 <- read.table(paste("./STB_netMHC4_pred/",allele,"_Sternberg2016_test", sep=""),header=F,stringsAsFactors = F)
  data2 <- data2[nchar(data2$V1) == 9,]
  df_tmp <- data.frame(x=as.factor(data2$V2), y=as.numeric(data1$nM))
  colnames(df_tmp)<-c("Elute","nM")
  df_tmp$Prediction <- 1-log10(df_tmp$nM)/log10(50000.0)

  # only use a subset of negative
  df_tmp <- df_tmp[1:(length(df_tmp[df_tmp$Elute == 1,])*50),]

  fnr1 <- dim(df_tmp[(df_tmp$Elute == 1) & (df_tmp$Prediction < 0.426),])[[1]]/dim(df_tmp[(df_tmp$Elute == 1),])[[1]]
  fpr1 <- dim(df_tmp[(df_tmp$Elute == 0) & (df_tmp$Prediction > 0.426),])[[1]]/dim(df_tmp[(df_tmp$Prediction > 0.426),])[[1]]

  data <- read.table(paste("./STB_netMHCpan4_pred/",allele,"_NetMHCpan.xls",sep=""),header=T,
                     sep="\t",skip = 1,stringsAsFactors = F)
  data <- data[nchar(data$Peptide)==9,]
  df_tmp <- data.frame(x=as.factor(data2$V2), y=as.numeric(data$Ave))
  colnames(df_tmp)<-c("Elute","Ave")

  df_tmp <- df_tmp[1:(length(df_tmp[df_tmp$Elute == 1,])*50),]

  fnr2 <- dim(df_tmp[(df_tmp$Elute == 1) & (df_tmp$Ave < 0.5),])[[1]]/dim(df_tmp[(df_tmp$Elute == 1),])[[1]]
  fpr2 <- dim(df_tmp[(df_tmp$Elute == 0) & (df_tmp$Ave > 0.5),])[[1]]/dim(df_tmp[(df_tmp$Ave > 0.5),])[[1]]

  data <- read.table(paste("./STB_MixMHCpred/",allele,".MixMHCpred",sep=""),stringsAsFactors = F,header=T,comment.char = "#",sep="\t")
  data <- data[nchar(data$Peptide)==9,]
  df_tmp <- data.frame(x=as.factor(data2$V2), y=as.numeric(data$Max_score))
  colnames(df_tmp) <- c("Elute","Prediction")

  df_tmp <- df_tmp[1:(length(df_tmp[df_tmp$Elute == 1,])*50),]

  df_tmp$Prediction <- (df_tmp$Prediction - min(df_tmp$Prediction))/(max(df_tmp$Prediction)-min(df_tmp$Prediction))
  fnr3 <- dim(df_tmp[(df_tmp$Elute == 1) & (df_tmp$Prediction < 0.75),])[[1]]/dim(df_tmp[(df_tmp$Elute == 1),])[[1]]
  fpr3 <- dim(df_tmp[(df_tmp$Elute == 0) & (df_tmp$Prediction > 0.75),])[[1]]/dim(df_tmp[(df_tmp$Prediction > 0.75),])[[1]]

  rateTmp <- c(fnr1,fnr2,fnr3,fpr1,fpr2,fpr3)
  methodTmp <- c("NetMHC4","NetMHCpan4","MixMHCpred","NetMHC4","NetMHCpan4","MixMHCpred")
  nameTmp <- c("FNr","FNr","FNr","FDr","FDr","FDr")

  FR_df_score[[i+length(DF_alleles)+length(AB_alleles)]]$Rate <- rateTmp
  FR_df_score[[i+length(DF_alleles)+length(AB_alleles)]]$Method <- methodTmp
  FR_df_score[[i+length(DF_alleles)+length(AB_alleles)]]$Metrics <- nameTmp
  FR_df_score[[i+length(DF_alleles)+length(AB_alleles)]]$Allele <- allele
  FR_df_score[[i+length(DF_alleles)+length(AB_alleles)]]$Dataset <- "Sternberg"
}
FR_all_score <- ldply(FR_df_score, data.frame)
FR_all_score$Endpoint <- "Pred.Score"

FR <- rbind(FR_all, FR_all_score)

## Boxplot ##
# jpeg(filename="box_r2_STBdata.jpeg", width=11, height=5, units="in", res=150)
# p <- ggplot(df_ST_all, aes(Allele,Prediction,fill=Class)) +
#   # geom_jitter(aes(Allele, Prediction, color=Class),df_DF_all, size=0.5) +
#   geom_boxplot(outlier.alpha = 0.2) +
#   scale_fill_manual(values = c("orange","cornflowerblue")) +
#   geom_hline(data=df_ST_all, aes(yintercept=Cutoff), colour="dimgrey", linetype="dashed") +
#   theme(axis.text.x = element_text(angle=65, vjust=0.6)) +
#   # scale_y_log10() +
#   facet_wrap(~Method) +
#   ggtitle("Sternberg Dataset") + labs(x="HLA Type", y="Prediction")
# print(p)

## Barplot ##
jpeg(filename="barplot_Fr_score_r2.jpeg",width = 11, height=5,units="in",res=300)
# FR_tmp = FR_all[FR_all$Endpoint == "Pred.Score",]
FR_tmp = FR_all_score
p <- ggplot(FR_tmp, aes(x=Allele,y=Method,fill=Rate)) +
  geom_tile(colour = "white") +
  facet_grid(Metrics~Dataset) +
  scale_fill_distiller(palette="Spectral") +
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) +
  ggtitle("False Classification Rates") + labs(x="HLA Type", y="Method")
print(p)
dev.off()

