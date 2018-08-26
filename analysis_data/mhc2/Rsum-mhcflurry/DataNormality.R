library(ROCR)

allele_list <- read.table("../CountLogStrip", sep="")
training_count <- read.table("Count_bd2013_alleles_for_benchmark", sep="")
df <- data.frame(allele_list[[1]])
auc <- list()
tr_count <- list()
test_count <- list()
skewness <- list()
binorm <- list()

for (j in 1:length(allele_list[[1]])){
  allele <- allele_list[[1]][j]
  file <- paste("../netmhc4/Rdata/",allele,".txt", sep="")
  if (file.exists(file)){
    datafile <- read.table(file, 
                           col.names=c("peptide","meas_nm","meas_bi",
                                       "meas_contin","predict","predict_rank"),sep="\t", 
                           header= FALSE, comment.char="#")
  }
  else{
    print(paste("file not exist for allele:", allele))
    auc[j] <- NULL
    tr_count[j] <- training_count[[1]][j]
    skewness[j] <- NULL
    next
  }
  data <- data.frame(datafile)
  headerline <- which (with(data,peptide=="peptide"))
  data <- data[-headerline,]
  
  auc[j] <- tryCatch(
    {
      pred <- prediction(as.numeric(data$predict),as.numeric(data$meas_bi))
      perf_auc <- performance(pred,"auc",fpr.stop=1.0)
      auc[j] <- perf_auc@y.values
    },
    error=function(cond){
      message(paste(cond, "\n", sep=""))
      return(NA)
    }
  )
  
  tr_count[j] <- training_count[[1]][j] 
  
  count <- c(0,0)

  for (i in 1:length(data$meas_bi)){
    if (as.numeric(as.character(data$meas_bi[i])) == 1){
      count[1] <- count[1] + 1
    } 
    else{
      count[2] <- count[2] + 1
    }
  }
  test_count[j] <- count[1] + count[2]
  skewness[j] <- count[1]/(count[1] + count[2])
  binorm_obj <- binom.test(count, p = 0.5,conf.level = 0.95)
  binorm[j] <- binorm_obj$p.value
}
df$auc <- auc
df$tr_count <- tr_count
df$test_count <- test_count
df$sk <- skewness
df$binorm <- binorm
print(as.matrix(df))
outputfile <- "test"
write.table(as.matrix(df), outputfile, col.names=c("Allele", "AUC", "Tr_count", "Test_count", "Skn", "Binorm_P"), 
            row.names=F, append=F, sep="\t")
