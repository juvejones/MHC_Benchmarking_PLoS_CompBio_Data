library(ROCR)
library(pROC)
library(psychometric)
library(RColorBrewer)

## Some constants
COLUMN_SEP = "\t"
## Plot individual alleles, need more colors
cols1 <- colorRampPalette(brewer.pal(n=9,name="Set1"))(45)

## Plot method summary, only ten colors needed
cols2 <- colorRampPalette(brewer.pal(n=7,name="Set3"))(7)

### Function to calculate auc and srcc for individual allele ###
allele_spec_plot <- function(methods, input_file){
  for (j in 1:length(methods)){
    method <- methods[j]
    output_file <- paste(method, "_auc_srcc.txt", sep="")
    dat <- read.table(input_file, sep="", comment.char="#")
    print(paste(method, "---->"))
    write("Allele data group to be plotted", stderr())
    print(dat[[1]])
    
    ## Initialize lists in the data frame
    df <- data.frame(dat[[1]])
    auc_row <- list()
    srcc_row <- list()
    ## Look for matching filenames, open files, and do necessary operations
    pdf(paste(method, "_allele_roc.pdf", sep=""), width=4,height=4)
    firstplot <- FALSE
    tag <- 0
    for (i in 1:length(dat[[1]])){
      print(paste(" Process Allele: ", dat[[1]][i], sep = ""))
      infile <- paste("../", method, "/Rdata/", dat[[1]][i], ".txt", sep = "")
      if (file.exists(infile)){
      	if (method == "mhcflurry"){
        	datafile <- read.table(infile,header=FALSE,sep=",",skip=1,
                              col.names=c("id","peptide","binding_core","meas_nm","meas_bi","meas_contin","predict","predict_rank"))
        }
        else {
        	datafile <- read.table(infile,header=TRUE,sep=COLUMN_SEP)
        }
        if (tag == 0){
          firstplot <- TRUE
        }
        tag <- tag + 1
      }
      else{
        print(paste("file not exist for:", dat[[1]][i]))
        auc_row[i] <- "NA"
        srcc_row[i] <- "NA"
        next
      }
     
      meas_bi <- list()
      ## Rectify meas_bi column of mhcflurry
      for (j in 1:length(datafile$meas_contin)){
        if (as.character(datafile$meas_contin[j]) > 0.438){
          meas_bi[j] = 1
        }
        else{
          meas_bi[j] = 0
        }
      } 
      ## Three of the methods use RANK as affinity meas and should be taken care of
      if (method == "consensus" || method == "comblib" || method == "tepitope"){
        pred <- prediction((1-datafile$predict_rank),datafile$meas_bi)
        srcc_row[i] <- cor(datafile$meas_contin,(1-datafile$predict_rank),method="spearman")
      }
      else if (method == "mhcflurry"){
        pred <- prediction(datafile$predict, as.numeric(meas_bi))
        srcc_row[i] <- cor(datafile$meas_contin,datafile$predict,method="spearman")
      }
      else{
        pred <- prediction(datafile$predict,datafile$meas_bi)
        srcc_row[i] <- cor(datafile$meas_contin,datafile$predict,method="spearman")
      }
      
      perf <- performance(pred,"tpr","fpr")
      if (firstplot){
        plot(perf,main="ROC Curves",col=cols1[i],
             xaxis.lwd=1.5,yaxis.lwd=1.5,cex.lab=1,xaxis.cex.axis=1,yaxis.cex.axis=1,add=FALSE)
        firstplot <- FALSE
      }
      else{
        plot(perf,main="ROC Curves",col=cols1[i],
             xaxis.lwd=1.5,yaxis.lwd=1.5,cex.lab=1,xaxis.cex.axis=1,yaxis.cex.axis=1,add=TRUE)
      }
      ##generate the random prediction line
      x_ran <- seq(0,1,length=10)
      y_ran <- seq(0,1,length=10)
      lines(x_ran,y_ran,col="gray50",lwd=1,lty=2)
      legend("bottomright",legend=lapply(dat[[1]],as.character),col=cols1,lwd=1.2,bty="n",cex=0.3,ncol=2)
      perf_auc <- performance(pred,"auc",fpr.stop=1.0)
      auc_row[i] <- perf_auc@y.values
    }
    df$auc <- auc_row
    df$srcc <- srcc_row
    write.table(as.matrix(df),output_file,col.names=F, row.names=T, append=F, sep=COLUMN_SEP)
    dev.off()
  }
  return("allele_spec_plot")
}

### Function to calculate auc and srcc for each method using summary.txt ###
### DEPRECATED ###
# all_allele_bi_eval <- function(methods){
#   input_file <- paste("../", method, "/Rdata/summary.txt", sep="")
# 	datafile <- read.table(input_file, col.names=c("peptide","binding_core","meas_nm","meas_bi","meas_contin","predict","predict_rank"), 
# 		sep=COLUMN_SEP, header= FALSE, comment.char="#")
# 	data <- data.frame(datafile)
# 	headerline <- which(with(data,peptide=="peptide"))
# 	data <- data[-headerline,]
#   
#   ##Here pROC package is used for getting ROC_CI
# 	rocobj <- roc(as.numeric(data$meas_bi),as.numeric(data$predict),ci=TRUE,of="auc")
# 	print(rocobj$auc)
# 	print(rocobj$ci)
# 	srcc <- cor(as.numeric(data$meas_contin),as.numeric(data$predict),method="spearman")
# 	srcc_err <- CIr(r=srcc,n=nrow(data),level=.95)
# 	sprintf("%.4f    %.4f", srcc, srcc_err)
# }

### Function to plot auc curves for methods ###
all_allele_bi_eval_plot <- function(methods){
  svg("methods_roc_r11.svg", width=3.5,height=3.5)
  
  for (i in 1:length(methods)){
    method <- methods[i]
    input_file <- paste("../", method, "/Rdata/summary.txt", sep="")
    if (method == "mhcflurry"){
    	datafile <- read.table(input_file, col.names=c("id","peptide","binding_core","meas_nm","meas_bi","meas_contin","predict","predict_rank"), 
                           sep=",", header= FALSE, comment.char="#")
    }
    else {
    	datafile <- read.table(input_file, col.names=c("peptide","binding_core","meas_nm","meas_bi","meas_contin","predict","predict_rank"), 
                           sep=COLUMN_SEP, header= FALSE, comment.char="#")
    }
    data <- data.frame(datafile)
    headerline <- which(with(data,peptide=="peptide"))
    data <- data[-headerline,]

    meas_bi <- list()
    ## Rectify meas_bi column of mhcflurry
    for (j in 1:length(data$meas_contin)){
    	if (as.numeric(as.character(data$meas_contin[j])) > 0.362){
            meas_bi[j] = 1
        }
        else{
            meas_bi[j] = 0
        }
    }
    ## Three of the methods use RANK as affinity meas and should be taken care of
    if (method == "consensus" || method == "comblib" || method == "tepitope"){
      pred <- prediction((1-as.numeric(data$predict_rank)), as.numeric(meas_bi))
    }
    else if (method == "mhcflurry"){
      pred <- prediction(as.numeric(data$predict), as.numeric(meas_bi))
    }
    else{
      pred <- prediction(as.numeric(data$predict), as.numeric(meas_bi))
    }
    
    perf <- performance(pred,"tpr","fpr")
    if (i == 1){
      plot(perf,main="ROC Curves",col=cols2[i],lwd=1.5,
           xaxis.lwd=1.5,yaxis.lwd=1.5,cex.lab=1,xaxis.cex.axis=1,yaxis.cex.axis=1,add=FALSE)
    }
    else{
      plot(perf,main="ROC Curves",col=cols2[i],lwd=1.5,
           xaxis.lwd=1.5,yaxis.lwd=1.5,cex.lab=1,xaxis.cex.axis=1,yaxis.cex.axis=1,add=TRUE)
    }
    
    ##generate tabulated report
    ## Three of the methods use RANK as affinity meas and should be taken care of
    #if (method == "consensus" || method == "comblib" || method == "tepitope"){
    #  rocobj <- roc(as.numeric(meas_bi),(1-as.numeric(data$predict_rank)),ci=TRUE,of="auc")
    #  srcc <- cor(as.numeric(data$meas_contin),(1-as.numeric(data$predict_rank)),method="spearman")
    #}
    #else if (method  == "mhcflurry"){
    #  rocobj <- roc(as.numeric(meas_bi),as.numeric(data$predict),ci=TRUE,of="auc")
    #  srcc <- cor(as.numeric(data$meas_contin),as.numeric(data$predict),method="spearman")
    #}
    #else{
    #  rocobj <- roc(as.numeric(meas_bi),as.numeric(data$predict),ci=TRUE,of="auc")
    #  srcc <- cor(as.numeric(data$meas_contin),as.numeric(data$predict),method="spearman")
    #}
    
    #print(method)
    #print(rocobj$auc) 
    #print(rocobj$ci) 
    #srcc_err <- CIr(r=srcc,n=nrow(data),level=.95)
    #print(paste(srcc, srcc_err, sep="\t"))
  }
  ##generate the random prediction line
  x_ran <- seq(0,1,length=10)
  y_ran <- seq(0,1,length=10)
  lines(x_ran,y_ran,col="gray50",lwd=1,lty=2)
  x2_ran <- seq(0,1,length=20)
  y2_ran <- seq(0.8,0.8,length=20)
  lines(x2_ran,y2_ran,col="black",lwd=1,lty=2)
  
  ##generate legend and finish plot
  legend("bottomright",legend=lapply(methods,as.character),col=cols2,lwd=1.5,bty="n",cex=0.8,ncol=1)
  dev.off()
  
  return("all_allele_bi_eval_plot")
}
## Read list of alleles from a file
args=(commandArgs(TRUE))
print(args)

if(length(args)<1){
	write("ERROR: Invalid arguments supplied.\n 
		Usage: '--args allele_file=\"allele.txt\"'", stderr())
	stop("No arguments supplied.")
} else {
	eval(parse(text=args[]))
}

#methods1 = c("NetMHCIIpan", "nn_align", "smm_align")
#methods2 = c("consensus", "comblib", "tepitope")
methods_all = c("NetMHCIIpan", "nn_align", "smm_align", "consensus", "comblib", "tepitope","mhcflurry")
#methods_all = c("mhcflurry")

## Functions return 'job' to reflect what has been done
#job1 <- allele_spec_plot(methods_all, allele_file)
job2 <- all_allele_bi_eval_plot(methods_all)

#print(paste("!Function:",job1, "Finished!"))
print(paste("!Function:",job2, "Finished!"))
