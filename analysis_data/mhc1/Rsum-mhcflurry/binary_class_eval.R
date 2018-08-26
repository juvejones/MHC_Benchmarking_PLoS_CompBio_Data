library(ROCR)
library(pROC)
library(psychometric)
library(RColorBrewer)

## Some constants
## Plot individual alleles, need more colors
#cols <- colorRampPalette(brewer.pal(n=9,name="Set1"))(36)

## Plot method summary, only ten colors needed
cols <- colorRampPalette(brewer.pal(n=11,name="Set3"))(11)

### Function to calculate auc and srcc for individual allele ###
allele_spec_plot <- function(input_file, output_file){
	dat <- read.table(input_file, sep="", comment.char="#")
	write("Allele data group to be plotted", stderr())
	print(dat[[1]])

## Initialize lists in the data frame
	df <- data.frame(dat[[1]])
	auc_row <- list()
	srcc_row <- list()

## Look for matching filenames, open files, and do necessary operations
  #pdf("mhcflurry_allele_roc.pdf", width=4,height=4)
	for (i in 1:length(dat[[1]])){
		print(paste(" Process Allele: ", dat[[1]][i], sep = ""))
		if (file.exists(paste("HLA-", dat[[1]][i], ".txt", sep = ""))){
			datafile <- read.table(paste("HLA-", dat[[1]][i], ".txt", sep = ""),header=FALSE, sep=",",
                                            col.names=c("idx","peptide","meas","pred"), comment.char="#")
		}
		else{
			next
		}
                data <- data.frame(datafile)
                #meas_bi <- list()
                #for (j in 1:length(data[[1]])){
                #  if (as.numeric(data$meas[[j]]) > 0.426){
                #    meas_bi[j] <- 1
                #  }
                #  else{
                #    meas_bi[j] <- 0
                #  }
                #}
                #data$meas_bi <- meas_bi
		#pred <- prediction(as.numeric(data$pred),as.numeric(data$meas_bi))
		#perf <- performance(pred,"tpr","fpr")
		#if (i == 1){
		#	plot(perf,main="ROC Curves",col=cols[i],
		#		xaxis.lwd=1.5,yaxis.lwd=1.5,cex.lab=1,xaxis.cex.axis=1,yaxis.cex.axis=1,add=FALSE)
		#}
		#else{
		#	plot(perf,main="ROC Curves",col=cols[i],
		#		xaxis.lwd=1.5,yaxis.lwd=1.5,cex.lab=1,xaxis.cex.axis=1,yaxis.cex.axis=1,add=TRUE)
		#}
	##generate the random prediction line
		#x_ran <- seq(0,1,length=10)
		#y_ran <- seq(0,1,length=10)
		#lines(x_ran,y_ran,col="gray50",lwd=1,lty=2)
		#legend("bottomright",legend=lapply(dat[[1]],as.character),col=cols,lwd=1.2,bty="n",cex=0.4,ncol=2)
		#perf_auc <- performance(pred,"auc",fpr.stop=1.0)
		#auc_row[i] <- perf_auc@y.values
		srcc_row[i] <- cor(data$meas,data$pred,method="spearman")
	}
	#df$auc <- auc_row
	df$srcc <- srcc_row
	write.table(as.matrix(df),output_file,col.names=F, row.names=T, append=F, sep="\t")
	#dev.off()
}

### Function to calculate auc and srcc for each method using summary.txt ###
all_allele_bi_eval <- function(input_file){
	#datafile <- read.table(input_file, col.names=c("peptide","meas_nm","meas_bi","meas_contin","predict","predict_rank"), 
	#	sep=COLUMN_SEP, header= FALSE, comment.char="#")
	datafile <- read.table(input_file, col.names=c("idx","peptide","meas","pred"),
               sep=",", header=FALSE, comment.char="#")
        data <- data.frame(datafile)
	#headerline <- which(with(data,peptide=="peptide"))
	#data <- data[-headerline,]
        ## Add missing binary labels to data
        meas_bi <- list()
        for (j in 1:length(data[[1]])){
          if (as.numeric(data$meas[[j]]) > 0.426){
            meas_bi[j] <- 1
          }
          else{
            meas_bi[j] <- 0
          } 
        }
        data$meas_bi <- meas_bi
  
  ##Here pROC package is used for getting ROC_CI
	rocobj <- roc(as.numeric(data$meas_bi),as.numeric(data$pred),ci=TRUE,of="auc")
	print(rocobj$auc)
	print(rocobj$ci)
	srcc <- cor(as.numeric(data$meas),as.numeric(data$pred),method="spearman")
	srcc_err <- CIr(r=srcc,n=nrow(data),level=.95)
        print(srcc) 
        print(srcc_err)
}

### Function to plot auc curves for methods ###
all_allele_bi_eval_plot <- function(methods){
  pdf("methods_roc_new_r1.pdf", width=9,height=4.5)
  par(mfrow=c(1,2))

  for (i in 1:length(methods)){
    method <- methods[i]
    if ((method == "mhcflurry") || (method == "mhcflurry_pan") || (method == "NetMHCpan4")){
      input_file <- paste("../", method, "/Rdata/summary.txt", sep="")
      datafile <- read.table(input_file, col.names=c("idx","peptide","meas","pred"),
                             sep=",", header= FALSE, comment.char="#")
      data <- data.frame(datafile)
      meas_bi <- list()
      for (j in 1:length(data[[1]])){
        if (as.numeric(data$meas[[j]]) > 0.426){
          meas_bi[j] <- 1
        }
        else{
          meas_bi[j] <- 0
        }
      }
      data$meas_bi <- meas_bi
      pred <- prediction(as.numeric(data$pred),as.numeric(data$meas_bi))
      perf <- performance(pred,"tpr","fpr")
    }
    else{
      input_file <- paste("../", method, "/Rdata/summary.txt", sep="")
      datafile <- read.table(input_file, col.names=c("peptide","meas_nm","meas_bi","meas_contin","predict","predict_rank"), 
                             sep="\t", header= FALSE, comment.char="#")
      data <- data.frame(datafile)
      headerline <- which(with(data,peptide=="peptide"))
      data <- data[-headerline,]
      pred <- prediction(as.numeric(data$predict),as.numeric(data$meas_bi))
      perf <- performance(pred,"tpr","fpr")
    }
    if (i == 1){
      plot(perf,main="ROC Curves",col=cols[i],lwd=1.5,
           xaxis.lwd=1.5,yaxis.lwd=1.5,cex.lab=1,xaxis.cex.axis=1,yaxis.cex.axis=1,add=FALSE)
    }
    else{
      plot(perf,main="ROC Curves",col=cols[i],lwd=1.5,
           xaxis.lwd=1.5,yaxis.lwd=1.5,cex.lab=1,xaxis.cex.axis=1,yaxis.cex.axis=1,add=TRUE)
    }
    if (i == 1){
      plot(perf,main="ROC Curves",col=cols[i],lwd=1.5, xlim=c(0.7,0.9),ylim=c(0.7,0.9),
           xaxis.lwd=1.5,yaxis.lwd=1.5,cex.lab=1,xaxis.cex.axis=1,yaxis.cex.axis=1,add=FALSE)
    }
    else{
      plot(perf,main="ROC Curves",col=cols[i],lwd=1.5, xlim=c(0.7,0.9),ylim=c(0.7,0.9),
           xaxis.lwd=1.5,yaxis.lwd=1.5,cex.lab=1,xaxis.cex.axis=1,yaxis.cex.axis=1,add=TRUE)
    }

  }
  ##generate the random prediction line
  x_ran <- seq(0,1,length=10)
  y_ran <- seq(0,1,length=10)
  lines(x_ran,y_ran,col="gray50",lwd=1,lty=2)
  
  ##generate legend and finish plot
  legend("bottomright",legend=lapply(methods,as.character),col=cols,lwd=1.5,bty="n",cex=0.7,ncol=1)
  dev.off()
  
  return("all_allele_bi_eval_plot")
}
## Read list of alleles from a file
args=(commandArgs(TRUE))
print(args)

if(length(args)<1){
	write("ERROR: Invalid arguments supplied.\n 
		Usage: '--args allele_file=\"allele.txt\" input_file=\"summary.txt\" output_file=\"allele_stat.txt\"'", stderr())
	stop("No arguments supplied.")
} else {
	eval(parse(text=args[]))
}

methods = c("smm", "smmpmbec", "ann", "NetMHC4", "PickPocket", "consensus", "NetMHCpan2.8", "NetMHCpan3", "NetMHCpan4", "NetMHCcons", "mhcflurry","mhcflurry_pan")

## Functions return 'job' to reflect what has been done
job <- allele_spec_plot(allele_file, output_file)
#job1 <- all_allele_bi_eval(input_file)
#job2 <- all_allele_bi_eval_plot(methods)

#print(paste("!Function:",job, "Finished!"))
#print(paste("!Function:",job1, "Finished!"))
