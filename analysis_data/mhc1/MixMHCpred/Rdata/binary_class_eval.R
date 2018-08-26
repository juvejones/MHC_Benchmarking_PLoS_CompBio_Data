library(ROCR)
library(pROC)
library(psychometric)
library(RColorBrewer)

## Some constants
## Plot individual alleles, need more colors
#cols <- colorRampPalette(brewer.pal(n=9,name="Set1"))(36)

## Plot method summary, only ten colors needed
cols <- colorRampPalette(brewer.pal(n=7,name="Set3"))(13)

## Define allele-specific threshold
repertoire <- read.table("../../Rsum-mhcflurry/Paul_2013_IEDB_epitope_repertoire",
                stringsAsFactors=F,header=T,sep="\t")
head(repertoire)
repertoire$AllelesL <- paste('HLA',repertoire$Alleles,sep="-")

### Function to calculate auc and srcc for individual allele ###
allele_spec_plot <- function(input_file, output_file){
	dat <- read.table(input_file, sep="", comment.char="#", stringsAsFactors=F)
	write("Allele data group to be plotted", stderr())
	print(dat[[1]])

## Initialize lists in the data frame
	df <- data.frame(dat[[1]])
	auc_row <- list()
	srcc_row <- list()

## Look for matching filenames, open files, and do necessary operations
	pdf("MixMHCpred_allele_roc.pdf", width=4,height=4)
	for (i in 1:length(dat[[1]])){
		print(paste(" Process Allele: ", dat[[1]][i], sep = ""))
		if (file.exists(paste("HLA-", dat[[1]][i], ".pred", sep = ""))){
			datafile <- read.table(paste("HLA-", dat[[1]][i], ".pred", sep = ""),header=TRUE, sep=",",
				stringsAsFactor=F, comment.char="#", skip=0)
		}
		else{
			next
		}
                data <- data.frame(datafile)
                #meas_bi <- list()
		#if (dat[[1]][i] %in% repertoire$Alleles){
                #  cutoff <- as.numeric(repertoire[repertoire$Alleles == dat[[1]][i],]$nM_75pct)
                #  cutoff <- 1 - log10(cutoff)/log10(50000.0)
                #}
                #else{
                #  cutoff <- 0.426
                #}

                #for (j in 1:length(data[[1]])){
                #  if (as.numeric(data$meas[[j]]) > cutoff){
                #    meas_bi[j] <- 1
                #  }
                #  else{
                #    meas_bi[j] <- 0
                #  }
                #}
                #data$meas_bi <- meas_bi
		pred <- prediction(as.numeric(as.vector(data$Max_score)),as.numeric(as.vector(data$meas_bi)))
		perf <- performance(pred,"tpr","fpr")
		if (i == 1){
			plot(perf,main="ROC Curves",col=cols[i],
				xaxis.lwd=1.5,yaxis.lwd=1.5,cex.lab=1,xaxis.cex.axis=1,yaxis.cex.axis=1,add=FALSE)
		}
		else{
			plot(perf,main="ROC Curves",col=cols[i],
				xaxis.lwd=1.5,yaxis.lwd=1.5,cex.lab=1,xaxis.cex.axis=1,yaxis.cex.axis=1,add=TRUE)
		}
	##generate the random prediction line
		x_ran <- seq(0,1,length=10)
		y_ran <- seq(0,1,length=10)
		lines(x_ran,y_ran,col="gray50",lwd=1,lty=2)
		legend("bottomright",legend=lapply(dat[[1]],as.character),col=cols,lwd=1.2,bty="n",cex=0.4,ncol=2)
		perf_auc <- performance(pred,"auc",fpr.stop=1.0)
		auc_row[i] <- perf_auc@y.values
		srcc_row[i] <- cor(1-log10(data$meas_contin)/log10(50000.0),data$Max_score,method="spearman")
	}
	df$auc <- auc_row
	df$srcc <- srcc_row
	write.table(as.matrix(df),output_file,col.names=F, row.names=T, append=F, sep="\t")
	dev.off()
}

### Function to calculate auc and srcc for each method using summary.txt ###
all_allele_bi_eval <- function(input_file){
	# IEDB tools block
	#datafile <- read.table(input_file, col.names=c("peptide","meas_nm","meas_bi","meas_contin","predict","predict_rank"), 
	#	sep=COLUMN_SEP, header= FALSE, comment.char="#")

	# mhcflurry/NetMHCpan4 block
	#datafile <- read.table(input_file, col.names=c("idx","peptide","meas","pred"),
        #       sep=",", header=FALSE, comment.char="#")
        #data <- data.frame(datafile)
	#headerline <- which(with(data,peptide=="peptide"))
	#data <- data[-headerline,]
        ## Add missing binary labels to data
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

	# MixMHCpred block
	datafile <- read.table(input_file, sep=",", header=TRUE, comment.char="#", stringsAsFactor=F)
	data <- as.data.frame(datafile)
	headerline <- which(with(data,Peptide=="Peptide"))
      	data <- data[-headerline,]
	data$meas_bi <- as.numeric(as.vector(data$`meas_bi`))
	data$pred <- as.numeric(as.vector(data$Max_score))
	data$meas <- 1-log10(as.numeric(as.vector(data$`meas_contin`)))/log10(50000.0)

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
  #svg("methods_roc_addMixpred_r11.svg", width=4.5,height=4.5)

  for (i in 1:length(methods)){
    method <- methods[i]
    if ((method == "mhcflurry") || (method == "mhcflurry_pan") || (method == "NetMHCpan4")){
      input_file <- paste("../../", method, "/Rdata/summary.txt", sep="")
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
    else if (method == "MixMHCpred") {
      input_file <- paste("../../", method, "/Rdata/summary.txt", sep="")
      datafile <- read.table(input_file, sep=",", header=TRUE, comment.char="#", stringsAsFactor=F)
      datafile <- as.data.frame(datafile)
      headerline <- which(with(datafile,Peptide=="Peptide"))
      datafile <- datafile[-headerline,]
      meas_bi <- as.numeric(as.vector(datafile$meas_bi))
      pred <- prediction(as.numeric(as.vector(datafile$Max_score)), meas_bi)
      perf <- performance(pred,"tpr","fpr")

      # Calculate R2
      score <- as.numeric(as.vector(datafile$Max_score))
      meas <- 1-log10(as.numeric(as.vector(datafile$`meas_contin`)))/log10(50000.0)
      lg <- lm(meas~score)
      print(summary(lg))
    }
    else{
      input_file <- paste("../../", method, "/Rdata/summary.txt", sep="")
      datafile <- read.table(input_file, col.names=c("peptide","meas_nm","meas_bi","meas_contin","predict","predict_rank"), 
                             sep="\t", header= FALSE, comment.char="#")
      data <- data.frame(datafile)
      headerline <- which(with(data,peptide=="peptide"))
      data <- data[-headerline,]
      pred <- prediction(as.numeric(data$predict),as.numeric(data$meas_bi))
      perf <- performance(pred,"tpr","fpr")
    }
    if (i == 1){
      plot(perf,main="ROC Curves",xaxis.lwd=1.5,yaxis.lwd=1.5,xaxis.cex.axis=1,yaxis.cex.axis=1,col=cols[i],lwd=1.5,cex.lab=1,add=F, ylim=c(0.0,1.0), xlim=c(0.0,1.0))
    }
    else {
      plot(perf,main="ROC Curves",xaxis.lwd=1.5,yaxis.lwd=1.5,xaxis.cex.axis=1,yaxis.cex.axis=1,col=cols[i],lwd=1.5,cex.lab=1,add=T, ylim=c(0.0,1.0), xlim=c(0.0,1.0))
    }
  }
  ##generate the random prediction line
  x_ran <- seq(0,1,length=10)
  y_ran <- seq(0,1,length=10)
  lines(x_ran,y_ran,col="gray50",lwd=1,lty=2)
  x2_ran <- seq(0,1,length=20)
  y2_ran <- seq(0.8,0.8,length=20)
  #lines(x2_ran,y2_ran,col="black",lwd=1,lty=2)  

  ##generate legend and finish plot
  legend("bottomright",legend=lapply(methods,as.character),col=cols,lwd=1.5,bty="n",cex=0.7,ncol=1)
  #dev.off()
  
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

methods = c("smm", "smmpmbec", "ann", "NetMHC4", "PickPocket", "consensus", "NetMHCpan2.8", "NetMHCpan3", "NetMHCpan4", "NetMHCcons", "mhcflurry","mhcflurry_pan","MixMHCpred")

## Functions return 'job' to reflect what has been done
#job <- allele_spec_plot(allele_file, output_file)
#job1 <- all_allele_bi_eval(input_file)
job2 <- all_allele_bi_eval_plot(methods)

#print(paste("!Function:",job, "Finished!"))
#print(paste("!Function:",job1, "Finished!"))
