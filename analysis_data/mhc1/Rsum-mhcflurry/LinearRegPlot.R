hlas <- read.table("../NetMHC4/Rdata/allele.txt",header=FALSE,comment.char="#")

hlas <- as.character(hlas$V1)
rs <- list()

for (hla in hlas){
	datafile <- read.table(paste("../NetMHC4/Rdata/",hla,".txt",sep=""), sep="\t", header= TRUE, comment.char="#")
	data <- data.frame(datafile)
	#headerline <- which (with(data,peptide=="peptide"))
	#data <- data[-headerline,]
	meas <- as.numeric(as.character(data$meas_contin))
	pred <- as.numeric(as.character(data$predict))
	lg <- lm(meas~pred)
	rs <- c(rs,summary(lg)$r.squared)
}
output = cbind(hlas, rs)
write.csv(output, file="./allele_rsquared2.txt",sep="\t",col.names=TRUE)
#plot(meas,pred,pch=1,cex=1.5,col="red",main="mhcflurry",xlab="1 - Ln(Meas)/Ln(50k)",ylab="1 - Ln(Pred)/Ln(50k)")
