library(ggplot2)
library(gridExtra)

lm_eqn <- function(m){
    eq <- substitute(~~italic(r)^2~"="~r2, list(r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));                 
}

theme_set(theme_bw(base_size = 18)) 
jpeg(filename="R2_50nM.jpeg", width=11, height=11, units="in", res=600)

datafile <- read.table("../NetMHC4/Rdata/summary.txt", col.names=c("peptide","meas_nm","meas_bi","meas_contin","predict","predict_rank"),sep="\t", header= FALSE, comment.char="#")
data <- data.frame(datafile)
headerline <- which (with(data,peptide=="peptide"))
data <- data[-headerline,]
data <- data[which (with(data, as.numeric(as.character(meas_contin))>0.6)),]

data$meas <- as.numeric(as.character(data$meas_contin))
data$pred <- as.numeric(as.character(data$predict))

m <- lm(pred ~ meas, data)
p1 <- ggplot(data, aes(meas, pred))
p1 <- p1 + geom_point(shape=21, fill="dodgerblue", alpha=0.1)
p1 <- p1 + geom_smooth(method = "lm", se=FALSE, color="red", size=2)
p1 <- p1 + geom_text(aes(x=0.7, y=1.0,label = lm_eqn(lm(pred ~ meas, data))), parse = TRUE)
p1 <- p1 + ylim(0.0,1.0)+ ylab("predicted affinity(1-log(IC50)/log(50k))") + xlab("\nmeasured affinity(1-log(IC50)/log(50k))")
p1 <- p1 + geom_hline(yintercept=0.426, linetype="dashed", color = "dimgrey", size=1.2)
p1 <- p1 + geom_hline(yintercept=0.638, linetype="dotted", color = "dimgrey", size=1.2)
p1 <- p1 + ggtitle("NetMHC4")
p1 <- p1 + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

datafile2 <- read.table("../mhcflurry/Rdata/summary.txt", col.names=c("id","peptide","meas_contin","predict"),sep=",", header= FALSE, comment.char="#")
data2 <- data.frame(datafile2)
# headerline2 <- which (with(data2,peptide=="peptide"))
# data2 <- data2[-headerline2,]
data2 <- data2[which (with(data2, as.numeric(as.character(meas_contin))>0.6)),]

data2$meas <- as.numeric(as.character(data2$meas_contin))
data2$pred <- as.numeric(as.character(data2$predict))

m <- lm(pred ~ meas, data2)
p2 <- ggplot(data2, aes(meas, pred))
p2 <- p2 + geom_point(shape=21, fill="dodgerblue", alpha=0.1)
p2 <- p2 + geom_smooth(method = "lm", se=FALSE, color="red", size=2)
p2 <- p2 + geom_text(aes(x=0.7, y=1.0,label = lm_eqn(lm(pred ~ meas, data2))), parse = TRUE)
p2 <- p2 + ylim(0.0,1.0)+ ylab("predicted affinity(1-log(IC50)/log(50k))") + xlab("\nmeasured affinity(1-log(IC50)/log(50k))")
p2 <- p2 + geom_hline(yintercept=0.426, linetype="dashed", color = "dimgrey", size=1.2)
p2 <- p2 + geom_hline(yintercept=0.638, linetype="dotted", color = "dimgrey", size=1.2)
p2 <- p2 + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p2 <- p2 + ggtitle("mhcflurry - MHC I")

datafile3 <- read.table("/home/zhaoweil/MHC2_IedbNewData_Results/project/nn_align/Rdata/summary.txt", col.names=c("peptide","binding_core","meas_nm","meas_bi","meas_contin","predict","predict_rank"),sep="\t", header= FALSE, comment.char="#")
data3 <- data.frame(datafile3)
headerline3 <- which (with(data3,peptide=="peptide"))
data3 <- data3[-headerline3,]
data3 <- data3[which (with(data3, as.numeric(as.character(meas_contin))>0.6)),]

data3$meas <- as.numeric(as.character(data3$meas_contin))
data3$pred <- as.numeric(as.character(data3$predict))

m <- lm(pred ~ meas, data3)
p3 <- ggplot(data3, aes(meas, pred))
p3 <- p3 + geom_point(shape=21, fill="dodgerblue", alpha=0.1)
p3 <- p3 + geom_smooth(method = "lm", se=FALSE, color="red", size=2)
p3 <- p3 + geom_text(aes(x=0.7, y=1.0,label = lm_eqn(lm(pred ~ meas, data3))), parse = TRUE)
p3 <- p3 + ylim(0.0,1.0)+ ylab("predicted affinity(1-log(IC50)/log(50k))") + xlab("\nmeasured affinity(1-log(IC50)/log(50k))")
p3 <- p3 + geom_hline(yintercept=0.426, linetype="dashed", color = "dimgrey", size=1.2)
p3 <- p3 + geom_hline(yintercept=0.638, linetype="dotted", color = "dimgrey", size=1.2)
p3 <- p3 + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p3 <- p3 + ggtitle("nn_align")

datafile4 <- read.table("/home/zhaoweil/MHC2_IedbNewData_Results/project/mhcflurry/Rdata/summary.txt", col.names=c("id","peptide","binding_core","meas_nm","meas_bi","meas_contin","predict","predict_rank"),sep=",", header= FALSE, comment.char="#")
data4 <- data.frame(datafile4)
headerline4 <- which (with(data4,peptide=="peptide"))
data4 <- data4[-headerline4,]
data4 <- data4[which (with(data4, as.numeric(as.character(meas_contin))>0.6)),]

data4$meas <- as.numeric(as.character(data4$meas_contin))
data4$pred <- as.numeric(as.character(data4$predict))

m <- lm(pred ~ meas, data4)
p4 <- ggplot(data4, aes(meas, pred))
p4 <- p4 + geom_point(shape=21, fill="dodgerblue", alpha=0.1)
p4 <- p4 + geom_smooth(method = "lm", se=FALSE, color="red", size=2)
p4 <- p4 + geom_text(aes(x=0.7, y=1.0,label = lm_eqn(lm(pred ~ meas, data4))), parse = TRUE)
p4 <- p4 + ylim(0.0,1.0)+ ylab("predicted affinity(1-log(IC50)/log(50k))") + xlab("\nmeasured affinity(1-log(IC50)/log(50k))")
p4 <- p4 + geom_hline(yintercept=0.426, linetype="dashed", color = "dimgrey", size=1.2)
p4 <- p4 + geom_hline(yintercept=0.638, linetype="dotted", color = "dimgrey", size=1.2)
p4 <- p4 + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p4 <- p4 + ggtitle("mhcflurry - MHC II")

g <- grid.arrange(arrangeGrob(p1, p2, p3, p4, nrow=2))
dev.off()

