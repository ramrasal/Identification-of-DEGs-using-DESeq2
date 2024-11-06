#import libraries
library(DESeq2)
library(ggplot2)

setwd("C:/Users/rajes/OneDrive/Desktop/RAMP/")

#input files in csv format
count<-read.csv('count_ens_after_conv.csv')
meta<-read.csv('meta.csv')

########################################################
##########Count data preprocessing######################
#count<-count[-c(2:6)]
#count
#write.csv(count,"count_new_ens.csv")

########################################################
dds<-DESeqDataSetFromMatrix(countData=count, 
                               colData=meta, 
                               design=~Condition, tidy = TRUE)
dds<-DESeq(dds)                
res<-results(dds)
res<-na.omit(res)
#########################################################
########VOLCANO PLOT#####################################
#########################################################
jpeg('volcano_ens.jpeg')
par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
topT <- as.data.frame(res)
with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", xlim=c(-10,10), cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value)))

with(subset(topT, padj<0.5 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=1.5))

with(subset(topT, padj<0.5 & log2FoldChange>2), points(log2FoldChange, -log10(padj), pch=20, col="blue", cex=1.5))

dev.off()

###########################################################
############Principal Component Analysis###################
jpeg('pca_ens.jpeg')
vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="Condition")
dev.off()
###########################################################
##########Save up and down regulated genes#################
up<-subset(topT, padj<0.5 & log2FoldChange>2)
down<-subset(topT, padj<0.5 & log2FoldChange<(-2))

write.csv(down,"down_ens.csv")
write.csv(up,"up_ens.csv")
###########################################################

