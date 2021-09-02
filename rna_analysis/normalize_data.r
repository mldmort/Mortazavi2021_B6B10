
Data_Dir = '../data/'

data = read.table("b6b10_counts_raw.txt",sep="\t",header=T,stringsAsFactors=F)
rownames(data)=data[,"X"]
data$X <- NULL
dim(data)
data[1:5, 1:5]

# keep only B6 and B10 strains
data = data[which(data['Group']=="B6" | data['Group']=="B10"),]
dim(data)

# omit 4 outliers from pca plots
drops = c("58828", "58852", "58956", "58921")
data = data[!(row.names(data) %in% drops),]

data.expr = data[1:nrow(data),3:ncol(data)]
data.expr[1:5, 1:5]

data.meta = data[1:nrow(data),1:2]
data.meta[1:5,1:2]

data.expr.t = t(data.expr)
data.expr.t[1:5, 1:5]

Groups = data.meta[,'Group']

library(edgeR)
d <- DGEList(counts=data.expr.t,group=factor(Groups))
dim(d)

### CPM > 1 for at least two samples
keep <- rowSums(cpm(d)>1) >= 2
d <- d[keep,]
dim(d)
d$samples$lib.size <- colSums(d$counts)
############################################

d <- calcNormFactors(d)
d1 <- estimateCommonDisp(d, verbose=T)
names(d1)
d2 <- estimateTagwiseDisp(d1)
names(d2)
plotBCV(d2)

data_ps_counts = t(d1$pseudo.counts)
data_ps_counts = cbind(data.meta, data_ps_counts)
data_ps_counts[1:5,1:5]
dim(data_ps_counts)

gene.sizes = read.table(paste(Data_Dir,"Mus_musculus.GRCm38.84.all_exons.merged.sizes.txt",sep=""), header=F)
row.names(gene.sizes) = gene.sizes$V1
dim(gene.sizes)
gene.sizes[1:5,]
gene.sizes = gene.sizes[rownames(d1$counts),]
dim(gene.sizes)
gene.sizes[1:5,]

data_rpkm = t(rpkm(d1, gene.length = gene.sizes[,2]))
data_rpkm = cbind(data.meta, data_rpkm)
data_rpkm[1:5,1:5]
dim(data_rpkm)

write.table(data_ps_counts,"ps_counts_b6b10_normalized_CPM1.txt",sep="\t",col.names=T,row.names=F,quote=F)
write.table(data_rpkm,"rpkm_b6b10_normalized_CPM1.txt",sep="\t",col.names=T,row.names=F,quote=F)

###### compute anova pval, fdr and DEG ######
data_ps_counts = read.table("ps_counts_b6b10_normalized_CPM1.txt",sep="\t",header=T,stringsAsFactors=F)
data_ps_b6 = data_ps_counts[which(data_ps_counts[,'Group']=="B6"),]
dim(data_ps_b6)

data_ps_b10 = data_ps_counts[which(data_ps_counts[,'Group']=="B10"),]
dim(data_ps_b10)

anova_pval_b6 = c()
for (i in 3:ncol(data_ps_b6)) {
	anova_pval_b6[i-2] = anova(lm(data_ps_b6[,i]~ as.factor(data_ps_b6[,"Strain"])))[1,5]
}

anova_pval_b10 = c()
for (i in 3:ncol(data_ps_b10)) {
	anova_pval_b10[i-2] = anova(lm(data_ps_b10[,i]~ as.factor(data_ps_b10[,"Strain"])))[1,5]
}

anova_pval_b6b10 = c()
for (i in 3:ncol(data_ps_counts)) {
	anova_pval_b6b10[i-2] = anova(lm(data_ps_counts[,i]~ as.factor(data_ps_counts[,"Strain"])))[1,5]
}

anova_fdr_b6 = p.adjust(anova_pval_b6, method ="fdr")
anova_fdr_b10 = p.adjust(anova_pval_b10, method ="fdr")
anova_fdr_b6b10 = p.adjust(anova_pval_b6b10, method ="fdr")

data_anova = data.frame(anova_pval_b6, row.names=colnames(data_ps_b6)[3:ncol(data_ps_b6)])
data_anova = cbind(data_anova,anova_fdr_b6)
data_anova = cbind(data_anova,anova_pval_b10)
data_anova = cbind(data_anova,anova_fdr_b10)
data_anova = cbind(data_anova,anova_pval_b6b10)
data_anova = cbind(data_anova,anova_fdr_b6b10)
data_anova[1:5,]

gene_names = read.table(paste(Data_Dir,"ensembl_symbol.txt",sep=""),header=T)
gene_names[1:5,]
data_anova = merge(gene_names, data_anova, by.x="Gene", by.y=0)
data_anova[1:5,]

write.table(data_anova,"anova_b6b10_normalized_CPM1.txt",sep="\t",col.names=T,row.names=F,quote=F)
###########################################################
