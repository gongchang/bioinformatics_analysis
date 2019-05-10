library(edgeR)
library(limma)
library(gplots)

#load in data
##the input data are count data that are usually either from Kallisto or Salmon
df=read.table("salmon_gene_counts.txt",sep="\t",header = TRUE, row.names = 1)
df2=read.table("samples_design.csv",sep=",",header=TRUE)


y=DGEList(df)

######y

#filtering to remove low counts
#a gene is only retained if it is expressed at a count-per-million (CPM) above 0.5 in at least 9 samples.
keep <- rowSums(cpm(y) > 0.5) >= 9
y <- y[keep,]

#cpm cutoff can be obtained this way:
#> cpm(10, mean(y$samples$lib.size))
#[,1]
#[1,] 0.3755316

######
#TMM normalization
y <- calcNormFactors(y)

#build the design matrix
condition=as.factor(df2$condition)
dose=as.numeric((df2$Dose))
set=as.factor(df2$set)
schedule=as.factor((df2$Schedule))

#help("model.matrix")
design <- model.matrix(~0+condition+dose+set+schedule)
#limma-voom
y.voom <- voom(y,design,plot=TRUE)
fit <- lmFit(y.voom, design)
fit <- eBayes(fit)

#contrast
Treat.vsControl <- makeContrasts(conditiontreatment-conditioncontrol,levels=design)
conFit <- contrasts.fit(fit,Treat.vsControl)

#adjust the p.value for differential expression (contrast=treatment-control)
#decrease expression in treatment
adjFit=topTable(conFit,number=nrow(conFit),adjust.method = "BH")
diff.expr=adjFit[((adjFit["logFC"]< (-1)) & (adjFit["adj.P.Val"]<0.05)),]


#adjust the p.value for dose
adjDose=topTable(fit,coef = "dose",number = nrow(fit),adjust.method = "BH")
dose.depend=adjDose[(adjDose["logFC"] < 0 & (adjDose["adj.P.Val"]<0.05)),]

#differential expressed genes (decrease in treatment) and dose-dependent
de.dose.depend=intersect(row.names(diff.expr),row.names(dose.depend))
m=y.voom$E[de.dose.depend,]  #the expression values

#differential expressed genes (increase in treatment) and dose-dependent
diff.expr.inc=adjFit[((adjFit["logFC"]> 1) & (adjFit["adj.P.Val"]<0.05)),]
dose.depend.inc=adjDose[(adjDose["logFC"] > 0 & (adjDose["adj.P.Val"]<0.05)),]
de.dose.depend.inc=intersect(row.names(diff.expr.inc),row.names(dose.depend.inc))

##modeling the 5-day data using the full model
df2_qdx5=df2[(df2$Schedule=="qdx5" & df2$Dose!=100),]
df_qdx5=df[,df2_qdx5$Id]

y_qdx5=DGEList(df_qdx5)
keep_qdx5 <- rowSums(cpm(y_qdx5) > 0.5) >= 9
y_qdx5 <- y_qdx5[keep_qdx5,]
y_qdx5 <- calcNormFactors(y_qdx5)

#build the design matrix
condition_qdx5=as.factor(df2_qdx5$condition)
dose_qdx5=as.numeric((df2_qdx5$Dose))
set_qdx5=as.factor(df2_qdx5$set)

design_qdx5 <- model.matrix(~0+condition_qdx5+dose_qdx5+set_qdx5)

y_qdx5.voom <- voom(y_qdx5,design_qdx5,plot=TRUE)
fit_qdx5 <- lmFit(y_qdx5.voom, design_qdx5)
fit_qdx5 <- eBayes(fit_qdx5)

#adjust the p.value for dose
adjDose_qdx5=topTable(fit_qdx5,coef = "dose_qdx5",number = nrow(fit_qdx5),adjust.method = "BH")


##get the differential expressed genes in qdx5
condition_qdx5=as.factor(df2_qdx5$condition)
dose_qdx5=as.numeric((df2_qdx5$Dose))
set_qdx5=as.factor(df2_qdx5$set)

design_qdx5_de <- model.matrix(~0+condition_qdx5+dose_qdx5+set_qdx5)

y_qdx5.voom.de <- voom(y_qdx5,design_qdx5_de,plot=TRUE)
fit_qdx5.de <- lmFit(y_qdx5.voom.de, design_qdx5_de)
fit_qdx5.de <- eBayes(fit_qdx5.de)

Treat.vsControl.qdx5 <- makeContrasts(condition_qdx5treatment-condition_qdx5control,levels=design_qdx5_de)
conFit.qdx5 <- contrasts.fit(fit_qdx5.de,Treat.vsControl.qdx5)
adjFit.qdx5=topTable(conFit.qdx5,number=nrow(conFit.qdx5),adjust.method = "BH")

diff.expr.qdx5=adjFit.qdx5[((adjFit.qdx5["logFC"]< (-1)) & (adjFit.qdx5["adj.P.Val"]<0.05)),]
diff.expr.qdx5.inc=adjFit.qdx5[((adjFit.qdx5["logFC"]> (1)) & (adjFit.qdx5["adj.P.Val"]<0.05)),]

#dose-dependent in qdx5 and decrase in treatment (in both qdx5 and all)
dose.depend.qdx5=adjDose_qdx5[(adjDose_qdx5["logFC"] < 0 & (adjDose_qdx5["adj.P.Val"]<0.05)),]
de.dose.depend.qdx5=intersect(row.names(diff.expr.qdx5),row.names(dose.depend.qdx5))
de.dose.depend.qdx5=intersect(de.dose.depend.qdx5,row.names(diff.expr))
length(de.dose.depend.qdx5)
m_qdx5=y_qdx5.voom$E[de.dose.depend.qdx5,]  #the expression values

#dose-dependent in qdx5 and increase in treatment (in both qdx5 and all)
dose.depend.qdx5.inc=adjDose_qdx5[(adjDose_qdx5["logFC"] > 0 & (adjDose_qdx5["adj.P.Val"]<0.05)),]
de.dose.depend.qdx5.inc=intersect(row.names(diff.expr.inc),row.names(dose.depend.qdx5.inc))
de.dose.depend.qdx5.inc=intersect(de.dose.depend.qdx5.inc,row.names(diff.expr.qdx5.inc))
length(de.dose.depend.qdx5.inc)

#zero rows if full model is used

########################
###########modeling the 28 day treatnebt with the full model
####qdx28

df=read.table("salmon_gene_counts.txt",sep="\t",header = TRUE, row.names = 1)
df2=read.table("samples_design.csv",sep=",",header=TRUE)

# ##
df2=df2[(df2$Schedule=="qdx29"),]
df=df[,as.vector(df2$Id)]

y=DGEList(df)

######y

#filtering to remove low counts
#a gene is only retained if it is expressed at a count-per-million (CPM) above 0.5 in at least 9 samples.
keep <- rowSums(cpm(y) > 0.5) >= 9
y <- y[keep,]

#cpm cutoff can be obtained this way:
#> cpm(10, mean(y$samples$lib.size))
#[,1]
#[1,] 0.3755316

######
#TMM normalization
y <- calcNormFactors(y)

#build the design matrix
condition=as.factor(df2$condition)
dose=as.numeric((df2$Dose))
set=as.factor(df2$set)
schedule=as.factor((df2$Schedule))
#help("model.matrix")
design <- model.matrix(~0+condition+dose+set+schedule)
#limma-voom
y.voom <- voom(y,design,plot=TRUE)
fit <- lmFit(y.voom, design)
fit <- eBayes(fit)
#contrast
Treat.vsControl <- makeContrasts(conditiontreatment-conditioncontrol,levels=design)
conFit <- contrasts.fit(fit,Treat.vsControl)
#adjust the p.value for differential expression (contrast=treatment-control)
#decrease expression in treatment
adjFit.qdx28=topTable(conFit,number=nrow(conFit),adjust.method = "BH")
diff.expr.qdx28=adjFit.qdx28[((adjFit.qdx28["logFC"]< (-1)) & (adjFit.qdx28["adj.P.Val"]<0.05)),]
########################

#####
genes.down.in.qdx5.and.all=intersect(rownames(diff.expr),rownames(diff.expr.qdx5))
genes.down=intersect(genes.down.in.qdx5.and.all,rownames(diff.expr.qdx28))
down.df=data.frame(diff.expr[genes.down,][,c("logFC","adj.P.Val")],diff.expr.qdx5[genes.down,][,c("logFC","adj.P.Val")],diff.expr.qdx28[genes.down,][,c("logFC","adj.P.Val")])
colnames(down.df)=c("logFC-all","adj.pval-all","logFC-qdx5","adj.pval-qdx5","logFC-qdx29","adj.pval-qdx29")
write.csv(down.df,file="GE0375_down_regulation.csv",quote = FALSE)
