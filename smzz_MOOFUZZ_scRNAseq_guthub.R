### Load data ###
data.parsed = read.table("GSE62270_data_counts_Whole_Organoid_Replicate_1.txt", sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors=F)
dim(data.parsed) ##[23630, 288]

### converts into boolean matrix ###
ycol <- "Nb Expressed Genes [count > 0]"
expressed.genes.per.sample = data.parsed
expressed.genes.per.sample[expressed.genes.per.sample > 0] = 1

### filtering ###
##step 1 filtering: cell filtering ###
coldata.sum = colSums(expressed.genes.per.sample)
length(coldata.sum) 
data.parsed.cellfilt <- data.parsed[,which(coldata.sum>2000)]
dim(data.parsed.cellfilt) 
  
expressed.genes.per.sample.cellfilt<-expressed.genes.per.sample[,which(coldata.sum>2000)]
dim(expressed.genes.per.sample.cellfilt) 

##step 2 filtering: gene filtering ###
rowdata.sum = rowSums(expressed.genes.per.sample.cellfilt)
length(rowdata.sum) 
data.parsed.cellgenefilt <- data.parsed.cellfilt[(rowdata.sum>3 & rowdata.sum<500),]
dim(data.parsed.cellgenefilt) 

###step 3 filtering: top 3000 highest variant genes filtering ###
gvar <- apply(data.parsed.cellgenefilt, MARGIN=1, FUN=var, na.rm=TRUE)
length(gvar)

gvar.sorted<-sort(gvar, decreasing =TRUE)
length(gvar.sorted)

data.parsed.cellgenefilt.gvarsorted <- data.parsed.cellgenefilt[match(names(gvar.sorted)[1:3000], rownames(data.parsed.cellgenefilt)),]
dim(data.parsed.cellgenefilt.gvarsorted)

###row-wise zero-mean normalization ###
data.parsed.cellgenefilt.gvarsorted_scaleddat <- t(scale(t(data.parsed.cellgenefilt.gvarsorted)))
apply(data.parsed.cellgenefilt.gvarsorted_scaleddat, 1, mean)
apply(data.parsed.cellgenefilt.gvarsorted_scaleddat, 1, sd) 

#######################Scnorm normalization ######################
library(devtools)
library(SCnorm)

#############Required inputs#################
ExampleSimSCData<-matrix(as.numeric(unlist(data.parsed.cellgenefilt.gvarsorted)), ncol = ncol(data.parsed.cellgenefilt.gvarsorted), byrow = FALSE)
mode(ExampleSimSCData)
dim(ExampleSimSCData)
rownames(ExampleSimSCData)<-rownames(data.parsed.cellgenefilt.gvarsorted)
colnames(ExampleSimSCData)<-colnames(data.parsed.cellgenefilt.gvarsorted)

sampsize<-ncol(ExampleSimSCData)
ExampleSimSCData[1:5,1:5]
str(ExampleSimSCData)
ExampleSimSCData <- SingleCellExperiment::SingleCellExperiment(assays = list('counts' = ExampleSimSCData))
ExampleSimSCData <- as(ExampleSimSCData, "SingleCellExperiment")
Conditions = rep(c(1), each= sampsize)
head(Conditions)

################SCnorm: Check count-depth relationship############
library("stats")
setEPS()
postscript("check_GSE62270datacount_wh.organoid.rep1_eval.eps")
par(mfrow=c(2,2))
countDeptEst <- plotCountDepth(Data = ExampleSimSCData, Conditions = Conditions,
                               FilterCellProportion = .1, NCores=3)
dev.off()
str(countDeptEst)
head(countDeptEst[[1]])


ExampleSimSCData = SingleCellExperiment::counts(ExampleSimSCData)
ExampleSimSCData.CPM <- t((t(ExampleSimSCData) / colSums(ExampleSimSCData)) *
                            mean(colSums(ExampleSimSCData)))
countDeptEst.CPM <- plotCountDepth(Data = ExampleSimSCData,
                                   NormalizedData = ExampleSimSCData.CPM,
                                   Conditions = Conditions,
                                   FilterCellProportion = .1, NCores=3)
str(countDeptEst.CPM)
head(countDeptEst.CPM[[1]])

#############SCnorm: Normalization###########
library("stats")

Conditions = rep(c(1), each= sampsize)
setEPS()
postscript("NormalizedData_GSE62270datacount_wh.organoid.rep1_k_evaluation.eps")
par(mfrow=c(2,2))
DataNorm <- SCnorm(Data = ExampleSimSCData,
                   Conditions = Conditions,
                   PrintProgressPlots = TRUE,
                   FilterCellNum = 10, K =1,
                   NCores=3, reportSF = TRUE)

dev.off()
DataNorm
NormalizedData <- SingleCellExperiment::normcounts(DataNorm)##final scmnormalized data
NormalizedData[1:5,1:5]
write.table(NormalizedData,file="NormalizedData_GSE62270datacount_wh.organoid.rep1.csv",sep=",",row.names=TRUE,col.names=FALSE)

##########Evaluate choice of K ##
setEPS()
postscript("MyNormalizedData_GSE62270datacount_wh.organoid.rep1_scnorm.eps")
par(mfrow=c(2,2))
countDeptEst.SCNORM <- plotCountDepth(Data = ExampleSimSCData,
                                      NormalizedData = NormalizedData,
                                      Conditions = Conditions,
                                      FilterCellProportion = .1, NCores=3)
dev.off()

###Partitioning Cluster Analysis Using Fuzzy C-Means###
library(ppclust)
library(factoextra)
library(cluster)
library(fclust)

NormalizedData

objfuncs<-matrix(,nrow=0, ncol=4)
colnames(objfuncs)<-c("Fuzzy Silhouette Index","Partition Entropy","Partition Coefficient","Modified Partition Coefficient")

seedval<-c(120,156,245,346,476,543,674,768,876)
  
for(i in 2:10)
{
print(paste0("i=",i))
set.seed(seedval[i-1])

v0<-i 
res.fcm <- fcm(t(NormalizedData), centers=v0, alginitv="hartiganwong", alginitu="imembrand", m=2)
res.fcm4 <- ppclust2(res.fcm, "fclust")
idxsf <- SIL.F(res.fcm4$Xca, res.fcm4$U, alpha=1)
idxpe <- PE(res.fcm4$U)
idxpc <- PC(res.fcm4$U)
idxmpc <- MPC(res.fcm4$U)

objfuncs<-rbind(objfuncs,c(idxsf,idxpe,idxpc,idxmpc))
}
objfuncs<-as.data.frame(objfuncs)

write.table(objfuncs,file="GSE62270_datcounts_Wh.Org.Rep1_fcm_cl2to10_objval.csv",sep=",",row.names=FALSE,col.names=TRUE)


#######topsis: multiple-criteria decision making ########
library("topsis")

weights<-c(1,1,1,1)
criteriaMinMax1<-c("+","-","+","+") ###c(idxsf,idxpe,idxpc,idxmpc)
result.topsis1<-topsis(as.matrix(objfuncs), weights, criteriaMinMax1)

temp23<-result.topsis1$score
cs<-seq(1, 9, by=1)
setEPS()
postscript("NormDat_GSE62270dtct_wh.organoid.rep1_topsis_optscore_eachcls_plot.eps")
par(mar=c(8.1,4.1,4.1,2.1))
names(temp23)<-paste0("#Cluster=",(cs+1))
original.parameters<- par( no.readonly = TRUE )##to disble x-axis position labels
par(xaxt="n")##to disble x-axis position labels
plot(temp23, type="b",col = "red", xlab="", ylab="TOPSIS optimal score", main = "TOPSIS optimal score for each case study")
axis(1, at=cs)
text(cs, par("usr")[3]-0.085, srt=90, pos=1, labels=names(temp23), xpd=TRUE)
dev.off()

rownames(result.topsis1)<-names(temp23)
write.table(result.topsis1,file="GSE62270_datcounts_Wh.Org.Rep1_fcm_cl2to10_toposis.rank.csv",sep=",",row.names=FALSE,col.names=TRUE)

########final fuzzy clustering with computed optimized (best) cluster number#####
v0<-(result.topsis1[which(result.topsis1$rank==1),]$alt.row+1)    #number of evaluated optimized clusters
v0
res.fcm <- fcm(t(NormalizedData), centers=v0, alginitv="hartiganwong", alginitu="imembrand", m=2)
res.fcm
as.data.frame(res.fcm$u)##fuzzy membership degree matrix
dim(res.fcm$u)
setEPS()
postscript("NormDat_GSE62270dtct_wh.organoid.rep1_fcm_fuzzmemval_boxplot.eps")
par(mfrow=c(2,2))
boxplot(res.fcm$u,notch=TRUE,col=(c("gold","darkgreen")),
           main="Fuzzy membership value of the cells", xlab="Cell Clusters", ylab="Fuzzy membership score")
dev.off()

res.fcm$v##cluster centroid after clustering
dim(res.fcm$v)

summary(res.fcm)
res.fcm$cluster##final cluster information of the cells
res.fcm$csize#the number of cells for each cluster 
res.fcm$iter#the number of iteration during clustering 

temp11<-res.fcm$cluster
names(temp11)<-colnames(NormalizedData)
write.table(temp11,file="GSE62270_datcounts_Wh.Org.Rep1_fcm_cl2to10_toposis.result.cellclsinfo.csv",sep=",",row.names=TRUE,col.names=FALSE)

setEPS()
postscript("NormDat_GSE62270dtct_wh.organoid.rep1_finalfcm_clusterplot.eps")
res.fcm2 <- ppclust2(res.fcm, "kmeans")
factoextra::fviz_cluster(res.fcm2, data = t(NormalizedData),
                         ellipse.type = "convex",
                         palette = "jco",
                         repel = TRUE)
dev.off()

res.fcm4 <- ppclust2(res.fcm, "fclust")
idxsf <- SIL.F(res.fcm4$Xca, res.fcm4$U, alpha=1)
idxpe <- PE(res.fcm4$U)
idxpc <- PC(res.fcm4$U)
idxmpc <- MPC(res.fcm4$U)
cat("Fuzzy Silhouette Index: ", idxsf)
cat("Partition Entropy: ", idxpe)
cat("Partition Coefficient: ", idxpc)
cat("Modified Partition Coefficient: ", idxmpc)


optimum.objval<-c(idxsf,idxpe,idxpc,idxmpc)
names(optimum.objval)<-c("Fuzzy Silhouette Index","Partition Entropy","Partition Coefficient","Modified Partition Coefficient")
write.table(t(as.matrix(optimum.objval)),file="GSE62270_datcounts_Wh.Org.Rep1_fcm_cl2to10_toposis.result.fianlclsvaldityind.csv",sep=",",row.names=FALSE,col.names=TRUE)



#################DEG analysis###############
#########Differential expression for each subtype vs rest with limma#################
temp11
cl

ord_cellclsinfo<-sort(temp11)
ord_NormData<-NormalizedData[,match(names(ord_cellclsinfo), colnames(NormalizedData))]

sample_clslabel_cl1vsrest_r <- rep(c("Cluster1","Cluster2"),c(res.fcm$csize[1], res.fcm$csize[2]))
sample_clslabel_cl1vsrest_r
design_cl1vsrest_r<-model.matrix(~ sample_clslabel_cl1vsrest_r- 1)##when 1st class-label is experimental, them use this design (i.e.,"-1")
head(design_cl1vsrest_r,n=15)
colnames(design_cl1vsrest_r) <- c("Cluster1","Cluster2")##use when u used "-1" in model.matrix

###zero-mean standardization####
ord_NormData_scaleddat <- t(scale(t(ord_NormData)))
apply(ord_NormData_scaleddat, 1, mean)
apply(ord_NormData_scaleddat, 1, sd) 

####omit the chromosome from genename of the data###
normdat_gnameonly<-sapply(strsplit(rownames(ord_NormData_scaleddat), "_"), head, 1)
rownames(ord_NormData_scaleddat)<-normdat_gnameonly

############Testing for differential expression for rna###
library(limma)
fit_cl1vsrest_r <- lmFit(ord_NormData_scaleddat,design_cl1vsrest_r)
names(fit_cl1vsrest_r)
contrast.matrix_cl1vsrest_r <- makeContrasts("Cluster1-Cluster2", levels = design_cl1vsrest_r)
contrast.matrix_cl1vsrest_r
fit_cl1vsrest_r <- contrasts.fit(fit_cl1vsrest_r, contrast.matrix_cl1vsrest_r)
fit_ebayes_cl1vsrest_r <- eBayes(fit_cl1vsrest_r)
dim(fit_ebayes_cl1vsrest_r)

options(digits=3)
top2_cl1vsrest_r <- topTable(fit_ebayes_cl1vsrest_r,number=Inf,adjust="bonferroni",sort.by="P")###use this if use "makeContrasts" function
sum(top2_cl1vsrest_r$adj.P.Val<0.05)

DEGlist_cl1vsrest_r<-top2_cl1vsrest_r[which(top2_cl1vsrest_r$adj.P.Val<0.05),]
DEGlist_srt_cl1vsrest_r<-DEGlist_cl1vsrest_r[order(DEGlist_cl1vsrest_r$adj.P.Val, decreasing = FALSE),]
fDEGlist_adjpvalFCfilt_cl1vsrest_r<-DEGlist_srt_cl1vsrest_r
nrow(fDEGlist_adjpvalFCfilt_cl1vsrest_r)

limmapvalsiggn_info<-cbind(as.matrix(rownames(fDEGlist_adjpvalFCfilt_cl1vsrest_r)),fDEGlist_adjpvalFCfilt_cl1vsrest_r[2:ncol(fDEGlist_adjpvalFCfilt_cl1vsrest_r)])
write.table(limmapvalsiggn_info,file="GSE62270_datcounts_Wh.Org.Rep1_optifcm.cl1vscl2cell.DEGinfo.csv",sep=",",row.names=FALSE,col.names=TRUE)


###adjusted p-value plot###
setEPS()
postscript("NormDat_GSE62270dtct_wh.organoid.rep1_optifcm_DEGrankwise_adjpval_plot.eps")
plot(limmapvalsiggn_info$adj.P.Val,xlab="Rank",ylab="Adjusted p-value",main="Rankwise adjusted p-value",col="red")
dev.off()

###############membership plot of first 50 cells after optimized fcm###
noex12<-50
cs12<-seq(1, noex12, by=1)
setEPS()
postscript("NormDat_GSE62270dtct_wh.organoid.rep1_optimizedfcm_firstfewcells_clustmemplot.eps")
par(mar=c(6.1,4.1,4.1,2.1))
plot(res.fcm$u[1:noex12,1],type = "o", frame = FALSE, pch = 19, 
     col = "red", main=paste0("Fuzzy membership values of first ",noex12," cells of two clusters"), xlab = "Cells", ylab = "Fuzzy membership value",ylim=c(0, 1),xaxt = 'n')

lines(res.fcm$u[1:noex12,2],col = "blue", type = "o", lty = 2)
axis(1, at=cs12,rownames(res.fcm$u)[1:noex12],labels = FALSE)
text(cs12, par("usr")[3]-0.045, srt=45, pos=1, labels=rownames(res.fcm$u)[1:noex12], xpd=TRUE,cex=0.5)
legend("topright", legend=c("Cluster 1", "Cluster 2"),
       col=c("red", "blue"), lty=1:2, cex=0.8)
dev.off()

