library("gplots")
library(RColorBrewer)
library(edgeR)
library(WGCNA)
library(flashClust)
library(DESeq2)

exprs <- as.matrix(read.table("RNA.seq.144.coga.inia.count.updated.header.txt", header=TRUE, sep=" ", row.names=1,as.is=TRUE))
cova = read.csv("coga.inia.detailed.pheno.04.12.17.csv", header=T)
rownames (cova) =cova$IID
table(rownames(cova)==colnames(exprs))
remove.samples <- rownames(cova)%in%c("X232","X459","X98","X442.1")
gExprOut <- exprs[,!remove.samples]
datTrait <- cova[!remove.samples,]
gExpr <- DGEList(gExprOut)
gExpr <- calcNormFactors(gExpr)
design <- model.matrix( ~ RNAsequencedby + RIN +PM.+ Gender +Age + Alc_status, datTrait)

#dds = DESeqDataSetFromMatrix( gExprOut, datTrait, design = ~ RNAsequencedby + RIN +PM.+ Gender +Age + alcohol_intake_gmsperday, datTrait)  ##for alcohol consumption deseq

vobjGenes <- voom(gExpr, design )
design <- model.matrix(~RNAsequencedby + RIN +PM.+ Gender +Age + Alc_status,  data=datTrait)
batch.design <- design[,1:2]
cov.design <- design[,-(1:2)]
vsd.removed=removeBatchEffect(vobjGenes, batch=datTrait$RNAsequencedby, covariates=cov.design)
write.table (vsd.removed, "batch.age.rin.sex.pm.alc.corrected.coga.inia.expression.txt", row.names=TRUE, quote=FALSE)

#sed 's/X214/id X214/' batch.age.rin.sex.pm.alc.corrected.coga.inia.expression.txt > temp
#mv temp batch.age.rin.sex.pm.alc.corrected.coga.inia.expression.txt

#head -n 1  batch.age.rin.sex.pm.alc.corrected.coga.inia.expression.txt | tr "\t" "\n" | sed 's/\./_/g' | sed 's/X//g' | awk 'FNR >1 {print $1,$1}' > coga.inia.final.subject.id.txt
##sample number 277_1 is coded as 277 and 89.61 as 89.610 in genotype file..make sure to change in "coga.inia.final.subject.id.txt"

#bash ~/join.unsorted -1 1 -2 1 ensembl.id.60k.genes.txt mart_export-2.txt | awk 'BEGIN {print "geneid chr s1 s2"} {print $1,"chr"$2,$3,$4}' | tr " " "\t" > ensembl.id.60k.gene.location.txt

#coga_data/COGA_INIA
#plink --bfile COGA_BB_Imputed --keep /hpc/users/kapoom02/work_dir/coga.inia.final.subject.id.txt --indiv-sort f /hpc/users/kapoom02/work_dir/coga.inia.final.subject.id.txt --make-bed --out /hpc/users/kapoom02/6/INIA_RNAseq/COGA_BB_Impute

#cd /hpc/users/kapoom02/6/INIA_RNAseq
#bash plink.recode.transpose.sh


exprs <- as.matrix(read.table("RNA.seq.144.coga.inia.count.updated.header.txt", header=TRUE, sep=" ", row.names=1,as.is=TRUE))
cova = read.csv("coga.inia.detailed.pheno.04.12.17.csv", header=T)
rownames (cova) =cova$IID
table(rownames(cova)==colnames(exprs))
remove.samples <- rownames(cova)%in%c("X232","X459","X98","X442.1")
gExprOut <- exprs[,!remove.samples]
datTrait <- cova[!remove.samples,]
gExpr <- DGEList(gExprOut)
gExpr <- calcNormFactors(gExpr)
keeps.na.out <- rownames(datTrait[!is.na(datTrait$alcohol_intake_gmsperday),])

remove.samples <- rownames(cova)%in%keeps.na.out
gExprOut <- exprs[,remove.samples]
datTrait <- cova[remove.samples,]
table(rownames(datTrait)==colnames(gExprOut))

datTrait$nBMI = as.numeric(datTrait$BMI)
keeps.na.out <- rownames(datTrait[!is.na(datTrait$nBMI),])
remove.samples <- rownames(datTrait)%in%keeps.na.out
gExprOut <- gExprOut[,remove.samples]
datTrait <- datTrait[remove.samples,]
table(rownames(datTrait)==colnames(gExprOut))

keeps.na.out <- rownames(datTrait[!is.na(datTrait$Liver_class),])
remove.samples <- rownames(datTrait)%in%keeps.na.out
gExprOut <- gExprOut[,remove.samples]
datTrait <- datTrait[remove.samples,]
table(rownames(datTrait)==colnames(gExprOut))

keeps.na.out <- rownames(datTrait[!is.na(datTrait$Total_drinking_yrs),])
remove.samples <- rownames(datTrait)%in%keeps.na.out
gExprOut <- gExprOut[,remove.samples]
datTrait <- datTrait[remove.samples,]
table(rownames(datTrait)==colnames(gExprOut))


gExpr <- DGEList(gExprOut)
gExpr <- calcNormFactors(gExprOut)

dds = DESeqDataSetFromMatrix( gExprOut, datTrait, design = ~ RNAsequencedby + RIN +PM.+ nBMI+ Brain_pH + Total_drinking_yrs + Frozentissue + Liver_class + AUDIT+ Agonal_phase + Alc_con)
dds<-DESeq(dds, parallel=TRUE)
res_AC <- results(dds)
write.table (res_AC, "diff.exp.alc.consumption.txt", quote=F)

dds = DESeqDataSetFromMatrix( gExprOut, datTrait, design = ~ RNAsequencedby + RIN +PM.+ Gender +Age + AUDIT)
dds<-DESeq(dds)
res_audit <- results(dds)
write.table (res_audit, "diff.exp.AUDIT.txt", quote=F)