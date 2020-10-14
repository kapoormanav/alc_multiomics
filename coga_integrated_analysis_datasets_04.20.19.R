## import_multiple_files_to_R
# Purpose: Import multiple files to the Global Environment in R

# set working directory
#setwd("/Users/manavkapoor/Box/alcohol.all.smr/")
#setwd("/Users/manavkapoor/Documents/smr_card_starnet/")

# list all csv files from the current directory
list.files(pattern=".smr$") # use the pattern argument to define a common pattern  for import files with regex. Here: .csv

# create a list from these files
#list.filenames<-list.files(pattern="CARDIOGENICS_MP.smr$") #combined_METAANALYSIS
list.filenames<-list.files(pattern=".smr$")
list.filenames

# create an empty list that will serve as a container to receive the incoming files
list.data<-list()

# create a loop to read in your data
for (i in 1:length(list.filenames))
{
  list.data[[i]]<-read.table(list.filenames[i], header=T)
}

# add the names of your data to the list

newnames = gsub("combined_", "", list.filenames)
newnames = gsub("chrall.CPRA_b37.", "", newnames)
newnames = gsub(".smr", "", newnames)



for (i in 1:length(newnames))
{
  names(list.data[[i]]) <- newnames[i]
}

## sort by P value
sorted.datset.list <- lapply(list.data, function(x) x[order(x$p_SMR), ])
substring.list <- lapply(seq_along(sorted.datset.list), 
                         function(i, x){
                           x[[i]]$Gene <- sub(":ILMN_.*", "", x[[i]]$Gene)
                           return (x[[i]])
                         }, sorted.datset.list
)
sorted.datset.list <- lapply(seq_along(substring.list), 
                             function(i, x){
                               x[[i]]$Gene <- sub(".*:", "", x[[i]]$Gene)
                               return (x[[i]])
                             }, substring.list
)
## Remove duplicates
deduplicated.datset.list <- lapply(sorted.datset.list, function(x) x[!duplicated(x$Gene), ])
##Limit to protein coding only genes
#deduplicated.datset.list <- lapply(deduplicated.datset.list, function(x) x[x$Gene %in% protein_coding_genes$Gene,])
##Calulate FDR
fdr.datset.list <- lapply(seq_along(deduplicated.datset.list), 
                          function(i, x){
                            x[[i]]$fdr <- p.adjust(x[[i]]$p_SMR, method = "fdr", n = length(x[[i]]$p_SMR))
                            return (x[[i]])
                          }, deduplicated.datset.list
)

microglia.expressed.genes <- read.table("~/Dropbox/microglia.expressed.genes.txt", quote="\"", comment.char="")
deduplicated.datset.list <- lapply(fdr.datset.list, function(x) x[x$Gene %in% microglia.expressed.genes$V1,])

###make plots

plot_list = list()
for (i in 1:length(fdr.datset.list))
{
  df1 = as.data.frame(fdr.datset.list[[i]])
  id.1 = list.filenames[i]
  df1 <- df1[c(5,6,7,19,17,17,3,22)]
  colnames(df1) <- c("snp", "chrom", "bp","pvalue", "beta","or","gene", "fdr")
  df1.sig <- subset(df1,fdr <= 0.1)
  fdr.cut = -log10(max(df1.sig$pvalue))
  p1 <- ggman(df1, snp = "gene", bp = "bp", chrom = "chrom", pvalue = "pvalue", sigLine = fdr.cut, relative.positions = TRUE, title = id.1) + theme_bw()
  p2 <- ggmanLabel((p1+scale_color_manual(values=c("darkblue","darkred"))) , labelDfm = df1.sig, snp = "gene", label = "gene")
  plot_list[[i]] = p2
  #dev.off()
  
}

plot_list1 = list()
for (i in 1:length(fdr.datset.list))
{
  df1 = as.data.frame(fdr.datset.list[[i]])
  id.1 = newnames[i]
  df1 <- df1[c(5,6,4,19,17,17,3,22)]
  colnames(df1) <- c("snp", "chrom", "bp","pvalue", "beta","or","gene", "fdr")
  df1.sig <- subset(df1,fdr <= 0.1)
  fdr.cut = -log10(max(df1.sig$pvalue))
  p1 <- ggman(df1, snp = "bp", bp = "bp", chrom = "chrom", pvalue = "pvalue", sigLine = fdr.cut, relative.positions = TRUE, title = id.1) + theme_bw()
  p2 <- ggmanLabel((p1+scale_color_manual(values=c("darkblue","darkred"))) , labelDfm = df1.sig, snp = "bp", label = "gene")
  plot_list1[[i]] = p2
  #dev.off()
  
}

plot_list3 = list()
for (i in 1:length(fdr.datset.list))
{
  df1 = as.data.frame(fdr.datset.list[[i]])
  id.1 = newnames[i]
  df1 <- df1[c(5,6,4,19,17,17,3,22)]
  colnames(df1) <- c("snp", "chrom", "bp","pvalue", "beta","or","gene", "fdr")
  df1.sig <- subset(df1,fdr <= 0.05)
  fdr.cut = -log10(max(df1.sig$pvalue))
  p1 <- ggman(df1, snp = "bp", bp = "bp", chrom = "chrom", pvalue = "pvalue", sigLine = fdr.cut, relative.positions = TRUE, title = id.1) + theme_bw()
  p2 <- ggmanLabel((p1+scale_color_manual(values=c("darkblue","darkred"))) , labelDfm = df1.sig, snp = "bp", label = "gene")
  plot_list3[[i]] = p2
  #dev.off()
  
}
names(plot_list) <- newnames
names(plot_list1) <- newnames
names(plot_list3) <- newnames
names(fdr.datset.list) <- newnames

pdf("~/Dropbox/smr.manhattan.plots.emma_sch_with_brain_meta.pdf", width = 14, height = 10)
for (i in 1:length(fdr.datset.list)) {
  print(plot_list[[i]])
}
dev.off()

names(plot_list) <- newnames

#########################annotations

###Keep three columns only before merging

fdr.datset.list3 <- lapply(fdr.datset.list, function(x) x[(names(x) %in% c("Gene", "p_SMR", "fdr", "b_SMR", "p_HEIDI", "p_GWAS", "topSNP", "topSNP_chr", "topSNP_bp"))])

###Merge all
#library (dplyr)
merged.mar.igap.mmeta = fdr.datset.list3 %>%
  Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Gene"), .)
v1=1
for (i in 1:18){    #######this is the number of files
  c1 = newnames[i]
  c0 = paste(c1, "gwas", sep = ".")
  c2 <- paste(c1, "beta", sep = ".")
  c3 <- paste(c1, "p", sep = ".")
  c4 <- paste(c1, "heidi", sep = ".")
  c5 <- paste(c1, "fdr", sep = ".")
  c6 <- paste(c1, "topSNP", sep = ".")
  c7 <- paste(c1, "topSNP_chr", sep = ".")
  c8 <- paste(c1, "topSNP_bp", sep = ".")
  v1 <- v1+1
  v2 <- v1+1
  v3 <- v2+1
  v4 <- v3+1
  v5 <- v4+1
  v6 <- v5+1
  v7 <- v6+1
  v8 <- v7+1
  colnames(merged.mar.igap.mmeta)[v1] <- c6
  colnames(merged.mar.igap.mmeta)[v2] <- c7
  colnames(merged.mar.igap.mmeta)[v3] <- c8
  colnames(merged.mar.igap.mmeta)[v4] <- c0
  colnames(merged.mar.igap.mmeta)[v5] <- c2
  colnames(merged.mar.igap.mmeta)[v6] <- c3
  colnames(merged.mar.igap.mmeta)[v7] <- c4
  colnames(merged.mar.igap.mmeta)[v8] <- c5

  v1=v8
}
##melt dataset for dado
#for_dado = melt(merged.mar.igap.mmeta)

#write.table(for_dado, "~/Dropbox/igap1.st1.marioni4.pgc.mono.meta.mfc.longformat.with.beta.txt", row.names=F, quote=F)

##Concatanate datasets and create unique chr.bp.gene identifiers for non-missing genes in all datasets
new = do.call("rbind", fdr.datset.list)
new <- new[c(3,4,2)]
new = new[!duplicated(new$Gene), ]
new$gene.id <- paste(new$ProbeChr, new$Probe_bp, new$Gene, sep=".")
save( merged.mar.igap.mmeta,plot_list,new, file="ukbb.alcohol.rda",compress=T)

#########
#library paths for Minerva
###
#> .libPaths()
#"/hpc/packages/minerva-common/rpackages/3.5.1/site-library"
#"/hpc/packages/minerva-common/rpackages/bioconductor/3.7"
# "/hpc/users/kapoom02/.Rlib"
#"/hpc/packages/minerva-common/R/3.5.1/lib64/R/library"

Results.coloc.M4.MP <- read.csv("~/Box/Alzhiemer/SMR_analysis/Results.coloc.marioni4_AD_2018----MP", sep="")
Results.coloc.M4.MP$PP.H3.pl.H4 <- Results.coloc.M4.MP$PP.H3 + Results.coloc.M4.MP$PP.H4
Results.coloc.M4.MP.V1 = Results.coloc.M4.MP[c(9,5,6,10)]

#) %>%
formatRound( c('PP.H3.pl.H4', 'PP.H4.abf', 'PP.H3.abf', 'beta', 'or'), 3) %>%
  formatSignif(c('pvalue', 'fdr'), 3)
df1$gene1 = sub(":ILMN_.*", "", df1$gene)
df1$gene1 = sub(".*:", "", df1$gene1)

merged.mar.igap.mmeta.s <- filter_at(merged.mar.igap.mmeta, vars(matches("Liu2019drnkwk.CARDIOGENICS_MP.fdr")), any_vars(. <= 0.05)) %>%
  select(ends_with("p"), starts_with("Gene"))

merged.mar.igap.mmeta.s <- filter_at(merged.mar.igap.mmeta.s, vars(-matches(".fdr"), -ends_with("p"),-starts_with("Gene"), -ends_with("beta")), all_vars(. > 0.05)) %>%
  select(ends_with("p"),starts_with("Gene"))

merged.mar.igap.mmeta.s <- merge(merged.mar.igap.mmeta.s, new, by="Gene", all.x=TRUE)
ndata = length(merged.mar.igap.mmeta.s)
start_n = ndata
end_n = ndata - 3
merged.mar.igap.mmeta.s <- merged.mar.igap.mmeta.s[c(start_n,2:end_n)] ## n+4, 2:n+1
#merged.mar.igap.mmeta.s <- merged.mar.igap.mmeta.s[c(6,2:3)]
melted.merged.mar.igap.mmeta.s <- melt(merged.mar.igap.mmeta.s)
colnames(melted.merged.mar.igap.mmeta.s)[3] <- "SMR_P"
colnames(melted.merged.mar.igap.mmeta.s)[2] <- "Dataset"
colnames(melted.merged.mar.igap.mmeta.s)[1] <- "Gene"
p1 <- ggplot(melted.merged.mar.igap.mmeta.s, aes(Dataset, Gene)) + geom_tile(aes(fill = -log10(SMR_P)), colour = "white") + scale_fill_gradientn(colours=c("white", "orange", "red"))+
  theme(axis.text=element_text(size=rel(0.75), colour="black"))

pdf("~/Dropbox/mvp.pgc.coga.smr.macrophage.heatmap.pdf", width = 12, height = 12)
print(p1)
dev.off()

##barplots

Heat_Map_2 <- read_excel("~/Downloads/Heat Map_2.xlsx", 
                         +     sheet = "33 vs 44M-")

library(readxl)
saima_data <- read_excel("~/Dropbox/saima_data2.xlsx", 
                          sheet = "Ko vs wt M+")
p1 = ggplot(saima_data,aes(x=saima_data$Pathway,y=NES,fill=FDR))+geom_bar(stat="identity")+coord_flip()+scale_fill_gradient(low = "red", high = "orange")+ ylim(-4.5, 4.5)

pdf("~/Dropbox/heatmap2_Ko_vs_wt_M+saima.pdf", width = 10, height = 12)
print(p1)
dev.off()

#AUD data heidi > 0.05, FDR < 20%
merged.mar.igap.mmeta.pvalues <- filter_at(merged.mar.igap.mmeta, vars(ends_with(".fdr"), -matches("Liu")), any_vars(. <= 0.20)) %>%
  dplyr::select(ends_with(".p"),ends_with(".heidi"), ends_with(".gwas"), -matches("CARDIOGENICS"),-matches("STARNET"),-matches("MESA"), -matches("Liu"), starts_with("Gene"))

merged.mar.igap.mmeta.pvalues <- filter_at(merged.mar.igap.mmeta.pvalues, vars(ends_with(".gwas")), any_vars(. <= 0.00005)) %>%
  dplyr::select(ends_with(".p"),ends_with(".heidi"), ends_with("topSNP"), ends_with("topSNP_chr"), ends_with("topSNP_bp"),-matches("CARDIOGENICS"),-matches("STARNET"),-matches("Liu"), starts_with("Gene"))

merged.mar.igap.mmeta.pvalues <- filter_at(merged.mar.igap.mmeta.pvalues, vars(ends_with(".heidi")), any_vars(. > 0.0)) %>%
  dplyr::select(ends_with(".p"), -matches("CARDIOGENICS"),-matches("STARNET"),-matches("MESA"), -matches("Liu"), starts_with("Gene"))

merged.mar.igap.mmeta.pvalues <- filter_at(merged.mar.igap.mmeta.pvalues, vars(ends_with(".p")), any_vars(. >0)) %>%
  dplyr::select(ends_with(".p"),-matches("CARDIOGENICS"),-matches("STARNET"),-matches("MESA"), -matches("Liu"), starts_with("Gene"))

merged.mar.igap.mmeta.pvalues <- merge(merged.mar.igap.mmeta.pvalues, new, by="Gene", all.x=TRUE)

merged_filtered_data <-apply(merged.mar.igap.mmeta.pvalues[,c(1:4)], MARGIN=2, function(x) {ifelse(x<=0.05, 1,0)})
merged_filtered_data = as.data.frame(merged_filtered_data)
merged_filtered_data = dplyr::select(merged_filtered_data, matches("METAANALYSIS_MVP"))
merged_filtered_data[is.na(merged_filtered_data)] <- 0

UpSetR::upset(merged_filtered_data, nsets = 8, nintersects = NA, keep.order = F, empty.intersections = "on", order.by = "freq", 
      queries = list(list(query = UpSetR::intersects, params = list("METAANALYSIS_MVP_AUD_PGC_COGA1.rsids.Brain_m_Meta.p", "METAANALYSIS_MVP_AUD_PGC_COGA1.rsids.FB_Brain_mQTL_rsid.p"), color = "red", active = T), 
                      list(query = UpSetR::intersects, params = list("METAANALYSIS_MVP_AUD_PGC_COGA1.rsids.Brain_m_Meta.p", "METAANALYSIS_MVP_AUD_PGC_COGA1.rsids.FB_Brain_mQTL_rsid.p", "METAANALYSIS_MVP_AUD_PGC_COGA1.rsids.eqtl_fetal_brain_brien_et_al_query_format_hg19.p", "METAANALYSIS_MVP_AUD_PGC_COGA1.tbl.meta_Psyc_COGA_ROSMAP.p"), color = "green", active = T), 
                        list(query = UpSetR::intersects, params = list("METAANALYSIS_MVP_AUD_PGC_COGA1.rsids.eqtl_fetal_brain_brien_et_al_query_format_hg19.p", "METAANALYSIS_MVP_AUD_PGC_COGA1.tbl.meta_Psyc_COGA_ROSMAP.p"), color = "orange", active = T)))

merged_filtered_data <-apply(merged.mar.igap.mmeta.pvalues[,c(1:4)], MARGIN=2, function(x) {ifelse(x<=0.05, 1,0)})
merged_filtered_data = as.data.frame(merged_filtered_data)
merged_filtered_data = dplyr::select(merged_filtered_data, matches("Liu2019"))
merged_filtered_data[is.na(merged_filtered_data)] <- 0

UpSetR::upset(merged_filtered_data, nsets = 8, nintersects = NA, keep.order = F, empty.intersections = "on", order.by = "freq", text.scale =1.7,
      queries = list(list(query = UpSetR::intersects, params = list("Liu2019drnkwk.rsid.Brain_m_Meta.p", "Liu2019drnkwk.FB_Brain_mQTL_rsid.p"), color = "red", active = T), 
                     list(query = UpSetR::intersects, params = list("Liu2019drnkwk.rsid.Brain_m_Meta.p", "Liu2019drnkwk.FB_Brain_mQTL_rsid.p", "Liu2019drnkwk.eqtl_fetal_brain_brien_et_al_query_format_hg19.p", "Liu2019drnkwk.meta_Psyc_COGA_ROSMAP.p"), color = "green", active = T), 
                     list(query = UpSetR::intersects, params = list("Liu2019drnkwk.eqtl_fetal_brain_brien_et_al_query_format_hg19.p", "Liu2019drnkwk.meta_Psyc_COGA_ROSMAP.p"), color = "orange", active = T)))

#DPW data heidi > 0.05, FDR < 20%
merged.mar.igap.mmeta.pvalues <- filter_at(merged.mar.igap.mmeta, vars(ends_with(".fdr"), -matches("MVP"), -matches("CARDIOGENICS"),-matches("STARNET"),-matches("MESA")), any_vars(. <= 0.20)) %>%
  dplyr::select(ends_with(".p"),ends_with(".heidi"), ends_with(".gwas"), -matches("CARDIOGENICS"),-matches("STARNET"),-matches("MESA"), -matches("MVP"), starts_with("Gene"))

merged.mar.igap.mmeta.pvalues <- filter_at(merged.mar.igap.mmeta.pvalues, vars(ends_with(".gwas"), -matches("MVP"), -matches("CARDIOGENICS"),-matches("STARNET"),-matches("MESA")), any_vars(. <= 0.00005)) %>%
  dplyr::select(ends_with(".p"),ends_with(".heidi"), -matches("CARDIOGENICS"),-matches("STARNET"),-matches("MESA"), -matches("MVP"), starts_with("Gene"))

merged.mar.igap.mmeta.pvalues <- filter_at(merged.mar.igap.mmeta.pvalues, vars(ends_with(".heidi")), any_vars(. > 0.05)) %>%
  dplyr::select(ends_with(".p"), -matches("CARDIOGENICS"),-matches("STARNET"),-matches("MESA"), -matches("MVP"), starts_with("Gene"))

merged.mar.igap.mmeta.pvalues <- filter_at(merged.mar.igap.mmeta.pvalues, vars(ends_with(".p")), any_vars(. >0)) %>%
  dplyr::select(ends_with(".p"),-matches("CARDIOGENICS"),-matches("STARNET"),-matches("MESA"), -matches("MVP"), starts_with("Gene"))

merged.mar.igap.mmeta.pvalues <- merge(merged.mar.igap.mmeta.pvalues, new, by="Gene", all.x=TRUE)

#### Select all SNPs heidi > 0.05; smr_fdr <20% ; GWAS_p < 5 x 10-5
merged.mar.igap.mmeta.pvalues <- filter_at(merged.mar.igap.mmeta, vars(ends_with(".fdr"), -matches("CARDIOGENICS"),-matches("STARNET"),-matches("MVP")), any_vars(. <= 0.20)) %>%
  dplyr::select(ends_with(".p"),ends_with(".heidi"), ends_with(".gwas"), ends_with("topSNP"), ends_with("topSNP_chr"), ends_with("topSNP_bp"), -matches("CARDIOGENICS"),-matches("STARNET"),-matches("MVP"), starts_with("Gene"))

merged.mar.igap.mmeta.pvalues <- filter_at(merged.mar.igap.mmeta.pvalues, vars(ends_with(".gwas"), -matches("CARDIOGENICS"),-matches("STARNET"),-matches("MVP")), any_vars(. <= 0.0005)) %>%
  dplyr::select(ends_with(".p"),ends_with(".heidi"), ends_with("topSNP"), ends_with("topSNP_chr"), ends_with("topSNP_bp"),-matches("CARDIOGENICS"),-matches("STARNET"),-matches("MVP"), starts_with("Gene"))

merged.mar.igap.mmeta.pvalues <- filter_at(merged.mar.igap.mmeta.pvalues, vars(ends_with(".heidi")), any_vars(. > 0.05)) %>%
  dplyr::select(ends_with(".p"), ends_with("topSNP"),ends_with("topSNP_chr"), ends_with("topSNP_bp"), -matches("CARDIOGENICS"),-matches("STARNET"),-matches("MVP"), starts_with("Gene"))

merged.mar.igap.mmeta.pvalues <- filter_at(merged.mar.igap.mmeta.pvalues, vars(ends_with(".p")), any_vars(. >0)) %>%
  dplyr::select(ends_with(".p"),ends_with("topSNP"), ends_with("topSNP_chr"), ends_with("topSNP_bp"),-matches("CARDIOGENICS"),-matches("STARNET"),-matches("MVP"), starts_with("Gene"))

merged.mar.igap.mmeta.pvalues <- merge(merged.mar.igap.mmeta.pvalues, new, by="Gene", all.x=TRUE)


AUD_fBrain_eQTL <- subset(merged.mar.igap.mmeta.pvalues, METAANALYSIS_MVP_AUD_PGC_COGA1.rsids.eqtl_fetal_brain_brien_et_al_query_format_hg19.p <= 0.05 ,
                  select=c(Gene))$Gene
AUD_fBrain_mQTL <- subset(merged.mar.igap.mmeta.pvalues, METAANALYSIS_MVP_AUD_PGC_COGA1.rsids.Brain_m_Meta.p <= 0.05 ,
                          select=c(Gene))$Gene
DPW_fBrain_eQTL <- subset(merged.mar.igap.mmeta.pvalues, Liu2019drnkwk.eqtl_fetal_brain_brien_et_al_query_format_hg19.p <= 0.05 ,
                          select=c(Gene))$Gene
DPW_fBrain_mQTL <- subset(merged.mar.igap.mmeta.pvalues, Liu2019drnkwk.FB_Brain_mQTL_rsid.p <= 0.05 ,
                          select=c(Gene))$Gene
group.venn(list(AUD_fBrain_eQTL=AUD_fBrain_eQTL, AUD_fBrain_mQTL=AUD_fBrain_mQTL, DPW_fBrain_eQTL = DPW_fBrain_eQTL, DPW_fBrain_mQTL = DPW_fBrain_mQTL), label=F, 
           fill = c("orange", "blue", "red", "green"),
           cat.pos = c(0, 0,0,0),
           lab.cex=0.2)
venn(list(AUD_fBrain_eQTL=AUD_fBrain_eQTL, AUD_fBrain_mQTL=AUD_fBrain_mQTL), fill = c("orange", "blue"))

##Emma's Manhattan plot

p1 <- ggman(df1, snp = "gene", bp = "bp", chrom = "chrom", pvalue = "pvalue", sigLine = fdr.cut, relative.positions = TRUE, title = "Schz + Alcohol with brain eQTLS") + theme_bw()
ggmanLabel((p1+scale_color_manual(values=c("coral4","coral3"))) , labelDfm = df1.sig, snp = "gene", label = "gene")
           

newnames = c("dpw_CARD_MC", " dpw_CARD_MP", " dpw_MESA_MC", " dpw_DLFPC_eMETA", "  dpw_DLFPC_mMETA", " dpw_STARNET_FC", " dpw_STARNET_MP", " dpw_FB-eQTL", " dpw_FB-mQTL", " AUD_DLFPC_mMETA", " AUD_FB-eQTL", " AUD_FB-mQTL", " AUD_CARD_MC", " AUD_CARD_MP", " AUD_MESA_MC", " AUD_DLFPC_eMETA", " AUD_STARNET_FC", " AUD_STARNET_MP")





