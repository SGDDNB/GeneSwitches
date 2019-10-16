## code to prepare `DATASET` dataset goes here

###pathways ----------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("fgsea")
# BiocManager::install("fgsea", version = "3.8")
library(fgsea)
pathways1 <- gmtPathways("~/Documents/geneswitch/pathways/msigdb/h.all.v6.2.symbols.gmt")
pathways2 <- gmtPathways("~/Documents/geneswitch/pathways/msigdb/c2.cp.kegg.v6.2.symbols.gmt")
pathways3 <- gmtPathways("~/Documents/geneswitch/pathways/msigdb/c5.all.v6.2.symbols.gmt")
msigdb_h_c2_c5 <- c(pathways1, pathways2, pathways3);length(pathways)
# save(pathways, file = "../shiny/pathways.RData")
usethis::use_data(msigdb_h_c2_c5, internal = FALSE)

###genelists ----------------------------------
load("~/Documents/general_scripts/genelists.RData")
table(genelists$genetype)
genelists <- genelists[genelists$genetype %in% c("TF", "surface"),];dim(genelists)
genelists[genelists$genetype == "surface", ]$genetype <- "Surface protein"
table(genelists$genetype)
# usethis::use_data(genelists, internal = TRUE, overwrite = TRUE)
genelists <- genelists[genelists$genetype == "TF",];dim(genelists)
genelists$genetype <- "TFs"
head(genelists)

gl <- read.table("~/Documents/general_scripts/CellSurfaceAtlas.txt", header = TRUE)
head(gl);class(gl);dim(gl)
# length(setdiff(genelists[genelists$genetype == "Surface protein",]$genenames, gl$entrez_symbol))
# length(setdiff(gl$entrez_symbol, genelists[genelists$genetype == "Surface protein",]$genenames))
gl2 <- data.frame(genenames = unique(gl$entrez_symbol), genetype = "Surface proteins")
rownames(gl2) <- gl2$genenames;dim(gl2)
head(gl2)
com <- intersect(rownames(gl2), rownames(genelists))
genelists[com,]$genetype <- "TF&Surface proteins"
gl2 <- gl2[!rownames(gl2) %in% genelists$genenames,];dim(gl2)
genelists <- rbind(gl2, genelists)
head(genelists);table(genelists$genetype)
gs_genelists <- genelists
colnames(gs_genelists)[2] <- "genetypes"
head(gs_genelists);table(gs_genelists$genetype)
usethis::use_data(gs_genelists, internal = FALSE, overwrite = TRUE)

genelists3 <- genelists
save(genelists3, file = "~/Documents/general_scripts/genelists3(new).RData")

gl <- read.table("~/Documents/general_scripts/ECM_proteins.txt", sep = "\t", header = TRUE)
head(gl);class(gl);dim(gl)
# length(setdiff(genelists[genelists$genetype == "Surface protein",]$genenames, gl$entrez_symbol))
# length(setdiff(gl$entrez_symbol, genelists[genelists$genetype == "Surface protein",]$genenames))
gl2 <- data.frame(genenames = unique(gl$entrez_symbol), genetype = "ECM")
rownames(gl2) <- gl2$genenames;dim(gl2)
head(gl2)
com <- intersect(rownames(gl2), rownames(genelists))
genelists[com,]
genelists <- genelists[!genelists$genenames %in% rownames(gl2),];dim(genelists)
genelists <- rbind(gl2, genelists)
head(genelists);table(genelists$genetype)

gl <- read.table("~/Documents/general_scripts/CytoKineRegistry.txt", sep = "\t", header = TRUE)
head(gl);class(gl);dim(gl)
# length(setdiff(genelists[genelists$genetype == "Surface protein",]$genenames, gl$entrez_symbol))
# length(setdiff(gl$entrez_symbol, genelists[genelists$genetype == "Surface protein",]$genenames))
gl2 <- data.frame(genenames = unique(gl$EntrezGeneSymbol), genetype = "Cytokine")
rownames(gl2) <- gl2$genenames;dim(gl2)
head(gl2)
com <- intersect(rownames(gl2), rownames(genelists));length(com)
genelists[com,]
genelists <- genelists[!genelists$genenames %in% rownames(gl2),];dim(genelists)
genelists <- rbind(gl2, genelists)
head(genelists);table(genelists$genetype)
genelists2 <- genelists
save(genelists2, file = "~/Documents/general_scripts/genelists2(new).RData")

###time monocle ----------------------------------
library(monocle)

seu3obj <- readRDS("~/Documents/geneswitch/cardiomyocytes/cardiacSeurat.rds")
allexpdata <- as.matrix(GetAssayData(object = seu3obj, slot = "counts"));dim(allexpdata)
cellinfo <- seu3obj@meta.data;dim(cellinfo)
all(rownames(cellinfo) == colnames(allexpdata))

#####subset cells
cellinfo_sub <- cellinfo[sample(nrow(cellinfo), 3000),];dim(cellinfo_sub)
table(cellinfo_sub$library)
head(cellinfo_sub)
expdata <- allexpdata[,rownames(cellinfo_sub)];dim(expdata)
all(rownames(cellinfo_sub) == colnames(expdata))

genenames <- as.data.frame(rownames(expdata));head(genenames)
colnames(genenames)[1] <- "gene_short_name"
rownames(genenames) <- genenames$gene_short_name;head(genenames)
cells <- cellinfo_sub
all(rownames(cells) == colnames(expdata))
all(rownames(genenames) == rownames(expdata))
###################

pd <- new("AnnotatedDataFrame",data=cells)
fd <- new("AnnotatedDataFrame",data=genenames)
# First create a CellDataSet from the relative expression levels
ORMM <- newCellDataSet(as.matrix(expdata), phenoData = pd, featureData =fd,
                       lowerDetectionLimit = 1, expressionFamily = negbinomial.size())

ORMM <- estimateSizeFactors(ORMM)
ORMM <- estimateDispersions(ORMM)
##Filtering low-quality cells
ORMM<- detectGenes(ORMM, min_expr = 0.1)
print(head(fData(ORMM)))
expressed_genes <- row.names(subset(fData(ORMM), num_cells_expressed >=10))
length(expressed_genes)
print(head(pData(ORMM)));dim(pData(ORMM))

##good to look at the distribution of mRNA totals across the cells
pData(ORMM)$Total_mRNAs <- Matrix::colSums(exprs(ORMM))

ORMM <- ORMM[,pData(ORMM)$Total_mRNAs < 1e6]

upper_bound <- 10^(mean(log10(pData(ORMM)$Total_mRNAs)) +
                     2*sd(log10(pData(ORMM)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(ORMM)$Total_mRNAs)) -
                     2*sd(log10(pData(ORMM)$Total_mRNAs)))
pdf(paste(pid,"_1total_mRNAs_distribution.pdf",sep = ""))
qplot(Total_mRNAs, data = pData(ORMM), color = Characteristics.inferred.lineage., geom = "density") +
  geom_vline(xintercept = lower_bound) +
  geom_vline(xintercept = upper_bound)
dev.off()
##removed the few cells with either very low mRNA recovery or far more mRNA that the typical cell
ORMM <- ORMM[,pData(ORMM)$Total_mRNAs > lower_bound &
               pData(ORMM)$Total_mRNAs < upper_bound]
ORMM <- detectGenes(ORMM, min_expr = 0.1)
# Log-transform each value in the expression matrix.
L <- log(exprs(ORMM[expressed_genes,]))
#install.packages("reshape2")
library(reshape2)
# Standardize each gene, so that they are all on the same scale,
# Then melt the data with plyr so we can plot it easily
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))

# Plot the distribution of the standardized gene expression values.
pdf(paste(pid,"_2standardized_gene_exp.pdf",sep = ""))
qplot(value, geom = "density", data = melted_dens_df) +
  stat_function(fun = dnorm, size = 0.5, color = 'red') +
  xlab("Standardized log(FPKM)") +
  ylab("Density")
dev.off()

#####Ordering based on genes that differ between clusters
disp_table <- dispersionTable(ORMM)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
ORMM <- setOrderingFilter(ORMM, unsup_clustering_genes$gene_id)
pdf(paste(pid,"_3ordering_genes.pdf",sep = ""))
plot_ordering_genes(ORMM)
dev.off()
pdf(paste(pid,"_4pc_variance.pdf",sep = ""))
plot_pc_variance_explained(ORMM, return_all = F) # norm_method='log'
dev.off()
ORMM <- reduceDimension(ORMM, max_components = 2, norm_method = 'log', num_dim = 5,
                        residualModelFormulaStr = "~num_genes_expressed",
                        reduction_method = 'tSNE', verbose = T)
ORMM <- clusterCells(ORMM, verbose = T)

##Constructing Single Cell Trajectories
diff_test_res <- differentialGeneTest(ORMM[expressed_genes,],
                                      fullModelFormulaStr = "~library",
                                      cores = 5, verbose = TRUE)
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
ORMM <- setOrderingFilter(ORMM, ordering_genes)
pdf(paste(pid,"_5ordering_genes_fortree.pdf",sep = ""))
plot_ordering_genes(ORMM)
dev.off()
ORMM <- reduceDimension(ORMM, max_components=2, method = 'DDRTree')
ORMM <- orderCells(ORMM)#, num_paths=2)
# save(ORMM, file = "/home/yiqun/Documents/geneswitch/cardiomyocytes/monocle2/ORMM(monocle_DDRTree).RData")
expressed_genes <- row.names(subset(fData(ORMM), num_cells_expressed >=10))
ORMM_filtered <- ORMM[expressed_genes,]
my_genes <- row.names(subset(fData(ORMM_filtered),
                             gene_short_name %in% c("POU5F1","SOX17","EOMES","ISL1","TNNI1","MYL7","THY1")))
cds_subset <- ORMM_filtered[my_genes,]
plot_cell_trajectory(ORMM, color_by = "Pseudotime", markers = "SOX17")


cardiac_monocle2 <- ORMM
cardiac_monocle2@reducedDimS[,1:5]
cardiac_monocle2@reducedDimS[1,] <- -(cardiac_monocle2@reducedDimS[1,])
cardiac_monocle2@reducedDimK[1,] <- -(cardiac_monocle2@reducedDimK[1,])
head(pData(cardiac_monocle2))
# colnames(pData(cardiac_monocle2))[7] <- "Clusters"
tiff("monocle_pseudotime.tiff", units="in", width=5, height=5, res=300)
plot_cell_trajectory(cardiac_monocle2, color_by = "Pseudotime", cell_size = 1)
dev.off()
plot_cell_trajectory(cardiac_monocle2, color_by = "library", cell_size = 1)

pdf(paste("monocle_","trajectory_cardiac",".pdf",sep = ""))
plot_cell_trajectory(cardiac_monocle2, color_by = "library")
plot_cell_trajectory(cardiac_monocle2, color_by = "RNA_snn_res.0.1")
plot_cell_trajectory(cardiac_monocle2, color_by = "Pseudotime")
# plot_genes_in_pseudotime(cds_subset, color_by = "library")
dev.off()
save(cardiac_monocle2, file = "cardiac_monocle2.RData")

mcells <- pData(cardiac_monocle2);head(mcells)
head(mcells);dim(mcells)

seu3obj <- readRDS("~/Documents/geneswitch/cardiomyocytes/cardiacSeurat.rds")
logexpdata <- as.matrix(GetAssayData(object = seu3obj, slot = "data")[,rownames(mcells)]);dim(logexpdata)
all(rownames(mcells) == colnames(logexpdata))
logexpdata[1:5,1:5]
save(logexpdata, file = "logexpdata.RData")

# load("~/Documents/geneswitch/cardiomyocytes/monocle2/cardiac_monocle2.RData")
# load("~/Documents/geneswitch/cardiomyocytes/monocle2/logexpdata.RData")
# usethis::use_data(logexpdata, cardiac_monocle2, genelists, internal = TRUE, overwrite = TRUE)
usethis::use_data(genelists, internal = TRUE, overwrite = TRUE)


# You can install the released version of GeneSwitches from [CRAN](https://CRAN.R-project.org) with:

#``` r
#install.packages("GeneSwitches")
#```
#You can install the development version from [GitHub](https://github.com/) with:
# user  system elapsed
# 205.137   5.573 210.613
# user  system elapsed
# 141.616   2.425 143.954
# user   system  elapsed
# 1956.396  242.979 2197.818
