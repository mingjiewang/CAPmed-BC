#!/usr/bin/env Rscript
#### Input args ####
library(jsonlite)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No JSON input provided.")
}
json_input <- args[1]
#json_input <- "GSEA.json"
print(json_input)
conf <- fromJSON(json_input)

#### pkgs and data ####
library(cowplot)
library(openxlsx)
# load data
ln = load("lncrna.Rdata")
ln

# read in conf
group1 = conf$Group1
group2 = conf$Group2
pval.cutoff = conf$DEG_p_cutoff
fc.cutoff = conf$DEG_FC_cutoff
term.number = conf$GSEA_term_number


##### DEG #####
library(limma)
#model design
Group = factor(target$Group1,levels=c(group1,group2))
design = model.matrix(~0+Group)
colnames(design) <- c("group1","group2")
design

#linear model fitness
data = datExp[,target$Group1 %in% c(group1, group2)]

fit <- lmFit(data, design)

#generate contrast matrix
contrast.matrix <- makeContrasts(group1 - group2, #1
                                 levels=design)
#constrast model fit 
fit2 <- contrasts.fit(fit, contrast.matrix)

#bayes model 
fit2 <- eBayes(fit2)

#get DEGs
p.adjust.cutoff = pval.cutoff
fold.change.cutoff = fc.cutoff
diff = topTable(fit2,adjust.method="fdr",
                p.value=p.adjust.cutoff,
                lfc=log(fold.change.cutoff,2),
                number=50000,sort.by = 'logFC')
diff$lncRNA_ID = rownames(diff)
diff$target_gene = anno$target_gene[match(diff$lncRNA_ID,anno$LncRNA_id)]
diff$Compare = paste(group1, "-", group2)

# filteration
library(stringr)
diff = diff[,c(7,8,1:6,9)]

##### GSEA analysis #####
# packages
pacman::p_load("enrichplot")
pacman::p_load("org.Hs.eg.db")
pacman::p_load("clusterProfiler")

# Load database
database <- org.Hs.egSYMBOL2EG
database <- as.list(database)

# match gene
gene.id.enrich <- database[match(unique(na.omit(diff$target_gene)), names(database))]
gene.id.enrich <- data.frame(GeneName = names(gene.id.enrich), ID = as.character(gene.id.enrich))
deg.file.enroll <- merge(gene.id.enrich, diff, by.x = "GeneName", by.y = "target_gene")
foldChange <- as.numeric(deg.file.enroll$logFC)
names(foldChange) <- deg.file.enroll$ID

foldChange <- foldChange[order(foldChange, decreasing = T)]
edo.gsea <- gseGO(foldChange, ont = "BP", OrgDb = org.Hs.eg.db,
                  keyType = "ENTREZID", pvalueCutoff = 1, pAdjustMethod = "BH")

# gene enriched gene name
id2gene <- function(id.list){
  tmp = bitr(strsplit(id.list,"/")[[1]],fromType = "ENTREZID",toType = "SYMBOL",OrgDb = "org.Hs.eg.db")
  return(paste(tmp$SYMBOL,collapse = "/"))
}

## output table
edo.gsea.out = as.data.frame(edo.gsea)[1:50,]
tmp = unlist(as.data.frame(sapply(edo.gsea.out$core_enrichment, id2gene))[1])
colnames(tmp) = NULL
edo.gsea.out$geneName = tmp
write.xlsx(edo.gsea.out,"GSEA.xlsx")

# Combine GSEA plot
view.id.list <- 1:term.number
p.combine = gseaplot2(edo.gsea, geneSetID = view.id.list,
                      pvalue_table = F)
p.combine


## output plots together
save_plot("GSEA.pdf",p.combine,base_width = 8, base_height = 8)

