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
ln = load("protome.Rdata")
ln

# read in conf
pval.cutoff = conf$DEG_p_cutoff
fc.cutoff = conf$DEG_FC_cutoff
term.number = conf$GSEA_term_number


##### DEG #####
fc.idx = (datExp$Fold.change > fc.cutoff) | (datExp$Fold.change < (1/fc.cutoff))
pval.idx = datExp$P.value < pval.cutoff

library(dplyr)
diff = datExp[fc.idx & pval.idx,] %>%
  na.omit() %>%
  arrange(desc(abs(log(Fold.change,2))))
diff$Protein.name = rownames(diff)


##### GSEA analysis #####
# packages
pacman::p_load("enrichplot")
pacman::p_load("org.Hs.eg.db")
pacman::p_load("clusterProfiler")

# Load database
database <- org.Hs.egSYMBOL2EG
database <- as.list(database)

# match gene
gene.id.enrich <- database[match(diff$Protein.name, names(database))]
gene.id.enrich <- data.frame(GeneName = names(gene.id.enrich), ID = as.character(gene.id.enrich))
gene.id.enrich <- na.omit(gene.id.enrich)
deg.file.enroll <- merge(gene.id.enrich, diff, by.x = "GeneName", by.y = "Protein.name")
foldChange <- as.numeric(log(deg.file.enroll$Fold.change,2))
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

