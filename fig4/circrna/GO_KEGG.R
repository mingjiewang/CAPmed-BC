#!/usr/bin/env Rscript
#### Input args ####
library(jsonlite)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No JSON input provided.")
}
json_input <- args[1]
# json_input <- "GO_KEGG.json"
print(json_input)
conf <- fromJSON(json_input)

#### pkgs and data ####
library(cowplot)
library(openxlsx)
# load data
ln = load("circrna.Rdata")
ln

# read in conf
group1 = conf$Group1
group2 = conf$Group2
pval.cutoff = conf$DEG_p_cutoff
fc.cutoff = conf$DEG_FC_cutoff
term.number = conf$GO_term_number


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
diff$circRNA_ID = rownames(diff)
diff$GeneName = anno$Gene.Symbol[match(diff$circRNA_ID,anno$Circ_id)]
diff$Compare = paste(group1, "-", group2)
diff = na.omit(diff)

# filteration
library(stringr)
diff = diff[,c(7,8,1:6,9)]
diff = diff[diff$GeneName != "---",]
diff = diff[!str_detect(diff$GeneName,"///"),]

##### GO analysis #####
# packages
pacman::p_load("enrichplot")
pacman::p_load("org.Hs.eg.db")
pacman::p_load("clusterProfiler")

# Load database
database <- org.Hs.egSYMBOL2EG
database <- as.list(database)


gene.list <- diff$GeneName %>%
  stringr::str_remove("^ +") %>%
  stringr::str_remove(" $+")

gene.list.enrich <- database[names(database) %in% gene.list]

# GO and KEGG enrichment
edo.bp <- enrichGO(as.character(unlist(gene.list.enrich)), OrgDb = org.Hs.eg.db,
                   ont = 'BP', pAdjustMethod = 'BH', 
                   pvalueCutoff = 1, qvalueCutoff = 1,
                   keyType = 'ENTREZID',readable = T)

edo.mf <- enrichGO(as.character(unlist(gene.list.enrich)), OrgDb = org.Hs.eg.db,
                   ont = 'MF', pAdjustMethod = 'BH', 
                   pvalueCutoff = 1, qvalueCutoff = 1,
                   keyType = 'ENTREZID',readable = T)

edo.cc <- enrichGO(as.character(unlist(gene.list.enrich)), OrgDb = org.Hs.eg.db,
                   ont = 'CC', pAdjustMethod = 'BH', 
                   pvalueCutoff = 1, qvalueCutoff = 1,
                   keyType = 'ENTREZID',readable = T)

edo.kegg <- enrichKEGG(as.character(unlist(gene.list.enrich)), organism = "hsa",
                       pAdjustMethod = 'BH', keyType="kegg",
                       pvalueCutoff = 1, qvalueCutoff = 1)

# gene enriched gene name
id2gene <- function(id.list){
  tmp = bitr(strsplit(id.list,"/")[[1]],fromType = "ENTREZID",toType = "SYMBOL",OrgDb = "org.Hs.eg.db")
  return(paste(tmp$SYMBOL,collapse = "/"))
}

## to excel


## bp out
edo.bp.out = as.data.frame(edo.bp)[1:50,]
edo.bp.out$Functional_enrichment_type = "Biological process"

## cc out
edo.cc.out = as.data.frame(edo.cc)[1:50,]
edo.cc.out$Functional_enrichment_type = "Cellular components"

## MF out
edo.mf.out = as.data.frame(edo.mf)[1:50,]
edo.mf.out$Functional_enrichment_type = "Molecular function"

## KEGG out
edo.kegg.out = as.data.frame(edo.kegg)[1:50,]
tmp = unlist(as.data.frame(sapply(edo.kegg.out$geneID, id2gene))[1])
colnames(tmp) = NULL
edo.kegg.out$geneName = tmp
edo.kegg.out = edo.kegg.out[,c(1:7,10,9)]
colnames(edo.kegg.out)[8] = "geneID"
edo.kegg.out$Functional_enrichment_type = "KEGG"

edo.out <- rbind(
  edo.bp.out,
  edo.cc.out,
  edo.mf.out,
  edo.kegg.out
)

edo.out <- edo.out[,c(10,1:9)]
write.xlsx(edo.out,"GO_KEGG.xlsx")

##Bubble plot 
# Packages
library(Hmisc)
library(ggplot2)
library(stringr)
library(cowplot)

# Functions to draw plots
DrawGOBubblePlot <- function(dat, category = "BP", top.number = 15, col="blue"){
  # Draw bubble plot using DAVID function enrichment results
  
  category = toupper(category)
  if (category == "BP"){
    main.title = "Biological Process"
  } else if (category == "CC"){
    main.title = "Cellular Components"
  } else if (category == "MF"){
    main.title = "Molecular Function"
  } else if (category == "KEGG"){
    main.title = "KEGG"
  } else {
    return("ERROR! Wrong input parameter [category].")
  }
  
  dat1 = dat[c(1:top.number),c(2,3,4,5,9)]
  dat1[,2] = str_remove(dat1[,2],"/.*")
  dat1[,3] = str_remove(dat1[,3],"/.*")
  dat1$Ratio = as.numeric(dat1$GeneRatio) / as.numeric(dat1$BgRatio)
  
  
  dat1$Description = capitalize(dat1$Description)
  dat1$Description = factor(dat1$Description,levels=dat1$Description[length(dat1$Description):1])
  dat1$pvalue = -log10(dat1$pvalue)
  
  p = ggplot(dat1,aes(Ratio,Description)) +
    geom_point(aes(size=Count,colour=pvalue)) +
    scale_colour_gradient(low=col,high="red") + 
    labs(colour=expression(-log[10]("P Value")),size="Gene counts",  
         x="Gene Ratio",y="",title=main.title) +
    theme_bw() +
    scale_x_continuous(limits = c(0,max(dat1$Ratio) * 1.2)) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 50))
  
  return(p)
}

# Read in data and generate the plots
#BP
p1 = DrawGOBubblePlot(edo.bp.out,"BP",10,"blue")

#MF
p2 = DrawGOBubblePlot(edo.mf.out,"MF",10,"blue")

#CC
p3 = DrawGOBubblePlot(edo.cc.out,"CC",10,"blue")

#KEGG
p4 = DrawGOBubblePlot(edo.kegg.out,"KEGG",10,"blue")

## output plots together
plot.out = plot_grid(p1,p2,p3,p4,labels=c("A","B","C","D"),align = "h",nrow=2)
save_plot("GO_KEGG.pdf",plot.out,base_width = 13, base_height = 9)

