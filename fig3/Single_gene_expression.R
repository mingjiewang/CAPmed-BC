#!/usr/bin/env Rscript
#### Input args ####
library(jsonlite)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No JSON input provided.")
}
json_input <- args[1]
# json_input <- "Single_gene_expression.json"
print(json_input)
conf <- fromJSON(json_input)

#### pkgs and data ####
library(ggpubr)
library(ggsci)
library(cowplot)
library(openxlsx)
# load data
ln = load("mrna.Rdata")
ln

# read in conf
groupCol = conf$ColPalette
gene.want = conf$Gene_Name
group1 = conf$Group1
group2 = conf$Group2
showAll = conf$ShowAll

#### analysis ####
gene.exp = data.frame(
  Sample = colnames(datExp),
  Expression = unlist(datExp[rownames(datExp) == gene.want,]),
  Group = as.factor(target$Group1))

gene.exp$Group = factor(gene.exp$Group,
                        levels=c("nonTNBC_Control",
                                 "nonTNBC_CAP_1h",
                                 "nonTNBC_CAP_8h",
                                 "TNBC_Control",
                                 "TNBC_CAP_1h",
                                 "TNBC_CAP_8h"))
if(showAll){
  gene.exp.want = gene.exp
}else{
  gene.exp.want = gene.exp[gene.exp$Group %in% c(group1, group2),]
}


##### plot #####
p1 <- ggviolin(gene.exp.want,x = "Group",y="Expression",color = NA,
         fill="Group",palette = groupCol,alpha = 0.3,
         add = "boxplot",add.params = list(color="black",size=0.05, alpha=0.7),
         title = paste(gene.want,"expression")) +
  theme_bw() +
  theme(title=element_text(size=15),
        legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  stat_compare_means(comparisons = list(c(group1, group2)),method = "t.test",label.y = max(gene.exp$Expression) * 1.2)
save_plot("GeneExpression.violinplot.pdf",p1,base_width = ifelse(showAll,7,5), base_height = 6)
write.xlsx(gene.exp,"GeneExpression.xlsx")
