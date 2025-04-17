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
ln = load("protome.Rdata")
ln

# read in conf
gene.want = conf$query_ID
groupCol = conf$ColPalette

#### analysis ####
gene.exp = data.frame(
  Sample = colnames(datExp)[1:6],
  Expression = unlist(datExp[rownames(datExp) == gene.want, 1:6]),
  Group = as.factor(target$Group)
)


gene.exp$Group = factor(gene.exp$Group,
                        levels=c("TNBC_Control",
                                 "TNBC_CAP"))

##### plot #####
p1 <- ggviolin(gene.exp,x = "Group",y="Expression",color = NA,
               fill="Group",palette = groupCol,alpha = 0.3,
               add = "boxplot",add.params = list(color="black",size=0.05, alpha=0.7),
               title = paste(gene.want,"expression")) +
  theme_bw() +
  theme(title=element_text(size=15),
        legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  annotate("segment",x=1,y=max(gene.exp$Expression)*1.2,xend=2,yend=max(gene.exp$Expression)*1.2)+
  annotate("text",x=1.5,y=max(gene.exp$Expression)*1.3,label=sprintf("%0.4f",datExp$P.value[rownames(datExp) == gene.want]))
p1
save_plot("GeneExpression.violinplot.pdf",p1,base_width = 5, base_height = 6)
write.xlsx(gene.exp,"GeneExpression.xlsx")
