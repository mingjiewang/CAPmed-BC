#!/usr/bin/env Rscript
#### Input args ####
library(jsonlite)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No JSON input provided.")
}
json_input <- args[1]
# json_input <- "Protein_modification_quantification.json"
print(json_input)
conf <- fromJSON(json_input)

#### pkgs and data ####
library(ggpubr)
library(ggsci)
library(cowplot)
library(openxlsx)
# load data
ln = load("lac.Rdata")
ln

# read in conf
site.want = conf$query_site
groupCol = conf$ColPalette

#### analysis ####
library(stringr)
gene.exp = data.frame(
  Sample = colnames(datExp)[1:9],
  Expression = unlist(datExp[rownames(datExp) == site.want, 1:9]),
  Group = as.factor(target$Group)
)


gene.exp$Group = factor(gene.exp$Group,
                        levels=c("TNBC_Control",
                                 "TNBC_CAP",
                                 "nonTNBC_CAP"))

##### plot #####
p1 <- ggviolin(gene.exp,x = "Group",y="Expression",color = NA,
               fill="Group",palette = groupCol,alpha = 0.3,
               add = "boxplot",add.params = list(color="black",size=0.05, alpha=0.7),
               title = paste(site.want,"quantification")) +
  theme_bw() +
  theme(title=element_text(size=15),
        legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  annotate("segment",x=1,y=max(gene.exp$Expression)*1.2,xend=2,yend=max(gene.exp$Expression)*1.2)+
  annotate("text",x=1.5,y=max(gene.exp$Expression)*1.3,label=sprintf("%0.4f",datExp$`P.value.TNBC_CAP-TNBC_Control`[rownames(datExp) == site.want]))+
  annotate("segment",x=2,y=max(gene.exp$Expression)*1.4,xend=3,yend=max(gene.exp$Expression)*1.4)+
  annotate("text",x=2.5,y=max(gene.exp$Expression)*1.5,label=sprintf("%0.4f",datExp$`P.value.nonTNBC_CAP-TNBC_CAP`[rownames(datExp) == site.want]))

p1
save_plot("Protein.site.modification.violinplot.pdf",p1,base_width = 6, base_height = 6)
write.xlsx(gene.exp,"Protein.site.modification.xlsx")
