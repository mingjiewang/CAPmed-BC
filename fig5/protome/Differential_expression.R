#!/usr/bin/env Rscript
#### Input args ####
library(jsonlite)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No JSON input provided.")
}
json_input <- args[1]
# json_input <- "Differential_expression.json"
print(json_input)
conf <- fromJSON(json_input)

#### pkgs and data ####
library(ggpubr)
# load data
ln = load("protome.Rdata")
ln

# read in conf
pval.cutoff = conf$pCutoff
fc.cutoff = conf$fcCutoff


#### analysis ####
fc.idx = (datExp$Fold.change > fc.cutoff) | (datExp$Fold.change < (1/fc.cutoff))
pval.idx = datExp$P.value < pval.cutoff

library(dplyr)
diff = datExp[fc.idx & pval.idx,] %>%
  na.omit() %>%
  arrange(desc(abs(log(Fold.change,2))))
diff$Protein.name = rownames(diff)
diff$Protein.accession = anno$Protein.accession[match(diff$Protein.name,anno$Gene.name)]
diff$Protein.description = anno$Protein.description[match(diff$Protein.name,anno$Gene.name)]
diff$Compare = "TNBC_CAP - TNBC_Control"
diff$logFC = log(diff$Fold.change,2)
diff.out = diff[,c(9,10,13,8,11,12)]
#View(head(diff.out))
library(openxlsx)
write.xlsx(diff.out,file = "DEP.xlsx")

##### volcanoplot #####
deg.data = na.omit(datExp)
deg.data$Symbol = rownames(deg.data)
deg.data$logP <- -log10(deg.data$P.value)
deg.data$logFC <- log(deg.data$Fold.change,2)

#p.adjust.cutoff = 0.01
p.cutoff = pval.cutoff
fold.change.cutoff = fc.cutoff

deg.data$Group = "Not significant"
deg.data$Group[(deg.data$P.value < p.cutoff) & (deg.data$logFC > log(fold.change.cutoff,2))] = "Up-regulated"
deg.data$Group[(deg.data$P.value < p.cutoff) & (deg.data$logFC < -log(fold.change.cutoff,2))] = "Down-regulated"
table(deg.data$Group)

p1 <- ggscatter(deg.data, x = "logFC", y = "logP",
                color = "Group", 
                palette = c("#2f5688", "#BBBBBB", "#CC0000"),
                size = 1,
                # label = deg.data$Label, 
                font.label = 8, 
                show.legend.text = F,
                repel = T,
                xlab = "log2(Fold Change)", 
                #ylab = "-log10(Adjust P-value)"
                ylab = "-log10(P-value)") + 
  theme_bw() + 
  geom_hline(yintercept = -log(p.cutoff,10), linetype="dashed") +
  #geom_hline(yintercept = -log(p.adjust.cutoff,10), linetype="dashed") +
  geom_vline(xintercept = c(log(fold.change.cutoff,2),-log(fold.change.cutoff,2)), linetype="dashed")


## output plot 
library(cowplot)
save_plot("DEP.volcano.pdf",p1,base_width = 8, base_height = 6)

##### heatmap #####
## draw heatmp
library(pheatmap)
annotation_col = as.data.frame(target$Group)
colnames(annotation_col) = "Group"
deg.exp =  diff[,c(1:6)]
rownames(annotation_col) = colnames(deg.exp)
color.key <- c("#3300CC", "#3399FF", "white", "#FF3333", "#CC0000")
sample.color.list = list(Group=c(get_palette("nejm",length(unique(annotation_col[,1])))))
names(sample.color.list$Group) = unique(annotation_col[,1])

# get expression values of DEGs
p2 = pheatmap(deg.exp, color = colorRampPalette(color.key)(50), 
              border_color = NA,
              annotation_col = annotation_col,
              labels_row = NULL,clustering_method = "ward.D2",
              show_rownames = F,show_colnames = T,fontsize_col = 5,
              annotation_colors = sample.color.list
)

save_plot("DEP.heatmap.pdf",p2,base_width = 6, base_height = 6)
