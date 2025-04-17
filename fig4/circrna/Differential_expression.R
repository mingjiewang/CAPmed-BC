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
ln = load("circrna.Rdata")
ln

# read in conf
group1 = conf$Group1
group2 = conf$Group2
pval.cutoff = conf$pCutoff
fc.cutoff = conf$fcCutoff


#### analysis ####
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
library(openxlsx)
write.xlsx(diff,"DEG.xlsx")

##### volcanoplot #####
deg.data = topTable(fit2,adjust.method="fdr",p.value=1,
                    lfc=log(1,2),number=50000,sort.by = 'logFC')
deg.data$Symbol = rownames(deg.data)
deg.data$logP <- -log10(deg.data$P.Value)

#p.adjust.cutoff = 0.01
p.cutoff = pval.cutoff
fold.change.cutoff = fc.cutoff

deg.data$Group = "Not significant"
deg.data$Group[which( (deg.data$P.Value < p.cutoff) & (deg.data$logFC > log(fold.change.cutoff,2)) )] = "Up-regulated"
deg.data$Group[which( (deg.data$P.Value < p.cutoff) & (deg.data$logFC < -log(fold.change.cutoff,2)) )] = "Down-regulated"

deg.data$Label = ""
deg.data <- deg.data[order(deg.data$P.Value), ]
up.genes <- head(deg.data$Symbol[which(deg.data$Group == "Up-regulated")], 10)
down.genes <- head(deg.data$Symbol[which(deg.data$Group == "Down-regulated")], 10)
deg.top10.genes <- c(as.character(up.genes), as.character(down.genes))
deg.data$Label[match(deg.top10.genes, deg.data$Symbol)] <- deg.top10.genes


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
save_plot("DEG.volcano.pdf",p1,base_width = 8, base_height = 6)

##### heatmap #####
## draw heatmp
library(pheatmap)
annotation_col = as.data.frame(target$Group1[target$Group1 %in% c(group1, group2)])
colnames(annotation_col) = "Group"
deg.exp =  data[match(diff$circRNA_ID,rownames(data)),]
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

save_plot("DEG.heatmap.pdf",p2,base_width = 6, base_height = 6)
