#!/usr/bin/env Rscript
#### Input args ####
library(jsonlite)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No JSON input provided.")
}
json_input <- args[1]
# json_input <- "PPI.json"
print(json_input)
conf <- fromJSON(json_input)

#### pkgs and data ####
library(cowplot)
library(openxlsx)
# load data
ln = load("pho.Rdata")
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
library(stringr)
diff$Protein.name = str_remove(rownames(diff),"_.*")


##### PPI #####
library(org.Hs.eg.db)
library(igraph)
library(ggraph)
library(clusterProfiler)
library(STRINGdb)

string_db <- STRINGdb$new(version="11", species = 9606, score_threshold = 400, input_directory="")

# covert to id
gene = unique(diff$Protein.name)
gene <- gene %>% bitr(
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = "org.Hs.eg.db",
  drop = T
)
## only top 200 selected
if(nrow(gene) > 200){
  gene = gene[1:200,]
}

data_mapped <- gene %>% string_db$map(my_data_frame_id_col_names = "ENTREZID",
                                      removeUnmappedRows = T)
# string_db$plot_network(data_mapped$STRING_id)

## only top 200 selected
data_links <- data_mapped$STRING_id %>% string_db$get_interactions()

# 转换stringID为Symbol，只取前两列和最后一列
links <- data_links %>%
  mutate(from = data_mapped[match(from, data_mapped$STRING_id), "SYMBOL"]) %>% 
  mutate(to = data_mapped[match(to, data_mapped$STRING_id), "SYMBOL"]) %>%  
  dplyr::select(from, to , last_col()) %>% 
  dplyr::rename(weight = combined_score)
openxlsx::write.xlsx(links,"PPI.network.xlsx")
# 节点数据
library(tidyverse)
nodes <- links %>% { data.frame(gene = c(.$from, .$to)) } %>% distinct()
# 创建网络图
# 根据links和nodes创建
net <- igraph::graph_from_data_frame(d=links,vertices=nodes,directed = F)
# 添加一些参数信息用于后续绘图
# V和E是igraph包的函数，分别用于修改网络图的节点（nodes）和连线(links)
igraph::V(net)$deg <- igraph::degree(net) # 每个节点连接的节点数
igraph::V(net)$size <- igraph::degree(net)/5 #
igraph::E(net)$width <- igraph::E(net)$weight/10

# 使用ggraph绘图
p1 <- ggraph(net,layout = "stress")+
  geom_edge_fan(aes(edge_width=width), color = "lightblue", show.legend = F)+
  geom_node_point(aes(size=size), color="orange", alpha=0.7)+
  geom_node_text(aes(filter=deg>quantile(igraph::V(net)$deg,0.9) , label=name), size = 5, repel = T)+
  scale_edge_width(range = c(0.2,1))+
  scale_size_continuous(range = c(1,10) )+
  guides(size=F)+
  theme_graph()
save_plot("PPI.network.pdf",p1,base_height = 8,base_width = 14)
