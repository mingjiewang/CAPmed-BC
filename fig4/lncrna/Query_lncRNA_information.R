#!/usr/bin/env Rscript
#### Input args ####
library(jsonlite)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No JSON input provided.")
}
json_input <- args[1]
# json_input <- "Query_circRNA_information.json"
print(json_input)
conf <- fromJSON(json_input)

#### pkgs and data ####
library(openxlsx)
# load data
ln = load("lncrna.Rdata")

# read in conf
query_type = conf$query_type
query_name = conf$query_name
show_seq = conf$show_lncRNA_seq


#### analysis ####
info.table = "No results match your inquery."
if(query_type == "lncRNA_ID"){
  info.table = anno[anno$LncRNA_id == query_name, ]
}else{
  info.table = na.omit(anno[anno$target_gene == query_name, ])
}


if(show_seq){
  ln2 = load("lncrna_seq.Rdata")
  
  getSeq <- function(seqid){
    paste0(lncrna.seq[[seqid]],collapse = "")
  }
  
  seqs <- (sapply(info.table$LncRNA_id,getSeq))
  
  info.table$Sequence = toupper(seqs)
  
}

write.xlsx(info.table,"Query_lncRNA_information.xlsx")
