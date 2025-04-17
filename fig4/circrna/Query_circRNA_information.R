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
ln = load("circrna.Rdata")

# read in conf
query_type = conf$query_type
query_name = conf$query_name
show_seq = conf$show_circRNA_seq


#### analysis ####
info.table = "No results match your inquery."
if(query_type == "circRNA_ID"){
  info.table = anno[anno$Circ_id == query_name, ]
}else{
  info.table = anno[anno$Gene.Symbol == query_name, ]
}

if(show_seq){
  ln2 = load("circrna_seq.Rdata")
  
  getSeq <- function(seqid){
    paste0(circrna.seq[[seqid]],collapse = "")
  }
  
  seqs <- (sapply(info.table$Circ_id,getSeq))
  
  info.table$Sequence = seqs
  
}

write.xlsx(info.table,"Query_circRNA_information.xlsx")
