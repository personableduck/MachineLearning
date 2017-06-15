# Combining Edge Files
edge_list <- list()
for (part in 1:114) {
  part_variable <- paste0("edge_part",part)
  tmp_part <- get(part_variable)
  edge_list <- append(edge_list,tmp_part)
  cat('Part',part,' Combining\n')
  if (part == 114) {
    cat('Complete\n')
  }
  rm(tmp_part)
}

for (part in 1:114) {
  part_variable <- paste0("edge_part",part)
  eval(parse(text=paste0("rm(",part_variable,")")))  
}
rm(part)
rm(part_variable)



# Creating Edge List Text File
sink(file = "c:/down/edge_list.txt", append = FALSE)
writeLines(unlist(lapply(edge_list, paste, collapse="\t")))
sink()

# Create Graph
g <- read.graph("C:/down/edge_list.txt", directed=TRUE, format="ncol")
length(E(g)) 
length(V(g))