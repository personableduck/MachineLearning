# Part 3


# Create a Graph
g <- read.graph("C:/down/edge_list.txt", directed=TRUE, format="ncol")
length(E(g)) 
length(V(g))

# Page Rank to the Graph
pr <- page.rank (g, algo = "prpack", directed = TRUE, damping = 0.85, weights = NULL, options = NULL)

# Finding Top 10 of PageRank
top10_act_list <- list()
for (i in 1:10) {
  # Get Vertex ID by ordering highest pagerank probability
  idx <- as.numeric(V(g)$name[which(pr$vector == sort(pr$vector, TRUE)[i])])
  # Find Name of ID
  top10_act_list[i] <- name_list[[idx]]
  cat(i,'\t',name_list[[idx]],'\t',sort(pr$vector, TRUE)[i],"\tMovies:",length(movie_list[[idx]]),'\n')
}



#Finding by name in PageRank
 
which(name_list == "Downey Jr., Robert")
movie_list[[which(name_list == "Downey Jr., Robert")]]
pr_idx <- which(attributes(pr$vector)$names == as.character(which(name_list == "Downey Jr., Robert")))
pr$vector[pr_idx]


which(name_list == "Diesel, Vin")
movie_list[[which(name_list == "Diesel, Vin")]]
pr_idx <- which(attributes(pr$vector)$names == as.character(which(name_list == "Diesel, Vin")))
pr$vector[pr_idx]




# Searching Actor/Actress Name by order of pagerank probability
for (i in 401:500) {
  idx <- as.numeric(V(g)$name[which(pr$vector == sort(pr$vector, TRUE)[i])])
  cat(i,'\t',name_list[[idx]],'\t',sort(pr$vector, TRUE)[i],"\tMovies:",length(movie_list[[idx]]),'\n')
}




# Significant Paring -> High Weight !! because Many Intersect Number
# Sort by Weight
# get Edge List to Matrix, Weight Round by 3
edge_mat <- cbind( get.edgelist(g) , round( E(g)$weight, 3 ))
sort_edge_mat <- edge_mat[order(edge_mat[,3],decreasing=TRUE),] 

for (i in 1:nrow(sort_edge_mat)) {
  if (sort_edge_mat[i,3] != 1) {
    cat(i)
    break
  }
}
length(which(E(g)$weight == 1))

# Top Paring Actors/Actresses
for (i in 1:100) {
  name1 <- name_list[[as.numeric(sort_edge_mat[i,1])]]
  name2 <- name_list[[as.numeric(sort_edge_mat[i,2])]]
  
  cat(i,' ',name1,'   ',name2,'\n')
}
