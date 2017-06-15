# Part 5

# Edge List Was too Big -> Sampling!!
sample_movie_edge_3000000 <- sample(m_edge_list,3000000)

sink(file = "c:/down/0612 Work/sample_edge_list_3000000.txt", append = FALSE)
writeLines(unlist(lapply(sample_movie_edge_3000000, paste, collapse="\t")))
sink()
rm(sample_movie_edge_3000000)


g <- read.graph("c:/down/0612 Work/sample_edge_list_3000000.txt", directed=FALSE, format="ncol")

g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = getIgraphOpt("edge.attr.comb"))
is.simple(g)
fg <- fastgreedy.community(g)

length(V(g))
length(E(g))

for (i in 1:length(genre_movie_list)) {
  genre_movie_list[[i]] <- strsplit(genre_movie_list[[i]],")")[[1]]
}


# ith Community members index return
#as.numeric(fg$names[which(fg$membership == i)])
#unlist(movie_index_list[as.numeric(fg$names[which(fg$membership == i)])])
#unlist(genre_number_list[which(unlist(movie_index_list[as.numeric(fg$names[which(fg$membership == i)])]) %in% genre_movie_list)])


genre_number_list <- vector("list", length(genre_list))
for (i in 1:length(genre_list)) {
  genre_number_list[[i]] <- which(genre_index_list %in% genre_list[[i]])
}


# Find communities which sizes are over 100
length(which(sizes(fg) > 100))
#which(sizes(fg) > 100)
sizes(fg)[which(sizes(fg) > 100)]

# Get community number from the Biggest size to lower
community_order_list <- list()
for (i in 1:length(which(sizes(fg) > 100))) {
  community_order_list[[i]] <- which(sizes(fg) == sort(sizes(fg)[which(sizes(fg) > 100)],decreasing=TRUE)[[i]])[[1]]
}


community_genre_list <- list()
for (i in 1:length(which(sizes(fg) > 100))) {
  community_genre_list[[i]] <- unlist(genre_number_list[which(unlist(movie_index_list[as.numeric(fg$names[which(fg$membership == community_order_list[[i]])])]) %in% genre_movie_list)])
}

for (i in 1:length(which(sizes(fg) > 100))) {
  hist(community_genre_list[[i]],29)
  hist(community_genre_list[[i]][which(community_genre_list[[i]] != 4)],29)
}

