# Part 6

nei_count <- 0
nei_index_list1 <- list()
target_list <- batman_actor_list
for (i in 1:length(movie_act_list)) {
  if (length(intersect(target_list,movie_act_list[[i]])) > 0) {
    nei_count <- nei_count + 1
    nei_index_list1[[nei_count]] <- i
    cat(i,' ',nei_count,'\n')
  }
}


nei_count <- 0
nei_index_list2 <- list()
target_list <- mission_actor_list
for (i in 1:length(movie_act_list)) {
  if (length(intersect(target_list,movie_act_list[[i]])) > 0) {
    nei_count <- nei_count + 1
    nei_index_list2[[nei_count]] <- i
    cat(i,' ',nei_count,'\n')
  }
}

nei_count <- 0
nei_index_list3 <- list()
target_list <- minions_actor_list
for (i in 1:length(movie_act_list)) {
  if (length(intersect(target_list,movie_act_list[[i]])) > 0) {
    nei_count <- nei_count + 1
    nei_index_list3[[nei_count]] <- i
    cat(i,' ',nei_count,'\n')
  }
}




save(nei_index_list1,nei_index_list2,nei_index_list3,file="c:/down/0612 Work/P6_nei.RData")



# Checking Neighbors of Batman, Mission, Minions in Sampled
length(intersect(as.numeric(V(g)$name),nei_index_list1))


length(intersect(as.numeric(V(g)),nei_index_list1))

length(intersect(as.numeric(V(g)$name),nei_index_list2))
length(intersect(as.numeric(V(g)$name),nei_index_list3))

# Movie Index
which(movie_index_list=="Batman v Superman: Dawn of Justice (2016")
which(movie_index_list=="Mission: Impossible - Rogue Nation (2015")
which(movie_index_list=="Minions (2015")

# Copy for new Graph
g2 <- g

# Get Neighbors and Weights
nei_batman_list <- unlist(intersect(as.numeric(V(g)$name),nei_index_list1))
nei_batman_w_list <- vector("integer", length(nei_batman_list))
for (i in 1:length(nei_batman_list)) {
  nei_actor_list <- movie_act_list[[nei_batman_list[[i]]]]
  nei_batman_w_list[[i]] <- length(intersect(nei_actor_list,batman_actor_list)) / (length(nei_actor_list) + length(batman_actor_list) - length(intersect(nei_actor_list,batman_actor_list)))
}

nei_mission_list <- unlist(intersect(as.numeric(V(g)$name),nei_index_list2))
nei_mission_w_list <- vector("integer", length(nei_mission_list))
for (i in 1:length(nei_mission_list)) {
  nei_actor_list <- movie_act_list[[nei_mission_list[[i]]]]
  nei_mission_w_list[[i]] <- length(intersect(nei_actor_list,mission_actor_list)) / (length(nei_actor_list) + length(mission_actor_list) - length(intersect(nei_actor_list,mission_actor_list)))
}

nei_minions_list <- unlist(intersect(as.numeric(V(g)$name),nei_index_list3))
nei_minions_w_list <- vector("integer", length(nei_minions_list))
for (i in 1:length(nei_minions_list)) {
  nei_actor_list <- movie_act_list[[nei_minions_list[[i]]]]
  nei_minions_w_list[[i]] <- length(intersect(nei_actor_list,minions_actor_list)) / (length(nei_actor_list) + length(minions_actor_list) - length(intersect(nei_actor_list,minions_actor_list)))
}

nei_list <- list(nei_batman_list,nei_mission_list,nei_minions_list)
nei_w_list <- list(nei_batman_w_list,nei_mission_w_list,nei_minions_w_list)


# Add Vertices
g2 <- g2 + vertices(as.character(which(movie_index_list=="Batman v Superman: Dawn of Justice (2016")),as.character(which(movie_index_list=="Mission: Impossible - Rogue Nation (2015")),as.character(which(movie_index_list=="Minions (2015")))

# Add Edges
for (i in 1:3) {
  new_vertex <- ""
  if (i == 1) {
    new_vertex <- as.character(which(movie_index_list=="Batman v Superman: Dawn of Justice (2016"))
  } else if (i == 2) {
    new_vertex <- as.character(which(movie_index_list=="Mission: Impossible - Rogue Nation (2015"))
  } else if (i == 3) {
    new_vertex <- as.character(which(movie_index_list=="Minions (2015"))
  }
  for (j in 1:length(nei_list[[i]])) {
    g2 <- g2 + edges(c(new_vertex, as.character(nei_list[[i]][j])), weight=nei_w_list[[i]][j])
  }    
}



# Return Top 5 nearest neighbors
for (i in 1:3) {
  for (j in 1:10) {
    if (length(which(nei_w_list[[i]] == sort(nei_w_list[[i]],decreasing=TRUE)[j])) == 1) {
      m_name <- movie_index_list[[nei_list[[i]][which(nei_w_list[[i]] == sort(nei_w_list[[i]],decreasing=TRUE)[j])]]]
      comm <- fg$membership[which(fg$name == as.character(nei_list[[i]][which(nei_w_list[[i]] == sort(nei_w_list[[i]],decreasing=TRUE)[j])]))]
      cat(j,' ',m_name,'  Community:',comm,'\n')
    } else if (length(which(nei_w_list[[i]] == sort(nei_w_list[[i]],decreasing=TRUE)[j])) > 1) {
      m_name <- movie_index_list[[nei_list[[i]][which(nei_w_list[[i]] == sort(nei_w_list[[i]],decreasing=TRUE)[j])][1]]]
      comm <- fg$membership[which(fg$name == as.character(nei_list[[i]][which(nei_w_list[[i]] == sort(nei_w_list[[i]],decreasing=TRUE)[j])][1]))]
      cat(j,' ',m_name,'  Community:',comm,'\n')
    }
  }
  cat('\n')
}






