# Part4 - 2 Create Act List
#
# Create Act Unlist Matching to Unlist Movie
movie_act_list <- movie_list
percent <- 0
total_number <- length(movie_act_list)
for (i in 1:length(movie_act_list)) {
  name <- name_list[[i]]
  for (j in 1:length(movie_act_list[[i]])) {
    movie_act_list[[i]][j] <- name
  }
  if (trunc(i/total_number*10000) > percent) {
    percent <- trunc(i/total_number*10000)
    cat((percent/100),'%','\n')
  }
}

unlist_movie_act <- unlist(movie_act_list)




