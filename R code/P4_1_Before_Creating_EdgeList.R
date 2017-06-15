# Part4 - Prepare for Creating Edge List by Movie

# Prepare Movie List for Part4
unlist_movie <- unlist(movie_list)
# Use Unique List to Index
movie_index_list <- unique(unlist_movie)
movie_index_list <- as.character(movie_index_list)
movie_index_list <- movie_index_list[c(-(which(movie_index_list=="")))]
sort_order <-  sort.list(movie_index_list, decreasing=FALSE, method=c("shell"))
movie_index_list <- as.list(rbind(movie_index_list)[,sort_order])

movie_index_data_frame <- as.data.frame(table(movie_index_list)) # Sorted
movie_index_list <- as.list(as.character(movie_index_data_frame[,1]))




duplicate_unlist_movie <- unlist_movie[duplicated(unlist_movie)]
unique_duplicate_movie <- unique(duplicate_unlist_movie)
count_duplicate_movie <- list()
total_number <- length(unique_duplicate_movie)

count_duplicate_movie <- as.data.frame(table(duplicate_unlist_movie))
movie_threshold <- 10
movie_nodes <- subset(count_duplicate_movie, Freq > movie_threshold) 
# 1st Row was 'NA' -> Remove
movie_nodes <- movie_nodes[c(-1),]



# Get Movie Genre File
genre_list <- list()
conn_genre=file("C:/project_2_data/movie_genre.txt",open="r")
genre_list <- strsplit(readLines(conn_genre),"\t\t")
close(conn_genre)

mat_genre <- do.call(rbind, genre_list)
genre_data_frame <- as.data.frame(mat_genre)

# Remove NA Row
genre_data_frame <- genre_data_frame[complete.cases(genre_data_frame),]
# movie_index_list <- as.list(as.character(genre_data_frame[,1]))
genre_list <- as.list(as.character(genre_data_frame[,2]))

# Get Unique Genre -> Use Index
genre_index_list <- unique(genre_list)



# Make Basic Movie List
p4_movie_list <- as.list(as.character(movie_nodes[,1]))
p4_movie_freq_list <-  as.list(as.character(movie_nodes[,2]))



#movie_index_list2 <- movie_index_list[1:100]
#movie_index_list2 <- sort(movie_index_list)



start_movie_index <- 1
j <- 1
if (exists("p4_movie_index_list")) {
  start_movie_index <- length(p4_movie_index_list) + 1
  j <- p4_movie_index_list[[length(p4_movie_index_list)]]
} else {
  p4_movie_index_list <- list()
}
which(movie_index_list == p4_movie_list[[1]])

percent <- 0
total_number <- length(p4_movie_list)

for (i in start_movie_index:length(p4_movie_list)) {
  while (TRUE) {
    if (p4_movie_list[[i]] == movie_index_list[[j]]) {
      p4_movie_list[[i]] <- j
      break
    } else {
      if (j == length(movie_index_list)) {
        cat('No Index Error ',i,'\n')
        break
      }
    }
    j <- j + 1
  }
  
  if (trunc(i/total_number*10000) > percent) {
    percent <- trunc(i/total_number*10000)
    cat((percent/100),'%','\n')
  } else {
    
  }
}
p4_movie_list[[1]]
