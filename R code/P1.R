actor_list <- list()
actress_list <- list()

conn_actor=file("C:/Users/1234/Dropbox/graph2016/Project2/data/actor_movies.txt",open="r")
actor_list <- strsplit(readLines(conn_actor),"\t\t")
close(conn_actor)

conn_actress=file("C:/Users/1234/Dropbox/graph2016/Project2/data/actress_movies.txt",open="r")
actress_list <- strsplit(readLines(conn_actress),"\t\t")
close(conn_actress)

total_list <- append(actor_list,actress_list)

name_list<- list()
movie_list <- list()

start <- 1
percent <- 0
if (length(name_list) == 0 && length(movie_list) == 0) {
  name_list <- list()
  movie_list <- list()
} else {
  name_length <- length(name_list)
  movie_length <- length(movie_list)
  if (movie_length == name_length) {
    start <- length(name_list)+1
  } else {
    name_list <- correction_func(name_list,1)
    movie_list <- correction_func(movie_list,2)
    name_length <- length(name_list)
    movie_length <- length(movie_list)
    start <- length(name_list)+1
  }
  
}
t_number <- length(total_list)
for (i in start:length(total_list)) {
  #name_list[[length(name_list)+1]] <- total_list[[i]][1]
  #movie_list[[length(movie_list)+1]] <- total_list[[i]][2:length(total_list[[i]])]
  
  name_list <- append(name_list,total_list[[i]][1])
  movie_list <- append(movie_list,list(total_list[[i]][2:length(total_list[[i]])]))
  
  if (trunc(i/t_number*10000) > percent) {
    percent <- trunc(i/t_number*10000)
    cat((percent/100),'%','\n')
  } else {
    
  }
  
}

name_list[[length(name_list)]] <- NULL
movie_list[[length(movie_list)]] <- NULL

correction_func <- function(list,flag) {
  if (flag == 1) {
    # Correct for Name
    percent <- 0
    t_number <- length(total_list)
    for (i in 1:length(total_list)) {
      
      if (identical(total_list[[i]][2:length(total_list[[i]])],list[[i]]) == 'FALSE') {
        first_part_tmp <- list[1:(i-1)]
        second_part_tmp <- list[i:length(list)]
        first_part_tmp <- list[1:(i-1)]
        second_part_tmp <- list[i:length(list)]
        first_part_tmp <- append(first_part_tmp,total_list[[i]][1])
        first_part_tmp <- append(first_part_tmp,second_part_tmp)
        list <- first_part_tmp
        rm(first_part_tmp)
        cat(i,' Error Correct\n')
      }
      if (trunc(i/t_number*10000) > percent) {
        percent <- trunc(i/t_number*10000)
        cat((percent/100),'%','\n')
      } else {
        
      }
    }
  } else if (flag == 2) {
    # Correct for Movie
    percent <- 0
    t_number <- length(total_list)
    for (i in 1:length(total_list)) {
      
      if (identical(total_list[[i]][2:length(total_list[[i]])],list[[i]]) == 'FALSE') {
        first_part_tmp <- list[1:(i-1)]
        second_part_tmp <- list[i:length(list)]
        first_part_tmp <- list[1:(i-1)]
        second_part_tmp <- list[i:length(list)]
        first_part_tmp <- append(first_part_tmp,list(total_list[[i]][2:length(total_list[[i]])]))
        first_part_tmp <- append(first_part_tmp,second_part_tmp)
        list <- first_part_tmp
        rm(first_part_tmp)
        cat(i,' Error Correct\n')
      }
      if (trunc(i/t_number*10000) > percent) {
        percent <- trunc(i/t_number*10000)
        cat((percent/100),'%','\n')
      } else {
        
      }
    }
  }
}



# Erase NA in movie_list
movie_list2 <- movie_list
if (exists(i_NA)) {
  continue_NA <- i_NA
} else {
  continue_NA <- 1
}
for (i_NA in continue_NA:length(movie_list2)) {
  if (length(which(is.na(movie_list2[[i_NA]]))) == 1) {
    movie_list2[[i_NA]] <- movie_list2[[i_NA]][-which(is.na(movie_list2[[i_NA]]))]
  } else {
    
  }
}
cat("Complete")


# Filtering Movie List -> Extract Only Title (Some includes Role Name after year)
total_count <- length(movie_list)
percent <- 0
start_row <- 1
start_col <- 1
if (exists("m_row") && exists("m_col")) {
  start_row <- m_row
  start_col <- m_col
}
for (m_row in 1:length(movie_list)) {
  for (m_col in 1:length(movie_list[[m_row]])) {
    movie_list[[m_row]][m_col] <- unlist(strsplit(movie_list[[m_row]][m_col], ")"))[1]
  }
  
  if (trunc(m_row/total_count*10000) > percent) {
    percent <- trunc(m_row/total_count*10000)
    cat((percent/100),'%','\n')
  } else {
    
  }
}

# Extract Actor/Actress Have 10 Movies
total_list <- append(actor_list,actress_list)
extracted_name_list <- list()
for (i in 1:length(total_list)) {
  if (length(total_list[[i]]) >= 11) {
    extracted_name_list[[length(extracted_name_list)+1]] <- total_list[[i]][1]
    #extracted_movie_list[[length(extracted_movie_list)+1]] <- total_list[[i]][2:length(total_list[[i]])]
  }  
}

extracted_movie_list <- vector("list", length(extracted_name_list))
percent <- 0
t_number <- length(extracted_name_list)
for (i in 1:length(extracted_name_list)) {
  extracted_movie_list[[i]] <- movie_list[[extracted_name_index_list[[i]]]]
  if (trunc(i/t_number*10000) > percent) {
    percent <- trunc(i/t_number*10000)
    cat((percent/100),'%','\n')
  }
}