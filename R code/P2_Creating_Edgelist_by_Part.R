# Edge List Creating

if (exists("part_variable") && exists("tmp_part")) {
  assign(part_variable,tmp_part)
}

count_all_function <- function(start_row,end_row,row_length) {
  if (end_row < start_row) {
    0
  } else {
    (end_row - start_row +1) * row_length - sum(start_row:end_row)
  }
  
}

r_intersect_function <- function(list1,list2) {
  #Reduce(intersect,list(list1,list2))
  intersect(list1,list2)
}

#####################################
# Set the Start Part and End Part
part_start <- 11
part_end <- 20
#####################################


for (part in part_start:part_end) {
  part_variable <- paste0("edge_part",part)
  part_variable_name <- part_variable
  if (exists(part_variable_name)) {
    part_start <- part
  } else {
    if (part == part_start) {
      part_start <- part
    } else {
      part_start <- part - 1
    }
    
    break
  }
}



for (part in part_start:part_end) {
  #part <- 1
  start_row <- ((part-1) * 1000) + 1
  if (part * 1000 > length(extracted_name_list)) {
    # Last Part
    end_row <- length(extracted_name_list)
  } else {
    end_row <- part * 1000
  }
  
  part_start_row <- start_row
  part_end_row <- end_row
  total_count <- count_all_function(part_start_row,part_end_row,length(extracted_name_list))
  
  
  
  
  part_variable <- paste0("edge_part",part)
  
  tmp_part <- list()
  
  continue <- 0
  
  
  part_variable_name <- part_variable
  
  
  
  if (exists(part_variable_name)) {
    
    continue <- 1
    
    tmp_part <- get(part_variable_name)
    start_row <- as.numeric(which(extracted_name_index_list==tmp_part[[length(tmp_part)]][1]))
    start_col <- as.numeric(which(extracted_name_index_list==tmp_part[[length(tmp_part)]][2]))
    
    if (start_row < start_col) {
      # Incomplete finish -> Erase it and continue
      tmp_part[[length(tmp_part)]] <- NULL
      start_row <- as.numeric(which(extracted_name_index_list==tmp_part[[length(tmp_part)]][2]))
      start_col <- as.numeric(which(extracted_name_index_list==tmp_part[[length(tmp_part)]][1]))
    } else {
      # Complete finish -> Continue to next
      start_row <- as.numeric(which(extracted_name_index_list==tmp_part[[length(tmp_part)]][2]))
      start_col <- as.numeric(which(extracted_name_index_list==tmp_part[[length(tmp_part)]][1]))
      if (start_col == length(extracted_name_list)) {
        # End of Row -> Go to Next Row
        start_row <- start_row + 1
        start_col <- start_row + 1
      } else {
        # Current Row Continue
        start_col <- start_col + 1
      }
    }
    
  } else {
    #assign(part_variable,list())
    start_row <- part_start_row
    start_col <- start_row + 1
  }
  
  if (continue == 1) {
    current_count <- count_all_function(part_start_row,(start_row-1),length(extracted_name_list)) + (start_col-start_row)
  } else {
    current_count <- 1
  }
  
  
  percent <- 0
  for (i in start_row:end_row) {
    from_index <- extracted_name_index_list[[i]]
    if (continue == 1) {
      # Continue Work
      for (j in start_col:length(extracted_name_list)) {
        intersect_list <- r_intersect_function(extracted_movie_list[[i]],extracted_movie_list[[j]])
        if (length(intersect_list) > 0) {
          to_index <- extracted_name_index_list[[j]]
          
          tmp_part <- append(tmp_part,list(c(from_index,to_index,length(intersect_list)/length(extracted_movie_list[[i]])),c(to_index,from_index,length(intersect_list)/length(extracted_movie_list[[j]]))))
        }
        #rm(intersect_list)
        current_count <- current_count + 1
      }
    } else {
      # Work to New Row
      for (j in (start_row+1):length(extracted_name_list)) {
        intersect_list <- r_intersect_function(extracted_movie_list[[i]],extracted_movie_list[[j]])
        if (length(intersect_list) > 0) {
          to_index <- extracted_name_index_list[[j]]
          
          tmp_part <- append(tmp_part,list(c(from_index,to_index,length(intersect_list)/length(extracted_movie_list[[i]])),c(to_index,from_index,length(intersect_list)/length(extracted_movie_list[[j]]))))
        }
        #rm(intersect_list)
        current_count <- current_count + 1
      }
    }
    
    if (trunc(current_count/total_count*10000) > percent) {
      percent <- trunc(current_count/total_count*10000)
      cat(part_variable,'  ',(percent/100),'%','\n')
    } else {
      
    }
  }
  
  
  
  # When Part is Complete Save to the Part List
  assign(part_variable,tmp_part)
  #rm(tmp_part)
}
