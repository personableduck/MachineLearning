# Part 4 - Creating Movie Edges
#
# Movie Edge List Creating

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
  intersect(list1,list2)
}

part_start <- 1
part_end <- 10

for (part in part_start:part_end) {
  part_variable <- paste0("movie_edge_part",part)
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
  if (part * 1000 > length(p4_movie_list)) {
    # Last Part
    end_row <- length(p4_movie_list)
  } else {
    end_row <- part * 1000
  }
  
  part_start_row <- start_row
  part_end_row <- end_row
  total_count <- count_all_function(part_start_row,part_end_row,length(p4_movie_list))
  
  
  
  
  part_variable <- paste0("movie_edge_part",part)
  
  tmp_part <- list()
  
  continue <- 0
  
  
  part_variable_name <- part_variable
  
  
  
  if (exists(part_variable_name)) {
    
    continue <- 1
    
    tmp_part <- get(part_variable_name)
    start_row <- as.numeric(which(p4_movie_index_list==tmp_part[[length(tmp_part)]][1]))
    start_col <- as.numeric(which(p4_movie_index_list==tmp_part[[length(tmp_part)]][2]))
    
    
  } else {
    #assign(part_variable,list())
    start_row <- part_start_row
    start_col <- start_row + 1
  }
  
  if (continue == 1) {
    current_count <- count_all_function(part_start_row,(start_row-1),length(p4_movie_list)) + (start_col-start_row)
  } else {
    current_count <- 1
  }
  
  
  percent <- 0
  for (i in start_row:end_row) {
    from_index <- p4_movie_index_list[[i]]
    if (continue == 1) {
      # Continue Work
      for (j in start_col:length(p4_movie_list)) {
        intersect_list <- r_intersect_function(p4_movie_act_list[[i]],p4_movie_act_list[[j]])
        if (length(intersect_list) > 0) {
          to_index <- p4_movie_index_list[[j]]
          
          
          tmp_part <- append(tmp_part,list(c(from_index,to_index,length(intersect_list)/(length(p4_movie_act_list[[i]])+length(p4_movie_act_list[[j]])-length(intersect_list)))))
        }
        current_count <- current_count + 1
      }
    } else {
      # Work to New Row
      for (j in (start_row+1):length(p4_movie_list)) {
        intersect_list <- r_intersect_function(p4_movie_act_list[[i]],p4_movie_act_list[[j]])
        if (length(intersect_list) > 0) {
          to_index <- p4_movie_index_list[[j]]
          
          tmp_part <- append(tmp_part,list(c(from_index,to_index,length(intersect_list)/(length(p4_movie_act_list[[i]])+length(p4_movie_act_list[[j]])-length(intersect_list)))))
        }
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
}

