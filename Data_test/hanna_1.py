# -*- coding: utf-8 -*-
"""
Created on Tue May 24 19:12:58 2016

@author: 1234
"""

import numpy as np
import csv

actor_list = []
actress_list = []

with open('C:/Users/Kenneth/Dropbox/graph2016/Project2/data/actor_movies.txt') as f:
	reader_actor = csv.reader(f, delimiter='\t')
 
actor_list = reader_actor.split()

with open('/Users/Kenneth/Dropbox/graph2016/Project2/data/actress_movies.txt') as j:
     reader_actress = csv.reader(j, delimiter='\t')
     
actress_list = reader_actor.split()     

total_list= actor_list.append(actress_list)

name_list = []
movie_list = []

start = 1
percent = 0

if len(name_list) == 0 & len(movie_list) == 0:
    name_list = []
    movie_list = []
else:
    name_length = len(name_list)
    movie_length = len(movie_list)

    if movie_length == name_length :
        start = len(name_list)+1
    else :
        name_list = correction_func(name_list,1)
        movie_list = correction_func(movie_list,2)
        name_length = len(name_list)
        movie_length = len(movie_list)
        start = len(name_list)+1


t_number = len(total_list)

for i in len(total_list):
    #name_list[[length(name_list)+1]] <- total_list[[i]][1]
    #movie_list[[length(movie_list)+1]] <- total_list[[i]][2:length(total_list[[i]])]
  
    name_list = name_list.append(total_list[[i]][1])
    movie_list = movie_list.append(list(total_list[[i]][2:len(total_list[[i]])]))
  
if trunc(i/t_number*10000) > percent:
    percent = trunc(i/t_number*10000)
    print ((percent/100),'%','\n')

      
name_list[[len(name_list)]] = None
movie_list[[len(movie_list)]] = None

#def 

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
