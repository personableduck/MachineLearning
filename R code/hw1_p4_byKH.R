#######################################################
# EE232E Homework 1
# 	Kyeong Ho Moon (104617357)
# 	Duckha Hwang   ()
# April 13, 2016
#######################################################


#######################################################
##### 0. load library
library(igraph)

#######################################################
##### 4(a)
n = 1000
fw.prob = 0.37
bw.factor = 0.32/0.37

g_4a <- forest.fire.game(n,fw.prob,bw.factor)
degdist_in <- degree.distribution(g_4a, mode="in")
degdist_out <- degree.distribution(g_4a, mode="out")
plot(degdist_in)
plot(degdist_out)

##### 4(b)
diameter(g_4a)

##### 4(c)
c <- fastgreedy.community(g_4a)
print(modularity(c))



