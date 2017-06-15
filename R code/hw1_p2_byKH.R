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
##### 2(a)
n = 1000
g_2a <- barabasi.game(n,directed=FALSE)
deg.dist = degree_distribution(g_2a)
plot(deg.dist, log='x')
print(diameter(g_2a, directed=FALSE, unconnected=is.connected(g_2a)))


##### 2(b)
is.connected(g_2a)
cl <- clusters(g_2a)
gccIndex = which.max(cl$csize)
nonGccNodes <- (1:vcount(g_2a))[cl$membership != gccIndex]
gcc <- delete.vertices(g_2a, nonGccNodes)

fg <- fastgreedy.community(gcc)
cmsize <- sizes(fg)
cmsize <- as.vector(sizes(fg))
print(modularity(g_2a,fg$membership))

##### 2(c)
n = 10000
g_2c <- barabasi.game(n,directed=FALSE)
is.connected(g_2c)
cl <- clusters(g_2c)
gccIndex = which.max(cl$csize)
nonGccNodes <- (1:vcount(g_2c))[cl$membership != gccIndex]
gcc <- delete.vertices(g_2c, nonGccNodes)
fg <- fastgreedy.community(gcc)
print(modularity(g_2c,fg$membership))

##### 2(d)



