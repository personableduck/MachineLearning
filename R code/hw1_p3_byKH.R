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
##### 3(a) and (b)
n = 1000

for (beta in 0:5) {
	g_3a <- aging.prefatt.game(n, pa.exp=1, aging.exp=beta, aging.bin=n, directed=FALSE)
	plot(degree.distribution(g_3a),log='xy')
	cat ("Press [enter] to continue")
	line <- readline()
	
	cl <- clusters(g_3a)
	gccIndex = which.max(cl$csize)
	nonGccNodes <- (1:vcount(g_3a))[cl$membership != gccIndex]
	gcc <- delete.vertices(g_3a, nonGccNodes)
	fg <- fastgreedy.community(gcc)
	cmsize <- sizes(fg)
	cmsize <- as.vector(sizes(fg))
	print(modularity(fg))
}