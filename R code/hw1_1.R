#######################################################
# EE232E Homework 1
# 	Kyeong Ho Moon (104617357)
# 	Duckha Hwang   ()
# April 13, 2016
#######################################################


#######################################################
##### 0. load data
library(igraph)

#######################################################
##### 1(a) and (b)
n = 1000 # number of nodes
p_vec = c(0.01, 0.05, 0.1) # vector of the probability for drawing an edge between two arbitrary vertices

for (p in p_vec) {
	g_1a <- erdos.renyi.game(n, p, type=c("gnp"), directed = FALSE)
	deg.dist = degree_distribution(g_1a)
	plot(deg.dist)
	print(is.connected(g_1a))
	print(diameter(g_1a, directed=FALSE,unconnected=))
	
	cat ("Press [enter] to continue")
	line <- readline()
	}

#######################################################
##### 1(c)
p = 0.01
count = 100
n_loop = 20

while (count > 0) {
	p = p - 0.001
	count = 0
	for (i in 1:n_loop) {
		g_1c <- erdos.renyi.game(n, p, type=c("gnp"), directed=FALSE)
		if (is.connected(g_1c)) {
			count = count + 1
		}
	}
	print(c(p, count))
}
print(p)



