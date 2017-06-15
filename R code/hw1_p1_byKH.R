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
	#g_1a <- erdos.renyi.game(n, p, type=c("gnp"), directed = FALSE)
	g_1a <- random.graph.game(n, p, directed=FALSE)
	print(is.connected(g_1a))
	
	deg.dist = degree_distribution(g_1a)
	plot(deg.dist, xlab="Degree", ylab="Density", main="Problem 1-(a)", col="blue", xlim=c(0,150), ylim=c(0,0.15))
	print(diameter(g_1a, directed=FALSE,unconnected=))
	
	cat ("Press [enter] to continue")
	line <- readline()
	}

#######################################################
##### 1(c)
n = 1000
n_loop = 1000
store_temp = c()

for (i in 1:n_loop) {
	p = 0.01
	flag = 0
	while (flag == 0) {
		p = p - 0.00001
		g_1c <- random.graph.game(n, p, directed=FALSE)
		if (is.connected(g_1c)==FALSE) {
			flag <- 1
			store_temp <- append(store_temp, p)
		}
	}
}
print(mean(store_temp))



