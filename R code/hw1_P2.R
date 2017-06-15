# Problem 2 (a)
# Generate Random Network for Power law degree distribution

g.big.ba = barabasi.game(1000,directed = FALSE)

fit_power_law = function(graph) {
  # calculate degree
  d = degree(graph, mode = "all")
  dd = degree.distribution(graph, mode = "all", cumulative = FALSE)
  degree = 1:max(d)
  density = dd[-1]
  # delete blank values
  nonzero.position = which(density != 0)
  density = density[nonzero.position]
  degree = degree[nonzero.position]
  reg = lm(log(density) ~ log(degree))
  cozf = coef(reg)
  power.law.fit = function(x) exp(cozf[[1]] + cozf[[2]] * log(x))
  alpha = -cozf[[2]]
  R.square = summary(reg)$r.squared
  print(paste("Alpha =", round(alpha, 3)))
  print(paste("R square =", round(R.square, 3)))
  # plot
  plot(density ~ degree, log = "xy", xlab = "Degree", ylab = "Density", 
       col = 1, main = "Degree Distribution")
  curve(power.law.fit, col = "red", add = T, n = length(d))
}


fit_power_law(g.big.ba)


#---------------------------------

# Problem 2 (a)
# Generate Random Network for Power law degree distribution
power = -3
numnodes = 1000  # Change to 10000
g<-barabasi.game(numnodes, power, directed = FALSE)
plot(degree.distribution(g),type="o")
plot(degree.distribution(g),log='xy',type="o")

# Calculate Diameter
diam_unconnected_true = diameter(g.big.ba, unconnected=TRUE)
diam_unconnected_false = diameter(g.big.ba, unconnected=FALSE)



# Problem 2 (b),(c)
# Clustering Graph
clusterlist <- clusters(g)
size <- clusterlist$csize
numcluster <- clusterlist$no

# Find Community structure of Overall Graph
communityobj <- fastgreedy.community(g, merges=TRUE, modularity=TRUE, membership=TRUE, weights=NULL)
# Modularity & Sizes of Communities
modularityval <- modularity(communityobj)
sizeofmembers <- sizes(communityobj)



# Find the GCC
GCCindex <- which.max(size)
nonGCCnodes<-(1:vcount(g))[clusterlist$membership!=GCCindex] #compares the membership vector with the index , then multiply TRUE with node list
GCCfinal <- delete.vertices(g,nonGCCnodes)

# Find Community structure of GCC
communityobjGCC = fastgreedy.community(GCCfinal, merges=TRUE, modularity=TRUE, membership=TRUE, weights=NULL)
# Modularity & Sizes of Communities
modularityvalGCC = modularity(communityobjGCC)
sizeofmembersGCC= sizes(communityobjGCC)

print(modularityval)
print(modularityvalGCC)
print(sizeofmembers)
print(sizeofmembersGCC)


# Problem 2 (d)
iter = 100
# Pick up random vertices iter times
randv<-igraph.sample(1, numnodes, iter)

# For each randv,pick up a random adjacent vertex
# So first calculate the degree of the random vertex , based on vertex iD
degrandv<-degree(g, randv)

# Find adjacent vertex IDs
adjacentv = g[[randv,]]     # or use adjcentv = get.adjlist(g)  #adjacentv is a vector 

# Find degree values of random neighbor nodes
nei_deg <- c(0)
for (i in 1:length(adjacentv)) {
	nei_deg[i] <- degree(g, adjacentv[[i]][sample(1:length(adjacentv[[i]][]),1)])
}

# Make Plot the Degree Distribution
h <- hist(nei_deg, breaks = 10, plot=FALSE)
h$counts = h$counts/sum(h$counts)
plot(h, ylab="degree.distribution")

