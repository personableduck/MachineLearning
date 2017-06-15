library ("igraph")

# Problem 1 (a)
p <- 0.01  # Change p value to 0.05, 0.1
g <- random.graph.game(1000, p, directed=F)
is.connected(g)
plot(degree.distribution(g), xlab="Degree",ylab="Density",main="problem 1-(a)", col="red",xlim=c(0,150),ylim=c(0,0.15))
par(new=TRUE)
p <- 0.05  # Change p value to 0.05, 0.1
g <- random.graph.game(1000, p, directed=F)
is.connected(g)
plot(degree.distribution(g), xlab="Degree",ylab="Density",main="problem 1-(a)", col="blue",xlim=c(0,150),ylim=c(0,0.15))
par(new=TRUE)
p <- 0.1  # Change p value to 0.05, 0.1
g <- random.graph.game(1000, p, directed=F)
is.connected(g)
plot(degree.distribution(g), xlab="Degree",ylab="Density",main="problem 1-(a)", col="green",xlim=c(0,150),ylim=c(0,0.15))

# Problem 1 (b)
p <- 0.005  # Change p value to 0.05, 0.1
g <- random.graph.game(1000, p, directed=F)

connectedindicator = is.connected(g) #check connected or not

diam_unconnected_false= diameter(g, unconnected=FALSE)
diam_unconnected_true = diameter(g, unconnected=TRUE)

# Problem 1 (c)
# Run (b) Several Times by changing 'p'
p <- 0.01  # Change p value to 0.05, 0.1
g <- random.graph.game(1000, p, directed=F)
is.connected(g)
plot(degree.distribution(g), xlab="Degree",ylab="Density",main="problem 1-(a)", col="red",xlim=c(0,150),ylim=c(0,0.15))

p = 0.005
count = 4

while (count > 0) {
  p = p + 0.001
  g_1c <- random.graph.game(1000, p, directed=F)
 if (is.connected(g_1c)) {
      count = count - 1
    }
  print(c(p, count))
}
print(p)





