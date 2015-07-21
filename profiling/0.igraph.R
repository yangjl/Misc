
# https://nsaunders.wordpress.com/2010/04/21/experiments-with-igraph/

library(igraph)
demo(package="igraph")



sgc <- spinglass.community(g)
V(g)$membership <- sgc$membership
# found 4 communities 0, 1, 2, 3
V(g) [ membership == 0 ]$color <- "cyan"
V(g) [ membership == 1 ]$color <- "green"
V(g) [ membership == 2 ]$color <- "blue"
V(g) [ membership == 3 ]$color <- "red"
V(g)$size <- 4
V(g) [ name == "the-life-scientists" ]$size <- 20
> png(filename = "tls.png", height = 800, width = 800)
> plot(g, layout=layout.fruchterman.reingold, vertex.color=V(g)$color, vertex.size = V(g)$size, vertex.label = NA, edge.arrow.size = 0.5)
> dev.off()


