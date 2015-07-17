x <- matrix(ncol=2, byrow=TRUE, c(
  c(0,0),
   c(1,0),
   c(0,1),
   c(-1,0),
   c(0,-1)
))

x <- x[sample(1:nrow(x)),]

euclidianDistance <- function(x,y)
{
  sqrt(sum((x-y)^2))
}

space <- lapply(apply(x, 1, as.list), as.numeric)

h <- hclust(dist(x), 'single')
print(h$merge)
plot(h)

hclust2(euclidianDistance, space)
