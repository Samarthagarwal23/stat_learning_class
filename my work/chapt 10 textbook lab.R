# chapt 10 textbook lab - p. 401 (pdf 416)

(states <- row.names(USArrests))
names(USArrests)

summary(USArrests)

apply(USArrests, 2, mean)
apply(USArrests, 2, var)

# since variances are so wide, would influance PCA - must scale!

pr.out <- prcomp(USArrests, scale = TRUE)

names(pr.out)
pr.out$center # means of variables before scaling
pr.out$scale # stdev prior to scaling

pr.out$rotation # each column gives loading
dim(pr.out$x) # columns are principal component score vectors

biplot(pr.out, scale = 0) # plot first 2 principal components

# reproduce fig 10.1
pr.out$rotation <- -pr.out$rotation
pr.out$x <- -pr.out$x
biplot(pr.out, scale = 0)

pr.out$sdev # stdev of each principal component

(pr.var <- pr.out$sdev^2) # variance

(pve <- pr.var/sum(pr.var)) # portion of variance explained by each pr. component

# plot portion of varience explained and cumultive PVE
plot(pve, xlab = "Prinicipal Component", 
     ylab = "Proportion of Variance Explained", ylim = c(0,1),
     type = 'b')

plot(cumsum(pve), xlab = "Prinicipal Component", 
     ylab = "Cumulative Proportion of Variance Explained", 
     ylim = c(0,1), type = 'b')

## clustering
set.seed(2)
x <- matrix(rnorm(50*2), ncol = 2)

x[1:25,1] <- x[1:25,1] + 3 #first 25 obs have a mean shift relative to the next 25
x[1:25,2] <- x[1:25,2] -4 #first 25 obs have a mean shift relative to the next 25

km.out <- kmeans(x, 2, nstart = 20)
km.out$cluster

plot(x, col = (km.out$cluster +1), 
     main = "K-means Clustering Results with K=2", xlab = "",
     ylab = "", pch = 20, cex = 2)

# what if we use 3 clusters?
set.seed(4)
km.out <- kmeans(x, 3, nstart = 20)
km.out

plot(x, col = (km.out$cluster +1), 
    main = "K-means Clustering Results with K=3", xlab = "",
    ylab = "", pch = 20, cex = 2)

# want to run with multiple initial cluster assignments
set.seed(3)
km.out <- kmeans(x, 3, nstart = 1)
km.out$tot.withinss

km.out <- kmeans(x, 3, nstart = 20)
km.out$tot.withinss

# "We strongly recommend always running K-means clustering with a large
# value of nstart, such as 20 or 50""

hc.complete <- hclust(dist(x), method = "complete")
hc.average <- hclust(dist(x), method = "average")
hc.single <- hclust(dist(x), method = "single")
par(mfrow = c(1,3))
plot(hc.complete, main = "Complete Linkage", xlab = "", sub = "", cex = .9)
plot(hc.average, main = "Average Linkage", xlab = "", sub = "", cex = .9)
plot(hc.single, main = "Single Linkage", xlab = "", sub = "", cex = .9)

cutree(hc.complete, 2)
cutree(hc.average, 2)
cutree(hc.single, 2) # there's one singleton

cutree(hc.single, 4) # better, but still two singletons

# scale first
xsc <- scale(x)
plot(hclust(dist(xsc), method = "complete"), 
     main = "Hierarchical Clustering with Scaled Features")

# correlation-based distance with as.dist() function
x <- matrix(rnorm(30*3), ncol = 3)
dd <- as.dist(1-cor(t(x)))
plot(hclust(dd, method = "complete"), 
     main = "Complete Linkage with Correlation-Based Distance",
     xlab = "", sub = "")

## Lab 3 - cancer microarray data - 6830 gene expression measurements on 64 cancer cell lines

library(ISLR)
nci.labs <- NCI60$labs #better to just use the data frame!
nci.data <- NCI60$data

# first, PCA

dim(nci.data)
nci.labs[1:4]
table(nci.labs)

pr.out <- prcomp(nci.data, scale = TRUE)
Cols <- function(vec){
  cols <- rainbow(length(unique(vec)))
  return(cols[as.numeric(as.factor(vec))])
}

par(mfrow = c(1,2))
plot(pr.out$x[,1:2], col = Cols(nci.labs), pch=19, 
     xlab = "Z1", ylab = "Z2")
plot(pr.out$x[,c(1,3)], col = Cols(nci.labs), pch=19, 
     xlab = "Z1", ylab = "Z3")

summary(pr.out)

plot(pr.out)

pve <- 100*pr.out$sdev^2/sum(pr.out$sdev^2)
par(mfrow=c(1,2))

plot(pve, type = "o", 
     ylab = "PVE", xlab = "Principal Component", col = "blue")

plot(cumsum(pve), type="o", ylab = "Cumulative PVE",
     xlab = "Principal Component", col="brown3")

# next, clustering
sd.data <- scale(nci.data)
par(mfrow = c(1,3))
data.dist <- dist(sd.data)
plot(hclust(data.dist), labels=nci.labs, 
     main = "Complete Linkage", xlab = "", 
     sub = "", ylab = "")

plot(hclust(data.dist, method = "average"), labels=nci.labs, 
     main = "Average Linkage", xlab = "", 
     sub = "", ylab = "")

plot(hclust(data.dist, method = "single"), labels=nci.labs, 
     main = "Single Linkage", xlab = "", 
     sub = "", ylab = "")

hc.out = hclust(dist(sd.data))
hc.clusters = cutree(hc.out, 4)
table(hc.clusters, nci.labs)

par(mfrow = c(1,1))
plot(hc.out, labels = nci.labs)
abline(h=139, col="red")

hc.out

set.seed(2)
km.out <- kmeans(sd.data, 4, nstart = 20)
km.clusters <- km.out$cluster
table(km.clusters, hc.clusters)

# cluster using first few principal component score vectors:
hc.out <- hclust(dist(pr.out$x[,1:5]))
plot(hc.out, labels = nci.labs,
     main = "Hier. Clust. on First Five Score Vectors")
table(cutree(hc.out, 4), nci.labs)

