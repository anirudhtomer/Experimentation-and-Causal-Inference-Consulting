
#First we do principal component analysis on the first data set. We have randomly checked 
#many other data sets and found that upto 5 principal components are useful to make clusters
# more than 75% variance accounted for by these 5 principal components
zhanmin = read.csv(file.choose(), header = T)[,-c(1:4)]
zhanmin_std = scale(zhanmin, center = T, scale = T)

princomp_zhanmin = princomp(zhanmin_std)
plot(princomp_zhanmin, type="l")
totalVariance = sum(princomp_zhanmin$sdev^2)
perVariance = (princomp_zhanmin$sdev^2 * 100)/totalVariance
print(round(perVariance))
round(cumsum(perVariance))

pcaScores = princomp_zhanmin$scores[,1:4]

#Now do k-means clustering. Honestly we don't expect more than 5-8 clusters
gradient_col = list(low = "steelblue", high = "white")
get_clust_tendency(pcaScores, n = 50, gradient = gradient_col)

fviz_nbclust(pcaScores, kmeans, method = "wss")
fviz_nbclust(pcaScores, pam, method = "wss") 
fviz_nbclust(pcaScores, hcut, method = "wss")    

fviz_nbclust(pcaScores, kmeans, method = "silhouette") 
fviz_nbclust(pcaScores, pam, method = "silhouette") 
fviz_nbclust(pcaScores, hcut, method = "silhouette")    

fviz_nbclust(pcaScores, kmeans, method = "gap_stat") 
fviz_nbclust(pcaScores, pam, method = "gap_stat") 
fviz_nbclust(pcaScores, hcut, method = "gap_stat")    

kmeansRes = kmeans(pcaScores, centers = 3, nstart = 50, iter.max = 1000)
fviz_cluster(kmeansRes, data=pcaScores, stand = FALSE, geom = "point", 
             choose.vars = c("Comp.2", "Comp.1"))

#For now I am choosing 2 clusters
zhanmin$cluster = kmeansRes$cluster
