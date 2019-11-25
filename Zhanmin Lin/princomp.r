install.packages("factoextra")
install.packages("cluster")
install.packages("ggplot2")
library(cluster)
library(factoextra)
library(ggplot2)

mainFolder = "C:/Users/838035/Documents/ErasmusMC_datasets/CPO/Zhanmin Lin/sGluR2/"
csvList = list.files(mainFolder, recursive = T)
csvCompletePath = paste(mainFolder, csvList, sep="")

ignoreCols = 1:4
expectedClusters = 2
fileIndex = 0

resTable = data.frame(index=numeric(), filepath=character(), medicine=character(), 
                      cluster1 = numeric(), cluster2=numeric())

levels(resTable$filepath) = csvCompletePath
levels(resTable$medicine) = c("wt","A", "B", "c", "d", "e")

##########
## Do this manually and check colors
#########
fileIndex = fileIndex + 1
csvPath = csvCompletePath[fileIndex]
ds = read.csv(csvPath, header = T)
        
std_ds = data.frame(scale(ds[,-ignoreCols], center = T, scale = T))
pcaScores = princomp(std_ds)$scores[,1:4]
    
kmeansRes = kmeans(pcaScores, centers = expectedClusters, nstart = 50, iter.max = 1000)
fviz_cluster(kmeansRes, data=pcaScores, stand = FALSE, geom = "point", 
             choose.vars = c("Comp.1", "Comp.2"))

newClusterNumber = sapply(kmeansRes$cluster, function(clusterNum){
    
    if(clusterNum==1){
        2
    }else if(clusterNum==2){
        1
    }
    
}, simplify = T)

resTable[fileIndex, ] = c(fileIndex, csvPath, "wt", table(newClusterNumber)/length(newClusterNumber))

ds$cluster = newClusterNumber
write.csv(ds, file = csvPath)

####################
write.csv(resTable, file=file.choose())
