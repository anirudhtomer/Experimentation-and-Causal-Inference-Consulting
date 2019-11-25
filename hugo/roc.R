#0 vs other

plotdf = data.frame(Marker2levels = Marker2levels, tpr = NA, fpr = NA)
Marker2levels = seq(min(ds.id$Marker.2), max(ds.id$Marker.2), by = 2)

plotdf$Marker2levels = Marker2levels

ds.id$hcVsOth = ifelse(ds.id$Category=="0", yes=0, 1)

i = 1

bestYouden = -Inf
bestCutoff = -Inf
bestTpr = NA
bestFpr = NA
bestFnr = NA
bestTnr = NA
for(marker2level in Marker2levels){
    ds.id$tempCat = ifelse(ds.id$Marker.2<marker2level, yes=1, 0)
    
    tp = sum(ds.id$tempCat ==0 & ds.id$hcVsOth == 0)
    fp = sum(ds.id$tempCat ==0 & ds.id$hcVsOth == 1)
    fn = sum(ds.id$tempCat ==1 & ds.id$hcVsOth == 0)
    tn = sum(ds.id$tempCat ==1 & ds.id$hcVsOth == 1)
    
    plotdf$tpr[i] = tp / (tp + fn)
    plotdf$fpr[i] = fp / (fp + tn)
    
    J = tp/(tp+fn)   + tn/(tn+fp) - 1
    #J = 2*tp/(fp + fn + 2*tp)
    if(J>bestYouden){
      bestYouden = J
      bestCutoff = marker2level
      bestTpr = plotdf$tpr[i]
      bestFpr = plotdf$fpr[i]
      bestFnr = fn / (fn + tp)
      bestTnr = tn / (tn + fp)
    }
    
    i = i+1
}

print(bestFnr)
print(bestFpr)
print(bestTnr)
print(bestTpr)

ggplot(data=plotdf) + 
  geom_line(aes(x=fpr, y=tpr)) + xlab("1 - Specificity") + ylab("Sensitivity") + 
  ggtitle("Healthy Controls vs. everyone else")

#2 vs other
Marker2levels = seq(min(ds.id$Marker.2), max(ds.id$Marker.2), by = 2)

ds.id$piVsOth = ifelse(ds.id$Category=="2", yes=0, 1)


i = 1

bestYouden = -Inf
bestCutoff = -Inf
bestTpr = NA
bestFpr = NA
bestFnr = NA
bestTnr = NA
for(marker2level in Marker2levels){
  ds.id$tempCat = ifelse(ds.id$Marker.2<marker2level, yes=0, 1)
  
  tp = sum(ds.id$tempCat ==0 & ds.id$piVsOth == 0)
  fp = sum(ds.id$tempCat ==0 & ds.id$piVsOth == 1)
  fn = sum(ds.id$tempCat ==1 & ds.id$piVsOth == 0)
  tn = sum(ds.id$tempCat ==1 & ds.id$piVsOth == 1)
  
  plotdf$tpr[i] = tp / (tp + fn)
  plotdf$fpr[i] = fp / (fp + tn)
  
  J = tp/(tp+fn)   + tn/(tn+fp) - 1
  #J = 2*tp/(fp + fn + 2*tp)
  if(J>bestYouden){
    bestYouden = J
    bestCutoff = marker2level
    bestTpr = plotdf$tpr[i]
    bestFpr = plotdf$fpr[i]
    bestFnr = fn / (fn + tp)
    bestTnr = tn / (tn + fp)
  }
  i = i+1
}

print(bestFnr)
print(bestFpr)
print(bestTnr)
print(bestTpr)


ggplot(data=ds.id) + geom_point(aes(x=Marker.2, y=y, color=Category)) +
  geom_vline(xintercept = 25.15) + geom_vline(xintercept = 97.15)
