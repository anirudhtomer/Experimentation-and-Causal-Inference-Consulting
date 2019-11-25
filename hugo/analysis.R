setwd("C:/Users/838035/Google Drive/PhD/src/cpo/hugo/")
ds = read.csv(file="ds.csv", header = T, sep = ";", dec = ",")
ds.id = ds[!duplicated(ds$Patient.Nr),]

ds.id$Marker.1 = c(by(data = ds$Marker.1, INDICES = ds$Patient.Nr, mean))
ds.id$Marker.2 = c(by(data = ds$Marker.2, INDICES = ds$Patient.Nr, mean))
ds.id$Marker.3 = c(by(data = ds$Marker.3, INDICES = ds$Patient.Nr, mean))

ds.id$Marker.2 = ds.id$Marker.1 + ds.id$Marker.2
ds.id$Marker.3 = ds.id$Marker.2 + ds.id$Marker.3

ds.id$Category = factor(sapply(ds.id$Category, function(cat){
  if(cat=="HC"){
    0
  }else if(cat=="PS"){
    1
  }else{
    2
  }
}), ordered = T)

a1.prop = vglm(Category~1+Marker.1, data=ds.id[-11,], family=cumulative(parallel=T))
a1 = vglm(Category~1+Marker.1, data=ds.id[-11,], family=cumulative(parallel=F), 
          control=vglm.control(maxit=2000000))


a2.prop = vglm(Category~Marker.2, data=ds.id[-11,], family=cumulative(parallel=T))
a2 = vglm(Category~Marker.2, data=ds.id[-11,], family=cumulative(parallel=F), 
          control=vglm.control(maxit=2000000))

a3.prop = vglm(Category~Marker.3, data=ds.id[-11,], family=cumulative(parallel=T))
a3 = vglm(Category~Marker.3, data=ds.id[-11,], family=cumulative(parallel=F),
          control=vglm.control(maxit=2000000), link=)

lrtest(a1.prop,a2.prop,a3.prop)

AIC(a1);AIC(a2);AIC(a3);AIC(a1.prop);AIC(a2.prop);AIC(a3.prop)

#######
a1 = polr(Category ~ Marker.1, data = ds.id[-11,], Hess=TRUE)
a2 = polr(Category ~ Marker.2, data = ds.id[-11,], Hess=TRUE)
a3 = polr(Category ~ Marker.3, data = ds.id[-11,], Hess=TRUE)
