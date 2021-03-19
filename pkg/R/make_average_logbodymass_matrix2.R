make_average_logbodymass_matrix2<-function(a){
 n=dim(a)
 res=matrix(0,(n[1]-1),n[2])
 for (i in 1:(n[2])){
  res[,i]=make_average_logbodymass(a[,i])
 }
res
}
