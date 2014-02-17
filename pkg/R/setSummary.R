setSummary <-
function(expr.data, gene.sets, sample.labels,is.log=TRUE){
  setSummary <- matrix(0,2,ncol(gene.sets))
    for(i in 1:ncol(gene.sets)){
      ind <- which(gene.sets[,i]==1)
      if(length(ind)==1){
        stop("Number of member genes just 1!")
      }
      expr <- expr.data[ind,]
      fc <- FC(expr, sample.labels)
      setSummary[1,i] <- mean(fc)
      setSummary[2,i] <- mean(abs(fc))
    }
  return(setSummary)
}
