toTable <-
function(mGSZobj, sample.perm.table = FALSE, method = c("mGSZ", "mGSA", "mAllez", "WRS", "SS", "SUM", "KS", "wKS") ,no.top.sets = 10){

  if(!mGSZobj$gene.perm.log & !sample.perm.table & !mGSZobj$other.methods){
    table <- mGSZobj$mGSZ[1:no.top.sets,]
  }
  
  else if(mGSZobj$gene.perm.log & !sample.perm.table & !mGSZobj$other.methods){
    table <- mGSZobj$mGSZ.gene.perm[1:no.top.sets,]
  }
  
  else if(mGSZobj$gene.perm.log & sample.perm.table & !mGSZobj$other.methods){
    table <- mGSZobj$mGSZ.sample.perm[1:no.top.sets,]
  }
  
  else if(!mGSZobj$gene.perm.log & !sample.perm.table & mGSZobj$other.methods){
    method <- match.arg(method)
    table <- mGSZobj[[method]][1:no.top.sets,]
  }
  
  else if(mGSZobj$gene.perm.log & !sample.perm.table & mGSZobj$other.methods){
    method <- match.arg(method)
    table <- mGSZobj$gene.perm[[method]][1:no.top.sets,]
  }
  
  else if(mGSZobj$gene.perm.log & sample.perm.table & mGSZobj$other.methods){
    method <- match.arg(method)
    table <- mGSZobj$sample.perm[[method]][1:no.top.sets,]
  }
  return(table)
}
