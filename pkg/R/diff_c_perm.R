diff_c_perm <-
function(data, po_correct, perm.number) {

# data: expression data matrix
# po_correct: Vector of response values (example: 1,2)
# perms.number: Number of sample permutations

    
    dime2 <- dim(data)
  pit <- length(po_correct)
  
  if( !( is.array(perm.number) || is.matrix(perm.number)) ) {
    
  all_perms <- replicate(perm.number,sample(po_correct,pit,replace=FALSE))}
  unique.perm <- unique(t(all_perms))
  if(dim(unique.perm)[1]<perm.number){
        cat("Number of unique permutations:", dim(unique.perm)[1])
        all_perms <- t(unique.perm)
  }
  dime  <- dim(all_perms)
  dime2 <- dim(data)
  diff_t = array(0, c(dime2[1], dime[2] + 1))
  diff_p = diff_t
  
### positive results ###

  length(po_correct)
  dim(data)
  mySet <- new("ExpressionSet",expr=as.matrix(data))
  design <- model.matrix(~po_correct)
  fit1 <- lmFit(mySet, design)
  fit2<-eBayes(fit1)
  diff_t[,1] <- fit2$t[,2]
  diff_p[,1] <- fit2$p.value[,2]
  
## permutation results ###

  for(k in 1:dime[2]){
    design <- model.matrix(~all_perms[,k])
    fit1 <- lmFit(mySet, design)
    fit2<-eBayes(fit1)
    diff_t[,k+1] <- fit2$t[,2]
    diff_p[,k+1] <- fit2$p.value[,2]
    }
  out <- list(p_values = diff_p, t_scores = diff_t, perms = unique.perm)
  }


################

diff.FC.perm <-
function(data, sample.labels,perm.number,is.log=TRUE){
	data <- as.matrix(data)
	pit <- length(sample.labels)
	all_perms <- replicate(perm.number,sample(c(1:pit),pit,replace=FALSE))
    unique.perm <- unique(t(all_perms))
    if(dim(unique.perm)[1]<perm.number){
        cat("Number of unique permutations:", dim(unique.perm)[1])
        all_perms <- t(unique.perm)
    }
	dime  <- dim(all_perms)
  	dime2 <- dim(data)
  	diff_FC <- array(0,c(dime2[1],dime[2]+1))
  	
  	#################################
  	## FC for positive data #########
  	#################################
  	

  	diff_FC[,1] <- FC(data,sample.labels)
  	
  	#################################
  	## FC for col permuted data #####
  	#################################
  	
  	for (i in 1:dime[2]){
  		
  		diff_FC[,i+1] <- FC(data[,all_perms[,i]],sample.labels)
  	}
  return(diff_FC)
}

###########

FC <-
function(expr.data,sample.labels,is.log=TRUE){
    # Sanity check
    if(length(unique(sample.labels))==1) {
        stop("No replicates to average!")
    }
    if(length(unique(sample.labels))>2) {
        stop("Fold change cannot be calculated for more than two groups!")}
    
    #Calculating averages
    dat<-matrix(nrow=nrow(expr.data), ncol=length(unique(sample.labels)), NA)
    
    index=1
    for( i in unique(sample.labels) ) {
        dat[,index]<-rowSums(data.frame(expr.data[,which(sample.labels==i)]))/ncol(data.frame(expr.data[,which(sample.labels==i)]))
        index=index+1
    }
    
    # Calculating the fold change

    FC <- dat[,1] / dat[,2]
        
    
    if (is.log==TRUE) {
    FC <- log2(FC)
    }

return(FC)
}