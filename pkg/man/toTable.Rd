\name{toTable}
\alias{toTable}

\title{Table with top gene sets}
\description{Table with top gene sets}
\usage{
toTable(mGSZobj, sample.perm.table = FALSE, method = c("mGSZ", "mGSA", "mAllez", "WRS", "SS", "SUM", "KS", "wKS"), no.top.sets = 10)
}

\arguments{
  \item{mGSZobj}{
mGSZ object
}

  \item{sample.perm.table}{
TRUE for table of top gene sets based on sample permutation when both gene and sample permutations were used.
}

  \item{method}{
Method for which table for top gene sets is required (Required only when other methods were used for the gene set analysis)
}
  \item{no.top.sets}{
Number of top gene sets in the table
}
}

\value{
A table with top gene sets
}
\references{
Mishra Pashupati, Toronen Petri, Leino Yrjo, Holm Liisa. Gene Set Analysis: Limitations in popular existing methods and proposed improvements (Not yet published)
http://ekhidna.biocenter.helsinki.fi/downloads/pashupati/mGSZ.html
}
\author{Pashupati Mishra, Petri Toronen}

\seealso{\code{\link{mGSZ}}}
\examples{
gene.names <- paste("g",1:1000, sep = "")

# create random gene expression data matrix

set.seed(100)
expr.data <- matrix(rnorm(1000*50),ncol=50)
rownames(expr.data) <- gene.names
b <- matrix(2*rnorm(2500),ncol=25)
ind <- sample(1:100,replace=FALSE)
expr.data[ind,26:50] <- expr.data[ind,26:50] + b

sample.labels <- rep(1:2,c(25,25))

# create random gene sets

gene.sets <- vector("list", 100)
for(i in 1:length(gene.sets)){
	gene.sets[[i]] <- sample(gene.names, size = 20)
}
names(gene.sets) <- paste("set", as.character(1:100), sep="")

mGSZ.obj <- mGSZ(expr.data, gene.sets, sample.labels, perm.number = 100)
top.mGSZ.sets <- toTable(mGSZ.obj, no.top.sets = 10) 
}
