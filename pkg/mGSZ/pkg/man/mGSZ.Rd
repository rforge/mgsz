\name{mGSZ}
\alias{mGSZ}

\title{Gene set analysis based on Gene Set Z-scoring function and asymptotic p-value}
\description{
Gene set analysis based on Gene Set Z scoring function and asymptotic p-value
}
\usage{
mGSZ(expr.data, gene.sets, sample.labels, flip.gene.sets = FALSE, select= 'T-score', is.log = TRUE, gene.perm = FALSE,  min.cl.sz=5, other.methods = FALSE, pre.var = 0, wgt1 = 0.2, wgt2 = 0.5, var.constant = 10, start.val=5, perm.number = 200)
}

\arguments{
  \item{expr.data}{
Gene expression data matrix (rows as genes and columns as samples)
}
  \item{gene.sets}{
Gene set data (dataframe/table/matrix/list)
}
  \item{sample.labels}{
Vector of response values (example:1,2)
}
  \item{flip.gene.sets}{
TRUE if gene set data is list with genes as list names
}
  \item{select}{
Gene level statistics (example: T-score/FC/P-value)
}
  \item{is.log}{
TRUE for log fold change as gene level statistics
}
  \item{gene.perm}{
TRUE for analysis with both gene and sample permutation data as the null distributions
}
  \item{min.cl.sz}{
Minimum size of gene sets (number of genes in a gene set) to be included in the analysis
}
  \item{other.methods}{
TRUE for gene set analysis with other methods (see the manuscript for details)
}
  \item{pre.var}{
Estimate of the variance associated with each observation
}
  \item{wgt1}{
Weight 1, parameter used to calculate the prior variance obtained with class size var.constant. This penalizes especially small classes and small subsets. Default is 0.2. Values around 0.1 - 0.5 are expected to be reasonable.
}
  \item{wgt2}{
Weight 2, parameter used to calculate the prior variance obtained with the same class size as that of the analyzed class. This penalizes small subsets from the gene list. Default is 0.5. Values around 0.3 and 0.5 are expected to be reasonable
}
  \item{var.constant}{
Size of the reference class used with wgt1. Default is 10
}
  \item{start.val}{
Number of first genes to be left out from the analysis
}
  \item{perm.number}{
Number of permutations for p-value calculation
}
}
\details{
A function for Gene set analysis based on Gene Set Z-scoring function and asymptotic p- value. It differs from GSZ (Toronen et al 2009) in that it implements asymptotic p-values instead of empirical p-values. Asymptotic p-values are based on fitting suitable distribution model to the permutation data. Unlike empirical p-values, the resolution of asymptotic p-values are independent of the number of permutations and hence requires consideralbly fewer permutations. In addition to GSZ, this function allows the users to carry out analysis with seven other scoring functions (visit http://ekhidna.biocenter.helsinki.fi/downloads/pashupati/mGSZ.html for a more detailed description) and compare the results.
}
\value{

 \item{mGSZ}{Dataframe with gene sets (in decreasing order based on the significance) reported by mGSZ method and their sizes, scores, p-values and gene set expression summary}
  \item{mGSA}{Dataframe with gene sets (in decreasing order based on the significance) reported by mGSA method and their sizes, scores, p-values and gene set expression summary}
 \item{mAllez}{Dataframe with gene sets (in decreasing order based on the significance) reported by mAllez method and their sizes, scores, p-values and gene set expression summary}
 \item{WRS}{Dataframe with gene sets (in decreasing order based on the significance) reported by WRS method and their sizes, scores, p-values and gene set expression summary}
 \item{SUM}{Dataframe with gene sets (in decreasing order based on the significance) reported by SUM method and their sizes, scores, p-values and gene set expression summary}
 \item{SS}{Dataframe with gene sets (in decreasing order based on the significance) reported by SS method and their sizes, scores, p-values and gene set expression summary}
 \item{KS}{Dataframe with gene sets (in decreasing order based on the significance) reported by KS method and their sizes, scores, p-values and gene set expression summary}
 \item{wKS}{Dataframe with gene sets (in decreasing order based on the significance) reported by wKS method and their sizes, scores, p-values and gene set expression summary}
 \item{sample.labels}{Vector of response values used}
 \item{perm.number}{Number of permutations used for p-value calculation}
 \item{expr.data}{For internal use}
 \item{gene.sets}{For internal use}
 \item{flip.gene.sets}{For internal use}
 \item{min.cl.sz}{For internal use}
 \item{other.methods}{For internal use}
 \item{pre.var}{For internal use}
 \item{wgt1}{For internal use}
 \item{wgt2}{For internal use}
 \item{var.constant}{For internal use}
 \item{start.val}{For internal use}
 \item{select}{For internal use}
 \item{is.log}{For internal use}
 \item{gene.perm.log}{For internal use}
}
\references{
Mishra Pashupati, Toronen Petri, Leino Yrjo, Holm Liisa. Gene Set Analysis: Limitations in popular existing methods and proposed improvements (Not yet published)
http://ekhidna.biocenter.helsinki.fi/downloads/pashupati/mGSZ.html

Toronen, P., Ojala, P. J., Marttinen, P., and Holm, L. (2009). Robust extraction of functional signals from gene set analysis using a generalized threshold free scoring function. BMC Bioinformatics, 10(1), 307.

}
\author{Pashupati Mishra, Petri Toronen}


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

# scoring function profile data across the ordered gene list for top 5 gene sets

data4plot <- StabPlotData(mGSZ.obj,rank.vector=c(1,2,3,4,5))

# profile plot for the top gene set
plotProfile(data4plot,1)  

# gene sets in a gmt format can be converted to mGSZ readable format as follows:
# gene.sets <- geneSetsList("gene.sets.gmt")

}
