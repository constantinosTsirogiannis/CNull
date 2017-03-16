\name{permutation.random.values.b}
\alias{permutation.random.values.b}
\title{Produces a set of random values for a beta diversity measure f on a given matrix, based on the permutation (SIM2) null model}
\description{
Given a matrix M, a beta diversity measure f and a number of repetitions k, the function produces k random values of f based on the permutation model. This is equivalent to shuffling M according to this model as many as k times , each time outputing the value of f only for a certain pair of rows (e.g. the two top ones) of the shuffled matrix. The output values can be used to determine the null distribution of f for a pair of rows in M. This distribution is the same for every pairs of rows in M. This is because the examined null model produces the same distribution for all row pairs of M; after shuffling M, each row has the same probability to store a given community as any other in the resulting matrix. 
}
\usage{
permutation.random.values.b(matrix,f,args,reps=1000)
}
\arguments{
  \item{matrix}{A matrix with integer values. The matrix should not contain any NA values.}
  \item{f}{A beta diversity function f. The interface of f should be such that f(matrix,args) returns a matrix V where the entry stored at the i-th row and j-th column of V is equal to the value of f when applied at the i-th and j-th row of the input matrix. To fit this interface, the user might have to develop f as a wrapper around an existing R function (see \bold{Examples}).}
  \item{args}{A list with extra arguments needed by f.}
  \item{reps}{The number of randomizations. This argument is optional and its default value is set to one thousand.}
}
\value{A vector of as many as reps elements. Stores the randomized values of f calculated based on the permutation null model.}
\references{
Gotelli, N. J., 2000. Null Model Analysis of Species Co-Occurrence Patterns. Ecology, 81(9), pp.2606-2621.

Tsirogiannis, C., A. Kalvisa, B. Sandel and T. Conradi. Column-Shuffling Null Models Are Simpler Than You Thought. To appear.
}
\author{Constantinos Tsirogiannis (tsirogiannis.c@gmail.com)}
\seealso{\code{\link{permutation.moments.b}}}
\examples{
#In the next example null-model calculations are
#performed using a function of phylogenetic diversity.
#Hence, we first load the required packages.
require(CNull)
require(ape)
require(PhyloMeasures)

#Load phylogenetic tree of bird families from package "ape"
data(bird.families, package = "ape")

#Create 100 random communities with 50 families each
comm = matrix(0,nrow = 100,ncol = length(bird.families$tip.label))
for(i in 1:nrow(comm)) {comm[i,sample(1:ncol(comm),50)] = 1}
colnames(comm) = bird.families$tip.label

#Set function f to be the Common Branch Length measure (CBL)
#as defined in the R package PhyloMeasures.
my.f = function(mt,args){ return (cbl.query(args[[1]],mt))}

# This function takes one extra argument, which is a phylogenetic tree.
# Hence, create a list whose only element is the desired tree.
arguments = list()
arguments[[1]] = bird.families

# Calculate 2000 randomized values of f on comm
# based on the permutation null model.
permutation.random.values.b(comm,f=my.f,args=arguments,reps=2000)
}

