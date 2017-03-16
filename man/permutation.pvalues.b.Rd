\name{permutation.pvalues.b}
\alias{permutation.pvalues.b}
\title{Produces the p-values of a beta diversity measure f on a given matrix, based on the permutation (SIM2) null model}
\description{
Given a matrix M and a beta diversity measure f, the function calculates the p-values of f based on the permutation null model.}
\usage{
permutation.pvalues.b(matrix, f, args, observed.vals, reps=1000)
}
\arguments{
  \item{matrix}{A matrix with integer values. The matrix should not contain any NA values.}
  \item{f}{A beta diversity function f. The interface of f should be such that f(matrix,args) returns a matrix V where the entry stored at the i-th row and j-th column of V is equal to the value of f when applied at the i-th and j-th row of the input matrix. To fit this interface, the user might have to develop f as a wrapper around an existing R function (see \bold{Examples}).}
  \item{args}{A list with extra arguments needed by f.}
  \item{observed.vals}{A set of pre-calculated values for which we want to compute their p-values.}
  \item{reps}{The number of randomizations. This argument is optional and its default value is set to one thousand.}
}
\value{A numeric vector that stores the p-values of f, calculated based on the permutation null model. The i-th element of the vector stores the p-value for the i-th element in vector observed.vals.}
\references{
Gotelli, N. J., 2000. Null Model Analysis of Species Co-Occurrence Patterns. Ecology, 81(9), pp.2606-2621.

Tsirogiannis, C., A. Kalvisa, B. Sandel and T. Conradi. Column-Shuffling Null Models Are Simpler Than You Thought. To appear.
}
\author{Constantinos Tsirogiannis (tsirogiannis.c@gmail.com)}
\seealso{\code{\link{permutation.communities.b}}}
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

#Compute the values of f for all pairs of observed communities in M.
#Turn the resulting matrix with the observed diversity values into a vector  
obs.v=my.f(comm,arguments)
vals = as.vector(t(obs.v))

# Calculate the p-values of f for the communities in comm 
# based on the permutation model, using 2000 Monte Carlo randomizations.
permutation.pvalues.b(comm,f=my.f,args=arguments,observed.vals=vals,reps=2000)
}

