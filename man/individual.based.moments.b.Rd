\name{individual.based.moments.b}
\alias{individual.based.moments.b}
\title{Given a matrix and a beta diversity measure f, calculates the mean and variance of f based on the individual-based null model (Stegen et al. 2013)}
\description{
Given a matrix M and a beta diversity measure f, the function computes the mean and variance of f between a pair of rows in M, when M is shuffled according to the individual-based model. The returned mean and variance is the same for every pair of rows in M. This is because the examined null model produces the same distribution for all pairs of rows in M; after shuffling M, each pair of rows has the same probability to store two specific communities as any other pair in the resulting matrix. 
}
\usage{
individual.based.moments.b(matrix,f,args,reps=1000)
}
\arguments{
  \item{matrix}{A matrix with integer values. The matrix should not contain any NA values.}
  \item{f}{A beta diversity function f. The interface of f should be such that f(matrix,args) returns a matrix V where the entry stored at the i-th row and j-th column of V is equal to the value of f when applied at the i-th and j-th row of the input matrix. To fit this interface, the user might have to develop f as a wrapper around an existing R function (see \bold{Examples}).}
  \item{args}{A list with extra arguments needed by f.}
  \item{reps}{The number of randomizations. This argument is optional and its default value is set to one thousand.}
}
\value{A list with two real numbers; the mean and the variance of f on the given matrix for the individual-based model.}
\references{
Stegen, J. C., Freestone, A. L., Crist, T. O., Anderson, M. J., Chase, J. M., Comita, L. S., Cornell, H. V., Davies, K. F., Harrison, S. P., Hurlbert, A. H., Inouye, B. D., Kraft, N. J. B., Myers, J. A., Sanders, N. J., Swenson, N. G., Vellend, M. (2013), Stochastic and Deterministic Drivers of Spatial and Temporal Turnover in Breeding Bird Communities. Global Ecology and Biogeography, 22: 202-212.

Tsirogiannis, C., A. Kalvisa, B. Sandel and T. Conradi. Column-Shuffling Null Models Are Simpler Than You Thought. To appear.
}
\author{Constantinos Tsirogiannis (tsirogiannis.c@gmail.com)}
\seealso{\code{\link{individual.based.pvalues.b}}}
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

# Calculate the mean and variance of f in the individual-based null
# model using 2000 Monte Carlo randomizations
individual.based.moments.b(comm,f=my.f,args=arguments,reps=2000)
}

