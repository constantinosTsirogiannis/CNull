\name{permutation.communities.b}
\alias{permutation.communities.b}
\title{Produces random pairs of communities from a given matrix, based on the permutation (SIM2) null model. These communities can be used for beta diversity computations}
\description{
Given a matrix M and a number of repetitions k, the function produces k pairs of random communities based on the permutation null model. This is equivalent to shuffling M according to this model as many as k times, each time outputing only two certain rows (e.g. the two top ones) of the shuffled matrix. A beta diversity measure f can be applied on the output communities to determine the null distribution of f for two rows in M. This distribution is the same for every pairs of rows in M. This is because the examined null model produces the same distribution for all pairs of rows in M; after shuffling M, each pair of rows has the same probability to store two specific communities as any other pair in the resulting matrix. 
}
\usage{
permutation.communities.b(matrix, reps=1000)
}
\arguments{
  \item{matrix}{A matrix with integer values. The matrix should not contain any NA values.}
  \item{reps}{The number of randomizations. This argument is optional and its default value is set to one thousand.}
}
\value{A matrix which stores pairs of communities generated based on the permutation model. The output matrix has two times reps rows, corresponding to as many as reps pairs of randomized communities. For every i in 1:reps, rows in the output matrix with indices 2i-1 and 2i designate a valid random pair of communities. The number of columns, and the names of the columns in the output matrix are the same as in the input matrix. The output matrix is \emph{not} a shuffled version of the input matrix.}
\references{
Gotelli, N. J., 2000. Null Model Analysis of Species Co-Occurrence Patterns. Ecology, 81(9), pp.2606-2621.

Tsirogiannis, C., A. Kalvisa, B. Sandel and T. Conradi. Column-Shuffling Null Models Are Simpler Than You Thought. To appear.
}
\author{Constantinos Tsirogiannis (tsirogiannis.c@gmail.com)}
\seealso{\code{\link{permutation.random.values.b}}}
\examples{
require(CNull)

#Create a random integer matrix
comm = matrix(sample(1:300),nrow=15,ncol=20)

#Use permutation model to produce 2000 random pairs of communities 
permutation.communities.b(comm,reps=2000)
}

