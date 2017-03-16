\name{individual.based.communities.a}
\alias{individual.based.communities.a}
\title{Produces a set of random communities from a given matrix, based on the individual-based null model (Stegen et al. 2013). These communities can be used for alpha diversity computations}
\description{
Given a matrix M and a number of repetitions k, the function produces k random communities based on the individual-based model. This is equivalent to shuffling M according to this model as many as k times , each time outputing only a certain row (e.g. the top one) of the shuffled matrix. An alpha diversity measure f can be applied on the output communities to determine the null distribution of f for a row in M. This distribution is the same for every row of M. This is because the examined null model produces the same distribution for all rows of M; after shuffling M, each row has the same probability to store a specific community as any other in the resulting matrix. 
}
\usage{
individual.based.communities.a(matrix, reps=1000)
}
\arguments{
  \item{matrix}{A matrix with integer values. The matrix should not contain any NA values.}
  \item{reps}{The number of randomizations. This argument is optional and its default value is set to one thousand.}
}
\value{A matrix which stores communities generated based on the individual-based model. Each row of the matrix corresponds to a different randomized community. The number of columns, and the names of the columns in the output matrix are the same as in the input matrix. The output matrix is \emph{not} a shuffled version of the input matrix.}
\references{
Stegen, J. C., Freestone, A. L., Crist, T. O., Anderson, M. J., Chase, J. M., Comita, L. S., Cornell, H. V., Davies, K. F., Harrison, S. P., Hurlbert, A. H., Inouye, B. D., Kraft, N. J. B., Myers, J. A., Sanders, N. J., Swenson, N. G., Vellend, M. (2013), Stochastic and Deterministic Drivers of Spatial and Temporal Turnover in Breeding Bird Communities. Global Ecology and Biogeography, 22:202-212.

Tsirogiannis, C., A. Kalvisa, B. Sandel and T. Conradi. Column-Shuffling Null Models Are Simpler Than You Thought. To appear.
}
\author{Constantinos Tsirogiannis (tsirogiannis.c@gmail.com)}
\seealso{\code{\link{individual.based.random.values.a}}}
\examples{
require(CNull)

#Create a random integer matrix
comm = matrix(sample(1:300),nrow=15,ncol=20)

#Use individual-based model to produce 2000 random communities 
individual.based.communities.a(comm,reps=2000)
}

