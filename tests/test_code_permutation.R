########################################################################################
##  Copyright (C) 2017,  Constantinos Tsirogiannis.  Email: tsirogiannis.c@gmail.com  ##
##                                                                                    ##
##  This file is part of CNull.                                                       ##
##                                                                                    ##
##  CNull is free software: you can redistribute it and/or modify                     ##
##  it under the terms of the GNU General Public License as published by              ##
##  the Free Software Foundation, either version 3 of the License, or                 ##
##  (at your option) any later version.                                               ##
##                                                                                    ##
##  CNull is distributed in the hope that it will be useful,                          ##
##  but WITHOUT ANY WARRANTY; without even the implied warranty of                    ##
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                     ##
##  GNU General Public License for more details.                                      ##
##                                                                                    ##
##  You should have received a copy of the GNU General Public License                 ##
##  along with CNull.  If not, see <http://www.gnu.org/licenses/>                     ##
########################################################################################

require(Rcpp)

source("standard_shuffling_methods.R")

#############################################################
#############################################################
################### Auxilliary functions ####################
#############################################################
#############################################################

myvar = function(vals)
{ return (mean(vals*vals) - (mean(vals)*mean(vals))) }

#############################################################
#############################################################
################# Alpha diversity functions #################
#############################################################
#############################################################

permutation.communities.a.tester = function(matrix,error, reps)
{
  samples = permutation.communities.a(matrix, reps)

  # Check that the dimensions of the output and its column names are correct.

  if(ncol(samples) != ncol(matrix) )
  {
    cat("
 The samples-matrix has a different number of columns than the input matrix.")
    cat("
 The number of columns of the input matrix is: " , ncol(matrix))
    cat("
 The number of columns of the sample-matrix is: " , ncol(samples) , "
")
    return(FALSE)
  }

  if(nrow(samples) != reps )
  {
    cat("
 The samples-matrix has a different number of rows than the requested repetitions.")
    cat("
 The number of requested repetitions is: " , reps)
    cat("
 The number of rows of the sample-matrix is: " , nrow(samples) , "
")
    return(FALSE)
  }

  cnm.input = colnames(matrix)
  cnm.samples = colnames(samples)

  for(i in 1:ncol(samples))
    if(cnm.input[i] != cnm.samples[i])
    {
      cat("
 The samples-matrix has at least one different column-name than the input matrix.")
      cat("
 The difference was observed at column #" , i)
      cat("
 The column name of the input matrix for this column is: " , cnm.input[i])
      cat("
 The column name of the samples matrix for this column is: " , cnm.samples[i] , "
")
      return(FALSE)
    }


  # Check if all output values are valid

  for( j in 1:ncol(samples) )
  {
    ubik = unique(matrix[,j])

    for( i in 1:nrow(samples) )
    {
      if( samples[i,j] %in% ubik == FALSE )
      {
        cat("
 There is at least one sample that has an element which is out of bounds.")
        cat("
 The discrepancy was found in element: ", i, " , ", j, " of the samples matrix.")
        cat("
 The value of samples matrix for this element is: " , samples[i,j])
        cat("
 For this column, the elements of the matrix are: " , matrix[,j] , "
")
        return(FALSE)
      }

    } #for( j in 1:nrow(samples) )

  } #for( j in 1:ncol(samples) )

  # Check if mean and variance of values in every column are close to the ideal ones

  for( j in 1:ncol(samples) )
  {
    mean.new = mean(samples[,j])
    var.new = myvar(samples[,j])

    mean.actual = mean(matrix[,j])
    var.actual = myvar(matrix[,j])

    abs.error.mean = abs(mean.actual-mean.new)
    
    rel.error.mean = 1.0
    
    if(mean.actual > 0.05)
      rel.error.mean = abs(mean.actual-mean.new)/abs(mean.actual)

    abs.error.sd = abs(sqrt(var.actual)-sqrt(var.new))

    rel.error.sd = 1.0

    if(var.actual > 0.05)
      rel.error.sd = abs(sqrt(var.actual)-sqrt(var.new))/abs(sqrt(var.actual)) 

    if(abs.error.mean > error && rel.error.mean > error  )
    {
      cat("
 The column mean computed by the alpha-diversity permutation communities function",
            "
 diverges a lot from the actual mean.")
      cat("
 The column where this was observed has index: ", j )
      cat("
 The elements of these column are: ", matrix[,j] )
      cat("
 The value returned by the new function is: " , mean.new)
      cat("
 The actual mean is: " , mean.actual)
      cat("
 The relative difference is (percent): " , (rel.error.mean*100) , "
")
      return(FALSE)
    }

    if(abs.error.sd > error && rel.error.sd > error  )
    {
      cat("
 The column st.deviation computed by the alpha-diversity permutation",
            "
 communities function diverges a lot from the actual deviation.")
      cat("
 The column where this was observed has index: ", j )
      cat("
 The value computed by the new function is: " , sqrt(var.new))
      cat("
 The actual deviation is: " , sqrt(var.actual))
      cat("
 The relative difference is (percent): " , (rel.error.sd*100) , "
")
      return(FALSE)
    } 

  } # for( j in 1:ncol(samples) )

  return (TRUE)

} # permutation.communities.a.tester = function(...)

permutation.random.values.a.tester = function(matrix, f, args, error=0.1, reps=1000)
{
  #res.standard = moments.permutation.shuffling.a.standard(matrix,f,args,reps)
  res.standard = moments.permutation.shuffling.a.single.row(matrix,f,args,reps)
  res.new = permutation.random.values.a(matrix,f,args,reps)
 
  if(length(res.new) != reps )
  {
    cat("
 The alpha-diversity permutation moments tester returns a wrong number of elements.")
    cat("
 The returned number of elements is: " , length(res.new) , "
")
    return(FALSE)
  }

  mean.standard = res.standard[[1]]
  var.standard = res.standard[[2]]

  mean.new = mean(res.new)
  var.new = var(res.new)

  error.mean = 0.0
  error.sd = 0.0

  if(mean.standard < 0.05)
    error.mean = abs(mean.standard-mean.new)
  else
    error.mean = abs(mean.standard-mean.new)/abs(mean.standard)

  if(var.standard < 0.05)
    error.sd = abs(sqrt(var.standard)-sqrt(var.new))
  else
    error.sd = abs(sqrt(var.standard)-sqrt(var.new))/abs(sqrt(var.standard))

  if(error.mean > error )
  {
    cat("
 The mean computed by the alpha-diversity permutation moments function diverges",
          "
 a lot from the value computed by the standard method.")
    cat("
 The value returned by the new function is: " , mean.new)
    cat("
 The value returned by the standard function is: " , mean.standard)
    cat("
 The relative difference is (percent): " , (error.mean*100) , "
")
    return(FALSE)
  }

  if(error.sd > error )
  {
    cat("
 The st.deviation computed by the alpha-diversity permutation moments function diverges",
          "
 a lot from the value computed by the standard method.")
    cat("
 The value computed by the new function is: " , sqrt(var.new))
    cat("
 The value computed by the standard function is: " , sqrt(var.standard))
    cat("
 The relative difference is (percent): " , (error.sd*100) , "
")
    return(FALSE)
  }

  return (TRUE)

} # permutation.random.values.a.tester = function(...)

permutation.moments.a.tester = function(matrix, f, args, error=0.1, reps=1000)
{
  #res.standard = moments.permutation.shuffling.a.standard(matrix,f,args,reps)
  res.standard = moments.permutation.shuffling.a.single.row(matrix,f,args,reps)
  res.new = permutation.moments.a(matrix,f,args,reps)
 
  if(length(res.new) != 2 )
  {
    cat("
 The alpha-diversity permutation moments tester returns a wrong number of elements.")
    cat("
 The returned number of elements is: " , length(res.new) , "
")
    return(FALSE)
  }

  mean.standard = res.standard[[1]]
  var.standard = res.standard[[2]]

  mean.new = res.new[[1]]
  var.new = res.new[[2]]

  error.mean = 0.0
  error.sd = 0.0

  if(mean.standard < 0.05)
    error.mean = abs(mean.standard-mean.new)
  else
    error.mean = abs(mean.standard-mean.new)/abs(mean.standard)

  if(var.standard < 0.05)
    error.sd = abs(sqrt(var.standard)-sqrt(var.new))
  else
    error.sd = abs(sqrt(var.standard)-sqrt(var.new))/abs(sqrt(var.standard))

  if(error.mean > error )
  {
    cat("
 The mean computed by the alpha-diversity permutation moments function diverges",
          "
 a lot from the value computed by the standard method.")
    cat("
 The value returned by the new function is: " , mean.new)
    cat("
 The value returned by the standard function is: " , mean.standard)
    cat("
 The relative difference is (percent): " , (error.mean*100) , "
")
    return(FALSE)
  }

  if(error.sd > error )
  {
    cat("
 The st.deviation computed by the alpha-diversity permutation moments function diverges",
          "
 a lot from the value computed by the standard method.")
    cat("
 The value computed by the new function is: " , sqrt(var.new))
    cat("
 The value computed by the standard function is: " , sqrt(var.standard))
    cat("
 The relative difference is (percent): " , (error.sd*100) , "
")
    return(FALSE)
  }

  return (TRUE)

} # permutation.moments.a.tester = function(...)


permutation.pvalues.a.tester = function(matrix, f, args, error= 0.01, reps=1000)
{
  pvalues.standard = pvalues.permutation.shuffling.a.standard(matrix,f,args,reps)
  #pvalues.standard = pvalues.permutation.shuffling.a.single.row(matrix,f,args,reps)
  pvalues.new = permutation.pvalues.a(matrix,f,args,reps)

  if(length(pvalues.new) != nrow(matrix) )
  {
    cat("
 There was an error in the output of the alpha-diversity permutation p-value function.")
    cat("
 The number of returned p-values is different from the number of rows in the input matrix.")
    cat("
 The number of rows in the input matrix is: " , nrow(matrix))
    cat("
 The number of returned p-values is: " , length(pvalues.new) , "
")
    return(FALSE)
  }

  for(i in 1:length(pvalues.new))
    if( pvalues.new[i] <= 0.0 || pvalues.new[i] > 1.0 )
    {
      cat("
 There was an error in the output of the alpha-diversity permutation p-value function.")
      cat("
 At least of the returned p-values is out of bounds.")
      cat("
 This is the value with index: " , i)
      cat("
 This p-value is equal to: " , pvalues.new[i] , "
")
      return(FALSE)
    }

  error.pvalues = vector(mode="numeric",length=nrow(matrix))
  
  for( i in 1:nrow(matrix) )
    error.pvalues[i] = abs(pvalues.standard[i]-pvalues.new[i])

  for(i in 1:length(pvalues.new))
    if( error.pvalues[i] > error )
    {
      cat("
 There was an error in the output of the alpha-diversity permutation p-value function.")
      cat("
 At least one of the p-values computed by the function diverges",
          "
 a lot from the value computed by the standard method.")
      cat("
 This is the value with index: " , i)
      cat("
 The value returned by the new function is: " , pvalues.new[i])
      cat("
 The value returned by the standard function is: " , pvalues.standard[i])
      cat("
 The absolute error is: " , error.pvalues[i], "
")
      return(FALSE)
    }

  return (TRUE)
  
} # permutation.pvalues.a.tester=function(...)

#####################
#####################
## Beta diversity  ##
#####################
#####################

permutation.communities.b.tester = function(matrix,error, reps)
{
  samples = permutation.communities.b(matrix, reps)

  if(ncol(samples) != ncol(matrix) )
  {
    cat("
 The samples-matrix has a different number of columns than the input matrix.")
    cat("
 The number of columns of the input matrix is: " , ncol(matrix) , "
")
    cat("
 The number of columns of the sample-matrix is: " , ncol(samples) , "
")
    return(FALSE)
  }

  if(nrow(samples) != (2*reps) )
  {
    cat("
 The samples-matrix has a different number of rows than TWO TIMES the requested repetitions.")
    cat("
 The number of requested repetitions is: " , reps , "
")
    cat("
 Times two that makes: " , (2*reps) , "
")
    cat("
 The number of rows of the sample-matrix is: " , nrow(samples) , "
")
    return(FALSE)
  }

  cnm.input = colnames(matrix)
  cnm.samples = colnames(samples)

  for(i in 1:ncol(samples))
    if(cnm.input[i] != cnm.samples[i])
    {
      cat("
 The samples-matrix has at least one different column-name than the input matrix.")
      cat("
 The difference was observed at column #" , i , "
")
      cat("
 The column name of the input matrix for this column is: " , cnm.input[i] , "
")
      cat("
 The column name of the samples matrix for this column is: " , cnm.samples[i] , "
")
      return(FALSE)
    }

  c.sums = colSums(matrix)


  for( j in 1:ncol(samples) )
  {
    ubik = unique(matrix[,j])

    for( i in seq(1,nrow(samples),by=2) )
    {
      if( samples[i,j] %in% ubik == FALSE || samples[i+1,j] %in% ubik == FALSE ||
          (samples[i,j] == samples[i+1,j] && length(which(matrix[,j]==samples[i,j])) == 1 ))
      {
        cat("
 There is at least one sample that has an invalid element.")
        cat("
 The discrepancy was found in element: ", i, " , ", j, " of the samples matrix.")
        cat("
 The value of the samples matrix for this element is: " , samples[i,j])
        cat("
 The value of the samples matrix for element at the next row is: " , samples[i+1,j])
        cat("
 For this column, the elements of the matrix are: " , matrix[,j] , "
")
        return(FALSE)
      }

    } #for( j in 1:nrow(samples) )

  } #for( j in 1:ncol(samples) )

  samples.a = samples[seq(1,nrow(samples)-1,2),]
  samples.b = samples[seq(2,nrow(samples),2),]

  for( j in 1:ncol(samples.a) )
  {
    mean.new = mean(samples.a[,j])
    var.new = myvar(samples.a[,j])

    mean.actual = mean(matrix[,j])
    var.actual = myvar(matrix[,j])

    abs.error.mean = abs(mean.actual-mean.new)    
    rel.error.mean = 1.0
    
    if(mean.actual > 0.05)
      rel.error.mean = abs(mean.actual-mean.new)/abs(mean.actual)

    abs.error.sd = abs(sqrt(var.actual)-sqrt(var.new))
    rel.error.sd = 1.0

    if(var.actual > 0.05)
      rel.error.sd = abs(sqrt(var.actual)-sqrt(var.new))/abs(sqrt(var.actual))

    if(abs.error.mean > error && rel.error.mean > error  )
    {
      cat("
 The column mean computed by the beta-diversity permutation communities function",
            "
 (single, interleaved output) diverges a lot from the actual mean.")
      cat("
 This divergence was observed on the samples of the odd output rows. ")
      cat("
 The column where this was observed has index: ", j )
      cat("
 The value returned by the new function is: " , mean.new)
      cat("
 The actual mean is: " , mean.actual)
      cat("
 The relative difference is (percent): " , (rel.error.mean*100) , "
")
      return(FALSE)
    }

    if(abs.error.sd > error && rel.error.sd > error  )
    {
      cat("
 The column st.deviation computed by the beta-diversity permutation communities function",
            "
 (single matrix, single output) diverges a lot from the actual deviation.")
      cat("
 This divergence was observed on the samples of the odd output rows. ")
      cat("
 The column where this was observed has index: ", j )
      cat("
 The value computed by the new function is: " , sqrt(var.new))
      cat("
 The actual deviation is: " , sqrt(var.actual))
      cat("
 The relative difference is (percent): " , (rel.error.sd*100) , "
")
      return(FALSE)
    } 

  } # for( j in 1:ncol(samples.a) )


  for( j in 1:ncol(samples.b) )
  {
    mean.new = mean(samples.b[,j])
    var.new = myvar(samples.b[,j])

    mean.actual = mean(matrix[,j])
    var.actual = myvar(matrix[,j])

    abs.error.mean = abs(mean.actual-mean.new)    
    rel.error.mean = 1.0
    
    if(mean.actual > 0.05)
      rel.error.mean = abs(mean.actual-mean.new)/abs(mean.actual)

    abs.error.sd = abs(sqrt(var.actual)-sqrt(var.new))
    rel.error.sd = 1.0

    if(var.actual > 0.05)
      rel.error.sd = abs(sqrt(var.actual)-sqrt(var.new))/abs(sqrt(var.actual))

    if(abs.error.mean > error && rel.error.mean > error  )
    {
      cat("
 The column mean computed by the beta-diversity permutation communities function",
            "
 (single, interleaved output) diverges a lot from the actual mean.")
      cat("
 This divergence was observed on the samples of the even output rows. ")
      cat("
 The column where this was observed has index: ", j )
      cat("
 The value returned by the new function is: " , mean.new)
      cat("
 The actual mean is: " , mean.actual)
      cat("
 The relative difference is (percent): " , (rel.error.mean*100) , "
")
      return(FALSE)
    }

    if(abs.error.sd > error && rel.error.sd > error  )
    {
      cat("
 The column st.deviation computed by the beta-diversity permutation communities function",
            "
 (single matrix, single output) diverges a lot from the actual deviation.")
      cat("
 This divergence was observed on the samples of the even output rows. ")
      cat("
 The column where this was observed has index: ", j )
      cat("
 The value computed by the new function is: " , sqrt(var.new))
      cat("
 The actual deviation is: " , sqrt(var.actual))
      cat("
 The relative difference is (percent): " , (rel.error.sd*100) , "
")
      return(FALSE)
    } 

  } # for( j in 1:ncol(samples.b) )


  for( j in 1:ncol(samples.a) )
  {
    mean.new = mean(samples.a[,j]+samples.b[,j])
    var.new = myvar(samples.a[,j]+samples.b[,j])

    mean.actual = mean(colSums(combn(matrix[,j],2)))
    var.actual = myvar(colSums(combn(matrix[,j],2)))

    abs.error.mean = abs(mean.actual-mean.new)    
    rel.error.mean = 1.0
    
    if(mean.actual > 0.05)
      rel.error.mean = abs(mean.actual-mean.new)/abs(mean.actual)

    abs.error.sd = abs(sqrt(var.actual)-sqrt(var.new))
    rel.error.sd = 1.0

    if(var.actual > 0.05)
      rel.error.sd = abs(sqrt(var.actual)-sqrt(var.new))/abs(sqrt(var.actual))

    if(abs.error.mean > error && rel.error.mean > error  )
    {
      cat("
 The summed column mean computed by the beta-diversity permutation communities function",
            "
 (single, interleaved) diverges a lot from the actual mean.")
      cat("
 This divergence was observed on the samples of the second output matrix. ")
      cat("
 The column where this was observed has index: ", j )
      cat("
 The value returned by the new function is: " , mean.new)
      cat("
 The actual mean is: " , mean.actual)
      cat("
 The relative difference is (percent): " , (rel.error.mean*100) , "
")
      return(FALSE)
    }

    if(abs.error.sd > error && rel.error.sd > error  )
    {
      cat("
 The summed column st.deviation computed by the beta-diversity permutation communities function",
            "
 (single, interleaved) diverges a lot from the actual deviation.")
      cat("
 This divergence was observed on the samples of the second output matrix. ")
      cat("
 The column where this was observed has index: ", j )
      cat("
 The value computed by the new function is: " , sqrt(var.new))
      cat("
 The actual deviation is: " , sqrt(var.actual))
      cat("
 The relative difference is (percent): " , (rel.error.sd*100) , "
")
      return(FALSE)
    } 

  } # for( j in 1:ncol(samples.a) )

  return (TRUE)

} # permutation.communities.b.tester = function(...)

permutation.random.values.b.tester = function(matrix, f, args, error=0.1, reps=1000)
{
  #res.standard = moments.permutation.shuffling.standard.b.single.all.pairs(matrix,f,args,reps)

  res.standard = moments.permutation.shuffling.single.pair.b.single.all.pairs(matrix,f,args,reps)
  res.new = permutation.random.values.b(matrix,f,args,reps)

  if(length(res.new) != reps )
  {
    cat("
 The beta-diversity permutation moments function (single, all pairs)",
        "
 returns a wrong number of elements.")
    cat("
 The returned number of elements is: " , length(res.new) , "
")
    return(FALSE)
  }

  mean.standard = res.standard[[1]]
  var.standard = res.standard[[2]]

  mean.new = mean(res.new)
  var.new = var(res.new)

  error.mean = 0.0
  error.sd = 0.0

  if(mean.standard < 0.05)
    error.mean = abs(mean.standard-mean.new)
  else
    error.mean = abs(mean.standard-mean.new)/abs(mean.standard)

  if(var.standard < 0.05)
    error.sd = abs(sqrt(var.standard)-sqrt(var.new))
  else
    error.sd = abs(sqrt(var.standard)-sqrt(var.new))/abs(sqrt(var.standard))

  if(error.mean > error )
  {
    cat("
 The mean computed by the beta-diversity (single, all pairs) permutation",
        "
 moments function diverges a lot from the value computed by the standard method.")
    cat("
 The value returned by the new function is: " , mean.new )
    cat("
 The value returned by the standard function is: " , mean.standard )
    cat("
 The relative difference is (percent): " , (error.mean*100) , "
")
    return(FALSE)
  }

  if(error.sd > error )
  {
    cat("
 The st.deviation computed by the beta-diversity (single, all pairs) permutation",
        "
 moments function diverges a lot from the value computed by the standard method.")
    cat("
 The value computed by the new function is: " , sqrt(var.new))
    cat("
 The value computed by the standard function is: " , sqrt(var.standard))
    cat("
 The relative difference is (percent): " , (error.sd*100) , "
")
    return(FALSE)
  }

  return (TRUE)

} # permutation.random.values.b.tester = function(...)

permutation.moments.b.tester = function(matrix, f, args, error=0.1, reps=1000)
{
  #res.standard = moments.permutation.shuffling.standard.b.single.all.pairs(matrix,f,args,reps)
  res.standard = moments.permutation.shuffling.single.pair.b.single.all.pairs(matrix,f,args,reps)
  res.new = permutation.moments.b(matrix,f,args,reps)
 
  if(length(res.new) != 2 )
  {
    cat("
 The beta-diversity permutation moments function (single, all pairs)",
        "
 returns a wrong number of elements.")
    cat("
 The returned number of elements is: " , length(res.new) , "
")
    return(FALSE)
  }

  mean.standard = res.standard[[1]]
  var.standard = res.standard[[2]]

  mean.new = res.new[[1]]
  var.new = res.new[[2]]

  error.mean = 0.0
  error.sd = 0.0

  if(mean.standard < 0.05)
    error.mean = abs(mean.standard-mean.new)
  else
    error.mean = abs(mean.standard-mean.new)/abs(mean.standard)

  if(var.standard < 0.05)
    error.sd = abs(sqrt(var.standard)-sqrt(var.new))
  else
    error.sd = abs(sqrt(var.standard)-sqrt(var.new))/abs(sqrt(var.standard))

  if(error.mean > error )
  {
    cat("
 The mean computed by the beta-diversity (single, all pairs) permutation",
        "
 moments function diverges a lot from the value computed by the standard method.")
    cat("
 The value returned by the new function is: " , mean.new )
    cat("
 The value returned by the standard function is: " , mean.standard )
    cat("
 The relative difference is (percent): " , (error.mean*100) , "
")
    return(FALSE)
  }

  if(error.sd > error )
  {
    cat("
 The st.deviation computed by the beta-diversity (single, all pairs) permutation",
        "
 moments function diverges a lot from the value computed by the standard method.")
    cat("
 The value computed by the new function is: " , sqrt(var.new))
    cat("
 The value computed by the standard function is: " , sqrt(var.standard))
    cat("
 The relative difference is (percent): " , (error.sd*100) , "
")
    return(FALSE)
  }

  return (TRUE)

} # permutation.moments.b.tester = function(...)

permutation.pvalues.b.tester = function(matrix, f, args, error= 0.01, reps=1000)
{
  #pvalues.standard = pvalues.permutation.shuffling.standard.b.single.all.pairs(matrix,f,args,reps)

  observed.vals = as.vector(f(matrix,args))

  pvalues.standard = pvalues.permutation.shuffling.single.pair.b.single.all.pairs(matrix,f,args,observed.vals,reps)
  pvalues.new = permutation.pvalues.b(matrix, f, args, observed.vals, reps)

  if(length(pvalues.new) != length(observed.vals) )
  {
    cat("
 There was an error in the output of the beta-diversity permutation p-value") 
    cat("
 (single, all pairs) function.")
    cat("
 The number of returned p-values is different from the number of rows in the pairs matrix.")
    cat("
 The number of rows in the pairs matrix is: " , length(observed.vals))
    cat("
 The number of returned p-values is: " , length(pvalues.new) , "
")
    return(FALSE)
  }

  for(i in 1:length(pvalues.new))
    if( pvalues.new[i] <= 0.0 || pvalues.new[i] > 1.0 )
    {
      cat("
 There was an error in the output of the beta-diversity permutation p-value") 
      cat("
 (single, all pairs) function.")
      cat("
 At least of the returned p-values is out of bounds.")
      cat("
 This is the value with index: " , i)
      cat("
 This p-value is equal to: " , pvalues.new[i] , "
")
      return(FALSE)
    }

  error.pvalues = vector(mode="numeric",length(observed.vals))
  
  for( i in 1:length(observed.vals) )
    error.pvalues[i] = abs(pvalues.standard[i]-pvalues.new[i])


  for(i in 1:length(pvalues.new))
    if( error.pvalues[i] > error )
    {
      cat("
 There was an error in the output of the beta-diversity permutation p-value") 
      cat("
 (single, all pairs) function.")
      cat("
 At least one of the p-values computed by the function diverges",
            "
 a lot from the value computed by the standard method.")
      cat("
 This is the value with index: " , i)
      cat("
 The value returned by the new function is: " , pvalues.new[i])
      cat("
 The value returned by the standard function is: " , pvalues.standard[i])
      cat("
 The absolute error is: " , error.pvalues[i], "
")
      return(FALSE)
    }

  return (TRUE)
  
} # permutation.pvalues.b.tester=function(...)
