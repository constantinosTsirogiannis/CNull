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

#####################################
#####################################
### Permutation shuffling methods ###
#####################################
#####################################

###############################
###############################
### Alpha diversity methods ###
###############################
###############################

communities.permutation.shuffling.standard.a=function(matrix,repetitions)
{
  for(i in 1:repetitions)
    temp.m = apply(matrix,2,sample)

} # communities.permutation.shuffling.standard.a(...)

moments.permutation.shuffling.standard.a=function(matrix,f,args,repetitions)
{
  r = nrow(matrix)
  s = ncol(matrix)

  sum.mean = rep(0.0,r)
  sum.var = rep(0.0,r)

  for(i in 1:repetitions)
  {
    temp.m = apply(matrix,2,sample)

    res = f(temp.m,args)        
    sum.mean = sum.mean+res
    sum.var = sum.var + (res*res)
  } 

  means = sum.mean/repetitions
  vars = (sum.var/repetitions) - (means*means)

  return(list(means,vars))

} # moments.permutation.shuffling.standard.a(...)

moments.permutation.shuffling.a.single.row=function(matrix,f,args,repetitions)
{
  mean = 0.0
  var = 0.0

  for(i in 1:repetitions)
  {
    temp.m = apply(matrix,2,sample)
    small.m = as.matrix(t(temp.m[1,]))

    res = f(small.m,args)        
    mean = mean+res[1]
    var = var + (res[1]*res[1])
  } 

  mean = mean/repetitions
  var = (var/repetitions) - (mean*mean)

  return(list(mean,var))

} # moments.permutation.shuffling.single.row.a(...)

pvalues.permutation.shuffling.a.standard=function(matrix,f,args,repetitions)
{
  r = nrow(matrix)
  s = ncol(matrix)

  vals = rep(1,r)
  
  actual.res = f(matrix,args)

  for(i in 1:repetitions)
  {
    temp.m = apply(matrix,2,sample)

    res = f(temp.m,args)        

    for(j in 1:r)
      if(actual.res[j] <= res[j])
        vals[j] = vals[j]+1
  } 

  vals = vals/(repetitions+1)

  return(vals)

} # pvalues.permutation.shuffling.standard.a(...)

pvalues.permutation.shuffling.a.single.row=function(matrix,f,args,repetitions)
{
  val = 1
  
  actual.res = f(as.matrix(t(matrix[1,])),args)

  for(i in 1:repetitions)
  {
    temp.m = apply(matrix,2,sample)
    small.m = as.matrix(t(temp.m[1,]))

    res = f(small.m,args)        

    if(actual.res[1] <= res[1])
      val=val+1
  } 

  val=val/(repetitions+1)

  return(val)

} # pvalues.permutation.shuffling.single.row.a(...)

##############################
##############################
### Beta diversity methods ###
##############################
##############################

############################################################
############################################################
## Beta-diversity functions, two matrices, specific pairs ##
############################################################
############################################################

communities.permutation.shuffling.standard.b.double=function(matrix.a, matrix.b,repetitions)
{
  for(i in 1:repetitions)
    temp.m = apply(matrix.a,2,sample)

  for(i in 1:repetitions)
    temp.m = apply(matrix.b,2,sample)

} # communities.permutation.shuffling.standard.b(...)

moments.permutation.shuffling.standard.b.double=function(matrix.a, matrix.b,pairs,f, args,repetitions)
{
  sum.mean = rep(0.0,nrow(pairs))
  sum.var = rep(0.0,nrow(pairs))

  for(i in 1:repetitions)
  {
    temp.m.a = apply(matrix.a,2,sample)
    temp.m.b = apply(matrix.b,2,sample)

    res = f(temp.m.a,temp.m.b,pairs,args)        
    sum.mean = sum.mean+res[1]
    sum.var = sum.var + (res[1]*res[1])
  } 

  means = sum.mean/repetitions
  vars = (sum.var/repetitions) - (means*means)

  return(list(means,vars))

} # moments.permutation.shuffling.standard.b.double(...)

moments.permutation.shuffling.single.pair.b.double=function(matrix.a, matrix.b,f,args,repetitions)
{
  mean = 0.0
  var = 0.0

  pairs = matrix(1, nrow=1, ncol=2)

  for(i in 1:repetitions)
  {
    temp.m.a = apply(matrix.a,2,sample)
    temp.m.b = apply(matrix.b,2,sample)
    small.m.a = as.matrix(t(temp.m.a[1,]))
    small.m.b = as.matrix(t(temp.m.b[1,]))

    res = f(small.m.a,small.m.b,pairs,args)        

    mean = mean+res[1]
    var = var + (res[1]*res[1])
  } 

  mean = mean/repetitions
  var = (var/repetitions) - (mean*mean)

  return(list(mean,var))

} # moments.permutation.shuffling.single.pair.b.double(...)


pvalues.permutation.shuffling.standard.b.double=function(matrix.a, matrix.b,pairs,f,args,repetitions)
{
  vals = rep(1,nrow(pairs))
  
  actual.res = f(matrix.a,matrix.b, pairs,args)

  for(i in 1:repetitions)
  {
    temp.m.a = apply(matrix.a,2,sample)
    temp.m.b = apply(matrix.b,2,sample)

    res = f(temp.m.a,temp.m.b,pairs,args)        

    for(j in 1:nrow(pairs))
      if(actual.res[j] <= res[j])
        vals[j] = vals[j]+1
  } 

  vals = vals/(repetitions+1)

  return(vals)

} # pvalues.permutation.shuffling.standard.b.double(...)


pvalues.permutation.shuffling.single.pair.b.double=function(matrix.a, matrix.b,pairs,f, args,repetitions)
{
  val = rep(1,nrow(pairs))
  
  small.a = as.matrix(t(matrix.a[1,]))
  small.b = as.matrix(t(matrix.b[1,]))

  actual.res = f(small.a,small.b,pairs,args)
  temp.pairs = matrix(1,nrow=1,ncol=2)

  for(i in 1:repetitions)
  {
    temp.m.a = apply(matrix.a,2,sample)
    temp.m.b = apply(matrix.b,2,sample)

    small.m.a = as.matrix(t(temp.m.a[1,]))
    small.m.b = as.matrix(t(temp.m.b[1,]))

    res = f(small.m.a, small.m.b, temp.pairs, args)        

    for(j in 1:length(actual.res))
      if(actual.res[j] <= res[1])
        val[j]=val[j]+1
  } 

  val=val/(repetitions+1)

  return(val)

} # pvalues.permutation.shuffling.single.pair.b.double(...)

#############################################################
#############################################################
## Beta-diversity functions, single matrix, specific pairs ##
#############################################################
#############################################################

#########################################################
## For the next function:                              ##
## Precondition: pairs[i,1] != pairs[i,2] for every i. ##
#########################################################

moments.permutation.shuffling.standard.b.single=function(matrix,pairs,f,args,repetitions)
{
  sum.mean = rep(0.0,nrow(pairs))
  sum.var = rep(0.0,nrow(pairs))

  for(i in 1:repetitions)
  {
    temp.m = apply(matrix,2,sample)

    res = f(temp.m,pairs,args)        
    sum.mean = sum.mean+res[1]
    sum.var = sum.var + (res[1]*res[1])
  } 

  means = sum.mean/repetitions
  vars = (sum.var/repetitions) - (means*means)

  return(list(means,vars))

} # moments.permutation.shuffling.standard.b.single(...)

moments.permutation.shuffling.single.pair.b.single=function(matrix,f,args,repetitions)
{
  mean = 0.0
  var = 0.0

  pairs = matrix(1:2, nrow=1, ncol=2)

  for(i in 1:repetitions)
  {
    temp.m = apply(matrix,2,sample)
    small.m = as.matrix(temp.m[1:2,])

    res = f(small.m,pairs,args)        
    mean = mean+res[1]
    var = var + (res[1]*res[1])
  } 

  mean = mean/repetitions
  var = (var/repetitions) - (mean*mean)

  return(list(mean,var))

} # moments.permutation.shuffling.single.pair.b.single(...)

#########################################################
## For the next function:                              ##
## Precondition: pairs[i,1] != pairs[i,2] for every i. ##
#########################################################

pvalues.permutation.shuffling.standard.b.single=function(matrix,pairs,f,args,repetitions)
{
  vals = rep(1,nrow(pairs))
  
  actual.res = f(matrix, pairs,args)

  for(i in 1:repetitions)
  {
    temp.m = apply(matrix,2,sample)

    res = f(temp.m,pairs,args)        

    for(j in 1:nrow(pairs))
      if(actual.res[j] <= res[j])
        vals[j] = vals[j]+1
  } 

  vals = vals/(repetitions+1)

  return(vals)

} # pvalues.permutation.shuffling.standard.b.single(...)


pvalues.permutation.shuffling.single.pair.b.single=function(matrix,pairs,f,args,repetitions)
{
  val = rep(1,nrow(pairs))
  
  actual.res = f(matrix,pairs,args)
  temp.pairs = matrix(1:2,nrow=1,ncol=2)

  for(i in 1:repetitions)
  {
    temp.m = apply(matrix,2,sample)

    small.m = as.matrix(temp.m[1:2,])
    res = f(small.m, temp.pairs, args)        

    for(j in 1:length(actual.res))
    if(actual.res[j] <= res[1])
      val[j]=val[j]+1
  } 

  val=val/(repetitions+1)

  return(val)

} # pvalues.permutation.shuffling.single.pair.b.single(...)


#######################################################
#######################################################
## Beta-diversity functions, two matrices, all pairs ##
#######################################################
#######################################################

moments.permutation.shuffling.standard.b.double.all.pairs=function(matrix.a, matrix.b,f,args,repetitions)
{
  sum.mean = rep(0.0,nrow(matrix.a)*nrow(matrix.b))
  sum.var = rep(0.0,nrow(matrix.a)*nrow(matrix.b))

  for(i in 1:repetitions)
  {
    temp.m.a = apply(matrix.a,2,sample)
    temp.m.b = apply(matrix.b,2,sample)

    res = f(temp.m.a,temp.m.b,args)        
    sum.mean = sum.mean+res
    sum.var = sum.var + (res*res)
  } 

  means = sum.mean/repetitions
  vars = (sum.var/repetitions) - (means*means)

  return(list(means[1],vars[1]))

} # moments.permutation.shuffling.standard.b.double.all.pairs(...)

moments.permutation.shuffling.single.pair.b.double.all.pairs=function(matrix.a, matrix.b,f,args,repetitions)
{
  mean = 0.0
  var = 0.0

  for(i in 1:repetitions)
  {
    temp.m.a = apply(matrix.a,2,sample)
    temp.m.b = apply(matrix.b,2,sample)
    small.m.a = t(as.matrix(temp.m.a[1,]))
    small.m.b = t(as.matrix(temp.m.b[1,]))

    res = f(small.m.a,small.m.b,args)        
    mean = mean+res[1]
    var = var + (res[1]*res[1])
  } 

  mean = mean/repetitions
  var = (var/repetitions) - (mean*mean)

  return(list(mean,var))

} # moments.permutation.shuffling.single.pair.b.double.all.pairs(...)


pvalues.permutation.shuffling.standard.b.double.all.pairs=function(matrix.a, matrix.b,f,args,repetitions)
{
  vals = rep(1,nrow(matrix.a)*nrow(matrix.b))
  
  actual.res = f(matrix.a,matrix.b,args)

  for(i in 1:repetitions)
  {
    temp.m.a = apply(matrix.a,2,sample)
    temp.m.b = apply(matrix.b,2,sample)

    res = f(temp.m.a,temp.m.b,args)        

    for(j in 1:(nrow(matrix.a)*nrow(matrix.b)))
      if(actual.res[j] <= res[j])
        vals[j] = vals[j]+1
  } 

  vals = vals/(repetitions+1)

  return(vals)

} # pvalues.permutation.shuffling.standard.b.double.all.pairs(...)


pvalues.permutation.shuffling.single.pair.b.double.all.pairs=function(matrix.a, matrix.b,f,args,observed.vals,repetitions)
{
  val = rep(1,length(observed.vals))
  
  for(i in 1:repetitions)
  {
    temp.m.a = apply(matrix.a,2,sample)
    temp.m.b = apply(matrix.b,2,sample)

    small.m.a = t(as.matrix(temp.m.a[1,]))
    small.m.b = t(as.matrix(temp.m.b[1,]))

    res = f(small.m.a, small.m.b, args)        

    for(j in 1:length(observed.vals))
      if(observed.vals[j] <= res[1])
        val[j]=val[j]+1
  } 

  val=val/(repetitions+1)

  return(val)

} # pvalues.permutation.shuffling.single.pair.b.double.all.pairs(...)

########################################################
########################################################
## Beta-diversity functions, single matrix, all pairs ##
########################################################
########################################################

#####################################################
# The next functions should work fine if the chosen #
# functin f returns an object of class dist.        #
#####################################################

moments.permutation.shuffling.standard.b.single.all.pairs=function(matrix,f, args,repetitions)
{
  sum.mean = 0.0
  sum.var = 0.0

  for(i in 1:repetitions)
  {
    temp.m = apply(matrix,2,sample)

    res = f(temp.m,args)        
    sum.mean = sum.mean+res[1,1]
    sum.var = sum.var + (res[1,1]*res[1,1])
  } 

  means = sum.mean/repetitions
  vars = (sum.var/repetitions) - (means*means)

  return(list(means,vars))

} # moments.permutation.shuffling.standard.b.single(...)

moments.permutation.shuffling.single.pair.b.single.all.pairs=function(mt,f, args,repetitions)
{
  sum.mean = 0.0
  sum.var = 0.0

  for(i in 1:repetitions)
  {
    temp.m = apply(mt,2,sample)

    res = f(as.matrix(temp.m[1:2,]),args)
        
    sum.mean = sum.mean+res[1,1]
    sum.var = sum.var + (res[1,1]*res[1,1])
  } 

  means = sum.mean/repetitions
  vars = (sum.var/repetitions) - (means*means)

  return(list(means,vars))

} # moments.permutation.shuffling.standard.b.single(...)

pvalues.permutation.shuffling.standard.b.single=function(matrix,f,args,repetitions)
{
  n = nrow(matrix)*(nrow(matrix)-1)/2

  vals = rep(1,n)
  
  actual.res = f(matrix, args)

  for(i in 1:repetitions)
  {
    temp.m = apply(matrix,2,sample)
    
    res = f(temp.m,args)        

    for(j in 1:n)
      if(actual.res[j] <= res[1,1])
        vals[j] = vals[j]+1
  } 

  vals = vals/(repetitions+1)

  return(vals)

} # pvalues.permutation.shuffling.standard.b.single(...)

pvalues.permutation.shuffling.single.pair.b.single.all.pairs=function(matrix,f,args,observed.vals,repetitions)
{
  vals = rep(1,length(observed.vals))

  for(i in 1:repetitions)
  {
    temp.m = apply(matrix,2,sample)
    
    res = f(as.matrix(temp.m[1:2,]),args)        

    for(j in 1:length(observed.vals))
      if(observed.vals[j] <= res[1,1])
        vals[j] = vals[j]+1
  } 

  vals = vals/(repetitions+1)

  return(vals)

} # pvalues.permutation.shuffling.single.pair.b.single(...)

##########################################
##########################################
### Individual-based shuffling methods ###
##########################################
##########################################

###############################
###############################
### Alpha diversity methods ###
###############################
###############################

individual.based.perturbation=function(m)
{
  csums = colSums(m)
  r = nrow(m)
  s = ncol(m)

  new.m = matrix(0, nrow = nrow(m), ncol = ncol(m) )

  for(j in seq(from=1,to=s, by=1))
    if(csums[j] > 0)
      for(k in seq(from=1,to=csums[j], by=1))
      {
        res = sample(1:r,1)
        new.m[res,j] = new.m[res,j]+1
      }

  return (new.m)

} # individual.based.perutrbation(m)


communities.individual.based.shuffling.standard.a=function(matrix,repetitions)
{
  csums = colSums(matrix)
  r = nrow(matrix)
  s = ncol(matrix)

  for(i in 1:repetitions)
    individual.based.perturbation(matrix)

} # communities.individual.based.shuffling.standard.a(...)


moments.individual.based.shuffling.a.standard=function(matrix,f,args,repetitions)
{
  csums = colSums(matrix)
  r = nrow(matrix)
  s = ncol(matrix)

  sum.mean = rep(0.0,r)
  sum.var = rep(0.0,r)

  for(i in 1:repetitions)
  {
    m = individual.based.perturbation(matrix)

    res = f(m,args)        
    sum.mean = sum.mean + res[1]
    sum.var = sum.var + (res[1]*res[1])

  } # for(i in 1:repetitions)

  means = sum.mean/repetitions
  vars = (sum.var/repetitions) - (means*means)

  return(list(means,vars))

} # moments.individual.based.shuffling.a(...)

moments.individual.based.shuffling.a.single.row=function(matrix,f,args,repetitions)
{
  csums = colSums(matrix)
  r = nrow(matrix)
  s = ncol(matrix)

  sum.mean = 0
  sum.var = 0

  for(i in 1:repetitions)
  {
    m = individual.based.perturbation(matrix)
    res = f(t(as.matrix((m[1,]))),args)
        
    sum.mean = sum.mean+res[1]
    sum.var = sum.var + (res[1]*res[1])

  } # for(i in 1:repetitions)

  mean = sum.mean/repetitions
  vr = (sum.var/repetitions) - (mean*mean)

  return(list(mean,vr))

} # moments.individual.based.shuffling.a.single.row(...)

pvalues.individual.based.shuffling.a.standard=function(matrix,f,args,repetitions)
{
  r = nrow(matrix)
  s = ncol(matrix)

  vals = rep(1,r)
  
  actual.res = f(matrix,args)

  for(i in 1:repetitions)
  {
    m = individual.based.perturbation(matrix)

    res = f(m,args)        

    for(j in 1:r)
      if(actual.res[j] <= res[j])
        vals[j] = vals[j]+1

  } # for(i in 1:repetitions)

  vals = vals/(repetitions+1)

  return(vals)

} # pvalues.individual.based.shuffling.a(...)


pvalues.individual.based.shuffling.a.single.row=function(matrix,f,args,repetitions)
{
  r = nrow(matrix)
  s = ncol(matrix)

  val = rep(1,r)
  
  actual.res = f(matrix,args)

  for(i in 1:repetitions)
  {
    m = individual.based.perturbation(matrix)

    res = f(as.matrix(t(m[1,])),args)        

    for(k in 1:r)    
      if(actual.res[k] <= res[1])
        val[k] = val[k]+1

  } # for(i in 1:repetitions)

  val = val/(repetitions+1)

  return(val)

} # pvalues.individual.based.shuffling.a.single.row(...)


##############################
##############################
### Beta diversity methods ###
##############################
##############################

############################################################
############################################################
## Beta-diversity functions, two matrices, specific pairs ##
############################################################
############################################################

communities.individual.based.shuffling.standard.b.double=function(matrix.a, matrix.b,repetitions)
{
  for(i in 1:repetitions)
    individual.based.perturbation(matrix.a)

  for(i in 1:repetitions)
    individual.based.perturbation(matrix.b)

} # communities.individual.based.shuffling.standard.b(...)

moments.individual.based.shuffling.standard.b.double=function(matrix.a, matrix.b,pairs,f, args,repetitions)
{
  sum.mean = rep(0.0,nrow(pairs))
  sum.var = rep(0.0,nrow(pairs))

  for(i in 1:repetitions)
  {
    temp.m.a = individual.based.perturbation(matrix.a)
    temp.m.b = individual.based.perturbation(matrix.b)

    res = f(temp.m.a,temp.m.b,pairs,args)        
    sum.mean = sum.mean+res[1]
    sum.var = sum.var + (res[1]*res[1])
  } 

  means = sum.mean/repetitions
  vars = (sum.var/repetitions) - (means*means)

  return(list(means,vars))

} # moments.individual.based.shuffling.standard.b.double(...)

moments.individual.based.shuffling.single.pair.b.double=function(matrix.a, matrix.b,f,args,repetitions)
{
  mean = 0.0
  var = 0.0

  pairs = matrix(1, nrow=1, ncol=2)

  for(i in 1:repetitions)
  {
    temp.m.a = individual.based.perturbation(matrix.a)
    temp.m.b = individual.based.perturbation(matrix.b)
    small.m.a = as.matrix(t(temp.m.a[1,]))
    small.m.b = as.matrix(t(temp.m.b[1,]))

    res = f(small.m.a,small.m.b,pairs,args)        
    mean = mean+res[1]
    var = var + (res[1]*res[1])
  } 

  mean = mean/repetitions
  var = (var/repetitions) - (mean*mean)

  return(list(mean,var))

} # moments.individual.based.shuffling.single.pair.b.double(...)


pvalues.individual.based.shuffling.standard.b.double=function(matrix.a, matrix.b,pairs,f,args,repetitions)
{
  vals = rep(1,nrow(pairs))
  
  actual.res = f(matrix.a,matrix.b, pairs,args)

  for(i in 1:repetitions)
  {
    temp.m.a = individual.based.perturbation(matrix.a)
    temp.m.b = individual.based.perturbation(matrix.b)

    res = f(temp.m.a,temp.m.b,pairs,args)        

    for(j in 1:nrow(pairs))
      if(actual.res[j] <= res[j])
        vals[j] = vals[j]+1
  } 

  vals = vals/(repetitions+1)

  return(vals)

} # pvalues.individual.based.shuffling.standard.b.double(...)


pvalues.individual.based.shuffling.single.pair.b.double=function(matrix.a, matrix.b,pairs,f,args,repetitions)
{
  val = rep(1,nrow(pairs))

  actual.res = f(matrix,pairs,args)
  temp.pairs = matrix(1:2,nrow=1,ncol=2)

  for(i in 1:repetitions)
  {
    temp.m.a = individual.based.perturbation(matrix.a)
    temp.m.b = individual.based.perturbation(matrix.b)

    small.m.a = as.matrix(t(temp.m.a[1,]))
    small.m.b = as.matrix(t(temp.m.b[1,]))

    res = f(small.m.a, small.m.b, temp.pairs, args)        

    for(j in 1:length(actual.res))
      if(actual.res[j] <= res[1])
        val[j]=val[j]+1
  } 

  val=val/(repetitions+1)

  return(val)

} # pvalues.individual.based.shuffling.single.pair.b.double(...)

#############################################################
#############################################################
## Beta-diversity functions, single matrix, specific pairs ##
#############################################################
#############################################################

#########################################################
## For the next function:                              ##
## Precondition: pairs[i,1] != pairs[i,2] for every i. ##
#########################################################

moments.individual.based.shuffling.standard.b.single=function(matrix,pairs,f,args,repetitions)
{
  sum.mean = rep(0.0,nrow(pairs))
  sum.var = rep(0.0,nrow(pairs))

  for(i in 1:repetitions)
  {
    temp.m = individual.based.perturbation(matrix)

    res = f(temp.m,pairs,args)        
    sum.mean = sum.mean+res[1]
    sum.var = sum.var + (res[1]*res[1])
  } 

  means = sum.mean/repetitions
  vars = (sum.var/repetitions) - (means*means)

  return(list(means,vars))

} # moments.individual.based.shuffling.standard.b.single(...)

moments.individual.based.shuffling.single.pair.b.single=function(matrix,f,args,repetitions)
{
  mean = 0.0
  var = 0.0

  pairs = matrix(1:2, nrow=1, ncol=2)

  for(i in 1:repetitions)
  {
    temp.m = individual.based.perturbation(matrix)
    small.m = as.matrix(temp.m[1:2,])

    res = f(small.m,pairs,args)        
    mean = mean+res[1]
    var = var + (res[1]*res[1])
  } 

  mean = mean/repetitions
  var = (var/repetitions) - (mean*mean)

  return(list(mean,var))

} # moments.individual.based.shuffling.single.pair.b.single(...)

#########################################################
## For the next function:                              ##
## Precondition: pairs[i,1] != pairs[i,2] for every i. ##
#########################################################

pvalues.individual.based.shuffling.standard.b.single=function(matrix,pairs,f,args,repetitions)
{
  vals = rep(1,nrow(pairs))
  
  actual.res = f(matrix, pairs,args)

  for(i in 1:repetitions)
  {
    temp.m = individual.based.perturbation(matrix)

    res = f(temp.m,pairs,args)        

    for(j in 1:nrow(pairs))
      if(actual.res[j] <= res[j])
        vals[j] = vals[j]+1
  } 

  vals = vals/(repetitions+1)

  return(vals)

} # pvalues.individual.based.shuffling.standard.b.single(...)


pvalues.individual.based.shuffling.single.pair.b.single=function(matrix,pairs,f,args,repetitions)
{
  val = rep(1,nrow(pairs))
  
  actual.res = f(matrix,pairs,args)
  temp.pairs = matrix(1:2,nrow=1,ncol=2)

  for(i in 1:repetitions)
  {
    temp.m = individual.based.perturbation(matrix)

    small.m = as.matrix(temp.m[1:2,])
    res = f(small.m, temp.pairs, args)        

    for(j in 1:length(actual.res))
      if(actual.res[j] <= res[1])
        val[j]=val[j]+1
  } 

  val=val/(repetitions+1)

  return(val)

} # pvalues.individual.based.shuffling.single.pair.b.single(...)


#######################################################
#######################################################
## Beta-diversity functions, two matrices, all pairs ##
#######################################################
#######################################################

moments.individual.based.shuffling.standard.b.double.all.pairs=function(matrix.a, matrix.b,f,args,repetitions)
{
  sum.mean = rep(0.0,nrow(matrix.a)*nrow(matrix.b))
  sum.var = rep(0.0,nrow(matrix.a)*nrow(matrix.b))

  for(i in 1:repetitions)
  {
    temp.m.a = individual.based.perturbation(matrix.a)
    temp.m.b = individual.based.perturbation(matrix.b)

    res = f(temp.m.a,temp.m.b,args)        
    sum.mean = sum.mean+res
    sum.var = sum.var + (res*res)
  } 

  means = sum.mean/repetitions
  vars = (sum.var/repetitions) - (means*means)

  return(list(means[1],vars[1]))

} # moments.individual.based.shuffling.standard.b.double.all.pairs(...)

moments.individual.based.shuffling.single.pair.b.double.all.pairs=function(matrix.a, matrix.b,f,args,repetitions)
{
  mean = 0.0
  var = 0.0

  for(i in 1:repetitions)
  {
    temp.m.a = individual.based.perturbation(matrix.a)
    temp.m.b = individual.based.perturbation(matrix.b)
    small.m.a = as.matrix(t(temp.m.a[1,]))
    small.m.b = as.matrix(t(temp.m.b[1,]))

    res = f(small.m.a,small.m.b,args)        
    mean = mean+res[1]
    var = var + (res[1]*res[1])
  } 

  mean = mean/repetitions
  var = (var/repetitions) - (mean*mean)

  return(list(mean,var))

} # moments.individual.based.shuffling.single.pair.b.double.all.pairs(...)


pvalues.individual.based.shuffling.standard.b.double.all.pairs=function(matrix.a, matrix.b,f,args,observed.vals,repetitions)
{
  vals = rep(1,length(observed.vals))

  for(i in 1:repetitions)
  {
    temp.m.a = individual.based.perturbation(matrix.a)
    temp.m.b = individual.based.perturbation(matrix.b)

    res = f(temp.m.a,temp.m.b,args)        

    for(j in 1:length(observed.vals))
      if(observed.vals[j] <= res[j])
        vals[j] = vals[j]+1
  } 

  vals = vals/(repetitions+1)

  return(vals)

} # pvalues.individual.based.shuffling.standard.b.double.all.pairs(...)


pvalues.individual.based.shuffling.single.pair.b.double.all.pairs=function(matrix.a, matrix.b,f,args,observed.vals,repetitions)
{
  val = rep(1,length(observed.vals))
  
  for(i in 1:repetitions)
  {
    temp.m.a = individual.based.perturbation(matrix.a)
    temp.m.b = individual.based.perturbation(matrix.b)

    small.m.a = as.matrix(t(temp.m.a[1,]))
    small.m.b = as.matrix(t(temp.m.b[1,]))

    res = f(small.m.a, small.m.b, args)        

    for(j in 1:length(observed.vals))
      if(observed.vals[j] <= res[1])
        val[j]=val[j]+1
  } 

  val=val/(repetitions+1)

  return(val)

} # pvalues.individual.based.shuffling.single.pair.b.double.all.pairs(...)

########################################################
########################################################
## Beta-diversity functions, single matrix, all pairs ##
########################################################
########################################################

#####################################################
# The next functions should work fine if the chosen #
# functin f returns an object of class dist.        #
#####################################################

moments.individual.based.shuffling.standard.b.single.all.pairs=function(matrix,f, args,repetitions)
{
  sum.mean= 0.0
  sum.var = 0.0

  for(i in 1:repetitions)
  {
    temp.m = individual.based.perturbation(matrix)

    res = f(temp.m,args)        
    sum.mean = sum.mean + res[1,1]
    sum.var = sum.var + (res[1,1]*res[1,1])
  } 

  means = sum.mean/repetitions
  vars = (sum.var/repetitions) - (means*means)

  return(list(means[1],vars[1]))

} # moments.individual.based.shuffling.standard.b.single(...)

moments.individual.based.shuffling.single.pair.b.single.all.pairs=function(matrix,f, args,repetitions)
{
  sum.mean = 0.0
  sum.var = 0.0

  for(i in 1:repetitions)
  {
    temp.m = individual.based.perturbation(matrix)

    res = f(as.matrix(temp.m[1:2,]),args)        
    sum.mean = sum.mean+res[1,1]
    sum.var = sum.var + (res[1,1]*res[1,1])
  } 

  means = sum.mean/repetitions
  vars = (sum.var/repetitions) - (means*means)

  return(list(means,vars))

} # moments.individual.based.shuffling.standard.b.single(...)

pvalues.individual.based.shuffling.standard.b.single.all.pairs=function(matrix,f,args,observed.vals,repetitions)
{
  vals = rep(1,length(observed.vals))
  
  for(i in 1:repetitions)
  {
    temp.m = individual.based.perturbation(matrix)
    
    res = f(temp.m,args)        

    for(j in 1:length(observed.vals))
      if(observed.vals[j] <= res[1,1])
        vals[j] = vals[j]+1
  } 

  vals = vals/(repetitions+1)

  return(vals)

} # pvalues.individual.based.shuffling.standard.b.single.all.pairs(...)

pvalues.individual.based.shuffling.single.pair.b.single.all.pairs=function(matrix,f,args,observed.vals,repetitions)
{
  vals = rep(1,length(observed.vals))
  
  for(i in 1:repetitions)
  {
    temp.m = individual.based.perturbation(matrix)
    
    res = f(as.matrix(temp.m[1:2,]),args)        

    for(j in 1:length(observed.vals))
      if(observed.vals[j] <= res[1,1])
        vals[j] = vals[j]+1
  } 

  vals = vals/(repetitions+1)

  return(vals)

} # pvalues.individual.based.shuffling.single.pair.b.single.all.pairs(...)
