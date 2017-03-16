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

require(CNull)

cat("
")

source("function_interfaces.R")
source("standard_shuffling_methods.R")
source("randomized_matrix_constructors.R")
source("test_code_individual_based.R")
source("test_code_permutation.R")

##################################################
##################################################
# Auxilliary function for adaptive recomputation #
##################################################
##################################################

recompute = function(f,reps,...)
{
  prw = 1

  for( i in 1:3)
  {    
    res = f(...,reps*prw)

    if(res==TRUE)
      return (TRUE)

    cat("               A discrepancy was observed, recomputing ...
")

    prw = prw*2

  } # for( i in 1:3)

  return (FALSE)

} # recompute = function(...)

###########################################
###########################################
# Test functions for permutation sampling #
###########################################
###########################################

source("test_code_permutation.R")

################################################
################################################
# Test functions for individual-based sampling #
################################################
################################################

source("test_code_individual_based.R")

#########################
# Create error log file #
#########################

if(file.exists("error_log.txt"))
  a=file.remove("error_log.txt")

b=file.create("error_log.txt")

#######################################################################
#######################################################################
#######################################################################
#######################  Main test suite code #########################
#######################################################################
#######################################################################
#######################################################################

cat("
")
cat("
")

cat("                     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
")
cat("                     ~ Constructing test matrices. ~ 
")
cat("                     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
")

original.matrix.a = construct.random.site.species.matrix(30,10,0.3,5)
original.matrix.b = construct.random.site.species.matrix(15,5,0.2,3)
original.matrix.c = construct.random.site.species.matrix(10,5,0.2,3)
reps=1000
func.a = species.richness.abundance.weighted.a
func.b = species.richness.abundance.weighted.b

pairs.a.to.a = matrix(0,nrow = nrow(original.matrix.a),ncol=2)
pairs.a.to.a[,1] = sample(1:nrow(original.matrix.a),size = nrow(original.matrix.a), replace = TRUE)
pairs.a.to.a[,2] = sample(1:nrow(original.matrix.a),size = nrow(original.matrix.a), replace = TRUE)

pairs.a.to.b = matrix(0,nrow = nrow(original.matrix.a),ncol=2)
pairs.a.to.b[,1] = sample(1:nrow(original.matrix.a),size = nrow(original.matrix.a), replace = TRUE)
pairs.a.to.b[,2] = sample(1:nrow(original.matrix.b),size = nrow(original.matrix.a), replace = TRUE)

pairs.b.to.b = matrix(0,nrow = nrow(original.matrix.b),ncol=2)
pairs.b.to.b[,1] = sample(1:nrow(original.matrix.b),size = nrow(original.matrix.b), replace = TRUE)
pairs.b.to.b[,2] = sample(1:nrow(original.matrix.b),size = nrow(original.matrix.b), replace = TRUE)

pairs.b.to.c = matrix(0,nrow = nrow(original.matrix.b),ncol=2)
pairs.b.to.c[,1] = sample(1:nrow(original.matrix.b),size = nrow(original.matrix.b), replace = TRUE)
pairs.b.to.c[,2] = sample(1:nrow(original.matrix.c),size = nrow(original.matrix.b), replace = TRUE)

args = list()
error=0.05
set.seed(7)

cat("                       ~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
")
cat("                       ~ Commencing experiments. ~ 
")
cat("                       ~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
")

cat("     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
")
cat("     ~ Executing experiments with permutation sampling functions. ~ 
")
cat("     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
")

cat("        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
")
cat("        ~ Executing experiments with alpha diversity functions. ~ 
")
cat("        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
")

res = recompute(permutation.communities.a.tester,10*reps,original.matrix.a,error)

if(res == FALSE)
  stop("Errors were detected in the tests, aborting...")

cat("                                    . 
")

res = recompute(permutation.random.values.a.tester,reps,original.matrix.a, func.a, args, error)

if(res == FALSE)
  stop("Errors were detected in the tests, aborting...")

cat("                                   . . 
")

res = recompute(permutation.moments.a.tester,reps,original.matrix.a, func.a, args, error)

if(res == FALSE)
  stop("Errors were detected in the tests, aborting...")

cat("                                  . . . 
")

res = recompute(permutation.pvalues.a.tester,reps,original.matrix.a, func.a, args, error)

if(res == FALSE)
  stop("Errors were detected in the tests, aborting...")


cat("         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
")
cat("         ~ Executing experiments with beta diversity functions. ~ 
")
cat("         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
")

res = recompute(permutation.communities.b.tester,10*reps,original.matrix.a,error)

if(res == FALSE)
  stop("Errors were detected in the tests, aborting...")

cat("                                    . 
")

res = recompute(permutation.random.values.b.tester, reps,original.matrix.a, func.b, args, error)

if(res == FALSE)
  stop("Errors were detected in the tests, aborting...")

cat("                                   . . 
")

res = recompute(permutation.moments.b.tester, reps,original.matrix.a, func.b, args, error)

if(res == FALSE)
  stop("Errors were detected in the tests, aborting...")

cat("                                  . . . 
")

res = recompute(permutation.pvalues.b.tester, reps,original.matrix.a, func.b, args, error)

if(res == FALSE)
  stop("Errors were detected in the tests, aborting...")

cat("   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
")
cat("   ~ Executing experiments with individual-based sampling functions. ~ 
")
cat("   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
")

cat("        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
")
cat("        ~ Executing experiments with alpha diversity functions. ~ 
")
cat("        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
")

res = recompute(individual.based.communities.a.tester,10*reps,original.matrix.a,error)

if(res == FALSE)
  stop("Errors were detected in the tests, aborting...")

cat("                                    . 
")


res = recompute(individual.based.random.values.a.tester, reps,original.matrix.b, func.a, args, error)

cat("                                   . . 
")

if(res == FALSE)
  stop("Errors were detected in the tests, aborting...")

res = recompute(individual.based.moments.a.tester, reps,original.matrix.b, func.a, args, error)

cat("                                  . . . 
")

if(res == FALSE)
  stop("Errors were detected in the tests, aborting...")

res = recompute(individual.based.pvalues.a.tester, reps,original.matrix.b, func.a, args, error)

if(res == FALSE)
  stop("Errors were detected in the tests, aborting...")

cat("         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
")
cat("         ~ Executing experiments with beta diversity functions. ~ 
")
cat("         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
")

res = recompute(individual.based.communities.b.tester,10*reps,original.matrix.a,error)

if(res == FALSE)
  stop("Errors were detected in the tests, aborting...")

cat("                                    . 
")

res = recompute(individual.based.random.values.b.tester,reps,original.matrix.b, func.b, args, error)

if(res == FALSE)
  stop("Errors were detected in the tests, aborting...")

cat("                                   . . 
")

res = recompute(individual.based.moments.b.tester,reps,original.matrix.b, func.b, args, error)

if(res == FALSE)
  stop("Errors were detected in the tests, aborting...")

cat("                                  . . . 
")

res = recompute(individual.based.pvalues.b.tester,reps,original.matrix.b, func.b,args,  error)

if(res == FALSE)
  stop("Errors were detected in the tests, aborting...")

cat("               ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
")
cat("               ~ All tests were executed successfully. ~ 
")
cat("               ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
 
")

#######################################################################
