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


############################
# Alpha diversity measures #
############################

species.richness.abundance.weighted.a = function(matrix,args)
{
  res <-rowSums(matrix)
  return(res)
}

###########################
# Beta diversity measures #
###########################

species.richness.abundance.weighted.b = function(mt,args)
{
  res.a <-rowSums(mt)
  res.c = matrix(0,nrow=nrow(mt),ncol=nrow(mt))

  for(i in 1:nrow(mt))
    for(j in 1:nrow(mt))
      res.c[i,j] = res.a[i]+res.a[j]

  return(res.c)
}
