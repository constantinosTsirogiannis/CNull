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

require(Matrix)

construct.random.site.species.matrix.fixed.row.size = function(n.rows,n.cols,row.density,
                                                               max.abundance, column.names)
{
  r.mat = matrix(0,nrow=n.rows,ncol=n.cols)
  #r.mat = rsparsematrix(nrow=n.rows, ncol=n.cols, density=row.density)

  row.entries = n.cols*row.density

  colnames(r.mat) = column.names

  for(i in 1:n.rows)
  {
    number.of.nonzero.entries = row.entries 
    entries = sort(sample(1:n.cols,number.of.nonzero.entries))
    
    for(j in 1:length(entries))
    {
      k = entries[j]
      r.mat[i,k] = sample(1:max.abundance,1)
    }

  }
  
  return(r.mat)

} # construct.random.site.species.matrix.fixed.row.size


construct.random.site.species.matrix.with.column.names = function(n.rows,n.cols,row.density,
                                                                  max.abundance, column.names)
{
  r.mat = rsparsematrix(nrow=n.rows, ncol=n.cols, density=row.density)
  colnames(r.mat) = column.names

  for(i in 1:n.rows)
    for(j in 1:n.cols)
      if( r.mat[i,j] !=0)
        r.mat[i,j] = sample(1:max.abundance,1)
  
  return(as.matrix(r.mat))

} # construct.random.site.species.matrix


construct.random.site.species.matrix = function(n.rows,n.cols,row.density,
                                                max.abundance)
{
  col.names = paste("Species", 1:n.cols)    
 
  return(construct.random.site.species.matrix.with.column.names(n.rows,n.cols,row.density,
                                                                 max.abundance, col.names))

} # construct.random.site.species.matrix

construct.random.coordinates.matrix = function(n.rows,n.cols, row.names)
{
  r.mat = matrix(runif(0.0,1.0),nrow=n.rows, ncol=n.cols)
  rownames(r.mat) = row.names
  
  return(r.mat)

} # construct.random.coordinates.matrix
