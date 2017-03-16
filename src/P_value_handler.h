////////////////////////////////////////////////////////////////////////////////////////
//  Copyright (C) 2017,  Constantinos Tsirogiannis.  Email: tsirogiannis.c@gmail.com  //
//                                                                                    //
//  This file is part of CNull.                                                       //
//                                                                                    //
//  CNull is free software: you can redistribute it and/or modify                     //
//  it under the terms of the GNU General Public License as published by              //
//  the Free Software Foundation, either version 3 of the License, or                 //
//  (at your option) any later version.                                               //
//                                                                                    //
//  CNull is distributed in the hope that it will be useful,                          //
//  but WITHOUT ANY WARRANTY; without even the implied warranty of                    //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                     //
//  GNU General Public License for more details.                                      //
//                                                                                    //
//  You should have received a copy of the GNU General Public License                 //
//  along with CNull.  If not, see <http://www.gnu.org/licenses/>                     //
////////////////////////////////////////////////////////////////////////////////////////

#ifndef P_VALUE_HANDLER_H
#define P_VALUE_HANDLER_H

#include<Rcpp.h>
#include<vector>
#include<algorithm>

using namespace Rcpp;

namespace SamplingFunctions{

template <typename KERNEL_TYPE>
class P_value_handler
{
 public:

  typedef KERNEL_TYPE                              Kernel;
  typedef typename Kernel::Number_type             Number_type;
  typedef typename Kernel::Numeric_traits          Numeric_traits;
  typedef typename Numeric_traits::Epsilon         Epsilon;

 private:

  struct Triplet
  {
    Triplet(){}

    Number_type val;
    int rank, rank_b;
  };

  class Is_greater_triplet
  {
   public:

    template<typename Triplet>
    bool operator()( const Triplet &t1, const Triplet &t2 ) const
    { 
      if(t1.val < t2.val)
        return true;

      if(t1.val > t2.val)
        return false;      

      if(t1.rank > t2.rank) // Order is reversed here, to make sure randomized values succeed actual values
        return true;

      if(t1.rank < t2.rank)
        return false; 

      if(t1.rank_b < t2.rank_b)
        return true;

      return false; 

    } // operator()( const Triplet &t1, const Triplet &t2 )

  }; // Is_greater_triplet

 public:

  NumericVector operator()(NumericVector &original_vals, NumericVector &randomized_vals)
  {
    std::vector<Triplet> triplets;
    Number_type  eps = Epsilon()();

    for(int i=0; i<original_vals.size(); i++)
    {
      Triplet tr;

      tr.val = original_vals(i);
      tr.rank = i;
      tr.rank_b = i;
      triplets.push_back(tr);

    } // for(int i=0; i<original_vals.size(); i++)

    for(int i=0; i<randomized_vals.size(); i++)
    {
      Triplet tr;

      tr.val = randomized_vals(i)+eps;
      tr.rank = -1;
      tr.rank_b = i;
      triplets.push_back(tr);

    } // for(int i=0; i<random_vals.size(); i++)

    std::sort(triplets.begin(),triplets.end(),Is_greater_triplet());

    NumericVector outvals(original_vals.size());

    int d = randomized_vals.size();
    int rand_count=d;

    for(int i=0; i<triplets.size(); i++)
      if(triplets[i].rank == -1)
        rand_count--;
      else
        outvals(triplets[i].rank) = Number_type(rand_count+1)/Number_type(d+1);

    return outvals;

  } // void operator()(NumericVector &original_vals, NumericVector &randomized_vals)

  NumericVector slow_operator(NumericVector &original_vals, NumericVector &randomized_vals)
  {
    NumericVector  output_vals(original_vals.size());
    Number_type    eps = Epsilon()();

    for(int i=0; i<output_vals.size(); i++)
    {
      output_vals(i) = 1;  

      for(int j=0; j<randomized_vals.size(); j++)
        if(original_vals(i) <= randomized_vals(j) + eps)
          output_vals(i)++;

      output_vals(i) = Number_type(output_vals(i))/Number_type(1+randomized_vals.size());
    }

    return output_vals;   

  } // slow_operator(NumericVector &original_vals, NumericVector &randomized_vals)

}; // P_value_handler

} // namespace SamplingFunctions


#endif //P_VALUE_HANDLER_H
