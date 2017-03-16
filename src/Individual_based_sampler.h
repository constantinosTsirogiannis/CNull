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

#ifndef INDIVIDUAL_BASED_SAMPLER_H
#define INDIVIDUAL_BASED_SAMPLER_H

#include<vector>
#include<Rcpp.h>

using namespace Rcpp;

namespace SamplingFunctions{

template <typename KERNEL_TYPE>
class Individual_based_sampler
{

 public:

  typedef KERNEL_TYPE                           Kernel;
  typedef typename Kernel::Number_type         Number_type;
  typedef typename Kernel::Binomial_generator  Binomial_generator;

 public:

  void alpha_diversity_operator(NumericMatrix &in_m, NumericMatrix &out_m)
  {
    std::vector<int>  col_sums;

    col_sums.assign(in_m.ncol(), 0);

    for(int j=0; j<in_m.ncol(); j++)
      for(int i=0; i<in_m.nrow(); i++)
        col_sums[j]+= in_m(i,j);

    Number_type p = Number_type(1.0)/Number_type(in_m.nrow());

    int repetitions = out_m.nrow();

    Binomial_generator gen;


    for(int j=0; j<in_m.ncol(); j++)    
    {  
      std::vector<int> random_ints;

      gen(col_sums[j],p,repetitions,std::back_inserter(random_ints));

      for(int i=0; i<out_m.nrow(); i++)
        out_m(i,j) = random_ints[i];

    } // for(int j=0; j<in_m.ncol(); j++) 

  } // void alpha_diversity_operator(NumericMatrix &in_m, NumericMatrix &out_m )


  void beta_diversity_operator( NumericMatrix &in_m, 
                                 NumericMatrix &out_a, 
                                 NumericMatrix &out_b)
  {
    std::vector<int>  col_sums;

    col_sums.assign(in_m.ncol(), 0);

    for(int j=0; j<in_m.ncol(); j++)
      for(int i=0; i<in_m.nrow(); i++)
        col_sums[j]+= in_m(i,j);

    Number_type p = Number_type(1.0)/Number_type(in_m.nrow());

    Binomial_generator gen;

    int repetitions = out_a.nrow();

    for(int j=0; j<in_m.ncol(); j++)    
    {  
      std::vector< std::pair<int,int> > random_ints;

      gen.two_draws_operator(col_sums[j],p,repetitions,std::back_inserter(random_ints));

      for(int i=0; i<repetitions; i++)
      {
        out_a(i,j) = random_ints[i].first;
        out_b(i,j) = random_ints[i].second;
      }

    } // for(int j=0; j<in_m.ncol(); j++) 
    
  } // void beta_diversity_operator(NumericMatrix &in_m, NumericMatrix &out_a, NumericMatrix &out_b )


  void beta_diversity_operator_interleaved( NumericMatrix & in_m, 
                                            NumericMatrix & out_m)
  {
    std::vector<int>  col_sums;

    col_sums.assign(in_m.ncol(), 0);

    for(int j=0; j<in_m.ncol(); j++)
      for(int i=0; i<in_m.nrow(); i++)
        col_sums[j]+= in_m(i,j);

    Number_type p = Number_type(1.0)/Number_type(in_m.nrow());

    Binomial_generator gen;

    int repetitions = out_m.nrow()/2;

    for(int j=0; j<in_m.ncol(); j++)    
    {  
      std::vector< std::pair<int,int> > random_ints;

      gen.two_draws_operator(col_sums[j],p,repetitions,std::back_inserter(random_ints));

      for(int i=0; i<repetitions; i++)
      {
        out_m(2*i,j) = random_ints[i].first;
        out_m((2*i)+1,j) = random_ints[i].second;
      }

    } // for(int j=0; j<in_m.ncol(); j++) 

  } // void beta_diversity_operator_interleaved(NumericMatrix &in_m, NumericMatrix &out_m )

}; // Individual_based_sampler

} // namespace SamplingFunctions


#endif //INDIVIDUAL_BASED_SAMPLER_H
