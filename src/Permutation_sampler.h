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

#ifndef PERMUTATION_SAMPLER_H
#define PERMUTATION_SAMPLER_H

#include<vector>
#include<Rcpp.h>

using namespace Rcpp;

namespace SamplingFunctions{

template <typename KERNEL_TYPE>
class Permutation_sampler
{

 public:

  typedef KERNEL_TYPE                          Kernel;
  typedef typename Kernel::Number_type         Number_type;
  typedef typename Kernel::Integer_generator  Integer_generator;

 public:

  void alpha_diversity_operator(NumericMatrix &in_m, NumericMatrix &out_m)
  {
    int repetitions = out_m.nrow();

    Integer_generator gen;

    for(int j=0; j<in_m.ncol(); j++)    
    {  
      std::vector<int> random_ints;

      gen(in_m.nrow(),repetitions,std::back_inserter(random_ints));

      for(int i=0; i<out_m.nrow(); i++)
        out_m(i,j) = in_m(random_ints[i],j);

    } // for(int j=0; j<in_m.ncol(); j++) 

  } // void alpha_diversity_operator(NumericMatrix &in_m, NumericMatrix &out_m )


  void beta_diversity_operator( NumericMatrix &in_m, 
                                NumericMatrix &out_a, 
                                NumericMatrix &out_b)
  {
    Integer_generator gen;

    int repetitions = out_a.nrow();

    for(int j=0; j<in_m.ncol(); j++)    
    {  
      std::vector< std::pair<int,int> > random_ints;

      gen.two_draws_operator(in_m.nrow(),repetitions,std::back_inserter(random_ints));

      for(int i=0; i<repetitions; i++)
      {
        out_a(i,j) = in_m(random_ints[i].first,j);
        out_b(i,j) = in_m(random_ints[i].second,j);
      }

    } // for(int j=0; j<in_m.ncol(); j++) 
    
  } // void beta_diversity_operator(NumericMatrix &in_m, NumericMatrix &out_a, NumericMatrix &out_b )


  void beta_diversity_operator_interleaved( NumericMatrix &in_m, 
                                            NumericMatrix &out_m)
  {
    Integer_generator gen;

    int repetitions = out_m.nrow()/2;

    for(int j=0; j<in_m.ncol(); j++)    
    {  
      std::vector< std::pair<int,int> > random_ints;

      gen.two_draws_operator(in_m.nrow(),repetitions,std::back_inserter(random_ints));

      for(int i=0; i<repetitions; i++)
      {
        out_m(2*i,j) = in_m(random_ints[i].first,j);
        out_m((2*i)+1,j) = in_m(random_ints[i].second,j);
      }

    } // for(int j=0; j<in_m.ncol(); j++) 
    
  } // void beta_diversity_operator_interleaved(NumericMatrix &in_m, NumericMatrix &out_m )

}; // Permutation_sampler

} // namespace SamplingFunctions


#endif //PERMUTATION_SAMPLER_H
