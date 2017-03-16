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

#ifndef SAMPLING_FUNCTIONS_KERNEL_H
#define SAMPLING_FUNCTIONS_KERNEL_H

#include "Numeric_traits_double.h"
#include "Random_variate_generators.h"
#include "Individual_based_sampler.h"
#include "Permutation_sampler.h"
#include "P_value_handler.h"

template< class NTS = typename SamplingFunctions::Numeric_traits_double>
class Sampling_functions_kernel
{
  public:

  //////////////////////////////
  // Kernel and numeric types //
  //////////////////////////////

  typedef NTS                                                    Numeric_traits;
  typedef Sampling_functions_kernel<Numeric_traits>              Self;
  typedef typename Numeric_traits::Number_type                   Number_type;

  //////////////////////////////
  // Random number generators //
  //////////////////////////////

  typedef typename SamplingFunctions::Binomial_generator<Self>   Binomial_generator;
  typedef typename SamplingFunctions::Integer_generator<Self>    Integer_generator;

  ///////////////////////////
  // Main sampling classes //
  ///////////////////////////

  typedef typename SamplingFunctions::Individual_based_sampler<Self>  Individual_based_sampler;
  typedef typename SamplingFunctions::Permutation_sampler<Self>       Permutation_sampler;

  /////////////////////////////////////////////
  // Class for efficient p-value computation //
  /////////////////////////////////////////////

  typedef typename SamplingFunctions::P_value_handler<Self>  P_value_handler; 

}; // class Sampling_functions_kernel

#endif //SAMPLING_FUNCTIONS_KERNEL_H
