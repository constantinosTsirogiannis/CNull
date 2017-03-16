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

#ifndef RANDOM_VARIATE_GENERATORS_H
#define RANDOM_VARIATE_GENERATORS_H

#include"Numeric_traits_double.h"
#include<Rcpp.h>
#include<vector>
#include<random>
#include<chrono>
#include<map>

namespace SamplingFunctions{

template <typename KERNEL_TYPE>
class Integer_generator
{
 public:

  typedef KERNEL_TYPE                    Kernel;
  typedef typename Kernel::Number_type  Number_type;
  
 public:

  Integer_generator():_generator(std::chrono::system_clock::now().time_since_epoch().count()){}

  Integer_generator(unsigned int seed):_generator(seed){}

  template< typename OutputIterator>
  void operator()(int max, int num, OutputIterator ot)
  {
    std::uniform_int_distribution<int> distribution(0,max-1);    

    for (int i=0; i<num; ++i)  
      *ot++ = distribution(_generator);

  } // void operator()(...)


  template< typename OutputIterator>
  void two_draws_operator(int max, int num, OutputIterator ot)
  {
    std::uniform_int_distribution<int> distribution(0,max-1);

    for(int i=0; i<num; i++)
    {
      int n_a = distribution(_generator),
          n_b = distribution(_generator);

      while(n_b == n_a)
        n_b = distribution(_generator);

      *ot++ = std::make_pair(n_a,n_b);
    }

  } // void two_draws_operator(...)

  std::default_random_engine _generator;
  
}; // class Integer_generator


template <typename KERNEL_TYPE>
class Binomial_generator
{

 public:

  typedef KERNEL_TYPE                    Kernel;
  typedef typename Kernel::Number_type  Number_type;
  

 public:

  Binomial_generator(){}

  template< typename OutputIterator>
  void operator()(int max, Number_type p, int num, OutputIterator ot)
  {
    Rcpp::NumericVector first_draw = Rcpp::rbinom(num,max,p);

    for (int i=0; i<first_draw.size(); ++i)  
      *ot++ = int(first_draw(i));

  } // void operator()(...)


  template< typename OutputIterator>
  void two_draws_operator(int max, Number_type p, int num, OutputIterator ot)
  {
    Rcpp::NumericVector first_draw = Rcpp::rbinom(num,max,Number_type(2.0)*p);

    for(int i=0; i<first_draw.size(); i++)
    {
      Rcpp::NumericVector second_draw =  Rcpp::rbinom(1,first_draw(i),Number_type(0.5));
      *ot++ = std::make_pair(second_draw(0), first_draw(i)-second_draw(0));
    }

  } // void two_draws_operator(...)


  template< typename OutputIterator>
  void two_draws_operator_version_II(int max, Number_type p, int num, OutputIterator ot)
  {
    Rcpp::NumericVector first_draw = Rcpp::rbinom(num,max,Number_type(2.0)*p);

    std::map<int,int> draws;

    for(int i=0; i<first_draw.size(); i++)
      if( draws.find(first_draw(i)) == draws.end() )
        draws[first_draw(i)] = 1;
      else
        draws[first_draw(i)]++;

    for(typename std::map<int,int>::iterator it=draws.begin(); it != draws.end(); it++)
    {
      int sum = it->first,
          count = it->second;

      Rcpp::NumericVector second_draw =  Rcpp::rbinom(count,sum,Number_type(0.5));

      for(int i=0; i<second_draw.size(); i++)
        *ot++ = std::make_pair(second_draw(i),sum-second_draw(i));
    }

  } // void two_draws_operator_version_II(...)
  
}; // class Binomial_generator

} // namespace SamplingFunctions

#endif // RANDOM_VARIATE_GENERATOR_H
