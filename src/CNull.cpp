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

#include<Rcpp.h>
#include<R_ext/Print.h>
#include<R.h>
#include<vector>
#include<limits>
#include"Numeric_traits_double.h"
#include"Sampling_functions_kernel.h"

typedef SamplingFunctions::Numeric_traits_double   Numeric_traits;
typedef Sampling_functions_kernel<Numeric_traits>  Kernel;
typedef Kernel::Number_type                        Number_type;
typedef Kernel::Individual_based_sampler           Individual_based_sampler;
typedef Kernel::Permutation_sampler                Permutation_sampler;
typedef Kernel::P_value_handler                    P_value_handler;

using namespace Rcpp;

bool check_NA_and_negative_values(NumericMatrix &in_m, bool allow_negative_values)
{
  bool contains_negative = false;

  for(int i=0; i<in_m.nrow(); i++)
    for(int j=0; j<in_m.ncol(); j++)
      if( in_m(i,j) == NA_REAL || in_m(i,j) == NA_INTEGER || std::isnan(in_m(i,j)) || in_m(i,j) < 0 )
      {
        if(in_m(i,j) == NA_REAL || in_m(i,j) == NA_INTEGER || std::isnan(in_m(i,j)))
        {
          Rcpp::Rcerr << std::endl << "Error: the input matrix contains NA values ... aborting." << std::endl; 

          return false;

        } // if(in_m(i,j) == NA_REAL)
        else if( allow_negative_values == false)
        {
          Rcpp::Rcerr << std::endl << "Error: the input matrix contains negative values ... aborting." << std::endl; 

          return false; 

        } // else if( allow_negative_values == false)
        else
          contains_negative = true;

      } // if( in_m(i,j == NA_REAL || in_m(i,j) < 0 )

 if(contains_negative == true)
   Rcpp::Rcerr << std::endl << "Warning: the input matrix contains negative values." << std::endl; 

 return true;
     
} // check_NA_and_negative_values(...)

// [[Rcpp::export]]
NumericMatrix  communities_permutation_sampling_alpha(NumericMatrix in_m, int repetitions)
{
  bool check = check_NA_and_negative_values(in_m, true);

  if(check==false)
  {
    NumericMatrix out_m(0,0);
    return out_m; 
  }

  Permutation_sampler ps;

  NumericMatrix out_m(repetitions,in_m.ncol());

  ps.alpha_diversity_operator(in_m,out_m);

  return out_m; 
}


// [[Rcpp::export]]
Rcpp::List communities_permutation_sampling_beta(NumericMatrix in_m, int repetitions)
{
  bool check = check_NA_and_negative_values(in_m, true);

  if(check==false)
  {
    NumericMatrix out_m(0,0);

    return Rcpp::List::create( Rcpp::Named("samples.a") = out_m,
                               Rcpp::Named("samples.b") = out_m); 
  }

  Permutation_sampler ps;

  NumericMatrix out_a(repetitions,in_m.ncol()),
                out_b(repetitions,in_m.ncol());

 ps.beta_diversity_operator(in_m,out_a, out_b);

  return Rcpp::List::create( Rcpp::Named("samples.a") = out_a,
                              Rcpp::Named("samples.b") = out_b);  
}


// [[Rcpp::export]]
NumericMatrix communities_permutation_sampling_beta_interleaved_matrices(NumericMatrix in_m, int repetitions)
{
  bool check = check_NA_and_negative_values(in_m, true);

  if(check==false)
  {
    NumericMatrix out_m(0,0);
    return out_m; 
  }

  Permutation_sampler ps;
  NumericMatrix out_m(2*repetitions,in_m.ncol());
  ps.beta_diversity_operator_interleaved(in_m,out_m);

  return out_m;
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
NumericMatrix  communities_individual_based_sampling_alpha(NumericMatrix in_m, int repetitions)
{
  bool check = check_NA_and_negative_values(in_m, false);

  if(check==false)
  {
    NumericMatrix out_m(0,0);
    return out_m; 
  }

  Individual_based_sampler ubs;

  NumericMatrix out_m(repetitions,in_m.ncol());

  ubs.alpha_diversity_operator(in_m,out_m);

  return out_m; 
}

// [[Rcpp::export]]
NumericVector  compute_pvalues(NumericVector original_vals, NumericVector randomised_vals)
{
  P_value_handler pvh;
  
  NumericVector out_v = pvh(original_vals,randomised_vals);

  return out_v; 
}

// [[Rcpp::export]]
Rcpp::List communities_individual_based_sampling_beta(NumericMatrix in_m, int repetitions)
{
  bool check = check_NA_and_negative_values(in_m, false);

  if(check==false)
  {
    NumericMatrix out_m(0,0);

    return Rcpp::List::create( Rcpp::Named("samples.a") = out_m,
                               Rcpp::Named("samples.b") = out_m); 
  }

  Individual_based_sampler ubs;

  NumericMatrix out_a(repetitions,in_m.ncol()),
                out_b(repetitions,in_m.ncol());

  ubs.beta_diversity_operator(in_m,out_a, out_b);

  return Rcpp::List::create( Rcpp::Named("samples.a") = out_a,
                              Rcpp::Named("samples.b") = out_b);  
}


// [[Rcpp::export]]
NumericMatrix communities_individual_based_sampling_beta_interleaved_matrices(NumericMatrix in_m, int repetitions)
{
  bool check = check_NA_and_negative_values(in_m, true);

  if(check==false)
  {
    NumericMatrix out_m(0,0);
    return out_m; 
  }

  Individual_based_sampler ubs;
  NumericMatrix out_m(2*repetitions,in_m.ncol());
  ubs.beta_diversity_operator_interleaved(in_m,out_m);

  return out_m;
}
