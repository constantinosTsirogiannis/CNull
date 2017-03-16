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

#ifndef NUMERIC_TRAITS_DOUBLE_H
#define NUMERIC_TRAITS_DOUBLE_H

#include<cmath>
#include "Numeric_traits_types/Protected_number_type.h" 

namespace SamplingFunctions 
{
  struct Numeric_traits_double
  {
    typedef double                                                      Number_type;
    typedef Numeric_traits_double                                        Self;

    class Is_exact
    {
      public:
 
        bool operator()(void)
        { return false;}
    };

    class To_double
    {
      public:

        double operator()(double x)
        { return x; }
    };

    class Power
    {
      public:

        double operator()(double x, int k)
        { return std::pow(x,k); }
    };

    class Ceiling
    {
      public:

        double operator()(double x)
        { return ceil(x); }
    };

    class Square_root
    {
      public:

        double operator()(double x)
        { return std::sqrt(x); }
    };

    class Absolute_value
    {
      public:
 
        double operator()(double x)
        { return std::abs(x);}
    };

    class Cosine
    {
      public:
 
        double operator()(double x)
        { return std::cos(x);}
    };

    class Sine
    {
      public:
 
        double operator()(double x)
        { return std::sin(x);}
    };

    class Epsilon
    {
      public:
 
        double operator()(void)
        { return double(0.01);}
    };


    /////////////////////////////////////////////////////////////////////
    // Type used for preserving accuracy in polynomial multiplications //
    /////////////////////////////////////////////////////////////////////

    typedef SamplingFunctions::Protected_number_type<Self>  Protected_number_type;

  }; // struct Numeric_traits_double

} //namespace SamplingFunctions

#endif // NUMERIC_TRAITS_DOUBLE_H
