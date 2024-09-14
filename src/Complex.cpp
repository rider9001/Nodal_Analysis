/// ------------------------------------------
/// @file Complex.cpp
///
/// @brief Source file for associated conversion functions
///
/// @note This file mainly used for collecting all types of complex
/// numbers and writing conversions between types.
/// ------------------------------------------

#include "../inc/Complex.h"

///--------------------------------------------------------
Complex_C_t polarToCart(const Complex_P_t& polar)
{
    return Complex_C_t{polar.real(), polar.imaginary()};
}

///--------------------------------------------------------
Complex_P_t cartToPolar(const Complex_C_t& cart)
{
    return Complex_P_t{cart.absolute(), cart.argument()};
}