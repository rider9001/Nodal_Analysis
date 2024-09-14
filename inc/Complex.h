/// ------------------------------------------
/// @file Complex.h
///
/// @brief Header file for all complex number structures and associated
/// functions
///
/// @note This file mainly used for collecting all types of complex
/// numbers and writing conversions between types.
/// ------------------------------------------
#pragma once

#include "Complex_C.h"
#include "Complex_P.h"

///--------------------------------------------------------
/// @brief Converts a polar complex number to cartesian form
///
/// @param polar polar complex number
///
/// @return cartesian complex number
Complex_C_t polarToCart(const Complex_P_t& polar);

///--------------------------------------------------------
/// @brief Converts a cartesian complex number to polar form
///
/// @param polar cartesian complex number
///
/// @return polar complex number
Complex_P_t cartToPolar(const Complex_C_t& cart);