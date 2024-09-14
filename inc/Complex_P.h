/// ------------------------------------------
/// @file Complex_P.h
///
/// @brief Header file for polar complex number structures and associated
/// functions
/// ------------------------------------------
#pragma once

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <string>

#include "Complex_C.h"

/// @brief Complex number structure, polar form
struct Complex_P_t
{
    /// @brief magnitude of polar complex
    double m_mag = 0;

    /// @brief argument/angle of polar complex, radians, should only ever be in the range [0, 2pi]
    double m_arg = 0;

    /// @brief Default constructor
    Complex_P_t()
    {
        m_mag = 0;
        m_arg = 0;
    };

    ///--------------------------------------------------------
    /// @brief Constructor
    ///
    /// @param mag magnitude of polar complex
    /// @param arg argument of polar complex
    Complex_P_t(const double& mag, const double& arg)
    {
        m_mag = mag;
        setArg(arg);
    };

    ///--------------------------------------------------------
    /// @brief Constructor, allows casting from real number
    ///
    /// @param mag magnitude
    Complex_P_t(const double& mag)
    {
        m_mag = mag;
        m_arg = 0;
    };

    ///--------------------------------------------------------
    /// @brief Setter function for the argument that bounds argument
    /// to [-pi, pi]
    ///
    /// @param arg to set m_arg to
    void setArg(const double& arg);

    ///--------------------------------------------------------
    /// @brief gets the real component of the polar complex
    ///
    /// @returns real component of the polar complex
    double real() const;

    ///--------------------------------------------------------
    /// @brief gets the imaginary component of the polar complex
    ///
    /// @returns imaginary component of the polar complex
    double imaginary() const;
};

///--------------------------------------------------------
/// @brief Overload of +, adds two complex numbers together
///
/// @param lcom left hand complex
/// @param rcom right hand complex
///
/// @return resulting complex number
Complex_P_t operator+(const Complex_P_t& lcom, const Complex_P_t& rcom);

///--------------------------------------------------------
/// @brief Overload of +, adds a complex and real number
///
/// @param lcom left hand complex
/// @param rreal right hand real
///
/// @return resulting complex number
Complex_P_t operator+(const Complex_P_t& lcom, const double& rreal);

///--------------------------------------------------------
/// @brief Overload of +, adds a real and complex number
///
/// @param lread left hand complex
/// @param rcom right hand real
///
/// @return resulting complex number
Complex_P_t operator+(const double& lreal, const Complex_P_t& rcom);

///--------------------------------------------------------
/// @brief Overload of +=, adds rcom to lcom
///
/// @param lcom left hand complex
/// @param rcom right hand complex
void operator+=(Complex_P_t& lcom, const Complex_P_t& rcom);

///--------------------------------------------------------
/// @brief Overload of +=, adds rreal to lcom
///
/// @param lcom left hand complex
/// @param rcom right hand real
void operator+=(Complex_P_t& lcom, const double& rreal);

///--------------------------------------------------------
/// @brief Overload of -, subtracts two complex numbers
///
/// @param lcom left hand complex
/// @param rcom right hand complex
///
/// @return resulting complex number
Complex_P_t operator-(const Complex_P_t& lcom, const Complex_P_t& rcom);

///--------------------------------------------------------
/// @brief Overload of -, subtracts a complex and real number
///
/// @param lcom left hand complex
/// @param rreal right hand real
///
/// @return resulting complex number
Complex_P_t operator-(const Complex_P_t& lcom, const double& rreal);

///--------------------------------------------------------
/// @brief Overload of -, subtracts real number and a complex
///
/// @param lreal left hand real
/// @param rcom right hand complex
///
/// @return resulting complex number
Complex_P_t operator-(const double& lreal, const Complex_P_t& rcom);

///--------------------------------------------------------
/// @brief Overload of -=, subtracts rcom from lcom
///
/// @param lcom left hand complex
/// @param rcom right hand complex
void operator-=(Complex_P_t& lcom, const Complex_P_t& rcom);

///--------------------------------------------------------
/// @brief Overload of +=, subtracts rreal from lcom
///
/// @param lcom left hand complex
/// @param rcom right hand real
void operator-=(Complex_P_t& lcom, const double& rreal);

///--------------------------------------------------------
/// @brief Overload of *, multiplies two complex numbers
///
/// @param lcom left hand complex
/// @param rcom right hand complex
///
/// @return resulting complex number
Complex_P_t operator*(const Complex_P_t& lcom, const Complex_P_t& rcom);

///--------------------------------------------------------
/// @brief Overload of *, multiplies a complex and real number
///
/// @param lcom left hand complex
/// @param rreal right hand real
///
/// @return resulting complex number
Complex_P_t operator*(const Complex_P_t& lcom, const double& rreal);

///--------------------------------------------------------
/// @brief Overload of *, multiplies a real number and a complex
///
/// @param lreal left hand real
/// @param rcom right hand complex
///
/// @return resulting complex number
Complex_P_t operator*(const double& lreal, const Complex_P_t& rcom);

///--------------------------------------------------------
/// @brief Overload of *=, multiplies lcom by rcom
///
/// @param lcom left hand complex
/// @param rcom right hand complex
void operator*=(Complex_P_t& lcom, const Complex_P_t& rcom);

///--------------------------------------------------------
/// @brief Overload of *=, multiplies lcom by rreal
///
/// @param lcom left hand complex
/// @param rreal right hand real
void operator*=(Complex_P_t& lcom, const double& rreal);

///--------------------------------------------------------
/// @brief Overload of /, divides two complex numbers
///
/// @param lcom left hand complex
/// @param rcom right hand complex
///
/// @return resulting complex number
Complex_P_t operator/(const Complex_P_t& lcom, const Complex_P_t& rcom);

///--------------------------------------------------------
/// @brief Overload of /, divides a complex by a real
///
/// @param lcom left hand complex
/// @param rreal right hand real
///
/// @return resulting complex number
Complex_P_t operator/(const Complex_P_t& lcom, const double& rreal);

///--------------------------------------------------------
/// @brief Overload of /, divides a real number by a complex
///
/// @param lreal left hand real
/// @param rcom right hand complex
///
/// @return resulting complex number
Complex_P_t operator/(const double& lreal, const Complex_P_t& rcom);

///--------------------------------------------------------
/// @brief Overload of /=, divides lcom by rcom
///
/// @param lcom left hand complex
/// @param rcom right hand complex
void operator/=(Complex_P_t& lcom, const Complex_P_t& rcom);

///--------------------------------------------------------
/// @brief Overload of /=, divides lcom by rreal
///
/// @param lcom left hand complex
/// @param rreal right hand real
void operator/=(Complex_P_t& lcom, const double& rreal);

///--------------------------------------------------------
/// @brief Overload of ==, are complex numbers equal?
///
/// @param lcom left hand complex
/// @param rcom right hand complex
///
/// @return equality boolean
bool operator==(const Complex_P_t& lcom, const Complex_P_t& rcom);

///--------------------------------------------------------
/// @brief Overload of ==, are complex and real equal?
///
/// @param lcom left hand complex
/// @param rreal right hand real
///
/// @return equality boolean
bool operator==(const Complex_P_t& lcom, const double& rreal);

///--------------------------------------------------------
/// @brief Overload of ==, is a complex and a real equal?
///
/// @param lreal left hand real
/// @param rcom right hand complex
///
/// @return equality boolean
bool operator==(const double& lreal, const Complex_P_t& rcom);

///--------------------------------------------------------
/// @brief Overload of !=, are complex numbers unequal?
///
/// @param lcom left hand complex
/// @param rcom right hand complex
///
/// @return inequality boolean
bool operator!=(const Complex_P_t& lcom, const Complex_P_t& rcom);

///--------------------------------------------------------
/// @brief Overload of !=, is a complex and a real unequal?
///
/// @param lcom left hand complex
/// @param rreal right hand real
///
/// @return equality boolean
bool operator!=(const Complex_P_t& lcom, const double& rreal);

///--------------------------------------------------------
/// @brief Overload of !=, is a complex and a real unequal?
///
/// @param lreal left hand real
/// @param rcom right hand complex
///
/// @return equality boolean
bool operator!=(const double& lreal, const Complex_P_t& rcom);

///--------------------------------------------------------
/// @brief Overload of <<, used to output string format of complex
///
/// @param os output stream
/// @param com complex number to string convert
///
/// @return output stream
std::ostream& operator<<(std::ostream& os, const Complex_P_t& com);