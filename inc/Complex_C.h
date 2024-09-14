/// ------------------------------------------
/// @file Complex_C.h
///
/// @brief Header file for cartesian complex number structures and associated
/// functions
/// ------------------------------------------
#pragma once

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <string>

/// @brief Complex number structure, cartesian form
struct Complex_C_t
{
    /// @brief real component of complex
    double m_real;

    /// @brief imaginary component of complex
    double m_imagine;

    ///--------------------------------------------------------
    /// @brief Default constructor
    Complex_C_t()
    {
        m_real = 0;
        m_imagine = 0;
    };

    ///--------------------------------------------------------
    /// @brief Constructor
    ///
    /// @param real real component
    /// @param imagine imaginary component
    Complex_C_t(const double& real, const double& imagine)
    {
        m_real = real;
        m_imagine = imagine;
    }

    ///--------------------------------------------------------
    /// @brief Cast constructor using a real number
    /// @param num real number to cast to complex
    Complex_C_t(const double& real)
    {
        m_real = real;
        m_imagine = 0;
    };

    ///--------------------------------------------------------
    /// @brief Find the conjugate of complex number
    ///
    /// @return conjugate of complex
    Complex_C_t conjugate() const;

    /// @brief Find the absolute value of the input com
    ///
    /// @return absolute value of complex number
    double absolute() const;

    ///--------------------------------------------------------
    /// @brief Find the argument of the input com
    ///
    /// @return argument in radians, [-pi, pi] range
    /// clockwise relative to 1+0i
    ///
    /// @note returns 0 if com = 0+0i
    double argument() const;
};

///--------------------------------------------------------
/// @brief Overload of +, adds two complex numbers together
///
/// @param lcom left hand complex
/// @param rcom right hand complex
///
/// @return resulting complex number
Complex_C_t operator+(const Complex_C_t& lcom, const Complex_C_t& rcom);

///--------------------------------------------------------
/// @brief Overload of +, adds a complex and real number
///
/// @param lcom left hand complex
/// @param rreal right hand real
///
/// @return resulting complex number
Complex_C_t operator+(const Complex_C_t& lcom, const double& rreal);

///--------------------------------------------------------
/// @brief Overload of +, adds a real and complex number
///
/// @param lread left hand complex
/// @param rcom right hand real
///
/// @return resulting complex number
Complex_C_t operator+(const double& lreal, const Complex_C_t& rcom);

///--------------------------------------------------------
/// @brief Overload of +=, adds rcom to lcom
///
/// @param lcom left hand complex
/// @param rcom right hand complex
void operator+=(Complex_C_t& lcom, const Complex_C_t& rcom);

///--------------------------------------------------------
/// @brief Overload of +=, adds rreal to lcom
///
/// @param lcom left hand complex
/// @param rcom right hand real
void operator+=(Complex_C_t& lcom, const double& rreal);

///--------------------------------------------------------
/// @brief Overload of -, subtracts two complex numbers
///
/// @param lcom left hand complex
/// @param rcom right hand complex
///
/// @return resulting complex number
Complex_C_t operator-(const Complex_C_t& lcom, const Complex_C_t& rcom);

///--------------------------------------------------------
/// @brief Overload of -, subtracts a complex and real number
///
/// @param lcom left hand complex
/// @param rreal right hand real
///
/// @return resulting complex number
Complex_C_t operator-(const Complex_C_t& lcom, const double& rreal);

///--------------------------------------------------------
/// @brief Overload of -, subtracts real number and a complex
///
/// @param lreal left hand real
/// @param rcom right hand complex
///
/// @return resulting complex number
Complex_C_t operator-(const double& lreal, const Complex_C_t& rcom);

///--------------------------------------------------------
/// @brief Overload of -=, subtracts rcom from lcom
///
/// @param lcom left hand complex
/// @param rcom right hand complex
void operator-=(Complex_C_t& lcom, const Complex_C_t& rcom);

///--------------------------------------------------------
/// @brief Overload of +=, subtracts rreal from lcom
///
/// @param lcom left hand complex
/// @param rcom right hand real
void operator-=(Complex_C_t& lcom, const double& rreal);

///--------------------------------------------------------
/// @brief Overload of *, multiplies two complex numbers
///
/// @param lcom left hand complex
/// @param rcom right hand complex
///
/// @return resulting complex number
Complex_C_t operator*(const Complex_C_t& lcom, const Complex_C_t& rcom);

///--------------------------------------------------------
/// @brief Overload of *, multiplies a complex and real number
///
/// @param lcom left hand complex
/// @param rreal right hand real
///
/// @return resulting complex number
Complex_C_t operator*(const Complex_C_t& lcom, const double& rreal);

///--------------------------------------------------------
/// @brief Overload of *, multiplies a real number and a complex
///
/// @param lreal left hand real
/// @param rcom right hand complex
///
/// @return resulting complex number
Complex_C_t operator*(const double& lreal, const Complex_C_t& rcom);

///--------------------------------------------------------
/// @brief Overload of *=, multiplies lcom by rcom
///
/// @param lcom left hand complex
/// @param rcom right hand complex
void operator*=(Complex_C_t& lcom, const Complex_C_t& rcom);

///--------------------------------------------------------
/// @brief Overload of *=, multiplies lcom by rreal
///
/// @param lcom left hand complex
/// @param rreal right hand real
void operator*=(Complex_C_t& lcom, const double& rreal);

///--------------------------------------------------------
/// @brief Overload of /, divides two complex numbers
///
/// @param lcom left hand complex
/// @param rcom right hand complex
///
/// @return resulting complex number
Complex_C_t operator/(const Complex_C_t& lcom, const Complex_C_t& rcom);

///--------------------------------------------------------
/// @brief Overload of /, divides a complex by a real
///
/// @param lcom left hand complex
/// @param rreal right hand real
///
/// @return resulting complex number
Complex_C_t operator/(const Complex_C_t& lcom, const double& rreal);

///--------------------------------------------------------
/// @brief Overload of /, divides a real number by a complex
///
/// @param lreal left hand real
/// @param rcom right hand complex
///
/// @return resulting complex number
Complex_C_t operator/(const double& lreal, const Complex_C_t& rcom);

///--------------------------------------------------------
/// @brief Overload of /=, divides lcom by rcom
///
/// @param lcom left hand complex
/// @param rcom right hand complex
void operator/=(Complex_C_t& lcom, const Complex_C_t& rcom);

///--------------------------------------------------------
/// @brief Overload of /=, divides lcom by rreal
///
/// @param lcom left hand complex
/// @param rreal right hand real
void operator/=(Complex_C_t& lcom, const double& rreal);

///--------------------------------------------------------
/// @brief Overload of ==, are complex numbers equal?
///
/// @param lcom left hand complex
/// @param rcom right hand complex
///
/// @return equality boolean
bool operator==(const Complex_C_t& lcom, const Complex_C_t& rcom);

///--------------------------------------------------------
/// @brief Overload of ==, are complex and real equal?
///
/// @param lcom left hand complex
/// @param rreal right hand real
///
/// @return equality boolean
bool operator==(const Complex_C_t& lcom, const double& rreal);

///--------------------------------------------------------
/// @brief Overload of ==, is a complex and a real equal?
///
/// @param lreal left hand real
/// @param rcom right hand complex
///
/// @return equality boolean
bool operator==(const double& lreal, const Complex_C_t& rcom);

///--------------------------------------------------------
/// @brief Overload of !=, are complex numbers unequal?
///
/// @param lcom left hand complex
/// @param rcom right hand complex
///
/// @return inequality boolean
bool operator!=(const Complex_C_t& lcom, const Complex_C_t& rcom);

///--------------------------------------------------------
/// @brief Overload of !=, is a complex and a real unequal?
///
/// @param lcom left hand complex
/// @param rreal right hand real
///
/// @return equality boolean
bool operator!=(const Complex_C_t& lcom, const double& rreal);

///--------------------------------------------------------
/// @brief Overload of !=, is a complex and a real unequal?
///
/// @param lreal left hand real
/// @param rcom right hand complex
///
/// @return equality boolean
bool operator!=(const double& lreal, const Complex_C_t& rcom);

///--------------------------------------------------------
/// @brief Overload of <<, used to output string format of complex
///
/// @param os output stream
/// @param com complex number to string convert
///
/// @return output stream
std::ostream& operator<<(std::ostream& os, const Complex_C_t& com);

///--------------------------------------------------------
/// @brief Raises eulers number by a complex value
///
/// @param com complex to raise e by
///
/// @return complex result of e^com
Complex_C_t raiseEComplex(const Complex_C_t& com);

///--------------------------------------------------------
/// @brief Raises a complex number by another complex number
///
/// @param base complex number to raise
/// @param raise power to raise the base to
///
/// @return resulting complex number
Complex_C_t powComplex(const Complex_C_t& base, const Complex_C_t& raise);
