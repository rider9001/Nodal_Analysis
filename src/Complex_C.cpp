/// ------------------------------------------
/// @file Complex_C.cpp
///
/// @brief Soource file for Complex_C_t number structure
/// ------------------------------------------

#include "../inc/Complex_C.h"

///--------------------------------------------------------
Complex_C_t operator+(const Complex_C_t& lcom, const Complex_C_t& rcom)
{
    return Complex_C_t{lcom.m_real + rcom.m_real, lcom.m_imagine + rcom.m_imagine};
}

///--------------------------------------------------------
Complex_C_t operator+(const Complex_C_t& lcom, const double& rreal)
{
    return Complex_C_t{lcom.m_real + rreal, lcom.m_imagine};
}

///--------------------------------------------------------
Complex_C_t operator+(const double& lreal, const Complex_C_t& rcom)
{
    // + is commutative so use the other arragement
    return rcom + lreal;
}

///--------------------------------------------------------
void operator+=(Complex_C_t& lcom, Complex_C_t const& rcom)
{
    lcom = lcom + rcom;
}

///--------------------------------------------------------
void operator+=(Complex_C_t& lcom, const double& rreal)
{
    lcom = lcom + rreal;
}

///--------------------------------------------------------
Complex_C_t operator-(const Complex_C_t& lcom, const Complex_C_t& rcom)
{
    return Complex_C_t{lcom.m_real - rcom.m_real, lcom.m_imagine - rcom.m_imagine};
}

///--------------------------------------------------------
Complex_C_t operator-(const Complex_C_t& lcom, const double& rreal)
{
    return Complex_C_t{lcom.m_real - rreal, lcom.m_imagine};
}

///--------------------------------------------------------
Complex_C_t operator-(const double& lreal, const Complex_C_t& rcom)
{
    return Complex_C_t{lreal - rcom.m_real, -rcom.m_imagine};
}

///--------------------------------------------------------
void operator-=(Complex_C_t& lcom, const Complex_C_t& rcom)
{
    lcom = lcom - rcom;
}

///--------------------------------------------------------
void operator-=(Complex_C_t& lcom, const double& rreal)
{
    lcom = lcom - rreal;
}

///--------------------------------------------------------
Complex_C_t operator*(const Complex_C_t& lcom, const Complex_C_t& rcom)
{
    return Complex_C_t{
        lcom.m_real * rcom.m_real - lcom.m_imagine * rcom.m_imagine,
        lcom.m_real * rcom.m_imagine + lcom.m_imagine * rcom.m_real
    };
}

///--------------------------------------------------------
Complex_C_t operator*(const Complex_C_t& lcom, const double& rreal)
{
    return Complex_C_t{lcom.m_real * rreal, lcom.m_imagine * rreal};
}

///--------------------------------------------------------
Complex_C_t operator*(const double& lreal, const Complex_C_t& rcom)
{
    // * is commutative so use the other arragement
    return rcom * lreal;
}

///--------------------------------------------------------
void operator*=(Complex_C_t& lcom, const Complex_C_t& rcom)
{
    lcom = lcom * rcom;
}

///--------------------------------------------------------
void operator*=(Complex_C_t& lcom, const double& rreal)
{
    lcom = lcom * rreal;
}

///--------------------------------------------------------
Complex_C_t operator/(const Complex_C_t& lcom, const Complex_C_t& rcom)
{
    double div = pow(rcom.m_real, 2) + pow(rcom.m_imagine, 2);
    return Complex_C_t{
        (lcom.m_real * rcom.m_real + lcom.m_imagine * rcom.m_imagine) / div,
        (lcom.m_imagine * rcom.m_real - lcom.m_real * rcom.m_imagine) / div
    };
}

///--------------------------------------------------------
Complex_C_t operator/(const Complex_C_t& lcom, const double& rreal)
{
    return Complex_C_t{lcom.m_real / rreal, lcom.m_imagine / rreal};
}

///--------------------------------------------------------
Complex_C_t operator/(const double& lreal, const Complex_C_t& rcom)
{
    double div = pow(rcom.m_real, 2) + pow(rcom.m_imagine, 2);
    return Complex_C_t{
        (lreal * rcom.m_real) / div,
        (-lreal * rcom.m_imagine) / div
    };
}

///--------------------------------------------------------
void operator/=(Complex_C_t& lcom, const Complex_C_t& rcom)
{
    lcom = lcom / rcom;
}

///--------------------------------------------------------
void operator/=(Complex_C_t& lcom, const double& rreal)
{
    lcom = lcom / rreal;
}

///--------------------------------------------------------
bool operator==(const Complex_C_t& lcom, const Complex_C_t& rcom)
{
    return (lcom.m_real == rcom.m_real) and (lcom.m_imagine == rcom.m_imagine);
}

///--------------------------------------------------------
bool operator==(const Complex_C_t& lcom, const double& rreal)
{
    return (lcom.m_real == rreal) and (lcom.m_imagine == 0);
}

///--------------------------------------------------------
bool operator==(const double& lreal, const Complex_C_t& rcom)
{
    return rcom == lreal;
}

///--------------------------------------------------------
bool operator!=(const Complex_C_t& lcom, const Complex_C_t& rcom)
{
    return !(lcom == rcom);
}

///--------------------------------------------------------
bool operator!=(const Complex_C_t& lcom, const double& rreal)
{
    return !(lcom == rreal);
}

///--------------------------------------------------------
bool operator!=(const double& lreal, const Complex_C_t& rcom)
{
    return !(lreal == rcom);
}

///--------------------------------------------------------
std::ostream& operator<<(std::ostream& os, Complex_C_t const& com)
{
    if (com.m_real < 0)
    {
        os << "-";
    }
    else
    {
        os << "+";
    }
    os << std::to_string(fabs(com.m_real));

    if (com.m_imagine < 0)
    {
        os << "-";
    }
    else
    {
        os << "+";
    }

    os << std::to_string(fabs(com.m_imagine)) << "i";
    return os;
}

///--------------------------------------------------------
Complex_C_t Complex_C_t::conjugate() const
{
    return Complex_C_t{m_real, -m_imagine};
}

///--------------------------------------------------------
double Complex_C_t::absolute() const
{
    return sqrt(pow(m_real, 2) + pow(m_imagine, 2));
}

///--------------------------------------------------------
double Complex_C_t::argument() const
{
    // Implemenation of atan that allows for +/- inf inputs
    // and scales output to [-pi, pi] range, relative to eastward 0 deg
    if (m_real == 0 and m_imagine == 0)
    {
        return 0;
    }

    double preRotation = 0;
    if (m_real < 0)
    {
        if (m_imagine > 0)
        {
            preRotation = M_PI_2;
        }
        else
        {
            preRotation = -M_PI_2;
        }
    }

    return atan(m_imagine / m_real) + preRotation;
}

///--------------------------------------------------------
Complex_C_t raiseEComplex(const Complex_C_t& com)
{
    // e^(b+ic) = (e^b)(e^(ic)) = (e^b)((cos c) + i(sin c))
    // e^(b+ic) = e^b * cos(c) + i * e^b * sin(c)
    double eb = exp(com.m_real);

    return Complex_C_t{
        eb * cos(com.m_imagine),
        eb * sin(com.m_imagine)
    };
}

///--------------------------------------------------------
Complex_C_t powComplex(const Complex_C_t& base, const Complex_C_t& raise)
{
    /*
    see here for explanation: https://math.stackexchange.com/q/476998
    (a+ib) ^ (c+id) = e^( ln(r)*(c+id) + iθ*(c+id) )
    r = abs(a+ib)
    θ = arg(a+ib)

    real = ln(r)*c - d*θ
    imagine = i( d*ln(r) + cθ )
    eRaiseComplex(real + imagine)
    */

    double logAbs = log(base.absolute());
    double arg = base.argument();

    return raiseEComplex(Complex_C_t{
    logAbs * raise.m_real - raise.m_imagine * arg,
    logAbs * raise.m_imagine + raise.m_real * arg
    });
}