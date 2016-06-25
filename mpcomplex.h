//
// Created by philipp on 24/06/16.
//

#ifndef SSH_CRAIG_COMPLEX_H
#define SSH_CRAIG_COMPLEX_H

#include <boost/multiprecision/gmp.hpp>

#ifndef PRECISION
#pragma message "PRECISION not defined - using 0"
#define PRECISION 0
#endif


using namespace std;
namespace mp = boost::multiprecision;

typedef mp::number<mp::gmp_float<PRECISION> > mpfloat_t;

/**
 * Implements a complex number type based on the boost GMP backend
 *
 * GMP Precision ist set to 1000 (meaning arbitrary precision). This can be overridden by
 * ```{c}
 * #define PRECISION {yourprecision}
 * ```
 * (before you `#include "complex.h"`) or by passing PRECISION to the preprocessor:
 * ```{bash}
 * gcc -D PRECISION=yourprecision
 * ```
 * in case you use gcc.
 */
//template <int PRECISION>
class mpcomplex {
public:
    //members
    mpfloat_t re;
    mpfloat_t im;

    //constructors
    mpcomplex(const mpfloat_t x, const mpfloat_t y);
    mpcomplex(const mpfloat_t x);
    mpcomplex();

    //deconstructor
    ~mpcomplex();

    //assign operators
    mpcomplex operator=(const mpcomplex a);
    mpcomplex operator=(const mpfloat_t a);

    //boolean operators

    // be compatible with std::complex
    mpfloat_t real() {
        return re;
    }
    mpfloat_t imag() {
        return im;
    }

    //complex conjugate
    mpcomplex cconj();
};

//aritmethic operators
mpcomplex operator+(const mpcomplex a, const mpcomplex b);
mpcomplex operator+(const mpcomplex a, const mpfloat_t b);

mpcomplex operator-(const mpcomplex a, const mpcomplex b);
mpcomplex operator-(const mpcomplex a, const mpfloat_t b);

mpcomplex operator*(const mpcomplex a, const mpcomplex b);
mpcomplex operator*(const mpcomplex a, const mpfloat_t b);

mpcomplex operator/(const mpcomplex a, const mpfloat_t b);
mpcomplex operator/(const mpcomplex a, const mpcomplex b);

//analytic functions
/*
 * TODO
 */
mpfloat_t abs(mpcomplex arg);
mpcomplex sin(mpcomplex arg);
mpcomplex cos(mpcomplex arg);
mpcomplex asin(mpcomplex arg);
mpcomplex acos(mpcomplex arg);
mpcomplex tan(mpcomplex arg);
mpcomplex atan(mpcomplex arg);
mpcomplex sinh(mpcomplex arg);
mpcomplex cosh(mpcomplex arg);
mpcomplex asinh(mpcomplex arg);
mpcomplex acosh(mpcomplex arg);
mpcomplex tanh(mpcomplex arg);
mpcomplex atanh(mpcomplex arg);
mpcomplex exp(mpcomplex arg);
mpcomplex log(mpcomplex arg);



#endif //SSH_CRAIG_COMPLEX_H
