//
// Created by philipp on 24/06/16.
//

#include "mpcomplex.h"

/**
* \param x real part
* \param y imaginary part
* Constructor with GMP numbers
*/
mpcomplex::mpcomplex(const mpfloat_t x, const mpfloat_t y) {
    re = x;
    im = y;
}
/**
* \param x real part
* Constructor with one GMP number
*/
mpcomplex::mpcomplex(const mpfloat_t x) {
    re = x;
    im = 0;
}

/**
 * Constructor with no argument. Will assign NULL.
 */
mpcomplex::mpcomplex() {
    re = NULL;
    im = NULL;
}

/**
 * Doconstructor to free unneeded memory
 */
mpcomplex::~mpcomplex() {
    mpf_clear((mpf_ptr &)re);
    mpf_clear((mpf_ptr &)im);
}
/**
 * Assign operator
 * \param a complex Number to be assigned
 */
mpcomplex mpcomplex::operator=(const mpcomplex a) {
    re = a.re;
    im = a.im;
    return *this;
}

/**
 * Assign operator
 * \param a real umber to be assigned
 */
mpcomplex mpcomplex::operator=(const mpfloat_t a) {
    re = (mpfloat_t) a;
    im = 0;
    return *this;
}

mpcomplex mpcomplex::cconj() {
    return mpcomplex(re,-im);
}

/**
 * Addition operator
 * Accepts mpcomplex, and mpfloat_t (which wraps the standard types)
 */
mpcomplex operator+(const mpcomplex a, const mpcomplex b) {
    return mpcomplex(b.re + a.re, b.im + a.im);
}
mpcomplex operator+(const mpcomplex a, const mpfloat_t b) {
    return a+mpcomplex(b);
}

/**
 * Substraction operator
 * Accepts mpcomplex, and mpfloat_t (which wraps the standard types)
 */
mpcomplex operator-(const mpcomplex a, const mpcomplex b) {
    return mpcomplex(b.re - a.re, b.im - a.im);
}
mpcomplex operator-(const mpcomplex a, const mpfloat_t b) {
    return a-mpcomplex(b);
}

/**
 * Multiplication operator
 * Accepts mpcomplex, and mpfloat_t (which wraps the standard types)
 */
mpcomplex operator*(const mpcomplex a, const mpcomplex b) {
    return mpcomplex((a.re*b.re - a.im*b.im), (a.re*b.im + a.im*b.re));
}
mpcomplex operator*(const mpcomplex a, const mpfloat_t b) {
    return a-mpcomplex(b);
}

/**
 * Division operator
 * Accepts mpcomplex, and mpfloat_t (which wraps the standard types)
 */
mpcomplex operator/(const mpcomplex a, const mpfloat_t b) {
    return mpcomplex(a.re/b,a.im/b);
}
mpcomplex operator/(const mpcomplex a, const mpcomplex b) {
    mpcomplex tmp=b;
    mpcomplex denominator_cmplx = b * tmp.cconj();
    mpfloat_t denominator = denominator_cmplx.re;
    delete &tmp;
    delete &denominator_cmplx;
    mpcomplex counter = a * b;
    return counter/denominator;
}
