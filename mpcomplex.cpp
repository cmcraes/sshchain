//
// Created by philipp on 24/06/16.
//

#include "mpcomplex.h"

/**
* Constructor with GMP numbers
* \param x real part
* \param y imaginary part
*/
mpcomplex::mpcomplex(const mpfloat_t x, const mpfloat_t y) {
    re = x;
    im = y;
}

/**
* Constructor with one GMP number
* \param x real part
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
    mpf_clear((mpf_ptr &) re);
    mpf_clear((mpf_ptr &) im);
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

/**
 * Conjugate the number
 * \returns the complex conjugate of the number it is called on
 */
mpcomplex mpcomplex::cconj() {
    return mpcomplex(re, -im);
}

/**
 * Conjugate the number
 * \param a complex number
 * \returns the complex conjugate of the input
 */
mpcomplex mpcomplex::cconj(mpcomplex c) {
    return mpcomplex(c.re, -c.im);
}

/**
 * Addition operator
 * for two mpcomplex objects
 * \memberof mpcomplex
 */
mpcomplex operator+(const mpcomplex a, const mpcomplex b) {
    return mpcomplex(b.re + a.re, b.im + a.im);
}

/**
 * Addition operator
 * for compatibiliy with other number types
 * \memberof mpcomplex
 */
mpcomplex operator+(const mpcomplex a, const mpfloat_t b) {
    return a + mpcomplex(b);
}

/**
 * Substraction operator
<<<<<<< HEAD
 * for two mpcomplex objects
 * Accepts mpcomplex, and mpfloat_t (which wraps the standard types)
 * \memberof mpcomplex
 */
mpcomplex operator-(const mpcomplex a, const mpcomplex b) {
    return mpcomplex(b.re - a.re, b.im - a.im);
}

/**
 * Substraction operator
 * for compatibiliy with other number types
 * \memberof mpcomplex
 */
mpcomplex operator-(const mpcomplex a, const mpfloat_t b) {
    return a - mpcomplex(b);
}

/**
 * Multiplication operator
 * for two mpcomplex objects
 * \memberof mpcomplex
 */
mpcomplex operator*(const mpcomplex a, const mpcomplex b) {
    return mpcomplex((a.re * b.re - a.im * b.im), (a.re * b.im + a.im * b.re));
}

/**
 * Multiplication operator
 * for compatibiliy with other number types
 * \memberof mpcomplex
 */
mpcomplex operator*(const mpcomplex a, const mpfloat_t b) {
    return a - mpcomplex(b);
}

/**
 * Division operator
 * for two mpcomplex objects
 * \memberof mpcomplex
 */
mpcomplex operator/(const mpcomplex a, const mpfloat_t b) {
    return mpcomplex(a.re / b, a.im / b);
}

/**
 * Division operator
 * for compatibiliy with other number types
 * \memberof mpcomplex
 */
mpcomplex operator/(const mpcomplex a, const mpcomplex b) {
    mpcomplex tmp = b;
    mpcomplex denominator_cmplx = b * tmp.cconj();
    mpfloat_t denominator = denominator_cmplx.re;
    delete &tmp;
    delete &denominator_cmplx;
    mpcomplex counter = a * b;
    return counter / denominator;
}
