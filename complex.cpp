//
// Created by philipp on 24/06/16.
//

#include "complex.h"

/**
     * \param x real part
     * \param y imaginary part
     * Constructor with mpf_class numbers
     */
Complex::Complex(mp::gmp_float<PRECISION> x, mp::gmp_float<PRECISION> y) {
    *real = x;
    *imaginary = y;
}


/**
 * \param x real part
 * Constructor to convert real numbers. Imaginary part is set to 0.
 */
Complex::Complex(mp::gmp_float<PRECISION> x) {
    *real = x;
    *imaginary = 0;
}

/**
 * Assign operator
 * \param a complex Number to be assigned
 */
Complex Complex::operator=(const Complex& a) {
    *real = a.real;
    *imaginary = a.imaginary;
}

/**
 * Silently cast real number to Complex type
 * \param a real Number to be assigned
 */
Complex Complex::operator=(const mp::gmp_float<PRECISION>& a){
    *real = a;
    *imaginary = 0;
}

/**
 * Addition operator
 * \param a complex Number to be assigned
 */
Complex Complex::operator+(const Complex *a) {
    real = *real + (mp::gmp_float &) a->real;
    //mpf_add (imaginary, imaginary, a.imaginary);
}
Complex operator+(const mp::gmp_float<PRECISION>& a) {}
Complex operator*(const Complex& a) {}
Complex operator*(const mp::gmp_float<PRECISION>& a) {}
