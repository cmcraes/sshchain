//
// Created by philipp on 24/06/16.
//

#ifndef SSH_CRAIG_COMPLEX_H
#define SSH_CRAIG_COMPLEX_H

#include <boost/multiprecision/gmp.hpp>
#include <mblas_gmp.h>
#include <mlapack_gmp.h>

#ifndef PRECISION
#define PRECISION 1000
#endif


using namespace std;
namespace mp = boost::multiprecision;

/**
 * Implements a complex number type based on the boost GMP backend
 */
//template <int PRECISION>
class Complex {
public:
    mp::gmp_float<PRECISION> *real;
    mp::gmp_float<PRECISION> *imaginary;

    Complex(mp::gmp_float<PRECISION> real, mp::gmp_float<PRECISION> imaginary);
    Complex(mp::gmp_float<PRECISION> real);

    Complex operator=(const Complex *a);
    Complex operator=(const mp::gmp_float<PRECISION>& a);
    Complex operator+(const Complex& a);
    Complex operator+(const mp::gmp_float<PRECISION>& a);
    Complex operator*(const Complex& a);
    Complex operator*(const mp::gmp_float<PRECISION>& a);



};


#endif //SSH_CRAIG_COMPLEX_H
