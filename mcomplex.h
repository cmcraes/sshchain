#include <iostream>

class mcomplex; 

mcomplex Conj(const mcomplex& z);
mpf_class Abs(const mcomplex& z);

class mcomplex{
public:
	mpf_class real;
	mpf_class imaginary;

	mcomplex(mpf_class x = 0, mpf_class y = 0) : real(x), imaginary(y) {}
	
	mcomplex(mpf_class val) : real(val), imaginary(0) {};	
	mcomplex(double val) : real(val), imaginary(0) {};	
	mcomplex(float val) : real(val), imaginary(0) {};	
	mcomplex(int val) : real(val), imaginary(0) {};	

    //assignment operator modifies object, therefore non-const
    mcomplex& operator=(const mcomplex& a){
        real = a.real;
        imaginary = a.imaginary;
        return *this;
    }
    mcomplex& operator=(const mpf_class& a){
        real = a;
        imaginary = 0.0;
        return *this;
    }
    mcomplex& operator=(const double& a){
        real = a;
        imaginary = 0.0;
        return *this;
    }
    mcomplex& operator=(const float& a){
        real = a;
        imaginary = 0.0;
        return *this;
    }
    mcomplex& operator=(const int& a){
        real = a;
        imaginary = 0.0;
        return *this;
    }
    //negation operator 
    mcomplex operator-(){
        return mcomplex(-real, -imaginary);
    }
    //add op. doesn't modify object therefore const
    mcomplex operator+(const mcomplex& a) const{
        return mcomplex(real + a.real, imaginary + a.imaginary);
    }
    mcomplex operator+(const mpf_class& a) const{
        return mcomplex(real + a, imaginary);
    }
    mcomplex operator+(const double& a) const{
        return mcomplex(real + a, imaginary);
    }
    mcomplex operator+(const float& a) const{
        return mcomplex(real + a, imaginary);
    }
    mcomplex operator+(const int& a) const{
        return mcomplex(real + a, imaginary);
    }
	//Subtraction
    mcomplex operator-(const mcomplex& a) const{
        return mcomplex(real - a.real, imaginary - a.imaginary);
    }
    mcomplex operator-(const mpf_class& a) const{
        return mcomplex(real - a, imaginary);
    }
    mcomplex operator-(const double& a) const{
        return mcomplex(real - a, imaginary);
    }
    mcomplex operator-(const float& a) const{
        return mcomplex(real - a, imaginary);
    }
    mcomplex operator-(const int& a) const{
        return mcomplex(real - a, imaginary);
    }
	//mult op. doesn't modify object therefore const
    mcomplex operator*(const mcomplex& a) const{
        return mcomplex((a.real*real - a.imaginary*imaginary), (a.real*imaginary + a.imaginary*real));
    }
    mcomplex operator*(const mpf_class& a) const{
        return mcomplex(real*a, imaginary*a);
    }
    mcomplex operator*(const double& a) const{
        return mcomplex(real*a, imaginary*a);
    }
    mcomplex operator*(const float& a) const{
        return mcomplex(real*a, imaginary*a);
    }
    mcomplex operator*(const int& a) const{
        return mcomplex(real*a, imaginary*a);
    }
    mcomplex operator/(const mpf_class& a) const{
        return mcomplex(real/a, imaginary/a);
    }
	mcomplex operator/(const mcomplex& a) const{
		mpf_class size = (a.real*a.real + a.imaginary*a.imaginary);
		return mcomplex ((real*a.real + imaginary*a.imaginary), (imaginary*a.real - real*a.imaginary))/size;
    }
	friend std::ostream& operator<<(std::ostream& os, const mcomplex& z);

};

//add ops
mcomplex operator+(const mpf_class& a, const mcomplex& z) {
    return mcomplex(z.real + a, z.imaginary);
}
mcomplex operator+(const double& a, const mcomplex& z) {
    return mcomplex(z.real + a, z.imaginary);
}
mcomplex operator+(const float& a, const mcomplex& z) {
    return mcomplex(z.real + a, z.imaginary);
}
mcomplex operator+(const int& a, const mcomplex& z) {
    return mcomplex(z.real + a, z.imaginary);
}
//mult ops
mcomplex operator*(const mpf_class& a, const mcomplex& z) {
    return mcomplex(z.real*a, z.imaginary*a);
}
mcomplex operator*(const double& a, const mcomplex& z) {
    return mcomplex(z.real*a, z.imaginary*a);
}
mcomplex operator*(const float& a, const mcomplex& z) {
    return mcomplex(z.real*a, z.imaginary*a);
}
mcomplex operator*(const int& a, const mcomplex& z) {
    return mcomplex(z.real*a, z.imaginary*a);
}

//define output of mcomplex 
std::ostream& operator<<(std::ostream& os, const mcomplex& z) {
    os << z.real << '+' << z.imaginary << 'i';
    return os;
}

mcomplex Conj(const mcomplex& z){
	return mcomplex(z.real, -z.imaginary);
}
mpf_class Abs(const mcomplex& z){
	return sqrt((z*Conj(z)).real);	
}
