#include <iostream>
#include <mblas_gmp.h>
#include <mlapack_gmp.h>
#include <ctime>
#include <iomanip>

#include "mpcomplex.h"

#define mint mpackint
#define N 100
#define N2 N*N
#define PREC 512
#define LDA N
using namespace std;

/*
 * TODO:: change to template or class
 * then you can use is also externally
 *
 * I would move this into an external class which also provides methods like abs, cconjugate, etc.
 */
struct mcomplex {
    mpf_class real;
    mpf_class imaginary;

    /**
     * \param x real part
     * \param y imaginary part
     * Constructor with mpf_class numbers
     */
    mcomplex(mpf_class x = 0, mpf_class y = 0) : real(x), imaginary(y) { }


    /**
     * \param x real part
     * Constructor to convert real numbers. Imaginary part id set to 0.
     * TODO change type to mpf_class
     */
    mcomplex(double val) : real(val), imaginary(0) { };


    // assignment operator modifies object, therefore non-const
    mcomplex &operator=(const mcomplex &a) {
        real = a.real;
        imaginary = a.imaginary;
        return *this;
    }

    // add op. doesn't modify object therefore const
    mcomplex operator+(const mcomplex &a) const {
        return mcomplex(real + a.real, imaginary + a.imaginary);
    }

    mcomplex operator+(const mpf_class &a) const {
        return mcomplex(real + a, imaginary);
    }

    // add op. doesn't modify object therefore const
    mcomplex operator*(const mcomplex &a) const {
        return mcomplex((a.real * real - a.imaginary * imaginary), (a.real * imaginary + a.imaginary * real));
    }

    mcomplex operator*(const mpf_class &a) const {
        return mcomplex(real * a, imaginary * a);
    }

};

/* Template for matrix multiplication. Alot of Annoying integer arithmetic, because were techinaclly multiplying
 * matrices but theyre stored as 1xN*N Vectors, Not NxN matrices*/
template<typename Mat>
Mat matMult(Mat A, Mat B, Mat C, bool Adiag, bool Bdiag) {

    if (!Adiag && !Bdiag) {//Neither diagonal
        for (int i = 0; i < N2; i++)
            for (int k = 0; k < N; k++)
                C[i] = C[i] + A[(i / N) * N + k] * B[(i % N) + k * N];
    }//if
    else if (!Adiag && Bdiag)//Right Diagonal
        for (int i = 0; i < N2; i++)
            C[i] = A[i] * B[(i % N) * (N + 1)];

    else if (Adiag && !Bdiag)//Left Diagonal
        for (int i = 0; i < N2; i++)
            C[i] = B[i] * A[(i / N) * (N + 1)];
    else {
        for (int i = 0; i < N2; i++) {//Both diagonal
            if (i % N == i / N)
                C[i] = A[i] * B[i];
            else
                C[i] = 0.0;
        }
    }
}

template<typename PrintMatrix>
PrintMatrix printmat(int n, int M, PrintMatrix A, int lda) {
    mpf_class mtmp;
    printf("\n");
    for (int i = 0; i < n; i++) {
        printf("[ ");
        for (int j = 0; j < M; j++) {
            mtmp = A[i + j * lda];
            gmp_printf("%8.6Fe", mtmp.get_mpf_t());
            if (j < M - 1)
                printf(", ");
        }
        printf("]; \n");
    }
    printf("\n");
}

void inverse(mpf_class *A, mpf_class *B) {
    /**
     * \param A the matrix to invert
     * \param B pointer to inverse
     *
     * TODO:: add check if A is invertible and return false if not
     */
    mint *IPIV = new mint[N];
    mint LWORK = -1;
    mpf_class *WORK = new mpf_class[1];
    mint INFO;

    //B becomes copy of A
    for (int j = 0; j < N2; j++)
        B[j] = A[j];

    //Workspace query
    Rgetri(N, A, N, IPIV, WORK, LWORK, &INFO);
    LWORK = int(WORK[0].get_d());
    delete[]WORK;
    WORK = new mpf_class[max(1, (int) LWORK)];

    Rgetrf(N, N, B, N, IPIV, &INFO);            //Finds the LU decomposition of B and Stores it in B
    Rgetri(N, B, N, IPIV, WORK, LWORK, &INFO);    //Finds Inverse of B from the LU decomp

    delete[]IPIV;
    delete[]WORK;
}

int main() {

    cout << fixed << setprecision(20);    //Output precision
    mpf_set_default_prec(PREC);            //Numeric precision
    const mpf_class H_BAR_INV = (9.4825217771312202812908507856980167254410816649238827 * pow(10, 33));

    //local variables
    mcomplex *E = new mcomplex[N2];            //Complex matrix with e^diagonal eigen elements
    mpf_class *A = new mpf_class[N2];        //Tridiagonal Hamiltonian
    mpf_class *B = new mpf_class[N2];        //Inverse EigenVector Matrix
    mpf_class *C = new mpf_class[N2];        //Final Result To be Held In C
    mpf_class *D = new mpf_class[N];        //Array of diagonal elements of A
    mpf_class *sD = new mpf_class[N - 1];        //Array of Subdiagonal Elements
    mpf_class *work = new mpf_class[2 * N - 2]; //WorkSpace
    mpf_class *val = new mpf_class[2];        //Array Of values to fill A with
    mpf_class delta;                        //Dimerization parameter
    mpf_class time;                            //Time evolution parameter
    mint info;                                //Variable to tell if an Mpack function converged or not
    bool index = false;                        //Index of val array, used for filling A

    //get parameters
    cout << "What is the value of the dimerization parameter? ";
    cin >> delta;

    cout << "Enter the time Evolution parameter: ";
    cin >> time;

    //Get System time for rough approxmation of Algorithm run time
    time_t result = std::time(NULL);
    cout << asctime(localtime(&result));

    //set matrix inputs
    val[0] = 1 + delta;
    val[1] = 1 - delta;

    //load matrix A with 0's exept on super and sub diagonal
    for (int i = 0; i < N2; i++) {
        if (i / N == i % N) {
            A[i] = 0;
            if (i != 0)
                A[i - 1] = val[index];
            index = !index;
            if (i != N2 - 1)
                A[i + 1] = val[index];
            ++i;
        }
        else
            A[i] = 0.0;
    }

    //get data for diagonal and subdiagonal entries
    //i/N gives row, i%N gives column
    for (int i = 0; i < N2; i++) {
        if (i / N == i % N) {
            D[i % N] = A[i];
            if (i != 0)
                sD[(i - 1) % N] = A[i - 1];
        }
    }

    //Rsteqr finds the Eigen values and Vectors of A,
    //on exit, A constains Matrix of Eiegnvectors, and D the array of eigenvalues
    Rsteqr("I", N, D, sD, A, N, work, &info);
    if (info == 0) {
        //Print out conrging Eigenvalues to see if they have the same magnitude at PREC
        //cout << D[N/2-1] << '\n' << D[N/2] << '\n';

        inverse(A, B);    //B becomes inverse of A

        //Fill E with e^(t*lambda/h_bar) & Zeros
        for (int i = 0; i < N2; i++) {
            if (i / N == i % N)
                E[i] = mcomplex(cos(time * D[i % N] * H_BAR_INV), -sin(time * D[i % N] * H_BAR_INV));
            else
                E[i] = 0;
        }

        cout << "Lowest Eigenvalue: " << A[0] << endl;
        //matMult(A, B, C, false, false);
        //printmat(N, N, E, N);
    }
    else
        cout << "Rsteqr did not converge! :(\n";
    //Free up workspace
    delete[]work;
    delete[]A;
    delete[]B;
    delete[]D;
    delete[]sD;
    delete[]val;

    //Print out time again to find Apprx Algorithm run time
    time_t resut = std::time(NULL);
    cout << asctime(localtime(&resut));
}