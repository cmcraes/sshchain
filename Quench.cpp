#include <iostream>
#include <mblas_gmp.h>
#include <mlapack_gmp.h>
#include <ctime>
#include <iomanip>
#include <fstream>
#include "mcomplex.h"

#define mint mpackint
#define N 12
#define N2 N*N
#define PREC 256
#define LDA N
using namespace std;

/** Loschmidt Echo Prgram
* 
* 
* Craig McRae
* 2016-07-03
*/

/*Template for matrix multiplication. Alot of Annoying integer arithmetic, because we're
 technically multiplying matrices but they're stored as 1xN*N Vectors, Not NxN matrices*/
template <typename Mat, typename matrix>
void matMult(matrix A, bool Adiag, Mat B, bool Bdiag, Mat C){

	if (!Adiag && !Bdiag)//Neither diagonal
	    for (int i = 0; i < N2; i++){
			for (int k = 0; k < N; k++){
				if (k == 0)	
					C[i] = A[(i/N)*N]*B[(i%N)];
				else
					C[i] = C[i] + A[(i/N)*N+k]*B[(i%N)+k*N];	
			}//for k
		}//for i
	else if (!Adiag && Bdiag)//Right Diagonal
		for (int i = 0; i < N2; i++)
			C[i] = A[i]*B[(i%N)*(N+1)];
	else if (Adiag && !Bdiag)//Left Diagonal
		for (int i = 0; i < N2; i++)
			C[i] = B[i]*A[(i/N)*(N+1)];
	else {
		for(int i = 0; i < N2; i++){//Both diagonal
			if (i%N == i/N)
				C[i] = A[i]*B[i];
			else
				C[i] = 0.0;
		}//for
	}//else
}//matMult

template <typename Mat, typename Matrix>
void matAdd(Matrix A, Mat B, Mat C){
	for (int i = 0; i < N2; i++)
		C[i] = A[i] + B[i];
}

template <typename Mat>
void matSub(Mat A, Mat B, Mat C){
	for (int i = 0; i < N2; i++)
		C[i] = A[i] - B[i];
}

template <typename PrintMatrix> 
void printMat(int n, int M, PrintMatrix A, int lda) {
	cout << '\n';
    	for (int j = 0; j < n; j++) {
		cout << "[ ";
		for (int i = 0; i < M; i++) {
		    cout << A[i + j * lda];
		    if (i < M - 1)
				cout << ", ";
		}
		cout << "]; \n";
	}
	cout << '\n';
}

void transpose(mpf_class* A){
	mpf_class temp;
	for (int i = 0; i < N2; i++){
		if (i/N > i%N){
			temp = A[i];
			A[i] = A[(i%N)*N + i/N];
			A[(i%N)*N + i/N] = temp;
		}
	}
}

void transpose(mcomplex* A){
	mcomplex temp;
	for (int i = 0; i < N2; i++){
		if (i/N > i%N){
			temp = A[i];
			A[i] = A[(i%N)*N + i/N];
			A[(i%N)*N + i/N] = temp;
		}
	}
}

void LU_decomp(mcomplex* a){
	for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
	        if (j >= i)
                for (int k = 0; k < i; k++)
                    a[j*N + i] = a[j*N + i] - a[j*N + k] * a[k*N + i];
        }
        for (int j = 0; j < N; j++) {
            if (j > i){
                a[i*N + j] = a[i*N + j] / a[i*N + i];
                for (int k = 0; k < i; k++)
                    a[i*N + j] = a[i*N + j] - ((a[i*N + k] * a[k*N + j]) / a[i*N + i]);
            }
        }
    }	
}//LU

mcomplex determinant(mcomplex* a){
	mcomplex det = mcomplex(1,0); 
	LU_decomp(a);	//Lu Decomp in-place of A
	
	//Determinant will be Product of diagonals
	for (int i = 0; i < N2; i+=(N+1))
		det = det * a[i];
	return det;
}

mpf_class determinant(mpf_class * A){
	mpf_class det = 1;
	mint *IPIV = new mint[N];
    	mint INFO;
	
    	Rgetrf(N, N, A, N, IPIV, &INFO);	//Finds the LU decomposition of A
	
	for (int i = 0; i < N2; i+=(N+1))
		det *= A[i];
	
    	delete[]IPIV;
	return det;
}

//get data for diagonal and subdiagonal entries
//i/N gives row, i%N gives column
void get_Diagonals(mpf_class* A, mpf_class* D, mpf_class* sD){
	for (int i = 0; i < N2; i++){
		if (i/N == i%N){
			if (D != NULL)
				D[i%N] = A[i];
			if (sD != NULL && i != 0)
				sD[(i-1)%N] = A[i-1];
		}
	}
}

void matrixExp(mpf_class* EigMat, mpf_class* EigVals, mpf_class time, mcomplex* workspace, mcomplex* ExpDiag){
	mcomplex Z;
	//Matrix of: e^(-i*t*lambda/h_bar) * Inverse EigenMatrix 	
	for (int i = 0; i < N; i++){
		Z = mcomplex(cos(time*EigVals[i]),-sin(time*EigVals[i]));
		for (int j = 0; j < N; j++)
			workspace[(i%N)*N + j] = Z*EigMat[(i%N)*N + j];
	}
	transpose(EigMat);
	//ExpDiag becomes: EigMat*e^(-i*t*lambda)*InvEigMat
	matMult(EigMat, false, workspace, false, ExpDiag); //ExpDiag now contains e^(-i*H_f*t) in original basis
}

int main(){
	mpf_set_default_prec(PREC);	//Numeric precision
	//local variables
	mcomplex* EigMat = new mcomplex[N2];	//EigenMat used to bring the ExpDiag into original magtrix
	mcomplex* ExpDiag = new mcomplex[N2];	//Complex matrix with e^diagonal eigen elements 
	mpf_class* R = new mpf_class[N2];		//Correlation Matrix
    	mpf_class* H_f = new mpf_class[N2];		//Tridiagonal Hamiltonians
	mpf_class* H_i = new mpf_class[N2];		//Entries are Conjugate of H_f
	mpf_class* D = new mpf_class[N];		//Array of diagonal elements of H
	mpf_class* sD = new mpf_class[N-1];		//Array of Subdiagonal Elements of H
    	mpf_class* work = new mpf_class[2*N-2]; //WorkSpace
	mpf_class* val = new mpf_class[2];		//Array Of values to fill H&E with
	mcomplex det;							//determinant
	mpf_class delta;						//Dimerization parameter
	mpf_class timeI, timeF, dt, curr;		//Time evolution parameters			
	mpf_class L = N;						//Chain Size
	mint* info = new mint[2];				//Variable to tell if an Mpack function converged or not
	ofstream output;
	bool index = true;						//Index of val array, used for filling A

	//get parameters
	cout << "What is the value of the dimerization parameter? ";
	cin >> delta;
	cout << "Enter the initial time: ";
	cin >> timeI;
	cout << "Enter the final time: ";
	cin >> timeF;
	cout << "Enter the change in time to calculate: ";
	cin >> dt;
	
	//Get System time for rough approixmation of Algorithm run time
	time_t result = std::time(NULL);
    	cout << asctime(localtime(&result));

	//set matrix inputs
	val[0] = 1 + delta;
	val[1] = 1 - delta;

	//load matrix H and E with 0's exept on super and sub diagonal 
	for (int j = 0; j < N2; j++){
		if (j/N == j%N){
			H_f[j] = 0.0;
			H_i[j] = 0.0;	
			if (j != 0){
				H_f[j-1] = val[index];
				H_i[j-1] = val[!index];
			}
			index = !index;
			if (j != N2-1){
				H_f[j+1] = val[index];
				H_i[j+1] = val[!index];
			}
			++j;
		}
		else{
			H_f[j] = 0.0;
			H_i[j] = 0.0;
		}
	}//for

	/*Rsteqr finds the Eigen values and Vectors of H. On exit H contains the matrix of 
	Eiegnvectors row-wise ascending W.R.T. eigenvalue, and D the array of eigenvalues*/
	get_Diagonals(H_i,D,sD);
    	Rsteqr("I", N, D, sD, H_i, N, work, &info[0]);
	
	get_Diagonals(H_f,D,sD);
    	Rsteqr("I", N, D, sD, H_f, N, work, &info[1]);

	if (info[0] == 0 && info[1] == 0) {//mlapack converged!
		//Create Correlation matrix R
		for (int i = 0; i < N2; i++)
			for (int j = 0; j < N/2; j++)
				R[i] = R[i] + H_i[j*N+i%N]*H_i[j*N+i/N];
		
		//H_i = Identity - R (Reuse H_i for space minimization & Preserving R)
		for (int i = 0; i < N2; i++){
			if (i%N == i/N)
				H_i[i] = 1.0 - R[i];
			else
				H_i[i] = -R[i];
		}
		output.open("output.txt");
		for (curr = timeI; curr < timeF; curr+=dt){
			matrixExp(H_f, D, curr, EigMat, ExpDiag); //ExpDiag now Contains e^(-i*H_final*t) in original basis
			matMult(R, false, ExpDiag, false, EigMat); //EigMat = R*e^(-i*H_final*t)
			matAdd(H_i, EigMat, EigMat);//EigMat = (I-R)+(R*e^(-i*H_final*t))
			det = determinant(EigMat); //Abs(Det(EigMat)) = Loschmidt Echo!!
			output << "(" << curr << "," << -log(Abs(det))/L << ")\n";
		}
		output.close();
	}
	else
		cout << "Rsteqr did not converge! :(\n";
	
	//Print out time again to find Apprx Algorithm run time
	time_t resut = std::time(NULL);
    	cout << asctime(localtime(&resut));

	//Free up workspace
    	delete[]H_f;
	delete[]H_i;
	delete[]D;
	delete[]sD;
    	delete[]work;
	delete[]val;
	delete[]EigMat;
	delete[]ExpDiag;
	delete[]R;	
}
