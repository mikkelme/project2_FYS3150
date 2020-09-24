//  Diagonalizing tridiagonal Toeplitz matrix  with Lapack functions
//  Compile as c++ -O3 -o Tridiag.exe TridiagToeplitz.cpp -larmadillo -llapack -lblas

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>

using namespace  std;
using namespace  arma;


int main(int argc, char* argv[]){

  int ;
  double eps, max_iterations, n;
  n = 
  max_iterations = n*n*n
  eps = 1e-8;
  // int iterations = 0;
  while(max_element(A)>eps && double iterations < max_iterations){
  	max = max_element(A);
  	rotate(A, R, k, l, n);
  	compute(A, R, k, l, n);
  	iterations++

  }	

  a_kl = max_element;
  tau = (A[l][l]-A[k][k])/(2*A[k][l])
  t1 = -tau + sqrt(1+tau*tau)
  t2 = -tau - sqrt(1+tau*tau)
  c = 

  return 0;
}
void rotate(double ** A, double ** R, int k, int l, int n){

}

void compute(double var1, var2, ){
	double tau, t, s, c;
	tau =
	if(A[k][l] != 0.0){
		t = 
	


	a_kk = A[k][k];
	a_ll = A[l][l];
	A[k][k] = c*c*a_kk -2.*c*s*A[k][l] + s*s*a_ll;
	A[l][l] = s*s*a_kk +2.*c*s*A[k][l] + c*c*a_ll;
	A[k][l] = 0.0;
	A[l][k] = 0.0;

	for(int i=0; i<nM i++){
		if(i != k && i != l){
			a_ik = A[i][k];
			a_il = A[i][l];
			A[i][k] = A[k][i] = c*a_ik - s*a_il;
			A[i][l] = A[l][i] = c*a_il + s*a_ik;
		}
	r_ik = R[i][k];
	r_il = R[i][l];
	R[i][k] = c*r_ik - s*r_il;
	R[i][l] = c*r_il + s*r_ik;
	}


}