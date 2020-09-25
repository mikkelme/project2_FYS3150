#include "JacobiMethod.h"
#include <iostream>
#include <cmath>
#include <armadillo>

using namespace  std;
using namespace  arma;


mat Jacobi::CreateMatrix(double Rmin, double Rmax,  int n){
  int     i, j;
  double  h, hh, d, a;

  //Intregration step length
  h = (Rmax - Rmin)/n; hh = h*h;
  d = 2.0 / (hh);   // DiagConst
  a =  -1.0 / (hh); //nondiagConst
  // Setting up tridiagonal matrix A
  mat A = zeros<mat>(n,n);
  A(0,0) = d;
  A(0,1) = a;
  for(i = 1; i < n-1; i++) {
    A(i,i-1) = a;
    A(i,i)   = d;
    A(i,i+1) = a;
  }
  A(n-1,n-2) = a;
  A(n-1,n-1) = d;
  return A;
}


void Jacobi::ShowMatrix(mat A){
  int ni = size(A)[0];
  int nj = size(A)[1];
  for (int i = 0; i < ni; i++){
    for (int j = 0; j < nj; j++){
        cout << setw(4) << A(i,j) << " ";
    }
    cout << endl;
  }
}


void Jacobi::MaxOffdiag(mat A, int *p, int*q, int n){
  double max;
  for (int i = 0; i < n; ++i){
    for (int j= i + 1; j < n; ++j){
      double aij = fabs(A(i,j));
      if (aij > max)
      {
        max = aij; *p = i; *q = j;
      }
    }
  }
}

void JacobiRotate(mat A, mat R, int k, int l, int n){
  double tau, tan, cos, sin;

  if (A(k,l) != 0.0){ //avoid divison by zero
    tau = (A(l,l) - A(k,k))/2*A(k,l);

    //Choose right sign for tau calculaiton
    if (tau >= 0){
      

    }

    tan = 1/(tau + sqrt(tau*tau + 1))
    cos = 1/sqrt(1 + tan*tan)
    sin = cos*tan

  }
  else {
    tan = 0.0
    cos = 1.0

  }





}

/*
int Jacobi::test(int a){
  b = 2;
  c = 3;
  d = 4;
  return a + b + c + d;
}
*/
