#include "JacobiMethod.h"
#include <iostream>
#include <cmath>
#include <armadillo>

using namespace  std;
using namespace  arma;


mat Jacobi::CreateMatrix(double Rmin, double Rmax, int n){
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


void Jacobi::MaxOffdiag(mat A, int& p, int& q, double& maxnondig, int n){
  double max;
  for (int i = 0; i < n; ++i){
    for (int j= i + 1; j < n; ++j){
      double aij = fabs(A(i,j));
      if (aij > max)
      {
        max = aij; p = i; q = j;
      }
    }
  }
  maxnondig = fabs(A(p,q));
}

void Jacobi::JacobiRotate(mat& A, mat& R, int k, int l, int n){
  double tau, tan, cos, sin;

  if (A(k,l) != 0.0){ //avoid divison by zero
    tau = (A(l,l) - A(k,k))/(2*A(k,l));

    //Choose right sign for tau calculaiton
    if (tau >= 0){
      tan = 1.0/(tau + sqrt(tau*tau + 1.0)); }
    else{
      tan = -1.0/(-tau + sqrt(tau*tau + 1.0)); }

    cos = 1.0/sqrt(1.0 + tan*tan);
    sin = cos*tan;
  }
  else {
    tan = 0.0;
    cos = 1.0;
  }

  //Perform similarity transformation A_new = (S^T)(A)(S)
  double a_kk_old = A(k,k);
  double a_ll_old = A(l,l);
  A(k,k) = a_kk_old*cos*cos - 2.0*A(k,l)*cos*sin + a_ll_old*sin*sin;
  A(l,l) = a_ll_old*cos*cos + 2.0*A(k,l)*cos*sin + a_kk_old*sin*sin;
  A(k,l) = 0.0; //Hard coding to zero
  A(l,k) = 0.0; //--------||--------
  for (int i = 0; i < n; ++i){
    if (i != k && i != l){
      double a_ik_old = A(i,k);
      double a_il_old = A(i,l);
      A(i,k) = a_ik_old*cos - a_il_old*sin;
      A(k,i) = A(i,k);
      A(i,l) = a_il_old*cos + a_ik_old*sin;
      A(l,i) = A(i,l);
    }
    //eigenvectors
    double r_ik_old = R(i,k);
    double r_il_old = R(i,l);

    R(i,k) = r_ik_old*cos - r_il_old*sin;
    R(i,l) = r_il_old*cos + r_ik_old*sin;
  }

return;

}
