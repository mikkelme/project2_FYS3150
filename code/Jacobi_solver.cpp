#include "JacobiMethod.h"
#include <iostream>
#include <cmath>
#include <armadillo>

using namespace  std;
using namespace  arma;


int main(int argc, char const *argv[]) {
  Jacobi my_functions;
  //Create Matrix
  double  Rmin, Rmax, h, hh, d, a;
  int n;
  if (argc > 1){
    n = atoi(argv[1]); }
  else {
    n = 20; }
  Rmin = 0; Rmax = 1;
  //Intregration step length
  h = (Rmax - Rmin)/n; hh = h*h;
  d = 2.0 / (hh);   // DiagConst
  a =  -1.0 / (hh); //nondiagConst

  //Create Matrix
  mat A = my_functions.CreateTridiagonal(d, a, n);
  mat R = eye<mat>(n,n); //For eigenvectors on each row

  double tol = 1.0E-10;
  int iter = 1;
  int maxiter = 5e4;
  double maxnondig;
  int p; int q;
  my_functions.MaxOffdiag(A, p, q, maxnondig, n);
  while (maxnondig > tol && iter <= maxiter){
    my_functions.JacobiRotate(A, R, p, q, n);
    my_functions.MaxOffdiag(A, p, q, maxnondig, n);
    iter++;
  }
  //Order result and show
  mat Eigval = my_functions.OrderEigenResults(A, R, n);

  my_functions.WriteIter(iter - 1);
  /*
  double pi = acos(-1.0);
  for(int i = 0; i < n; i++) {
    double Exact = d+2*a*cos((i+1)*pi/(n+1));
    cout << Eigval(i) << " " << Exact << endl;
  }
  */





return 0;
}











//
