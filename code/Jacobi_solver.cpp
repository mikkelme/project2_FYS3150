#include "JacobiMethod.h"
#include <iostream>
#include <cmath>
#include <armadillo>

using namespace  std;
using namespace  arma;


/*
int main(int argc, char const *argv[]) {
  //Parameters
  int a = 3;
  int output;
  Jacobi my_solver;
  output = my_solver.test(a);
  cout << output << endl;
  return 0;
}
*/

int main(int argc, char const *argv[]) {
  Jacobi my_functions;
  //Create Matrix
  double Rmin = 0;
  double Rmax = 1;
  int n = 20;
  mat A = my_functions.CreateMatrix(Rmin, Rmax, n);
  mat R = eye<mat>(n,n); //For eigenvectors on each row
  //my_functions.ShowMatrix(A);

  double tol = 1.0E-10;
  int iter = 0;
  int maxiter = 1e3;
  double maxnondig;
  int p; int q;
  my_functions.MaxOffdiag(A, p, q, maxnondig, n);
  while (maxnondig > tol && iter <= maxiter){
    my_functions.JacobiRotate(A, R, p, q, n);
    my_functions.MaxOffdiag(A, p, q, maxnondig, n);
    iter++;
  }

  //Order result and show
  //my_functions.OrderEigenvalues(A)
  cout << "iter = " << iter << endl;
  for (int i = 0; i < n; i++){
    cout << A(i,i) << endl;
  }

return 0;
}











//
