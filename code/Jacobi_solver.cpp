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
  int n = 10;
  mat A = my_functions.CreateMatrix(Rmin, Rmax, n);
  //my_functions.ShowMatrix(A);

  double tol = 1.0E-10;
  int iter = 0;
  int maxiter = 10;
  int maxnondig = 1;
  int p; int q;
  while (maxnondig > tol && iter <= maxiter){
    my_functions.MaxOffdiag(A, &p, &q, n);
    cout << p <<" "<< q << " "<< A(p,q) << endl;
  //  JacobiRotate(A, R, p, q, n);
  //  maxnondig = A(&p, &q)
    iter++;
  }


}











//
