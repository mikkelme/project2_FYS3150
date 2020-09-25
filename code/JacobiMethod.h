#ifndef JACOBIMETHOD_H
#define JACOBIMETHOD_H

#include <armadillo>
using namespace  arma;

class Jacobi {
private:



public:
  mat CreateMatrix(double Rmin, double Rmax, int n);
  void ShowMatrix(mat A);
  void MaxOffdiag(mat A, int *p, int*q, int n);
  void JacobiRotate(mat A, mat R, int k, int l, int n)


};


#endif
