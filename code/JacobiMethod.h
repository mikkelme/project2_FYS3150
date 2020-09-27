#ifndef JACOBIMETHOD_H
#define JACOBIMETHOD_H

#include <armadillo>
using namespace  arma;

class Jacobi {
private:



public:
  //void Task(int, int);
  mat CreateTridiagonal(double d, double a, int n, int arg, double rho_max, double w);
  void ShowMatrix(mat A);
  void MaxOffdiag(mat A, int& p, int& q, double& maxnondig, int n);
  void JacobiRotate(mat& A, mat& R, int k, int l, int n);
  mat OrderEigenResults(mat& A, mat&R, int n);
  void WriteIter(int iter);
  void WriteMeanError(mat& Eigval,double d, double a, int n);
  void WriteEig(mat& Eigval, mat& R, int n);
  void WriteTime(double timeused);
  void OrthTest(double tol);
  void EigValTest(double tol);
  void MaxOffTest(double tol);

  
};


#endif
