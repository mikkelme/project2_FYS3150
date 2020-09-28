#ifndef SOLVER_FUNCTIONS_H
#define SOLVER_FUNCTIONS_H

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
  vec OrderEigenResults(mat& A, mat&R, int n);
  vec CalculateExact(double d, double a, int n);
  void WriteIter(int iter);
  void WriteMeanError(vec& Eigval, vec& Exact, int n);
  void PrintResults(vec& Eigval, vec& Exact, mat& R, int n, bool jacobi_solve, bool armadillo_solve, int arg);
  void WriteTime(double timeused);
  void OrthTest(double tol);
  void EigValTest(double tol);
  void MaxOffTest(double tol);


};


#endif
