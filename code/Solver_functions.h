#ifndef SOLVER_FUNCTIONS_H
#define SOLVER_FUNCTIONS_H

#include <armadillo>
using namespace std;
using namespace  arma;

class Jacobi {
private:



public:
  //void Task(int, int);
  mat CreateTridiagonal(double d, double a, int n, int arg, double rho_max, double w);
  void ShowMatrix(mat A);
  void MaxOffdiag(mat A, int& p, int& q, double& maxnondig, int n);
  void JacobiRotate(mat& A, mat& Eigvec, int k, int l, int n);
  vec OrderEigenResults(mat& A, mat&Eigvec, int n);
  pair<vec, mat> CalculateExact(double d, double a, int n);
  void WriteIter(int iter);
  void WriteMeanError(vec& Eigval, vec& Exact, int n);
  void WriteTime(double timeused);
  void WriteEig(mat& Eigval, mat& Eigvec, int n);
  void PrintResults(vec& Eigval, mat& Eigvec, vec& Exact_eigval, mat& Exact_eigvec, int n, bool jacobi_solve, bool armadillo_solve, int arg);
  void OrthTest(double tol);
  void EigValTest(double tol);
  void EigvecTest(double tol, vec& Eigval, mat& Eigvec, mat& A_original);
  void MaxOffTest(double tol);

};


#endif
