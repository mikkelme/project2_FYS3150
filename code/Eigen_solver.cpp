#include "Solver_functions.h"
#include <iostream>
#include <cmath>
#include <armadillo>
#include "time.h"

using namespace  std;
using namespace  arma;


int main(int argc, char *argv[]) {
  Jacobi my_functions; //Name of class
  int n, arg, iter, maxiter, p, q;
  double  Rmin, Rmax, h, hh, d, a, rho_max, w, tol, maxnondig;

  //Initialize problem
  n = atoi(argv[1]);
  Rmin = 0; Rmax = 1;

  h = (Rmax - Rmin)/n; hh = h*h; //Intregration step length
  d = 2.0 / (hh);   // DiagConst
  a =  -1.0 / (hh); //nondiagConst

  arg = atoi(argv[2]); //arg == 1: one electron, arg == 2: two electrons, arg == else: no potential
  if(arg == 1 || arg == 2){
    rho_max = atof(argv[3]);
    w = atof(argv[4]);
  }
  else{
    rho_max = 0.;
    w = 0.;
  }

  //Create Matrix
  mat A = my_functions.CreateTridiagonal(d, a, n, arg, rho_max, w);
  mat Eigvec = eye<mat>(n,n); //For eigenvectors on each colum
  vec Eigval(n); //Vector for eigenvalues

  //Run tests
  bool run_test = false;
  if (run_test){
    double test_tol = 1.0E-10;
    my_functions.OrthTest(test_tol);
    my_functions.EigValTest(test_tol);
    my_functions.MaxOffTest(test_tol);
  }

  //Select eigenvalue solver
  bool jacobi_solve = false;
  bool armadillo_solve = true;

  clock_t start, finish; //For timing

  //Jacobi Method
  if (jacobi_solve){
    tol = 1.0E-10;
    iter = 0;
    maxiter = 5e4;
    start = clock(); //Start timer
    my_functions.MaxOffdiag(A, p, q, maxnondig, n);
    while (maxnondig > tol && iter <= maxiter){
      my_functions.JacobiRotate(A, Eigvec, p, q, n);
      my_functions.MaxOffdiag(A, p, q, maxnondig, n);
      iter++;
    }
    finish = clock(); //End timer
    Eigval = my_functions.OrderEigenResults(A, Eigvec, n);
  }

  //Armdadillo eig_sym
  else if (armadillo_solve){
    start = clock(); //Start timer
    eig_sym(Eigval, Eigvec, A);
    finish = clock(); //End timer
    iter += 1;
  }
  else {
    cout << "Program failure: Must choose a solver method internally in Eigen_solver.cpp" << endl;;
    return 1;
  }
  double timeused = (double) (finish - start)/(CLOCKS_PER_SEC );

  //Calculate Analytical eigenvalues and eigenvectors
  pair<vec, mat> pau = my_functions.CalculateExact(d,a,n);
  vec Exact_eigval = pau.first;
  mat Exact_eigvec = pau.second;

  //Eigenvector test
  if (run_test){
    double test_tol = 1.0E-10;
    mat A_original = my_functions.CreateTridiagonal(d, a, n, arg, rho_max, w);
    my_functions.EigvecTest(test_tol, Eigval, Eigvec, A_original);
  }


  //Write results
  my_functions.WriteIter(iter);
  my_functions.WriteMeanError(Eigval, Exact_eigval, n);
  my_functions.WriteTime(timeused);
  my_functions.WriteEig(Eigval, Eigvec, n);

  //Print results (Only for n ≤ 10)
  //my_functions.PrintResults(Eigval, Eigvec, Exact_eigval, Exact_eigvec, n, jacobi_solve, armadillo_solve, arg);



  return 0;
}
