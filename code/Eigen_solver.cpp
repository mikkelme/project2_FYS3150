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
  mat R = eye<mat>(n,n); //For eigenvectors on each row

  //Run tests
  bool run_test = false;
  if (run_test){
    my_functions.OrthTest(tol);
    my_functions.EigValTest(tol);
    my_functions.MaxOffTest(tol);
    cout << "All tests passed." << endl;
  }

  //Select eigenvalue solver
  bool jacobi_solve = true;
  bool armadillo_solve = false;

  clock_t start, finish; //For timing
  vec Eigval(n); //Vector for eigenvalues

  //Jacobi Method
  if (jacobi_solve){
    tol = 1.0E-10;
    iter = 0;
    maxiter = 5e4;
    start = clock(); //Start timer
    my_functions.MaxOffdiag(A, p, q, maxnondig, n);
    while (maxnondig > tol && iter <= maxiter){
      my_functions.JacobiRotate(A, R, p, q, n);
      my_functions.MaxOffdiag(A, p, q, maxnondig, n);
      iter++;
    }
    finish = clock(); //End timer
    Eigval = my_functions.OrderEigenResults(A, R, n);
  }

  //Armdadillo eig_sym
  else if (armadillo_solve){
    cout << "armadillo" << endl;
    start = clock(); //Start timer
    eig_sym(Eigval, A);
    finish = clock(); //End timer
    iter += 1;
  }
  double timeused = (double) (finish - start)/(CLOCKS_PER_SEC );

  //Calculate Analytical eigenvalues
  vec Exact = my_functions.CalculateExact(d,a,n);

  //Write results
  my_functions.WriteIter(iter);
  my_functions.WriteMeanError(Eigval, Exact, n);
  my_functions.WriteTime(timeused);

  //Print results
  my_functions.PrintResults(Eigval, Exact, R, n, jacobi_solve, armadillo_solve, arg);
return 0;
}


/*
all: compile execute

compile:
	c++ -o Jacobi_solver.exe Jacobi_solver.cpp JacobiMethod.cpp -larmadillo

execute:
	./Jacobi_solver.exe
*/





//
