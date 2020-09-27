#include "JacobiMethod.h"
#include <iostream>
#include <cmath>
#include <armadillo>
#include "time.h"

using namespace  std;
using namespace  arma;


int main(int argc, char *argv[]) {
  Jacobi my_functions;
  //Create Matrix
  int n, arg;
  double  Rmin, Rmax, h, hh, d, a, rho_max, w;
  
  n = atoi(argv[1]);
  Rmin = 0; Rmax = 1;
  //Intregration step length
  h = (Rmax - Rmin)/n; hh = h*h;
  d = 2.0 / (hh);   // DiagConst
  a =  -1.0 / (hh); //nondiagConst
  arg = atoi(argv[2]);
  // if arg == 1, one electron potential, if arg == 2, two electron potential
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

  double tol = 1.0E-10;
  int iter = 1;
  int maxiter = 5e4;
  double maxnondig;
  int p; int q;

  my_functions.OrthTest(tol);
  my_functions.EigValTest(tol);
  my_functions.MaxOffTest(tol);
  cout << "All tests passed." << endl;

  clock_t start, finish;
  start = clock();
  
  my_functions.MaxOffdiag(A, p, q, maxnondig, n);
  while (maxnondig > tol && iter <= maxiter){
    my_functions.JacobiRotate(A, R, p, q, n);
    my_functions.MaxOffdiag(A, p, q, maxnondig, n);
    iter++;
  }

  finish = clock();
  double timeused = (double) (finish - start)/(CLOCKS_PER_SEC );

  //Order result and write to files
  mat Eigval = my_functions.OrderEigenResults(A, R, n);
  my_functions.WriteIter(iter - 1);
  my_functions.WriteMeanError(Eigval, d, a, n);
  my_functions.WriteEig(Eigval, R, n);
  my_functions.WriteTime(timeused);


return 0;
}











//
