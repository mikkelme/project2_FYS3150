//  Diagonalizing tridiagonal Toeplitz matrix  with Lapack functions
//  Compile as c++ -O3 -o Tridiag.exe TridiagToeplitz.cpp -larmadillo -llapack -lblas

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>

using namespace  std;
using namespace  arma;


int main(int argc, char* argv[]){

  int     i, j, N;
  double  Rmin, Rmax, h, hh, d, a;
  Rmin = 0; Rmax = 1; N = 20;

  //Intregration step length
  h = (Rmax - Rmin)/N; hh = h*h;
  d = 2.0 / (hh);   // DiagConst
  a =  -1.0 / (hh); //NondiagConst

  // Setting up tridiagonal matrix A
  mat A = zeros<mat>(N,N);
  A(0,0) = d;
  A(0,1) = a;
  for(i = 1; i < N-1; i++) {
    A(i,i-1) = a;
    A(i,i)   = d;
    A(i,i+1) = a;
  }
  A(N-1,N-2) = a;
  A(N-1,N-1) = d;

  // Diagonalize and obtain eigenvalues
  vec Eigval(N);
  eig_sym(Eigval, A);

  // Print numerical vs analytical results
  double pi = acos(-1.0);
  cout << "RESULTS:" << endl;
  cout << setiosflags(ios::showpoint | ios::uppercase);
  cout <<"Number of Eigenvalues = " << setw(15) << N << endl;
  cout << "Difference: abs(analytical - numerical)" << endl;
  for(int i = 0; i < N; i++) {
    double Exact = d+2*a*cos((i+1)*pi/(N+1));
    cout << setw(15) << setprecision(8) << fabs(Eigval[i]-Exact) << endl;
  }

  return 0;
}  //  end of main function
