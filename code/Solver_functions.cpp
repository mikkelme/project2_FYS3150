#include "Solver_functions.h"
#include <iostream>
#include <cmath>
#include <iomanip>
#include <armadillo>
#include <cassert>
#include <vector>

using namespace  std;
using namespace  arma;

mat Jacobi::CreateTridiagonal(double d, double a, int n, int arg, double rho_max, double w){
  //Create tridiagonal matr ix
  int     i, j;
  double rho[n];

  //rho[0] = 0.0;
  rho[n] = rho_max;
  double h = rho[n]/n;
  for(i=0; i<n; i++){
  	rho[i] = i*h;
    rho[i] = rho[i]*rho[i];
  }

  if (arg != 1 || arg != 2){
    for(i = 0; i<n-1; i++){
      rho[i] = 0;
    }
  }
  if (arg == 1){
    for(i = 0; i<n-1; i++){
      rho[i] = i*h;
      rho[i] = rho[i]*rho[i];
    }
  }
  if (arg == 2){
    h = (rho[n]-rho[0])/n;
    for(i = 0; i<n-1; i++){
      rho[i] = i*h;
      rho[i] = w*w*rho[i]*rho[i] + 1/rho[i];
    }
  }

  mat A = zeros<mat>(n,n);
  A(0,0) = d;
  A(0,1) = a;
  for(i = 1; i < n-1; i++) {
    A(i,i-1) = a;
    A(i,i)   = d+rho[i];
    A(i,i+1) = a;
  }
  A(n-1,n-2) = a;
  A(n-1,n-1) = d;
  return A;
}
void Jacobi::ShowMatrix(mat A){
  //Prints a matrix
  int ni = size(A)[0];
  int nj = size(A)[1];
  for (int i = 0; i < ni; i++){
    for (int j = 0; j < nj; j++){
        cout << setw(4) << A(i,j) << " ";
    }
    cout << endl;
  }
}
void Jacobi::MaxOffdiag(mat A, int& p, int& q, double& maxnondig, int n){
  //Find max off-diagonal element
  //Update indexes p and q and return absolute value
  double max;
  for (int i = 0; i < n; ++i){
    for (int j= i + 1; j < n; ++j){
      double aij = fabs(A(i,j));
      if (aij > max)
      {
        max = aij; p = i; q = j;
      }
    }
  }
  maxnondig = fabs(A(p,q));
}
void Jacobi::JacobiRotate(mat& A, mat& Eigvec, int k, int l, int n){
  //Perform jacobi raotation
  double tau, tan, cos, sin;

  if (A(k,l) != 0.0){ //avoid divison by zero
    tau = (A(l,l) - A(k,k))/(2*A(k,l));

    //Choose right sign for tau calculaiton
    if (tau >= 0){
      tan = 1.0/(tau + sqrt(tau*tau + 1.0)); }
    else{
      tan = -1.0/(-tau + sqrt(tau*tau + 1.0)); }

    cos = 1.0/sqrt(1.0 + tan*tan);
    sin = cos*tan;
  }
  else {
    tan = 0.0;
    cos = 1.0;
  }

  //Perform similarity transformation A_new = (S^T)(A)(S)
  double a_kk_old = A(k,k);
  double a_ll_old = A(l,l);
  A(k,k) = a_kk_old*cos*cos - 2.0*A(k,l)*cos*sin + a_ll_old*sin*sin;
  A(l,l) = a_ll_old*cos*cos + 2.0*A(k,l)*cos*sin + a_kk_old*sin*sin;
  A(k,l) = 0.0; //Hard coding to zero
  A(l,k) = 0.0; //--------||--------
  for (int i = 0; i < n; ++i){
    if (i != k && i != l){
      double a_ik_old = A(i,k);
      double a_il_old = A(i,l);
      A(i,k) = a_ik_old*cos - a_il_old*sin;
      A(k,i) = A(i,k);
      A(i,l) = a_il_old*cos + a_ik_old*sin;
      A(l,i) = A(i,l);
    }
    //eigenvectors
    double r_ik_old = Eigvec(i,k);
    double r_il_old = Eigvec(i,l);

    Eigvec(i,k) = r_ik_old*cos - r_il_old*sin;
    Eigvec(i,l) = r_il_old*cos + r_ik_old*sin;
  }
return;
}
vec Jacobi::OrderEigenResults(mat& A, mat& Eigvec, int n){
  //Return sorted array for the eigenvalues
  //& sort the eigenvector-matrix accordingly
  vec Eigval = zeros<mat>(n);
  for (int i = 0; i < n; i++){
    Eigval(i) = A(i,i);
  }
  uvec idx = sort_index(Eigval);
  Eigvec = Eigvec.rows(idx(span(0, n-1)));
  return sort(Eigval);

}
pair<vec, mat> Jacobi::CalculateExact(double d, double a, int n){
  //Calculate analytical eigenvalues
  vec Exact_eigval = zeros<mat>(n);
  mat Exact_eigvec = zeros <mat>(n,n);

  double pi = acos(-1.0);
  for(int i = 0; i < n; i++){
    Exact_eigval(i) = d+2*a*cos((i+1)*pi/(n+1));
    for (int j = 0; j < n; j++){
      Exact_eigvec(i,j) = sin((j+1)*(i+1)*pi/(n+1));
    }
  }
  return make_pair(Exact_eigval, Exact_eigvec);
}
void Jacobi::WriteIter(int iter){
  ofstream ofile;
  string output_file = "iter_dim.txt";
  ofile.open(output_file, ios::out | ios::app);
  if (ofile.fail()){
    throw ios_base::failure(strerror(errno));
  }
  ofile.exceptions(ofile.exceptions() | ios::failbit | ifstream::badbit);

  ofile << iter << endl;
  ofile.close();
}
void Jacobi::WriteMeanError(vec& Eigval, vec& Exact_eigval, int n){
  ofstream ofile;
  string output_file = "MeanError.txt";
  double MeanError;

  ofile.open(output_file, ios::out | ios::app);
  if (ofile.fail()){
    throw ios_base::failure(strerror(errno));
  }
  ofile.exceptions(ofile.exceptions() | ios::failbit | ifstream::badbit);
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  for(int i = 0; i < n; i++) {
    MeanError += fabs(Exact_eigval(i) - Eigval(i));
  }
  MeanError /= n;
  ofile << setprecision(8) << MeanError <<endl;
  ofile.close();
}
void Jacobi::PrintResults(vec& Eigval, mat& Eigvec, vec& Exact_eigval, mat& Exact_eigvec, int n, bool jacobi_solve, bool armadillo_solve, int arg){
  string solver; string potential;
  if (jacobi_solve)         {solver = "Solver = Jacobi Method";}
  else if (armadillo_solve) {solver = "Solver = Armadillo eig_sym";}
  if (arg == 1)             {potential = "Potential = one electron";}
  else if (arg == 2)        {potential = "Potential = two electrons";}
  else                      {potential = "Potential = None ";}
  cout << "---------------------------" << endl;
  cout << "Program details:" << endl;
  cout << solver << endl;
  cout << "Dimensions (n) = "<< n << endl;
  cout << potential << endl;
  cout << "---------------------------\n" << endl;

  cout << setw(11) << "Eigenvalues";
  cout << setw(20) << "Numerical";
  cout << setw(20) << "Analytical";
  cout << setw(20) << "Absolute Error" << endl;
  for (int i = 0; i < n; i++){
    cout << setw(11) << i+1;
    cout << setw(20) << Eigval(i);
    cout << setw(20) << Exact_eigval(i);
    cout << setw(20) << fabs(Exact_eigval(i)-Eigval(i)) << endl;
  }
  cout << "\nNumerical eigenvectors:" << endl;
  for (int i = 0; i < n; i++){
    cout << i+1 << "):  [ ";
    for (int j = 0; j < n; j++){
      cout << Eigvec(i,j) << " ";
    }
    cout << "]"<< endl;
  }
  cout << endl;

  cout << "\nAnalytical eigenvectors:" << endl;
  for (int i = 0; i < n; i++){
    cout << i+1 << "):  [ ";
    for (int j = 0; j < n; j++){
      cout << Exact_eigvec(i,j) << " ";
    }
    cout << "]"<< endl;
  }
  cout << endl;
}
void Jacobi::WriteTime(double timeused){
  ofstream ofile;
  string output_file = "Timeused.txt";

  ofile.open(output_file, ios::out | ios::app);
  if (ofile.fail()){
    throw ios_base::failure(strerror(errno));
  }
  ofile.exceptions(ofile.exceptions() | ios::failbit | ifstream::badbit);
  ofile << setiosflags(ios::showpoint | ios::uppercase);

  ofile << setprecision(8) << timeused <<endl;
  ofile.close();
}
void Jacobi::OrthTest(double tol){
	cout << "Testing orthogonality..." << endl;

	int n = 4;
	mat A(n,n, fill::randu); // Random 4x4 matrix
	mat Eigvec = eye<mat>(n,n);

	int p, q;
	int iter = 1;
	int maxiter = 1000;
	double maxnondig;

	Jacobi::MaxOffdiag(A, p, q, maxnondig, n);
	while (maxnondig > tol && iter <= maxiter){
	    Jacobi::JacobiRotate(A, Eigvec, p, q, n);
	    Jacobi::MaxOffdiag(A, p, q, maxnondig, n);
	    iter++;
	}

	mat RT = Eigvec.t();
	mat I = eye<mat>(n,n);

	bool test = false;
	mat testMatrix = abs(RT*Eigvec-I);
	if (all(all(testMatrix < tol))){
		test = true;
	}

	assert(test);
	cout << "Test passed successufully." << endl;
}
void Jacobi::EigValTest(double tol){
	cout << "Testing eigenvalues..." << endl;

	int n = 4;
	mat A = {{2, 1, 0, 0}, // Symmetric 4x4 tridiagonal matrix
			 {1, 2, 1, 0},
			 {0, 1, 2, 1},
			 {0, 0, 1, 2}};
	mat Eigvec = eye<mat>(n,n);

	int p, q;
	int iter = 1;
	int maxiter = 1000;
	double maxnondig;

	vec lambda(n);
	eig_sym(lambda,A);

	Jacobi::MaxOffdiag(A, p, q, maxnondig, n);
	while (maxnondig > tol && iter <= maxiter){
	    Jacobi::JacobiRotate(A, Eigvec, p, q, n);
	    Jacobi::MaxOffdiag(A, p, q, maxnondig, n);
	    iter++;
	}
	vec rotate_lambda = {A(0,0), A(1,1), A(2,2), A(3,3)};

	bool test = false;
	if (all(sort(lambda)- sort(rotate_lambda) < tol)){
		test = true;
	}

	assert(test);
	cout << "Test passed successufully." << endl;
}
void Jacobi::MaxOffTest(double tol){
	cout << "Testing max off-diagonal..." << endl;

	int n = 4;
	mat A = {{1, 2, 5, 3}, // Some matrix where an off-diagonal element is the largest = 5	int p, q;
			 {2, 3, 4, 1},
			 {1, 0, 2, 3},
			 {3, -2, 3, 2}};
	mat Eigvec = eye<mat>(n,n);


  int p; int q;
  double maxnondig;

	Jacobi::MaxOffdiag(A, p, q, maxnondig, n);
  int arma_max = A.max();

	bool test = false;
	if (maxnondig == 5 && maxnondig == arma_max){
		test = true;
	}

	assert(test);
	cout << "Test passed successufully." << endl;
}
