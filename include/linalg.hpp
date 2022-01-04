#ifndef MYVECTOR_H
#define MYVECTOR_H

#ifndef MYMATRIX_H
#define MYMATRIX_H

#ifndef MYBANDEDMATRIX_H
#define MYBANDEDMATRIX_H

#include <iostream>
#include <vector>

// Vector class.
class MyVector {
public:
  // Three constructors...
  MyVector();
  explicit MyVector(int n);
  MyVector(int n, double x);

  // Get the size of the vector.
  int Size() const;

  // Access elements.
  double &operator[](int index);
  double operator[](int index) const;

  // Norms.
  double LInfNorm() const;
  double L2Norm() const;

  // Make normal.
  void MakeNormal(double mu, double sigma);

private:
  std::vector<double> v;
};
#endif

// Matrix class.
class MyMatrix {
public:
  MyMatrix();
  MyMatrix(int m, int n, double x = 0.0);

  // Sets all matrix entries to an input double.
  MyMatrix &operator=(double x);

  // Access matrix elements.
  double &operator()(int i, int j);
  double operator()(int i, int j) const;

  // Get number of rows and columns.
  int nRows() const;
  int nCols() const;

  // Transpose.
  MyMatrix Transpose() const;

  // Fill with random normally distributed entries.
  void MakeNormal(double mu, double sigma);

private:
  int rows, cols;
  std::vector<double> A;
};

#endif

class MyBandedMatrix {
public:
  MyBandedMatrix();
  MyBandedMatrix(int n, int m, int lband, int rband, double x = 0);

  // access element [rvalue]
  double operator()(int i, int j) const;
  // access element [lvalue]
  double &operator()(int i, int j);

  int nRows() const;
  int nCols() const;

  int Bands() const;  // total number of bands
  int LBands() const; // number of left bands
  int RBands() const; // number of right bands

private:
  int rows, cols;
  std::vector<double> A;
  int l, r;     // these are the number of left/right diagonals
  double z = 0; // needed because &operator() must return an lvalue
};

#endif

/////////////////////////
// Functions.
////////////////////////

// Overload the << operator to display a MyVector.
std::ostream &operator<<(std::ostream &os, const MyVector &V);
// Overload the << operator to display a MyMatrix.
std::ostream &operator<<(std::ostream &os, const MyMatrix &M);
// Overload the << operator to display a MyBandedMatrix.
std::ostream &operator<<(std::ostream &output, const MyBandedMatrix &banded);

// Dot product of two MyVectors.
double dot(const MyVector &lhs, const MyVector &rhs);

// Generate a Poisson matrix.
MyMatrix Poisson(int n);
// Generate a banded Posson matrix.
MyBandedMatrix BandedPoisson(int n);

// Gradient descent.
int GradientDescent(const MyMatrix &A, const MyVector &b, MyVector &x,
                    int maxIterations, double tol);

// Conjugate gradient.
MyVector ConjugateGradient(const MyMatrix &A, MyVector &x, const MyVector &b,
                           int maxIterations, double tol);

// scalar * MyVector.
inline MyVector operator*(const double &lhs, const MyVector &rhs) {
  MyVector temp(rhs);
  for (int i = 0; i < temp.Size(); i++) {
    temp[i] *= lhs;
  }
  return temp;
}

// MyVector * scalar.
inline MyVector operator*(const MyVector &lhs, const double &rhs) {
  MyVector temp(lhs);
  for (int i = 0; i < temp.Size(); i++) {
    temp[i] *= rhs;
  }
  return temp;
}

// MyVector + MyVector.
inline MyVector operator+(const MyVector &lhs, const MyVector &rhs) {
  if (lhs.Size() != rhs.Size()) {
    std::cout << "Error! Those two vectors aren't the same size!" << std::endl;
    exit(1);
  }
  MyVector temp(rhs);
  for (int i = 0; i < temp.Size(); i++) {
    temp[i] += lhs[i];
  }
  return temp;
}

// MyVector - MyVector.
inline MyVector operator-(const MyVector &lhs, const MyVector &rhs) {
  if (lhs.Size() != rhs.Size()) {
    std::cout << "Error! Those two vectors aren't the same size!" << std::endl;
    exit(1);
  }
  MyVector temp(lhs);
  for (int i = 0; i < temp.Size(); i++) {
    temp[i] -= rhs[i];
  }
  return temp;
}

// Overload * for MyVector * MyMatrix.
inline MyVector operator*(const MyMatrix &A, const MyVector &x) {
  if (A.nCols() != x.Size()) {
    std::cout
        << "That MyMatrix and that MyVector are incompatible! Check the sizes!"
        << std::endl;
    exit(1);
  }

  MyVector u(A.nRows());
  for (int i = 0; i < A.nRows(); i++) {
    for (int j = 0; j < A.nCols(); j++) {
      u[i] += A(i, j) * x[j];
    }
  }
  return u;
}

// Overload * operator for MyMatrix * MyMatrix (i.e., matrix multiplication).
// This is the simplest possible way to do it; there are better algorithms but
// they're more difficult to code.
inline MyMatrix operator*(const MyMatrix &A, const MyMatrix &B) {
  if (A.nCols() != B.nRows()) {
    std::cout << "Error! Those two matrices can't be multiplied!" << std::endl;
    exit(1);
  }

  MyMatrix C(A.nRows(), B.nCols()); // C has the same number of rows as A and
                                    // the same number of columns as B.

  for (int i = 0; i < C.nRows(); i++) {
    for (int j = 0; j < C.nCols(); j++) {
      for (int k = 0; k < B.nRows(); k++) {
        C(i, j) += A(i, k) * B(k, j);
      }
    }
  }
  return C;
}

// MyBandedMatrix * MyVector.
inline MyVector operator*(const MyBandedMatrix &B, const MyVector &x) {
  if (B.nCols() != x.Size()) {
    std::cout << "That MyBandedMatrix and that MyVector are incompatible! "
                 "Check the sizes!"
              << std::endl;
    exit(1);
  }

  MyVector y(B.nRows(), 0);
  for (int i = 0; i < y.Size(); i++) {
    for (int j = 0; j < B.nCols(); j++) {
      if (j < i - B.LBands() || j > i + B.RBands()) {
        continue;
      } // we ignore all the zero elements outside the band
      y[i] += B(i, j) * x[j];
    }
  }
  return y;
}
