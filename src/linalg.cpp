#include "linalg.hpp"
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <random>

/*////////////////////////////////////////////////////////////////////
   MyVector class.
 */////////////////////////////////////////////////////////////////////

// Constructors.
MyVector::MyVector(int n, double x) : v(n, x) {
}
MyVector::MyVector() {
}
MyVector::MyVector(int n) : v(n) {
}

// Size of vector.
int MyVector::Size() const {
        return v.size();
}

// Access elements.
double &MyVector::operator[](int index) {
        if (index < 0 || index > v.size() - 1) {
                std::cout << "Check your vector indexing!" << std::endl;
                exit(1);
        }
        return v[index];
}

double MyVector::operator[](int index) const {
        if (index < 0 || index > v.size() - 1) {
                std::cout << "Check your vector indexing!" << std::endl;
                exit(1);
        }
        return v[index];
}

// L-infinity norm.
double MyVector::LInfNorm() const {
        double mx = 0.0;
        for (double const &val : v) {
                mx = std::max(std::abs(val), mx);
        }
        return mx;
}
// L2 norm.
double MyVector::L2Norm() const {
        double sum = 0;
        for (double const &val : v) {
                sum += val * val;
        }
        return std::sqrt(sum);
}

void MyVector::MakeNormal(double mu, double sigma) {
        std::random_device rnd;
        std::mt19937 mersenne_engine(rnd());
        std::normal_distribution<double> dist(mu, sigma);
        auto gen = std::bind(dist, mersenne_engine);
        generate(begin(v), end(v), gen); // Fill v with normally distributed numbers
}

// Overload the << operator to display a MyVector
std::ostream &operator<<(std::ostream &os, const MyVector &V) {
        int n = V.Size();
        os << "(";
        for (int i = 0; i < n - 1; i++) {
                os << V[i] << ", ";
        }
        os << V[n - 1] << ")" << std::endl;
        return os;
}

/*////////////////////////////////////////////////////////////////////
   MyMatrix class.
 */////////////////////////////////////////////////////////////////////

// Constructors.
MyMatrix::MyMatrix() : rows(0), cols(0) {
}
MyMatrix::MyMatrix(int m, int n, double x) : rows(m), cols(n), A(m * n, x) {
}

// Set all elements to input value.
MyMatrix &MyMatrix::operator=(double x) {
        for (int i = 0; i < rows * cols; i++) {
                A[i] = x;
        }
        return *this;
}

double MyMatrix::operator()(int i, int j) const {
        if (i < 0 || i > rows - 1 || j < 0 || j > cols - 1) {
                std::cout << "Check your matrix indexing!" << std::endl;
                exit(1);
        }
        return A[j + i * cols];
}

double &MyMatrix::operator()(int i, int j) {
        if (i < 0 || i > rows - 1 || j < 0 || j > cols - 1) {
                std::cout << "Check your matrix indexing!" << std::endl;
                exit(1);
        }
        return A[j + i * cols];
}

int MyMatrix::nRows() const {
        return rows;
}
int MyMatrix::nCols() const {
        return cols;
}

MyMatrix MyMatrix::Transpose() const {
        MyMatrix B(cols, rows);
        for (int i = 0; i < cols; i++) {
                for (int j = 0; j < rows; j++) {
                        B(i, j) = A[i + j * cols];
                }
        }
        return B;
}

void MyMatrix::MakeNormal(double mu, double sigma) {
        std::random_device rnd;
        std::mt19937 mersenne_engine(rnd());
        std::normal_distribution<double> dist(mu, sigma);
        auto gen = std::bind(dist, mersenne_engine);
        generate(begin(A), end(A), gen); // Fill A with normally distributed numbers.
}

// Overload the << operator to display a MyMatrix
std::ostream &operator<<(std::ostream &os, const MyMatrix &M) {
        for (int i = 0; i < M.nRows(); i++) {
                for (int j = 0; j < M.nCols(); j++) {
                        os.width(15);
                        os << M(i, j);
                }
                os << std::endl;
        }
        return os;
}

/*////////////////////////////////////////////////////////////////////
   MyBandedMatrix class.
*/////////////////////////////////////////////////////////////////////

// Constructors.
MyBandedMatrix::MyBandedMatrix() : rows(0), cols(0) {
}
MyBandedMatrix::MyBandedMatrix(int n, int m, int lband, int rband, double x)
        : rows(n), cols(m), A(n * (lband + rband + 1), x), l(lband), r(rband) {
}

// Access element (rvalue).
double MyBandedMatrix::operator()(int i, int j) const {
        if (i < 0 || i > rows - 1 || j < 0 || j > cols - 1) {
                std::cout << "Check your matrix indexing!" << std::endl;
                exit(1);
        }
        if (j < i - l || j > i + r) {
                return 0; // returns zero if index is not in the band
        }
        return A[(j + l - i) + i * (l + r + 1)];
}

// Access element (lvalue).
double &MyBandedMatrix::operator()(int i, int j) {
        if (i < 0 || i > rows - 1 || j < 0 || j > cols - 1) {
                std::cout << "Check your matrix indexing!" << std::endl;
                exit(1);
        }
        if (j < i - l || j > i + r) {
                return z; // returns z = 0 if index is not in the band
        }
        return A[(j + l - i) + i * (l + r + 1)];
}

int MyBandedMatrix::nRows() const {
        return rows;
}
int MyBandedMatrix::nCols() const {
        return cols;
}
int MyBandedMatrix::Bands() const {
        return r + l + 1;
}                                                       // total number of bands
int MyBandedMatrix::LBands() const {
        return l;
}                                                       // number of left bands
int MyBandedMatrix::RBands() const {
        return r;
}                                                       // number of right bands

// Overload the << operator to display a MBandedMatrix
std::ostream &operator<<(std::ostream &output, const MyBandedMatrix &banded) {
        for (int i = 0; i < banded.nRows();
             i++) { // calculate position of lower and upper band
                int jmin = std::max(std::min(i - banded.LBands(), banded.nCols()), 0);
                int jmax = std::min(i + banded.RBands() + 1, banded.nCols());
                output << "( ";
                for (int j = 0; j < jmin; j++)
                        output << 0 << "\t ";
                for (int j = jmin; j < jmax; j++)
                        output << banded(i, j) << "\t ";
                for (int j = jmax; j < banded.nCols(); j++)
                        output << 0 << "\t ";
                output << ")\n";
        }
        return output;
}

/*////////////////////////////////////////////////////////////////////
   Functions.
 */////////////////////////////////////////////////////////////////////

// Dot product of two MyVectors.
double dot(const MyVector &lhs, const MyVector &rhs) {
        if (lhs.Size() != rhs.Size()) {
                std::cout << "Error! Those two vectors aren't the same size!" << std::endl;
                exit(1);
        }
        double d = 0;
        for (int i = 0; i < lhs.Size(); i++) {
                d += (lhs[i] * rhs[i]);
        }
        return d;
}

// Generate a Poisson matrix.
MyMatrix Poisson(int n) {
        MyMatrix P(n, n); // note all entries are initialised to zero.
        for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                        if (i == j) {
                                P(i, j) = 2;
                        } else if (std::abs(i - j) == 1) {
                                P(i, j) = -1;
                        }
                }
        }
        return P;
}

// Banded Poisson matrix.
MyBandedMatrix BandedPoisson(int n) {
        MyBandedMatrix B(n, n, 1, 1); // all entries initialised to zero
        for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                        if (i == j) {
                                B(i, j) = 2;
                        } else if (std::abs(i - j) == 1) {
                                B(i, j) = -1;
                        }
                }
        }
        return B;
}

// Gradient descent.
int GradientDescent(const MyMatrix &A, const MyVector &b, MyVector &x,
                    int maxIterations = 1000, double tol = 1e-6) {
        int iterations = -1; // Returns -1 if method does not converge.
        MyMatrix T = A.Transpose();
        MyVector r = T * (b - A * x);

        for (int i = 0; i < maxIterations; i++) {
                double alpha =
                        dot(r, r) /
                        dot(r, T * A * r); // dot product rather than vector transpose.
                x = x + alpha * r;
                r = r - alpha * (T * A * r); // no operator for double * matrix.
                if (r.L2Norm() < tol) {
                        iterations = i + 1;
                        break;
                }
        }
        return iterations;
}

// Conjugate gradient method.
MyVector ConjugateGradient(const MyMatrix& A, MyVector &x, const MyVector& b, int maxIterations = 1000, double tol = 1e-6){

        MyVector r = b - A * x;
        MyVector p = r;

        for (int i = 0; i < maxIterations; i++) {
                double alpha = dot(r, r) / dot(p, A * p);
                x = x + alpha * p; // Note += operator not overloaded.
                MyVector r1 = r; // Needed to calculate beta.
                r = r - alpha * (A * p);
                if (r.L2Norm() < tol) { break;}
                double beta = dot(r, r) / dot(r1, r1);
                p = r + beta * p;
        }
        return x;
}
