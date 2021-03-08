# Table of Contents
- [Table of Contents](#table-of-contents)
- [Intro](#intro)
- [Storage](#storage)
  - [Vector](#vector)
  - [General](#general)
  - [BGeneral](#bgeneral)
  - [Triangle](#triangle)
  - [BTriangle](#btriangle)
  - [PTriangle](#ptriangle)
  - [Symmetric](#symmetric)
  - [BSymmetric](#bsymmetric)
  - [PSymmetric](#psymmetric)
  - [Hermitian](#hermitian)
  - [BHermitian](#bhermitian)
  - [PHermitian](#phermitian)
  - [Posdef](#posdef)
- [Operations](#operations)
  - [Givens](#givens)
  - [Rotate](#rotate)
  - [Swap](#swap)
  - [Scale](#scale)
  - [Copy](#copy)
  - [Axpy](#axpy)
  - [Norm2](#norm2)
  - [ASum](#asum)
  - [IAMax](#iamax)
  - [ExDot](#exdot)
  - [Dot](#dot)
  - [Solve](#solve)
  - [Update](#update)
  - [Inverse](#inverse)
  - [Determinant](#determinant)
  - [LUFact](#lufact)
  - [Cholesky](#cholesky)
  - [QRFact](#qrfact)
  - [Eigen](#eigen)
  - [Schur](#schur)
  - [LSquares](#lsquares)
  - [SVD](#svd)
  - [Rank](#rank)
- [Examples](#examples)

# Intro
You can include ```BLASW``` with

``` c++
#include <blasw/blasw.h>
```

The API is in ```Blasw``` namespace and consists of some wrappers around different types of matrix storage and functions that take these wrappers in and do some operation. The right way to use this documentation is to first choose the operation you want to perform and then find it in [the operations section](#operations). When you find it you can see what data types and matrix storages it supports. Then you look for a way to wrap your data pointer in [the storage scheme section](#storage). You can look at [the examples sections](#examples) to get familiar with the API.

# Storage
Each class that will be introduced here has some functions that create instance of the class. They usually optionally take stride as argument. If stride is not given, Then padding will be zero. You can see technical details of each storage type [here](https://www.netlib.org/lapack/lug/node121.html).

## Vector
1D vector.

+ ```Blasw::vec(Type *data, Size size, Size stride=0)```: wraps vector.

## General
General row or column major matrix.

+ ```Blasw::rmat(Type *data, Size rows, Size cols, Size stride=0)```: wraps row major general matrix.
+ ```Blasw::cmat(Type *data, Size rows, Size cols, Size stride=0)```: wraps column major general matrix.

## BGeneral
General banded row or column major matrix.

+ ```Blasw::rbmat(Type *data, Size rows, Size cols, Size sub, Size super, Size stride=0)```: wraps banded row major general matrix.
+ ```Blasw::cbmat(Type *data, Size rows, Size cols, Size sub, Size super, Size stride=0)```: wraps banded column major general matrix.

## Triangle
Upper or lower triangular row or column major matrix.

+ ```Blasw::rupper(Type *data, Size size, Size stride=0)```: wraps upper row major matrix.
+ ```Blasw::rlower(Type *data, Size size, Size stride=0)```: wraps lower row major matrix.
+ ```Blasw::cupper(Type *data, Size size, Size stride=0)```: wraps upper column major matrix.
+ ```Blasw::clower(Type *data, Size size, Size stride=0)```: wraps lower column major matrix.

## BTriangle
Upper or lower triangular banded row or column major matrix.

+ ```Blasw::rbupper(Type *data, Size size, Size stride=0)```: wraps upper banded row major matrix.
+ ```Blasw::rblower(Type *data, Size size, Size stride=0)```: wraps lower banded row major matrix.
+ ```Blasw::cbupper(Type *data, Size size, Size stride=0)```: wraps upper banded column major matrix.
+ ```Blasw::cblower(Type *data, Size size, Size stride=0)```: wraps lower banded column major matrix.

## PTriangle
Upper or lower triangular banded row or column major matrix.

+ ```Blasw::rpupper(Type *data, Size size, Size stride=0)```: wraps upper packed row major matrix.
+ ```Blasw::rplower(Type *data, Size size, Size stride=0)```: wraps lower packed row major matrix.
+ ```Blasw::cpupper(Type *data, Size size, Size stride=0)```: wraps upper packed column major matrix.
+ ```Blasw::cplower(Type *data, Size size, Size stride=0)```: wraps lower packed column major matrix.

## Symmetric
Upper or lower triangular row or column major symmetric matrix.

+ ```Blasw::rusym(Type *data, Size size, Size stride=0)```: wraps upper row major symmetric matrix.
+ ```Blasw::rlsym(Type *data, Size size, Size stride=0)```: wraps lower row major symmetric matrix.
+ ```Blasw::cusym(Type *data, Size size, Size stride=0)```: wraps upper column major symmetric matrix.
+ ```Blasw::clsym(Type *data, Size size, Size stride=0)```: wraps lower column major symmetric matrix.

## BSymmetric
Upper or lower triangular banded row or column major symmetric matrix.

+ ```Blasw::rbusym(Type *data, Size size, Size stride=0)```: wraps upper banded row major symmetric matrix.
+ ```Blasw::rblsym(Type *data, Size size, Size stride=0)```: wraps lower banded row major symmetric matrix.
+ ```Blasw::cbusym(Type *data, Size size, Size stride=0)```: wraps upper banded column major symmetric matrix.
+ ```Blasw::cblsym(Type *data, Size size, Size stride=0)```: wraps lower banded column major symmetric matrix.

## PSymmetric
Upper or lower triangular banded row or column major symmetric matrix.

+ ```Blasw::rpusym(Type *data, Size size, Size stride=0)```: wraps upper packed row major symmetric matrix.
+ ```Blasw::rplsym(Type *data, Size size, Size stride=0)```: wraps lower packed row major symmetric matrix.
+ ```Blasw::cpusym(Type *data, Size size, Size stride=0)```: wraps upper packed column major symmetric matrix.
+ ```Blasw::cplsym(Type *data, Size size, Size stride=0)```: wraps lower packed column major symmetric matrix.

## Hermitian
Upper or lower triangular row or column major hermitian matrix.

+ ```Blasw::ruherm(Type *data, Size size, Size stride=0)```: wraps upper row major hermitian matrix.
+ ```Blasw::rlherm(Type *data, Size size, Size stride=0)```: wraps lower row major hermitian matrix.
+ ```Blasw::cuherm(Type *data, Size size, Size stride=0)```: wraps upper column major hermitian matrix.
+ ```Blasw::clherm(Type *data, Size size, Size stride=0)```: wraps lower column major hermitian matrix.

## BHermitian
Upper or lower triangular banded row or column major hermitian matrix.

+ ```Blasw::rbuherm(Type *data, Size size, Size stride=0)```: wraps upper banded row major hermitian matrix.
+ ```Blasw::rblherm(Type *data, Size size, Size stride=0)```: wraps lower banded row major hermitian matrix.
+ ```Blasw::cbuherm(Type *data, Size size, Size stride=0)```: wraps upper banded column major hermitian matrix.
+ ```Blasw::cblherm(Type *data, Size size, Size stride=0)```: wraps lower banded column major hermitian matrix.

## PHermitian
Upper or lower triangular banded row or column major hermitian matrix.

+ ```Blasw::rpuherm(Type *data, Size size, Size stride=0)```: wraps upper packed row major hermitian matrix.
+ ```Blasw::rplherm(Type *data, Size size, Size stride=0)```: wraps lower packed row major hermitian matrix.
+ ```Blasw::cpuherm(Type *data, Size size, Size stride=0)```: wraps upper packed column major hermitian matrix.
+ ```Blasw::cplherm(Type *data, Size size, Size stride=0)```: wraps lower packed column major hermitian matrix.

## Posdef
Row or major positive definite matrix.

+ ```Blasw::rupod(Type *data, Size size, Size stride=0)```: wraps upper row major positive definite matrix.
+ ```Blasw::rlpod(Type *data, Size size, Size stride=0)```: wraps lower row major positive definite matrix.
+ ```Blasw::cupod(Type *data, Size size, Size stride=0)```: wraps upper column major positive definite matrix.
+ ```Blasw::clpod(Type *data, Size size, Size stride=0)```: wraps lower column major positive definite matrix.

# Operations
The above classes are template classes so they accept any data type but the actual operations in ```BLAS``` and ```LAPACK``` support single precision, double precision and complex single precision and complex double precision data types. Type names used in the operations:

+ F: ```float```
+ D: ```double```
+ CF: ```std::complex<float>```
+ DF: ```std::complex<double>```

For example ```<F, D>``` means that data type of variable can be float or double. If you want to know what the actual operations do, check out the [documentation of ```BLAS```](http://www.netlib.org/blas/#_blas_routines).

## Givens
Setup givens or modified givens rotation, uses ```rotg, rotmg```.

+ ```Blasw::givens(<F, D> &A, <F, D> &B, <F, D> &C, <F, D> &S)```
+ ```Blasw::givens(<F, D> &D1, <F, D> &D2, <F, D> &B1, <F, D> &B2, <F, D> P[5])```

## Rotate
Apply givens or modified givens rotation, uses ```rotm```.

+ ```Blasw::rotate(Vector<F, D> X, Vector<F, D> Y, <F, D> C, <F, D> S)```
+ ```Blasw::rotate(Vector<F, D> X, Vector<F, D> Y, <F, D> P[5])```

## Swap
Swap two vectors, uses ```swap```.

+ ```Blasw::swap(Vector<F, D, CF, DF> X, Vector<F, D, CF, DF> Y)```

## Scale
X = A * X, uses ```scal```.

+ ```Blasw::scale(Vector<F, D, CF, DF> X, <F, D> alpha)```

## Copy
Copies one vector into another, uses ```copy```.

+ ```Blasw::copy(Vector<F, D, CF, DF> X, Vector<F, D, CF, DF> Y)```

## Axpy
Y = A * X + Y, uses```axpy```.

+ ```Blasw::axpy(Vector<F, D, CF, DF> X, Vector<F, D, CF, DF> Y, <F, D, CF, DF> A)```

## Norm2
Euclidean norm, uses ```nrm2```.

+ ```Blasw::norm2(Vector<F, D, CF, CD> X)```

## ASum
Sum of absolute values of vector, uses ```asum```.

+ ```Blasw::asum(Vector<F, D, CF, CD> X)```

## IAMax
 index of max abs values of vector, uses ```iamax```.

+ ```Blasw::iamax(Vector<F, D, CF, CD> X)```

## ExDot
Vector dot product with extended precision accumulation, uses ```dsdot```.

+ ```Blasw::dot(Vector<F> X, Vector<F> Y)```
+ ```Blasw::dot(Vector<F> X, Vector<F> Y, F alpha)```

## Dot
Matrices and vectors multiplication.

+ ```Blasw::dot(Vector<F, D> X, Vector<F, D> Y)```
+ ```Blasw::dot(Vector<CF, CD> X, Vector<CF, CD> Y, bool conjugate=true)```
+ ```Blasw::dot(General<F, D, CF, CD> A, Vector<F, D, CF, CD> X, Vector<F, D, CF, CD> Y, <F, D, CF, CD> alpha, <F, D, CF, CD> beta)```
+ ```Blasw::dot(BGeneral<F, D, CF, CD> A, Vector<F, D, CF, CD> X, Vector<F, D, CF, CD> Y, <F, D, CF, CD> alpha, <F, D, CF, CD> beta)```
+ ```Blasw::dot(Symmetric<F, D> A, Vector<F, D> X, Vector<F, D> Y, <F, D> alpha, <F, D> beta)```
+ ```Blasw::dot(Hermitian<CF, CD> A, Vector<CF, CD> X, Vector<CF, CD> Y, <CF, CD> alpha, <CF, CD> beta)```
+ ```Blasw::dot(BSymmetric<F, D> A, Vector<F, D> X, Vector<F, D> Y, <F, D> alpha, <F, D> beta)```
+ ```Blasw::dot(BHermitian<CF, CD> A, Vector<CF, CD> X, Vector<CF, CD> Y, <CF, CD> alpha, <CF, CD> beta)```
+ ```Blasw::dot(PSymmetric<F, D> A, Vector<F, D> X, Vector<F, D> Y, <F, D> alpha, <F, D> beta)```
+ ```Blasw::dot(PHermitian<CF, CD> A, Vector<CF, CD> X, Vector<CF, CD> Y, <CF, CD> alpha, <CF, CD> beta)```
+ ```Blasw::dot(Triangle<F, D, CF, CD> X, Vector<F, D, CF, CD> Y)```
+ ```Blasw::dot(BTriangle<F, D, CF, CD> X, Vector<F, D, CF, CD> Y)```
+ ```Blasw::dot(PTriangle<F, D, CF, CD> X, Vector<F, D, CF, CD> Y)```
+ ```Blasw::dot(General<F, D, CF, CD> A, General<F, D, CF, CD> B, General<F, D, CF, CD> C, <F, D, CF, CD> alpha, <F, D, CF, CD> beta)```
+ ```Blasw::dot(Symmetric<F, D, CF, CD> A, General<F, D, CF, CD> B, General<F, D, CF, CD> C, <F, D, CF, CD> alpha, <F, D, CF, CD> beta)```
+ ```Blasw::dot(General<F, D, CF, CD> A, Symmetric<F, D, CF, CD> B, General<F, D, CF, CD> C, <F, D, CF, CD> alpha, <F, D, CF, CD> beta)```
+ ```Blasw::dot(Hermitian<F, D, CF, CD> A, General<F, D, CF, CD> B, General<F, D, CF, CD> C, <F, D, CF, CD> alpha, <F, D, CF, CD> beta)```
+ ```Blasw::dot(General<F, D, CF, CD> A, Hermitian<F, D, CF, CD> B, General<F, D, CF, CD> C, <F, D, CF, CD> alpha, <F, D, CF, CD> beta)```
+ ```Blasw::dot(Triangle<F, D, CF, CD> A, General<F, D, CF, CD> B, <F, D, CF, CD> alpha)```

## Solve
Solve triangular matrix problems.

+ ```Blasw::solve(Triangle<F, D, CF, CD> X, Vector<F, D, CF, CD> Y)```
+ ```Blasw::solve(BTriangle<F, D, CF, CD> X, Vector<F, D, CF, CD> Y)```
+ ```Blasw::solve(PTriangle<F, D, CF, CD> X, Vector<F, D, CF, CD> Y)```
+ ```Blasw::solve(Triangle<F, D, CF, CD> X, General<F, D, CF, CD> B, <F, D, CF, CD> alpha)```
+ ```Blasw::solve(General<F, D, CF, CD> X, Triangle<F, D, CF, CD> B, <F, D, CF, CD> alpha)```
+ ```Blasw::solve(General<F, D, CF, CD> A, General<F, D, CF, CD> B)```
+ ```Blasw::solve(Symmetric<F, D, CF, CD> A, General<F, D, CF, CD> B)```
+ ```Blasw::solve(Hermitian<CF, CD> A, General<CF, CD> B)```

## Update
Perform rank 1 or 2 update

+ ```Blasw::update(Vector<F, D> X, Vector<F, D> Y, General<F, D> A, <F, D> alpha)```
+ ```Blasw::update(Vector<CF, DF> X, Vector<CF, DF> Y, General<CF, DF> A, <CF, CD> alpha, bool conj = true)```
+ ```Blasw::update(Vector<F, D> X, Symmetric<F, D> A, <F, D> alpha)```
+ ```Blasw::update(Vector<CF, CD> X, Hermitian<CF, CD> A, <F, D> alpha)```
+ ```Blasw::update(Vector<F, D> X, PSymmetric<F, D> A, <F, D> alpha)```
+ ```Blasw::update(Vector<CF, CD> X, PHermitian<CF, CD> A, <F, D> alpha)```
+ ```Blasw::update(Vector<F, D> X, Vector<F, D> Y, Symmetric<F, D> A, <F, D> alpha)```
+ ```Blasw::update(Vector<CF, CD> X, Vector<CF, CD> Y, Hermitian<CF, CD> A, <CF, CD> alpha)```
+ ```Blasw::update(General<F, D, CF, CD> A, Symmetric<F, D, CF, CD> C, <F, D, CF, CD> alpha, <F, D, CF, CD> beta)```
+ ```Blasw::update(General<CF, CD> A, Hermitian<CF, CD> C, <F, D> alpha, <F, D> beta)```
+ ```Blasw::update(General<F, D, CF, CD> A, General<F, D, CF, CD> B, Symmetric<F, D, CF, CD> C, <F, D, CF, CD> alpha, <F, D, CF, CD> beta)```
+ ```Blasw::update(General<CF, CD> A, General<CF, CD> B, Hermitian<T> C, <CF, CD> alpha, <F, D> beta)```

## Inverse
Inverse of matrix

+ ```Blasw::inverse(General<F, D, CF, CD> A)```
+ ```Blasw::inverse(Symmetric<F, D, CF, CD> A)```

## Determinant
Determinant of matrix

+ ```Blasw::determinant(General<F, D, CF, CD> A)```
+ ```Blasw::determinant(Symmetric<F, D, CF, CD> A)```
+ ```Blasw::determinant(Hermitian<CF, CD> A)```

## LUFact
Lower upper factorization

+ ```Blasw::lufact(General<F, D, CF, CD> A, Vector<int> P)```
+ ```Blasw::lufact(Symmetric<F, D, CF, CD> A, Vector<int> P)```
+ ```Blasw::lufact(Hermitian<CF, CD> A, Vector<int> P)```

## Cholesky
cholesky factorization

+ ```Blasw::cholesky(Posdef<F, D, CF, CD> A)```

## QRFact
QR Factorization

+ ```Blasw::qrfact(General<F, D, CF, CD> A, Vector<F, D, CF, CD> T)```

## Eigen
Eigen decomposition

+ ```Blasw::eigen(General<F, D> A, Vector<CF, CD> E, General<F, D> L, General<F, D> R)```
+ ```Blasw::eigen(General<CF, CD> A, Vector<CF, CD> E, General<CF, CD> L, General<CF, CD> R)```
+ ```Blasw::eigen(Symmetric<F, D> A, Vector<F, D> W, bool vectors)```
+ ```Blasw::eigen(Hermitian<CF, CD> A, Vector<F, D> W, bool vectors)```

## Schur
Schur decomposition

+ ```Blasw::schur(General<F, D> A, Vector<CF, CD> E, General<F, D> V)```
+ ```Blasw::schur(General<CF, CD> A, Vector<CF, CD> E, General<CF, CD> V)```

## LSquares
Least Squares solution

+ ```Blasw::lsquares(General<F, D, CF, CD> A, General<F, D, CF, CD> B)```

## SVD
SVD decomposition

+ ```Blasw::svd(General<F, D, CF, CD> A, Vector<F, D> S, General<F, D, CF, CD> U, General<F, D, CF, CD> VT)```

## Rank
Matrix rank

+ ```Blasw::rank(General<F, D, CF, CD> A, F, D, CF, CD epsilon)```

# Examples
Axpy Example:

``` c++
Blasw::Size N = 100;
auto X = new float[N];
auto Y = new float[N];
float alpha = 2;

for (Blasw::Size i = 0; i < N; ++i) X[i] = 1;
Blasw::axpy(Blasw::vec(X, N), Blasw::vec(Y, N), alpha);
```

General Matrix Multiplication Example:

``` c++
Blasw::Size X = 10, Y = 100, Z = 20;
auto A = new float[X * Y];
auto B = new float[Y * Z];
auto C = new float[X * Z];
float alpha = 1, beta = 0;

for (Blasw::Size i = 0; i < X * Y; ++i) A[i] = 1;
for (Blasw::Size i = 0; i < Y * Z; ++i) A[i] = 2;
Blasw::dot(Blasw::rmat(A, X, Y), Blasw::rmat(B, Y, Z), Blasw::rmat(C, X, Z), alpha, beta);
```
