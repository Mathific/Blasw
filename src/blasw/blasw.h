/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2021, Shahriar Rezghi <shahriar25.ss@gmail.com>
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#pragma once

#include <complex>

#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>

#include <blasw/config.h>

#ifdef BLASW_CBLAS_MKL
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif

#ifdef BLASW_LAPACKE_FOUND
#ifdef BLASW_LAPACKE_MKL
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif
#endif

#include <cstdlib>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>

#define CHECK(expr) \
    if (!static_cast<bool>(expr)) throw std::runtime_error(std::string("Expression \"") + #expr + "\" has failed!");

namespace Blasw
{
using Size = int;
using Index = CBLAS_INDEX;

enum class Major
{
    Row = CblasRowMajor,
    Col = CblasColMajor,
};
enum class State
{
    None = CblasNoTrans,
    Trans = CblasTrans,
    ConjTrans = CblasConjTrans,
};
enum class Triangular
{
    Upper = CblasUpper,
    Lower = CblasLower,
};
enum class Diagonal
{
    Unit = CblasUnit,
    NonUnit = CblasNonUnit,
};

template <typename T>
struct From
{
    using Type = T;
};
template <typename T>
struct From<std::complex<T>>
{
    using Type = T;
};
template <typename T>
struct To
{
    using Type = std::complex<T>;
};
template <typename T>
struct To<std::complex<T>>
{
    using Type = std::complex<T>;
};

namespace Impl
{
inline CBLAS_ORDER blcvt(Major major) { return CBLAS_ORDER(major); }
inline CBLAS_TRANSPOSE blcvt(State state) { return CBLAS_TRANSPOSE(state); }
inline CBLAS_UPLO blcvt(Triangular tri) { return CBLAS_UPLO(tri); }
inline CBLAS_DIAG blcvt(Diagonal diag) { return CBLAS_DIAG(diag); }
}  // namespace Impl

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

#define STATE(name)               \
    name<T> &trans()              \
    {                             \
        state = State::Trans;     \
        return *this;             \
    }                             \
    name<T> &adjoint()            \
    {                             \
        state = State::ConjTrans; \
        return *this;             \
    }

template <typename T>
struct Vector
{
    T *data;
    Size size, stride;

    Vector() {}
    Vector(T *data, Size size, Size stride) : data(data), size(size) { this->stride = (stride == 0 ? 1 : stride); }
};

template <typename T>
Vector<T> vec(T *data, Size size, Size stride = 0)
{
    return Vector<T>(data, size, stride);
}

template <typename T>
struct General
{
    using Type = T;

    T *data;
    Size rows, cols, stride;
    Major major;
    State state;
    STATE(General);

    General() : data(nullptr), rows(0), cols(0), stride(0), major(Major::Row), state(State::None) {}
    General(T *data, Size rows, Size cols, Size stride, Major major, State state)
        : data(data), rows(rows), cols(cols), major(major), state(state)
    {
        this->stride = (stride == 0 ? (major == Major::Row ? cols : rows) : stride);
    }

    bool _trans() const { return state == State::Trans || state == State::ConjTrans; }
    Size _rows() const { return _trans() ? cols : rows; }
    Size _cols() const { return _trans() ? rows : cols; }
};

#define GENERAL(name, major)                                                    \
    template <typename T>                                                       \
    General<T> name(T *data, Size rows, Size cols, Size stride = 0)             \
    {                                                                           \
        return General<T>(data, rows, cols, stride, Major::major, State::None); \
    }
GENERAL(rmat, Row) GENERAL(cmat, Col);

template <typename T>
struct BGeneral
{
    using Type = T;

    T *data;
    Size rows, cols, stride;
    Major major;
    State state;
    Size sub, super;
    STATE(BGeneral);

    BGeneral() : data(nullptr), rows(0), cols(0), stride(0), major(Major::Row), state(State::None), sub(0), super(0) {}
    BGeneral(T *data, Size rows, Size cols, Size stride, Major major, State state, Size sub, Size super)
        : data(data), rows(rows), cols(cols), major(major), state(state), sub(sub), super(super)
    {
        this->stride = (stride == 0 ? (major == Major::Row ? cols : rows) : stride);
    }

    bool _trans() const { return state == State::Trans || state == State::ConjTrans; }
    Size _rows() const { return _trans() ? cols : rows; }
    Size _cols() const { return _trans() ? rows : cols; }
};

#define BGENERAL(name, major, str)                                                           \
    template <typename T>                                                                    \
    BGeneral<T> name(T *data, Size rows, Size cols, Size sub, Size super, Size stride = 0)   \
    {                                                                                        \
        return BGeneral<T>(data, rows, cols, stride, Major::major, State::None, sub, super); \
    }
BGENERAL(rbmat, Row, cols) BGENERAL(cbmat, Col, rows);

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <typename T>
struct Triangle
{
    using Type = T;

    T *data;
    Size size, stride;
    Major major;
    State state;
    Triangular tri;
    Diagonal diag;
    STATE(Triangle);

    Triangle()
        : data(nullptr),
          size(0),
          stride(0),
          major(Major::Row),
          state(State::None),
          tri(Triangular::Upper),
          diag(Diagonal::NonUnit)
    {
    }
    Triangle(T *data, Size size, Size stride, Major major, State state, Triangular tri, Diagonal diag)
        : data(data), size(size), major(major), state(state), tri(tri), diag(diag)
    {
        this->stride = (stride == 0 ? size : stride);
    }
    Triangle<T> &unit()
    {
        diag = Diagonal::Unit;
        return *this;
    }
};

#define TRIANGLE(name, major, tri)                                                          \
    template <typename T>                                                                   \
    Triangle<T> name(T *data, Size size, Size stride = 0)                                   \
    {                                                                                       \
        return Triangle<T>(data, size, stride, Major::major, State::None, Triangular::tri); \
    }
TRIANGLE(rupper, Row, Upper) TRIANGLE(cupper, Col, Upper);
TRIANGLE(rlower, Row, Lower) TRIANGLE(clower, Col, Lower);

template <typename T>
struct BTriangle
{
    using Type = T;

    T *data;
    Size size, stride;
    Major major;
    State state;
    Triangular tri;
    Diagonal diag;
    Size super;
    STATE(BTriangle);

    BTriangle()
        : data(nullptr),
          size(0),
          stride(0),
          major(Major::Row),
          state(State::None),
          tri(Triangular::Upper),
          diag(Diagonal::NonUnit),
          super(0)
    {
    }
    BTriangle(T *data, Size size, Size stride, Major major, State state, Triangular tri, Diagonal diag, Size super)
        : data(data), size(size), major(major), state(state), tri(tri), diag(diag), super(super)
    {
        this->stride = (stride == 0 ? size : stride);
    }

    BTriangle<T> &unit()
    {
        diag = Diagonal::Unit;
        return *this;
    }
};

#define BTRIANGLE(name, major, tri)                                                                 \
    template <typename T>                                                                           \
    BTriangle<T> name(T *data, Size size, Size super, Size stride = 0)                              \
    {                                                                                               \
        return BTriangle<T>(data, size, stride, Major::major, State::None, Triangular::tri, super); \
    }
BTRIANGLE(rbupper, Row, Upper) BTRIANGLE(cbupper, Col, Upper);
BTRIANGLE(rblower, Row, Lower) BTRIANGLE(cblower, Col, Lower);

template <typename T>
struct PTriangle
{
    using Type = T;

    T *data;
    Size size;
    Major major;
    State state;
    Triangular tri;
    Diagonal diag;
    STATE(PTriangle);

    PTriangle()
        : data(nullptr), size(0), major(Major::Row), state(State::None), tri(Triangular::Upper), diag(Diagonal::NonUnit)
    {
    }
    PTriangle(T *data, Size size, Major major, State state, Triangular tri, Diagonal diag)
        : data(data), size(size), major(major), state(state), tri(tri), diag(diag)
    {
    }

    PTriangle<T> &unit()
    {
        diag = Diagonal::Unit;
        return *this;
    }
};

#define PTRIANGLE(name, major, tri)                                                  \
    template <typename T>                                                            \
    PTriangle<T> name(T *data, Size size)                                            \
    {                                                                                \
        return PTriangle<T>(data, size, Major::major, State::None, Triangular::tri); \
    }
PTRIANGLE(rpupper, Row, Upper) PTRIANGLE(cpupper, Col, Upper);
PTRIANGLE(rplower, Row, Lower) PTRIANGLE(cplower, Col, Lower);

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <typename T>
struct Symmetric
{
    using Type = T;

    T *data;
    Size size, stride;
    Major major;
    Triangular tri;

    Symmetric() : data(nullptr), size(0), stride(0), major(Major::Row), tri(Triangular::Upper) {}
    Symmetric(T *data, Size size, Size stride, Major major, Triangular tri)
        : data(data), size(size), major(major), tri(tri)
    {
        this->stride = (stride == 0 ? size : stride);
    }
};

#define SYMMETRIC(name, major, tri)                                             \
    template <typename T>                                                       \
    Symmetric<T> name(T *data, Size size, Size stride = 0)                      \
    {                                                                           \
        return Symmetric<T>(data, size, stride, Major::major, Triangular::tri); \
    }
SYMMETRIC(rusym, Row, Upper) SYMMETRIC(cusym, Col, Upper);
SYMMETRIC(rlsym, Row, Lower) SYMMETRIC(clsym, Col, Lower);

template <typename T>
struct BSymmetric
{
    using Type = T;

    T *data;
    Size size, stride;
    Major major;
    Triangular tri;
    Size super;

    BSymmetric() : data(nullptr), size(0), stride(0), major(Major::Row), tri(Triangular::Upper), super(0) {}
    BSymmetric(T *data, Size size, Size stride, Major major, Triangular tri, Size super)
        : data(data), size(size), major(major), tri(tri), super(super)
    {
        this->stride = (stride == 0 ? size : stride);
    }
};

#define BSYMMETRIC(name, major, tri)                                                    \
    template <typename T>                                                               \
    BSymmetric<T> name(T *data, Size size, Size super, Size stride = 0)                 \
    {                                                                                   \
        return BSymmetric<T>(data, size, stride, Major::major, Triangular::tri, super); \
    }
BSYMMETRIC(rbusym, Row, Upper) BSYMMETRIC(cbusym, Col, Upper);
BSYMMETRIC(rblsym, Row, Lower) BSYMMETRIC(cblsym, Col, Lower);

template <typename T>
struct PSymmetric
{
    using Type = T;

    T *data;
    Size size;
    Major major;
    Triangular tri;

    PSymmetric() : data(nullptr), size(0), major(Major::Row), tri(Triangular::Upper) {}
    PSymmetric(T *data, Size size, Major major, Triangular tri) : data(data), size(size), major(major), tri(tri) {}
};

#define PSYMMETRIC(name, major, tri)                                     \
    template <typename T>                                                \
    PSymmetric<T> name(T *data, Size size)                               \
    {                                                                    \
        return PSymmetric<T>(data, size, Major::major, Triangular::tri); \
    }
PSYMMETRIC(rpusym, Row, Upper) PSYMMETRIC(cpusym, Col, Upper);
PSYMMETRIC(rplsym, Row, Lower) PSYMMETRIC(cplsym, Col, Lower);

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <typename T>
struct Hermitian
{
    using Type = T;

    T *data;
    Size size, stride;
    Major major;
    Triangular tri;

    Hermitian() : data(nullptr), size(0), stride(0), major(Major::Row), tri(Triangular::Upper) {}
    Hermitian(T *data, Size size, Size stride, Major major, Triangular tri)
        : data(data), major(major), size(size), tri(tri)
    {
        this->stride = (stride == 0 ? size : stride);
    }
};

#define HERMITIAN(name, major, tri)                                             \
    template <typename T>                                                       \
    Hermitian<T> name(T *data, Size size, Size stride = 0)                      \
    {                                                                           \
        return Hermitian<T>(data, size, stride, Major::major, Triangular::tri); \
    }
HERMITIAN(ruherm, Row, Upper) HERMITIAN(cuherm, Col, Upper);
HERMITIAN(rlherm, Row, Lower) HERMITIAN(clherm, Col, Lower);

template <typename T>
struct BHermitian
{
    using Type = T;

    T *data;
    Size size, stride;
    Major major;
    Triangular tri;
    Size super;

    BHermitian() : data(nullptr), size(0), stride(0), major(Major::Row), tri(Triangular::Upper), super(0) {}
    BHermitian(T *data, Size size, Size stride, Major major, Triangular tri, Size super)
        : data(data), size(size), major(major), tri(tri), super(super)
    {
        this->stride = (stride == 0 ? size : stride);
    }
};

#define BHERMITIAN(name, major, tri)                                                    \
    template <typename T>                                                               \
    BHermitian<T> name(T *data, Size size, Size super, Size stride = 0)                 \
    {                                                                                   \
        return BHermitian<T>(data, size, stride, Major::major, Triangular::tri, super); \
    }
BHERMITIAN(rbuherm, Row, Upper) BHERMITIAN(cbuherm, Col, Upper);
BHERMITIAN(rblherm, Row, Lower) BHERMITIAN(cblherm, Col, Lower);

template <typename T>
struct PHermitian
{
    using Type = T;

    T *data;
    Size size;
    Major major;
    Triangular tri;

    PHermitian() : data(nullptr), size(0), major(Major::Row), tri(Triangular::Upper) {}
    PHermitian(T *data, Size size, Major major, Triangular tri) : data(data), size(size), major(major), tri(tri) {}
};

#define PHERMITIAN(name, major, tri)                                     \
    template <typename T>                                                \
    PHermitian<T> name(T *data, Size size)                               \
    {                                                                    \
        return PHermitian<T>(data, size, Major::major, Triangular::tri); \
    }
PHERMITIAN(rpusym, Row, Upper) PHERMITIAN(cpusym, Col, Upper);
PHERMITIAN(rplsym, Row, Lower) PHERMITIAN(cplsym, Col, Lower);

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <typename T>
struct Posdef
{
    using Type = T;

    T *data;
    Size size, stride;
    Major major;
    Triangular tri;

    Posdef() : data(nullptr), size(0), stride(0), major(Major::Row), tri(Triangular::Upper) {}
    Posdef(T *data, Size size, Size stride, Major major, Triangular tri)
        : data(data), size(size), major(major), tri(tri)
    {
        this->stride = (stride == 0 ? size : stride);
    }
};

#define POSDEF(name, major, tri)                                             \
    template <typename T>                                                    \
    Posdef<T> name(T *data, Size size, Size stride = 0)                      \
    {                                                                        \
        return Posdef<T>(data, size, stride, Major::major, Triangular::tri); \
    }
POSDEF(rupod, Row, Upper) POSDEF(cupod, Col, Upper);
POSDEF(rlpod, Row, Lower) POSDEF(clpod, Col, Lower);

////////////////////////////////////////////////////////////////////////////////
/// BLAS LEVEL 1 ///////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

#define REPEAT0(func, name1, name2) \
    func(cblas_##name1, float, );   \
    func(cblas_##name2, double, );

#define REPEAT1(func, name1, name2)                      \
    func(cblas_##name1, std::complex<float>, (void *)&); \
    func(cblas_##name2, std::complex<double>, (void *)&);

#define REPEAT(func, name1, name2, name3, name4) \
    REPEAT0(func, name1, name2)                  \
    REPEAT1(func, name3, name4)

#define ROTG(F, T, R) \
    inline void givens(T &a, T &b, T &c, T &s) { F(&a, &b, &c, &s); }

#define ROTMG(F, T, R) \
    inline void givens(T &d1, T &d2, T &b1, T &b2, T P[5]) { F(&d1, &d2, &b1, b2, P); }

REPEAT0(ROTG, srotg, drotg)
REPEAT0(ROTMG, srotmg, drotmg)

#define ROT(F, T, R)                                                               \
    inline void rotate(Vector<T> X, Vector<T> Y, From<T>::Type c, From<T>::Type s) \
    {                                                                              \
        CHECK(X.size == Y.size);                                                   \
        return F(X.size, X.data, X.stride, Y.data, Y.stride, c, s);                \
    }

REPEAT0(ROT, srot, drot)

#define ROTM(F, T, R)                                                \
    inline void rotate(Vector<T> X, Vector<T> Y, From<T>::Type P[5]) \
    {                                                                \
        CHECK(X.size == Y.size);                                     \
        return F(X.size, X.data, X.stride, Y.data, Y.stride, P);     \
    }

REPEAT0(ROTM, srotm, drotm)

#define SWAP(F, T, R)                                         \
    inline void swap(Vector<T> X, Vector<T> Y)                \
    {                                                         \
        CHECK(X.size == Y.size);                              \
        return F(X.size, X.data, X.stride, Y.data, Y.stride); \
    }

REPEAT(SWAP, sswap, dswap, cswap, zswap)

#define SCAL(F, T, R) \
    inline void scale(Vector<T> X, T alpha) { return F(X.size, R alpha, X.data, X.stride); }

REPEAT(SCAL, sscal, dscal, cscal, zscal)

#define SSCAL(F, T, R) \
    inline void scale(Vector<T> X, From<T>::Type alpha) { return F(X.size, alpha, X.data, X.stride); }

REPEAT1(SSCAL, csscal, zdscal)

#define COPY(F, T, R)                                         \
    inline void copy(Vector<T> X, Vector<T> Y)                \
    {                                                         \
        CHECK(X.size == Y.size);                              \
        return F(X.size, X.data, X.stride, Y.data, Y.stride); \
    }

REPEAT(COPY, scopy, dcopy, ccopy, zcopy)

#define AXPY(F, T, R)                                                  \
    inline void axpy(Vector<T> X, Vector<T> Y, T alpha)                \
    {                                                                  \
        CHECK(X.size == Y.size);                                       \
        return F(X.size, R alpha, X.data, X.stride, Y.data, Y.stride); \
    }

REPEAT(AXPY, saxpy, daxpy, caxpy, zaxpy)

#define DOT(F, T, R)                                          \
    inline T dot(Vector<T> X, Vector<T> Y)                    \
    {                                                         \
        CHECK(X.size == Y.size);                              \
        return F(X.size, X.data, X.stride, Y.data, Y.stride); \
    }

REPEAT0(DOT, sdot, ddot)

#define DOTU(F, T, R)                                                       \
    inline T dot(Vector<T> X, Vector<T> Y, bool conj = true)                \
    {                                                                       \
        CHECK(X.size == Y.size);                                            \
        T result;                                                           \
        if (conj)                                                           \
            F##c_sub(X.size, X.data, X.stride, Y.data, Y.stride, R result); \
        else                                                                \
            F##u_sub(X.size, X.data, X.stride, Y.data, Y.stride, R result); \
        return result;                                                      \
    }

REPEAT1(DOTU, cdot, zdot)

inline float exdot(Vector<float> X, Vector<float> Y, float beta)
{
    CHECK(X.size == Y.size);
    return cblas_sdsdot(X.size, beta, X.data, X.stride, Y.data, Y.stride);
}
inline float exdot(Vector<float> X, Vector<float> Y)
{
    CHECK(X.size == Y.size);
    return cblas_dsdot(X.size, X.data, X.stride, Y.data, Y.stride);
}

#define NRM2(F, T, R) \
    inline From<T>::Type norm2(Vector<T> X) { return F(X.size, X.data, X.stride); }

REPEAT(NRM2, snrm2, dnrm2, scnrm2, dznrm2)

#define ASUM(F, T, R) \
    inline From<T>::Type asum(Vector<T> X) { return F(X.size, X.data, X.stride); }

REPEAT(ASUM, sasum, dasum, scasum, dzasum)

#define IAMAX(F, T, R) \
    inline Index iamax(Vector<T> X) { return F(X.size, X.data, X.stride); }

REPEAT(IAMAX, isamax, idamax, icamax, izamax)

////////////////////////////////////////////////////////////////////////////////
/// BLAS LEVEL 2 ///////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

#define GEMV(F, T, R)                                                                                              \
    inline void dot(General<T> A, Vector<T> X, Vector<T> Y, T alpha, T beta)                                       \
    {                                                                                                              \
        CHECK(A._cols() == X.size) CHECK(A._rows() == Y.size);                                                     \
        F(Impl::blcvt(A.major), Impl::blcvt(A.state), A.rows, A.cols, R alpha, A.data, A.stride, X.data, X.stride, \
          R beta, Y.data, Y.stride);                                                                               \
    }

REPEAT(GEMV, sgemv, dgemv, cgemv, zgemv)

#define GBMV(F, T, R)                                                                                            \
    inline void dot(BGeneral<T> A, Vector<T> X, Vector<T> Y, T alpha, T beta)                                    \
    {                                                                                                            \
        CHECK(A._cols() == X.size) CHECK(A._rows() == Y.size);                                                   \
        F(Impl::blcvt(A.major), Impl::blcvt(A.state), A.rows, A.cols, A.sub, A.super, R alpha, A.data, A.stride, \
          X.data, X.stride, R beta, Y.data, Y.stride);                                                           \
    }

REPEAT(GBMV, sgbmv, dgbmv, cgbmv, zgbmv)

#define SYMV(F, T, R)                                                                                            \
    inline void dot(Symmetric<T> A, Vector<T> X, Vector<T> Y, T alpha, T beta)                                   \
    {                                                                                                            \
        CHECK(A.size == X.size) CHECK(A.size == Y.size);                                                         \
        F(Impl::blcvt(A.major), Impl::blcvt(A.tri), A.size, R alpha, A.data, A.stride, X.data, X.stride, R beta, \
          Y.data, Y.stride);                                                                                     \
    }

#define HEMV(F, T, R)                                                                                            \
    inline void dot(Hermitian<T> A, Vector<T> X, Vector<T> Y, T alpha, T beta)                                   \
    {                                                                                                            \
        CHECK(A.size == X.size) CHECK(A.size == Y.size);                                                         \
        F(Impl::blcvt(A.major), Impl::blcvt(A.tri), A.size, R alpha, A.data, A.stride, X.data, X.stride, R beta, \
          Y.data, Y.stride);                                                                                     \
    }

REPEAT0(SYMV, ssymv, dsymv)
REPEAT1(HEMV, chemv, zhemv)

#define SBMV(F, T, R)                                                                                             \
    inline void dot(BSymmetric<T> A, Vector<T> X, Vector<T> Y, T alpha, T beta)                                   \
    {                                                                                                             \
        CHECK(A.size == X.size) CHECK(A.size == Y.size);                                                          \
        F(Impl::blcvt(A.major), Impl::blcvt(A.tri), A.size, A.super, R alpha, A.data, A.stride, X.data, X.stride, \
          R beta, Y.data, Y.stride);                                                                              \
    }

#define HBMV(F, T, R)                                                                                             \
    inline void dot(BHermitian<T> A, Vector<T> X, Vector<T> Y, T alpha, T beta)                                   \
    {                                                                                                             \
        CHECK(A.size == X.size) CHECK(A.size == Y.size);                                                          \
        F(Impl::blcvt(A.major), Impl::blcvt(A.tri), A.size, A.super, R alpha, A.data, A.stride, X.data, X.stride, \
          R beta, Y.data, Y.stride);                                                                              \
    }

REPEAT0(SBMV, ssbmv, dsbmv)
REPEAT1(HBMV, chbmv, zhbmv)

#define SPMV(F, T, R)                                                                                          \
    inline void dot(PSymmetric<T> A, Vector<T> X, Vector<T> Y, T alpha, T beta)                                \
    {                                                                                                          \
        CHECK(A.size == X.size) CHECK(A.size == Y.size);                                                       \
        F(Impl::blcvt(A.major), Impl::blcvt(A.tri), A.size, R alpha, A.data, X.data, X.stride, R beta, Y.data, \
          Y.stride);                                                                                           \
    }

#define HPMV(F, T, R)                                                                                          \
    inline void dot(PHermitian<T> A, Vector<T> X, Vector<T> Y, T alpha, T beta)                                \
    {                                                                                                          \
        CHECK(A.size == X.size) CHECK(A.size == Y.size);                                                       \
        F(Impl::blcvt(A.major), Impl::blcvt(A.tri), A.size, R alpha, A.data, X.data, X.stride, R beta, Y.data, \
          Y.stride);                                                                                           \
    }

REPEAT0(SPMV, sspmv, dspmv)
REPEAT1(HPMV, chpmv, zhpmv)

#define TRMV(F, T, R)                                                                                          \
    inline void dot(Triangle<T> A, Vector<T> X)                                                                \
    {                                                                                                          \
        CHECK(A.size == X.size);                                                                               \
        F(Impl::blcvt(A.major), Impl::blcvt(A.tri), Impl::blcvt(A.state), Impl::blcvt(A.diag), A.size, A.data, \
          A.stride, X.data, X.stride);                                                                         \
    }

REPEAT(TRMV, strmv, dtrmv, ctrmv, ztrmv)

#define TBMV(F, T, R)                                                                                           \
    inline void dot(BTriangle<T> A, Vector<T> X)                                                                \
    {                                                                                                           \
        CHECK(A.size == X.size);                                                                                \
        F(Impl::blcvt(A.major), Impl::blcvt(A.tri), Impl::blcvt(A.state), Impl::blcvt(A.diag), A.size, A.super, \
          A.data, A.stride, X.data, X.stride);                                                                  \
    }

REPEAT(TBMV, stbmv, dtbmv, ctbmv, ztbmv)

#define TPMV(F, T, R)                                                                                                  \
    inline void dot(PTriangle<T> A, Vector<T> X)                                                                       \
    {                                                                                                                  \
        CHECK(A.size == X.size);                                                                                       \
        F(Impl::blcvt(A.major), Impl::blcvt(A.tri), Impl::blcvt(A.state), Impl::blcvt(A.diag), A.size, A.data, X.data, \
          X.stride);                                                                                                   \
    }

REPEAT(TPMV, stpmv, dtpmv, ctpmv, ztpmv)

#define TRSV(F, T, R)                                                                                          \
    inline void solve(Triangle<T> A, Vector<T> X)                                                              \
    {                                                                                                          \
        CHECK(A.size == X.size);                                                                               \
        F(Impl::blcvt(A.major), Impl::blcvt(A.tri), Impl::blcvt(A.state), Impl::blcvt(A.diag), A.size, A.data, \
          A.stride, X.data, X.stride);                                                                         \
    }

REPEAT(TRSV, strsv, dtrsv, ctrsv, ztrsv)

#define TBSV(F, T, R)                                                                                           \
    inline void solve(BTriangle<T> A, Vector<T> X)                                                              \
    {                                                                                                           \
        CHECK(A.size == X.size);                                                                                \
        F(Impl::blcvt(A.major), Impl::blcvt(A.tri), Impl::blcvt(A.state), Impl::blcvt(A.diag), A.size, A.super, \
          A.data, A.stride, X.data, X.stride);                                                                  \
    }

REPEAT(TBSV, stbsv, dtbsv, ctbsv, ztbsv)

#define TPSV(F, T, R)                                                                                                  \
    inline void solve(PTriangle<T> A, Vector<T> X)                                                                     \
    {                                                                                                                  \
        CHECK(A.size == X.size);                                                                                       \
        F(Impl::blcvt(A.major), Impl::blcvt(A.tri), Impl::blcvt(A.state), Impl::blcvt(A.diag), A.size, A.data, X.data, \
          X.stride);                                                                                                   \
    }

REPEAT(TPSV, stpsv, dtpsv, ctpsv, ztpsv)

#define GER(F, T, R)                                                                                            \
    inline void update(Vector<T> X, Vector<T> Y, General<T> A, T alpha)                                         \
    {                                                                                                           \
        CHECK(X.size == A.rows) CHECK(Y.size == A.cols) CHECK(A.state == State::None);                          \
        F(Impl::blcvt(A.major), X.size, Y.size, R alpha, X.data, X.stride, Y.data, Y.stride, A.data, A.stride); \
    }

REPEAT0(GER, sger, dger)

#define GERU(F, T, R)                                                                                                  \
    inline void update(Vector<T> X, Vector<T> Y, General<T> A, T alpha, bool conj = true)                              \
    {                                                                                                                  \
        CHECK(X.size == A.rows) CHECK(Y.size == A.cols) CHECK(A.state == State::None);                                 \
        if (conj)                                                                                                      \
            F##c(Impl::blcvt(A.major), X.size, Y.size, R alpha, X.data, X.stride, Y.data, Y.stride, A.data, A.stride); \
        else                                                                                                           \
            F##u(Impl::blcvt(A.major), X.size, Y.size, R alpha, X.data, X.stride, Y.data, Y.stride, A.data, A.stride); \
    }

REPEAT1(GERU, cger, zger)

#define SYR(F, T, R)                                                                                      \
    inline void update(Vector<T> X, Symmetric<T> A, T alpha)                                              \
    {                                                                                                     \
        CHECK(A.size == X.size);                                                                          \
        F(Impl::blcvt(A.major), Impl::blcvt(A.tri), A.size, R alpha, X.data, X.stride, A.data, A.stride); \
    }

#define HER(F, T, R)                                                                                    \
    inline void update(Vector<T> X, Hermitian<T> A, From<T>::Type alpha)                                \
    {                                                                                                   \
        CHECK(A.size == X.size);                                                                        \
        F(Impl::blcvt(A.major), Impl::blcvt(A.tri), A.size, alpha, X.data, X.stride, A.data, A.stride); \
    }

REPEAT0(SYR, ssyr, dsyr)
REPEAT1(HER, cher, zher)

#define SPR(F, T, R)                                                                            \
    inline void update(Vector<T> X, PSymmetric<T> A, T alpha)                                   \
    {                                                                                           \
        CHECK(A.size == X.size);                                                                \
        F(Impl::blcvt(A.major), Impl::blcvt(A.tri), A.size, R alpha, X.data, X.stride, A.data); \
    }

#define HPR(F, T, R)                                                                          \
    inline void update(Vector<T> X, PHermitian<T> A, From<T>::Type alpha)                     \
    {                                                                                         \
        CHECK(A.size == X.size);                                                              \
        F(Impl::blcvt(A.major), Impl::blcvt(A.tri), A.size, alpha, X.data, X.stride, A.data); \
    }

REPEAT0(SPR, sspr, dspr)
REPEAT1(HPR, chpr, zhpr)

#define SYR2(F, T, R)                                                                                            \
    inline void update(Vector<T> X, Vector<T> Y, Symmetric<T> A, T alpha)                                        \
    {                                                                                                            \
        CHECK(A.size == X.size) CHECK(X.size == Y.size);                                                         \
        F(Impl::blcvt(A.major), Impl::blcvt(A.tri), A.size, R alpha, X.data, X.stride, Y.data, Y.stride, A.data, \
          A.stride);                                                                                             \
    }

#define HER2(F, T, R)                                                                                            \
    inline void update(Vector<T> X, Vector<T> Y, Hermitian<T> A, T alpha)                                        \
    {                                                                                                            \
        CHECK(A.size == X.size) CHECK(X.size == Y.size);                                                         \
        F(Impl::blcvt(A.major), Impl::blcvt(A.tri), A.size, R alpha, X.data, X.stride, Y.data, Y.stride, A.data, \
          A.stride);                                                                                             \
    }

REPEAT0(SYR2, ssyr2, dsyr2)
REPEAT1(HER2, cher2, zher2)

#define SPR2(F, T, R)                                                                                             \
    inline void update(Vector<T> X, Vector<T> Y, PSymmetric<T> A, T alpha)                                        \
    {                                                                                                             \
        CHECK(A.size == X.size) CHECK(X.size == Y.size);                                                          \
        F(Impl::blcvt(A.major), Impl::blcvt(A.tri), A.size, R alpha, X.data, X.stride, Y.data, Y.stride, A.data); \
    }

#define HPR2(F, T, R)                                                                                             \
    inline void update(Vector<T> X, Vector<T> Y, PHermitian<T> A, T alpha)                                        \
    {                                                                                                             \
        CHECK(A.size == X.size) CHECK(X.size == Y.size);                                                          \
        F(Impl::blcvt(A.major), Impl::blcvt(A.tri), A.size, R alpha, X.data, X.stride, Y.data, Y.stride, A.data); \
    }

REPEAT0(SPR2, sspr2, dspr2)
REPEAT1(HPR2, chpr2, zhpr2)

////////////////////////////////////////////////////////////////////////////////
/// BLAS LEVEL 3 ///////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

#define GEMM(F, T, R)                                                                                           \
    inline void dot(General<T> A, General<T> B, General<T> C, T alpha, T beta)                                  \
    {                                                                                                           \
        CHECK(A.major == B.major) CHECK(A.major == C.major) CHECK(C.state == State::None);                      \
        CHECK(A._cols() == B._rows()) CHECK(A._rows() == C.rows) CHECK(B._cols() == C.cols);                    \
        F(Impl::blcvt(A.major), Impl::blcvt(A.state), Impl::blcvt(B.state), C.rows, C.cols, A._cols(), R alpha, \
          A.data, A.stride, B.data, B.stride, R beta, C.data, C.stride);                                        \
    }

REPEAT(GEMM, sgemm, dgemm, cgemm, zgemm)

#define SYMM(F, T, R)                                                                                              \
    inline void dot(Symmetric<T> A, General<T> B, General<T> C, T alpha, T beta)                                   \
    {                                                                                                              \
        CHECK(A.major == B.major) CHECK(A.major == C.major);                                                       \
        CHECK(B.state == State::None) CHECK(C.state == State::None);                                               \
        CHECK(A.size == B.rows) CHECK(A.size == C.rows) CHECK(B.cols == C.cols);                                   \
        F(Impl::blcvt(A.major), CblasLeft, Impl::blcvt(A.tri), C.rows, C.cols, R alpha, A.data, A.stride, B.data,  \
          B.stride, R beta, C.data, C.stride);                                                                     \
    }                                                                                                              \
    inline void dot(General<T> B, Symmetric<T> A, General<T> C, T alpha, T beta)                                   \
    {                                                                                                              \
        CHECK(A.major == B.major) CHECK(A.major == C.major);                                                       \
        CHECK(B.state == State::None) CHECK(C.state == State::None);                                               \
        CHECK(A.size == B.cols) CHECK(A.size == C.cols) CHECK(B.rows == C.rows);                                   \
        F(Impl::blcvt(A.major), CblasRight, Impl::blcvt(A.tri), C.rows, C.cols, R alpha, A.data, A.stride, B.data, \
          B.stride, R beta, C.data, C.stride);                                                                     \
    }

REPEAT(SYMM, ssymm, dsymm, csymm, zsymm)

#define HEMM(F, T, R)                                                                                              \
    inline void dot(Hermitian<T> A, General<T> B, General<T> C, T alpha, T beta)                                   \
    {                                                                                                              \
        CHECK(A.major == B.major) CHECK(A.major == C.major);                                                       \
        CHECK(B.state == State::None) CHECK(C.state == State::None);                                               \
        CHECK(A.size == B.rows) CHECK(A.size == C.rows) CHECK(B.cols == C.cols);                                   \
        F(Impl::blcvt(A.major), CblasLeft, Impl::blcvt(A.tri), C.rows, C.cols, R alpha, A.data, A.stride, B.data,  \
          B.stride, R beta, C.data, C.stride);                                                                     \
    }                                                                                                              \
    inline void dot(General<T> B, Hermitian<T> A, General<T> C, T alpha, T beta)                                   \
    {                                                                                                              \
        CHECK(A.major == B.major) CHECK(A.major == C.major);                                                       \
        CHECK(B.state == State::None) CHECK(C.state == State::None);                                               \
        CHECK(A.size == B.cols) CHECK(A.size == C.cols) CHECK(B.rows == C.rows);                                   \
        F(Impl::blcvt(A.major), CblasRight, Impl::blcvt(A.tri), C.rows, C.cols, R alpha, A.data, A.stride, B.data, \
          B.stride, R beta, C.data, C.stride);                                                                     \
    }

REPEAT1(HEMM, chemm, zhemm)

#define TRMM(F, T, R)                                                                                              \
    inline void dot(Triangle<T> A, General<T> B, T alpha)                                                          \
    {                                                                                                              \
        CHECK(A.major == B.major) CHECK(A.size == B.rows) CHECK(B.state == State::None);                           \
        F(Impl::blcvt(A.major), CblasLeft, Impl::blcvt(A.tri), Impl::blcvt(A.state), Impl::blcvt(A.diag), B.rows,  \
          B.cols, R alpha, A.data, A.stride, B.data, B.stride);                                                    \
    }                                                                                                              \
    inline void dot(General<T> B, Triangle<T> A, T alpha)                                                          \
    {                                                                                                              \
        CHECK(A.major == B.major) CHECK(A.size == B.cols) CHECK(B.state == State::None);                           \
        F(Impl::blcvt(A.major), CblasRight, Impl::blcvt(A.tri), Impl::blcvt(A.state), Impl::blcvt(A.diag), B.rows, \
          B.cols, R alpha, A.data, A.stride, B.data, B.stride);                                                    \
    }

REPEAT(TRMM, strmm, dtrmm, ctrmm, ztrmm)

#define TRSM(F, T, R)                                                                                              \
    inline void solve(Triangle<T> A, General<T> B, T alpha)                                                        \
    {                                                                                                              \
        CHECK(A.major == B.major) CHECK(A.size == B.rows) CHECK(B.state == State::None);                           \
        F(Impl::blcvt(A.major), CblasLeft, Impl::blcvt(A.tri), Impl::blcvt(A.state), Impl::blcvt(A.diag), B.rows,  \
          B.cols, R alpha, A.data, A.stride, B.data, B.stride);                                                    \
    }                                                                                                              \
    inline void solve(General<T> B, Triangle<T> A, T alpha)                                                        \
    {                                                                                                              \
        CHECK(A.major == B.major) CHECK(A.size == B.cols) CHECK(B.state == State::None);                           \
        F(Impl::blcvt(A.major), CblasRight, Impl::blcvt(A.tri), Impl::blcvt(A.state), Impl::blcvt(A.diag), B.rows, \
          B.cols, R alpha, A.data, A.stride, B.data, B.stride);                                                    \
    }

REPEAT(TRSM, strsm, dtrsm, ctrsm, ztrsm)

#define SYRK(F, T, R)                                                                                            \
    inline void update(General<T> A, Symmetric<T> C, T alpha, T beta)                                            \
    {                                                                                                            \
        CHECK(A.major == C.major) CHECK(A._rows() == C.size);                                                    \
        F(Impl::blcvt(A.major), Impl::blcvt(C.tri), Impl::blcvt(A.state), A._rows(), A._cols(), R alpha, A.data, \
          A.stride, R beta, C.data, C.stride);                                                                   \
    }

REPEAT(SYRK, ssyrk, dsyrk, csyrk, zsyrk)

#define HERK(F, T, R)                                                                                          \
    inline void update(General<T> A, Hermitian<T> C, From<T>::Type alpha, From<T>::Type beta)                  \
    {                                                                                                          \
        CHECK(A.major == C.major) CHECK(A._rows() == C.size);                                                  \
        F(Impl::blcvt(A.major), Impl::blcvt(C.tri), Impl::blcvt(A.state), A._rows(), A._cols(), alpha, A.data, \
          A.stride, beta, C.data, C.stride);                                                                   \
    }

REPEAT1(HERK, cherk, zherk)

#define SYR2K(F, T, R)                                                                                           \
    inline void update(General<T> A, General<T> B, Symmetric<T> C, T alpha, T beta)                              \
    {                                                                                                            \
        CHECK(A.major == B.major) CHECK(A.major == C.major) CHECK(A.state == B.state);                           \
        CHECK(A._rows() == C.size) CHECK(A._rows() == B._rows()) CHECK(A._cols() == B._cols());                  \
        F(Impl::blcvt(A.major), Impl::blcvt(C.tri), Impl::blcvt(A.state), A._rows(), B._cols(), R alpha, A.data, \
          A.stride, B.data, B.stride, R beta, C.data, C.stride);                                                 \
    }

REPEAT(SYR2K, ssyr2k, dsyr2k, csyr2k, zsyr2k)

#define HER2K(F, T, R)                                                                                           \
    inline void update(General<T> A, General<T> B, Hermitian<T> C, T alpha, From<T>::Type beta)                  \
    {                                                                                                            \
        CHECK(A.major == B.major) CHECK(A.major == C.major) CHECK(A.state == B.state);                           \
        CHECK(A._rows() == C.size) CHECK(A._rows() == B._rows()) CHECK(A._cols() == B._cols());                  \
        F(Impl::blcvt(A.major), Impl::blcvt(C.tri), Impl::blcvt(A.state), A._rows(), B._cols(), R alpha, A.data, \
          A.stride, B.data, B.stride, beta, C.data, C.stride);                                                   \
    }

REPEAT1(HER2K, cher2k, zher2k)

////////////////////////////////////////////////////////////////////////////////
/// LAPACK /////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

#ifdef BLASW_LAPACKE_FOUND

#undef REPEAT0
#undef REPEAT1
#undef REPEAT

#define REPEAT0(func, name1, name2) \
    func(LAPACKE_##name1, float, ); \
    func(LAPACKE_##name2, double, );

#define REPEAT1(func, name1, name2)                        \
    func(LAPACKE_##name1, std::complex<float>, (void *)&); \
    func(LAPACKE_##name2, std::complex<double>, (void *)&);

#define REPEAT(func, name1, name2, name3, name4) \
    REPEAT0(func, name1, name2)                  \
    REPEAT1(func, name3, name4)

namespace Impl
{
inline int lpcvt(Major major) { return major == Major::Row ? LAPACK_ROW_MAJOR : LAPACK_ROW_MAJOR; }

inline char lpcvt(State state)
{
    if (state == State::None)
        return 'N';
    else if (state == State::Trans)
        return 'T';
    else if (state == State::ConjTrans)
        return 'C';
    else
        return '\0';
}

inline char lpcvt(Triangular tri) { return tri == Triangular::Upper ? 'U' : 'L'; }

template <typename T>
std::unique_ptr<T> alloc(Size size)
{
    auto P = aligned_alloc(64, size * sizeof(T));
    return std::unique_ptr<T>((T *)P);
}
}  // namespace Impl

#define GETRI(F, T, R)                                                                                    \
    inline bool inverse(General<T> A)                                                                     \
    {                                                                                                     \
        CHECK(A.rows == A.cols) CHECK(A.state == State::None);                                            \
        auto P = Impl::alloc<int>(A.rows);                                                                \
        if (F##getrf(Impl::lpcvt(A.major), A.rows, A.cols, A.data, A.stride, P.get()) != 0) return false; \
        return F##getri(Impl::lpcvt(A.major), A.rows, A.data, A.stride, P.get()) == 0;                    \
    }

REPEAT(GETRI, s, d, c, z)

#define SYTRI(F, T, R)                                                                                                \
    inline bool inverse(Symmetric<T> A)                                                                               \
    {                                                                                                                 \
        auto P = Impl::alloc<int>(A.size);                                                                            \
        if (F##sytrf(Impl::lpcvt(A.major), Impl::lpcvt(A.tri), A.size, A.data, A.stride, P.get()) != 0) return false; \
        return F##sytri(Impl::lpcvt(A.major), Impl::lpcvt(A.tri), A.size, A.data, A.stride, P.get()) == 0;            \
    }

REPEAT(SYTRI, s, d, c, z)

#define GEDTR(F, T, R)                                                                                 \
    inline T determinant(General<T> A)                                                                 \
    {                                                                                                  \
        CHECK(A.rows == A.cols) CHECK(A.state == State::None);                                         \
        auto _P = Impl::alloc<int>(A.rows);                                                            \
        if (F##getrf(Impl::lpcvt(A.major), A.rows, A.cols, A.data, A.stride, _P.get()) != 0) return 0; \
                                                                                                       \
        T det = 1;                                                                                     \
        auto P = _P.get();                                                                             \
        for (Size i = 0; i < A.rows; ++i) det *= A.data[i * A.stride + i] * T(P[i] == i ? 1 : -1);     \
        return det;                                                                                    \
    }

REPEAT(GEDTR, s, d, c, z)

#define SYDTR(F, T, R)                                                                                             \
    inline T determinant(Symmetric<T> A)                                                                           \
    {                                                                                                              \
        auto _P = Impl::alloc<int>(A.size);                                                                        \
        if (F##sytrf(Impl::lpcvt(A.major), Impl::lpcvt(A.tri), A.size, A.data, A.stride, _P.get()) != 0) return 0; \
                                                                                                                   \
        T det = 1;                                                                                                 \
        auto P = _P.get();                                                                                         \
        for (Size i = 0; i < A.size; ++i) det *= A.data[i * A.stride + i] * T(P[i] == i ? 1 : -1);                 \
        return det;                                                                                                \
    }

REPEAT(SYDTR, s, d, c, z)

#define HEDTR(F, T, R)                                                                                             \
    inline T determinant(Hermitian<T> A)                                                                           \
    {                                                                                                              \
        auto _P = Impl::alloc<int>(A.size);                                                                        \
        if (F##hetrf(Impl::lpcvt(A.major), Impl::lpcvt(A.tri), A.size, A.data, A.stride, _P.get()) != 0) return 0; \
                                                                                                                   \
        T det = 1;                                                                                                 \
        auto P = _P.get();                                                                                         \
        for (Size i = 0; i < A.size; ++i) det *= A.data[i * A.stride + i] * T(P[i] == i ? 1 : -1);                 \
        return det;                                                                                                \
    }

REPEAT1(HEDTR, c, z)

#define GETRF(F, T, R)                                                                                \
    inline bool lufact(General<T> A, Vector<int> P)                                                   \
    {                                                                                                 \
        CHECK(A.state == State::None) CHECK(std::min(A.rows, A.cols) == P.size) CHECK(P.stride == 1); \
        return F(Impl::lpcvt(A.major), A.rows, A.cols, A.data, A.stride, P.data) == 0;                \
    }

REPEAT(GETRF, sgetrf, dgetrf, cgetrf, zgetrf)

#define SYTRF(F, T, R)                                                                             \
    inline bool lufact(Symmetric<T> A, Vector<int> P)                                              \
    {                                                                                              \
        CHECK(A.size == P.size) CHECK(P.stride == 1);                                              \
        return F(Impl::lpcvt(A.major), Impl::lpcvt(A.tri), A.size, A.data, A.stride, P.data) == 0; \
    }

REPEAT(SYTRF, ssytrf, dsytrf, csytrf, zsytrf)

#define HETRF(F, T, R)                                                                             \
    inline bool lufact(Hermitian<T> A, Vector<int> P)                                              \
    {                                                                                              \
        CHECK(A.size == P.size) CHECK(P.stride == 1);                                              \
        return F(Impl::lpcvt(A.major), Impl::lpcvt(A.tri), A.size, A.data, A.stride, P.data) == 0; \
    }

REPEAT1(HETRF, chetrf, zhetrf)

#define POTRF(F, T, R)                                                                     \
    inline bool cholesky(Posdef<T> A)                                                      \
    {                                                                                      \
        return F(Impl::lpcvt(A.major), Impl::lpcvt(A.tri), A.size, A.data, A.stride) == 0; \
    }

REPEAT(POTRF, spotrf, dpotrf, cpotrf, zpotrf)

#define GEQRF(F, TYPE, R)                                                                             \
    inline bool qrfact(General<TYPE> A, Vector<TYPE> T)                                               \
    {                                                                                                 \
        CHECK(A.state == State::None) CHECK(std::min(A.rows, A.cols) == T.size) CHECK(T.stride == 1); \
        return F(Impl::lpcvt(A.major), A.rows, A.cols, A.data, A.stride, T.data) == 0;                \
    }

REPEAT(GEQRF, sgeqrf, dgeqrf, cgeqrf, zgeqrf)

#define GEEV(F, T, REF)                                                                                           \
    inline bool eigen(General<T> A, Vector<std::complex<T>> E, General<T> L, General<T> R)                        \
    {                                                                                                             \
        CHECK(A.major == L.major) CHECK(A.major == R.major);                                                      \
        CHECK(A.state == State::None) CHECK(L.state == State::None) CHECK(R.state == State::None);                \
        CHECK(A.rows == A.cols) CHECK(L.rows == L.cols) CHECK(R.rows == R.cols);                                  \
        CHECK(A.rows == E.size) CHECK(E.stride == 1);                                                             \
                                                                                                                  \
        if (L.data == nullptr) L.stride = A.rows;                                                                 \
        if (R.data == nullptr) R.stride = A.rows;                                                                 \
        auto _WR = Impl::alloc<T>(A.rows), _WI = Impl::alloc<T>(A.rows);                                          \
        if (F(Impl::lpcvt(A.major), L.data == nullptr ? 'N' : 'V', R.data == nullptr ? 'N' : 'V', A.rows, A.data, \
              A.stride, _WR.get(), _WI.get(), L.data, L.stride, R.data, R.stride) != 0)                           \
            return false;                                                                                         \
                                                                                                                  \
        auto WR = _WR.get(), WI = _WI.get();                                                                      \
        for (Size i = 0; i < A.rows; ++i) E.data[i].real(WR[i]), E.data[i].imag(WI[i]);                           \
        return true;                                                                                              \
    }

#define GEEVC(F, T, REF)                                                                                             \
    inline bool eigen(General<T> A, Vector<T> E, General<T> L, General<T> R)                                         \
    {                                                                                                                \
        CHECK(A.major == L.major) CHECK(A.major == R.major);                                                         \
        CHECK(A.state == State::None) CHECK(L.state == State::None) CHECK(R.state == State::None);                   \
        CHECK(A.rows == A.cols) CHECK(L.rows == L.cols) CHECK(R.rows == R.cols);                                     \
        CHECK(A.rows == E.size) CHECK(E.stride == 1);                                                                \
                                                                                                                     \
        if (L.data == nullptr) L.stride = A.rows;                                                                    \
        if (R.data == nullptr) R.stride = A.rows;                                                                    \
        return F(Impl::lpcvt(A.major), L.data == nullptr ? 'N' : 'V', R.data == nullptr ? 'N' : 'V', A.rows, A.data, \
                 A.stride, E.data, L.data, L.stride, R.data, R.stride) == 0;                                         \
    }

REPEAT0(GEEV, sgeev, dgeev)
REPEAT1(GEEVC, cgeev, zgeev)

#define SYEV(F, T, R)                                                                                 \
    inline bool eigen(Symmetric<T> A, Vector<T> E, bool vectors)                                      \
    {                                                                                                 \
        CHECK(A.size == E.size) CHECK(E.stride == 1);                                                 \
        auto V = vectors ? 'V' : 'N';                                                                 \
        return F(Impl::lpcvt(A.major), V, Impl::lpcvt(A.tri), A.size, A.data, A.stride, E.data) == 0; \
    }

#define HEEV(F, T, R)                                                                                 \
    inline bool eigen(Hermitian<T> A, Vector<From<T>::Type> E, bool vectors)                          \
    {                                                                                                 \
        CHECK(A.size == E.size) CHECK(E.stride == 1);                                                 \
        auto V = vectors ? 'V' : 'N';                                                                 \
        return F(Impl::lpcvt(A.major), V, Impl::lpcvt(A.tri), A.size, A.data, A.stride, E.data) == 0; \
    }

REPEAT0(SYEV, ssyev, dsyev)
REPEAT1(HEEV, cheev, zheev)

#define GEES(F, T, R)                                                                                             \
    inline bool schur(General<T> A, Vector<std::complex<T>> E, General<T> V)                                      \
    {                                                                                                             \
        CHECK(A.major == V.major) CHECK(A.state == State::None) CHECK(V.state == State::None);                    \
        CHECK(A.rows == A.cols) CHECK(A.rows == E.size) CHECK(E.stride == 1) CHECK(V.rows == V.cols);             \
                                                                                                                  \
        int sdim = 0;                                                                                             \
        if (V.data == nullptr) V.stride = A.rows;                                                                 \
        auto _WR = Impl::alloc<T>(A.rows), _WI = Impl::alloc<T>(A.rows);                                          \
        if (F(Impl::lpcvt(A.major), V.data == nullptr ? 'N' : 'V', 'N', nullptr, A.rows, A.data, A.stride, &sdim, \
              _WR.get(), _WI.get(), V.data, V.stride) != 0)                                                       \
            return false;                                                                                         \
                                                                                                                  \
        auto WR = _WR.get(), WI = _WI.get();                                                                      \
        for (Size i = 0; i < A.rows; ++i) E.data[i].real(WR[i]), E.data[i].imag(WI[i]);                           \
        return true;                                                                                              \
    }

#define GEESC(F, T, R)                                                                                               \
    inline bool schur(General<T> A, Vector<T> E, General<T> V)                                                       \
    {                                                                                                                \
        CHECK(A.major == V.major) CHECK(A.state == State::None) CHECK(V.state == State::None);                       \
        CHECK(A.rows == A.cols) CHECK(A.rows == E.size) CHECK(E.stride == 1) CHECK(V.rows == V.cols);                \
                                                                                                                     \
        int sdim = 0;                                                                                                \
        if (V.data == nullptr) V.stride = A.rows;                                                                    \
        return F(Impl::lpcvt(A.major), V.data == nullptr ? 'N' : 'V', 'N', nullptr, A.rows, A.data, A.stride, &sdim, \
                 E.data, V.data, V.stride) == 0;                                                                     \
    }

REPEAT0(GEES, sgees, dgees)
REPEAT1(GEESC, cgees, zgees)

#define GESV(F, T, R)                                                                                     \
    inline bool solve(General<T> A, General<T> B)                                                         \
    {                                                                                                     \
        CHECK(A.major == B.major) CHECK(A.state == State::None) CHECK(B.state == State::None);            \
        CHECK(A.rows == A.cols) CHECK(A.rows == B.rows);                                                  \
        auto P = Impl::alloc<int>(A.rows);                                                                \
        return F(Impl::lpcvt(A.major), A.rows, B.cols, A.data, A.stride, P.get(), B.data, B.stride) == 0; \
    }

REPEAT(GESV, sgesv, dgesv, cgesv, zgesv)

#define SYSV(F, T, R)                                                                                         \
    inline bool solve(Symmetric<T> A, General<T> B)                                                           \
    {                                                                                                         \
        CHECK(A.major == B.major) CHECK(B.state == State::None) CHECK(A.size == B.rows);                      \
        auto P = Impl::alloc<int>(A.size);                                                                    \
        return F(Impl::lpcvt(A.major), Impl::lpcvt(A.tri), A.size, B.cols, A.data, A.stride, P.get(), B.data, \
                 B.stride) == 0;                                                                              \
    }

REPEAT(SYSV, ssysv, dsysv, csysv, zsysv)

#define HESV(F, T, R)                                                                                         \
    inline bool solve(Hermitian<T> A, General<T> B)                                                           \
    {                                                                                                         \
        CHECK(A.major == B.major) CHECK(B.state == State::None) CHECK(A.size == B.rows);                      \
        auto P = Impl::alloc<int>(A.size);                                                                    \
        return F(Impl::lpcvt(A.major), Impl::lpcvt(A.tri), A.size, B.cols, A.data, A.stride, P.get(), B.data, \
                 B.stride) == 0;                                                                              \
    }

REPEAT1(HESV, chesv, zhesv)

#define GELS(F, T, R)                                                                                                \
    inline bool lsquares(General<T> A, General<T> B)                                                                 \
    {                                                                                                                \
        CHECK(A.major == B.major) CHECK(B.state == State::None) CHECK(A._rows() == B.rows);                          \
        return F(Impl::lpcvt(A.major), Impl::lpcvt(A.state), A._rows(), A._cols(), B.cols, A.data, A.stride, B.data, \
                 B.stride) == 0;                                                                                     \
    }

REPEAT(GELS, sgels, dgels, cgels, zgels)

#define GESVD(F, T, R)                                                                                                \
    inline bool svd(General<T> A, Vector<From<T>::Type> S, General<T> U, General<T> VT)                               \
    {                                                                                                                 \
        CHECK(A.major == U.major) CHECK(A.major == VT.major);                                                         \
        CHECK(A.state == State::None) CHECK(U.state == State::None) CHECK(VT.state == State::None);                   \
        CHECK(std::min(A.rows, A.cols) == S.size) CHECK(S.stride == 1);                                               \
                                                                                                                      \
        if (U.data != nullptr) CHECK(A.rows == U.rows) CHECK(U.rows == U.cols);                                       \
        if (VT.data != nullptr) CHECK(A.cols == VT.rows) CHECK(VT.rows == VT.cols);                                   \
                                                                                                                      \
        if (U.data == nullptr) U.stride = A.rows;                                                                     \
        if (VT.data == nullptr) VT.stride = A.cols;                                                                   \
        auto SB = Impl::alloc<From<T>::Type>(std::min(A.rows, A.cols));                                               \
        return F(Impl::lpcvt(A.major), U.data == nullptr ? 'N' : 'A', VT.data == nullptr ? 'N' : 'A', A.rows, A.cols, \
                 A.data, A.stride, S.data, U.data, U.stride, VT.data, VT.stride, SB.get()) == 0;                      \
    }

REPEAT(GESVD, sgesvd, dgesvd, cgesvd, zgesvd)

#define GERK(F, T, R)                                                                              \
    inline Size rank(General<T> A, From<T>::Type epsilon)                                          \
    {                                                                                              \
        CHECK(A.state == State::None);                                                             \
        auto SB = Impl::alloc<From<T>::Type>(std::min(A.rows, A.cols));                            \
        auto _S = Impl::alloc<From<T>::Type>(std::min(A.rows, A.cols));                            \
        if (F(Impl::lpcvt(A.major), 'N', 'N', A.rows, A.cols, A.data, A.stride, _S.get(), nullptr, \
              std::max(A.rows, A.cols), nullptr, std::max(A.rows, A.cols), SB.get()) != 0)         \
            return 0;                                                                              \
                                                                                                   \
        Size rank = 0;                                                                             \
        auto S = _S.get();                                                                         \
        for (Size i = 0; i < std::min(A.rows, A.cols); ++i)                                        \
            if (std::abs(S[i]) >= std::abs(epsilon)) ++rank;                                       \
        return rank;                                                                               \
    }

REPEAT(GERK, sgesvd, dgesvd, cgesvd, zgesvd)

#undef GETRI
#undef SYTRI
#undef GEDTR
#undef SYDTR
#undef HEDTR
#undef GETRF
#undef SYTRF
#undef HETRF
#undef POTRF
#undef GEQRF
#undef GEEV
#undef GEEVC
#undef SYEV
#undef HEEV
#undef GEES
#undef GEESC
#undef GESV
#undef SYSV
#undef HESV
#undef GELS
#undef GESVD
#undef GERK

#endif
}  // namespace Blasw

#undef CHECK
#undef REPEAT0
#undef REPEAT1
#undef REPEAT
#undef STATE
#undef GENERAL
#undef BGENERAL
#undef TRIANGLE
#undef BTRIANGLE
#undef PTRIANGLE
#undef SYMMETRIC
#undef BSYMMETRIC
#undef PSYMMETRIC
#undef HERMITIAN
#undef BHERMITIAN
#undef PHERMITIAN
#undef POSDEF
#undef ROTG
#undef ROTMG
#undef ROT
#undef ROTM
#undef SWAP
#undef SCAL
#undef SSCAL
#undef COPY
#undef AXPY
#undef DOT
#undef DOTU
#undef NRM2
#undef ASUM
#undef IAMAX
#undef GEMV
#undef GBMV
#undef SYMV
#undef HEMV
#undef SBMV
#undef HBMV
#undef SPMV
#undef HPMV
#undef TRMV
#undef TBMV
#undef TPMV
#undef TRSV
#undef TBSV
#undef TPSV
#undef GER
#undef GERU
#undef SYR
#undef HER
#undef SPR
#undef HPR
#undef SYR2
#undef HER2
#undef SPR2
#undef HPR2
#undef GEMM
#undef SYMM
#undef HEMM
#undef TRMM
#undef TRSM
#undef SYRK
#undef HERK
#undef SYR2K
#undef HER2K
