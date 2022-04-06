# <a href="https://spectralib.org"><img src="https://spectralib.org/img/logo.png" width="200px" /></a>

[![Build Status](https://travis-ci.org/yixuan/spectra.svg?branch=master)](https://travis-ci.org/yixuan/spectra) ![Basic CI](https://github.com/yixuan/spectra/workflows/Basic%20CI/badge.svg) [![codecov](https://codecov.io/gh/yixuan/spectra/branch/master/graph/badge.svg)](https://codecov.io/gh/yixuan/spectra)

> **NOTE**: Spectra 1.0.0 is released, with a lot of
> API-breaking changes. Please see the [migration guide](MIGRATION.md)
> for a smooth transition to the new version.

> **NOTE**: If you are interested in the future development of Spectra, please join
> [this thread](https://github.com/yixuan/spectra/issues/92) to share your comments and suggestions.

[**Spectra**](https://spectralib.org) stands for **Sp**arse **E**igenvalue **C**omputation **T**oolkit
as a **R**edesigned **A**RPACK. It is a C++ library for large scale eigenvalue
problems, built on top of [Eigen](http://eigen.tuxfamily.org),
an open source linear algebra library.

**Spectra** is implemented as a header-only C++ library, whose only dependency,
**Eigen**, is also header-only. Hence **Spectra** can be easily embedded in
C++ projects that require calculating eigenvalues of large matrices.

## Relation to ARPACK

[ARPACK](http://www.caam.rice.edu/software/ARPACK/) is a software package written in
FORTRAN for solving large scale eigenvalue problems. The development of
**Spectra** is much inspired by ARPACK, and as the full name indicates,
**Spectra** is a redesign of the ARPACK library using the C++ language.

In fact, **Spectra** is based on the algorithm described in the
[ARPACK Users' Guide](http://www.caam.rice.edu/software/ARPACK/UG/ug.html),
the implicitly restarted Arnoldi/Lanczos method. However,
it does not use the ARPACK code, and it is **NOT** a clone of ARPACK for C++.
In short, **Spectra** implements the major algorithms in ARPACK,
but **Spectra** provides a completely different interface, and it does not
depend on ARPACK.

## Common Usage

**Spectra** is designed to calculate a specified number (`k`) of eigenvalues
of a large square matrix (`A`). Usually `k` is much smaller than the size of the matrix
(`n`), so that only a few eigenvalues and eigenvectors are computed, which
in general is more efficient than calculating the whole spectral decomposition.
Users can choose eigenvalue selection rules to pick the eigenvalues of interest,
such as the largest `k` eigenvalues, or eigenvalues with largest real parts, etc.

To use the eigen solvers in this library, the user does not need to directly
provide the whole matrix, but instead, the algorithm only requires certain operations
defined on `A`. In the basic setting, it is simply the matrix-vector
multiplication. Therefore, if the matrix-vector product `A * x` can be computed
efficiently, which is the case when `A` is sparse, **Spectra**
will be very powerful for large scale eigenvalue problems.

There are two major steps to use the **Spectra** library:

1. Define a class that implements a certain matrix operation, for example the
matrix-vector multiplication `y = A * x` or the shift-solve operation
`y = inv(A - σ * I) * x`. **Spectra** has defined a number of
helper classes to quickly create such operations from a matrix object.
See the documentation of
[DenseGenMatProd](https://spectralib.org/doc/classSpectra_1_1DenseGenMatProd.html),
[DenseSymShiftSolve](https://spectralib.org/doc/classSpectra_1_1DenseSymShiftSolve.html), etc.
2. Create an object of one of the eigen solver classes, for example
[SymEigsSolver](https://spectralib.org/doc/classSpectra_1_1SymEigsSolver.html)
for symmetric matrices, and
[GenEigsSolver](https://spectralib.org/doc/classSpectra_1_1GenEigsSolver.html)
for general matrices. Member functions
of this object can then be called to conduct the computation and retrieve the
eigenvalues and/or eigenvectors.

Below is a list of the available eigen solvers in **Spectra**:

- [SymEigsSolver](https://spectralib.org/doc/classSpectra_1_1SymEigsSolver.html):
For real symmetric matrices
- [GenEigsSolver](https://spectralib.org/doc/classSpectra_1_1GenEigsSolver.html):
For general real matrices
- [SymEigsShiftSolver](https://spectralib.org/doc/classSpectra_1_1SymEigsShiftSolver.html):
For real symmetric matrices using the shift-and-invert mode
- [GenEigsRealShiftSolver](https://spectralib.org/doc/classSpectra_1_1GenEigsRealShiftSolver.html):
For general real matrices using the shift-and-invert mode,
with a real-valued shift
- [GenEigsComplexShiftSolver](https://spectralib.org/doc/classSpectra_1_1GenEigsComplexShiftSolver.html):
For general real matrices using the shift-and-invert mode,
with a complex-valued shift
- [SymGEigsSolver](https://spectralib.org/doc/classSpectra_1_1SymGEigsSolver.html):
For generalized eigen solver with real symmetric matrices
- [SymGEigsShiftSolver](https://spectralib.org/doc/classSpectra_1_1SymGEigsShiftSolver.html):
For generalized eigen solver with real symmetric matrices, using the shift-and-invert mode
- [DavidsonSymEigsSolver](https://spectralib.org/doc/classSpectra_1_1DavidsonSymEigsSolver.html):
Jacobi-Davidson eigen solver for real symmetric matrices, with the DPR correction method

## Examples

Below is an example that demonstrates the use of the eigen solver for symmetric
matrices.

```cpp
#include <Eigen/Core>
#include <Spectra/SymEigsSolver.h>
// <Spectra/MatOp/DenseSymMatProd.h> is implicitly included
#include <iostream>

using namespace Spectra;

int main()
{
    // We are going to calculate the eigenvalues of M
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(10, 10);
    Eigen::MatrixXd M = A + A.transpose();

    // Construct matrix operation object using the wrapper class DenseSymMatProd
    DenseSymMatProd<double> op(M);

    // Construct eigen solver object, requesting the largest three eigenvalues
    SymEigsSolver<DenseSymMatProd<double>> eigs(op, 3, 6);

    // Initialize and compute
    eigs.init();
    int nconv = eigs.compute(SortRule::LargestAlge);

    // Retrieve results
    Eigen::VectorXd evalues;
    if(eigs.info() == CompInfo::Successful)
        evalues = eigs.eigenvalues();

    std::cout << "Eigenvalues found:\n" << evalues << std::endl;

    return 0;
}
```

Sparse matrix is supported via classes such as `SparseGenMatProd`
and `SparseSymMatProd`.

```cpp
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <iostream>

using namespace Spectra;

int main()
{
    // A band matrix with 1 on the main diagonal, 2 on the below-main subdiagonal,
    // and 3 on the above-main subdiagonal
    const int n = 10;
    Eigen::SparseMatrix<double> M(n, n);
    M.reserve(Eigen::VectorXi::Constant(n, 3));
    for(int i = 0; i < n; i++)
    {
        M.insert(i, i) = 1.0;
        if(i > 0)
            M.insert(i - 1, i) = 3.0;
        if(i < n - 1)
            M.insert(i + 1, i) = 2.0;
    }

    // Construct matrix operation object using the wrapper class SparseGenMatProd
    SparseGenMatProd<double> op(M);

    // Construct eigen solver object, requesting the largest three eigenvalues
    GenEigsSolver<SparseGenMatProd<double>> eigs(op, 3, 6);

    // Initialize and compute
    eigs.init();
    int nconv = eigs.compute(SortRule::LargestMagn);

    // Retrieve results
    Eigen::VectorXcd evalues;
    if(eigs.info() == CompInfo::Successful)
        evalues = eigs.eigenvalues();

    std::cout << "Eigenvalues found:\n" << evalues << std::endl;

    return 0;
}
```

And here is an example for user-supplied matrix operation class.

```cpp
#include <Eigen/Core>
#include <Spectra/SymEigsSolver.h>
#include <iostream>

using namespace Spectra;

// M = diag(1, 2, ..., 10)
class MyDiagonalTen
{
public:
    using Scalar = double;  // A typedef named "Scalar" is required
    int rows() const { return 10; }
    int cols() const { return 10; }
    // y_out = M * x_in
    void perform_op(const double *x_in, double *y_out) const
    {
        for(int i = 0; i < rows(); i++)
        {
            y_out[i] = x_in[i] * (i + 1);
        }
    }
};

int main()
{
    MyDiagonalTen op;
    SymEigsSolver<MyDiagonalTen> eigs(op, 3, 6);
    eigs.init();
    eigs.compute(SortRule::LargestAlge);
    if(eigs.info() == CompInfo::Successful)
    {
        Eigen::VectorXd evalues = eigs.eigenvalues();
        std::cout << "Eigenvalues found:\n" << evalues << std::endl;
    }

    return 0;
}
```

## Shift-and-invert Mode

When it is needed to find eigenvalues that are closest to a number `σ`,
for example to find the smallest eigenvalues of a positive definite matrix
(in which case `σ = 0`), it is advised to use the shift-and-invert mode
of eigen solvers.

In the shift-and-invert mode, selection rules are applied to `1/(λ - σ)`
rather than `λ`, where `λ` are eigenvalues of `A`.
To use this mode, users need to define the shift-solve matrix operation. See
the documentation of
[SymEigsShiftSolver](https://spectralib.org/doc/classSpectra_1_1SymEigsShiftSolver.html)
for details.

## Documentation

The [API reference](https://spectralib.org/doc/) page contains the documentation
of **Spectra** generated by [Doxygen](http://www.doxygen.org/),
including all the background knowledge, example code and class APIs.

More information can be found in the project page [https://spectralib.org](https://spectralib.org).

## Installation

An optional CMake installation is supported, if you have CMake with at least v3.10 installed. You can install the headers using the following commands:

```bash
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX='intended installation directory' -DCMAKE_PREFIX_PATH='path where the installation of Eigen3 can be found' -DBUILD_TESTS=TRUE
make all && make tests && make install
```

By installing **Spectra** in this way, you also create a CMake target `Spectra::Spectra` that can be used in subsequent build procedures for other programs.

## License

**Spectra** is an open source project licensed under
[MPL2](https://www.mozilla.org/MPL/2.0/), the same license used by **Eigen**.
