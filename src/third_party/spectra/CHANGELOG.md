## [1.0.1] - 2022-04-06
### Added
- Added SIMD support for `UpperHessenbergSchur`. This should accelerate
  general eigen solvers such as `GenEigsSolver`
- Added test code for `UpperHessenbergSchur`

### Changed
- Fixed several bugs in the examples caused by the `const` keyword,
  reported by [@alexpghayes](https://github.com/alexpghayes)
  ([#135](https://github.com/yixuan/spectra/issues/135),
  [#137](https://github.com/yixuan/spectra/pull/137))
- Updated the included [Catch2](https://github.com/catchorg/Catch2) to v2.13.8


## [1.0.0] - 2021-07-01
### Added
- Added version macros `SPECTRA_MAJOR_VERSION`, `SPECTRA_MINOR_VERSION`,
  `SPECTRA_PATCH_VERSION`, and `SPECTRA_VERSION` that are included by all eigen solvers
- Added the wrapper class `SparseGenComplexShiftSolve` for eigen solver with complex shifts
- Added the `SymGEigsShiftSolver` class for symmetric generalized eigen solver with real shifts
- Added the wrapper class `SymShiftInvert` that can be used with `SymGEigsShiftSolver`
- Added test code for symmetric generalized eigen solver with real shifts
- Added an internal class `UpperHessenbergSchur` to compute the Schur decomposition of
  upper Hessenberg matrices more efficiently
- Added a `Flags` template parameter to every matrix operation class
  (e.g. `DenseCholesky` and `DenseSymMatProd`), whose possible values are `Eigen::ColMajor`
  and `Eigen::RowMajor`. This parameter allows these wrapper classes to handle row-major matrices.
  If the input matrix is inconsistent with the `Flags` parameter (e.g., if `Flags` is
  `Eigen::ColMajor` but the input matrix is row-major), a compiler error will occur
- Added the member function `info()` and convergence tests to `SparseRegularInverse`,
  suggested by [@Spammed](https://github.com/Spammed) ([#111](https://github.com/yixuan/spectra/issues/111))
- Added symmetric Davidson eigen solver `DavidsonSymEigsSolver`, written by Felipe Zapata,
  Nicolas Renaud, Victor Azizi, Pablo Lopez-Tarifa, and Jens Wehner from the Netherlands eScience Center
- Extended matrix operations in `DenseGenMatProd`, `DenseSymMatProd`, `SparseGenMatProd`, and
  `SparseSymMatProd` to handle matrix-matrix products and coefficient-wise accessors

### Changed
- **API change**: Spectra now requires C++11
- **API change**: All enumerations have been converted to enum classes
  (e.g. `LARGEST_MAGN` is now `SortRule::LargestMagn`)
- **API change**: Selection rules are no longer template parameters. They are now
  specified in the `compute()` member function as arguments
- **API change**: The `Scalar` template parameter has been removed from eigen solvers.
  Instead, matrix operation classes now need to define a public type named `Scalar`
- **API change**: Constructors of solvers now request references of matrix operators
  instead of pointers
- Clang-Format now uses the C++11 standard to format code
- Updated documentation to reflect the new API
- Many internal changes to make use of C++11 features
- Added a `SPECTRA_` prefix to each header guard to prevent potential name clash
- Changed the default value of the `Flags` template parameter that exists in various
  class templates from `0` to the more readable constant `Eigen::ColMajor`
- Renamed the function `mat_prod` to `perform_op` in the `SparseRegularInverse` wrapper class.
  This makes the API more consistent when implementing new generalized eigen solvers
- Improved the precision of `UpperHessenbergQR` and `TridiagQR` by computing the
  Givens rotations in a more stable way
- Added a deflation test to `TridiagQR` to accelerate the convergence of eigen solvers
- Improved the precision of `TridiagQR::matrix_QtHQ()` by directly applying rotations
  to the original input matrix
- Improved the precision of `DoubleShiftQR` by computing the Householder reflectors
  in a more stable way
- Improved the deflation test in `DoubleShiftQR`
- More careful computation of residual vectors in the `Lanczos` process
- Initial vectors in the `Lanczos` and `Arnoldi` processes are now forced to be in the
  range of the `A` matrix
- More sensible test for orthogonality in generating new random vectors in the
  `Lanczos` and `Arnoldi` processes
- In symmetric eigen solvers large shifts are applied first to increase precision
- Updated the included [Catch2](https://github.com/catchorg/Catch2) to v2.13.6


## [0.9.0] - 2020-05-19
### Added
- Added support for CMake build, contributed by
  [Guillaume Acke](https://github.com/guacke) and [Jens Wehner](https://github.com/JensWehner)
  ([#70](https://github.com/yixuan/spectra/pull/70), [#88](https://github.com/yixuan/spectra/pull/88))
- Spectra can now be installed via [conda-forge](https://github.com/conda-forge/spectralib-feedstock),
  thanks to [Guillaume Acke](https://github.com/guacke) and [Julien Schueller](https://github.com/jschueller)
  ([#81](https://github.com/yixuan/spectra/pull/81), [#85](https://github.com/yixuan/spectra/pull/85))
- The source code of Spectra is now formatted using
  [Clang-Format](https://clang.llvm.org/docs/ClangFormat.html), suggested by
  [Jens Wehner](https://github.com/JensWehner)

### Changed
- Fixed a compiler warning caused by unused parameter, contributed by
  [Julien Schueller](https://github.com/jschueller) ([#80](https://github.com/yixuan/spectra/pull/80))
- Changed the implementation of `BKLDLT` solver to improve precision in some tests


## [0.8.1] - 2019-06-05
### Changed
- Fixed a bug in `BKLDLT` in which a wrong type was used, thanks to
  [@jdbancal](https://github.com/jdbancal) for the issue
  [#64](https://github.com/yixuan/spectra/pull/64)
- Fixed a bug in `BKLDLT` that caused segmentation fault in some edge
  cases, also reported by [@jdbancal](https://github.com/jdbancal) in issue
  [#66](https://github.com/yixuan/spectra/issues/66)
- The type `Eigen::Index` is now globally used for indices and sizes, in order to
  handle potentially large matrices. This was suggested by
  [Yuan Yao](https://github.com/y-yao) in issue
  [#19](https://github.com/yixuan/spectra/issues/19)


## [0.8.0] - 2019-04-03
### Added
- Added a `BKLDLT` class that implements the Bunch-Kaufman LDLT decomposition
  for symmetric indefinite matrices. According to the Eigen documentation,
  currently `Eigen::LDLT` cannot handle some special indefinite matrices such
  as `[0, 1; 1, 0]`, but `BKLDLT` is applicable to any symmetric matrices as
  long as it is not singular. LDLT decomposition is used in shift-and-invert
  solvers (see below)
- Added a unit test for `BKLDLT`

### Changed
- `DenseSymShiftSolve` now uses the newly added `BKLDLT` class to do the
  decomposition. This change broadens the class of matrices that
  `DenseSymShiftSolve` can handle, reduces memory use, and should also improve
  the numerical stability of the solver
- Replaced `Eigen::SimplicialLDLT` with `Eigen::SparseLU` in the `SparseSymShiftSolve`
  class, as some edge-case indefinite matrices may break `Eigen::SimplicialLDLT`
- `SparseSymShiftSolve` and `SparseGenRealShiftSolve` will throw an error if
  the factorization failed, for example, on singular matrices
- Fixed a missing `#include` in `DenseCholesky.h`, thanks to
  [Lennart Trunk](https://github.com/TheScarfix) for the issue
  [#59](https://github.com/yixuan/spectra/issues/59)
- Fixed errors in examples ([#60](https://github.com/yixuan/spectra/issues/60)),
  thanks to [@linuxfreebird](https://github.com/linuxfreebird)
- Updated the included [Catch2](https://github.com/catchorg/Catch2) to v2.7.0


## [0.7.0] - 2019-01-10
### Added
- Added a directory `contrib` to include code contributed by users. It is not
  formally a part of the Spectra library, but it may contain useful solvers
  and applications based on Spectra. Code in `contrib` may not be fully tested,
  so please use with caution. Feedback and report of issues are always welcome
- Added an eigen solver `LOBPCGSolver` in the `contrib` directory using the
  [LOBPCG](https://en.wikipedia.org/wiki/LOBPCG) algorithm,
  contributed by [Anna Araslanova](https://github.com/AnnaAraslanova)
- Added a partial SVD solver `PartialSVDSolver` in the `contrib` directory
- Added two internal classes `Arnoldi` and `Lanczos` to compute the Arnoldi/Lanczos
  factorization in eigen solvers
- Added a few other internal classes to refactor the eigen solver classes (see below)

### Changed
- **API change**: Spectra now requires Eigen >= 3.3
- **API change**: The library header files are moved into a directory
  named `Spectra`. Hence the recommended include directive would look like
  `#include <Spectra/SymEigsSolver.h>`
- All eigen solvers have been refactored using a cleaner class hierarchy.
  It may potentially make the implementation of new eigen solvers easier,
  especially for generalized eigen problems
- The matrix operation classes (e.g. `DenseSymMatProd` and `SparseSymMatProd`)
  are now internally using an
  [Eigen::Ref](https://eigen.tuxfamily.org/dox/classEigen_1_1Ref.html) object
  to wrap the user matrices, thanks to
  [Dario Mangoni](https://github.com/dariomangoni) who raised this issue in
  [#16](https://github.com/yixuan/spectra/issues/16)
- Fixed inappropriate range of random numbers in the tests
- Updated the included [Catch2](https://github.com/catchorg/Catch2) to v2.4.2


## [0.6.2] - 2018-05-22
### Changed
- Fixed regressions in v0.6.0 on some edge cases
- Improved the accuracy of restarting processes in `SymEigsSolver` and `GenEigsSolver`
- Updated the included [Catch2](https://github.com/catchorg/Catch2) to v2.2.2
- Code and documentation cleanup


## [0.6.1] - 2018-03-03
### Changed
- Fixed a bug of uninitialized memory
- Updated the included [Catch2](https://github.com/catchorg/Catch2) to v2.1.2


## [0.6.0] - 2018-03-03
### Added
- Added virtual destructors to the `SymEigsSolver` and `UpperHessenbergQR` classes
  to fix compiler warnings, by [Julian Kent](https://github.com/jkflying)
- Added a `NUMERICAL_ISSUE` entry to the `COMPUTATION_INFO` enumeration to indicate
  the status of Cholesky decomposition
- Added the `info()` member function to `DenseCholesky` and `SparseCholesky` to
  report the status of the decomposition
- Added a missing `#include` item in `SparseCholesky.h`, thanks to
  [Maxim Torgonsky](https://github.com/kriolog)
- Added a `TypeTraits` class to retrieve additional numeric limits of scalar value
  types

### Changed
- Documentation updates
- Updated the project URL to [https://spectralib.org](https://spectralib.org)
- Some internal improvements, such as pre-allocating vectors in loops, and changing
  return type to reference, thanks to
  [Angelos Mantzaflaris](https://github.com/filiatra)
- Improved the accuracy of symmetric and general eigen solvers
- Reduced the memory use of `UpperHessenbergQR` and `TridiagQR` decompositions
- Updated the included [Catch2](https://github.com/catchorg/Catch2) to v2.0.1
- Updated the testing code using the new API of Catch2
- Updated Travis CI script


## [0.5.0] - 2017-02-05
### Added
- Added the generalized eigen solver `SymGEigsSolver` in the regular inverse mode
- Added the wrapper class `SparseRegularInverse` that can be used with
  `SymGEigsSolver` in the regular inverse mode
- Added test code for generalized eigen solver in the regular inverse mode

### Changed
- Improved the numerical precision and stability of some internal linear
  algebra classes, including `TridiagEigen`, `UpperHessenbergEigen`, and
  `DoubleShiftQR`
- **API change**: The `x_in` argument in matrix operation functions, e.g.
  `perform_op()`, is now labelled to be constant
- Fixed a [bug](https://github.com/yixuan/spectra/issues/15) that
  `GenEigsComplexShiftSolver` gave wrong results when transforming back the
  eigenvalues, discovered by [@jdbancal](https://github.com/jdbancal)
- Updated included [Catch](https://github.com/philsquared/Catch) to v1.7.0
- Documentation improvement


## [0.4.0] - 2016-11-14
### Added
- Added an `Uplo` template parameter to the `DenseSymShiftSolve` class
- Added the generalized eigen solver `SymGEigsSolver` in the Cholesky mode
- Added the wrapper classes `DenseCholesky` and `SparseCholesky` that can be
  used with `SymGEigsSolver` in the Cholesky mode
- Added test code for generalized eigen solver in the Cholesky mode

### Changed
- Updated included [Catch](https://github.com/philsquared/Catch) to v1.5.7
- Improved documentation
- Updated Travis CI script
- Allowing basic math functions such as `abs()` and `sqrt()` to be overloaded
  (avoid using `std::abs` and `std::sqrt` directly), thanks to
  [@jdbancal](https://github.com/jdbancal). This makes it possible to use
  user-defined float number types with Spectra
- Replaced other `std` functions by their Eigen counterparts, for example using
  `Eigen::NumTraits<Scalar>::epsilon()` to substitute
  `std::numeric_limits<Scalar>::epsilon()`
- Improved the numerical stability of several operations, e.g. the function
  `hypot(x, y)` is used to compute `sqrt(x^2 + y^2)`
- More careful use of "approximate zero" constants
- Fixed an out-of-bound [bug](https://github.com/yixuan/spectra/issues/14)
  detected by [@jdbancal](https://github.com/jdbancal)


## [0.3.0] - 2016-07-03
### Added
- Added the wrapper classes `SparseSymMatProd` and `SparseSymShiftSolve`
  for sparse symmetric matrices
- Added the wrapper class `SparseGenRealShiftSolve` for general sparse matrices
- Added tests for sparse matrices
- Using Travis CI for automatic unit test

### Changed
- Updated included [Catch](https://github.com/philsquared/Catch) to v1.5.6
- **API change**: Each eigen solver was moved to its own header file.
  For example to use `SymEigsShiftSolver` one needs to include
  `<SymEigsShiftSolver.h>`
- Header files for internal use were relocated


## [0.2.0] - 2016-02-28
### Added
- Benchmark script now outputs number of matrix operations
- Added this change log
- Added a simple built-in random number generator, so that the algorithm
  was made to be deterministic
- Added the wrapper class `DenseSymMatProd` for symmetric matrices

### Changed
- Improved Arnoldi factorization
  - Iteratively corrects orthogonality
  - Creates new residual vector when invariant subspace is found
  - Stability for matrices with repeated eigenvalues is greatly improved
- Adjusted deflation tolerance in double shift QR
- Updated result analyzer
- Updated included [Catch](https://github.com/philsquared/Catch) to v1.3.4
- Updated copyright information
- **API change**: Default operator of `SymEigsSolver` was changed from
  `DenseGenMatProd` to `DenseSymMatProd`


## [0.1.0] - 2015-12-19
### Added
- Initial release of Spectra
