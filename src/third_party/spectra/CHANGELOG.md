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
