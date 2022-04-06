## [0.9.0] → [1.0.0]

Spectra 1.0.0 introduces a lot of API-breaking changes, but the migration
should be straightforward following the guide below.

### Toolchain

Spectra 1.0.0 requires a compiler supporting the C++11 standard.
Any modern C++ compiler should already have this.

### Matrix Operation Classes

In most cases you do not need to change anything for the code involving
built-in matrix operation classes such as `DenseSymMatProd` and
`SparseGenMatProd`. However, if you have defined your own class, you
need to add a public type definition named `Scalar`, as the example
below shows. The type `Scalar` indicates the element type of the matrix.

```cpp
// A user-defined matrix operation class
// representing the matrix A=diag(1, 2, ..., 10)
class MyDiagonalTen
{
public:
    // The line below is new
    using Scalar = double;

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
```

### Eigen Solvers

The biggest change happens in the eigen solvers:

1. The first template parameter `Scalar` has been removed.
2. The second template parameter, `SelectionRule`, has been changed to
   a runtime parameter `selection` in the `compute()` member function.
3. In the constructor, matrix operation objects are now passed as
   references instead of pointers.
4. All enumerations have been converted to enum classes (see the
   conversion table below).

Below shows the one-to-one conversion of the code reflecting the
changes above:

[0.9.0] code:
```cpp
// Construct matrix operation object using the wrapper class
DenseSymMatProd<double> op(M);

// Construct eigen solver object, requesting the largest three eigenvalues
SymEigsSolver< double, LARGEST_ALGE, DenseSymMatProd<double> > eigs(&op, 3, 6);

// Initialize, and compute with at most 1000 iterations
eigs.init();
int nconv = eigs.compute(1000);

// Retrieve results
Eigen::VectorXd evalues;
if(eigs.info() == SUCCESSFUL)
    evalues = eigs.eigenvalues();
```

[1.0.0] code:
```cpp
// Construct matrix operation object using the wrapper class
DenseSymMatProd<double> op(M);

// Construct eigen solver object, requesting the largest three eigenvalues
SymEigsSolver<DenseSymMatProd<double>> eigs(op, 3, 6);

// Initialize, and compute with at most 1000 iterations
eigs.init();
int nconv = eigs.compute(SortRule::LargestAlge, 1000);

// Retrieve results
Eigen::VectorXd evalues;
if(eigs.info() == CompInfo::Successful)
    evalues = eigs.eigenvalues();
```

### Enumerations

| [0.9.0]                 | [1.0.0]                     |
|-------------------------|-----------------------------|
| `SUCCESSFUL`            | `CompInfo::Successful`      |
| `NOT_COMPUTED`          | `CompInfo::NotComputed`     |
| `NOT_CONVERGING`        | `CompInfo::NotConverging`   |
| `NUMERICAL_ISSUE`       | `CompInfo::NumericalIssue`  |
| `LARGEST_MAGN`          | `SortRule::LargestMagn`     |
| `LARGEST_REAL`          | `SortRule::LargestReal`     |
| `LARGEST_IMAG`          | `SortRule::LargestImag`     |
| `LARGEST_ALGE`          | `SortRule::LargestAlge`     |
| `SMALLEST_MAGN`         | `SortRule::SmallestMagn`    |
| `SMALLEST_REAL`         | `SortRule::SmallestReal`    |
| `SMALLEST_IMAG`         | `SortRule::SmallestImag`    |
| `SMALLEST_ALGE`         | `SortRule::SmallestAlge`    |
| `BOTH_ENDS`             | `SortRule::BothEnds`        |
| `GEIGS_CHOLESKY`        | `GEigsMode::Cholesky`       |
| `GEIGS_REGULAR_INVERSE` | `GEigsMode::RegularInverse` |
| `GEIGS_SHIFT_INVERT`    | `GEigsMode::ShiftInvert`    |
| `GEIGS_BUCKLING`        | `GEigsMode::Buckling`       |
| `GEIGS_CAYLEY`          | `GEigsMode::Cayley`         |
