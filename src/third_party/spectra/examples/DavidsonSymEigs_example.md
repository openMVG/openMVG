This is an example of how to use the Jacobi-Davidson Symmetric Eigenvalue Solver with DPR correction method. This test can also be found as a full file in the [example/DavidsonSymEigs_example.cpp](example/DavidsonSymEigs_example.cpp) file and can be compiled with cmake and run afterwards

```bash
mkdir build && cd build && cmake ../
make DavidsonSymEigs_example
./example/DavidsonSymEigs_example
```

Suppose we want to find the 2 eigenpairs with the Largest value from a 1000x1000 Matrix A, then we could use this solver to quickly find them.


- First we have to construct the matrix

`Note: The Matrix has to be diagonally dominant otherwise the method will not converge`

```cpp

#include <Eigen/Dense>
#include <Spectra/DavidsonSymEigsSolver.h>
#include <Spectra/MatOp/DenseSymMatProd.h> 
#include <iostream>

using namespace Spectra;

int main()
{
Eigen::Index n = 1000;
    Eigen::MatrixXd mat = 0.03 * Eigen::MatrixXd::Random(n, n);
    Eigen::MatrixXd mat1 = mat + mat.transpose();
    for (Eigen::Index i=0; i<n; i++) {
        mat1(i,i) += i+1;
    }
```

- Then we have to construct a Matrix Product operation, which is provided by Spectra for Dense Symmetric Eigen matrices. 

`Note: For the solver only a Matrix product operation is required, thus you can specify a custom one without underlying matrix if you wish`

```cpp

DenseSymMatProd<double> op(mat); // Create the Matrix Product operation
```

- Afterwards the solver can be constructed. Both the operator and the desired number of eigen values must be specfied in the constructor. While their defaults values should be adequate for most situations, several internal parameters of the solver, most notably the maximum size of the search space and the number of correction vectors to append to the search size at each iteration, can be tuned :

```cpp
Eigen::Index num_of_eigenvalues = 5;
DavidsonSymEigsSolver<DenseSymMatProd<double>> solver(op_dense, num_of_eigenvalues);  //Create Solver
```

- This solver can then be executed through the compute method, where we also specify which EigenPairs we want through the [Sortrule enum](https://spectralib.org/doc/selectionrule_8h_source). The maximum number of iterations of the solver as well as the convergence criteria for the 
norm of the residues can also be specified in the call of the `compute()` method. `compute()` returns the number of converged eigenvalues.

```cpp

    Eigen::Index iterations = 100;
    double tolerance = 1e-3;
    int nconv = solver.compute(SortRule::LargestAlge, iterations, tolerance);

     // Retrieve results
    Eigen::VectorXd evalues;
    if (solver.info() == CompInfo::Successful){
        evalues = solver.eigenvalues();

    std::cout <<nconv<< " Eigenvalues found:\n"
              << evalues << std::endl;
    }else{
        std::cout <<"Calculation failed"<<std::endl;
    }
    return 0;
}
```

- It is also possible to provide a staring values for the eigenvectors. This can be done with the `compute_with_guess` method, which takes an additional `Eigen::Matrix` as an input. The guess is expected to be normalized.
