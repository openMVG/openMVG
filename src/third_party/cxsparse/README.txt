CXSparse: Copyright (c) 2006-2012, Timothy A. Davis.
http://www.suitesparse.com

Derived from CSparse.  Conversion originally by David Bateman, Motorola,
and then modified by Tim Davis.  ANSI C99 is required, with support for
the _Complex data type.
(if you use a C++ compiler, the C++ complex type is used instead).

CXSparse is a version of CSparse that operates on both real and complex
matrices, using either int or SuiteSparse_long integers.  A SuiteSparse_long is
normally just a long on most platforms, but becomes __int64 on WIN64.  It now
includes a MATLAB interface, enabling the use of CXSparse functions on both
32-bit and 64-bit platforms.

To install for use in MATLAB, simply type "cs_install" in the MATLAB Command
Window, while in the CXSparse/MATLAB directory.  (NOTE: Windows users cannot
use the "lcc" command; run "mex -setup" first, and select a different
compiler).  If you use the Unix "make" command in that directory instead and
are using a 64-bit platform, then you must edit the CXSparse/MATLAB/Makefile
first.  Refer to the instructions in that file.

Refer to "Direct Methods for Sparse Linear Systems," Timothy A. Davis,
SIAM, Philadelphia, 2006.  No detailed user guide is included in this
package; the user guide is the book itself.

To compile the C-library (./Source), C demo programs (./Demo) just type "make"
in this directory.  To run the exhaustive statement coverage tests, type
"make" in the Tcov directory; the Tcov tests assume you are using Linux.  To
remove all files not in the original distribution, type "make distclean".
I recommend that you use a different level of
optimization than "cc -O", which was chosen so that the Makefile is portable.
See Source/Makefile.

If your C compiler does not support the ANSI C99 complex type, the
#include <complex.h> statement will fail.  If this happens, compile the
code with the -DNCOMPLEX flag (the MATLAB cs_install will do this for you).

This package is backward compatible with CSparse.  That is, user code that
uses CSparse may switch to using CXSparse without any changes to the user code.
Each CXSparse function has a generic version with the same name as the CSparse
function, and four type-specific versions.  For example:

    cs_add      same as cs_add_di by default, but can be changed to use
                SuiteSparse_long
                integers if user code is compiled with -DCS_LONG, and/or can
                be changed to operate on complex matrices with -DCS_COMPLEX.

    cs_di_add   double/int version of cs_add
    cs_dl_add   double/SuiteSparse_long version of cs_add
    cs_ci_add   complex/int version of cs_add
    cs_cl_add   complex/SuiteSparse_long version of cs_add

The sparse matrix data structures are treated in the same way:  cs, css,
csn, and csd become cs_di, cs_dis, cs_din, and cs_did for the double/int case,
cs_cl, cs_cls, cs_cln, and cs_cld for the complex/SuiteSparse_long case, and so
on.

See cs_demo.c for a type-generic user program, and cs_cl_demo.c for a
type-specific version of the same program (complex/SuiteSparse_long).

Several macros are available in CXSparse (but not in CSparse) to allow user
code to be written in a type-generic manner:

    CS_INT      int by default, SuiteSparse_long if -DCS_LONG compiler flag
                is used
    CS_ENTRY    double by default, double complex if -DCS_COMPLEX flag is used.
    CS_ID       "%d" or "%ld", for printf and scanf of the CS_INT type.
    CS_INT_MAX  INT_MAX or LONG_MAX, the largest possible value of CS_INT.
    CS_REAL(x)  x or creal(x)
    CS_IMAG(x)  0 or cimag(x)
    CS_CONJ(x)  x or conj(x)
    CS_ABS(x)   fabs(x) or cabs(x)

Even the name of the include file (cs.h) is the same.  To use CXSparse instead
of CSparse, simply compile with -ICXSparse/Source instead of -ICSparse/Source,
and link against libcxsparse.a instead of the CSparse libcsparse.a library.

To determine at compile time if CXSparse or CSparse is being used:

    #ifdef CXSPARSE
        CXSparse is in use.  The generic functions equivalent to CSparse may
        be used (cs_add, etc).  These generic functions can use different
        types, depending on the -DCS_LONG and -DCS_COMPLEX compile flags, with
        the default being double/int.  The type-specific functions and data
        types (cs_di_add, cs_di, CS_INT, etc.) can be used.
    #else
        CSparse is in use.  Only the generic functions "cs_add", etc., are
        available, and they are of type double/int.
    #endif

See cs.h for the prototypes of each function, and the book "Direct Methods
for Sparse Linear Systems" for full documentation of CSparse and CXSparse.

Other changes from CSparse:  cs_transpose performs the complex conjugate
transpose if values>0 (C=A'), the pattern-only transpose if values=0
(C=spones(A') in MATLAB), and the array transpose if values<0 (C=A.' in
MATLAB notation).  A set of four conversion routines are included in CXSparse,
to convert real matrices to/from complex matrices.
The Householder reflection constructed by cs_house.c also differs slightly, to
accomodate both the real and complex cases properly.

CXSparse is generated automatically from CSparse.  Refer to
http://www.suitesparse.com for details.

--------------------------------------------------------------------------------
Contents:
--------------------------------------------------------------------------------

Demo/           demo C programs that use CXSparse
Doc/            license and change log
Makefile        Makefile for the whole package
MATLAB/         MATLAB interface, demos, and tests for CXSparse
Matrix/         sample matrices (with extra complex matrices for CXSparse)
README.txt      this file
Source/         primary CXSparse source files
Tcov/           CXSparse tests

--------------------------------------------------------------------------------
./Doc:          license and change log
--------------------------------------------------------------------------------

ChangeLog       changes in CSparse since first release
lesser.txt      the GNU LGPL
License.txt     license (GNU LGPL)

--------------------------------------------------------------------------------
./Source:       Primary source code for CXSparse
--------------------------------------------------------------------------------

cs_add.c        add sparse matrices
cs_amd.c        approximate minimum degree
cs_chol.c       sparse Cholesky
cs_cholsol.c    x=A\b using sparse Cholesky
cs_compress.c   convert a compress form to compressed-column form
cs_counts.c     column counts for Cholesky and QR
cs_convert.c    convert real to complex and complex to real (not in CSparse)
cs_cumsum.c     cumulative sum
cs_dfs.c        depth-first-search
cs_dmperm.c     Dulmage-Mendelsohn permutation
cs_droptol.c    drop small entries from a sparse matrix
cs_dropzeros.c  drop zeros from a sparse matrix
cs_dupl.c       remove (and sum) duplicates
cs_entry.c      add an entry to a triplet matrix
cs_ereach.c     nonzero pattern of Cholesky L(k,:) from etree and triu(A(:,k))
cs_etree.c      find elimination tree
cs_fkeep.c      drop entries from a sparse matrix
cs_gaxpy.c      sparse matrix times dense matrix
cs.h            include file for CXSparse
cs_happly.c     apply Householder reflection
cs_house.c      Householder reflection (*** NOTE: different algo. from CSparse)
cs_ipvec.c      x(p)=b
cs_leaf.c       determine if j is a leaf of the skeleton matrix and find lca
cs_load.c       load a sparse matrix from a file
cs_lsolve.c     x=L\b
cs_ltsolve.c    x=L'\b
cs_lu.c         sparse LU factorization
cs_lusol.c      x=A\b using sparse LU factorization
cs_malloc.c     memory manager
cs_maxtrans.c   maximum transveral (permutation for zero-free diagonal)
cs_multiply.c   sparse matrix multiply
cs_norm.c       sparse matrix norm
cs_permute.c    permute a sparse matrix
cs_pinv.c       invert a permutation vector
cs_post.c       postorder an elimination tree
cs_print.c      print a sparse matrix
cs_pvec.c       x=b(p)
cs_qr.c         sparse QR
cs_qrsol.c      solve a least-squares problem
cs_randperm.c   random permutation
cs_reach.c      find nonzero pattern of x=L\b for sparse L and b
cs_scatter.c    scatter a sparse vector
cs_scc.c        strongly-connected components
cs_schol.c      symbolic Cholesky
cs_spsolve.c    x=Z\b where Z, x, and b are sparse, and Z upper/lower triangular
cs_sqr.c        symbolic QR (also can be used for LU)
cs_symperm.c    symmetric permutation of a sparse matrix
cs_tdfs.c       depth-first-search of a tree
cs_transpose.c  transpose a sparse matrix
cs_updown.c     sparse rank-1 Cholesky update/downate
cs_usolve.c     x=U\b
cs_util.c       various utilities (allocate/free matrices, workspace, etc)
cs_utsolve.c    x=U'\b
Makefile        Makefile for CXSparse
README.txt      README file for CXSparse

--------------------------------------------------------------------------------
./Demo:         C program demos
--------------------------------------------------------------------------------

cs_ci_demo1.c   complex/int version of cs_demo1.c
cs_ci_demo2.c   complex/int version of cs_demo2.c
cs_ci_demo3.c   complex/int version of cs_demo3.c
cs_ci_demo.c    complex/int version of cs_demo.c
cs_ci_demo.h    complex/int version of cs_demo.h

cs_cl_demo1.c   complex/SuiteSparse_long version of cs_demo1.c
cs_cl_demo2.c   complex/SuiteSparse_long version of cs_demo2.c
cs_cl_demo3.c   complex/SuiteSparse_long version of cs_demo3.c
cs_cl_demo.c    complex/SuiteSparse_long version of cs_demo.c
cs_cl_demo.h    complex/SuiteSparse_long version of cs_demo.h

cs_demo1.c      read a matrix from a file and perform basic matrix operations
cs_demo2.c      read a matrix from a file and solve a linear system
cs_demo3.c      read a matrix, solve a linear system, update/downdate
cs_demo.c       support routines for cs_demo*.c
cs_demo.h       include file for demo programs

cs_demo.out     output of "make", which runs the demos on some matrices

cs_di_demo1.c   double/int version of cs_demo1.c
cs_di_demo2.c   double/int version of cs_demo2.c
cs_di_demo3.c   double/int version of cs_demo3.c
cs_di_demo.c    double/int version of cs_demo.c
cs_di_demo.h    double/int version of cs_demo.h

cs_dl_demo1.c   double/SuiteSparse_long version of cs_demo1.c
cs_dl_demo2.c   double/SuiteSparse_long version of cs_demo2.c
cs_dl_demo3.c   double/SuiteSparse_long version of cs_demo3.c
cs_dl_demo.c    double/SuiteSparse_long version of cs_demo.c
cs_dl_demo.h    double/SuiteSparse_long version of cs_demo.h

cs_idemo.c      convert real matrices to/from complex (int version)
cs_ldemo.c      convert real matrices to/from complex (SuiteSparse_long version)

Makefile        Makefile for Demo programs
readhb.f        read a Rutherford-Boeing matrix (real matrices only)
README.txt      Demo README file

--------------------------------------------------------------------------------
./MATLAB:       MATLAB interface, demos, and tests
--------------------------------------------------------------------------------

cs_install.m    MATLAB function for compiling and installing CSparse for MATLAB
CSparse/        MATLAB interface for CSparse
Demo/           MATLAB demos for CSparse
Makefile        MATLAB interface Makefile
README.txt      MATLAB README file
Test/           MATLAB test for CSparse, and "textbook" routines
UFget/          MATLAB interface to UF Sparse Matrix Collection


--------------------------------------------------------------------------------
./MATLAB/CSparse:   MATLAB interface for CSparse
--------------------------------------------------------------------------------

Contents.m          Contents of MATLAB interface to CSparse
cs_add.m            add two sparse matrices
cs_add_mex.c
cs_amd.m            approximate minimum degree
cs_amd_mex.c
cs_chol.m           sparse Cholesky
cs_chol_mex.c
cs_cholsol.m        x=A\b using a sparse Cholesky
cs_cholsol_mex.c
cs_counts.m         column counts for Cholesky or QR (like "symbfact" in MATLAB)
cs_counts_mex.c
cs_dmperm.m         Dulmage-Mendelsohn permutation
cs_dmperm_mex.c
cs_dmsol.m          x=A\b using dmperm
cs_dmspy.m          plot a picture of a dmperm-permuted matrix
cs_droptol.m        drop small entries
cs_droptol_mex.c
cs_esep.m           find edge separator
cs_etree.m          compute elimination tree
cs_etree_mex.c
cs_gaxpy.m          sparse matrix times dense vector
cs_gaxpy_mex.c
cs_lsolve.m         x=L\b where L is lower triangular
cs_lsolve_mex.c
cs_ltsolve.m        x=L'\b where L is lower triangular
cs_ltsolve_mex.c
cs_lu.m             sparse LU factorization
cs_lu_mex.c
cs_lusol.m          x=A\b using sparse LU factorization
cs_lusol_mex.c
cs_make.m           compiles CSparse for use in MATLAB
cs_mex.c            support routines for CSparse mexFunctions
cs_mex.h
cs_multiply.m       sparse matrix multiply
cs_multiply_mex.c
cs_must_compile.m   determine if a source file needs to be compiled with mex
cs_nd.m             nested dissection
cs_nsep.m           find node separator
cs_permute.m        permute a sparse matrix
cs_permute_mex.c
cs_print.m          print a sparse matrix
cs_print_mex.c
cs_qleft.m          apply Householder vectors to the left
cs_qright.m         apply Householder vectors to the right
cs_qr.m             sparse QR factorization
cs_qr_mex.c
cs_qrsol.m          solve a sparse least squares problem
cs_qrsol_mex.c
cs_randperm.m       randdom permutation
cs_randperm_mex.c
cs_scc.m            strongly-connected components
cs_scc_mex.c
cs_sep.m            convert an edge separator into a node separator
cs_sparse.m         convert a triplet form matrix to a compress-column form
cs_sparse_mex.c
cs_symperm.m        symmetric permutation of a sparse matrix
cs_symperm_mex.c
cs_sqr.m            symbolic QR ordering and analysis
cs_sqr_mex.c
cs_thumb_mex.c      compute small "thumbnail" of a sparse matrix (for cspy).
cs_transpose.m      transpose a sparse matrix
cs_transpose_mex.c
cs_updown.m         sparse Cholesky update/downdate
cs_updown_mex.c
cs_usolve.m         x=U\b where U is upper triangular 
cs_usolve_mex.c
cs_utsolve.m        x=U'\b where U is upper triangular 
cs_utsolve_mex.c
cspy.m              a color "spy"
Makefile            Makefile for CSparse MATLAB interface
README.txt          README file for CSparse MATLAB interface


--------------------------------------------------------------------------------
./MATLAB/Demo:      MATLAB demos for CSparse
--------------------------------------------------------------------------------

Contents.m          Contents of MATLAB demo for CSparse
cs_demo.m           run all MATLAB demos for CSparse
cs_demo1.m          MATLAB version of Demo/cs_demo1.c
cs_demo2.m          MATLAB version of Demo/cs_demo2.c
cs_demo3.m          MATLAB version of Demo/cs_demo3.c
private/            private functions for MATLAB demos
README.txt          README file for CSparse MATLAB demo


--------------------------------------------------------------------------------
./MATLAB/Demo/private: private functions for MATLAB demos
--------------------------------------------------------------------------------

demo2.m             demo 2
demo3.m             demo 3
ex_1.m              example 1
ex2.m               example 2
ex3.m               example 3
frand.m             generate a random finite-element matrix
get_problem.m       get a matrix
is_sym.m            determine if a matrix is symmetric
mesh2d1.m           construct a 2D mesh (method 1)
mesh2d2.m           construct a 2D mesh (method 2)
mesh3d1.m           construct a 3D mesh (method 1)
mesh3d2.m           construct a 3D mesh (method 2)
print_order.m       print the ordering method used
resid.m             compute residual
rhs.m               create right-hand-side


--------------------------------------------------------------------------------
./MATLAB/Test:      Extensive test of CSparse, in MATLAB
--------------------------------------------------------------------------------

Makefile            Makefile for MATLAB Test directory
README.txt          README file for MATLAB/Test
Contents.m          Contents of MATLAB/Test, "textbook" files only

chol_downdate.m     downdate a Cholesky factorization.
chol_left.m         left-looking Cholesky factorization.
chol_left2.m        left-looking Cholesky factorization, more details.
chol_right.m        right-looking Cholesky factorization.
chol_super.m        left-looking "supernodal" Cholesky factorization.
chol_up.m           up-looking Cholesky factorization.
chol_update.m       update a Cholesky factorization.
chol_updown.m       update or downdate a Cholesky factorization.
cond1est.m          1-norm condition estimate.
cs_fiedler.m        the Fiedler vector of a connected graph.
givens2.m           find a Givens rotation.
house.m             find a Householder reflection.
lu_left.m           left-looking LU factorization.
lu_right.m          right-looking LU factorization.
lu_rightp.m         right-looking LU factorization, with partial pivoting.
lu_rightpr.m        recursive right-looking LU, with partial pivoting.
lu_rightr.m         recursive right-looking LU.
norm1est.m          1-norm estimate.
qr_givens.m         Givens-rotation QR factorization.
qr_givens_full.m    Givens-rotation QR factorization, for full matrices.
qr_left.m           left-looking Householder QR factorization.
qr_right.m          right-looking Householder QR factorization.
cs_fiedler.m        Fiedler vector

cs_frand.m          generate a random finite-element matrix
cs_frand_mex.c
cs_ipvec.m          x(p)=b
cs_ipvec_mex.c
cs_maxtransr.m      recursive maximum matching algorithm
cs_maxtransr_mex.c
cs_pvec.m           x=b(p)
cs_pvec_mex.c       interface for cs_pvec
cs_reach.m          non-recursive reach (interface to CSparse cs_reach)
cs_reach_mex.c      non-recursive x=spones(L\sparse(b))
cs_reachr.m         recursive reach (interface to CSparse cs_reachr)
cs_reachr_mex.c
cs_rowcnt.m         row counts for sparse Cholesky
cs_rowcnt_mex.c     row counts for sparse Cholesky
cs_sparse2.m        same as cs_sparse, to test cs_entry function
cs_sparse2_mex.c    like cs_sparse, but for testing cs_entry

cs_test_make.m      compiles MATLAB tests

check_if_same.m     check if two inputs are identical or not
choldn.m            Cholesky downdate
cholup.m            Cholesky update, using Given's rotations
cholupdown.m        Cholesky update/downdate (Bischof, Pan, and Tang method)
cs_q1.m             construct Q from Householder vectors
cs_test_make.m      compiles the CSparse, Demo, and Test mexFunctions.
dmperm_test.m       test cs_dmperm
chol_example.m      simple Cholesky factorization example
etree_sample.m      construct a sample etree and symbolic factorization
gqr3.m              QR factorization, based on Givens rotations
happly.m            apply Householder reflection to a vector
hmake1.m            construct a Householder reflection
mynormest1.m        estimate norm(A,1), using LU factorization (L*U = P*A*Q).
myqr.m              QR factorization using Householder reflections
another_colormap.m  try another color map
cspy_test.m         test cspy and cs_dmspy
qr2.m               QR factorization based on Householder reflections
sample_colormap.m   try a colormap for use in cspy
signum.m            compute and display the sign of a column vector x
sqr_example.m       test cs_sqr
dmspy_test.m        test cspy, cs_dmspy, and cs_dmperm
test_qr.m           test various QR factorization methods
test_randperms.m    test random permutations
testh.m             test Householder reflections
test_qr1.m          test QR factorizations
test_qrsol.m        test cs_qrsol
test_sep.m          test cs_sep, and compare with Gilbert's meshpart vtxsep
testall.m           test all CSparse functions (run tests 1 to 28 below)
test1.m             test cs_transpose
test2.m             test cs_sparse
test3.m             test cs_lsolve, cs_ltsolve, cs_usolve, cs_chol
test4.m             test cs_multiply
test5.m             test cs_add
test6.m             test cs_reach, cs_reachr, cs_lsolve, cs_usolve
test7.m             test cs_lu
test8.m             test cs_cholsol, cs_lusol
test9.m             test cs_qr
test10.m            test cs_qr
test11.m            test cs_rowcnt
test12.m            test cs_qr and compare with svd
test13.m            test cs_counts, cs_etree
test14.m            test cs_droptol
test15.m            test cs_amd
test16.m            test cs_amd
test17.m            test cs_qr, cs_qright, cs_q1, cs_qrleft, cs_qrsol
test18.m            test iterative refinement after backslash
test19.m            test cs_dmperm, cs_maxtransr, cs_dmspy, cs_scc
test20.m            test cholupdown
test21.m            test cs_updown
test22.m            test cond1est
test23.m            test cs_dmspy
test24.m            test cs_fielder
test25.m            test cs_nd
test26.m            test cs_dmsol and cs_dmspy
test27.m            test cs_qr, cs_utsolve, cs_qrsol
test28.m            test cs_randperm, cs_dmperm


--------------------------------------------------------------------------------
./MATLAB/UFget:     MATLAB interface for the UF Sparse Matrix Collection
--------------------------------------------------------------------------------

Contents.m          Contents of UFget
mat/                default directory where downloaded matrices will be put
README.txt          README file for UFget
UFget_defaults.m    default parameter settings
UFget_example.m     example of use
UFget_install.m     installs UFget temporarily (for current session)
UFget_java.class    read a url and load it in into MATLAB (compiled Java code)
UFget_java.java     read a url and load it in into MATLAB (Java source code)
UFget_lookup.m      look up a matrix in the index
UFget.m             UFget itself (primary user interface)
UFweb.m             open url for a matrix or collection
mat/UF_Index.mat    index of matrices in UF Sparse Matrix Collection


--------------------------------------------------------------------------------
./Matrix:           Sample matrices, most from Rutherford/Boeing collection
--------------------------------------------------------------------------------

ash219              overdetermined pattern of Holland survey.  Ashkenazi, 1974.
bcsstk01            stiffness matrix for small generalized eigenvalue problem
bcsstk16            stiffness matrix, Corp of Engineers dam
fs_183_1            unsymmetric facsimile convergence matrix
lp_afiro            NETLIB afiro linear programming problem
mbeacxc             US economy, 1972.  Dan Szyld, while at NYU
t1                  small example used in Chapter 2
west0067            Cavett problem with 5 components (chemical eng., Westerberg)

c_mbeacxc           complex version of mbeacxc
c_west0067          complex version of west0067
mhd1280b            Alfven spectra in magnetohydrodynamics (complex)
neumann             complex matrix
qc324               model of H+ in an electromagnetic field (complex)
t2                  small complex matrix
t3                  small complex matrix
t4                  small complex matrix
c4                  small complex matrix
young1c             aeronautical problem (complex matrix)

--------------------------------------------------------------------------------
./Tcov:             Exhaustive test coverage of CXSparse
--------------------------------------------------------------------------------

covall              same as covall.linux
covall.linux        find coverage (Linux)
covall.sol          find coverage (Solaris)
cov.awk             coverage summary
cover               print uncovered lines
covs                print uncovered lines
cstcov_malloc_test.c    malloc test
cstcov_malloc_test.h
cstcov_test.c       main program for Tcov tests
gcovs               run gcov (Linux)
Makefile            Makefile for Tcov tests
nil                 an empty matrix
zero                a 1-by-1 zero matrix
czero               a 1-by-1 complex zero matrix
README.txt          README file for Tcov directory


--------------------------------------------------------------------------------
Change Log:
--------------------------------------------------------------------------------

Refer to CSparse for changes in CSparse, which are immediately propagated
into CXSparse (those Change Log entries are not repeated here).

Jun 1, 2012.  version 3.1.0

    * now based on CSparse v3.1.0
    * This version of CXSparse changes the 'long' integer from UF_long to
        cs_long_t.  UF_long is still available to user codes, however, so this
        change is backward compatible with user codes.  Future user codes
        should use cs_long_t instead of UF_long.
    * changed unsigned integer in cs_amd.c to signed, for hash code.
    * in Source, only changes are to cs_demo*.c, cs_print.c

Nov 1, 2007.  version 2.2.1

    CXSparse/MATLAB/Test ported to Windows


May 31, 2007.  version 2.2.0

    * back-port to MATLAB 7.2 and earlier (which does not have mwIndex).

    * more graceful failure in cs_make when attempting complex matrix support
        (Windows, in particular)

    * correction to CXSparse/Demo/Makefile

    * added sizeof(CS_INT) printout to cs_idemo.c, cs_ldemo.c

Mar 14, 2007.  Version 2.1.0.

    * MATLAB interface added for CXSparse.

    * cs_complex_t type added (a #define for "double _Complex", which is the
        complex type used in CXSparse 2.0.x).  When compiling with a C++ 
        compiler, the std::compex<double> type is used for the complex case.

    * bug fix in complex sparse Cholesky (cs_chol.c).

    * bug fix in complex sparse Cholesky update/downdate (cs_updown.c).

    * bug fix in cs_symperm for the complex case.

    * "beta" changed from complex to real, in sparse QR (cs_house.c,
        cs_happly.c, cs_qr.c), (a performance/memory improvement, not a
        bug fix).  Similar change to "nz2" in cs_cumsum.c.

May 5, 2006.  Version 2.0.1 released.

    * long changed to UF_long, dependency in ../UFconfig/UFconfig.h added.
        "UF_long" is a #define'd term in UFconfig.h.  It is normally defined
        as "long", but can be redefined as something else if desired.
        On Windows-64, it becomes __int64.

Mar 6, 2006

    "double complex" changed to "double _Complex", to avoid conflicts when
    CXSparse is compiled with a C++ compiler.  Other minor changes to cs.h.
