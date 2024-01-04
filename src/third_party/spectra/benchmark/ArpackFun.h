#ifndef SPECTRA_ARPACKFUN_H
#define SPECTRA_ARPACKFUN_H

#define F77_CALL(x) x##_
#define F77_NAME(x) F77_CALL(x)

enum BMAT
{
    BMAT_I = 0,
    BMAT_G
};

enum WHICH
{
    WHICH_LM = 0,
    WHICH_SM,
    WHICH_LR,
    WHICH_SR,
    WHICH_LI,
    WHICH_SI,
    WHICH_LA,
    WHICH_SA,
    WHICH_BE
};

enum HOWMNY
{
    HOWMNY_A = 0,
    HOWMNY_P,
    HOWMNY_S
};

extern "C" {

// ARPACK Fortran functions
void F77_NAME(dsaupdwr)(int *ido, int *bmati, int *n, int *whichi,
                        int *nev, double *tol, double *resid,
                        int *ncv, double *v, int *ldv,
                        int *iparam, int *ipntr, double *workd,
                        double *workl, int *lworkl, int *info);

void F77_NAME(dseupdwr)(int *rvec, int *howmnyi, int *select, double *d,
                        double *z, int *ldz, double *sigma, int *bmati,
                        int *n, int *whichi, int *nev, double *tol,
                        double *resid, int *ncv, double *v, int *ldv,
                        int *iparam, int *ipntr, double *workd, double *workl,
                        int *lworkl, int *info);

void F77_NAME(dnaupdwr)(int *ido, int *bmati, int *n, int *whichi,
                        int *nev, double *tol, double *resid,
                        int *ncv, double *v, int *ldv,
                        int *iparam, int *ipntr, double *workd,
                        double *workl, int *lworkl, int *info);

void F77_NAME(dneupdwr)(int *rvec, int *howmnyi, int *select, double *dr, double *di,
                        double *z, int *ldz, double *sigmar, double *sigmai, double *workev,
                        int *bmati, int *n, int *whichi, int *nev, double *tol,
                        double *resid, int *ncv, double *v, int *ldv, int *iparam,
                        int *ipntr, double *workd, double *workl, int *lworkl, int *info);

}  // extern "C"

// Map char *which to enum type
// WHICH_LM is the default if unusual case happens
inline int whichenum(char *which)
{
    switch (which[0])
    {
        case 'L':
            switch (which[1])
            {
                case 'M':
                    return (int) WHICH_LM;
                case 'R':
                    return (int) WHICH_LR;
                case 'I':
                    return (int) WHICH_LI;
                case 'A':
                    return (int) WHICH_LA;
                default:
                    return (int) WHICH_LM;
            }
        case 'S':
            switch (which[1])
            {
                case 'M':
                    return (int) WHICH_SM;
                case 'R':
                    return (int) WHICH_SR;
                case 'I':
                    return (int) WHICH_SI;
                case 'A':
                    return (int) WHICH_SA;
                default:
                    return (int) WHICH_LM;
            }
        case 'B':
            if (which[1] == 'E')
                return (int) WHICH_BE;
            else
                return (int) WHICH_LM;
        default:
            return (int) WHICH_LM;
    }

    return (int) WHICH_LM;
}

// C++ Wrapper of the functions above
inline void saupd(int &ido, char bmat, int n, char *which,
                  int nev, double &tol, double resid[],
                  int ncv, double v[], int ldv,
                  int iparam[], int ipntr[], double workd[],
                  double workl[], int lworkl, int &info)
{
    int bmati = (bmat == 'G') ? BMAT_G : BMAT_I;
    int whichi = whichenum(which);

    F77_CALL(dsaupdwr)
    (&ido, &bmati, &n, &whichi,
     &nev, &tol, resid,
     &ncv, v, &ldv,
     iparam, ipntr, workd,
     workl, &lworkl, &info);
}

inline void seupd(bool rvec, char howmny, double d[],
                  double z[], int ldz, double sigma, char bmat,
                  int n, char *which, int nev, double tol,
                  double resid[], int ncv, double v[], int ldv,
                  int iparam[], int ipntr[], double workd[], double workl[],
                  int lworkl, int &info)
{
    int rvec_pass = (int) rvec;
    int *select_pass = new int[ncv];
    double *z_pass = (z == NULL) ? v : z;
    int howmnyi = (howmny == 'P') ? HOWMNY_P : ((howmny == 'S') ? HOWMNY_S : HOWMNY_A);
    int bmati = (bmat == 'G') ? BMAT_G : BMAT_I;
    int whichi = whichenum(which);

    F77_CALL(dseupdwr)
    (&rvec_pass, &howmnyi, select_pass, d,
     z_pass, &ldz, &sigma, &bmati,
     &n, &whichi, &nev, &tol,
     resid, &ncv, v, &ldv,
     iparam, ipntr, workd, workl,
     &lworkl, &info);

    delete[] select_pass;
}

inline void naupd(int &ido, char bmat, int n, char *which,
                  int nev, double &tol, double resid[],
                  int ncv, double v[], int ldv,
                  int iparam[], int ipntr[], double workd[],
                  double workl[], int lworkl, int &info)
{
    int bmati = (bmat == 'G') ? BMAT_G : BMAT_I;
    int whichi = whichenum(which);

    F77_CALL(dnaupdwr)
    (&ido, &bmati, &n, &whichi,
     &nev, &tol, resid,
     &ncv, v, &ldv,
     iparam, ipntr, workd,
     workl, &lworkl, &info);
}

inline void neupd(bool rvec, char howmny, double dr[], double di[],
                  double z[], int ldz, double sigmar, double sigmai, double workev[],
                  char bmat, int n, char *which, int nev, double tol,
                  double resid[], int ncv, double v[], int ldv, int iparam[],
                  int ipntr[], double workd[], double workl[], int lworkl, int &info)
{
    int rvec_pass = (int) rvec;
    int *select_pass = new int[ncv];
    double *z_pass = (z == NULL) ? v : z;
    int howmnyi = (howmny == 'P') ? HOWMNY_P : ((howmny == 'S') ? HOWMNY_S : HOWMNY_A);
    int bmati = (bmat == 'G') ? BMAT_G : BMAT_I;
    int whichi = whichenum(which);

    F77_CALL(dneupdwr)
    (&rvec_pass, &howmnyi, select_pass, dr, di,
     z_pass, &ldz, &sigmar, &sigmai, workev,
     &bmati, &n, &whichi, &nev, &tol,
     resid, &ncv, v, &ldv, iparam,
     ipntr, workd, workl, &lworkl, &info);

    delete[] select_pass;
}

#endif  // SPECTRA_ARPACKFUN_H
