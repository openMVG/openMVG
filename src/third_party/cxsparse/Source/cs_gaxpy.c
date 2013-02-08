#include "cs.h"
/* y = A*x+y */
CS_INT cs_gaxpy (const cs *A, const CS_ENTRY *x, CS_ENTRY *y)
{
    CS_INT p, j, n, *Ap, *Ai ;
    CS_ENTRY *Ax ;
    if (!CS_CSC (A) || !x || !y) return (0) ;       /* check inputs */
    n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    for (j = 0 ; j < n ; j++)
    {
        for (p = Ap [j] ; p < Ap [j+1] ; p++)
        {
            y [Ai [p]] += Ax [p] * x [j] ;
        }
    }
    return (1) ;
}
