#include "cs.h"
/* solve Ux=b where x and b are dense.  x=b on input, solution on output. */
CS_INT cs_usolve (const cs *U, CS_ENTRY *x)
{
    CS_INT p, j, n, *Up, *Ui ;
    CS_ENTRY *Ux ;
    if (!CS_CSC (U) || !x) return (0) ;                     /* check inputs */
    n = U->n ; Up = U->p ; Ui = U->i ; Ux = U->x ;
    for (j = n-1 ; j >= 0 ; j--)
    {
        x [j] /= Ux [Up [j+1]-1] ;
        for (p = Up [j] ; p < Up [j+1]-1 ; p++)
        {
            x [Ui [p]] -= Ux [p] * x [j] ;
        }
    }
    return (1) ;
}
