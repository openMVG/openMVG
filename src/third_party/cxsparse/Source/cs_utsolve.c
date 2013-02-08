#include "cs.h"
/* solve U'x=b where x and b are dense.  x=b on input, solution on output. */
CS_INT cs_utsolve (const cs *U, CS_ENTRY *x)
{
    CS_INT p, j, n, *Up, *Ui ;
    CS_ENTRY *Ux ;
    if (!CS_CSC (U) || !x) return (0) ;                     /* check inputs */
    n = U->n ; Up = U->p ; Ui = U->i ; Ux = U->x ;
    for (j = 0 ; j < n ; j++)
    {
        for (p = Up [j] ; p < Up [j+1]-1 ; p++)
        {
            x [j] -= CS_CONJ (Ux [p]) * x [Ui [p]] ;
        }
        x [j] /= CS_CONJ (Ux [Up [j+1]-1]) ;
    }
    return (1) ;
}
