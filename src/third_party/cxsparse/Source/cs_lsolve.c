#include "cs.h"
/* solve Lx=b where x and b are dense.  x=b on input, solution on output. */
CS_INT cs_lsolve (const cs *L, CS_ENTRY *x)
{
    CS_INT p, j, n, *Lp, *Li ;
    CS_ENTRY *Lx ;
    if (!CS_CSC (L) || !x) return (0) ;                     /* check inputs */
    n = L->n ; Lp = L->p ; Li = L->i ; Lx = L->x ;
    for (j = 0 ; j < n ; j++)
    {
        x [j] /= Lx [Lp [j]] ;
        for (p = Lp [j]+1 ; p < Lp [j+1] ; p++)
        {
            x [Li [p]] -= Lx [p] * x [j] ;
        }
    }
    return (1) ;
}
