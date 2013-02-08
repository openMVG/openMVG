#include "cs.h"
/* solve L'x=b where x and b are dense.  x=b on input, solution on output. */
CS_INT cs_ltsolve (const cs *L, CS_ENTRY *x)
{
    CS_INT p, j, n, *Lp, *Li ;
    CS_ENTRY *Lx ;
    if (!CS_CSC (L) || !x) return (0) ;                     /* check inputs */
    n = L->n ; Lp = L->p ; Li = L->i ; Lx = L->x ;
    for (j = n-1 ; j >= 0 ; j--)
    {
        for (p = Lp [j]+1 ; p < Lp [j+1] ; p++)
        {
            x [j] -= CS_CONJ (Lx [p]) * x [Li [p]] ;
        }
        x [j] /= CS_CONJ (Lx [Lp [j]]) ;
    }
    return (1) ;
}
