#include "cs.h"
/* x(p) = b, for dense vectors x and b; p=NULL denotes identity */
CS_INT cs_ipvec (const CS_INT *p, const CS_ENTRY *b, CS_ENTRY *x, CS_INT n)
{
    CS_INT k ;
    if (!x || !b) return (0) ;                              /* check inputs */
    for (k = 0 ; k < n ; k++) x [p ? p [k] : k] = b [k] ;
    return (1) ;
}
