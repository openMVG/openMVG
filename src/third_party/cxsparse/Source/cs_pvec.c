#include "cs.h"
/* x = b(p), for dense vectors x and b; p=NULL denotes identity */
CS_INT cs_pvec (const CS_INT *p, const CS_ENTRY *b, CS_ENTRY *x, CS_INT n)
{
    CS_INT k ;
    if (!x || !b) return (0) ;                              /* check inputs */
    for (k = 0 ; k < n ; k++) x [k] = b [p ? p [k] : k] ;
    return (1) ;
}
