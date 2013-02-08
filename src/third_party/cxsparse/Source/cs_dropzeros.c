#include "cs.h"
static CS_INT cs_nonzero (CS_INT i, CS_INT j, CS_ENTRY aij, void *other)
{
    return (aij != 0) ;
}
CS_INT cs_dropzeros (cs *A)
{
    return (cs_fkeep (A, &cs_nonzero, NULL)) ;  /* keep all nonzero entries */
} 
