void linscan_query(UINT32 *counter, UINT32 *res, UINT8 *codes, UINT8 *queries, 
                   int N, UINT32 NQ, int B, unsigned int K, int dim1codes, 
                   int dim1queries);

void linscan_AH_query(REAL *dists, UINT32 *res, UINT8 *codes, REAL *centers,
                      REAL *queries, int N, UINT32 NQ, int B, int K,
                      int dim1codes, int dim1queries, int subdim);
