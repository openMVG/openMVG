#include "cudaSift.h"
#include "cudautils.h"

//================= Device matching functions =====================//

__global__ void MatchSiftPoints(SiftPoint *sift1, SiftPoint *sift2, float *corrData, int numPts1, int numPts2)
{
  __shared__ float siftPoint[128];
  __shared__ float sums[16*16];
  const int tx = threadIdx.x;
  const int ty = threadIdx.y;
  const int p1 = blockIdx.x;
  const int p2 = blockIdx.y*16 + ty;
  const float *ptr1 = sift1[p1].data;
  const float *ptr2 = sift2[p2].data;
  const int i = 16*ty + tx;
  if (ty<8)
    siftPoint[i] = ptr1[i];
  __syncthreads();
  float sum = 0.0f;
  if (p2<numPts2)
    for (int j=0;j<8;j++)
      sum += siftPoint[16*j+tx] * ptr2[16*j+tx];
  sums[i] = sum;
  __syncthreads();
  if (tx<8)
    sums[i] += sums[i+8];
  __syncthreads();
  if (tx<4)
    sums[i] += sums[i+4];
  __syncthreads();
  if (ty==0) {
    sum = sums[16*tx+0] + sums[16*tx+1] + sums[16*tx+2] + sums[16*tx+3];
    corrData[p1*gridDim.y*16 + blockIdx.y*16 + tx] = sum;
  }
  __syncthreads();
}


__global__ void MatchSiftPoints2(SiftPoint *sift1, SiftPoint *sift2, float *corrData, int numPts1, int numPts2)
{
  __shared__ float siftPoints1[16*128];
  __shared__ float siftPoints2[16*128];
  const int tx = threadIdx.x;
  const int ty = threadIdx.y;
  const float *ptr1 = sift1[min(numPts1-1,blockIdx.x*16 + ty)].data;
  const float *ptr2 = sift2[min(numPts2-1,blockIdx.y*16 + ty)].data;
  for (int i=0;i<8;i++) {
    siftPoints1[128*ty+16*i+tx] = ptr1[16*i+tx];
    siftPoints2[128*ty+16*i+tx] = ptr2[16*i+tx];
  }
  __syncthreads();
  const int p1 = blockIdx.x*16 + ty;
  const int p2 = blockIdx.y*16 + tx;
  const float *pt1 = &siftPoints1[ty*128];
  const float *pt2 = &siftPoints2[tx*128];
  float sum = 0.0f;
  for (int i=0;i<128;i++) {
    int itx = (i + tx)&127; // avoid bank conflicts
    sum += pt1[itx]*pt2[itx];
  }
  if (p1<numPts1)
    corrData[p1*gridDim.y*16 + p2] = (p2<numPts2 ? sum : -1.0f);
}

__global__ void FindMaxCorr(float *corrData, SiftPoint *sift1, SiftPoint *sift2, int numPts1, int corrWidth, int siftSize)
{
  __shared__ float maxScore[16*16];
  __shared__ float maxScor2[16*16];
  __shared__ int maxIndex[16*16];
  const int tx = threadIdx.x;
  const int ty = threadIdx.y;
  const int idx = ty*16 + tx;
  int p1 = blockIdx.x*16 + threadIdx.y;
  p1 = (p1>=numPts1 ? numPts1-1 : p1);
  maxScore[idx] = -1.0f;
  maxScor2[idx] = -1.0f;
  maxIndex[idx] = -1;
  __syncthreads();
  float *corrs = &corrData[p1*corrWidth];
  for (int i=tx;i<corrWidth;i+=16) {
    float val = corrs[i];
    if (val>maxScore[idx]) {
      maxScor2[idx] = maxScore[idx];
      maxScore[idx] = val;
      maxIndex[idx] = i;
    } else if (val>maxScor2[idx])
      maxScor2[idx] = val;
  }
  //if (p1==1)
  //  printf("tx = %d, score = %.2f, scor2 = %.2f, index = %d\n", 
  //	   tx, maxScore[idx], maxScor2[idx], maxIndex[idx]);
  __syncthreads();
  for (int len=8;len>0;len/=2) {
    if (tx<8) {
      float val = maxScore[idx+len];
      int i = maxIndex[idx+len];
      if (val>maxScore[idx]) {
	maxScor2[idx] = maxScore[idx];
	maxScore[idx] = val;
	maxIndex[idx] = i;
      } else if (val>maxScor2[idx])
	maxScor2[idx] = val;
      float va2 = maxScor2[idx+len];
      if (va2>maxScor2[idx])
	maxScor2[idx] = va2;
    }
    __syncthreads();
    //if (p1==1 && tx<len) 
    //  printf("tx = %d, score = %.2f, scor2 = %.2f, index = %d\n", 
    //	     tx, maxScore[idx], maxScor2[idx], maxIndex[idx]);
  }
  if (tx==6)
    sift1[p1].score = maxScore[ty*16];
  if (tx==7)
    sift1[p1].ambiguity = maxScor2[ty*16] / (maxScore[ty*16] + 1e-6);
  if (tx==8)
    sift1[p1].match = maxIndex[ty*16];
  if (tx==9)
    sift1[p1].match_xpos = sift2[maxIndex[ty*16]].xpos;
  if (tx==10)
    sift1[p1].match_ypos = sift2[maxIndex[ty*16]].ypos;
  __syncthreads();
  //if (tx==0)
  //  printf("index = %d/%d, score = %.2f, ambiguity = %.2f, match = %d\n", 
  //	p1, numPts1, sift1[p1].score, sift1[p1].ambiguity, sift1[p1].match);
}

template <int size>
__device__ void InvertMatrix(float elem[size][size], float res[size][size]) 
{  
  int indx[size];
  float b[size];
  float vv[size];
  for (int i=0;i<size;i++)
    indx[i] = 0;
  int imax = 0;
  float d = 1.0;
  for (int i=0;i<size;i++) { // find biggest element for each row
    float big = 0.0;
    for (int j=0;j<size;j++) {
      float temp = fabs(elem[i][j]); 
      if (temp>big) 
	big = temp;
    }
    if (big>0.0)
      vv[i] = 1.0/big;
    else
      vv[i] = 1e16;
  }
  for (int j=0;j<size;j++) { 
    for (int i=0;i<j;i++) { // i<j
      float sum = elem[i][j]; // i<j (lower left)
      for (int k=0;k<i;k++) // k<i<j
	sum -= elem[i][k]*elem[k][j]; // i>k (upper right), k<j (lower left)
      elem[i][j] = sum; // i<j (lower left)
    }
    float big = 0.0;
    for (int i=j;i<size;i++) { // i>=j
      float sum = elem[i][j]; // i>=j (upper right)
      for (int k=0;k<j;k++) // k<j<=i
	sum -= elem[i][k]*elem[k][j]; // i>k (upper right), k<j (lower left)
      elem[i][j] = sum; // i>=j (upper right)
      float dum = vv[i]*fabs(sum);
      if (dum>=big) {
	big = dum;
	imax = i;  
      }
    }
    if (j!=imax) { // imax>j
      for (int k=0;k<size;k++) {
	float dum = elem[imax][k]; // upper right and lower left
	elem[imax][k] = elem[j][k];
	elem[j][k] = dum;
      }
      d = -d;
      vv[imax] = vv[j];
    }
    indx[j] = imax;
    if (elem[j][j]==0.0)  // j==j (upper right)
      elem[j][j] = 1e-16;
    if (j!=(size-1)) {
      float dum = 1.0/elem[j][j];
      for (int i=j+1;i<size;i++) // i>j
	elem[i][j] *= dum; // i>j (upper right)
    }
  }
  for (int j=0;j<size;j++) {
    for (int k=0;k<size;k++) 
      b[k] = 0.0;  
    b[j] = 1.0;
    int ii = -1;
    for (int i=0;i<size;i++) {
      int ip = indx[i];
      float sum = b[ip];
      b[ip] = b[i];
      if (ii!=-1)
	for (int j=ii;j<i;j++) 
	  sum -= elem[i][j]*b[j]; // i>j (upper right)
      else if (sum!=0.0)
        ii = i;
      b[i] = sum;
    }
    for (int i=size-1;i>=0;i--) {
      float sum = b[i];
      for (int j=i+1;j<size;j++) 
	sum -= elem[i][j]*b[j]; // i<j (lower left)
      b[i] = sum/elem[i][i]; // i==i (upper right)
    }
    for (int i=0;i<size;i++)
      res[i][j] = b[i];
  }
}

__global__ void ComputeHomographies(float *coord, int *randPts, float *homo, 
  int numPts) 
{
  float a[8][8], ia[8][8];
  float b[8]; 
  const int bx = blockIdx.x;
  const int tx = threadIdx.x;
  const int idx = blockDim.x*bx + tx;
  const int numLoops = blockDim.x*gridDim.x;
  for (int i=0;i<4;i++) {
    int pt = randPts[i*numLoops+idx];
    float x1 = coord[pt+0*numPts];
    float y1 = coord[pt+1*numPts];
    float x2 = coord[pt+2*numPts];
    float y2 = coord[pt+3*numPts];
    float *row1 = a[2*i+0];
    row1[0] = x1;
    row1[1] = y1;
    row1[2] = 1.0;
    row1[3] = row1[4] = row1[5] = 0.0;
    row1[6] = -x2*x1;
    row1[7] = -x2*y1;
    float *row2 = a[2*i+1];
    row2[0] = row2[1] = row2[2] = 0.0;
    row2[3] = x1;
    row2[4] = y1;
    row2[5] = 1.0;
    row2[6] = -y2*x1;
    row2[7] = -y2*y1;
    b[2*i+0] = x2;
    b[2*i+1] = y2;
  }
  InvertMatrix<8>(a, ia);
  __syncthreads();
  for (int j=0;j<8;j++) {
    float sum = 0.0f;
    for (int i=0;i<8;i++) 
      sum += ia[j][i]*b[i];
    homo[j*numLoops+idx] = sum;
  }
  __syncthreads();
}

#define TESTHOMO_TESTS 16 // number of tests per block,  alt. 32, 32
#define TESTHOMO_LOOPS 16 // number of loops per block,  alt.  8, 16 

__global__ void TestHomographies(float *d_coord, float *d_homo, 
  int *d_counts, int numPts, float thresh2)
{
  __shared__ float homo[8*TESTHOMO_LOOPS];
  __shared__ int cnts[TESTHOMO_TESTS*TESTHOMO_LOOPS];
  const int tx = threadIdx.x;
  const int ty = threadIdx.y;
  const int idx = blockIdx.y*blockDim.y + tx;
  const int numLoops = blockDim.y*gridDim.y;
  if (ty<8 && tx<TESTHOMO_LOOPS)
    homo[tx*8+ty] = d_homo[idx+ty*numLoops];
  __syncthreads();
  float a[8];
  for (int i=0;i<8;i++) 
    a[i] = homo[ty*8+i];
  int cnt = 0;
  for (int i=tx;i<numPts;i+=TESTHOMO_TESTS) {
    float x1 = d_coord[i+0*numPts];
    float y1 = d_coord[i+1*numPts];
    float x2 = d_coord[i+2*numPts];
    float y2 = d_coord[i+3*numPts];
    float nomx = __fmul_rz(a[0],x1) + __fmul_rz(a[1],y1) + a[2];
    float nomy = __fmul_rz(a[3],x1) + __fmul_rz(a[4],y1) + a[5];
    float deno = __fmul_rz(a[6],x1) + __fmul_rz(a[7],y1) + 1.0f;
    float errx = __fmul_rz(x2,deno) - nomx;
    float erry = __fmul_rz(y2,deno) - nomy;
    float err2 = __fmul_rz(errx,errx) + __fmul_rz(erry,erry);
    if (err2<__fmul_rz(thresh2,__fmul_rz(deno,deno)))
      cnt ++;
  }
  int kty = TESTHOMO_TESTS*ty;
  cnts[kty + tx] = cnt;
  __syncthreads();
  int len = TESTHOMO_TESTS/2;
  while (len>0) {
    if (tx<len)
      cnts[kty + tx] += cnts[kty + tx + len];
    len /= 2;
    __syncthreads();
  }
  if (tx<TESTHOMO_LOOPS && ty==0)
    d_counts[idx] = cnts[TESTHOMO_TESTS*tx];
  __syncthreads();
}

//================= Host matching functions =====================//

double FindHomography(SiftData &data, float *homography, int *numMatches, int numLoops, float minScore, float maxAmbiguity, float thresh)
{
  *numMatches = 0;
  homography[0] = homography[4] = homography[8] = 1.0f;
  homography[1] = homography[2] = homography[3] = 0.0f;
  homography[5] = homography[6] = homography[7] = 0.0f;
#ifdef MANAGEDMEM
  SiftPoint *d_sift = data.m_data;
#else
  if (data.d_data==NULL)
    return 0.0f;
  SiftPoint *d_sift = data.d_data;
#endif
  TimerGPU timer(0);
  numLoops = iDivUp(numLoops,16)*16;
  int numPts = data.numPts;
  if (numPts<8)
    return 0.0f;
  int numPtsUp = iDivUp(numPts, 16)*16;
  float *d_coord, *d_homo;
  int *d_randPts, *h_randPts;
  int randSize = 4*sizeof(int)*numLoops;
  int szFl = sizeof(float);
  int szPt = sizeof(SiftPoint);
  safeCall(cudaMalloc((void **)&d_coord, 4*sizeof(float)*numPtsUp));
  safeCall(cudaMalloc((void **)&d_randPts, randSize));
  safeCall(cudaMalloc((void **)&d_homo, 8*sizeof(float)*numLoops));
  h_randPts = (int*)malloc(randSize);
  float *h_scores = (float *)malloc(sizeof(float)*numPtsUp);
  float *h_ambiguities = (float *)malloc(sizeof(float)*numPtsUp);
  safeCall(cudaMemcpy2D(h_scores, szFl, &d_sift[0].score, szPt, szFl, numPts, cudaMemcpyDeviceToHost));
  safeCall(cudaMemcpy2D(h_ambiguities, szFl, &d_sift[0].ambiguity, szPt, szFl, numPts, cudaMemcpyDeviceToHost));
  int *validPts = (int *)malloc(sizeof(int)*numPts);
  int numValid = 0;
  for (int i=0;i<numPts;i++) {
    if (h_scores[i]>minScore && h_ambiguities[i]<maxAmbiguity)
      validPts[numValid++] = i;
  }
  free(h_scores);
  free(h_ambiguities);
  if (numValid>=8) {
    for (int i=0;i<numLoops;i++) {
      int p1 = rand() % numValid;
      int p2 = rand() % numValid;
      int p3 = rand() % numValid;
      int p4 = rand() % numValid;
      while (p2==p1) p2 = rand() % numValid;
      while (p3==p1 || p3==p2) p3 = rand() % numValid;
      while (p4==p1 || p4==p2 || p4==p3) p4 = rand() % numValid;
      h_randPts[i+0*numLoops] = validPts[p1];
      h_randPts[i+1*numLoops] = validPts[p2];
      h_randPts[i+2*numLoops] = validPts[p3];
      h_randPts[i+3*numLoops] = validPts[p4];
    }
    safeCall(cudaMemcpy(d_randPts, h_randPts, randSize, cudaMemcpyHostToDevice));
    safeCall(cudaMemcpy2D(&d_coord[0*numPtsUp], szFl, &d_sift[0].xpos, szPt, szFl, numPts, cudaMemcpyDeviceToDevice));
    safeCall(cudaMemcpy2D(&d_coord[1*numPtsUp], szFl, &d_sift[0].ypos, szPt, szFl, numPts, cudaMemcpyDeviceToDevice));
    safeCall(cudaMemcpy2D(&d_coord[2*numPtsUp], szFl, &d_sift[0].match_xpos, szPt, szFl, numPts, cudaMemcpyDeviceToDevice));
    safeCall(cudaMemcpy2D(&d_coord[3*numPtsUp], szFl, &d_sift[0].match_ypos, szPt, szFl, numPts, cudaMemcpyDeviceToDevice));
    ComputeHomographies<<<numLoops/16, 16>>>(d_coord, d_randPts, d_homo, numPtsUp);
    safeCall(cudaThreadSynchronize());
    checkMsg("ComputeHomographies() execution failed\n");
    dim3 blocks(1, numLoops/TESTHOMO_LOOPS);
    dim3 threads(TESTHOMO_TESTS, TESTHOMO_LOOPS);
    TestHomographies<<<blocks, threads>>>(d_coord, d_homo, d_randPts, numPtsUp, thresh*thresh);
    safeCall(cudaThreadSynchronize());
    checkMsg("TestHomographies() execution failed\n");
    safeCall(cudaMemcpy(h_randPts, d_randPts, sizeof(int)*numLoops, cudaMemcpyDeviceToHost));
    int maxIndex = -1, maxCount = -1;
    for (int i=0;i<numLoops;i++) 
      if (h_randPts[i]>maxCount) {
	maxCount = h_randPts[i];
	maxIndex = i;
      }
    *numMatches = maxCount;
    safeCall(cudaMemcpy2D(homography, szFl, &d_homo[maxIndex], sizeof(float)*numLoops, szFl, 8, cudaMemcpyDeviceToHost));
  }
  free(validPts);
  free(h_randPts);
  safeCall(cudaFree(d_homo));
  safeCall(cudaFree(d_randPts));
  safeCall(cudaFree(d_coord));
  double gpuTime = timer.read();
#ifdef VERBOSE
  printf("FindHomography time =         %.2f ms\n", gpuTime);
#endif
  return gpuTime;
}

double MatchSiftData(SiftData &data1, SiftData &data2)
{
  TimerGPU timer(0);
  int numPts1 = data1.numPts;
  int numPts2 = data2.numPts;
  if (!numPts1 || !numPts2) 
    return 0.0;
#ifdef MANAGEDMEM
  SiftPoint *sift1 = data1.m_data;
  SiftPoint *sift2 = data2.m_data;
#else
  if (data1.d_data==NULL || data2.d_data==NULL)
    return 0.0f;
  SiftPoint *sift1 = data1.d_data;
  SiftPoint *sift2 = data2.d_data;
#endif
  
  float *d_corrData; 
  int corrWidth = iDivUp(numPts2, 16)*16;
  int corrSize = sizeof(float)*numPts1*corrWidth;
  safeCall(cudaMalloc((void **)&d_corrData, corrSize));
#if 0
  dim3 blocks1(numPts1, iDivUp(numPts2, 16));
  dim3 threads1(16, 16); // each block: 1 points x 16 points
  MatchSiftPoints<<<blocks1, threads1>>>(sift1, sift2, d_corrData, numPts1, numPts2);
#else
  dim3 blocks(iDivUp(numPts1,16), iDivUp(numPts2, 16));
  dim3 threads(16, 16); // each block: 1 points x 16 points
  MatchSiftPoints2<<<blocks, threads>>>(sift1, sift2, d_corrData, numPts1, numPts2);
#endif
  safeCall(cudaThreadSynchronize());
  dim3 blocksMax(iDivUp(numPts1, 16));
  dim3 threadsMax(16, 16);
  FindMaxCorr<<<blocksMax, threadsMax>>>(d_corrData, sift1, sift2, numPts1, corrWidth, sizeof(SiftPoint));
  safeCall(cudaThreadSynchronize());
  checkMsg("MatchSiftPoints() execution failed\n");
  safeCall(cudaFree(d_corrData));
  if (data1.h_data!=NULL) {
    float *h_ptr = &data1.h_data[0].score;
    float *d_ptr = &data1.d_data[0].score;
    safeCall(cudaMemcpy2D(h_ptr, sizeof(SiftPoint), d_ptr, sizeof(SiftPoint), 5*sizeof(float), data1.numPts, cudaMemcpyDeviceToHost));
  }

  double gpuTime = timer.read();
#ifdef VERBOSE
  printf("MatchSiftData time =          %.2f ms\n", gpuTime);
#endif
  return gpuTime;
}		 
  
