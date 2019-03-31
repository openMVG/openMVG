//********************************************************//
// CUDA SIFT extractor by Marten Bjorkman aka Celebrandil //
//********************************************************//  

#include "cudautils.h"
#include "cudaSiftD.h"
#include "cudaSift.h"

///////////////////////////////////////////////////////////////////////////////
// Kernel configuration
///////////////////////////////////////////////////////////////////////////////

__constant__ int d_MaxNumPoints;
__device__ unsigned int d_PointCounter[8*2+1];
__constant__ float d_ScaleDownKernel[5]; 
__constant__ float d_LowPassKernel[2*LOWPASS_R+1]; 
__constant__ float d_LaplaceKernel[8*12*16]; 

///////////////////////////////////////////////////////////////////////////////
// Lowpass filter and subsample image
///////////////////////////////////////////////////////////////////////////////
__global__ void ScaleDownDenseShift(float *d_Result, float *d_Data, int width, int pitch, int height, int newpitch)
{
#define BW (SCALEDOWN_W+4)
#define BH (SCALEDOWN_H+4)
#define W2 (SCALEDOWN_W/2)
#define H2 (SCALEDOWN_H/2)
  __shared__ float brows[BH*BW];
  const int tx = threadIdx.x;
  const int ty = threadIdx.y;
  const int xp = blockIdx.x*SCALEDOWN_W + tx;
  const int yp = blockIdx.y*SCALEDOWN_H + ty;
  const float k0 = d_ScaleDownKernel[0];
  const float k1 = d_ScaleDownKernel[1];
  const float k2 = d_ScaleDownKernel[2];
  const int xl = min(width-1,  max(0, xp-2));
  const int yl = min(height-1, max(0, yp-2));
  if (xp<(width+4) && yp<(height+4)) {
    float v = d_Data[yl*pitch + xl];
    brows[BW*ty + tx]  = k0*(v + ShiftDown(v, 4)) + k1*(ShiftDown(v, 1) + ShiftDown(v, 3)) + k2*ShiftDown(v, 2);
  }
  __syncthreads();
  const int xs = blockIdx.x*W2 + tx;
  const int ys = blockIdx.y*H2 + ty;
  if (tx<W2 && ty<H2 && xs<(width/2) && ys<(height/2)) {
    float *ptr = &brows[BW*(ty*2) + (tx*2)];
    d_Result[ys*newpitch + xs] = k0*(ptr[0] + ptr[4*BW]) + k1*(ptr[1*BW] + ptr[3*BW]) + k2*ptr[2*BW];
  } 
}

__global__ void ScaleDownDense(float *d_Result, float *d_Data, int width, int pitch, int height, int newpitch)
{
#define BW (SCALEDOWN_W+4)
#define BH (SCALEDOWN_H+4)
#define W2 (SCALEDOWN_W/2)
#define H2 (SCALEDOWN_H/2)
  __shared__ float irows[BH*BW]; 
  __shared__ float brows[BH*W2];
  const int tx = threadIdx.x;
  const int ty = threadIdx.y;
  const int xp = blockIdx.x*SCALEDOWN_W + tx;
  const int yp = blockIdx.y*SCALEDOWN_H + ty;
  const int xl = min(width-1,  max(0, xp-2));
  const int yl = min(height-1, max(0, yp-2));
  const float k0 = d_ScaleDownKernel[0];
  const float k1 = d_ScaleDownKernel[1];
  const float k2 = d_ScaleDownKernel[2];
  if (xp<(width+4) && yp<(height+4))
    irows[BW*ty + tx] = d_Data[yl*pitch + xl];
  __syncthreads();
  if (yp<(height+4) && tx<W2) {
    float *ptr = &irows[BW*ty + 2*tx];
    brows[W2*ty + tx] = k0*(ptr[0] + ptr[4]) + k1*(ptr[1] + ptr[3]) + k2*ptr[2];
  }
  __syncthreads();
  const int xs = blockIdx.x*W2 + tx;
  const int ys = blockIdx.y*H2 + ty;
  if (tx<W2 && ty<H2 && xs<(width/2) && ys<(height/2)) {
    float *ptr = &brows[W2*(ty*2) + tx];
    d_Result[ys*newpitch + xs] = k0*(ptr[0] + ptr[4*W2]) + k1*(ptr[1*W2] + ptr[3*W2]) + k2*ptr[2*W2];
  } 
}

__global__ void ScaleDown(float *d_Result, float *d_Data, int width, int pitch, int height, int newpitch)
{
  __shared__ float inrow[SCALEDOWN_W+4]; 
  __shared__ float brow[5*(SCALEDOWN_W/2)];
  __shared__ int yRead[SCALEDOWN_H+4];
  __shared__ int yWrite[SCALEDOWN_H+4];
  #define dx2 (SCALEDOWN_W/2)
  const int tx = threadIdx.x;
  const int tx0 = tx + 0*dx2;
  const int tx1 = tx + 1*dx2;
  const int tx2 = tx + 2*dx2;
  const int tx3 = tx + 3*dx2;
  const int tx4 = tx + 4*dx2;
  const int xStart = blockIdx.x*SCALEDOWN_W;
  const int yStart = blockIdx.y*SCALEDOWN_H;
  const int xWrite = xStart/2 + tx;
  float k0 = d_ScaleDownKernel[0];
  float k1 = d_ScaleDownKernel[1];
  float k2 = d_ScaleDownKernel[2];
  if (tx<SCALEDOWN_H+4) {
    int y = yStart + tx - 2; 
    y = (y<0 ? 0 : y);
    y = (y>=height ? height-1 : y);
    yRead[tx] = y*pitch;
    yWrite[tx] = (yStart + tx - 4)/2 * newpitch;
  }
  __syncthreads();
  int xRead = xStart + tx - 2;
  xRead = (xRead<0 ? 0 : xRead);
  xRead = (xRead>=width ? width-1 : xRead);

  int maxtx = min(dx2, width/2 - xStart/2);
  for (int dy=0;dy<SCALEDOWN_H+4;dy+=5) {
    {
      inrow[tx] = d_Data[yRead[dy+0] + xRead];
      __syncthreads();
      if (tx<maxtx) {
	brow[tx4] = k0*(inrow[2*tx]+inrow[2*tx+4]) + k1*(inrow[2*tx+1]+inrow[2*tx+3]) + k2*inrow[2*tx+2];
	if (dy>=4 && !(dy&1))
	  d_Result[yWrite[dy+0] + xWrite] = k2*brow[tx2] + k0*(brow[tx0]+brow[tx4]) + k1*(brow[tx1]+brow[tx3]);
      }
      __syncthreads();
    }
    if (dy<(SCALEDOWN_H+3)) {
      inrow[tx] = d_Data[yRead[dy+1] + xRead];
      __syncthreads();
      if (tx<maxtx) {
	brow[tx0] = k0*(inrow[2*tx]+inrow[2*tx+4]) + k1*(inrow[2*tx+1]+inrow[2*tx+3]) + k2*inrow[2*tx+2];
	if (dy>=3 && (dy&1))
	  d_Result[yWrite[dy+1] + xWrite] = k2*brow[tx3] + k0*(brow[tx1]+brow[tx0]) + k1*(brow[tx2]+brow[tx4]);
      }
      __syncthreads();
    }
    if (dy<(SCALEDOWN_H+2)) {
      inrow[tx] = d_Data[yRead[dy+2] + xRead];
      __syncthreads();
      if (tx<maxtx) {
	brow[tx1] = k0*(inrow[2*tx]+inrow[2*tx+4]) + k1*(inrow[2*tx+1]+inrow[2*tx+3]) + k2*inrow[2*tx+2];
	if (dy>=2 && !(dy&1))
	  d_Result[yWrite[dy+2] + xWrite] = k2*brow[tx4] + k0*(brow[tx2]+brow[tx1]) + k1*(brow[tx3]+brow[tx0]);
      }
      __syncthreads();
    }
    if (dy<(SCALEDOWN_H+1)) {
      inrow[tx] = d_Data[yRead[dy+3] + xRead];
      __syncthreads();
      if (tx<maxtx) {
	brow[tx2] = k0*(inrow[2*tx]+inrow[2*tx+4]) + k1*(inrow[2*tx+1]+inrow[2*tx+3]) + k2*inrow[2*tx+2];
	if (dy>=1 && (dy&1))
	  d_Result[yWrite[dy+3] + xWrite] = k2*brow[tx0] + k0*(brow[tx3]+brow[tx2]) + k1*(brow[tx4]+brow[tx1]);
      }
      __syncthreads();
    }
    if (dy<SCALEDOWN_H) {
      inrow[tx] = d_Data[yRead[dy+4] + xRead];
      __syncthreads();
      if (tx<dx2 && xWrite<width/2) {
	brow[tx3] = k0*(inrow[2*tx]+inrow[2*tx+4]) + k1*(inrow[2*tx+1]+inrow[2*tx+3]) + k2*inrow[2*tx+2];
	if (!(dy&1))
	  d_Result[yWrite[dy+4] + xWrite] = k2*brow[tx1] + k0*(brow[tx4]+brow[tx3]) + k1*(brow[tx0]+brow[tx2]);
      }
      __syncthreads();
    }
  }
}

__global__ void ScaleUp(float *d_Result, float *d_Data, int width, int pitch, int height, int newpitch)
{
  #undef BW
  #undef BH
  #define BW (SCALEUP_W/2 + 2)
  #define BH (SCALEUP_H/2 + 2)
  __shared__ float buffer[BW*BH];
  const int tx = threadIdx.x;
  const int ty = threadIdx.y;
  if (tx<BW && ty<BH) {
    int x = min(max(blockIdx.x*(SCALEUP_W/2) + tx - 1, 0), width-1);
    int y = min(max(blockIdx.y*(SCALEUP_H/2) + ty - 1, 0), height-1);
    buffer[ty*BW + tx] = d_Data[y*pitch + x];
  }
  __syncthreads();
  int x = blockIdx.x*SCALEUP_W + tx;
  int y = blockIdx.y*SCALEUP_H + ty;
  if (x<2*width && y<2*height) {
    int bx = (tx + 1)/2;
    int by = (ty + 1)/2;
    int bp = by*BW + bx;
    float wx = 0.25f + (tx&1)*0.50f;
    float wy = 0.25f + (ty&1)*0.50f;
    d_Result[y*newpitch + x] = wy*(wx*buffer[bp] + (1.0f-wx)*buffer[bp+1]) +
      (1.0f-wy)*(wx*buffer[bp+BW] + (1.0f-wx)*buffer[bp+BW+1]);
  }
}

__global__ void ExtractSiftDescriptors(cudaTextureObject_t texObj, SiftPoint *d_sift, int fstPts, float subsampling)
{
  __shared__ float gauss[16];
  __shared__ float buffer[128];
  __shared__ float sums[4];

  const int tx = threadIdx.x; // 0 -> 16
  const int ty = threadIdx.y; // 0 -> 8
  const int idx = ty*16 + tx;
  const int bx = blockIdx.x + fstPts;  // 0 -> numPts
  if (ty==0)
    gauss[tx] = exp(-(tx-7.5f)*(tx-7.5f)/128.0f);
  buffer[idx] = 0.0f;
  __syncthreads();

  // Compute angles and gradients
  float theta = 2.0f*3.1415f/360.0f*d_sift[bx].orientation;
  float sina = sinf(theta);           // cosa -sina
  float cosa = cosf(theta);           // sina  cosa
  float scale = 12.0f/16.0f*d_sift[bx].scale;
  float ssina = scale*sina; 
  float scosa = scale*cosa;

  for (int y=ty;y<16;y+=8) {
    float xpos = d_sift[bx].xpos + (tx-7.5f)*scosa - (y-7.5f)*ssina + 0.5f;
    float ypos = d_sift[bx].ypos + (tx-7.5f)*ssina + (y-7.5f)*scosa + 0.5f;
    float dx = tex2D<float>(texObj, xpos+cosa, ypos+sina) - 
      tex2D<float>(texObj, xpos-cosa, ypos-sina);
    float dy = tex2D<float>(texObj, xpos-sina, ypos+cosa) - 
      tex2D<float>(texObj, xpos+sina, ypos-cosa);
    float grad = gauss[y]*gauss[tx] * sqrtf(dx*dx + dy*dy);
    float angf = 4.0f/3.1415f*atan2f(dy, dx) + 4.0f;
    
    int hori = (tx + 2)/4 - 1;      // Convert from (tx,y,angle) to bins      
    float horf = (tx - 1.5f)/4.0f - hori;
    float ihorf = 1.0f - horf;           
    int veri = (y + 2)/4 - 1;
    float verf = (y - 1.5f)/4.0f - veri;
    float iverf = 1.0f - verf;
    int angi = angf;
    int angp = (angi<7 ? angi+1 : 0);
    angf -= angi;
    float iangf = 1.0f - angf;
    
    int hist = 8*(4*veri + hori);   // Each gradient measure is interpolated 
    int p1 = angi + hist;           // in angles, xpos and ypos -> 8 stores
    int p2 = angp + hist;
    if (tx>=2) { 
      float grad1 = ihorf*grad;
      if (y>=2) {   // Upper left
        float grad2 = iverf*grad1;
	atomicAdd(buffer + p1, iangf*grad2);
	atomicAdd(buffer + p2,  angf*grad2);
      }
      if (y<=13) {  // Lower left
        float grad2 = verf*grad1;
	atomicAdd(buffer + p1+32, iangf*grad2); 
	atomicAdd(buffer + p2+32,  angf*grad2);
      }
    }
    if (tx<=13) { 
      float grad1 = horf*grad;
      if (y>=2) {    // Upper right
        float grad2 = iverf*grad1;
	atomicAdd(buffer + p1+8, iangf*grad2);
	atomicAdd(buffer + p2+8,  angf*grad2);
      }
      if (y<=13) {   // Lower right
        float grad2 = verf*grad1;
	atomicAdd(buffer + p1+40, iangf*grad2);
	atomicAdd(buffer + p2+40,  angf*grad2);
      }
    }
  }
  __syncthreads();

  // Normalize twice and suppress peaks first time
  float sum = buffer[idx]*buffer[idx];
  for (int i=16;i>0;i/=2)
    sum += ShiftDown(sum, i);
  if ((idx&31)==0)
    sums[idx/32] = sum;
  __syncthreads();
  float tsum1 = sums[0] + sums[1] + sums[2] + sums[3]; 
  tsum1 = min(buffer[idx] * rsqrtf(tsum1), 0.2f);
  
  sum = tsum1*tsum1; 
  for (int i=16;i>0;i/=2)
    sum += ShiftDown(sum, i);
  if ((idx&31)==0)
    sums[idx/32] = sum;
  __syncthreads();

  float tsum2 = sums[0] + sums[1] + sums[2] + sums[3];
  float *desc = d_sift[bx].data;
  desc[idx] = tsum1 * rsqrtf(tsum2);
  if (idx==0) {
    d_sift[bx].xpos *= subsampling;
    d_sift[bx].ypos *= subsampling;
    d_sift[bx].scale *= subsampling;
  }
}
 

__global__ void ExtractSiftDescriptorsCONST(cudaTextureObject_t texObj, SiftPoint *d_sift, float subsampling, int octave)
{
  __shared__ float gauss[16];
  __shared__ float buffer[128];
  __shared__ float sums[4];

  const int tx = threadIdx.x; // 0 -> 16
  const int ty = threadIdx.y; // 0 -> 8
  const int idx = ty*16 + tx;
  if (ty==0)
    gauss[tx] = exp(-(tx-7.5f)*(tx-7.5f)/128.0f);

  int fstPts = min(d_PointCounter[2*octave-1], d_MaxNumPoints);
  int totPts = min(d_PointCounter[2*octave+1], d_MaxNumPoints);
  //if (tx==0 && ty==0)
  //  printf("%d %d %d %d\n", octave, fstPts, min(d_PointCounter[2*octave], d_MaxNumPoints), totPts); 
  for (int bx = blockIdx.x + fstPts; bx < totPts; bx += gridDim.x) {
    
    buffer[idx] = 0.0f;
    __syncthreads();

    // Compute angles and gradients
    float theta = 2.0f*3.1415f/360.0f*d_sift[bx].orientation;
    float sina = sinf(theta);           // cosa -sina
    float cosa = cosf(theta);           // sina  cosa
    float scale = 12.0f/16.0f*d_sift[bx].scale;
    float ssina = scale*sina; 
    float scosa = scale*cosa;
    
    for (int y=ty;y<16;y+=8) {
      float xpos = d_sift[bx].xpos + (tx-7.5f)*scosa - (y-7.5f)*ssina + 0.5f; 
      float ypos = d_sift[bx].ypos + (tx-7.5f)*ssina + (y-7.5f)*scosa + 0.5f;
      float dx = tex2D<float>(texObj, xpos+cosa, ypos+sina) - 
	tex2D<float>(texObj, xpos-cosa, ypos-sina);
      float dy = tex2D<float>(texObj, xpos-sina, ypos+cosa) - 
	tex2D<float>(texObj, xpos+sina, ypos-cosa);
      float grad = gauss[y]*gauss[tx] * sqrtf(dx*dx + dy*dy);
      float angf = 4.0f/3.1415f*atan2f(dy, dx) + 4.0f;
      
      int hori = (tx + 2)/4 - 1;      // Convert from (tx,y,angle) to bins      
      float horf = (tx - 1.5f)/4.0f - hori;
      float ihorf = 1.0f - horf;           
      int veri = (y + 2)/4 - 1;
      float verf = (y - 1.5f)/4.0f - veri;
      float iverf = 1.0f - verf;
      int angi = angf;
      int angp = (angi<7 ? angi+1 : 0);
      angf -= angi;
      float iangf = 1.0f - angf;
      
      int hist = 8*(4*veri + hori);   // Each gradient measure is interpolated 
      int p1 = angi + hist;           // in angles, xpos and ypos -> 8 stores
      int p2 = angp + hist;
      if (tx>=2) { 
	float grad1 = ihorf*grad;
	if (y>=2) {   // Upper left
	  float grad2 = iverf*grad1;
	  atomicAdd(buffer + p1, iangf*grad2);
	  atomicAdd(buffer + p2,  angf*grad2);
	}
	if (y<=13) {  // Lower left
	  float grad2 = verf*grad1;
	  atomicAdd(buffer + p1+32, iangf*grad2); 
	  atomicAdd(buffer + p2+32,  angf*grad2);
	}
      }
      if (tx<=13) { 
	float grad1 = horf*grad;
	if (y>=2) {    // Upper right
	  float grad2 = iverf*grad1;
	  atomicAdd(buffer + p1+8, iangf*grad2);
	  atomicAdd(buffer + p2+8,  angf*grad2);
	}
	if (y<=13) {   // Lower right
	  float grad2 = verf*grad1;
	  atomicAdd(buffer + p1+40, iangf*grad2);
	  atomicAdd(buffer + p2+40,  angf*grad2);
	}
      }
    }
    __syncthreads();
    
    // Normalize twice and suppress peaks first time
    float sum = buffer[idx]*buffer[idx];
    for (int i=16;i>0;i/=2)
      sum += ShiftDown(sum, i);
    if ((idx&31)==0)
      sums[idx/32] = sum;
    __syncthreads();
    float tsum1 = sums[0] + sums[1] + sums[2] + sums[3]; 
    tsum1 = min(buffer[idx] * rsqrtf(tsum1), 0.2f);
     
    sum = tsum1*tsum1; 
    for (int i=16;i>0;i/=2)
      sum += ShiftDown(sum, i);
    if ((idx&31)==0)
      sums[idx/32] = sum;
    __syncthreads();
    
    float tsum2 = sums[0] + sums[1] + sums[2] + sums[3];
    float *desc = d_sift[bx].data;
    desc[idx] = tsum1 * rsqrtf(tsum2);
    if (idx==0) {
      d_sift[bx].xpos *= subsampling;
      d_sift[bx].ypos *= subsampling;
      d_sift[bx].scale *= subsampling;
    }
    __syncthreads();
  }
}
 

__global__ void ExtractSiftDescriptorsOld(cudaTextureObject_t texObj, SiftPoint *d_sift, int fstPts, float subsampling)
{
  __shared__ float gauss[16];
  __shared__ float buffer[128];
  __shared__ float sums[128];

  const int tx = threadIdx.x; // 0 -> 16
  const int ty = threadIdx.y; // 0 -> 8
  const int idx = ty*16 + tx;
  const int bx = blockIdx.x + fstPts;  // 0 -> numPts
  if (ty==0)
    gauss[tx] = exp(-(tx-7.5f)*(tx-7.5f)/128.0f);
  buffer[idx] = 0.0f;
  __syncthreads();

  // Compute angles and gradients
  float theta = 2.0f*3.1415f/360.0f*d_sift[bx].orientation;
  float sina = sinf(theta);           // cosa -sina
  float cosa = cosf(theta);           // sina  cosa
  float scale = 12.0f/16.0f*d_sift[bx].scale;
  float ssina = scale*sina; 
  float scosa = scale*cosa;

  for (int y=ty;y<16;y+=8) {
    float xpos = d_sift[bx].xpos + (tx-7.5f)*scosa - (y-7.5f)*ssina + 0.5f;
    float ypos = d_sift[bx].ypos + (tx-7.5f)*ssina + (y-7.5f)*scosa + 0.5f;
    float dx = tex2D<float>(texObj, xpos+cosa, ypos+sina) - 
      tex2D<float>(texObj, xpos-cosa, ypos-sina);
    float dy = tex2D<float>(texObj, xpos-sina, ypos+cosa) - 
      tex2D<float>(texObj, xpos+sina, ypos-cosa);
    float grad = gauss[y]*gauss[tx] * sqrtf(dx*dx + dy*dy);
    float angf = 4.0f/3.1415f*atan2f(dy, dx) + 4.0f;
    
    int hori = (tx + 2)/4 - 1;      // Convert from (tx,y,angle) to bins      
    float horf = (tx - 1.5f)/4.0f - hori;  
    float ihorf = 1.0f - horf;           
    int veri = (y + 2)/4 - 1;
    float verf = (y - 1.5f)/4.0f - veri;
    float iverf = 1.0f - verf;
    int angi = angf;
    int angp = (angi<7 ? angi+1 : 0);
    angf -= angi;
    float iangf = 1.0f - angf;
    
    int hist = 8*(4*veri + hori);   // Each gradient measure is interpolated 
    int p1 = angi + hist;           // in angles, xpos and ypos -> 8 stores
    int p2 = angp + hist;
    if (tx>=2) { 
      float grad1 = ihorf*grad;
      if (y>=2) {   // Upper left
        float grad2 = iverf*grad1;
	atomicAdd(buffer + p1, iangf*grad2);
	atomicAdd(buffer + p2,  angf*grad2);
      }
      if (y<=13) {  // Lower left
        float grad2 = verf*grad1;
	atomicAdd(buffer + p1+32, iangf*grad2); 
	atomicAdd(buffer + p2+32,  angf*grad2);
      }
    }
    if (tx<=13) { 
      float grad1 = horf*grad;
      if (y>=2) {    // Upper right
        float grad2 = iverf*grad1;
	atomicAdd(buffer + p1+8, iangf*grad2);
	atomicAdd(buffer + p2+8,  angf*grad2);
      }
      if (y<=13) {   // Lower right
        float grad2 = verf*grad1;
	atomicAdd(buffer + p1+40, iangf*grad2);
	atomicAdd(buffer + p2+40,  angf*grad2);
      }
    }
  }
  __syncthreads();

  // Normalize twice and suppress peaks first time
  if (idx<64)
    sums[idx] = buffer[idx]*buffer[idx] + buffer[idx+64]*buffer[idx+64];
  __syncthreads();      
  if (idx<32) sums[idx] = sums[idx] + sums[idx+32];
  __syncthreads();      
  if (idx<16) sums[idx] = sums[idx] + sums[idx+16];
  __syncthreads();      
  if (idx<8)  sums[idx] = sums[idx] + sums[idx+8];
  __syncthreads();      
  if (idx<4)  sums[idx] = sums[idx] + sums[idx+4];
  __syncthreads();      
  float tsum1 = sums[0] + sums[1] + sums[2] + sums[3]; 
  buffer[idx] = buffer[idx] * rsqrtf(tsum1);

  if (buffer[idx]>0.2f)
    buffer[idx] = 0.2f;
  __syncthreads();
  if (idx<64)
    sums[idx] = buffer[idx]*buffer[idx] + buffer[idx+64]*buffer[idx+64];
  __syncthreads();      
  if (idx<32) sums[idx] = sums[idx] + sums[idx+32];
  __syncthreads();      
  if (idx<16) sums[idx] = sums[idx] + sums[idx+16];
  __syncthreads();      
  if (idx<8)  sums[idx] = sums[idx] + sums[idx+8];
  __syncthreads();      
  if (idx<4)  sums[idx] = sums[idx] + sums[idx+4];
  __syncthreads();      
  float tsum2 = sums[0] + sums[1] + sums[2] + sums[3]; 

  float *desc = d_sift[bx].data;
  desc[idx] = buffer[idx] * rsqrtf(tsum2);
  if (idx==0) {
    d_sift[bx].xpos *= subsampling;
    d_sift[bx].ypos *= subsampling;
    d_sift[bx].scale *= subsampling;
  }
}


__device__ void ExtractSiftDescriptor(cudaTextureObject_t texObj, SiftPoint *d_sift, float subsampling, int octave, int bx)
{
  __shared__ float gauss[16];
  __shared__ float buffer[128];
  __shared__ float sums[4];

  const int idx = threadIdx.x;
  const int tx = idx & 15; // 0 -> 16
  const int ty = idx / 16; // 0 -> 8
  if (ty==0)
    gauss[tx] = exp(-(tx-7.5f)*(tx-7.5f)/128.0f);
  buffer[idx] = 0.0f;
  __syncthreads();

  // Compute angles and gradients
  float theta = 2.0f*3.1415f/360.0f*d_sift[bx].orientation;
  float sina = sinf(theta);           // cosa -sina
  float cosa = cosf(theta);           // sina  cosa
  float scale = 12.0f/16.0f*d_sift[bx].scale;
  float ssina = scale*sina; 
  float scosa = scale*cosa;
  
  for (int y=ty;y<16;y+=8) {
    float xpos = d_sift[bx].xpos + (tx-7.5f)*scosa - (y-7.5f)*ssina + 0.5f;
    float ypos = d_sift[bx].ypos + (tx-7.5f)*ssina + (y-7.5f)*scosa + 0.5f;
    float dx = tex2D<float>(texObj, xpos+cosa, ypos+sina) - 
      tex2D<float>(texObj, xpos-cosa, ypos-sina);
    float dy = tex2D<float>(texObj, xpos-sina, ypos+cosa) - 
      tex2D<float>(texObj, xpos+sina, ypos-cosa);
    float grad = gauss[y]*gauss[tx] * sqrtf(dx*dx + dy*dy);
    float angf = 4.0f/3.1415f*atan2f(dy, dx) + 4.0f;
    
    int hori = (tx + 2)/4 - 1;      // Convert from (tx,y,angle) to bins      
    float horf = (tx - 1.5f)/4.0f - hori;
    float ihorf = 1.0f - horf;           
    int veri = (y + 2)/4 - 1;
    float verf = (y - 1.5f)/4.0f - veri;
    float iverf = 1.0f - verf;
    int angi = angf;
    int angp = (angi<7 ? angi+1 : 0);
    angf -= angi;
    float iangf = 1.0f - angf;
    
    int hist = 8*(4*veri + hori);   // Each gradient measure is interpolated 
    int p1 = angi + hist;           // in angles, xpos and ypos -> 8 stores
    int p2 = angp + hist;
    if (tx>=2) { 
      float grad1 = ihorf*grad;
      if (y>=2) {   // Upper left
	float grad2 = iverf*grad1;
	atomicAdd(buffer + p1, iangf*grad2);
	atomicAdd(buffer + p2,  angf*grad2);
      }
      if (y<=13) {  // Lower left
	float grad2 = verf*grad1;
	atomicAdd(buffer + p1+32, iangf*grad2); 
	atomicAdd(buffer + p2+32,  angf*grad2);
      }
    }
    if (tx<=13) { 
      float grad1 = horf*grad;
      if (y>=2) {    // Upper right
	float grad2 = iverf*grad1;
	atomicAdd(buffer + p1+8, iangf*grad2);
	atomicAdd(buffer + p2+8,  angf*grad2);
      }
      if (y<=13) {   // Lower right
	float grad2 = verf*grad1;
	atomicAdd(buffer + p1+40, iangf*grad2);
	atomicAdd(buffer + p2+40,  angf*grad2);
      }
    }
  }
  __syncthreads();
    
  // Normalize twice and suppress peaks first time
  float sum = buffer[idx]*buffer[idx];
  for (int i=16;i>0;i/=2)
    sum += ShiftDown(sum, i);
  if ((idx&31)==0)
    sums[idx/32] = sum;
  __syncthreads();
  float tsum1 = sums[0] + sums[1] + sums[2] + sums[3]; 
  tsum1 = min(buffer[idx] * rsqrtf(tsum1), 0.2f);
  
  sum = tsum1*tsum1; 
  for (int i=16;i>0;i/=2)
    sum += ShiftDown(sum, i);
  if ((idx&31)==0)
    sums[idx/32] = sum;
  __syncthreads();
  
  float tsum2 = sums[0] + sums[1] + sums[2] + sums[3];
  float *desc = d_sift[bx].data;
  desc[idx] = tsum1 * rsqrtf(tsum2);
  if (idx==0) {
    d_sift[bx].xpos *= subsampling;
    d_sift[bx].ypos *= subsampling;
    d_sift[bx].scale *= subsampling;
  }
  __syncthreads();
}


__global__ void RescalePositions(SiftPoint *d_sift, int numPts, float scale)
{
  int num = blockIdx.x*blockDim.x + threadIdx.x;
  if (num<numPts) {
    d_sift[num].xpos *= scale;
    d_sift[num].ypos *= scale;
    d_sift[num].scale *= scale;
  }
}


__global__ void ComputeOrientations(cudaTextureObject_t texObj, SiftPoint *d_Sift, int fstPts)
{
  __shared__ float hist[64];
  __shared__ float gauss[11];
  const int tx = threadIdx.x;
  const int bx = blockIdx.x + fstPts;
  float i2sigma2 = -1.0f/(4.5f*d_Sift[bx].scale*d_Sift[bx].scale);
  if (tx<11) 
    gauss[tx] = exp(i2sigma2*(tx-5)*(tx-5));
  if (tx<64)
    hist[tx] = 0.0f;
  __syncthreads();
  float xp = d_Sift[bx].xpos - 4.5f;
  float yp = d_Sift[bx].ypos - 4.5f;
  int yd = tx/11;
  int xd = tx - yd*11;
  float xf = xp + xd;
  float yf = yp + yd;
  if (yd<11) {
    float dx = tex2D<float>(texObj, xf+1.0, yf) - tex2D<float>(texObj, xf-1.0, yf); 
    float dy = tex2D<float>(texObj, xf, yf+1.0) - tex2D<float>(texObj, xf, yf-1.0); 
    int bin = 16.0f*atan2f(dy, dx)/3.1416f + 16.5f;
    if (bin>31)
      bin = 0;
    float grad = sqrtf(dx*dx + dy*dy);
    atomicAdd(&hist[bin], grad*gauss[xd]*gauss[yd]);
  }
  __syncthreads();
  int x1m = (tx>=1 ? tx-1 : tx+31);
  int x1p = (tx<=30 ? tx+1 : tx-31);
  if (tx<32) {
    int x2m = (tx>=2 ? tx-2 : tx+30);
    int x2p = (tx<=29 ? tx+2 : tx-30);
    hist[tx+32] = 6.0f*hist[tx] + 4.0f*(hist[x1m] + hist[x1p]) + (hist[x2m] + hist[x2p]);
  }
  __syncthreads();
  if (tx<32) {
    float v = hist[32+tx];
    hist[tx] = (v>hist[32+x1m] && v>=hist[32+x1p] ? v : 0.0f);
  }
  __syncthreads();
  if (tx==0) {
    float maxval1 = 0.0;
    float maxval2 = 0.0;
    int i1 = -1;
    int i2 = -1;
    for (int i=0;i<32;i++) {
      float v = hist[i];
      if (v>maxval1) {
	maxval2 = maxval1;
	maxval1 = v;
	i2 = i1;
	i1 = i;
      } else if (v>maxval2) {
	maxval2 = v;
	i2 = i;
      }
    }
    float val1 = hist[32+((i1+1)&31)];
    float val2 = hist[32+((i1+31)&31)];
    float peak = i1 + 0.5f*(val1-val2) / (2.0f*maxval1-val1-val2);
    d_Sift[bx].orientation = 11.25f*(peak<0.0f ? peak+32.0f : peak);
    if (maxval2>0.8f*maxval1) {
      float val1 = hist[32+((i2+1)&31)];
      float val2 = hist[32+((i2+31)&31)];
      float peak = i2 + 0.5f*(val1-val2) / (2.0f*maxval2-val1-val2);
      unsigned int idx = atomicInc(d_PointCounter, 0x7fffffff);
      if (idx<d_MaxNumPoints) {
	d_Sift[idx].xpos = d_Sift[bx].xpos;
	d_Sift[idx].ypos = d_Sift[bx].ypos;
	d_Sift[idx].scale = d_Sift[bx].scale;
	d_Sift[idx].sharpness = d_Sift[bx].sharpness;
	d_Sift[idx].edgeness = d_Sift[bx].edgeness;
	d_Sift[idx].orientation = 11.25f*(peak<0.0f ? peak+32.0f : peak);;
	d_Sift[idx].subsampling = d_Sift[bx].subsampling;
      }
    } 
  }
} 

// With constant number of blocks
__global__ void ComputeOrientationsCONST(cudaTextureObject_t texObj, SiftPoint *d_Sift, int octave)
{
  __shared__ float hist[64];
  __shared__ float gauss[11];
  const int tx = threadIdx.x;
  
  int fstPts = min(d_PointCounter[2*octave-1], d_MaxNumPoints);
  int totPts = min(d_PointCounter[2*octave+0], d_MaxNumPoints);  
  for (int bx = blockIdx.x + fstPts; bx < totPts; bx += gridDim.x) {
 
    float i2sigma2 = -1.0f/(4.5f*d_Sift[bx].scale*d_Sift[bx].scale);
    if (tx<11) 
      gauss[tx] = exp(i2sigma2*(tx-5)*(tx-5));
    if (tx<64)
      hist[tx] = 0.0f;
    __syncthreads();
    float xp = d_Sift[bx].xpos - 4.5f;
    float yp = d_Sift[bx].ypos - 4.5f;
    int yd = tx/11;
    int xd = tx - yd*11;
    float xf = xp + xd;
    float yf = yp + yd;
    if (yd<11) {
      float dx = tex2D<float>(texObj, xf+1.0, yf) - tex2D<float>(texObj, xf-1.0, yf); 
      float dy = tex2D<float>(texObj, xf, yf+1.0) - tex2D<float>(texObj, xf, yf-1.0); 
      int bin = 16.0f*atan2f(dy, dx)/3.1416f + 16.5f;
      if (bin>31)
	bin = 0;
      float grad = sqrtf(dx*dx + dy*dy);
      atomicAdd(&hist[bin], grad*gauss[xd]*gauss[yd]);
    }
    __syncthreads();
    int x1m = (tx>=1 ? tx-1 : tx+31);
    int x1p = (tx<=30 ? tx+1 : tx-31);
    if (tx<32) {
      int x2m = (tx>=2 ? tx-2 : tx+30);
      int x2p = (tx<=29 ? tx+2 : tx-30);
      hist[tx+32] = 6.0f*hist[tx] + 4.0f*(hist[x1m] + hist[x1p]) + (hist[x2m] + hist[x2p]);
    }
    __syncthreads();
    if (tx<32) {
      float v = hist[32+tx];
      hist[tx] = (v>hist[32+x1m] && v>=hist[32+x1p] ? v : 0.0f);
    }
    __syncthreads();
    if (tx==0) {
      float maxval1 = 0.0;
      float maxval2 = 0.0;
      int i1 = -1;
      int i2 = -1;
      for (int i=0;i<32;i++) {
	float v = hist[i];
	if (v>maxval1) {
	  maxval2 = maxval1;
	  maxval1 = v;
	  i2 = i1;
	  i1 = i;
	} else if (v>maxval2) {
	  maxval2 = v;
	  i2 = i;
	}
      }
      float val1 = hist[32+((i1+1)&31)];
      float val2 = hist[32+((i1+31)&31)];
      float peak = i1 + 0.5f*(val1-val2) / (2.0f*maxval1-val1-val2);
      d_Sift[bx].orientation = 11.25f*(peak<0.0f ? peak+32.0f : peak);
      atomicMax(&d_PointCounter[2*octave+1], d_PointCounter[2*octave+0]); 
      if (maxval2>0.8f*maxval1) {
	float val1 = hist[32+((i2+1)&31)];
	float val2 = hist[32+((i2+31)&31)];
	float peak = i2 + 0.5f*(val1-val2) / (2.0f*maxval2-val1-val2);
	unsigned int idx = atomicInc(&d_PointCounter[2*octave+1], 0x7fffffff);
	if (idx<d_MaxNumPoints) {
	  d_Sift[idx].xpos = d_Sift[bx].xpos;
	  d_Sift[idx].ypos = d_Sift[bx].ypos;
	  d_Sift[idx].scale = d_Sift[bx].scale;
	  d_Sift[idx].sharpness = d_Sift[bx].sharpness;
	  d_Sift[idx].edgeness = d_Sift[bx].edgeness;
	  d_Sift[idx].orientation = 11.25f*(peak<0.0f ? peak+32.0f : peak);;
	  d_Sift[idx].subsampling = d_Sift[bx].subsampling;
	}
      }
    }
    __syncthreads();
  }
} 

// With constant number of blocks
__global__ void OrientAndExtractCONST(cudaTextureObject_t texObj, SiftPoint *d_Sift, float subsampling, int octave)
{
  __shared__ float hist[64];
  __shared__ float gauss[11];
  __shared__ unsigned int idx; //%%%%
  const int tx = threadIdx.x;
  
  int fstPts = min(d_PointCounter[2*octave-1], d_MaxNumPoints);
  int totPts = min(d_PointCounter[2*octave+0], d_MaxNumPoints);  
  for (int bx = blockIdx.x + fstPts; bx < totPts; bx += gridDim.x) {
 
    float i2sigma2 = -1.0f/(4.5f*d_Sift[bx].scale*d_Sift[bx].scale);
    if (tx<11) 
      gauss[tx] = exp(i2sigma2*(tx-5)*(tx-5));
    if (tx<64)
      hist[tx] = 0.0f;
    __syncthreads();
    float xp = d_Sift[bx].xpos - 4.5f;
    float yp = d_Sift[bx].ypos - 4.5f;
    int yd = tx/11;
    int xd = tx - yd*11;
    float xf = xp + xd;
    float yf = yp + yd;
    if (yd<11) {
      float dx = tex2D<float>(texObj, xf+1.0, yf) - tex2D<float>(texObj, xf-1.0, yf); 
      float dy = tex2D<float>(texObj, xf, yf+1.0) - tex2D<float>(texObj, xf, yf-1.0); 
      int bin = 16.0f*atan2f(dy, dx)/3.1416f + 16.5f;
      if (bin>31)
	bin = 0;
      float grad = sqrtf(dx*dx + dy*dy);
      atomicAdd(&hist[bin], grad*gauss[xd]*gauss[yd]);
    }
    __syncthreads();
    int x1m = (tx>=1 ? tx-1 : tx+31);
    int x1p = (tx<=30 ? tx+1 : tx-31);
    if (tx<32) {
      int x2m = (tx>=2 ? tx-2 : tx+30);
      int x2p = (tx<=29 ? tx+2 : tx-30);
      hist[tx+32] = 6.0f*hist[tx] + 4.0f*(hist[x1m] + hist[x1p]) + (hist[x2m] + hist[x2p]);
    }
    __syncthreads();
    if (tx<32) {
      float v = hist[32+tx];
      hist[tx] = (v>hist[32+x1m] && v>=hist[32+x1p] ? v : 0.0f);
    }
    __syncthreads();
    if (tx==0) {
      float maxval1 = 0.0;
      float maxval2 = 0.0;
      int i1 = -1;
      int i2 = -1;
      for (int i=0;i<32;i++) {
	float v = hist[i];
	if (v>maxval1) {
	  maxval2 = maxval1;
	  maxval1 = v;
	  i2 = i1;
	  i1 = i;
	} else if (v>maxval2) {
	  maxval2 = v;
	  i2 = i;
	}
      }
      float val1 = hist[32+((i1+1)&31)];
      float val2 = hist[32+((i1+31)&31)];
      float peak = i1 + 0.5f*(val1-val2) / (2.0f*maxval1-val1-val2);
      d_Sift[bx].orientation = 11.25f*(peak<0.0f ? peak+32.0f : peak);
      idx = 0xffffffff; //%%%%
      atomicMax(&d_PointCounter[2*octave+1], d_PointCounter[2*octave+0]); 
      if (maxval2>0.8f*maxval1) {
	float val1 = hist[32+((i2+1)&31)];
	float val2 = hist[32+((i2+31)&31)];
	float peak = i2 + 0.5f*(val1-val2) / (2.0f*maxval2-val1-val2);
	idx = atomicInc(&d_PointCounter[2*octave+1], 0x7fffffff); //%%%%
	if (idx<d_MaxNumPoints) {
	  d_Sift[idx].xpos = d_Sift[bx].xpos;
	  d_Sift[idx].ypos = d_Sift[bx].ypos;
	  d_Sift[idx].scale = d_Sift[bx].scale;
	  d_Sift[idx].sharpness = d_Sift[bx].sharpness;
	  d_Sift[idx].edgeness = d_Sift[bx].edgeness;
	  d_Sift[idx].orientation = 11.25f*(peak<0.0f ? peak+32.0f : peak);;
	  d_Sift[idx].subsampling = d_Sift[bx].subsampling;
	}
      }
    }
    __syncthreads();
    ExtractSiftDescriptor(texObj, d_Sift, subsampling, octave, bx); //%%%%
    if (idx<d_MaxNumPoints) //%%%%
      ExtractSiftDescriptor(texObj, d_Sift, subsampling, octave, idx); //%%%%
  }
} 


///////////////////////////////////////////////////////////////////////////////
// Subtract two images (multi-scale version)
///////////////////////////////////////////////////////////////////////////////
  
__global__ void FindPointsMultiTest(float *d_Data0, SiftPoint *d_Sift, int width, int pitch, int height, float subsampling, float lowestScale, float thresh, float factor, float edgeLimit, int octave)
{
  #define MEMWID (MINMAX_W + 2)
  __shared__ unsigned int cnt;
  __shared__ unsigned short points[3*MEMWID];

  if (blockIdx.x==0 && blockIdx.y==0 && threadIdx.x==0 && threadIdx.y==0) {
    atomicMax(&d_PointCounter[2*octave+0], d_PointCounter[2*octave-1]); 
    atomicMax(&d_PointCounter[2*octave+1], d_PointCounter[2*octave-1]);
  }
  int tx = threadIdx.x;
  int ty = threadIdx.y;
  if (tx==0 && ty==0)
    cnt = 0; 
  __syncthreads();

  int ypos = MINMAX_H*blockIdx.y + ty;
  if (ypos>=height)
    return;
  int block = blockIdx.x/NUM_SCALES; 
  int scale = blockIdx.x - NUM_SCALES*block;
  int minx = block*MINMAX_W;
  int maxx = min(minx + MINMAX_W, width);
  int xpos = minx + tx;
  int size = pitch*height;
  int ptr = size*scale + max(min(xpos-1, width-1), 0);

  float maxv = fabs(d_Data0[ptr + ypos*pitch + 1*size]);
  maxv = fmaxf(maxv, ShiftDown(maxv, 16, MINMAX_W));
  maxv = fmaxf(maxv, ShiftDown(maxv, 8, MINMAX_W));
  maxv = fmaxf(maxv, ShiftDown(maxv, 4, MINMAX_W));
  maxv = fmaxf(maxv, ShiftDown(maxv, 2, MINMAX_W));
  maxv = fmaxf(maxv, ShiftDown(maxv, 1, MINMAX_W));
    
  if (Shuffle(maxv, 0)>thresh) {
    int yptr1 = ptr + ypos*pitch;
    int yptr0 = ptr + max(0,ypos-1)*pitch;
    int yptr2 = ptr + min(height-1,ypos+1)*pitch;
    float d20 = d_Data0[yptr0 + 1*size];
    float d21 = d_Data0[yptr1 + 1*size];
    float d22 = d_Data0[yptr2 + 1*size];
    float d31 = d_Data0[yptr1 + 2*size];
    float d11 = d_Data0[yptr1];
    
    float d10 = d_Data0[yptr0];
    float d12 = d_Data0[yptr2];
    float ymin1 = fminf(fminf(d10, d11), d12);
    float ymax1 = fmaxf(fmaxf(d10, d11), d12);
    float d30 = d_Data0[yptr0 + 2*size];
    float d32 = d_Data0[yptr2 + 2*size]; 
    float ymin3 = fminf(fminf(d30, d31), d32);
    float ymax3 = fmaxf(fmaxf(d30, d31), d32);
    float ymin2 = fminf(fminf(ymin1, fminf(fminf(d20, d22), d21)), ymin3);
    float ymax2 = fmaxf(fmaxf(ymax1, fmaxf(fmaxf(d20, d22), d21)), ymax3);
    
    float nmin2 = fminf(ShiftUp(ymin2, 1), ShiftDown(ymin2, 1));
    float nmax2 = fmaxf(ShiftUp(ymax2, 1), ShiftDown(ymax2, 1));
    if (tx>0 && tx<MINMAX_W+1 && xpos<=maxx) {
      if (d21<-thresh) {
	float minv = fminf(fminf(nmin2, ymin1), ymin3);
	minv = fminf(fminf(minv, d20), d22);
	if (d21<minv) { 
	  int pos = atomicInc(&cnt, MEMWID-1);
	  points[3*pos+0] = xpos - 1;
	  points[3*pos+1] = ypos;
	  points[3*pos+2] = scale;
	}
      } 
      if (d21>thresh) {
	float maxv = fmaxf(fmaxf(nmax2, ymax1), ymax3);
	maxv = fmaxf(fmaxf(maxv, d20), d22);
	if (d21>maxv) { 
	  int pos = atomicInc(&cnt, MEMWID-1);
	  points[3*pos+0] = xpos - 1;
	  points[3*pos+1] = ypos;
	  points[3*pos+2] = scale;
	}
      }
    }
  }
  __syncthreads();
  if (ty==0 && tx<cnt) {
    int xpos = points[3*tx+0];
    int ypos = points[3*tx+1];
    int scale = points[3*tx+2];
    int ptr = xpos + (ypos + (scale+1)*height)*pitch;
    float val = d_Data0[ptr];
    float *data1 = &d_Data0[ptr];
    float dxx = 2.0f*val - data1[-1] - data1[1];
    float dyy = 2.0f*val - data1[-pitch] - data1[pitch];
    float dxy = 0.25f*(data1[+pitch+1] + data1[-pitch-1] - data1[-pitch+1] - data1[+pitch-1]);
    float tra = dxx + dyy;
    float det = dxx*dyy - dxy*dxy;
    if (tra*tra<edgeLimit*det) {
      float edge = __fdividef(tra*tra, det);
      float dx = 0.5f*(data1[1] - data1[-1]);
      float dy = 0.5f*(data1[pitch] - data1[-pitch]); 
      float *data0 = d_Data0 + ptr - height*pitch;
      float *data2 = d_Data0 + ptr + height*pitch;
      float ds = 0.5f*(data0[0] - data2[0]); 
      float dss = 2.0f*val - data2[0] - data0[0];
      float dxs = 0.25f*(data2[1] + data0[-1] - data0[1] - data2[-1]);
      float dys = 0.25f*(data2[pitch] + data0[-pitch] - data2[-pitch] - data0[pitch]);
      float idxx = dyy*dss - dys*dys;
      float idxy = dys*dxs - dxy*dss;   
      float idxs = dxy*dys - dyy*dxs;
      float idet = __fdividef(1.0f, idxx*dxx + idxy*dxy + idxs*dxs);
      float idyy = dxx*dss - dxs*dxs;
      float idys = dxy*dxs - dxx*dys;
      float idss = dxx*dyy - dxy*dxy;
      float pdx = idet*(idxx*dx + idxy*dy + idxs*ds);
      float pdy = idet*(idxy*dx + idyy*dy + idys*ds);
      float pds = idet*(idxs*dx + idys*dy + idss*ds);
      if (pdx<-0.5f || pdx>0.5f || pdy<-0.5f || pdy>0.5f || pds<-0.5f || pds>0.5f) {
	pdx = __fdividef(dx, dxx);
	pdy = __fdividef(dy, dyy);
	pds = __fdividef(ds, dss);
      }
      float dval = 0.5f*(dx*pdx + dy*pdy + ds*pds);
      int maxPts = d_MaxNumPoints;
      float sc = powf(2.0f, (float)scale/NUM_SCALES) * exp2f(pds*factor);
      if (sc>=lowestScale) {
	unsigned int idx = atomicInc(&d_PointCounter[2*octave+0], 0x7fffffff);
	idx = (idx>=maxPts ? maxPts-1 : idx);
	d_Sift[idx].xpos = xpos + pdx;
	d_Sift[idx].ypos = ypos + pdy;
	d_Sift[idx].scale = sc;
	d_Sift[idx].sharpness = val + dval;
	d_Sift[idx].edgeness = edge;
	d_Sift[idx].subsampling = subsampling;
      }
    }
  }
}

__global__ void FindPointsMulti(float *d_Data0, SiftPoint *d_Sift, int width, int pitch, int height, float subsampling, float lowestScale, float thresh, float factor, float edgeLimit, int octave)
{
  #define MEMWID (MINMAX_W + 2)
  __shared__ unsigned int cnt;
  __shared__ unsigned short points[3*MEMWID];

  if (blockIdx.x==0 && blockIdx.y==0 && threadIdx.x==0) {
    atomicMax(&d_PointCounter[2*octave+0], d_PointCounter[2*octave-1]);
    atomicMax(&d_PointCounter[2*octave+1], d_PointCounter[2*octave-1]);
  }
  int tx = threadIdx.x;
  int block = blockIdx.x/NUM_SCALES; 
  int scale = blockIdx.x - NUM_SCALES*block;
  int minx = block*MINMAX_W;
  int maxx = min(minx + MINMAX_W, width);
  int xpos = minx + tx;
  int size = pitch*height;
  int ptr = size*scale + max(min(xpos-1, width-1), 0);

  int yloops = min(height - MINMAX_H*blockIdx.y, MINMAX_H);
  float maxv = 0.0f;
  for (int y=0;y<yloops;y++) {
    int ypos = MINMAX_H*blockIdx.y + y;
    int yptr1 = ptr + ypos*pitch;
    float val = d_Data0[yptr1 + 1*size];
    maxv = fmaxf(maxv, fabs(val));
  }
  maxv = fmaxf(maxv, ShiftDown(maxv, 16, MINMAX_W));
  maxv = fmaxf(maxv, ShiftDown(maxv, 8, MINMAX_W));
  maxv = fmaxf(maxv, ShiftDown(maxv, 4, MINMAX_W));
  maxv = fmaxf(maxv, ShiftDown(maxv, 2, MINMAX_W));
  maxv = fmaxf(maxv, ShiftDown(maxv, 1, MINMAX_W));
  if (Shuffle(maxv, 0)<=thresh)
    return;
  
  if (tx==0)
    cnt = 0; 
  __syncthreads();

  for (int y=0;y<yloops;y++) {

    int ypos = MINMAX_H*blockIdx.y + y;
    int yptr1 = ptr + ypos*pitch;
    int yptr0 = ptr + max(0,ypos-1)*pitch;
    int yptr2 = ptr + min(height-1,ypos+1)*pitch;
    float d20 = d_Data0[yptr0 + 1*size];
    float d21 = d_Data0[yptr1 + 1*size];
    float d22 = d_Data0[yptr2 + 1*size];
    float d31 = d_Data0[yptr1 + 2*size];
    float d11 = d_Data0[yptr1];
    
    float d10 = d_Data0[yptr0];
    float d12 = d_Data0[yptr2];
    float ymin1 = fminf(fminf(d10, d11), d12);
    float ymax1 = fmaxf(fmaxf(d10, d11), d12);
    float d30 = d_Data0[yptr0 + 2*size];
    float d32 = d_Data0[yptr2 + 2*size]; 
    float ymin3 = fminf(fminf(d30, d31), d32);
    float ymax3 = fmaxf(fmaxf(d30, d31), d32);
    float ymin2 = fminf(fminf(ymin1, fminf(fminf(d20, d22), d21)), ymin3);
    float ymax2 = fmaxf(fmaxf(ymax1, fmaxf(fmaxf(d20, d22), d21)), ymax3);
    
    float nmin2 = fminf(ShiftUp(ymin2, 1), ShiftDown(ymin2, 1));
    float nmax2 = fmaxf(ShiftUp(ymax2, 1), ShiftDown(ymax2, 1));
    if (tx>0 && tx<MINMAX_W+1 && xpos<=maxx) {
      if (d21<-thresh) {
	float minv = fminf(fminf(nmin2, ymin1), ymin3);
	minv = fminf(fminf(minv, d20), d22);
	if (d21<minv) { 
	  int pos = atomicInc(&cnt, MEMWID-1);
	  points[3*pos+0] = xpos - 1;
	  points[3*pos+1] = ypos;
	  points[3*pos+2] = scale;
	}
      } 
      if (d21>thresh) {
	float maxv = fmaxf(fmaxf(nmax2, ymax1), ymax3);
	maxv = fmaxf(fmaxf(maxv, d20), d22);
	if (d21>maxv) { 
	  int pos = atomicInc(&cnt, MEMWID-1);
	  points[3*pos+0] = xpos - 1;
	  points[3*pos+1] = ypos;
	  points[3*pos+2] = scale;
	}
      }
    }
  }
  if (tx<cnt) {
    int xpos = points[3*tx+0];
    int ypos = points[3*tx+1];
    int scale = points[3*tx+2];
    int ptr = xpos + (ypos + (scale+1)*height)*pitch;
    float val = d_Data0[ptr];
    float *data1 = &d_Data0[ptr];
    float dxx = 2.0f*val - data1[-1] - data1[1];
    float dyy = 2.0f*val - data1[-pitch] - data1[pitch];
    float dxy = 0.25f*(data1[+pitch+1] + data1[-pitch-1] - data1[-pitch+1] - data1[+pitch-1]);
    float tra = dxx + dyy;
    float det = dxx*dyy - dxy*dxy;
    if (tra*tra<edgeLimit*det) {
      float edge = __fdividef(tra*tra, det);
      float dx = 0.5f*(data1[1] - data1[-1]);
      float dy = 0.5f*(data1[pitch] - data1[-pitch]); 
      float *data0 = d_Data0 + ptr - height*pitch;
      float *data2 = d_Data0 + ptr + height*pitch;
      float ds = 0.5f*(data0[0] - data2[0]); 
      float dss = 2.0f*val - data2[0] - data0[0];
      float dxs = 0.25f*(data2[1] + data0[-1] - data0[1] - data2[-1]);
      float dys = 0.25f*(data2[pitch] + data0[-pitch] - data2[-pitch] - data0[pitch]);
      float idxx = dyy*dss - dys*dys;
      float idxy = dys*dxs - dxy*dss;   
      float idxs = dxy*dys - dyy*dxs;
      float idet = __fdividef(1.0f, idxx*dxx + idxy*dxy + idxs*dxs);
      float idyy = dxx*dss - dxs*dxs;
      float idys = dxy*dxs - dxx*dys;
      float idss = dxx*dyy - dxy*dxy;
      float pdx = idet*(idxx*dx + idxy*dy + idxs*ds);
      float pdy = idet*(idxy*dx + idyy*dy + idys*ds);
      float pds = idet*(idxs*dx + idys*dy + idss*ds);
      if (pdx<-0.5f || pdx>0.5f || pdy<-0.5f || pdy>0.5f || pds<-0.5f || pds>0.5f) {
	pdx = __fdividef(dx, dxx);
	pdy = __fdividef(dy, dyy);
	pds = __fdividef(ds, dss);
      }
      float dval = 0.5f*(dx*pdx + dy*pdy + ds*pds);
      int maxPts = d_MaxNumPoints;
      float sc = powf(2.0f, (float)scale/NUM_SCALES) * exp2f(pds*factor);
      if (sc>=lowestScale) {
	atomicMax(&d_PointCounter[2*octave+0], d_PointCounter[2*octave-1]); 
	unsigned int idx = atomicInc(&d_PointCounter[2*octave+0], 0x7fffffff);
	idx = (idx>=maxPts ? maxPts-1 : idx);
	d_Sift[idx].xpos = xpos + pdx;
	d_Sift[idx].ypos = ypos + pdy;
	d_Sift[idx].scale = sc;
	d_Sift[idx].sharpness = val + dval;
	d_Sift[idx].edgeness = edge;
	d_Sift[idx].subsampling = subsampling;
      }
    }
  }
}


__global__ void FindPointsMultiOld(float *d_Data0, SiftPoint *d_Sift, int width, int pitch, int height, float subsampling, float lowestScale, float thresh, float factor, float edgeLimit, int octave)
{
  #define MEMWID (MINMAX_W + 2)
  __shared__ float ymin1[MEMWID], ymin2[MEMWID], ymin3[MEMWID];
  __shared__ float ymax1[MEMWID], ymax2[MEMWID], ymax3[MEMWID];
  __shared__ unsigned int cnt;
  __shared__ unsigned short points[3*MEMWID];

  if (blockIdx.x==0 && blockIdx.y==0 && threadIdx.x==0) {
    atomicMax(&d_PointCounter[2*octave+0], d_PointCounter[2*octave-1]); 
    atomicMax(&d_PointCounter[2*octave+1], d_PointCounter[2*octave-1]);
  }
  int tx = threadIdx.x;
  int block = blockIdx.x/NUM_SCALES; 
  int scale = blockIdx.x - NUM_SCALES*block;
  int minx = block*MINMAX_W;
  int maxx = min(minx + MINMAX_W, width);
  int xpos = minx + tx;
  int size = pitch*height;
  int ptr = size*scale + max(min(xpos-1, width-1), 0);

  int yloops = min(height - MINMAX_H*blockIdx.y, MINMAX_H);
  float maxv = 0.0f;
  for (int y=0;y<yloops;y++) {
    int ypos = MINMAX_H*blockIdx.y + y;
    int yptr1 = ptr + ypos*pitch;
    float val = d_Data0[yptr1 + 1*size];
    maxv = fmaxf(maxv, fabs(val));
  }
  maxv = fmaxf(maxv, ShiftDown(maxv, 16, MINMAX_W));
  maxv = fmaxf(maxv, ShiftDown(maxv, 8, MINMAX_W));
  maxv = fmaxf(maxv, ShiftDown(maxv, 4, MINMAX_W));
  maxv = fmaxf(maxv, ShiftDown(maxv, 2, MINMAX_W));
  maxv = fmaxf(maxv, ShiftDown(maxv, 1, MINMAX_W));
  if (Shuffle(maxv, 0)<=thresh)
    return;
  
  if (tx==0)
    cnt = 0; 
  __syncthreads();

  for (int y=0;y<yloops;y++) {

    int ypos = MINMAX_H*blockIdx.y + y;
    int yptr1 = ptr + ypos*pitch;
    int yptr0 = ptr + max(0,ypos-1)*pitch;
    int yptr2 = ptr + min(height-1,ypos+1)*pitch;
    float d20 = d_Data0[yptr0 + 1*size];
    float d21 = d_Data0[yptr1 + 1*size];
    float d22 = d_Data0[yptr2 + 1*size];
    float d31 = d_Data0[yptr1 + 2*size];
    float d11 = d_Data0[yptr1];

    float d10 = d_Data0[yptr0];
    float d12 = d_Data0[yptr2];
    ymin1[tx] = fminf(fminf(d10, d11), d12);
    ymax1[tx] = fmaxf(fmaxf(d10, d11), d12);
    float d30 = d_Data0[yptr0 + 2*size];
    float d32 = d_Data0[yptr2 + 2*size]; 
    ymin3[tx] = fminf(fminf(d30, d31), d32);
    ymax3[tx] = fmaxf(fmaxf(d30, d31), d32);
    ymin2[tx] = fminf(fminf(ymin1[tx], fminf(fminf(d20, d22), d21)), ymin3[tx]);
    ymax2[tx] = fmaxf(fmaxf(ymax1[tx], fmaxf(fmaxf(d20, d22), d21)), ymax3[tx]);
    
    __syncthreads(); 

    if (tx>0 && tx<MINMAX_W+1 && xpos<=maxx) {
      if (d21<-thresh) {
	float minv = fminf(fminf(fminf(ymin2[tx-1], ymin2[tx+1]), ymin1[tx]), ymin3[tx]);
	minv = fminf(fminf(minv, d20), d22);
	if (d21<minv) { 
	  int pos = atomicInc(&cnt, MEMWID-1);
	  points[3*pos+0] = xpos - 1;
	  points[3*pos+1] = ypos;
	  points[3*pos+2] = scale;
	}
      } 
      if (d21>thresh) {
	float maxv = fmaxf(fmaxf(fmaxf(ymax2[tx-1], ymax2[tx+1]), ymax1[tx]), ymax3[tx]);
	maxv = fmaxf(fmaxf(maxv, d20), d22);
	if (d21>maxv) { 
	  int pos = atomicInc(&cnt, MEMWID-1);
	  points[3*pos+0] = xpos - 1;
	  points[3*pos+1] = ypos;
	  points[3*pos+2] = scale;
	}
      }
    }
    __syncthreads();
  }
  if (tx<cnt) {
    int xpos = points[3*tx+0];
    int ypos = points[3*tx+1];
    int scale = points[3*tx+2];
    int ptr = xpos + (ypos + (scale+1)*height)*pitch;
    float val = d_Data0[ptr];
    float *data1 = &d_Data0[ptr];
    float dxx = 2.0f*val - data1[-1] - data1[1];
    float dyy = 2.0f*val - data1[-pitch] - data1[pitch];
    float dxy = 0.25f*(data1[+pitch+1] + data1[-pitch-1] - data1[-pitch+1] - data1[+pitch-1]);
    float tra = dxx + dyy;
    float det = dxx*dyy - dxy*dxy;
    if (tra*tra<edgeLimit*det) {
      float edge = __fdividef(tra*tra, det);
      float dx = 0.5f*(data1[1] - data1[-1]);
      float dy = 0.5f*(data1[pitch] - data1[-pitch]); 
      float *data0 = d_Data0 + ptr - height*pitch;
      float *data2 = d_Data0 + ptr + height*pitch;
      float ds = 0.5f*(data0[0] - data2[0]); 
      float dss = 2.0f*val - data2[0] - data0[0];
      float dxs = 0.25f*(data2[1] + data0[-1] - data0[1] - data2[-1]);
      float dys = 0.25f*(data2[pitch] + data0[-pitch] - data2[-pitch] - data0[pitch]);
      float idxx = dyy*dss - dys*dys;
      float idxy = dys*dxs - dxy*dss;   
      float idxs = dxy*dys - dyy*dxs;
      float idet = __fdividef(1.0f, idxx*dxx + idxy*dxy + idxs*dxs);
      float idyy = dxx*dss - dxs*dxs;
      float idys = dxy*dxs - dxx*dys;
      float idss = dxx*dyy - dxy*dxy;
      float pdx = idet*(idxx*dx + idxy*dy + idxs*ds);
      float pdy = idet*(idxy*dx + idyy*dy + idys*ds);
      float pds = idet*(idxs*dx + idys*dy + idss*ds);
      if (pdx<-0.5f || pdx>0.5f || pdy<-0.5f || pdy>0.5f || pds<-0.5f || pds>0.5f) {
	pdx = __fdividef(dx, dxx);
	pdy = __fdividef(dy, dyy);
	pds = __fdividef(ds, dss);
      }
      float dval = 0.5f*(dx*pdx + dy*pdy + ds*pds);
      int maxPts = d_MaxNumPoints;
      float sc = powf(2.0f, (float)scale/NUM_SCALES) * exp2f(pds*factor);
      if (sc>=lowestScale) {
	unsigned int idx = atomicInc(&d_PointCounter[2*octave+0], 0x7fffffff);
	idx = (idx>=maxPts ? maxPts-1 : idx);
	d_Sift[idx].xpos = xpos + pdx;
	d_Sift[idx].ypos = ypos + pdy;
	d_Sift[idx].scale = sc;
	d_Sift[idx].sharpness = val + dval;
	d_Sift[idx].edgeness = edge;
	d_Sift[idx].subsampling = subsampling;
      }
    }
  }
}


__global__ void LaplaceMultiTex(cudaTextureObject_t texObj, float *d_Result, int width, int pitch, int height, int octave)
{
  __shared__ float data1[(LAPLACE_W + 2*LAPLACE_R)*LAPLACE_S];
  __shared__ float data2[LAPLACE_W*LAPLACE_S];
  const int tx = threadIdx.x;
  const int xp = blockIdx.x*LAPLACE_W + tx;
  const int yp = blockIdx.y;
  const int scale = threadIdx.y;
  float *kernel = d_LaplaceKernel + octave*12*16 + scale*16;
  float *sdata1 = data1 + (LAPLACE_W + 2*LAPLACE_R)*scale; 
  float x = xp-3.5;
  float y = yp+0.5;
  sdata1[tx] = kernel[4]*tex2D<float>(texObj, x, y) + 
    kernel[3]*(tex2D<float>(texObj, x, y-1.0) + tex2D<float>(texObj, x, y+1.0)) + 
    kernel[2]*(tex2D<float>(texObj, x, y-2.0) + tex2D<float>(texObj, x, y+2.0)) + 
    kernel[1]*(tex2D<float>(texObj, x, y-3.0) + tex2D<float>(texObj, x, y+3.0)) + 
    kernel[0]*(tex2D<float>(texObj, x, y-4.0) + tex2D<float>(texObj, x, y+4.0));
  __syncthreads();
  float *sdata2 = data2 + LAPLACE_W*scale; 
  if (tx<LAPLACE_W) {
    sdata2[tx] = kernel[4]*sdata1[tx+4] + 
      kernel[3]*(sdata1[tx+3] + sdata1[tx+5]) + 
      kernel[2]*(sdata1[tx+2] + sdata1[tx+6]) + 
      kernel[1]*(sdata1[tx+1] + sdata1[tx+7]) + 
      kernel[0]*(sdata1[tx+0] + sdata1[tx+8]);
  }
  __syncthreads(); 
  if (tx<LAPLACE_W && scale<LAPLACE_S-1 && xp<width) 
    d_Result[scale*height*pitch + yp*pitch + xp] = sdata2[tx] - sdata2[tx+LAPLACE_W];
}


__global__ void LaplaceMultiMem(float *d_Image, float *d_Result, int width, int pitch, int height, int octave)
{
  __shared__ float buff[(LAPLACE_W + 2*LAPLACE_R)*LAPLACE_S];
  const int tx = threadIdx.x;
  const int xp = blockIdx.x*LAPLACE_W + tx;
  const int yp = blockIdx.y;
  float *data = d_Image + max(min(xp - 4, width-1), 0);
  float temp[9], kern[LAPLACE_S][LAPLACE_R+1];
  if (xp<(width + 2*LAPLACE_R)) {
    for (int i=0;i<4;i++)
      temp[i] = data[max(0, min(yp+i-4, height-1))*pitch];
    for (int i=4;i<8+1;i++)
      temp[i] = data[min(yp+i-4, height-1)*pitch];
    for (int scale=0;scale<LAPLACE_S;scale++) {
      float *buf = buff + (LAPLACE_W + 2*LAPLACE_R)*scale; 
      float *kernel = d_LaplaceKernel + octave*12*16 + scale*16; 
      for (int i=0;i<LAPLACE_R+1;i++)
	kern[scale][i] = kernel[i];
      buf[tx] = kern[scale][4]*temp[4] +
	kern[scale][3]*(temp[3] + temp[5]) + kern[scale][2]*(temp[2] + temp[6]) + 
	kern[scale][1]*(temp[1] + temp[7]) + kern[scale][0]*(temp[0] + temp[8]);
    }
  }
  __syncthreads();
  if (tx<LAPLACE_W && xp<width) {
    float oldRes = kern[0][4]*buff[tx+4] + 
      kern[0][3]*(buff[tx+3] + buff[tx+5]) + kern[0][2]*(buff[tx+2] + buff[tx+6]) + 
      kern[0][1]*(buff[tx+1] + buff[tx+7]) + kern[0][0]*(buff[tx+0] + buff[tx+8]);
    for (int scale=1;scale<LAPLACE_S;scale++) {
      float *buf = buff + (LAPLACE_W + 2*LAPLACE_R)*scale; 
      float res = kern[scale][4]*buf[tx+4] + 
	kern[scale][3]*(buf[tx+3] + buf[tx+5]) + kern[scale][2]*(buf[tx+2] + buf[tx+6]) + 
	kern[scale][1]*(buf[tx+1] + buf[tx+7]) + kern[scale][0]*(buf[tx+0] + buf[tx+8]);
      d_Result[(scale-1)*height*pitch + yp*pitch + xp] = res - oldRes;
      oldRes = res;
    }
  }
}

__global__ void LaplaceMultiMemTest(float *d_Image, float *d_Result, int width, int pitch, int height, int octave)
{
  __shared__ float data1[(LAPLACE_W + 2*LAPLACE_R)*LAPLACE_S];
  __shared__ float data2[LAPLACE_W*LAPLACE_S];
  const int tx = threadIdx.x;
  const int xp = blockIdx.x*LAPLACE_W + tx;
  const int yp = LAPLACE_H*blockIdx.y;
  const int scale = threadIdx.y;
  float *kernel = d_LaplaceKernel + octave*12*16 + scale*16; 
  float *sdata1 = data1 + (LAPLACE_W + 2*LAPLACE_R)*scale; 
  float *data = d_Image + max(min(xp - 4, width-1), 0);
  int h = height-1;
  float temp[8+LAPLACE_H], kern[LAPLACE_R+1];
  for (int i=0;i<4;i++)
    temp[i] = data[max(0, min(yp+i-4, h))*pitch];
  for (int i=4;i<8+LAPLACE_H;i++)
    temp[i] = data[min(yp+i-4, h)*pitch];
  for (int i=0;i<LAPLACE_R+1;i++)
    kern[i] = kernel[i];
  for (int j=0;j<LAPLACE_H;j++) {
    sdata1[tx] = kern[4]*temp[4+j] +
      kern[3]*(temp[3+j] + temp[5+j]) + kern[2]*(temp[2+j] + temp[6+j]) + 
      kern[1]*(temp[1+j] + temp[7+j]) + kern[0]*(temp[0+j] + temp[8+j]);
    __syncthreads();
    float *sdata2 = data2 + LAPLACE_W*scale; 
    if (tx<LAPLACE_W) {
      sdata2[tx] = kern[4]*sdata1[tx+4] + 
	kern[3]*(sdata1[tx+3] + sdata1[tx+5]) + kern[2]*(sdata1[tx+2] + sdata1[tx+6]) + 
	kern[1]*(sdata1[tx+1] + sdata1[tx+7]) + kern[0]*(sdata1[tx+0] + sdata1[tx+8]);
    }
    __syncthreads(); 
    if (tx<LAPLACE_W && scale<LAPLACE_S-1 && xp<width && (yp+j)<height) 
      d_Result[scale*height*pitch + (yp+j)*pitch + xp] = sdata2[tx] - sdata2[tx+LAPLACE_W];
  }
}

__global__ void LaplaceMultiMemOld(float *d_Image, float *d_Result, int width, int pitch, int height, int octave)
{
  __shared__ float data1[(LAPLACE_W + 2*LAPLACE_R)*LAPLACE_S];
  __shared__ float data2[LAPLACE_W*LAPLACE_S];
  const int tx = threadIdx.x;
  const int xp = blockIdx.x*LAPLACE_W + tx;
  const int yp = blockIdx.y;
  const int scale = threadIdx.y;
  float *kernel = d_LaplaceKernel + octave*12*16 + scale*16; 
  float *sdata1 = data1 + (LAPLACE_W + 2*LAPLACE_R)*scale; 
  float *data = d_Image + max(min(xp - 4, width-1), 0);
  int h = height-1;
  sdata1[tx] = kernel[4]*data[min(yp, h)*pitch] +
    kernel[3]*(data[max(0, min(yp-1, h))*pitch] + data[min(yp+1, h)*pitch]) + 
    kernel[2]*(data[max(0, min(yp-2, h))*pitch] + data[min(yp+2, h)*pitch]) + 
    kernel[1]*(data[max(0, min(yp-3, h))*pitch] + data[min(yp+3, h)*pitch]) + 
    kernel[0]*(data[max(0, min(yp-4, h))*pitch] + data[min(yp+4, h)*pitch]);
  __syncthreads();
  float *sdata2 = data2 + LAPLACE_W*scale; 
  if (tx<LAPLACE_W) {
    sdata2[tx] = kernel[4]*sdata1[tx+4] + 
      kernel[3]*(sdata1[tx+3] + sdata1[tx+5]) + kernel[2]*(sdata1[tx+2] + sdata1[tx+6]) + 
      kernel[1]*(sdata1[tx+1] + sdata1[tx+7]) + kernel[0]*(sdata1[tx+0] + sdata1[tx+8]);
  }
  __syncthreads(); 
  if (tx<LAPLACE_W && scale<LAPLACE_S-1 && xp<width) 
    d_Result[scale*height*pitch + yp*pitch + xp] = sdata2[tx] - sdata2[tx+LAPLACE_W];
}

__global__ void LowPass(float *d_Image, float *d_Result, int width, int pitch, int height)
{
  __shared__ float buffer[(LOWPASS_W + 2*LOWPASS_R)*LOWPASS_H];
  const int tx = threadIdx.x;
  const int ty = threadIdx.y;
  const int xp = blockIdx.x*LOWPASS_W + tx;
  const int yp = blockIdx.y*LOWPASS_H + ty;
  float *kernel = d_LowPassKernel;
  float *data = d_Image + max(min(xp - 4, width-1), 0);
  float *buff = buffer + ty*(LOWPASS_W + 2*LOWPASS_R);
  int h = height-1;
  if (yp<height) 
    buff[tx] = kernel[4]*data[min(yp, h)*pitch] +
      kernel[3]*(data[max(0, min(yp-1, h))*pitch] + data[min(yp+1, h)*pitch]) + 
      kernel[2]*(data[max(0, min(yp-2, h))*pitch] + data[min(yp+2, h)*pitch]) + 
      kernel[1]*(data[max(0, min(yp-3, h))*pitch] + data[min(yp+3, h)*pitch]) + 
      kernel[0]*(data[max(0, min(yp-4, h))*pitch] + data[min(yp+4, h)*pitch]);
  __syncthreads();
  if (tx<LOWPASS_W && xp<width && yp<height)
    d_Result[yp*pitch + xp] = kernel[4]*buff[tx+4] + 
      kernel[3]*(buff[tx+3] + buff[tx+5]) + kernel[2]*(buff[tx+2] + buff[tx+6]) + 
      kernel[1]*(buff[tx+1] + buff[tx+7]) + kernel[0]*(buff[tx+0] + buff[tx+8]);
}

__global__ void LowPassBlockOld(float *d_Image, float *d_Result, int width, int pitch, int height)
{
  __shared__ float xrows[16][32];          
  const int tx = threadIdx.x;
  const int ty = threadIdx.y;
  const int xp = blockIdx.x*LOWPASS_W + tx;
  const int yp = blockIdx.y*LOWPASS_H + ty;
  const int N = 16;
  float *k = d_LowPassKernel;
  int xl = max(min(xp - 4, width-1), 0);
  for (int l=-8;l<=LOWPASS_H;l+=4) {
    if (l<LOWPASS_H) {
      int yl = max(min(yp + l + 4, height-1), 0);
      float val = d_Image[yl*pitch + xl];
      xrows[(l + 8 + ty)%N][tx] = k[4]*ShiftDown(val, 4) +
	k[3]*(ShiftDown(val, 5) + ShiftDown(val, 3)) +
	k[2]*(ShiftDown(val, 6) + ShiftDown(val, 2)) +
	k[1]*(ShiftDown(val, 7) + ShiftDown(val, 1)) +
	k[0]*(ShiftDown(val, 8) + val);
    }
    if (l>=4) {
      int ys = yp + l - 4;
      if (xp<width && ys<height && tx<LOWPASS_W)
	d_Result[ys*pitch + xp] = k[4]*xrows[(l + 0 + ty)%N][tx] +
	     k[3]*(xrows[(l - 1 + ty)%N][tx] + xrows[(l + 1 + ty)%N][tx]) +
	     k[2]*(xrows[(l - 2 + ty)%N][tx] + xrows[(l + 2 + ty)%N][tx]) +
	     k[1]*(xrows[(l - 3 + ty)%N][tx] + xrows[(l + 3 + ty)%N][tx]) +
	     k[0]*(xrows[(l - 4 + ty)%N][tx] + xrows[(l + 4 + ty)%N][tx]);
    }
    if (l>=0)
      __syncthreads();
  }
}

__global__ void LowPassBlock(float *d_Image, float *d_Result, int width, int pitch, int height)
{
  __shared__ float xrows[16][32];          
  const int tx = threadIdx.x;
  const int ty = threadIdx.y;
  const int xp = blockIdx.x*LOWPASS_W + tx;
  const int yp = blockIdx.y*LOWPASS_H + ty;
  const int N = 16;
  float *k = d_LowPassKernel;
  int xl = max(min(xp - 4, width-1), 0);
#pragma unroll
  for (int l=-8;l<4;l+=4) {
    int ly = l + ty;
    int yl = max(min(yp + l + 4, height-1), 0);
    float val = d_Image[yl*pitch + xl];
    val = k[4]*ShiftDown(val, 4) +
      k[3]*(ShiftDown(val, 5) + ShiftDown(val, 3)) +
      k[2]*(ShiftDown(val, 6) + ShiftDown(val, 2)) +
      k[1]*(ShiftDown(val, 7) + ShiftDown(val, 1)) +
      k[0]*(ShiftDown(val, 8) + val);
    xrows[ly + 8][tx] = val;
  }
  __syncthreads();
#pragma unroll
  for (int l=4;l<LOWPASS_H;l+=4) {
    int ly = l + ty;
    int yl = min(yp + l + 4, height-1);
    float val = d_Image[yl*pitch + xl];
    val = k[4]*ShiftDown(val, 4) +
      k[3]*(ShiftDown(val, 5) + ShiftDown(val, 3)) +
      k[2]*(ShiftDown(val, 6) + ShiftDown(val, 2)) +
      k[1]*(ShiftDown(val, 7) + ShiftDown(val, 1)) +
      k[0]*(ShiftDown(val, 8) + val);
    xrows[(ly + 8)%N][tx] = val;
    int ys = yp + l - 4;
    if (xp<width && ys<height && tx<LOWPASS_W)
      d_Result[ys*pitch + xp] = k[4]*xrows[(ly + 0)%N][tx] +
		       k[3]*(xrows[(ly - 1)%N][tx] + xrows[(ly + 1)%N][tx]) +
		       k[2]*(xrows[(ly - 2)%N][tx] + xrows[(ly + 2)%N][tx]) +
		       k[1]*(xrows[(ly - 3)%N][tx] + xrows[(ly + 3)%N][tx]) +
		       k[0]*(xrows[(ly - 4)%N][tx] + xrows[(ly + 4)%N][tx]);
    __syncthreads();
  }
  int ly = LOWPASS_H + ty;
  int ys = yp + LOWPASS_H - 4;
  if (xp<width && ys<height && tx<LOWPASS_W)
    d_Result[ys*pitch + xp] = k[4]*xrows[(ly + 0)%N][tx] +
		     k[3]*(xrows[(ly - 1)%N][tx] + xrows[(ly + 1)%N][tx]) +
		     k[2]*(xrows[(ly - 2)%N][tx] + xrows[(ly + 2)%N][tx]) +
		     k[1]*(xrows[(ly - 3)%N][tx] + xrows[(ly + 3)%N][tx]) +
		     k[0]*(xrows[(ly - 4)%N][tx] + xrows[(ly + 4)%N][tx]);
}

