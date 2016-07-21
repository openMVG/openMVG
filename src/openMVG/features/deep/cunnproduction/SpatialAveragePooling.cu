#include <THC/THC.h>


/*
 * Description:
 *    this function avg-pools an input 3D tensor along dimensions 1 and 2
 *    3D input, 3D output
 */
__global__ void subsample(float *input, float *output, 
                          int input_n, int input_h, int input_w, int output_h, int output_w,
                          int kH, int kW, int dH, int dW)
{
  // iterators
  int xx, yy;

  // compute offsets based on thread/block ID
  int o = blockIdx.x;
  int i = o;

  int xx_start = threadIdx.x;
  int xx_end = output_w;
  int xx_step = blockDim.x;

  int yy_start = blockDim.y*blockIdx.y + threadIdx.y;
  int yy_end = output_h;
  int yy_step = blockDim.y*gridDim.y;

  // select input/output plane
  output = output + o*output_w*output_h;
  input = input + i*input_w*input_h;

  // For all output pixels...
  for(yy = yy_start; yy < yy_end; yy+=yy_step) {
    for(xx = xx_start; xx < xx_end; xx+=xx_step) {
      // Compute the mean of the input image...
      float *ptr_input = input + yy*dH*input_w + xx*dW;
      float *ptr_output = output + yy*output_w + xx;
      float sum = 0;
      int nElements = 0;
      int kx, ky;
      for(ky = 0; ky < kH; ky++) {
        for(kx = 0; kx < kW; kx++) {
          if((xx*dW+kx < input_w) & (yy*dH+ky < input_h)) {
            sum += ptr_input[kx];
            nElements++;
          }
        }
        ptr_input += input_w; // next input line
      }
      // Update output
      *ptr_output = sum / float(nElements);
    }
  }
}

extern "C"
void cunnrelease_SpatialAveragePooling(THCState* state, THCudaTensor* input,
    THCudaTensor* output, int kW, int kH, int dW, int dH, bool is_ceil)
{
  if (input->nDimension == 3) {
    long nInputCols = input->size[2];
    long nInputRows = input->size[1];
    long nOutputRows, nOutputCols;
    if(is_ceil)
    {
      nOutputCols = ceil(float(nInputCols - kW) / float(dW) + 1);
      nOutputRows = ceil(float(nInputRows - kH) / float(dH) + 1);
    }
    else
    {
      nOutputCols = floor(float(nInputCols - kW) / float(dW) + 1);
      nOutputRows = floor(float(nInputRows - kH) / float(dH) + 1);
    }
    long nInputPlane = input->size[0];

    input = THCudaTensor_newContiguous(state, input);
    float* input_data = THCudaTensor_data(state, input);

    THCudaTensor_resize3d(state, output, nInputPlane, nOutputRows, nOutputCols);
    float* output_data = THCudaTensor_data(state, output);

    // cuda blocks & threads:
    int yblocks = (int)(16L / nInputPlane);
    yblocks = yblocks < 1 ? 1 : yblocks;
    dim3 blocks(nInputPlane,yblocks);
    dim3 threads(32,8);

    // run subsample kernel
    subsample <<<blocks, threads, 0, THCState_getCurrentStream(state)>>>
      (input_data, output_data,
      nInputPlane, nInputRows, nInputCols, nOutputRows, nOutputCols, kH, kW, dH, dW);
  } else {
    long nInputCols = input->size[3];
    long nInputRows = input->size[2];
    long nbatch = input->size[0];
    long nOutputRows, nOutputCols;
    if(is_ceil)
    {
      nOutputCols = ceil(float(nInputCols - kW) / float(dW) + 1);
      nOutputRows = ceil(float(nInputRows - kH) / float(dH) + 1);
    }
    else
    {
      nOutputCols = floor(float(nInputCols - kW) / float(dW) + 1);
      nOutputRows = floor(float(nInputRows - kH) / float(dH) + 1);
    }
    long nInputPlane = input->size[1];

    input = THCudaTensor_newContiguous(state, input);
    float* input_data = THCudaTensor_data(state, input);

    THCudaTensor_resize4d(state, output, nbatch, nInputPlane, nOutputRows, nOutputCols);
    float* output_data = THCudaTensor_data(state, output);

    // cuda blocks & threads:
    int yblocks = (int)(16L / nInputPlane);
    yblocks = yblocks < 1 ? 1 : yblocks;
    dim3 blocks(nInputPlane*nbatch,yblocks);
    dim3 threads(32,8);

    // run subsample kernel
    subsample <<<blocks, threads, 0, THCState_getCurrentStream(state)>>>
      (input_data, output_data,
      nInputPlane, nInputRows, nInputCols, nOutputRows, nOutputCols, kH, kW, dH, dW);
  }

  // clean
  THCudaTensor_free(state, input);

  // check for errors
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    printf("error in SpatialSubsampling.updateOutput: %s\n", cudaGetErrorString(err));
    THError("aborting");
  }
}

