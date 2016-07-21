#include <THC/THC.h>


/*
 * Description:
 *    this function maxpools an input 4D tensor along dimensions 2 and 3
 *    4D input, 4D output, 4D argmax x and y 
 */
__global__ void maxpool(float *input, float *output,
                        int input_n, int input_h, int input_w, int output_h, int output_w,
                        int kH, int kW, int dH, int dW)
{
  // iterators
  int xx, yy;

  // compute offsets based on thread/block ID
  int o = blockIdx.x;
  int i = o;
  //int k = blockIdx.x % input_n;

  int xx_start = threadIdx.x;
  int xx_end = output_w;
  const int xx_step = blockDim.x;

  int yy_start = blockDim.y*blockIdx.y + threadIdx.y;
  int yy_end = output_h;
  const int yy_step = blockDim.y*gridDim.y;

  // select input/output plane
  output = output + o*output_w*output_h;
  input = input + i*input_w*input_h;

  // For all output pixels...
  for(yy = yy_start; yy < yy_end; yy+=yy_step) {
    for(xx = xx_start; xx < xx_end; xx+=xx_step) {
      // Compute the mean of the input image...
      float *ptr_input = input + yy*dH*input_w + xx*dW;
      float *ptr_output = output + yy*output_w + xx;
      float max = -FLT_MAX;
      int kx, ky;
      for(ky = 0; ky < kH; ky++) {
        for(kx = 0; kx < kW; kx++) {
          if((xx*dW+kx < input_w) & (yy*dH+ky < input_h))
	    max = fmaxf(max, ptr_input[kx]);
        }
        ptr_input += input_w; // next input line
      }
      // Update output and argmax
      *ptr_output = max;
    }
  }
}


extern "C"
void cunnrelease_SpatialMaxPooling(THCState* state, THCudaTensor* input, 
    THCudaTensor* output, int kW, int kH, int dW, int dH, bool is_ceil)
{
  float *output_data;
  float *input_data;

  //luaL_argcheck(L, input->nDimension == 3 || input->nDimension == 4, 2, "3D or 4D (batch) tensor expected");

  if (input->nDimension == 3) {
    long nInputCols = input->size[2];
    long nInputRows = input->size[1];
    long nInputPlane = input->size[0];
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

    //luaL_argcheck(L, nInputCols >= kW && nInputRows >= kH, 2, "input image smaller than kernel size");

    input = THCudaTensor_newContiguous(state, input);
    input_data = THCudaTensor_data(state, input);

    THCudaTensor_resize3d(state, output, nInputPlane, nOutputRows, nOutputCols);
    
    output_data = THCudaTensor_data(state, output);

    // cuda blocks & threads:
    int yblocks = (int)(16L / nInputPlane);
    yblocks = yblocks < 1 ? 1 : yblocks;
    dim3 blocks(nInputPlane,yblocks);
    dim3 threads(32,8);

    // run maxpool kernel
    maxpool <<<blocks, threads, 0, THCState_getCurrentStream(state)>>>
      (input_data, output_data, 
      nInputPlane, nInputRows, nInputCols, nOutputRows, nOutputCols, kH, kW, dH, dW);
  } else {
    long nInputCols = input->size[3];
    long nInputRows = input->size[2];
    long nInputPlane = input->size[1];
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

    //luaL_argcheck(L, nInputCols >= kW && nInputRows >= kH, 2, "input image smaller than kernel size");

    input = THCudaTensor_newContiguous(state, input);
    input_data = THCudaTensor_data(state, input);

    THCudaTensor_resize4d(state, output, nbatch, nInputPlane, nOutputRows, nOutputCols);

    output_data = THCudaTensor_data(state, output);

    // cuda blocks & threads:
    int yblocks = (int)(16L / nInputPlane);
    yblocks = yblocks < 1 ? 1 : yblocks;
    dim3 blocks(nInputPlane*nbatch,yblocks);
    dim3 threads(32,8);

    // run maxpool kernel
    maxpool <<<blocks, threads, 0, THCState_getCurrentStream(state)>>>
      (input_data, output_data,
      nInputPlane, nInputRows, nInputCols, nOutputRows, nOutputCols, kH, kW, dH, dW);
  }

  // clean
  THCudaTensor_free(state, input);

  // check for errors
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    printf("error in SpatialMaxsampling.updateOutput: %s\n", cudaGetErrorString(err));
    THError("aborting");
  }
}

