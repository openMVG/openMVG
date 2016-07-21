#include <THC/THC.h>

extern "C"
void cunnrelease_Linear(THCState *state,
    THCudaTensor *input,
    THCudaTensor *output,
    THCudaTensor *weight,
    THCudaTensor *bias,
    THCudaTensor *buffer)
{
  long nOutputPlane = weight->size[0];
  long bs = input->size[0];

  if(THCudaTensor_nDimension(state, buffer) == 0 || buffer->size[0] != bs)
  {
    THCudaTensor_resize1d(state, buffer, bs);
    THCudaTensor_fill(state, buffer, 1);
  }

  THCudaTensor_resize2d(state, output, bs, nOutputPlane);

  weight = THCudaTensor_newTranspose(state, weight, 0, 1);
  THCudaTensor_addmm(state, output, 0, output, 1, input, weight);
  THCudaTensor_addr(state, output, 1, output, 1, buffer, bias); 
  THCudaTensor_free(state, weight);
}
