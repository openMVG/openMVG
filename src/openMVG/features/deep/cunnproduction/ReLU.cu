#include <THC/THC.h>
#include <THC/THCApply.cuh>

struct ReLUUpdateOutputIP {
  __device__ __forceinline__ void operator()(float* x) {
    *x = (*x > 0.f) ? *x : 0.f;
  }
};

extern "C"
void cunnrelease_ReLUIP(THCState *state, THCudaTensor *input)
{
  THC_pointwiseApply1(state, input, ReLUUpdateOutputIP());
}
