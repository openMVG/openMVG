#include "cunn.h"

#include <sstream>

extern "C" {
void cunnrelease_SpatialConvolution(THCState *state,
    THCudaTensor *input,
    THCudaTensor *weight,
    THCudaTensor *bias,
    THCudaTensor *columns,
    THCudaTensor *ones,
    THCudaTensor *output,
    int nInputPlane, int nOutputPlane, int kW, int kH, int dW, int dH, int padding);

void cunnrelease_SpatialMaxPooling(THCState* state,
    THCudaTensor* input, 
    THCudaTensor* output,
    int kW, int kH, int dW, int dH, bool is_ceil);

void cunnrelease_SpatialAveragePooling(THCState* state,
    THCudaTensor* input,
    THCudaTensor* output,
    int kW, int kH, int dW, int dH, bool is_ceil);

void cunnrelease_Linear(THCState *state,
    THCudaTensor *input,
    THCudaTensor *output,
    THCudaTensor *weight,
    THCudaTensor *bias,
    THCudaTensor *buffer);

void cunnrelease_ReLUIP(THCState *state,
    THCudaTensor *input);
}

namespace cunn {

Module::Module(THCState *state) : state(state) {
  output = THCudaTensor_new(state);
}

Module::~Module() {
  THCudaTensor_free(state, output);
}


SpatialConvolutionMM::SpatialConvolutionMM(THCState *state,
    int nInputPlane, int nOutputPlane, int kW, int kH, int dW, int dH, int padding) :
    	Module(state), nInputPlane(nInputPlane), nOutputPlane(nOutputPlane), kW(kW), kH(kH), dW(dW), dH(dH), padding(padding)  {
  weight = THCudaTensor_newWithSize2d(state, nOutputPlane, nInputPlane*kW*kH);
  bias = THCudaTensor_newWithSize1d(state, nOutputPlane);
  finput = THCudaTensor_new(state);
  fgradinput = THCudaTensor_new(state);
}
/*
SpatialConvolutionMM::SpatialConvolutionMM(THCState *state,
    int nInputPlane, int nOutputPlane, int kW, int kH, int dW, int dH, int padding, int samp) :
    	Module(state), nInputPlane(nInputPlane), nOutputPlane(nOutputPlane), kW(kW), kH(kH), dW(dW), dH(dH), padding(padding)  {
  weight = THCudaTensor_newWithSize2d(state, nOutputPlane, nInputPlane*kW*kH);
  bias = THCudaTensor_newWithSize1d(state, nOutputPlane);
  finput = THCudaTensor_new(state);
  fgradinput = THCudaTensor_new(state);
}
*/
THCudaTensor*
SpatialConvolutionMM::forward(THCudaTensor *input)
{
  cunnrelease_SpatialConvolution(state, input, weight, bias, finput, fgradinput, output,
      nInputPlane, nOutputPlane, kW, kH, dW, dH, padding);
  return output;
}

SpatialConvolutionMM::~SpatialConvolutionMM()
{
  THCudaTensor_free(state, weight);
  THCudaTensor_free(state, bias);
  THCudaTensor_free(state, finput);
  THCudaTensor_free(state, fgradinput);
}

SpatialMaxPooling::SpatialMaxPooling(THCState *state, int kW, int kH, int dW, int dH, bool is_ceil) :
  Module(state), kW(kW), kH(kH), dW(dW), dH(dH), is_ceil(is_ceil) {}

SpatialMaxPooling::~SpatialMaxPooling() {}

THCudaTensor*
SpatialMaxPooling::forward(THCudaTensor *input)
{
  cunnrelease_SpatialMaxPooling(state, input, output, kW, kH, dW, dH, is_ceil);
  return output;
}

SpatialAveragePooling::SpatialAveragePooling(THCState *state, int kW, int kH, int dW, int dH, bool is_ceil) :
  Module(state), kW(kW), kH(kH), dW(dW), dH(dH), is_ceil(is_ceil) {}

SpatialAveragePooling::~SpatialAveragePooling() {}

THCudaTensor*
SpatialAveragePooling::forward(THCudaTensor *input)
{
  cunnrelease_SpatialAveragePooling(state, input, output, kW, kH, dW, dH, is_ceil);
  return output;
}


ReLU::ReLU(THCState *state) : Module(state) {}

ReLU::~ReLU() {}

THCudaTensor*
ReLU::forward(THCudaTensor *input)
{
  cunnrelease_ReLUIP(state, input);
  return input;
}


Linear::Linear(THCState *state, int nInputPlane, int nOutputPlane) :
  Module(state), nInputPlane(nInputPlane), nOutputPlane(nOutputPlane)
{
  weight = THCudaTensor_newWithSize2d(state, nOutputPlane, nInputPlane);
  bias = THCudaTensor_newWithSize1d(state, nOutputPlane);
  buffer = THCudaTensor_new(state);
}

Linear::~Linear()
{
  THCudaTensor_free(state, weight);
  THCudaTensor_free(state, bias);
  THCudaTensor_free(state, buffer);
}

THCudaTensor*
Linear::forward(THCudaTensor *input)
{
  cunnrelease_Linear(state, input, output, weight, bias, buffer);
  return output;
}


Sequential::Sequential(THCState *state) : Module(state) {}

Sequential::~Sequential() {}

void
Sequential::add(Module::Ptr module)
{
  modules.push_back(module);
}

THCudaTensor*
Sequential::forward(THCudaTensor* input)
{
  THCudaTensor* out = input;
  for(auto& it: modules)
    out = it->forward(out);
  return out;
}

const std::string
Sequential::tostring() const
{
  std::stringstream s;
  s << "nn.Sequential {\n";
  int counter = 1;
  for(auto &it: modules)
    s << "  (" << counter++ << ") " <<  it->tostring() << std::endl;
  s << "}\n";
  return s.str();
}


Concat::Concat(THCState *state, int dim) : Module(state), dim(dim) {}

Concat::~Concat() {}

void
Concat::add(Module::Ptr module)
{
  modules.push_back(module);
}

THCudaTensor*
Concat::forward(THCudaTensor *input)
{
  std::vector<THCudaTensor*> outputs;
  for(auto &it: modules)
    outputs.push_back(it->forward(input));

  std::vector<long> incs(outputs.size(), 0);
  long total_elems = 0;
  long total_dim_elems = 0;
  for(auto &jt: outputs)
  {
    total_dim_elems += jt->size[dim];
    total_elems += THCudaTensor_nElement(state, jt);
  }
  for(size_t i=1; i<outputs.size(); ++i)
    incs[i] = incs[i-1] + outputs[i]->size[dim];
  
  if(total_elems != THCudaTensor_nElement(state, output))
  {
    THLongStorage* size = THCudaTensor_newSizeOf(state, outputs.front());
    THLongStorage_data(size)[dim] = total_dim_elems;
    THCudaTensor_resize(state, output, size, NULL);
    THLongStorage_free(size);
  }

  for(size_t i=0; i<outputs.size(); ++i)
  {
    THCudaTensor* cur = THCudaTensor_newNarrow(state, output, dim, incs[i], outputs[i]->size[dim]);
    THCudaTensor_copy(state, cur, outputs[i]);
    THCudaTensor_free(state, cur);
  }
  return output;
}

const std::string
Concat::tostring() const
{
  std::stringstream s;
  s << "nn.Concat(" << dim << ") {\n";
  int counter = 1;
  for(auto &it: modules)
    s << "  (" << counter++ << ") " <<  it->tostring() << std::endl;
  s << "}\n";
  return s.str();
}


Parallel::Parallel(THCState *state, int input_dim, int output_dim) : 
  Module(state), input_dim(input_dim), output_dim(output_dim) {}

Parallel::~Parallel() {}

void
Parallel::add(Module::Ptr module)
{
  modules.push_back(module);
}

THCudaTensor*
Parallel::forward(THCudaTensor* input)
{
  std::vector<THCudaTensor*> outputs;

  int counter = 0;
  for(auto &it: modules)
  {
    THCudaTensor *in = THCudaTensor_newSelect(state, input, input_dim, counter++);
    outputs.push_back(it->forward(in));
    THCudaTensor_free(state, in);
  }

  std::vector<long> incs(outputs.size(), 0);
  long total_elems = 0;
  long total_dim_elems = 0;
  for(auto &jt: outputs)
  {
    total_dim_elems += jt->size[output_dim];
    total_elems += THCudaTensor_nElement(state, jt);
  }
  for(size_t i=1; i<outputs.size(); ++i)
    incs[i] = incs[i-1] + outputs[i]->size[output_dim];
  
  if(total_elems != THCudaTensor_nElement(state, output))
  {
    THLongStorage* size = THCudaTensor_newSizeOf(state, outputs.front());
    THLongStorage_data(size)[output_dim] = total_dim_elems;
    THCudaTensor_resize(state, output, size, NULL);
    THLongStorage_free(size);
  }

  for(size_t i=0; i<outputs.size(); ++i)
  {
    THCudaTensor* cur = THCudaTensor_newNarrow(state, output, output_dim, incs[i], outputs[i]->size[output_dim]);
    THCudaTensor_copy(state, cur, outputs[i]);
    THCudaTensor_free(state, cur);
  }
  return output;
}

const std::string
Parallel::tostring() const
{
  std::stringstream s;
  s << "nn.Parallel(" << input_dim << ", " << output_dim << ") {\n";
  int counter = 1;
  for(auto &it: modules)
    s << "  (" << counter++ << ") " <<  it->tostring() << std::endl;
  s << "}\n";
  return s.str();
}

Reshape::Reshape(THCState *state, const std::vector<size_t>& sizes) : Module(state), sizes(sizes) {}

Reshape::~Reshape() {}

THCudaTensor*
Reshape::forward(THCudaTensor* input)
{
  //size_t ndim = THCudaTensor_nDimension(state, input);
  // support only one case for now
 /* 
  THLongStorage *size = THCudaTensor_newSizeOf(state, input);
  THLongStorage_resize(size, sizes.size() + 1);
  long* data = THLongStorage_data(size);
  for(size_t i=1; i < sizes.size() + 1; ++i)
    data[i] = sizes[i-1];
  
  //THCudaTensor_resize2d(state, output, input->size[0], sizes[0]);
  //THCudaTensor_copy(state, output, input);
  THCudaTensor_resize(state, output, size, NULL);
  THCudaTensor_copy(state, output, input);
  */
  if(sizes.size() == 3)
    THCudaTensor_resize4d(state, output, input->size[0], sizes[0], sizes[1], sizes[2]);
  else if(sizes.size() == 1)
    THCudaTensor_resize2d(state, output, input->size[0], sizes[0]);
  THCudaTensor_copy(state, output, input);
  return output; 
}

} // namespace cunn
