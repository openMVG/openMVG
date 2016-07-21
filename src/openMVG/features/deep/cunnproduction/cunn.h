#ifndef CUNN_H
#define CUNN_H

#include <THC/THC.h>
#include <memory>
#include <vector>
#include <string>

namespace cunn {

/*
 * Abstract class as nn.Module
 */
class Module {
public:
  typedef std::shared_ptr<Module> Ptr;

  Module(THCState *state);
  ~Module();

  virtual THCudaTensor* forward(THCudaTensor *input) = 0;
  virtual const std::string tostring() const { return std::string("name not defined"); }
  THCState *state;
  THCudaTensor *output;
};

/*
 * nn.Sequential
 */
class Sequential : public Module {
public:
  typedef std::shared_ptr<Sequential> Ptr;

  Sequential(THCState* state);
  ~Sequential();

  inline Module::Ptr get(int i) const { return modules[i];  }
  void add(Module::Ptr module);

  THCudaTensor* forward(THCudaTensor* input);
  const std::string tostring() const;

  std::vector<Module::Ptr> modules;
};

/*
 * nn.Concat
 */
class Concat : public Module {
public:
  typedef std::shared_ptr<Concat> Ptr;

  Concat(THCState *state, int dim);
  ~Concat();

  inline Module::Ptr get(int i) const { return modules[i];  }
  void add(Module::Ptr module);

  THCudaTensor* forward(THCudaTensor* input);
  const std::string tostring() const;

  std::vector<Module::Ptr> modules;
  int dim;
};

/*
 * nn.Parallel
 */
class Parallel : public Module {
public:
  typedef std::shared_ptr<Parallel> Ptr;

  Parallel(THCState *state, int input_dim, int output_dim);
  ~Parallel();

  inline Module::Ptr get(int i) const { return modules[i];  }
  void add(Module::Ptr module);

  THCudaTensor* forward(THCudaTensor* input);
  const std::string tostring() const;

  std::vector<Module::Ptr> modules;
  int input_dim, output_dim;
};

/*
 * nn.SpatialConvolutionMM
 */
class SpatialConvolutionMM : public Module {
public:
  SpatialConvolutionMM(THCState *state, int nInputPlane, int nOutputPlane, int kW, int kH, int dW = 1, int dH = 1, int padding = 0);
  // Do NOT use this method. Just some strange linker errors between externals and this.
  //SpatialConvolutionMM(THCState *state, int nInputPlane, int nOutputPlane, int kW, int kH, int dW = 1, int dH = 1, int padding = 0, int samp = 0);
 
  ~SpatialConvolutionMM();

  THCudaTensor* forward(THCudaTensor *input);
  inline const std::string tostring() const { return std::string("cunn.SpatialConvolutionMM"); }

  THCudaTensor *weight, *bias;
  THCudaTensor *finput, *fgradinput;
  int nInputPlane, nOutputPlane, kW, kH, dW, dH, padding;
};

/*
 * nn.SpatialMaxPooling
 */
class SpatialMaxPooling : public Module {
public:
  SpatialMaxPooling(THCState *state, int kW, int kH, int dW, int dH, bool is_ceil = false);
  ~SpatialMaxPooling();

  THCudaTensor* forward(THCudaTensor *input);
  inline const std::string tostring() const { return std::string("cunn.SpatialMaxPooling"); }

  int kW, kH, dW, dH;
  bool is_ceil;
};

/*
 * nn.SpatialAveragePooling
 */
class SpatialAveragePooling : public Module {
public:
  SpatialAveragePooling(THCState *state, int kW, int kH, int dW, int dH, bool is_ceil = false);
  ~SpatialAveragePooling();

  THCudaTensor* forward(THCudaTensor *input);
  inline const std::string tostring() const { return std::string("cunn.SpatialAveragePooling"); }

  int kW, kH, dW, dH;
  bool is_ceil;
};

/*
 * nn.ReLU
 */
class ReLU : public Module {
public:
  ReLU(THCState *state);
  ~ReLU();

  THCudaTensor* forward(THCudaTensor *input);
  inline const std::string tostring() const { return std::string("cunn.ReLU"); }
};

/*
 * nn.Linear
 */
class Linear : public Module {
public:
  typedef std::shared_ptr<Linear> Ptr;

  Linear(THCState *state, int nInputPlane, int nOutputPlane);
  ~Linear();

  THCudaTensor* forward(THCudaTensor *input);
  inline const std::string tostring() const { return std::string("cunn.Linear"); }

  THCudaTensor *weight, *bias, *buffer;
  int nOutputPlane, nInputPlane;
};

/*
 * nn.Reshape
 */
class Reshape : public Module {
public:
  Reshape(THCState *state, const std::vector<size_t>& sizes);
  ~Reshape();

  THCudaTensor* forward(THCudaTensor* input);
  inline const std::string tostring() const { return std::string("cunn.Reshape"); }

  std::vector<size_t> sizes;
};

}
#endif // CUNN_H
