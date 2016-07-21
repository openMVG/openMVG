CUNN-PRODUCTION
=====

A simple CUDA/C++ wrapper over some cunn modules to enable easy embedding of Torch7 networks in C++ projects. Only depends on TH and THC.

The following modules supported for now:
```
nn.Sequential
nn.Concat
nn.Parallel
nn.SpatialConvolutionMM
nn.SpatialMaxPooling
nn.SpatialAveragePooling
nn.ReLU
nn.Linear
```

I tried to follow Torch7 nn architecture. There is an abstract cunn::Module with forward() function and all the modules are inherited from it, so syntax has to be familiar to torch users:

```c++
auto net = std::make_shared<cunn::Sequential>();
net->add(std::make_shared<cunn::SpatialConvolutionMM(state,3,96,11,11,3,3)>);
net->add(std::make_shared<cunn::ReLU>(state));
net->add(std::make_shared<cunn::SpatialMaxPooling>(state, 2,2,2,2));
net->add(std::make_shared<cunn::SpatialConvolutionMM>(state,96,192,5,5));
net->add(std::make_shared<cunn::ReLU>(state));

THCudaTensor *input = THCudaTensor_newWithSize4d(state, 1,3,224,224)
THCudaTensor *output = net->forward(input);
```

Only CUDA 7 supported.
