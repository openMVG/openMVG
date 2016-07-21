// Copyright 2015 Sergey Zagoruyko, Nikos Komodakis
// sergey.zagoruyko@imagine.enpc.fr, nikos.komodakis@enpc.fr
// Ecole des Ponts ParisTech, Universite Paris-Est, IMAGINE
//
// The software is free to use only for non-commercial purposes.
// IF YOU WOULD LIKE TO USE IT FOR COMMERCIAL PURPOSES, PLEASE CONTACT
// Prof. Nikos Komodakis (nikos.komodakis@enpc.fr)

#include "loader.h"

#include <iostream>
#include <cunn.h>
#include <sstream>

#define CHECK_EQ(CONDITION) \
  do { \
    if( !(CONDITION) ) \
      printf("%d: assertion failed\n", __LINE__); \
  } while (0)

// Additional simple module for iris cropping
namespace cunn {

class CentralCrop : public Module {
public:
  CentralCrop(THCState *state, int dimx, int dimy) : Module(state), dimx(dimx), dimy(dimy) {}
  ~CentralCrop() {}

  THCudaTensor* forward(THCudaTensor *input)
  {
    THCudaTensor* crop1 = THCudaTensor_newNarrow(state, input, 2, dimx/4, dimx/2);
    THCudaTensor* crop2 = THCudaTensor_newNarrow(state, crop1, 3, dimy/4, dimy/2);
    THCudaTensor_resizeAs(state, output, crop2);
    THCudaTensor_copy(state, output, crop2);
    THCudaTensor_free(state, crop2);
    THCudaTensor_free(state, crop1);
    return output;
  }

  int dimx, dimy;
};

}

// read weights for one module from file
template<class T>
void readParameters(THCState *state, T* casted, FILE* f)
{
  THFloatTensor *weight = THFloatTensor_newWithSize1d(THCudaTensor_nElement(state, casted->weight));
  THFloatTensor *bias = THFloatTensor_newWithSize1d(THCudaTensor_nElement(state, casted->bias));
  int w_numel = 0, b_numel = 0;
  fread(&w_numel, sizeof(int), 1, f);
  fread(&b_numel, sizeof(int), 1, f);

  CHECK_EQ(THFloatTensor_nElement(weight) == w_numel);
  CHECK_EQ(THFloatTensor_nElement(bias) == b_numel);
  //std::cout << THFloatTensor_nElement(weight) << " " << w_numel << std::endl;
  //std::cout << THFloatTensor_nElement(bias) << " " << b_numel << std::endl;

  fread(THFloatTensor_data(weight), sizeof(float)*THFloatTensor_nElement(weight), 1, f);
  fread(THFloatTensor_data(bias), sizeof(float)*THFloatTensor_nElement(bias), 1, f);
  
  THCudaTensor_copyFloat(state, casted->weight, weight);
  THCudaTensor_copyFloat(state, casted->bias, bias);
  THFloatTensor_free(weight);
  THFloatTensor_free(bias);
}

cunn::Concat::Ptr createSiam2streamBranch(THCState *state)
{
  std::vector<size_t> size = {1, 64, 64};
  cunn::Sequential::Ptr sub_branch1_1 = std::make_shared<cunn::Sequential>(state);
  sub_branch1_1->add(std::make_shared<cunn::Reshape>(state, size));
  sub_branch1_1->add(std::make_shared<cunn::CentralCrop>(state, 64, 64));
  sub_branch1_1->add(std::make_shared<cunn::SpatialConvolutionMM>(state, 1,96, 4,4, 2,2));
  sub_branch1_1->add(std::make_shared<cunn::ReLU>(state));
  sub_branch1_1->add(std::make_shared<cunn::SpatialMaxPooling>(state, 2,2,2,2));
  sub_branch1_1->add(std::make_shared<cunn::SpatialConvolutionMM>(state, 96,192, 3,3));
  sub_branch1_1->add(std::make_shared<cunn::ReLU>(state));
  sub_branch1_1->add(std::make_shared<cunn::SpatialConvolutionMM>(state, 192,256, 3,3));
  sub_branch1_1->add(std::make_shared<cunn::ReLU>(state));
  sub_branch1_1->add(std::make_shared<cunn::SpatialConvolutionMM>(state, 256,256, 3,3));
  sub_branch1_1->add(std::make_shared<cunn::Reshape>(state, std::vector<size_t>(1,256)));

  cunn::Sequential::Ptr sub_branch1_2 = std::make_shared<cunn::Sequential>(state);
  sub_branch1_2->add(std::make_shared<cunn::Reshape>(state, size));
  sub_branch1_2->add(std::make_shared<cunn::SpatialAveragePooling>(state, 2,2,2,2));
  sub_branch1_2->add(std::make_shared<cunn::SpatialConvolutionMM>(state, 1,96, 4,4, 2,2));
  sub_branch1_2->add(std::make_shared<cunn::ReLU>(state));
  sub_branch1_2->add(std::make_shared<cunn::SpatialMaxPooling>(state, 2,2,2,2));
  sub_branch1_2->add(std::make_shared<cunn::SpatialConvolutionMM>(state, 96,192, 3,3));
  sub_branch1_2->add(std::make_shared<cunn::ReLU>(state));
  sub_branch1_2->add(std::make_shared<cunn::SpatialConvolutionMM>(state, 192,256, 3,3));
  sub_branch1_2->add(std::make_shared<cunn::ReLU>(state));
  sub_branch1_2->add(std::make_shared<cunn::SpatialConvolutionMM>(state, 256,256, 3,3));
  sub_branch1_2->add(std::make_shared<cunn::Reshape>(state, std::vector<size_t>(1,256)));

  cunn::Concat::Ptr sub_branch1 = std::make_shared<cunn::Concat>(state, 1);

  sub_branch1->add(sub_branch1_1);
  sub_branch1->add(sub_branch1_2);

  return sub_branch1;
}

// load the full network
cunn::Sequential::Ptr
loadNetwork(THCState* state, const char* filename)
{
  cunn::Sequential::Ptr net = std::make_shared<cunn::Sequential>(state);

  FILE *f = fopen(filename, "rb");
  if(f == NULL)
  {
    std::stringstream s;
    s << "file " << filename << " not found";
    THError(s.str().c_str());
    return net;
  }

  std::string net_type;
  {
    int type_size;
    fread(&type_size, sizeof(int), 1, f);
    std::vector<unsigned char> type(type_size);
    fread(type.data(), sizeof(char) * type_size, 1, f);
    net_type = std::string(type.begin(), type.end());
  }
  //std::cout << "Reading file with: " << net_type << std::endl;
  
  if(net_type == "2ch")
  {
    net->add(std::make_shared<cunn::SpatialConvolutionMM>(state, 2,96, 7,7, 3,3));
    net->add(std::make_shared<cunn::ReLU>(state));
    net->add(std::make_shared<cunn::SpatialMaxPooling>(state, 2,2, 2,2));
    net->add(std::make_shared<cunn::SpatialConvolutionMM>(state, 96, 192, 5,5));
    net->add(std::make_shared<cunn::ReLU>(state));
    net->add(std::make_shared<cunn::SpatialMaxPooling>(state, 2,2, 2,2));
    net->add(std::make_shared<cunn::SpatialConvolutionMM>(state, 192, 256, 3,3));
    net->add(std::make_shared<cunn::ReLU>(state));
    net->add(std::make_shared<cunn::Reshape>(state, std::vector<size_t>(1,256)));
    net->add(std::make_shared<cunn::Linear>(state, 256, 1));

    for(int i=0; i<4; ++i)
    {
      int net_j = 0;
      fread(&net_j, sizeof(int), 1, f);
      CHECK_EQ(net_j >= 1);
      //std::cout << "Reading weights of layer " << net_j << std::endl;
      auto module = net->get(net_j-1);
      if(i!=3)
	readParameters(state, (cunn::SpatialConvolutionMM*)module.get(), f);
      else
	readParameters(state, (cunn::Linear*)module.get(), f);
    }
  }
  else if(net_type == "2chdeep")
  {
    net->add(std::make_shared<cunn::SpatialConvolutionMM>(state, 2,96, 4,4, 3,3));
    net->add(std::make_shared<cunn::ReLU>(state));
    net->add(std::make_shared<cunn::SpatialConvolutionMM>(state, 96,96, 3,3));
    net->add(std::make_shared<cunn::ReLU>(state));
    net->add(std::make_shared<cunn::SpatialConvolutionMM>(state, 96,96, 3,3));
    net->add(std::make_shared<cunn::ReLU>(state));
    net->add(std::make_shared<cunn::SpatialConvolutionMM>(state, 96,96, 3,3));
    net->add(std::make_shared<cunn::ReLU>(state));

    net->add(std::make_shared<cunn::SpatialMaxPooling>(state, 2,2, 2,2, true));

    net->add(std::make_shared<cunn::SpatialConvolutionMM>(state, 96,192, 3,3));
    net->add(std::make_shared<cunn::ReLU>(state));
    net->add(std::make_shared<cunn::SpatialConvolutionMM>(state, 192,192, 3,3));
    net->add(std::make_shared<cunn::ReLU>(state));
    net->add(std::make_shared<cunn::SpatialConvolutionMM>(state, 192,192, 3,3));
    net->add(std::make_shared<cunn::ReLU>(state));

    net->add(std::make_shared<cunn::Reshape>(state, std::vector<size_t>(1,192*2*2)));
    net->add(std::make_shared<cunn::Linear>(state, 192*2*2, 1));

    for(int i=0; i<8; ++i)
    {
      int net_j = 0;
      fread(&net_j, sizeof(int), 1, f);
      CHECK_EQ(net_j >= 1);
      //std::cout << "Reading weights of layer " << net_j << std::endl;
      auto module = net->get(net_j-1);
      if(i!=7)
	readParameters(state, (cunn::SpatialConvolutionMM*)module.get(), f);
      else
	readParameters(state, (cunn::Linear*)module.get(), f);
    }
  }
  else if(net_type == "2ch2stream")
  {
    int featureOut = 192*2*2;

    cunn::Sequential::Ptr branch1 = std::make_shared<cunn::Sequential>(state);
    branch1->add(std::make_shared<cunn::SpatialAveragePooling>(state, 2,2,2,2));
    branch1->add(std::make_shared<cunn::SpatialConvolutionMM>(state, 2, 96, 5,5));
    branch1->add(std::make_shared<cunn::ReLU>(state));
    branch1->add(std::make_shared<cunn::SpatialMaxPooling>(state, 2,2, 2,2));
    branch1->add(std::make_shared<cunn::SpatialConvolutionMM>(state, 96, 96, 3,3));
    branch1->add(std::make_shared<cunn::ReLU>(state));
    branch1->add(std::make_shared<cunn::SpatialMaxPooling>(state, 2,2, 2,2));
    branch1->add(std::make_shared<cunn::SpatialConvolutionMM>(state, 96, 192, 3,3));
    branch1->add(std::make_shared<cunn::ReLU>(state));
    branch1->add(std::make_shared<cunn::SpatialConvolutionMM>(state, 192, 192, 3,3));
    branch1->add(std::make_shared<cunn::ReLU>(state));
    branch1->add(std::make_shared<cunn::Reshape>(state, std::vector<size_t>(1,featureOut)));

    cunn::Sequential::Ptr branch2 = std::make_shared<cunn::Sequential>(state);
    branch2->add(std::make_shared<cunn::CentralCrop>(state, 64, 64));
    branch2->add(std::make_shared<cunn::SpatialConvolutionMM>(state, 2, 96, 5,5));
    branch2->add(std::make_shared<cunn::ReLU>(state));
    branch2->add(std::make_shared<cunn::SpatialMaxPooling>(state, 2,2, 2,2));
    branch2->add(std::make_shared<cunn::SpatialConvolutionMM>(state, 96, 96, 3,3));
    branch2->add(std::make_shared<cunn::ReLU>(state));
    branch2->add(std::make_shared<cunn::SpatialMaxPooling>(state, 2,2, 2,2));
    branch2->add(std::make_shared<cunn::SpatialConvolutionMM>(state, 96, 192, 3,3));
    branch2->add(std::make_shared<cunn::ReLU>(state));
    branch2->add(std::make_shared<cunn::SpatialConvolutionMM>(state, 192, 192, 3,3));
    branch2->add(std::make_shared<cunn::ReLU>(state));
    branch2->add(std::make_shared<cunn::Reshape>(state, std::vector<size_t>(1,featureOut)));

    cunn::Concat::Ptr concat = std::make_shared<cunn::Concat>(state, 1);
    concat->add(branch1);
    concat->add(branch2);

    net->add(concat);

    net->add(std::make_shared<cunn::Reshape>(state, std::vector<size_t>(1,featureOut*2)));
    net->add(std::make_shared<cunn::Linear>(state, featureOut*2, featureOut));
    net->add(std::make_shared<cunn::ReLU>(state));
    net->add(std::make_shared<cunn::Linear>(state, featureOut, 1));

    for(int i=0; i<4; ++i)
    {
      int net_j = 0;
      fread(&net_j, sizeof(int), 1, f);
      //std::cout << "Reading weights of layer " << net_j << std::endl;
      auto module = branch1->get(net_j-1);
      readParameters(state, (cunn::SpatialConvolutionMM*)module.get(), f);
    }

    for(int i=0; i<4; ++i)
    {
      int net_j = 0;
      fread(&net_j, sizeof(int), 1, f);
      //std::cout << "Reading weights of layer " << net_j << std::endl;
      auto module = branch2->get(net_j-1);
      readParameters(state, (cunn::SpatialConvolutionMM*)module.get(), f);
    }

    int net_j = 0;
    fread(&net_j, sizeof(int), 1, f);
    readParameters(state, (cunn::Linear*)net->get(2).get(), f);
    fread(&net_j, sizeof(int), 1, f);
    readParameters(state, (cunn::Linear*)net->get(4).get(), f);
  }
  else if(net_type == "siam_desc")
  {
    net->add(std::make_shared<cunn::SpatialConvolutionMM>(state, 1,96, 7,7, 3,3));
    net->add(std::make_shared<cunn::ReLU>(state));
    net->add(std::make_shared<cunn::SpatialMaxPooling>(state, 2,2,2,2));
    net->add(std::make_shared<cunn::SpatialConvolutionMM>(state, 96,192, 5,5));
    net->add(std::make_shared<cunn::ReLU>(state));
    net->add(std::make_shared<cunn::SpatialMaxPooling>(state, 2,2,2,2));
    net->add(std::make_shared<cunn::SpatialConvolutionMM>(state, 192,256, 3,3));
    net->add(std::make_shared<cunn::ReLU>(state));
    net->add(std::make_shared<cunn::Reshape>(state, std::vector<size_t>(1,256)));

    for(int i=0; i<3; ++i)
    {
      int net_j = 0;
      fread(&net_j, sizeof(int), 1, f);
      //std::cout << "Reading weights of layer " << net_j << std::endl;
      auto module = net->get(net_j-2);
      readParameters(state, (cunn::SpatialConvolutionMM*)module.get(), f);
    }
  }
  else if(net_type == "siam")
  {
    std::vector<size_t> size = {1, 64, 64};
    cunn::Sequential::Ptr branch1 = std::make_shared<cunn::Sequential>(state);
    branch1->add(std::make_shared<cunn::Reshape>(state, size));
    branch1->add(std::make_shared<cunn::SpatialConvolutionMM>(state, 1,96, 7,7, 3,3));
    branch1->add(std::make_shared<cunn::ReLU>(state));
    branch1->add(std::make_shared<cunn::SpatialMaxPooling>(state, 2,2,2,2));
    branch1->add(std::make_shared<cunn::SpatialConvolutionMM>(state, 96,192, 5,5));
    branch1->add(std::make_shared<cunn::ReLU>(state));
    branch1->add(std::make_shared<cunn::SpatialMaxPooling>(state, 2,2,2,2));
    branch1->add(std::make_shared<cunn::SpatialConvolutionMM>(state, 192,256, 3,3));
    branch1->add(std::make_shared<cunn::ReLU>(state));
    branch1->add(std::make_shared<cunn::Reshape>(state, std::vector<size_t>(1,256)));

    cunn::Sequential::Ptr branch2 = std::make_shared<cunn::Sequential>(state);
    branch2->add(std::make_shared<cunn::Reshape>(state, size));
    branch2->add(std::make_shared<cunn::SpatialConvolutionMM>(state, 1,96, 7,7, 3,3));
    branch2->add(std::make_shared<cunn::ReLU>(state));
    branch2->add(std::make_shared<cunn::SpatialMaxPooling>(state, 2,2,2,2));
    branch2->add(std::make_shared<cunn::SpatialConvolutionMM>(state, 96,192, 5,5));
    branch2->add(std::make_shared<cunn::ReLU>(state));
    branch2->add(std::make_shared<cunn::SpatialMaxPooling>(state, 2,2,2,2));
    branch2->add(std::make_shared<cunn::SpatialConvolutionMM>(state, 192,256, 3,3));
    branch2->add(std::make_shared<cunn::ReLU>(state));
    branch2->add(std::make_shared<cunn::Reshape>(state, std::vector<size_t>(1,256)));

    cunn::Parallel::Ptr branches = std::make_shared<cunn::Parallel>(state, 1, 1);
    branches->add(branch1);
    branches->add(branch2);

    net->add(branches);
    net->add(std::make_shared<cunn::Linear>(state, 256*2, 256*2));
    net->add(std::make_shared<cunn::ReLU>(state));
    net->add(std::make_shared<cunn::Linear>(state, 256*2, 1));

    for(int i=0; i<3; ++i)
    {
      int net_j = 0;
      fread(&net_j, sizeof(int), 1, f);
      //std::cout << "Reading weights of layer " << net_j << std::endl;
      auto casted1 = (cunn::SpatialConvolutionMM*)(branch1->get(net_j-1)).get();
      auto casted2 = (cunn::SpatialConvolutionMM*)(branch2->get(net_j-1)).get();
      readParameters(state, casted1, f);
      THCudaTensor_copy(state, casted2->weight, casted1->weight);
      THCudaTensor_copy(state, casted2->bias, casted1->bias);
    }
    int net_j = 0;
    fread(&net_j, sizeof(int), 1, f);
    readParameters(state, (cunn::Linear*)net->get(1).get(), f);
    fread(&net_j, sizeof(int), 1, f);
    readParameters(state, (cunn::Linear*)net->get(3).get(), f);
  }
  else if(net_type == "siam_decision")
  {
    net->add(std::make_shared<cunn::Linear>(state, 256*2, 256*2));
    net->add(std::make_shared<cunn::ReLU>(state));
    net->add(std::make_shared<cunn::Linear>(state, 256*2, 1));

    int net_j = 0;
    fread(&net_j, sizeof(int), 1, f);
    readParameters(state, (cunn::Linear*)net->get(0).get(), f);
    fread(&net_j, sizeof(int), 1, f);
    readParameters(state, (cunn::Linear*)net->get(2).get(), f);
  }
  else if(net_type == "siam2stream")
  {
    cunn::Parallel::Ptr branches = std::make_shared<cunn::Parallel>(state, 1,1);
    auto branch1 = createSiam2streamBranch(state);
    auto branch2 = createSiam2streamBranch(state);
    branches->add(branch1);
    branches->add(branch2);
    auto sub_branch1_1 = (cunn::Sequential*)branch1->get(0).get();
    auto sub_branch1_2 = (cunn::Sequential*)branch1->get(1).get();
    auto sub_branch2_1 = (cunn::Sequential*)branch2->get(0).get();
    auto sub_branch2_2 = (cunn::Sequential*)branch2->get(1).get();

    net->add(branches);
    net->add(std::make_shared<cunn::Linear>(state, 1024, 1024));
    net->add(std::make_shared<cunn::ReLU>(state));
    net->add(std::make_shared<cunn::Linear>(state, 1024, 1));

    for(int i=0; i<4; ++i)
    {
      int net_j = 0;
      fread(&net_j, sizeof(int), 1, f);
      auto casted1 = (cunn::SpatialConvolutionMM*)(sub_branch1_1->get(net_j-1)).get();
      auto casted2 = (cunn::SpatialConvolutionMM*)(sub_branch2_1->get(net_j-1)).get();
      readParameters(state, casted1, f);
      THCudaTensor_copy(state, casted2->weight, casted1->weight);
      THCudaTensor_copy(state, casted2->bias, casted1->bias);
    }
    for(int i=0; i<4; ++i)
    {
      int net_j = 0;
      fread(&net_j, sizeof(int), 1, f);
      auto casted1 = (cunn::SpatialConvolutionMM*)(sub_branch1_2->get(net_j-1)).get();
      auto casted2 = (cunn::SpatialConvolutionMM*)(sub_branch2_2->get(net_j-1)).get();
      readParameters(state, casted1, f);
      THCudaTensor_copy(state, casted2->weight, casted1->weight);
      THCudaTensor_copy(state, casted2->bias, casted1->bias);
    }
    int net_j = 0;
    fread(&net_j, sizeof(int), 1, f);
    readParameters(state, (cunn::Linear*)net->get(1).get(), f);
    fread(&net_j, sizeof(int), 1, f);
    readParameters(state, (cunn::Linear*)net->get(3).get(), f);
  }
  else if(net_type == "siam2stream_decision")
  {
    net->add(std::make_shared<cunn::Linear>(state, 1024, 1024));
    net->add(std::make_shared<cunn::ReLU>(state));
    net->add(std::make_shared<cunn::Linear>(state, 1024, 1));

    int net_j = 0;
    fread(&net_j, sizeof(int), 1, f);
    readParameters(state, (cunn::Linear*)net->get(0).get(), f);
    fread(&net_j, sizeof(int), 1, f);
    readParameters(state, (cunn::Linear*)net->get(2).get(), f);
  }
  else if(net_type == "siam2stream_desc")
  {
    auto branch = createSiam2streamBranch(state);
    auto sub_branch1 = (cunn::Sequential*)branch->get(0).get();
    auto sub_branch2 = (cunn::Sequential*)branch->get(1).get();
    net->add(branch);

    for(int i=0; i<4; ++i)
    {
      int net_j = 0;
      fread(&net_j, sizeof(int), 1, f);
      auto casted = (cunn::SpatialConvolutionMM*)(sub_branch1->get(net_j-1)).get();
      readParameters(state, casted, f);
    }
    for(int i=0; i<4; ++i)
    {
      int net_j = 0;
      fread(&net_j, sizeof(int), 1, f);
      auto casted = (cunn::SpatialConvolutionMM*)(sub_branch2->get(net_j-1)).get();
      readParameters(state, casted, f);
    }
  }

  //std::cout << net->tostring() << std::endl;

  fclose(f);
  return net;
}

