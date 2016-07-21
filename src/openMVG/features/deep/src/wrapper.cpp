#include "loader.h"
#include "wrapper.h"

struct Impl {
  cunn::Sequential::Ptr net;
  THCState *state;
};

Network::Network(const char* filename) {
  THCState *state = (THCState*)malloc(sizeof(THCState));
  THCudaInit(state);

  ptr = std::make_shared<Impl>();
  ptr->state = state;
  init(filename);
}

Network::~Network() {
  THCState* state = ptr->state;
  ptr.reset();
  THCudaShutdown(state);
}

void Network::init(const char* filename) {
  ptr->net = loadNetwork(ptr->state, filename);
}

std::string Network::tostring() const {
  return ptr->net->tostring();
}

void Network::setDevice(int i) const {
  cudaSetDevice(i);
}

void Network::reset() {
  ptr->net.reset();
}

void forwardAnyNet(THCState *state, cunn::Module* net,
    THFloatTensor *input, THFloatTensor *output) {
  THLongStorage *input_size = THFloatTensor_newSizeOf(input);
  THCudaTensor *input_cuda = THCudaTensor_newWithSize(state, input_size, NULL);
  THCudaTensor_copyFloat(state, input_cuda, input);

  THCudaTensor* output_cuda = net->forward(input_cuda);

  THLongStorage *output_size = THCudaTensor_newSizeOf(state, output_cuda);
  THFloatTensor_resize(output, output_size, NULL);
  THFloatTensor_copyCuda(state, output, output_cuda);

  THCudaTensor_free(state, input_cuda);
  THLongStorage_free(input_size);
  THLongStorage_free(output_size);
}

void Network::forward(THFloatTensor *input, THFloatTensor *output) {
  forwardAnyNet(ptr->state, (cunn::Module*)ptr->net.get(), input, output);
}
