#include <TH/TH.h>
#include <memory>

struct Impl;

struct Network {
  Network(const char* filename);
  ~Network();
  void init(const char* filename);
  std::string tostring() const;
  void reset();
  void setDevice(int i) const;
  void forward(THFloatTensor *input, THFloatTensor *output);

  std::shared_ptr<Impl> ptr;
};

