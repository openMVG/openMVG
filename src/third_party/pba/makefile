#To use AVX instead of SSE, you need to add -DCPUPBA_USE_AVX to CFLAGS
#################################################################

CUDA_INSTALL_PATH = /usr/local/cuda
CUDA_LIB_PATH = $(CUDA_INSTALL_PATH)/lib64
ifdef CUDA_BIN_PATH
	NVCC = $(CUDA_BIN_PATH)/nvcc
else
	NVCC = $(CUDA_INSTALL_PATH)/bin/nvcc
endif
ifndef CUDA_INC_PATH
	CUDA_INC_PATH = $(CUDA_INSTALL_PATH)/include
endif
#################################################################

# detect OS
OSUPPER = $(shell uname -s 2>/dev/null | tr [:lower:] [:upper:])
OSLOWER = $(shell uname -s 2>/dev/null | tr [:upper:] [:lower:])
DARWIN = $(strip $(findstring DARWIN, $(OSUPPER)))


SHELL = /bin/sh
BIN_DIR = ./bin
SRC_PBA = ./src/pba
SRC_DRIVER = ./src/driver
OUT_PBA = ./bin/out

CC = g++
CFLAGS =  -I$(CUDA_INC_PATH) -L$(CUDA_LIB_PATH) -fPIC -L/usr/lib64 -L/usr/lib -L$(BIN_DIR) -O2 -Wall -Wno-deprecated -pthread  -march=native -mfpmath=sse -fpermissive
NVCC_FLAGS = -I$(CUDA_INC_PATH) -O2 -Xcompiler -fPIC


# siftgpu header files
_HEADER_PBA = pba.h ConfigBA.h CuTexImage.h DataInterface.h ProgramCU.h SparseBundleCU.h
_HEADER_PBA_LIB = pba.h 
_OBJ_PBA = pba.o CuTexImage.o ConfigBA.o SparseBundleCU.o SparseBundleCPU.o

all: makepath pba driver 
 
#the dependencies of SiftGPU library 
DEPS_PBA = $(patsubst %, $(SRC_PBA)/%, $(_HEADER_PBA))

#rules for the rest of the object files
$(OUT_PBA)/%.o: $(SRC_PBA)/%.cpp $(DEPS_PBA) 
	$(CC) -o $@ $< $(CFLAGS) -c 

#build rule for CUDA 
$(OUT_PBA)/ProgramCU.o: $(SRC_PBA)/ProgramCU.cu $(DEPS_PBA)
	$(NVCC) $(NVCC_FLAGS) -o $@ $< -c
_OBJ_PBA += ProgramCU.o

OBJ_PBA  = $(patsubst %,$(OUT_PBA)/%,$(_OBJ_PBA))
LIBS_DRIVER = $(BIN_DIR)/libpba.a $(LIBS_PBA) 

pba: makepath $(OBJ_PBA)
	ar rcs $(BIN_DIR)/libpba.a $(OBJ_PBA)
	$(CC) -o $(BIN_DIR)/libpba.so $(OBJ_PBA) $(LIBS_PBA) $(CFLAGS) -lcudart -shared -fPIC
 
driver: makepath 
	$(CC) -o $(BIN_DIR)/driver $(SRC_DRIVER)/driver.cpp -lpba -lcudart $(CFLAGS) 
	
makepath:
	mkdir -p $(OUT_PBA)
	mkdir -p $(BIN_DIR) 
 
clean:
	rm -f $(OUT_PBA)/*.o
	rm -f $(BIN_DIR)/libpba.a
	rm -f $(BIN_DIR)/libpba.so
	rm -f $(BIN_DIR)/driver


