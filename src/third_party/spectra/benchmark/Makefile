CXX = g++
CXXFLAGS = -Wall -O2
CPPFLAGS = -I../include -I.
FC = gfortran
FFLAGS = -O3
LDFLAGS =
LIBS = -larpack -llapack -lblas -lgfortran

HEADERS = $(wildcard ../include/Spectra/MatOp/*.h) \
	$(wildcard ../include/Spectra/LinAlg/*.h) \
	$(wildcard ../include/Spectra/Util/*.h) \
	$(wildcard ../include/Spectra/*.h)
OBJS = F77.o Cpp.o wrapper.o

.PHONY: all clean

all: benchmark.out

benchmark.out: $(HEADERS) $(OBJS) main.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) main.cpp $(OBJS) -o benchmark.out $(LDFLAGS) $(LIBS)

%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	-rm *.out *.o

