CXX = g++
NVCC = nvcc -ccbin g++
CFLAGS = -std=c++11 -O2 

all: CrossbarSHD cxcompare

CrossbarSHD: crossbar_shd.cu compare.o reference.o commandline.o
	$(NVCC) $(CFLAGS) -o CrossbarSHD crossbar_shd.cu compare.o reference.o commandline.o 

compare.o: compare.h compare.cpp
	$(CXX) $(CFLAGS) -c compare.cpp 

reference.o: reference.h reference.cpp
	$(CXX) $(CFLAGS) -c reference.cpp

commandline.o: commandline.h commandline.cpp
	$(CXX) $(CFLAGS) -c commandline.cpp

cxcompare: crossbar_compare.cpp
	$(CXX) $(CFLAGS) crossbar_compare.cpp -o cxcompare

clean:
	rm CrossbarSHD cxcompare *.o
