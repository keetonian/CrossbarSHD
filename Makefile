CXX = g++
CFLAGS = -std=c++11 -O2 -Wall 
THREAD = -pthread -lpthread 

all: CrossbarSHD cxcompare

CrossbarSHD: crossbar_shd.cpp ctpl.h compare.o reference.o commandline.o
	$(CXX) $(CFLAGS) $(THREAD) -o CrossbarSHD crossbar_shd.cpp compare.o reference.o commandline.o ctpl.h

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
