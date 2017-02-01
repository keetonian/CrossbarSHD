CXX = g++
CFLAGS = -std=c++11 -O2 -Wall 
THREAD = -pthread -lpthread 

all: CrossbarSHD cxcompare

CrossbarSHD: crossbar_shd.cpp ctpl.h
	$(CXX) $(CFLAGS) $(THREAD) -o CrossbarSHD crossbar_shd.cpp ctpl.h

cxcompare: crossbar_compare.cpp
	$(CXX) $(CFLAGS) crossbar_compare.cpp -o cxcompare

clean:
	rm CrossbarSHD cxcompare
