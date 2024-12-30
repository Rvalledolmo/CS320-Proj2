CXX = g++
CXXFLAGS = -Wall -g 

all: cache-sim

cache-sim: cache.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

clean:
	rm -f cache-sim *output.txt