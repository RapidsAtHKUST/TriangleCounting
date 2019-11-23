CXX = g++
RM = rm
CXXFLAGS = -O3 -fopenmp -lpthread -std=c++14 -march=native

.PHONY: clean

tc: log.o triangle_counting_main.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^

log.o: util/log.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $^

clean: 
	$(RM) -rf *.o tc
