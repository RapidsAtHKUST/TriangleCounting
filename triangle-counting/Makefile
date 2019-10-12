CXX = g++
RM = rm
CXXFLAGS = -O3 -fopenmp -lpthread -std=c++11 -march=native

.PHONY: clean

tc: log.o search_util.o triangle_counting_main.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^

search_util.o: util/search_util.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $^

log.o: util/log.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $^

clean: 
	$(RM) -rf *.o tc
