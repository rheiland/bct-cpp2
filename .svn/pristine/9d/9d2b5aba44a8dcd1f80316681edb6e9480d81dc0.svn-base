.PHONY: all clean

all: example_float example example_long_double

example_float: example_float.cpp
	$(CXX) $(CXXFLAGS) -lbct_float -lgsl -lgslcblas -o $@ $<

example: example.cpp
	$(CXX) $(CXXFLAGS) -lbct -lgsl -lgslcblas -o $@ $<

example_long_double: example_long_double.cpp
	$(CXX) $(CXXFLAGS) -lbct_long_double -lgsl -lgslcblas -o $@ $<

clean:
	-rm example_float example example_long_double
