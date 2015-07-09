CXX=g++
CXXDEF=
CXXFLAGS=-O3 


DEPS=help.h
OBJ= main.o help.o

%.o: %.cpp $(DEPS)
	$(CXX) ${CXXDEF} -c -o $@ $< $(CXXFLAGS)
test:$(OBJ)
	$(CXX)  -o $@ $^ $(CXXFLAGS)
.PHONY: clean

clean:
	rm *.o  test
