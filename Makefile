CXX=g++
CXXDEF=
CXXFLAGS=-O3 


DEPS=       help.h system.h pairing.h infinite.h
OBJ= main.o help.o system.o pairing.o infinite.o

%.o: %.cpp $(DEPS)
	$(CXX) ${CXXDEF} -c -o $@ $< $(CXXFLAGS)
test:$(OBJ)
	$(CXX)  -o $@ $^ $(CXXFLAGS)
.PHONY: clean

clean:
	rm *.o  test
