# This file is part of recombsim

all: Makefile mtrand.cpp mtrand.h recombsim.cpp 
	g++ -O3 mtrand.cpp recombsim.cpp -lbpp-core -lbpp-seq -g -o recombsim

clean:
	rm recombsim *~
