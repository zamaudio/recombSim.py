all: Makefile mtrand.cpp mtrand.h recombsim.cpp 
	g++ mtrand.cpp recombsim.cpp -lbpp-core -lbpp-seq -g -o recombsim

clean:
	rm recombsim *.o *~
