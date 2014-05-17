all: Makefile mtrand.cpp mtrand.h recombsim.cpp 
	g++ mtrand.cpp recombsim.cpp -o recombsim

clean:
	rm recombsim *.o *~
