# Makefile for building the C language shared library for the CalcMatrixAvg demonstration package.
C = gcc
OPTS = -c -fPIC -I/usr/include/R 
LOADER = gcc

OBJECTS = src/image.false.nearest.o

 ImageFalseNearest.so: $(OBJECTS)
	R CMD SHLIB -o ImageFalseNearest.so $(OBJECTS) -lm

.c.o: ; $(C) $(OPTS) -c $<

clean:
	find -name "*~" -exec rm -rf {} \;
	rm -rf *.o *.so
