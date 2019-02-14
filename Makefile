CC=mpic++
all: Fractal.cpp
	$(CC) -o fractal_parallel Fractal.cpp

clean: 
	rm fractal_parallel mandelbrot_par.tga