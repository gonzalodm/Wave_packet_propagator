CC := nvcc
CFLAGS := -O3 --gpu-architecture=sm_50

objects    = particle.o wp_subs.o wp_cuda.o wp_main.o
executable = wp_program.out

$(executable): $(objects)
			$(CC) -o $@ $^ ${CFLAGS}

particle.o: particle.cpp particle.h
			$(CC) -o $@ -c $< ${CFLAGS}

wp_subs.o: wp_subs.cpp wp_subs.h
			$(CC) -o $@ -c $< ${CFLAGS}

wp_cuda.o: wp_cuda.cu wp_cuda.h
			$(CC) -o $@ -c $< ${CFLAGS}

wp_main.o: wp_main.cpp
			$(CC) -o $@ -c $< ${CFLAGS}

.PHONY: clean
clean:
	      rm $(objects) $(executable)
