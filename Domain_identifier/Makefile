CC = mpicc

#CFLAGS = -c -O3 -I/home/$(USER)/local/include/ -I/usr/include/ -I/usr/include/hdf5/serial/
CFLAGS = -c -g -I/home/$(USER)/local/include/ -I/usr/include/ -I/usr/include/hdf5/serial/

LFLAGS = -lm -L/home/$(USER)/local/lib -Wl,"-R /home/$(USER)/local/lib"  -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lgsl -lgslcblas -lhdf5

PROGRAM = Domain_identifier


$(PROGRAM): 
	$(CC) $(CFLAGS) $@.c -o $@.o
	$(CC) $@.o $(LFLAGS) -o $@

clean:
	rm -rf $(PROGRAM)
	rm -rf *~
	rm -rf *.out
	rm -rf *#
	rm -rf *.o
	rm -rf *.a	
	rm -rf *.so

plot:
	gnuplot script.gpl

# run:
# 	$(PROGRAM)
