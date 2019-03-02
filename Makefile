CPP = mpicc
CPPFLAGS = -std=c++0x -O3 -D _MPI

qvoter_parallel: main.o network_sim.o network_sim_parallel.o
	${CPP} ${CPPFLAGS} main.o network_sim.o network_sim_parallel.o -o qvoter_parallel

network_sim_parallel.o: network_sim_parallel.cpp
	${CPP} ${CPPFLAGS} network_sim_parallel.cpp -c

network_sim.o: network_sim.cpp
	${CPP} ${CPPFLAGS} network_sim.cpp -c

main.o: main.cpp
	${CPP} ${CPPFLAGS} main.cpp -c
	
clean:
	rm -rf *.o
	rm qvoter_parallel
	clear