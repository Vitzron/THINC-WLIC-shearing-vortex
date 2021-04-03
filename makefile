CXX = g++
CXXFLAGS = -std=c++11 -Wall -Ofast
LIBS = -I /usr/local/include/eigen3/Eigen 			# change the path to Eigen directory
PROM = main
OBJS = Grid.o Grid_Utilities.o Eigen_Utilities.o THINC_WLIC.o main.o

${PROM}: ${OBJS}
	${CXX} -o ${PROM} ${OBJS} ${LIBS}
clean:
	rm -f *.o *.dat ${PROM}