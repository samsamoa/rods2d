# Basic Makefile with variables


#CC = /opt/intel/Compiler/11.0/074/bin/intel64/icpc
#CFLAGS = -I../Dynamics -O3 -axHost -ipo -opt-class-analysis -fp-model fast
CC = g++
CFLAGS = -I/Users/samsamoa/Work-Hagan/Dynamics/stable/src -Wall -O3

COMPILE = $(CC) $(CFLAGS) -c

all: rods2d clump

clump: clumps.o Rod.o
	cd build; $(CC) -o clump clumps.o Rod.o

fluctuate: densityFluctuations.o Rod.o
	cd build; $(CC) -o fluctuate densityFluctuations.o Rod.o

rods2d: main.o RodSimulation.o Rand.o Rod.o InteractionPair.o MotorGrid.o
	cd build; $(CC) -o bd main.o RodSimulation.o Rand.o Rod.o InteractionPair.o MotorGrid.o


clumps.o: clumps.cpp 
	$(COMPILE) -o build/clumps.o clumps.cpp

densityFluctuations.o: densityFluctuations.cpp 
	$(COMPILE) -o build/densityFluctuations.o densityFluctuations.cpp

main.o: main.cpp
	$(COMPILE) -o build/main.o main.cpp


RodSimulation.o: RodSimulation.cpp
	$(COMPILE) -o build/RodSimulation.o RodSimulation.cpp


Rod.o: Rod.cpp
	$(COMPILE) -o build/Rod.o Rod.cpp

Rand.o: ../../Dynamics/stable/src/Rand.cpp
	$(COMPILE) -o build/Rand.o ../../Dynamics/stable/src/Rand.cpp

InteractionPair.o: InteractionPair.cpp
	$(COMPILE) -o build/InteractionPair.o InteractionPair.cpp

MotorGrid.o: MotorGrid.cpp
	$(COMPILE) -o build/MotorGrid.o MotorGrid.cpp
