CXX=g++
CXXFLAGS= -Wall -pedantic -O3
LDLIBS =  -lblas -fopenmp -lboost_program_options## -O0 -g -o perf

#Files to be included
OBJS 	 = main.o ShallowWater.o 
HEADERS  = ShallowWater.h
TARGET   = main

#Setting out test cases
ARGS1=--dt 0.1 --T 80.0 --Nx 100 --Ny 100 --ic 1 
ARGS2=--dt 0.1 --T 80.0 --Nx 100 --Ny 100 --ic 2 
ARGS3=--dt 0.1 --T 80.0 --Nx 100 --Ny 100 --ic 3
ARGS4=--dt 0.1 --T 80.0 --Nx 100 --Ny 100 --ic 4

# Number of Threads used
np = 9

# 1 for Loop Based Calculations, 2 for BLAS Based Calculations
calc = 1

default: $(TARGET)
all: $(TARGET)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $< $(LDLIBS)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS) 


test1: $(TARGET)
	./$(TARGET) $(ARGS1) --np $(np) --calc $(calc)

test2: $(TARGET)
	./$(TARGET) $(ARGS2) --np $(np) --calc $(calc)

test3: $(TARGET)
	./$(TARGET) $(ARGS3) --np $(np) --calc $(calc)

test4: $(TARGET)
	./$(TARGET) $(ARGS4) --np $(np) --calc $(calc)

#removing files created in directory
clean:
	rm -rf $(TARGET) *.o

##to plot run command gnuplot myplot.gnu