CXX = g++

# Compiler flags
CXXFLAGS = -Wall -Wextra -std=c++11 -O3

# Libraries
LIBS = -lm

# Source files
SRC = main.cpp

# Executable name
EXEC = a.out

# Python script for visualization
VIS_SCRIPT = visualize.py

# Target for compiling the .cpp file
all: $(EXEC)

$(EXEC): $(SRC)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

# Target for running the Python visualization script
visualize:
	python3 $(VIS_SCRIPT)

# Target for cleaning up generated files
clean:
	rm -f $(EXEC)

.PHONY: all clean visualize