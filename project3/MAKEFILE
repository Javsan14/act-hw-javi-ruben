# Compiler
FC = gfortran

# Compiler flags
FFLAGS = -Wall -Wextra -O2

# Source file
SRC = src/main.f90

# Executable name
EXE = md

# Default target
all: $(EXE)

# Rule to compile the Fortran source file into an executable
$(EXE): $(SRC)
	$(FC) $(FFLAGS) -o $(EXE) $(SRC)

# Clean up build files
clean:
	rm -f $(EXE) *.mod *.o

# Phony targets (not actual files)
.PHONY: all clean