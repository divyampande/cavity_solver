# Makefile for Lid-Driven Cavity Solver

# Compiler and flags
FC = gfortran
FFLAGS = -O3 -Wall -std=f2008
DEBUG_FLAGS = -g -fcheck=all -fbacktrace

# Directories
SRC_DIR = src
BUILD_DIR = build
RESULTS_DIR = results

# Source files (order matters for module dependencies)
SOURCES = $(SRC_DIR)/input.f90 \
          $(SRC_DIR)/boundary.f90 \
          $(SRC_DIR)/solver.f90 \
          $(SRC_DIR)/main.f90

# Object files
OBJECTS = $(patsubst $(SRC_DIR)/%.f90,$(BUILD_DIR)/%.o,$(SOURCES))

# Executable
TARGET = cavity_solver

# Default target
all: $(BUILD_DIR) $(RESULTS_DIR) $(TARGET)

# Create build directory
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Create results directory
$(RESULTS_DIR):
	mkdir -p $(RESULTS_DIR)

# Link executable
$(TARGET): $(OBJECTS)
	$(FC) $(FFLAGS) -o $@ $^
	@echo "Build successful: $(TARGET)"

# Compile object files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f90 | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -c $< -o $@

# Debug build
debug: FFLAGS = $(DEBUG_FLAGS)
debug: clean all

# Run simulation
INPUT ?= 
run: $(TARGET) | $(RESULTS_DIR)
	./$(TARGET) $(INPUT)

# Clean build artifacts
clean:
	rm -rf $(BUILD_DIR) $(TARGET)

# Clean everything including results
cleanall: clean
	rm -rf $(RESULTS_DIR)/*.dat

# Install dependencies (for reference)
.PHONY: help
help:
	@echo "Lid-Driven Cavity Solver"
	@echo "========================"
	@echo "Targets:"
	@echo "  make          - Build the solver"
	@echo "  make run      - Build and run simulation"
	@echo "  make debug    - Build with debug flags"
	@echo "  make clean    - Remove build artifacts"
	@echo "  make cleanall - Remove build artifacts and results"
	@echo ""
	@echo "Requirements:"
	@echo "  - gfortran compiler"
	@echo "  - Python 3 with numpy, matplotlib for visualization"

.PHONY: all clean cleanall debug run help
