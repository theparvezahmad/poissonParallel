TARGET := openmp_cpp.x
BUILD_DIR :=
SRC_DIR := ../src/

SRCS := $(SRC_DIR)main_openmp.cpp

OBJS := ${SRCS:$(SRC_DIR)%.cpp=$(BUILD_DIR)%.o}

FC := icc -diag-disable=10441 -qopenmp#fortran compiler
FFLAGS :=# -CB -traceback# -xCORE-AVX2# -qopenmp #-CB -traceback
# FFLAGS :=-qmkl -g -check all -fpe0 -warn -traceback -debug extended
# FC := gfortran
# FFLAGS := -g# -fopenmp #-CB -traceback
# LIBS := -llapack -lblas
# LIBS= -L/usr/local/lib/ -llapack -lblas

$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS) $(LIBS)

$(BUILD_DIR)%.o: $(SRC_DIR)%.cpp #wildcard rule, creation of *.o depends on *.f90
	$(FC) $(FFLAGS) -o $@ -c $< $(LIBS)

.PHONY: clean
clean:
	@rm -f -v $(BUILD_DIR)*.mod *.o *.x
