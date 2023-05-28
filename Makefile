.PHONY: all clean

CCXX := icpc
CCXXFLAGS := -std=c++17 -Wall -O2 -O3 -DNDEBUG -funroll-loops -ffast-math -lstdc++ 

OTHERFLAGS := -lm -O2  -no-multibyte-chars -fopenmp -std=c++17

EIGEN := /data/finite/sazizi/new_non_adiabatic/forThesisi/Eigen/Dense

GSL := -lgsl -lgslcblas

OPENMP := -fopenmp

BUILDDIR := ./build/

$(shell mkdir -p $(BUILDDIR))
$(shell mkdir -p basissets)

all: main 

main: $(addprefix $(BUILDDIR), $(patsubst %.cpp, %.o, $(wildcard *.cpp)))
	$(CCXX) $(CXXFLAGS) $^ -o $@ $(GSL) $(OPENMP)

$(BUILDDIR)%.o: %.cpp
	@echo ">> (g++) Compiling $<...";
	@$(CCXX) $(CXXFLAGS) -c -MD $< -o $@ $(OTHERFLAGS)


clean:
	@rm -r $(BUILDDIR)
	@rm -f main
#######icpc -Wall -O2 -std=c++17 -O3 -DNDEBUG -funroll-loops -ffast-math -lstdc++ -fopenmp -lgsl -lgslcblas 2d_coupled_tise.cpp -lm -no-multibyte-chars 

