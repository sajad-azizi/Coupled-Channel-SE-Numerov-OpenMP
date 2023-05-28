#pragma once

#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <fstream>
#include <iomanip>
#include <functional>
#include <algorithm>
#include <complex>
#include <omp.h>

using std::cout;
using std::endl;
typedef std::complex<double> dcompx;

constexpr dcompx I{0.0, 1.0};
constexpr dcompx zero{0.0,0.0};

#include "/data/finite/sazizi/new_non_adiabatic/forThesisi/Eigen/Dense"
#include "./Parameters.hpp"

class Potentials{
    
    public:
        explicit Potentials(Parameters *parameters);//int N_grid, int channels, double dr, int Nroots,  int NTHREADS, int external_parameter);
        
        double Potential(double r, double phi);
        void potMatrixElement();
        inline int delta(int k, int l){return k == l;}
        
        
        //angle part
        //double vphi_max = 2.0*M_PI;
        
        
        std::vector<Eigen::MatrixXcd> pot_component;
        
    private:
        
        Parameters *parameters;
        
        int external_parameter;
        
        int channels;
        int N_grid;
        double dr;
        
        //the number of OpenMP threads
        int NTHREADS;
        
        // number of roots of Gauss-Legender integration
        int Nroots;
        double vphi_max;//angle part
};
