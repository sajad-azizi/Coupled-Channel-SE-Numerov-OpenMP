#pragma once
#include <iostream>

using std::cout;
using std::endl;

#include "./Wavefunctions.hpp"

class DipoleMat{
    
    public:
        explicit DipoleMat(Wavefunctions *wavefunctions, Parameters *parameters);
        void calculate_complex_dipole_matrix_element_ingoingBC(Eigen::MatrixXcd &A, Eigen::MatrixXcd &B, const double energy);
        double dipole_function(const double phi);
        void real_dipole_matrix_element(std::vector<dcompx> &d_real);
        
        //double vphi_max = 2.0*M_PI;
        
    private:
        
        Wavefunctions *wavefunctions;
        Parameters *parameters;
        
        int N_grid;
        int channels;
        int l_max;
        double dr;
        int NTHREADS;
        int N_phi;
        double dphi;
        double vphi_max;
        int Nroots;
    
};
