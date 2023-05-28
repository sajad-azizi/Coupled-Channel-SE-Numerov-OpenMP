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

#include "/data/finite/sazizi/new_non_adiabatic/forThesisi/Eigen/Dense"
#include <gsl/gsl_sf_bessel.h>

#include "./Equations.hpp"


class Wavefunctions{
    
    public:
        
        explicit Wavefunctions(Equations *equations, Parameters *parameters);
        
        void calculate_eigenfunction(const double Energy, const int i_match);
        void calculate_eigenfunction_continuum(const double E);
        void Normalization(std::vector<Eigen::VectorXcd> &eigfunc);
        void calculate_channel_wavefunction(const double Energy);
        
        void calculate_final_continumm_states(Eigen::MatrixXcd &A, Eigen::MatrixXcd &B, const double energy);
        void calculate_A_B_matrices( Eigen::MatrixXcd &A, Eigen::MatrixXcd &B, const double Energy);
        std::pair<double, double> getCoefficents(std::vector<double> &data, double k, int loc_m);
        
        void getBessel_zeroes(std::vector<double> &bessel_zeros, double k, int m);
        double riccati_Jn(int n,double x);
        double riccati_Yn(int n,double x);
        double spherical_jn(int n,double x);
        double spherical_yn(int n,double x);
        
        
        std::vector<Eigen::VectorXcd> eigfunc;//bound state eigenfunction
        std::vector<Eigen::MatrixXcd> scattering_eigenfunc;//scattering eigenfunction
        
        
    private:
        
        Equations *equations;
        Parameters *parameters;
        double Energy;
        int i_match;
        int N_grid;
        int channels;
        int l_max;
        double dr;
        int NTHREADS;
        int divide;
        
        //diagonalization 
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es;
};
