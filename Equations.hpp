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

#include "./Potentials.hpp"

class Equations{
    
    public:
        explicit Equations(Potentials *potentials, Parameters *parameters);
        
        void propagateBackward(double Energy, int i_match, Eigen::MatrixXcd &resRmp1, bool save);
        void propagateForward(double Energy, int i_match, Eigen::MatrixXcd &resRm, bool save);
        void proper_initialization_R(double Energy, Eigen::MatrixXcd &resRinv);
        std::pair<int, double> OutwardNodeCounting(double Energy);
        dcompx imaginary_riccati_Jn(int n, dcompx x);
        
        // saved vectors
        std::vector<Eigen::MatrixXcd> Rinv_vector;
        std::vector<Eigen::MatrixXcd> Rinv_vector_back;
        std::vector<Eigen::MatrixXcd> Winv_vector;
        
    private:
        
        Potentials *potentials;
        Parameters *parameters;

        int N_grid;
        int channels;
        double dr;
        int l_max;// = int((channels-1)/2);
        int i_match;
        double Energy;
        
        //the numer of correction points `p' of the ratio matrix R
        int p; // min value = 1
        
        //renormalized Numerov Method functions
        Eigen::MatrixXcd Rinv;
        Eigen::MatrixXcd In;
        Eigen::MatrixXcd R;
        Eigen::MatrixXcd Wmat;
        Eigen::MatrixXcd U;
        
        //diagonalization 
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es;
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es1;
        
};
