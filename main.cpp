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
#define dmax(a,b) ((a>b)?a:b)
#define dmin(a,b) ((a<b)?a:b)
#define pi_ 3.1415926535897932384626433

typedef std::complex<double> dcompx;


#include "/data/finite/sazizi/new_non_adiabatic/forThesisi/Eigen/Dense"
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>


#include "./Potentials.hpp"
#include "./Equations.hpp"
#include "./Eigenvalues.hpp"
#include "./Wavefunctions.hpp"
#include "./DipoleMat.hpp"

//#include "./Parameters.hpp"

int main(){

    //space mech
    double r_max = 100.0;
    double dr = 0.01;
    int N_grid = int(r_max/dr)+1;
    
    //angle part
    double vphi_max = 2.0*M_PI;
    double phi_m = -M_PI;
    double dphi = 0.001;
    int N_phi = int((vphi_max)/dphi)+1;
    
    int l_max = 2; //angular momentum
    int channels = 2*l_max+1;
    
    //energy space     
    double Emin = -4.7;
    double Emax = 2.0;
    
    double Energy = 1.0;
    
    // number of roots of Gauss-Legender integration
    int Nroots = 1000;
    int external_parameter = 0;

    //the number of OpenMP threads
    int NTHREADS = omp_get_max_threads();
    
    //border between interaction and asymptotic region
    int divide = 8;
    
    //the numer of correction points `p' of the ratio matrix R
    int p = 9; // min value = 1
    
    
    cout << "Nthreads: " << NTHREADS << endl;
    std::cout << std::fixed << std::setprecision(14);
    std::cerr << std::fixed << std::setprecision(14);
    
    
    Parameters parameters(N_grid, channels, dr, Nroots, NTHREADS, N_phi, dphi, Emin,  Emax, divide, p, external_parameter);
    
    
    Potentials potentials(&parameters);// N_grid, channels, dr, Nroots, NTHREADS, external_parameter );
    
    potentials.potMatrixElement();
    cout << "pot done! \n";
    
    // renormalized functions
   Equations equations(&potentials, &parameters);
    //int node_counting;
    //double node_pos;
   // std::tie(node_counting, node_pos) = equations.OutwardNodeCounting(-0.35947853840436);
    
    // ground state finder
    Eigenvalues eigenvalues(&equations, &parameters);
    
    eigenvalues.groundstate_finder();
    /*for(int m = 0; m < channels; m++){
        std::ofstream pot_out("pot_anglar_mom_"+std::to_string(m)+".dat");
        for(int ir = 1; ir  < N_grid; ir++){
            pot_out << ir*dr <<"\t"<<  real(potentials.pot_component[ir](m,m)) <<"\t" << imag(potentials.pot_component[ir](m,m)) << endl;
        }
        pot_out.close();
    }*/
    
    /*cout  << "i_match " << eigenvalues.i_match << "\t" << eigenvalues.gsEnergy << endl;
    Wavefunctions wavefunctions(&equations, &parameters);
    wavefunctions.calculate_eigenfunction(eigenvalues.gsEnergy, eigenvalues.i_match);
    wavefunctions.Normalization(wavefunctions.eigfunc);
    
    
    std::ofstream bout;
    for(int m = 0; m < channels; m++){
        bout.open("basissets/psi_"+std::to_string(m)+".dat");
        for(int k = 0; k < N_grid; ++k){
            bout<<k*dr<<"\t"<<real(wavefunctions.eigfunc[k](m))<<"\t"<<imag(wavefunctions.eigfunc[k](m))<<std::endl;
        }
        bout.close();
    }
    
    Eigen::MatrixXcd In = Eigen::MatrixXcd::Identity(channels,channels);
    Eigen::MatrixXcd A = Eigen::MatrixXcd::Zero(channels,channels);
    Eigen::MatrixXcd B = Eigen::MatrixXcd::Zero(channels,channels);
    Eigen::MatrixXcd K = Eigen::MatrixXcd::Zero(channels,channels);
    
    
    wavefunctions.calculate_channel_wavefunction(Energy);
    wavefunctions.calculate_A_B_matrices(A,B,Energy);
    cout <<"\n A:\n"<<A<<endl;
    //cout <<"\n A.inverse():\n"<<A.inverse()<<endl;
    K = B*A.inverse();
    Eigen::MatrixXcd S = (In+I*K)*(In-I*K).inverse();
   // cout <<"\n S:\n"<<S<<endl;
    
    Eigen::MatrixXcd unitary = S.conjugate().transpose()*S;
    std::ofstream Sout("sk_matrix.dat");
    for(int j = 0; j < channels; j++){
        dcompx sum_uni = 0.0;
        for(int m = 0; m < channels; m++){
            sum_uni += unitary(j,m);
            Sout<<real(S(j,m))<<"\t"<<imag(S(j,m))<<"\t"<<real(K(j,m))<<"\t"<<imag(K(j,m))<<endl;;
        }
        cout<<-l_max+j<<" (~1) "<<abs(sum_uni)*abs(sum_uni)<<" (~0) "<<1.0-abs(sum_uni)*abs(sum_uni)<<" diag abs_S^2: "<<abs(S(j,j))*abs(S(j,j))<<" diag S-Mat/2: "<<std::arg(S(j,j))/2.0<<endl;
    }
    Sout.close();
    
    DipoleMat dipoleMat(&wavefunctions, &parameters);
    dipoleMat.calculate_complex_dipole_matrix_element_ingoingBC(A, B,Energy);
    
    cout << "printing continumm final wavefunction ...\n";
    //wavefunctions.calculate_final_continumm_states(A, B, Energy);
    cout << "running is done.\n"<< endl;
    
    */
    
    return 0;
}
