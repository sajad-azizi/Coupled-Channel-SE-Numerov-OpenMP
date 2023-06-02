#include "./DipoleMat.hpp"

DipoleMat::DipoleMat(Wavefunctions *wavefunctions, Parameters *parameters):wavefunctions(wavefunctions),parameters(parameters){
    
    this->N_grid = parameters->N_grid;
    this->channels = parameters->channels;
    this->dr = parameters->dr;
    this->l_max = int((parameters->channels-1)/2);
    this->NTHREADS = parameters->NTHREADS;
    this->N_phi = parameters->N_phi;
    this->dphi = parameters->dphi;
    this->vphi_max = 2.0*M_PI;//(parameters->N_phi-1)*parameters->dphi;
    this->Nroots = parameters->Nroots;
    
}


double DipoleMat::dipole_function(double phi){
    return cos(phi);
}


//================calculate scattering state====================================
void DipoleMat::calculate_complex_dipole_matrix_element_ingoingBC( Eigen::MatrixXcd &A, Eigen::MatrixXcd &B, const double energy){

    
    std::ofstream oout("real_dipole.dat");
    std::vector<dcompx> d_real(channels,0.0);
    real_dipole_matrix_element(d_real);
    for(int beta = 0; beta < channels; beta++){
        oout << -l_max + beta <<"\t"<<real(d_real[beta])<<"\t"<<imag(d_real[beta])<<endl;
    }
    oout.close();
    cout <<"real dipole done!\n";
    
    
    double k_wave = sqrt(2.0*energy);
    
    Eigen::MatrixXcd AiB_inv = ( (A-I*B).inverse() ).transpose().conjugate();
    
    std::vector<dcompx> d_complex(channels,0.0);
    
    std::ofstream reduced_dipole_out("reduced_dipole.dat");
    for(int m = 0; m < channels; m++){
        int loc_m = -l_max+m;
        dcompx sum_beta = 0.0;
        for(int beta = 0; beta < channels; beta++){
            sum_beta += AiB_inv(m,beta)*d_real[beta];
        }
        d_complex[m] = sum_beta;
        reduced_dipole_out<<loc_m<<"\t"<<real(d_complex[m])<<"\t"<<imag(d_complex[m])<<"\t"<<abs(d_complex[m])*abs(d_complex[m])<<"\t"<<std::arg(d_complex[m])<<endl;
    }
    reduced_dipole_out.close();
    
    cout << " complex matrix element is done \n";
    
    
    d_real.clear();
    
    
    std::vector<Eigen::VectorXcd> d_S;
    d_S.resize(N_phi);
    for(int k = 0; k < N_phi; ++k)
        d_S[k] = Eigen::VectorXcd::Zero(channels);
    
    
    for(int i = 0; i < N_phi; i++){
        double phi = i*dphi;
        for(int m = 0; m < channels; m++){
            int loc_m = -l_max+m;
            d_S[i](m)= 2.0*std::pow(-I,loc_m)/sqrt(k_wave)*d_complex[m]*exp(I*double(loc_m)*phi);
        }
    }

}
void DipoleMat::real_dipole_matrix_element(std::vector<dcompx> &d_real){
    
    std::vector<double> roots(Nroots,0.0);
    std::vector<double> weights(Nroots,0.0);
    std::ifstream fin("/data/finite/sazizi/new_non_adiabatic/forThesisi/roots_legendre_"+std::to_string(Nroots)+".dat");
    if(!fin.is_open()){std::cerr<<"requested file does not exist! :( \n"; exit(0);}
    for(int i = 0; i < Nroots; i++){
        fin >> roots[i] >> weights[i];
    }
    
    
    std::function<dcompx(double,int) > func;
    func = [=](double vphi,int m){ return ( cos(m*vphi)+I*sin(m*vphi) )*dipole_function(vphi);};
    
    Eigen::MatrixXcd angular_integration = Eigen::MatrixXcd::Zero(channels,channels);
    for(int m = 0; m < channels; m ++){
        int loc_m = -l_max+m;
        for(int k = 0; k < channels; k++){
            int loc_k = -l_max+k;
            double sum_phi_re = 0.0;
            double sum_phi_im = 0.0;
            #pragma omp parallel for default(shared) reduction(+:sum_phi_re,sum_phi_im)
            for(int j = 0; j < roots.size(); j++){
                double tt = 0.5*(vphi_max-0.0) * roots[j] + (0.0 + vphi_max)/2.0;
                sum_phi_re += 0.5*(vphi_max-0.0)*( cos((loc_k-loc_m)*tt)*dipole_function(tt) )*(weights[j])*1.0/(2.0*M_PI);
                sum_phi_im += 0.5*(vphi_max-0.0)*( sin((loc_k-loc_m)*tt)*dipole_function(tt) )*(weights[j])*1.0/(2.0*M_PI);
            }
            angular_integration(m,k)= sum_phi_re+I*sum_phi_im;
            //cout << loc_m <<"\t"<<loc_k<<"\t"<<angular_integration(m,k)<<endl;
        }
    }
    
    
    std::function<dcompx(double,int, int, int) > radial_func;
    radial_func = [=](double r,int m, int k, int beta){
        int ir = int(r/dr);
        return  angular_integration(m,k)*dr*std::conj( wavefunctions->scattering_eigenfunc[ir](m,beta) )*r*wavefunctions->eigfunc[ir](k);
    };
    
    
    for(int beta = 0; beta < channels; beta++){
        dcompx realDipol_sum = 0.0;
        //do nothong
        for(int m = 0; m < channels; m ++){
            int loc_m = -l_max+m;
            for(int k = 0; k < channels; k++){
                int loc_k = -l_max+k;
                
                double sum_radial_re = 0.0;
                double sum_radial_im = 0.0;
                #pragma omp parallel for default(shared) reduction(+:sum_radial_re,sum_radial_im)
                for(int ir = 1; ir < N_grid-1; ir++){
                    double r = ir*dr;
                    if(ir % 2 == 0){
                        sum_radial_re += 2.0*real( dr*std::conj( wavefunctions->scattering_eigenfunc[ir](m,beta) )*r*wavefunctions->eigfunc[ir](k) );
                        sum_radial_im += 2.0*imag( dr*std::conj( wavefunctions->scattering_eigenfunc[ir](m,beta) )*r*wavefunctions->eigfunc[ir](k) );
                    }
                    else{
                        sum_radial_re += 4.0*real( dr*std::conj( wavefunctions->scattering_eigenfunc[ir](m,beta) )*r*wavefunctions->eigfunc[ir](k) );
                        sum_radial_im += 4.0*imag( dr*std::conj( wavefunctions->scattering_eigenfunc[ir](m,beta) )*r*wavefunctions->eigfunc[ir](k) );
                    }
                }
                realDipol_sum += angular_integration(m,k)*(\
                0.0+dr*std::conj(wavefunctions->scattering_eigenfunc[N_grid-1](m,beta))*(N_grid-1)*dr*wavefunctions->eigfunc[N_grid-1](k)+sum_radial_re + I*sum_radial_im )/3.0;
            }
        }
        
        d_real[beta] = realDipol_sum;
    }
    
    angular_integration.resize(0,0);
    roots.clear();
    weights.clear();

}
