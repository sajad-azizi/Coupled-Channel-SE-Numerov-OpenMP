#include "./Equations.hpp"

Equations::Equations(Potentials *potentials, Parameters *parameters): potentials(potentials), parameters(parameters){
    
    this->N_grid = parameters->N_grid;
    this->channels = parameters->channels;
    this->dr = parameters->dr;
    this->l_max = int((parameters->channels-1)/2);
    this->p = parameters->p;
    
    Rinv_vector.resize(N_grid);
    for(int k = 0; k < N_grid; ++k)
        Rinv_vector[k] = Eigen::MatrixXcd::Zero(channels, channels);
    Rinv_vector_back.resize(N_grid);
    for(int k = 0; k < N_grid; ++k)
        Rinv_vector_back[k] = Eigen::MatrixXcd::Zero(channels, channels);
    Winv_vector.resize(N_grid);
    for(int k = 0; k < N_grid; ++k)
        Winv_vector[k] = Eigen::MatrixXcd::Zero(channels, channels);
    
    
     Rinv = Eigen::MatrixXcd::Zero(channels, channels);
     In = Eigen::MatrixXcd::Identity(channels, channels);
     R = Eigen::MatrixXcd::Zero(channels, channels);
     Wmat = Eigen::MatrixXcd::Zero(channels, channels);
     U = Eigen::MatrixXcd::Zero(channels, channels);
    
}



void Equations::proper_initialization_R(double Energy, Eigen::MatrixXcd &resRinv){
    
    
    double hh12 = dr*dr/12.0;
    
    
    
    // correction of m=0 for few step of R matrix
    // we use J0(kr) for the wavefunction
    dcompx k_wave = 0.0;
    
    for(int i = 1; i < p; i++){
        
        Wmat = In + hh12*(Energy*In - potentials->pot_component[i])*2.0;
        U = 12.0*Wmat.inverse() - 10.0*In;
        R = U - Rinv;
        
        double r = i*dr, rp1 = (i+1)*dr;
        
        double loc_energy = Energy - potentials->Potential(r,0);
        if( loc_energy < 0.) k_wave = I*sqrt( 2.0*abs(loc_energy) );
        else k_wave = sqrt( 2.0*abs(loc_energy) );
        
        dcompx F1  = Wmat(l_max,l_max) * imaginary_riccati_Jn(0,  k_wave*r)/sqrt(2.0*M_PI);
       /// cout << i << "\t" << F1 << "\t" << loc_energy << "\t" << k_wave << "\t" << imaginary_riccati_Jn(0,  k_wave*r) << "\t" << Wmat(l_max,l_max) << "\t" << l_max << endl;
        Wmat(l_max,l_max)  = In(l_max,l_max)  + hh12*(Energy*In(l_max,l_max)  - potentials->pot_component[i+1](l_max,l_max) )*2.0;
        
        loc_energy = Energy - potentials->Potential(rp1,0);
        if(loc_energy < 0.) k_wave = I*sqrt( 2.0*abs(loc_energy) );
        else k_wave = sqrt( 2.0*abs(loc_energy) );
        
        dcompx F2  = Wmat(l_max,l_max) * imaginary_riccati_Jn(0,  k_wave*rp1)/sqrt(2.0*M_PI);
        R(l_max,l_max) = F2/F1;
        
        
        Rinv = R.inverse();
        // reserve them for calculating the eigenfunctions
        if(1){
            Rinv_vector[i] = Rinv;
            Winv_vector[i] = Wmat.inverse();
        }
        
    }
    resRinv = Rinv;
    
    
   // cout << resRinv << endl;
    
    
    
}

std::pair<int, double> Equations::OutwardNodeCounting(double Energy){
    
    double hh12 = dr*dr/12.0;
    
    
    proper_initialization_R(Energy, Rinv);

    //std::ofstream fout("R_foreward.dat");
    
    int node_c = 0;
    double node_pos = 0;

    for(int i = p; i < N_grid; i++){
        
        Wmat = In + hh12*(Energy*In - potentials->pot_component[i])*2.0;
        U = 12.0*Wmat.inverse() - 10.0*In;
        R = U - Rinv;
        
        Rinv = R.inverse();
        es.compute(R, Eigen::EigenvaluesOnly);
        es1.compute(Wmat, Eigen::EigenvaluesOnly);
        //fout << real(R.determinant()) << endl;
        //if(real(R.determinant()) > 0.0){
        for(int m = 0; m < channels; m++){
            if((es1.eigenvalues().real()(m)) > 0.0){
                //cout << "requarments are  full filled!\n";
                if((es.eigenvalues().real()(m)) < 0.0){
                    node_c ++;
                    node_pos = i*dr;
                    //cout <<"eigVal: " << m<< "\t" << real(R(m,m)) << "\t" << node_c<<endl;
                }
            }
        }
    }

    return std::make_pair(node_c, node_pos);

}

void Equations::propagateBackward(double Energy, int i_match, Eigen::MatrixXcd &resRmp1, bool save){
    
    
    double hh12 = dr*dr/12.0;
    
    //std::ofstream fout("R_backward.dat");

    for(int i = N_grid-1-1; i > i_match; i--){
        
        Wmat = In + hh12*(Energy*In - potentials->pot_component[i])*2.0;
        U = 12.0*Wmat.inverse() - 10.0*In;
        R = U - Rinv;
        
        Rinv = R.inverse();
        
        if(save){
            
            Rinv_vector[i] = Rinv;
            Winv_vector[i] = Wmat.inverse();
            //if(i==i_match+1)cout<<Rinv_vector[i]<<endl;
        }
        Rinv_vector_back[i] = Rinv;
        //fout<<i*dr<<"\t"<<R.determinant().real()<<"\t"<<R.determinant().imag()<<"\t"<<abs(R.determinant())<<"\t"<<R.determinant().real()<<"\t"<<Rinv.determinant().imag()<<"\t"<<abs(Rinv.determinant())<<endl;
    }//fout.close();
    
    resRmp1 = R;
}

void Equations::propagateForward(double Energy, int i_match, Eigen::MatrixXcd &resRm, bool save){
    
    
    double hh12 = dr*dr/12.0;
    
    proper_initialization_R(Energy,Rinv);
    
    
    for(int i = p; i <= i_match; i++){
        
        Wmat = In + hh12*( Energy*In - potentials->pot_component[i] )*2.0;
        U = 12.0*Wmat.inverse() - 10.0*In;
        R = U - Rinv;
        
        Rinv = R.inverse();
        
        if(save){
            Rinv_vector[i] = Rinv;
            Winv_vector[i] = Wmat.inverse();
        }
    }
    
    resRm = R;
}


//imaginary argument of bessel function
dcompx Equations::imaginary_riccati_Jn(int n, dcompx x){
    
    if(real(x) == 0.){
    
        if(n < 0){
            return pow(-1,abs(n))*I*sqrt(M_PI*abs(x)/2.0)*exp(I*M_PI*double(n)/2.0)*std::cyl_bessel_i(abs(n), abs(x));
        }
        else{
            return I*sqrt(M_PI*abs(x)/2.0)*exp(I*M_PI*double(n)/2.0)*std::cyl_bessel_i(abs(n), abs(x));
        }
    }
    else if(imag(x) == 0.){
        if(n<0)return pow(-1,abs(n))*sqrt(M_PI*abs(x)/2.0)*std::cyl_bessel_j(abs(n), abs(x));
        else return sqrt(M_PI*abs(x)/2.0)*std::cyl_bessel_j(n, abs(x));
    }
    else{
        std::cerr << "argument of bessel function is complex, it is not pure imag or real \n";
        exit(0);
        return 0.; // to shut up compiler warnings :(
    }
    
}






