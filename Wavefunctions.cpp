#include "./Wavefunctions.hpp"

Wavefunctions::Wavefunctions(Equations *equations, Parameters *parameters): equations(equations){
    
    this->N_grid = parameters->N_grid;
    this->channels = parameters->channels;
    this->dr = parameters->dr;
    this->l_max = int((parameters->channels-1)/2);
    this->NTHREADS = parameters->NTHREADS;
    this->divide = parameters->divide;
    
    eigfunc.resize(N_grid);
    for(int k = 0; k < N_grid; ++k)
        eigfunc[k] = Eigen::VectorXcd::Zero(channels); 
    
    scattering_eigenfunc.resize(N_grid);
    for(int k = 0; k < N_grid; ++k)
        scattering_eigenfunc[k] = Eigen::MatrixXcd::Zero(channels,channels);
    
}




void Wavefunctions::calculate_final_continumm_states(Eigen::MatrixXcd &A, Eigen::MatrixXcd &B, const double energy){
    

    Eigen::MatrixXcd AiB_inv = ( (A-I*B).inverse() );
    
    

    std::vector<Eigen::MatrixXcd> scattering_eigenfunc_energyNorm;//scattering eigenfunction
    scattering_eigenfunc_energyNorm.resize(N_grid);
    for(int k = 0; k < N_grid; ++k)
        scattering_eigenfunc_energyNorm[k] = Eigen::MatrixXcd::Zero(channels,channels);
    
    
    for(int k = 0; k < N_grid; ++k)
        scattering_eigenfunc_energyNorm[k] = scattering_eigenfunc[k] * AiB_inv;
    
    std::ofstream mout;
    for(int m = 0; m < channels; m++){
        for(int beta = 0; beta < channels; beta++){
            mout.open("final_psi_"+std::to_string(m)+"_"+std::to_string(beta)+".dat");
            for(int k = 0; k < int(N_grid/divide); k+=10){
                mout<<real( scattering_eigenfunc_energyNorm[k](m,beta) )<<"\t"<<imag( scattering_eigenfunc_energyNorm[k](m,beta) )<<endl;
            }
            mout.close();
        }
    }
    
}


void Wavefunctions::calculate_eigenfunction(double E, int i_match){
    
    std::cout << std::fixed << std::setprecision(14);
    std::cerr << std::fixed << std::setprecision(14);
    //std::cerr << "(calculate_eigenfunction) Eigenvalue: " << E << "i match: "<< i_match << std::endl;
    
    Eigen::MatrixXcd Rm, Rmp1;
    
    equations->propagateForward(E,i_match,Rm,true);
    equations->propagateBackward(E,i_match,Rmp1,true);
    
    
    Eigen::MatrixXcd diff = Rmp1.inverse()-Rm;
    
    //cout<<Rinv_vector[i_match-1]<<endl;
    Eigen::VectorXd Veig = Eigen::VectorXd::Zero(channels);
    es.compute(diff);
    for(int m = 0; m < channels; m++){
        Veig(m) = abs( es.eigenvalues().real()(m));
    }
    //cout << Veig << endl;
    Eigen::VectorXd::Index jm,jrm;
    double minV_eig = Veig.minCoeff(&jm);
    
    //cout<<"(index,min,max) of pot at inf: " << jm << "\t" << minV_eig<< "\t" << Veig.maxCoeff() << endl;
    //cout << es.eigenvectors().col(jm)<< endl; 
    //cout << "eigenvalues diff: " << endl << Veig << endl;
    //cout << "min eigenvectors: " << endl << es.eigenvectors().col(jm) << endl;
    //cout << " eigenvectors: " << endl << es.eigenvectors()<< endl;
    
    
    Eigen::FullPivLU<Eigen::MatrixXcd> lu(diff);
    lu = lu.setThreshold(1.0e-2);
    Eigen::MatrixXcd ker = lu.kernel();
    //cout << "kernel: " << endl << ker << endl;
    
    

    Eigen::VectorXcd fn_prev, fn, fn_next, psin;
    Eigen::VectorXcd fn_matching_point;
    
    
    fn_matching_point =  es.eigenvectors().col(jm);//ker
    //cout << "ker:" << ker << endl;
    fn_next = fn_matching_point;
    fn_prev = fn_matching_point;

    
    
    for(int k = i_match - 1; k > 0; --k,fn_next = fn){
        //cout << "k: " << k << "; Rinv: " << endl << Rinv_vector[k] << endl;
        fn = equations->Rinv_vector[k] * fn_next;
        psin = equations->Winv_vector[k] * fn;
        eigfunc[k] = psin;

        //cout << "eigfunc[" << k << "]: " << endl << eigfunc[k] << endl;
        //fn_next = fn;
    }
    
    for(int k = i_match + 1; k < N_grid-1; ++k, fn_prev = fn){
        //cout << "k: " << k << "; Rinv: " << endl << Rinv_vector[k] << endl;
        fn = equations->Rinv_vector[k] * fn_prev;
        psin = equations->Winv_vector[k] * fn;
        eigfunc[k] = psin;
        //cout << Rinv_vector[k] << endl;
        //cout << "eigfunc[" << k << "]: " << endl << eigfunc[k] << endl;	
    }

    psin = equations->Winv_vector[i_match] * fn_matching_point;
    eigfunc[i_match] = psin;

    cout << "eigfunc[" << i_match << "]: " << std::endl << eigfunc[i_match] <<endl;

    //return eigfunc;
}



void Wavefunctions::calculate_eigenfunction_continuum(double E){
    
    std::cout << std::fixed << std::setprecision(14);
    std::cerr << std::fixed << std::setprecision(14);
    //std::cerr << "(calculate_eigenfunction) Eigenvalue: " << E << "i match: "<< i_match << std::endl;
    
    Eigen::MatrixXcd Rm, Rmp1;
    
    equations->propagateForward(E,N_grid-1,Rm,true);
    
    Eigen::MatrixXcd diff = Rm;
    
    //cout<<Rinv_vector[i_match-1]<<endl;
    Eigen::VectorXd Veig = Eigen::VectorXd::Zero(channels);
    es.compute(diff);
    for(int m = 0; m < channels; m++){
        Veig(m) = abs( es.eigenvalues().real()(m));
    }
    //cout << Veig << endl;
    Eigen::VectorXd::Index jm,jrm;
    double minV_eig = Veig.minCoeff(&jm);
    
    //cout<<"(index,min,max) of pot at inf: " << jm << "\t" << minV_eig<< "\t" << Veig.maxCoeff() << endl;
    //cout << es.eigenvectors().col(jm)<< endl; 
    //cout << "eigenvalues diff: " << endl << Veig << endl;
    //cout << "min eigenvectors: " << endl << es.eigenvectors().col(jm) << endl;
    //cout << " eigenvectors: " << endl << es.eigenvectors()<< endl;
    
    
    Eigen::FullPivLU<Eigen::MatrixXcd> lu(diff);
    lu = lu.setThreshold(1.0e-2);
    Eigen::MatrixXcd ker = lu.kernel();
    //cout << "kernel: " << endl << ker << endl;
    
    

    Eigen::VectorXcd fn_prev, fn, fn_next, psin;
    Eigen::VectorXcd fn_matching_point;
    
    
    fn_matching_point =  es.eigenvectors().col(jm);//ker
    //cout << "ker:" << ker << endl;
    fn_next = fn_matching_point;
    fn_prev = fn_matching_point;

    
    
    for(int k = N_grid-1 - 1; k > 0; --k,fn_next = fn){
        //cout << "k: " << k << "; Rinv: " << endl << Rinv_vector[k] << endl;
        fn = equations->Rinv_vector[k] * fn_next;
        psin = equations->Winv_vector[k] * fn;
        eigfunc[k] = psin;

        //cout << "eigfunc[" << k << "]: " << endl << eigfunc[k] << endl;
        //fn_next = fn;
    }

    psin = equations->Winv_vector[N_grid-1] * fn_matching_point;
    eigfunc[N_grid-1] = psin;

    cout << "eigfunc[" << N_grid-1 << "]: " << std::endl << eigfunc[N_grid-1] <<endl;

    //return eigfunc;
}


void Wavefunctions::Normalization(std::vector<Eigen::VectorXcd> &eigfunc){
    
    double sum = 0.;
    for(int l = 0; l < channels; l++){
        for(int i = 0; i < N_grid; i++){
            sum += abs(eigfunc[i](l))*abs(eigfunc[i](l))*dr;
        }
    }
    sum = sqrt(sum);
    for(int l = 0; l < channels; l++){
        if(sum!=0){
            for(int i = 0; i < N_grid; i++){
                eigfunc[i](l) = eigfunc[i](l)/sum;;
            }
        }
    }
}



//wavefunction for each open channel (j)
void Wavefunctions::calculate_channel_wavefunction(const double Energy){
    
    //std::cout << std::fixed << std::setprecision(14);
    //std::cerr << std::fixed << std::setprecision(14);
    //std::cerr << "(calculate_eigenfunction) Eigenvalue: " << Energy << std::endl;

    Eigen::MatrixXcd Rm;
    equations->propagateForward(Energy,N_grid-1,Rm,true);
    
    //cout << "R->inf \n" << Rm << endl;

    
    std::vector<Eigen::VectorXcd> eigfunc_(N_grid);
    for(int k = 0; k < N_grid; ++k)
        eigfunc_[k] = Eigen::VectorXcd::Zero(channels);

    Eigen::VectorXcd fn_matching_point(channels);
    
    Eigen::VectorXcd fn_prev, fn, fn_next, psin;

    //cout << "win : \n" <<  equations->Winv_vector[N_grid-1].inverse() << endl;
    //cout << "win_row j : \n" << equations->Winv_vector[N_grid-1].inverse().row(0) << endl;

    for(int j = 0; j < channels; j++){

        fn_matching_point = equations->Winv_vector[N_grid-1].inverse().row(j); //or .col(j)
        //far from the origin the potential is vanished and just centrifugal barrier remained and that is diagonal
        //therefore each row (column) represeneds as an open channel which are indepenednt from each other

        fn_next = fn_matching_point;

        //std::ofstream fout("data.txt");
        for(int k = N_grid-1 - 1; k > 0; --k, fn_next = fn){
            //cout << "k: " << k << "; Rinv: " << endl << equations->Rinv_vector[k] << endl;
            fn = equations->Rinv_vector[k] * fn_next;
            psin = equations->Winv_vector[k] * fn;
            eigfunc_[k] = psin;
            //fout << k*dr << "\t" << real(eigfunc_[k](l_max))<< "\t" << imag(eigfunc_[k](l_max)) << endl;	
        }/*
        
        
        if(abs(1.0 - Rm(j,j)) > 0.2 ){
            cout << "1.0 - Rm(j,j) " << 1.0 - Rm(j,j) << endl;
            for(int  ir = 0; ir < N_grid; ir++){
                for(int m = 0; m < channels; m++){
                    scattering_eigenfunc_save[ir](m,j) = eigfunc_[ir](m);
                }
            }
        }*/
        
        
        //each column has an indepenednt solution labeled by j
        for(int  ir = 0; ir < N_grid; ir++){
            for(int m = 0; m < channels; m++){
                scattering_eigenfunc[ir](m,j) = eigfunc_[ir](m);
            }
        }
        //cout << "psi at inf :" << j << endl << eigfunc_[N_grid-2] << endl;

    }
    //cout << N_grid <<"\t" << channels << endl;
    //cout << "psi at inf :"<< endl << eigfunc_[N_grid-2] << endl;
    
    eigfunc_.clear();

}







void Wavefunctions::calculate_A_B_matrices( Eigen::MatrixXcd &A, Eigen::MatrixXcd &B, const double energy){
    
    double k = sqrt(2.0*energy);
    Eigen::MatrixXcd In = Eigen::MatrixXcd::Identity(channels, channels);
    
    
    //method 1
    std::vector<double> psi2d_real(N_grid,0.0);
    std::vector<double> psi2d_imag(N_grid,0.0);
    for(int j = 0; j < channels; j++){
        for(int m = 0; m < channels; m++){
            int loc_m = -l_max+m;
            
            for(int i = 0; i < N_grid; i++){
                double r = i*dr;
                psi2d_real[i] = real( scattering_eigenfunc[i](m,j) );
                psi2d_imag[i] = imag( scattering_eigenfunc[i](m,j) );
            }
            
            double a = 0.,b = 0.,c=0,d=0;
            std::tie(a,b) = getCoefficents(psi2d_real,k,loc_m);
            std::tie(c,d) = getCoefficents(psi2d_imag,k,loc_m);
            
            A(m,j) = a+I*c;
            B(m,j) = b+I*d;
            //cout << -l_max+j <<"\t"<<loc_m<<"\t"<<a<<", "<<b<<endl;
        }
    }
    

    //cout << "\n A_j: \n" << A <<endl;
    //cout << "\n B_j: \n" << B <<endl;
    
    /*Eigen::MatrixXcd K = B*A.inverse();
    cout << "\n K_re :\n"<<K<<endl;
    Eigen::MatrixXcd S = (In+I*K).inverse()*(In-I*K);
    //cout <<"\n S:\n"<<S<<endl;
    
    Eigen::MatrixXcd unitary = S.conjugate().transpose()*S;
    for(int j = 0; j < channels; j++){
        dcompx sum_uni = 0.0;
        for(int m = 0; m < channels; m++){
            sum_uni += unitary(j,m);
        }
        cout<<-l_max+j<<" (~1) "<<abs(sum_uni)*abs(sum_uni)<<" (~0) "<<1.0-abs(sum_uni)*abs(sum_uni)<<" abs_S^2: "<<abs(S(j,j))*abs(S(j,j))<<" diag S-Mat/2: "<<std::arg(S(j,j))/2.0<<endl;
    }*/
    
}


std::pair<double, double> Wavefunctions::getCoefficents(std::vector<double> &data, double k, int loc_m){
    
        int i_save = N_grid-6;
        std::vector<double> wf_f;

        
        /*std::vector<double> bessel_zeros;
        getBessel_zeroes(bessel_zeros,k,loc_m);

        //double x_alp = 4.0*finding_max_min(data) - 3.0*r_max;
        //i_save =bessel_zeros[bessel_zeros.size()-4];// int(x_alp/dr);
        int i_count = 2;
        while(1){
            i_save = bessel_zeros[bessel_zeros.size()-i_count];
            if(i_save < N_grid - 4){
            break;
            }
            i_count ++;
        }
        cout << "i_save: "<<i_save<<"\t"<< "size: "<<N_grid-1 - i_save<<"\n";
        bessel_zeros.clear();*/
        
        double zeros_bessel;
        int oo=1;
        while(1){
            zeros_bessel  = gsl_sf_bessel_zero_Jnu(abs(loc_m),oo);
            if(zeros_bessel > k*(N_grid-1)*dr){
                break;
            }
            oo++;
        }
        int i_count = 2;
        while(1){
            //i_save = bessel_zeros[bessel_zeros.size()-i_count];
	        i_save = int(gsl_sf_bessel_zero_Jnu(abs(loc_m),oo-i_count)/dr);
            if(i_save < (N_grid - 4) ){
                break;
            }
            i_count ++;
        }
	    //i_save = int(gsl_sf_bessel_zero_Jnu(abs(loc_m),oo-2)/dr);
	    //cout << "i_save: "<<i_save<<"\t"<< "size: "<<N_grid-1 - i_save<<"\n";

        //int i_save_r = i_save;
        //cout << "i_save: "<<i_save << " corss_x : " << i_save*dr  << endl;
        for(int i = i_save; i < N_grid-5; i++){
            wf_f.push_back(data[i]);
        }
        
        
        
        /*double sumA = 0,sumB = 0,sumC = 0,sumD = 0,sumE = 0;
        for(int i = 0; i < wf_f.size(); i++){
            double r_loc = (i_save+i)*dr;
            sumA += wf_f[i]*riccati_Jn(loc_m,k*r_loc);
            sumB += riccati_Jn(loc_m,k*r_loc)*riccati_Jn(loc_m,k*r_loc);
            sumC += riccati_Yn(loc_m,k*r_loc)*riccati_Jn(loc_m,k*r_loc);
            sumE += wf_f[i]*riccati_Yn(loc_m,k*r_loc);
            sumD += riccati_Yn(loc_m,k*r_loc)*riccati_Yn(loc_m,k*r_loc);
        }
        */
        double sumA = 0,sumB = 0,sumC = 0,sumD = 0,sumE = 0;
        #pragma omp parallel for default(shared) reduction(+:sumA)
        for(int i = 0; i < wf_f.size(); i++){
            double r_loc = (i_save+i)*dr;
            sumA += wf_f[i]*riccati_Jn(loc_m,k*r_loc);
        }
        #pragma omp parallel for default(shared) reduction(+:sumB)
        for(int i = 0; i < wf_f.size(); i++){
            double r_loc = (i_save+i)*dr;
            sumB += riccati_Jn(loc_m,k*r_loc)*riccati_Jn(loc_m,k*r_loc);
        }
        #pragma omp parallel for default(shared) reduction(+:sumC)
        for(int i = 0; i < wf_f.size(); i++){
            double r_loc = (i_save+i)*dr;
            sumC += riccati_Yn(loc_m,k*r_loc)*riccati_Jn(loc_m,k*r_loc);
        }
        #pragma omp parallel for default(shared) reduction(+:sumE)
        for(int i = 0; i < wf_f.size(); i++){
            double r_loc = (i_save+i)*dr;
            sumE += wf_f[i]*riccati_Yn(loc_m,k*r_loc);
        }
        #pragma omp parallel for default(shared) reduction(+:sumD)
        for(int i = 0; i < wf_f.size(); i++){
            double r_loc = (i_save+i)*dr;
            sumD += riccati_Yn(loc_m,k*r_loc)*riccati_Yn(loc_m,k*r_loc);
        }
        
        //cout << "coeffs: " << sumA << "\t" <<sumB<<"\t"<< sumC << "\t" <<sumE << "\t" <<sumD << endl;
        double determ;
        if(wf_f.size() == 0){cout<<" Achtung!! size is zero!!\n";exit(0);}
        determ = sumB*sumD-sumC*sumC;
        double a_ = (sumA*sumD-sumE*sumC)/determ;
        double b_ = (-sumA*sumC+sumE*sumB)/determ;

        //cout << "m: " << loc_m <<" a: "<<a_<<" b: "<<b_  <<" determ:"<<determ<<"\n";
        
        wf_f.clear();
        
        return std::make_pair(a_,b_);
}




void Wavefunctions::getBessel_zeroes(std::vector<double> &bessel_zeros, double k, int m){
    
    int zeros_count = 0;
    for(int i = 1; i < N_grid; i++){
        if( riccati_Yn(m, k*i*dr) != copysign(riccati_Yn(m, k*i*dr),riccati_Yn(m, k*(i-1)*dr))  ){
            bessel_zeros.push_back( i   );
            //cout << zeros_count << "\t" << i*dr << endl;
            //zeros_count++;
        }
    }
}
double Wavefunctions::riccati_Jn(int n,double x){
    if(n<0)return pow(-1,abs(n))*sqrt(M_PI*x/2.0)*std::cyl_bessel_j(abs(n), x);
    else return sqrt(M_PI*x/2.0)*std::cyl_bessel_j(n, x);
}
double Wavefunctions::riccati_Yn(int n,double x){
    if(n<0)return pow(-1,abs(n))*sqrt(M_PI*x/2.0)*std::cyl_neumann(abs(n), x);
    else return sqrt(M_PI*x/2.0)*std::cyl_neumann(n, x);
}
double Wavefunctions::spherical_jn(int n,double x){
    if(n<0)return pow(-1,abs(n))*std::sph_bessel(abs(n), x);
    else return std::sph_bessel(n, x);
}
double Wavefunctions::spherical_yn(int n,double x){
    if(n<0)return pow(-1,abs(n))*std::sph_neumann(abs(n), x);
    else return std::sph_neumann(n, x);
}
