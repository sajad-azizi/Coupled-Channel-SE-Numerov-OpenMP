#include "./Potentials.hpp"

Potentials::Potentials(Parameters *parameters){//;//(int N_grid, int channels, double dr, int Nroots, int NTHREADS, int external_parameter){
    
    this->N_grid = parameters->N_grid;
    this->channels = parameters->channels;
    this->dr = parameters->dr;
    this->Nroots = parameters->Nroots;
    this->NTHREADS = parameters->NTHREADS;
    this->external_parameter = parameters->external_parameter;
    this->vphi_max = 2.0*M_PI;//(parameters->N_phi-1)*parameters->dphi;
    
    
    pot_component.resize(N_grid);
    for(int k = 0; k < N_grid; ++k)
        pot_component[k] = Eigen::MatrixXcd::Zero(channels, channels);
    
}

double Potentials::Potential(double r, double phi){

    
    int iAA = external_parameter;
    
    double x = r*cos(phi);
    double y = r*sin(phi);

    double v0 = -1.0/r;
    double v1 = -exp(-r*r/3.0);
    double v2 = -exp(-0.2*r)/r;
    double v3 = 0.5*r*r;
    double alp = 2.1325;
    double v4 = -(1.00+exp(-alp*r))/r;
    double v5 = 0.0;

    double v6= -exp(-(x*x/3.0+y*y/27.0));
    double v7  = -exp(-((x+3.)*(x+3.)+(y-2.)*(y-2.))/9.0)\
                    -exp(-((x-3.)*(x-3.)+(y+2.)*(y+2.))/9.0)\
                        -exp(-((x+1.)*(x+1.)+(y+3.)*(y+3.))/9.0);
                        
                        
    double w = 1.5;
    double v8 = 0.5*(x*x+w*w*y*y); //E = (nx+1/2)+w(ny+1/2)->E0 = 0.5+0.5*w

    //double a,b,c;
    //a=b=1;
   // c = 4.0;

    //double v9 = a*x*x+b*y*y+c*x*y;

    double R = 2.0;
    //cout <<"R: "<<R<<endl;
    double v10 = -1.0/sqrt((x-R/2.0)*(x-R/2.0)+y*y+0.6384)-1.0/sqrt((x+R/2.0)*(x+R/2.0)+(y)*(y)+0.6384);//-1.0/sqrt((y+R)*(y+R)+x*x+0.6384);

    double v11 = -1.0/sqrt((x)*(x)+(y-3)*(y-3)+0.6384)\
                        -1.0/sqrt((x)*(x)+(y+3.)*(y+3.0)+0.6384)\
                            -1.0/sqrt((x-3)*(x-3)+(y)*(y)+0.6384)\
                                -1.0/sqrt((x+3)*(x+3)+(y)*(y)+0.6384);
    double v12 = -1.0/((r+R/2.0))-1.0/((r-R))-1.0/((r+2.0*r*cos(phi)*R/2.0))-1.0/sqrt(r*r+r*cos(phi)*R+r*sin(phi)*R+0.001);
    
    double theta = M_PI/2.0;
    double VV0 = 1.0;
    
    double v13 = VV0*(r*r+1-3.0*(r*cos(phi)*sin(theta) + cos(theta)))/pow(r*r+1,5.0/2.0);
    
    double v14 = -1.0/sqrt((x-3)*(x-3)+(y-3)*(y-3)+0.6384)\
                        -1.0/sqrt((x+3)*(x+3)+(y+3.)*(y+3.0)+0.6384);
    
    double v15 = -1.0/sqrt((x)*(x)+y*y+0.6384);
    
    double v16 = -1.0/sqrt((x-R/2.0)*(x-R/2.0)+y*y+0.6384)-1.0/sqrt((x+R/2.0)*(x+R/2.0)+(y)*(y)+0.6384)-1.0/sqrt((y+R)*(y+R)+x*x+0.6384);
    
    
    double v17 = -exp( -0.3*( (x-R/2.0)*(x-R/2.0) + (y)*(y) ) )\
                    -exp( -0.3*( (x+R/2.0)*(x+R/2.0) + (y)*(y) ) );
    
    double v18 = -exp(-r*r/3);

    double v19 = -exp( -( (x-R/2.0)*(x-R/2.0) + (y)*(y) )/9. )\
                    -exp( -( (x+R/2.0)*(x+R/2.0) + (y)*(y) )/9. );


    double v20 = -exp( -0.3*( (x-R/2.0)*(x-R/2.0) + (y)*(y) ) )\
                    -exp( -0.3*( (x+5*R)*(x+5*R) + (y)*(y) ) )\
                          -exp( -0.3*( (x)*(x) + (y)*(y) ) );

    double b = -0;
    double a = 0;//-1.5;
    double c = 0;//-1.5;
    double v21 = -exp( -0.3*( (x-R/2.0)*(x-R/2.0) + (y)*(y) ) ) * ( (y/a)*(y/a) + (y/b)*(y/b) );
    double v22 = -0.7*exp(-0.1*(2*(x*x)+0.8*(y*y)+2.4*x*y));//-exp(-0.9*((x-a)*(x-a)+(y-b)*(y-b)));//-0.5*exp(-0.9*((x+c)*(x+c)+(y+b)*(y+b)));


    double vv0 = -1.0;
    double theta_ = R*M_PI;//1.6;// M_PI/8.0;
    double v23 = vv0*(x*x+y*y - 2*x*y + 3.0*( r*cos(phi)*sin(theta_) + cos(theta_) )*(r*cos(phi)*sin(theta_) + cos(theta_) ) )/pow(r*r+1.0,2.5);


    double a0 = 1.;
    double a1 = 1.;
    double a2 = 0.;
    double a3 = -2.0+R;

    double epsi = 0.5;

    double v25 = -1.0*exp(-(a0*x*x+a1*y*y+a2*x+a3*x*y))/pow( (r+epsi),2.5 );
    
    double v26 = -1.0/(r*r*r);

    double RR = 50, aa = 0.01, cc = 15.0;
    double v27 = - 1.0/(1.0 + exp((r - RR)/aa));
    
    
    double v28 = - 2.0*exp(-(r-cc)*(r-cc));
    
    double v29 = -exp( -( (x-R/2.0)*(x-R/2.0) + (y)*(y) )/3. )\
                    -exp( -( (x+R/2.0)*(x+R/2.0) + (y)*(y) )/3. );

    double x00 = 1.0;
    double v30 = -exp( -( (x-x00)*(x-x00) + (y)*(y) )/3.0 );
    
    double v31 = -exp( -( (x-x00)*(x-x00) + (y)*(y) )/3. ) -exp( -( (x-x00)*(x-x00) + (y-x00)*(y-x00) )/3. ) ;
    
    
    double v32 = -exp( -( (x-x00)*(x-x00) + (y)*(y) )/3. ) -exp( -( (x-x00)*(x-x00) + (y-x00)*(y-x00) )/3. )\
                    -exp(-r*r/3)-exp( -( (x)*(x) + (y+2*x00)*(y+2*x00) )/3. )-exp( -( (x+2*x00)*(x+2*x00) + (y-x00)*(y-x00) )/3. );
    
    return v18;

}




void Potentials::potMatrixElement(){
    
    double r_max = (N_grid-1)*dr;
    int l_max = int((channels-1)/2);
    
    std::vector<double> roots(Nroots,0.0);
    std::vector<double> weights(Nroots,0.0);
    std::ifstream fin("./roots_legendre_"+std::to_string(int(Nroots))+".dat");
    if(!fin.is_open()){std::cerr<<"requested file does not exist! :( \n"; exit(0);}
    for(int i = 0; i < Nroots; i++){
        fin >> roots[i] >> weights[i];
    }

    std::function<dcompx(double,double,int) > func;
    func = [=](double vphi, double r, int m){ return ( cos(m*vphi)+I*sin(m*vphi) )*Potential(r,vphi);};
    
    std::ifstream pot_in("../pot_component_"+std::to_string(int(r_max))+"_"+std::to_string(int(channels))+".dat");
    if(!pot_in.is_open()){
        cout << "pot_component has not already been saved, try calculating ...\n";
        // print potential
        /*std::ofstream p_out("potential.dat");
        for(int i = 0; i < N_grid; i+=200){
            for(int j = 0; j < N_phi; j++){
                p_out << i*dr <<"\t"<<j*dphi<< "\t" << i*dr*cos(j*dphi) <<"\t"<<i*dr*sin(j*dphi)<<"\t"<<Potential(i*dr,j*dphi)<<endl;
            }
        }
        p_out.close();*/

        std::vector<std::vector<dcompx> > dipolePot(channels*channels,std::vector<dcompx>(N_grid,0.0));

        int local_NTHREADS = int(NTHREADS);
        int CHUNK = int(N_grid/local_NTHREADS);
        omp_set_num_threads(local_NTHREADS);
        #pragma omp parallel for default(shared) schedule(static, CHUNK)
        for(int ir = 1; ir < N_grid; ir++){
            double r = ir*dr;
            //do nothong
            for(int m = 0; m < channels; m ++){
                int loc_m = -l_max+m;
                for(int k = 0; k < channels; k++){
                    int loc_k = -l_max+k;
                    
                    dcompx sum_phi = 0.0;
                    for(int j = 0; j < roots.size(); j++){
                        double tt = 0.5*(vphi_max-0.0) * roots[j] + (0.0 + vphi_max)/2.0;
                        sum_phi += 0.5*(vphi_max-0.0)*(func(tt,r,loc_m-loc_k))*(weights[j]);
                    }
                    dipolePot[m*channels+k][ir]=( (loc_m*loc_m-1.0/4.0)/(2.0*r*r) )*delta(loc_m,loc_k)+sum_phi/(2.0*M_PI);
                    //pot_component[ir](m,k)=( (loc_m*loc_m-1.0/4.0)/(2.0*r*r) )*delta(loc_m,loc_k)+sum_phi/(2.0*M_PI);
                }
            }
        }
        
        cout << "done 1\n";


        std::ofstream pot_out("pot_component_"+std::to_string(int(r_max))+"_"+std::to_string(int(channels))+".dat");
        //omp_set_num_threads(local_NTHREADS);
        //#pragma omp parallel for default(shared) schedule(static, CHUNK)
        for(int ir = 0; ir < N_grid; ir++){
            for(int m = 0; m < channels; m ++){
                // to lift the degenracy
                if(m > 0){
                    for(int k = 0; k < channels; k++){
                        pot_component[ir](m,k) = dipolePot[m*channels+k][ir] ;//+ 1e-8*(sin(ir*dr/r_max * M_PI))*(sin(ir*dr/r_max * M_PI));
                        pot_out << real(pot_component[ir](m,k)) <<"\t"<<imag(pot_component[ir](m,k)) << endl;
                    }
                }
                else{
                    for(int k = 0; k < channels; k++){
                        pot_component[ir](m,k) = dipolePot[m*channels+k][ir] ;//- 1e-10*(sin(ir*dr/r_max * M_PI))*(sin(ir*dr/r_max * M_PI));
                        pot_out << real(pot_component[ir](m,k)) <<"\t"<<imag(pot_component[ir](m,k)) << endl;
                    }
                }
            }
        }pot_out.close();
        
        cout << "done 2\n";
        
        dipolePot.clear();
    }
    else{
        cout << "pot_component has already been saved, try reading the existed file ...\n";
        double in_re, in_im;
        for(int ir = 0; ir < N_grid; ir++){
            for(int m = 0; m < channels; m ++){
                for(int k = 0; k < channels; k++){
                    pot_in >> in_re >> in_im;
                    pot_component[ir](m,k) = in_re + I*in_im;
                }
            }
        }
    }
    
    
    //cout << "\n pot: \n" << pot_component[100] << endl;
    

    cout << "done 3\n";
}


