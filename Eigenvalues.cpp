#include "./Eigenvalues.hpp"

Eigenvalues::Eigenvalues(Equations *equations, Parameters *parameters  ): equations(equations){
    
    this->Emin = parameters->Emin;
    this->Emax = parameters->Emax;
    this->N_grid = parameters->N_grid;
    
}


void Eigenvalues::groundstate_finder(){ //bisection algorithm
    
    double e_L = Emin;
    double e_H = Emax;
    
    //int i_match;
    
    int desire_node = 1;
     
    std::ofstream eout("ground_energy.dat");eout << std::fixed << std::setprecision(15);

    int oo = 0;

    //e_L = Emin;
    e_H = Emax;
    double e_tr = e_L;
    int node_counting = 0;
    double node_pos = 0;
    while(1){
        std::tie(node_counting, node_pos) = equations->OutwardNodeCounting(e_tr);
        if(node_counting < desire_node){
            e_L = e_tr;
        }
        else{
            e_H = e_tr;
        }
        
        cout  << e_tr << "\t" << e_L << "\t" << e_H << "\t" << node_pos << "\t" << node_counting << "\t" << desire_node << endl;
        if(abs(e_H-e_L) < 1e-12){
            
            i_match = matching_point_finder(e_tr);
            gsEnergy = e_tr;
            cout<<"done ! energy: "<<e_tr<<" i_m_:"<<i_match<<" nod pos: "<<node_pos<<" count node: "<<node_counting<<" desire node: "<<desire_node<<endl;
            eout << e_tr <<"\t" << i_match << endl;
            break;
        }
        e_tr = (e_L + e_H)/2.0;
    }
    e_L = e_tr;
    eout.close();

}




void Eigenvalues::Boundstates_finder(){
            
    double e_L = Emin;
    double e_H = Emax;
    int i_match;
    
    int desire_node = 7;
    
    std::ofstream psout; 
    std::ofstream eout("eigenvalues.dat");eout << std::fixed << std::setprecision(15);
    
    
    int node_max = 0;
    double node_pos_max = 0;
    std::tie(node_max, node_pos_max) = equations->OutwardNodeCounting(Emax);
    
    
    int node_min = 0;
    double node_pos_min = 0;
    std::tie(node_min, node_pos_min) = equations->OutwardNodeCounting(Emin);
    
    
    int oo = 0;
    for(desire_node = 1; desire_node < node_max; desire_node++){
        //e_L = Emin;
        e_H = Emax;
        double e_tr = e_L;
        int node_counting = 0;
        double node_pos = 0;
        while(1){
            std::tie(node_counting, node_pos) = equations->OutwardNodeCounting(e_tr);
            if(node_counting < desire_node){
                e_L = e_tr;
            }
            else{
                e_H = e_tr;
            }
            
            //cout  << e_tr << "\t" << e_L << "\t" << e_H << "\t" << node_pos << "\t" << node_counting << "\t" << desire_node << endl;
            if(abs(e_H-e_L) < 1e-12){
                
                if(e_tr < 0){
                    i_match = matching_point_finder(e_tr);
                }
                else{
                    i_match = N_grid-1;
                }
                
                cout << "done ! e: " << e_tr << "\t" << i_match << "\t" << node_pos << "\t" << node_counting << "\t" << desire_node << endl;
                eout << oo <<"\t"<< e_tr <<"\t" << i_match << endl;
                oo++;
                break;
            }
            e_tr = (e_L + e_H)/2.0;
        }
        e_L = e_tr;
    }
    
    eout.close();
}


int Eigenvalues::matching_point_finder(double Energy){
    
    int i_match;
    
    Eigen::MatrixXcd Rm, Rmp1;
    equations->propagateForward(Energy,N_grid-1,Rm,true);
    equations->propagateBackward(Energy,0,Rmp1,false);
    
    double per_deter = 0.0;
    //std::ofstream ffout("matching_point.dat");
    for(int i = 1; i < N_grid-2; i++){
        dcompx deter = (equations->Rinv_vector[i].inverse() - equations->Rinv_vector_back[i]).determinant();
        if(abs(deter) > per_deter && i > 5){
            i_match = i;
            cout << "i_match : " << i_match << endl;
            break;
        }
        //ffout << i*dr << "\t" << abs(deter) << "\t" << abs((Rinv_vector[i].inverse()).determinant()) <<"\t" << abs(( Rinv_vector_back[i]).determinant()) << endl;
        per_deter = abs(deter);
    }
    
    return i_match;
}




