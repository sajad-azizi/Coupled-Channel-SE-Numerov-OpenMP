#pragma once



#include "./Equations.hpp"

class Eigenvalues{
    
    public:
        
        explicit Eigenvalues( Equations *equations, Parameters *parameters);
        void Boundstates_finder();
        void groundstate_finder();
        int matching_point_finder(double Energy);
        
        double gsEnergy;
        int i_match;
        
    private:
    
        Equations *equations;
        Parameters *parameters;
        //energy space
        double Emin;
        double Emax;
        double Energy;
        int N_grid;
        int channels;
    
    
};
