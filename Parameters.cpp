#include "./Parameters.hpp"

Parameters::Parameters(int N_grid, int channels, double dr, int Nroots, int NTHREADS, int N_phi, double dphi, double Emin, double Emax, int divide, int p, int external_parameter):
N_grid(N_grid), channels(channels), dr(dr), Nroots(Nroots), NTHREADS(NTHREADS), N_phi(N_phi), dphi(dphi) , Emin(Emin), Emax(Emax) , divide(divide), p(p), external_parameter(external_parameter){
    this->N_grid = N_grid;
    this->channels = channels;
    this->dr = dr;
    this->Nroots = Nroots;
    this->NTHREADS = NTHREADS;
    this->external_parameter = external_parameter;
    this->Emin = Emin;
    this->Emax = Emax;
    this->dphi = dphi;
    this->N_phi = N_phi;
    this->divide = divide;
    this->p = p;
}
