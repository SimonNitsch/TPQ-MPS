#include "TPQ-MPS/TPQ-MPS.h"
#include "itensor/all.h"
#include <vector>
#include <array>
#include <iostream>
#include <ios>
#include "TDVP/tdvp.h"
#include <complex>

 
int main(){
    std::ios_base::sync_with_stdio(false);
    std::string filename = "Heisenberg_2D";

    int N = 8;
    int M = 6;
    double J = 1;
    double beta = 0.05;
    double K = 1;
    int auxiliaries = 5;
    itensor::SpinHalf sites;

    int TimeSteps = 200;
    int Evols = 20;
    int random_states = 1;

    auto H_U = TPQ_MPS::Create_Heisenberg_Model_2D(N,M,J,beta,sites,auxiliaries);
    itensor::PrintData(H_U[0]);
    std::vector<std::array<double,2>> Energies = TPQ_MPS::Calculate_Energies_WI(TimeSteps,Evols,H_U,sites,random_states);
    TPQ_MPS::Save_Data(filename,Energies);

    for (auto& e : Energies){
        std::cout << e[0] << " " << e[1] << "\n";
    }



}







