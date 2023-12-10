#include "TPQ-MPS/TPQ-MPS.h"
#include "itensor/all.h"
#include <vector>
#include <array>
#include <iostream>
#include <ios>

 
int main(){
    std::ios_base::sync_with_stdio(false);
    std::string filename = "Heisenberg_1D";

    int N = 21;
    double J = 1;
    double beta = 0.05;
    int auxiliaries = 5;
    itensor::SpinHalf sites;

    int TimeSteps = 200;
    int Evols = 20;

    auto H_U = TPQ_MPS::Create_Heisenberg_Model_1D(N,J,beta,sites,auxiliaries);
    std::vector<std::array<double,2>> Energies = TPQ_MPS::Calculate_Energies(TimeSteps,Evols,H_U,sites);
    TPQ_MPS::Save_Data(filename,Energies);

    for (auto& e : Energies){
        std::cout << e[0]/20 << " " << e[1]/20 << "\n";
    }

}







