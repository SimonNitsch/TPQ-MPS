#include "TPQ-MPS-v2-1/TPQ-MPS.h"
#include <iostream>
#include <array>
#include <vector>
#include "itensor/all.h"
#include "TDVP/tdvp.h"



int main(){
    itensor::SpinHalf sites;
    int N = 8;
    int M = 6;
    std::array<double,3> K = {1,1,1};
    std::vector<double> intervals = {5,10,20};

    int TimeSteps = 100;
    int Evols = 40;
    int auxiliaries = 7;
    std::string filename = "E1";

    auto H = TPQ_MPS::Create_Kitaev_Triangular_Model_2D(N,M,K,sites,auxiliaries);
    std::vector<std::array<double,2>> Energies = TPQ_MPS::Calculate_Energies_TDVP(TimeSteps,intervals,Evols,H,sites);
    TPQ_MPS::Save_Data_txt(filename,Energies);


}




