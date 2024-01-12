#include "TPQ-MPS-v2/TPQ-MPS.h"
#include "itensor/all.h"
#include <iostream>
#include <array>



int main(){
    itensor::SpinHalf sites;
    int N = 8;
    int M = 6;
    std::array<double,3> K = {1,1,1};
    std::array<double,3> h = {1,1,1};
    int auxiliaries = 7;
    std::string filename = "E1";

    double beta = 0.05;
    int TimeSteps = 120;
    int Evols = 50;
    
    auto H = TPQ_MPS::Create_Kitaev_Honeycomb_Model_2D(N,M,K,h,sites,auxiliaries);
    auto Energies = TPQ_MPS::Calculate_Energies_TDVP(TimeSteps,Evols,beta,H,sites);
    TPQ_MPS::Save_Data_txt(filename,Energies);

}





