#include "TPQ-MPS-v2/TPQ-MPS.h"
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
    //std::array<double,3> h = {1,1,1};
    int auxiliaries = 7;
    std::string filename = "GSE";

    auto H = TPQ_MPS::Create_Kitaev_Triangular_Model_2D(N,M,K,sites,auxiliaries);
    double E = TPQ_MPS::DMRG(H,sites);
    TPQ_MPS::Save_Data_txt(filename,E);


}




