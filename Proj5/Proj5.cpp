#include "itensor/all.h"
#include "TPQ-MPS/TPQ-MPS.h"



int main(){
    itensor::SpinHalf sites;    
    int N = 8;
    int M = 6;
    int aux = 7;
    std::array<double,3> K = {1,1,1};
    std::array<double,3> h = {1,1,1};

    int tsteps = 60;
    double beta = 0.05;
    int evols = 40;

    //std::string filename = "K_Honey_E";

    auto H = TPQ_MPS::Create_Kitaev_Honeycomb_Model_2D(N,M,K,h,sites,aux);
    /*
    double E = TPQ_MPS::DMRG(H,sites);

    std::string filename = "GSE";
    TPQ_MPS::Save_Data_txt(filename,E);
    */

   //std::vector<std::array<double,2>> Energies = TPQ_MPS::Calculate_Energies_TDVP(tsteps,evols,beta,H,sites);
   //TPQ_MPS::Save_Data_txt(filename,Energies);

}



