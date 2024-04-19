#include <iostream>
#include <array>
#include <vector>
#include "itensor/all.h"
#include "TDVP/tdvp.h"
#include "TPQ-MPS-v3-09/main.h"
#include "TPQ-MPS-v3-07/main.h"
#include <ios>
#include <chrono>

/*
int main(){
    auto sites = itensor::SpinHalf(10,{"ConserveQNs=",false});
    auto a = itensor::randomMPS(sites,10);
    itensor::PrintData(a);
}
*/




int main(){
    std::ios_base::sync_with_stdio(false);
    TPQ_MPS::Hamiltonian H_Details;
    H_Details.set("K",4./3.);
    TPQ_MPS_old::Hamiltonian H_Details_old;
    H_Details_old.set("K",4./3.);
    //H_Details.set("hz",0.15);
    //H_Details.set("J",1);

    int LX = 4;
    int LY = 6;
    int auxiliaries = 2;

    itensor::CustomSpin sites;
    itensor::SpinHalf sites_old;
    int spin = 1;
    std::string filename = "aaaaaaaaaa";

    int TimeSteps = 100;
    int Evols = 30;
    std::vector<double> intervals = {5,20,35,40};
    int init_states = 32;
    int max_states = 128;

    auto t = std::chrono::system_clock::now();
    std::chrono::time_point<std::chrono::system_clock,std::chrono::duration<double>> tmin;
    auto t0 = std::chrono::duration<double>(t-tmin);
    long long t00 = static_cast<long long>(t0.count());
    long long t000 = t00; //- (t00 / 10000) * 10000;
    std::cout << t000 << "\n" << std::flush;
    std::cout << t000/10000000 * 10000000 << std::flush;
    


    auto Model = TPQ_MPS::Kitaev_Model(LX,LY,H_Details,spin,auxiliaries,"Honeycomb");
    auto Model2 = TPQ_MPS_old::Kitaev_Model(LX,LY,H_Details_old,sites_old,auxiliaries,"Honeycomb");
    //itensor::PrintData(Model.H2);
    //itensor::PrintData(Model.H2);
    itensor::PrintData(Model.H_flux);
    itensor::PrintData(Model2.H_flux);




/*
    try{
        Model.Time_Evolution(TimeSteps,intervals,Evols,false,max_states,init_states,"OneSite");
    }
    catch(std::exception& e){
        std::cout << e.what();
    }
  */  
    
    //Model.Save(filename);



}

