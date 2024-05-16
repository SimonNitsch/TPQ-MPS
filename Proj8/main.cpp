#include <iostream>
#include <array>
#include <vector>
#include "itensor/all.h"
#include "TDVP/tdvp.h"
#include "TPQ-MPS-v3-10/main.h"
#include "TPQ-MPS-v3-07/main.h"
#include <ios>
#include <chrono>
#include "Eigen/Dense"
#include "unsupported/Eigen/MatrixFunctions"

/*
int main(){
    Eigen::MatrixX<itensor::Cplx> a = Eigen::MatrixXcd::Zero(3,3);
    a(0,0) = -1. * M_PI * itensor::Cplx_i;
    a(2,2) = M_PI * itensor::Cplx_i;
    a(1,2) = 0.5 * itensor::Cplx_1;

    Eigen::MatrixX<itensor::Cplx> c = a.exp();
    std::cout << c << "\n";

    for (int i = 0; i != 3; i++){
        for (int j = 0; j != 3; j++){
            auto b = c(i,j);
            b += itensor::Cplx_1;
            std::cout << std::real(b) << std::imag(b) << " ";
        }
        std::cout << "\n";
    }
    std::cout << c << "\n";
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
    std::cout << t000/10000000 * 10000000 << "\n" << std::flush;
    


    auto Model = TPQ_MPS::Kitaev_Model(LX,LY,H_Details,spin,auxiliaries,"Honeycomb");
    auto Model2 = TPQ_MPS_old::Kitaev_Model(LX,LY,H_Details_old,sites_old,auxiliaries,"Honeycomb");
    //itensor::PrintData(Model.H2);
    //itensor::PrintData(Model.H2);
    //itensor::PrintData(Model.H_flux);
    
    auto H02 = Model.H0(10);
    //itensor::PrintData(H02);
    auto H02ind = itensor::inds(H02);
    auto it = itensor::ITensor(H02ind);
    
    /*
    for (auto I : itensor::iterInds(H02)){
        std::cout << "a " << std::flush;
        for (auto& i : I){
            std::cout << itensor::val(i) << " ";
        }
        auto e = itensor::eltC(H02,I);
        std::cout << e << "\n";
    }
    auto psi2 = itensor::MPS(Model.sitestan);
    auto psi3 = itensor::MPO(Model.sitestan);

    auto ind1 = rightLinkIndex(Model.H0,1);
    auto ind2 = siteInds(psi3,1)(1);
    auto ind3 = siteInds(psi3,1)(2);
    auto h1 = ITensor(ind1,ind2,ind3);
    psi3.set(1,h1);
    
    //itensor::PrintData(it);
    itensor::PrintData(itensor::siteInds(psi3,2));
    */

    auto Hexp = itensor::toExpH(Model.ampo,-1);
    itensor::println(Hexp);
    itensor::println(Model.H0);

    auto H0tan = Model.mpo_to_tanmpo(Model.H0);
    auto Hmpstan = Model.mpo_to_tanmps(Hexp);

    for (int i = 0; i != 50; i++){
        auto rmps = itensor::randomMPS(Model.sites);
        auto rmpstan = Model.mps_to_tanmps(rmps);

        std::cout << itensor::innerC(rmps,Hexp,Model.H0,rmps) << " " << std::flush;
        std::cout << itensor::innerC(Hmpstan,H0tan,rmpstan) << "\n" << std::flush;
    }






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

