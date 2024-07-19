#include <iostream>
#include <array>
#include <vector>
#include "itensor/all.h"
#include "TDVP-MPS/main.h"
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
    TDVP_MPS::Hamiltonian H_Details;
    H_Details.set("K",1.);
    TPQ_MPS_old::Hamiltonian H_Details_old;
    H_Details_old.set("K",4./3.);
    H_Details.set("hz",0.1);
    //H_Details.set("J",1);

    char env[] = "MKL_NUM_THREADS=1";
    putenv(env);

    int LX = 2;
    int LY = 2;
    int auxiliaries = 1;
    int sec_auxiliaries = 1;

    itensor::CustomSpin sites;
    itensor::SpinHalf sites_old;
    int spin = 1;
    std::string filename = "aaaaaaaaaa";

    std::vector<int> timesteps= {50,100};
    int Evols = 2;
    std::vector<double> intervals = {1,9};
    int init_states = 32;
    int max_states = 128;

    auto t = std::chrono::system_clock::now();
    std::chrono::time_point<std::chrono::system_clock,std::chrono::duration<double>> tmin;
    auto t0 = std::chrono::duration<double>(t-tmin);
    long long t00 = static_cast<long long>(t0.count());
    long long t000 = t00; //- (t00 / 10000) * 10000;
    std::cout << t000 << "\n" << std::flush;
    std::cout << t000/10000000 * 10000000 << "\n" << std::flush;
    


    auto Model = TDVP_MPS::Kitaev_Model(LX,LY,H_Details,spin,auxiliaries,sec_auxiliaries,"Triangular");
    //Model.TPQ_MPS(timesteps,intervals,Evols,256,32,"TwoSite",0.005,"z");
    Model.TPQ_MPS(timesteps,intervals,Evols);
    Model.Save("benis3");
    //auto Model2 = TPQ_MPS_old::Kitaev_Model(LX,LY,H_Details_old,sites_old,auxiliaries,"Honeycomb");
    //itensor::PrintData(Model.H2);
    //itensor::PrintData(Model.H2);
    //itensor::PrintData(Model.H_flux);

    
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

    auto Szampo = AutoMPO(Model.sites);
    Szampo += "Sz",3;
    MPO Sz = toMPO(Szampo);

    
    auto Hexphalf = itensor::toExpH(Model.ampo,0.1);
    auto Hexphalfinv = itensor::toExpH(Model.ampo,-0.1);

    //itensor::println(Hexp);
    //itensor::println(Model.H0);

  
    
    //itensor::println(Hexphalf);
    //itensor::println(Hexphalfinv);
    //auto rmps = itensor::randomMPS(Model.sites,256); 
    //itensor::println(rmps);

    //Model.Tan_Evolution(TimeSteps,intervals);

    

    /*
    auto Flux12 = Model.honeycomb_flux_operator_half(LX,LY,auxiliaries);
    auto Fluxex = Model.honeycomb_flux_operator(LX,LY,auxiliaries);
    auto Hexp = itensor::toExpH(Model.ampo,itensor::Cplx_1);

    for (int i = 0; i != 50; i++){
        auto rmps = itensor::randomMPS(Model.sites,64);

        std::cout << itensor::innerC(rmps,Flux12,Hexp,rmps) << " ";
        std::cout << itensor::innerC(rmps,Fluxex,Hexp,rmps) << "\n";
    }
    */


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

