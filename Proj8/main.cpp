#include <iostream>
#include <array>
#include <vector>
#include "itensor/all.h"
#include "TPQ-MPS-v3-14/main.h"
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
    H_Details.set("K",-1.);
    TPQ_MPS_old::Hamiltonian H_Details_old;
    H_Details_old.set("K",4./3.);
    //H_Details.set("hz",0.15);
    //H_Details.set("J",1);

    int LX = 3;
    int LY = 2;
    int auxiliaries = 0;

    itensor::CustomSpin sites;
    itensor::SpinHalf sites_old;
    int spin = 1;
    std::string filename = "aaaaaaaaaa";

    int TimeSteps = 50;
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
    


    auto Model = TPQ_MPS::Kitaev_Model(LX,LY,H_Details,spin,auxiliaries,"Triangular");
    //Model.Tan_Evolution(TimeSteps,intervals,256,1,5);
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

    auto Hexphalf_tan = Model.mpo_to_tanmpo(Hexphalf);
    auto Hexphalfinv_tan = Model.mpo_to_tanmpo(Hexphalfinv);
    //itensor::println(Hexp);
    //itensor::println(Model.H0);

    auto H0tan = Model.mpo_to_tanmpo(Model.H0);
    auto Hmpstan = Model.mpo_to_tanmps(Hexphalf);
    auto Hmpstaninv = Model.mpo_to_tanmps(Hexphalfinv);

    MPO Hcalc3 = nmultMPO(Hexphalfinv,prime(Sz));
    Hcalc3 = nmultMPO(Hcalc3,prime(Hexphalf));
    Hcalc3 = nmultMPO(Hcalc3,prime(Hexphalf));
    Hcalc3 = nmultMPO(Hcalc3,prime(Sz));
    Hcalc3 = nmultMPO(Hcalc3,prime(Hexphalfinv));

    MPO Hcalc3inv = nmultMPO(Hexphalf,prime(Sz));
    Hcalc3inv = nmultMPO(Hcalc3inv,prime(Hexphalfinv));
    Hcalc3inv = nmultMPO(Hcalc3inv,prime(Hexphalfinv));
    Hcalc3inv = nmultMPO(Hcalc3inv,prime(Sz));
    Hcalc3inv = nmultMPO(Hcalc3inv,prime(Hexphalf));

    MPS Hcalc3mps = applyMPO(H0tan,Hmpstaninv);
    Hcalc3mps = applyMPO(Hexphalf_tan,Hcalc3mps);
    Hcalc3mps = applyMPO(Hexphalf_tan,Hcalc3mps);
    Hcalc3mps = applyMPO(H0tan,Hcalc3mps);
    Hcalc3mps = applyMPO(Hexphalfinv_tan,Hcalc3mps);
    MPO Hcalc32 = Model.tanmps_to_mpo(Hcalc3mps);

    MPO Hcalc4 = nmultMPO(Model.H0,prime(Hexphalf));
    Hcalc4 = nmultMPO(Hcalc4,prime(Hexphalf));
    Hcalc4 = nmultMPO(Hcalc4,prime(Model.H0));

    MPO Hcalc5 = itensor::nmultMPO(Hexphalf_tan,prime(H0tan));
    Hcalc5 = nmultMPO(Hcalc5,prime(Hexphalf_tan));

    MPO Hcalc6 = nmultMPO(H0tan,prime(Hexphalf_tan));
    Hcalc6 = nmultMPO(Hcalc6,prime(Hexphalf_tan));
    Hcalc6 = nmultMPO(Hcalc6,prime(H0tan));

    Hmpstan.orthogonalize();
    //Hmpstan.normalize();
    Hmpstaninv.orthogonalize();
    Hmpstaninv.normalize();
    //itensor::println(Hcalc);


    MPO Hcalc4t = nmultMPO(Model.H0,prime(Hexphalf));
    //Hcalc4t.noPrime();
    MPO Hcalc4td = Hcalc4t;


    MPO Hcalc4tan = Model.mpo_to_tanmpo(Hcalc4t);




    double E1 = 0.;
    double E2 = 0.;
    double Ec1 = 0.;
    double Ec2 = 0.;
    double Enorm = 0.;

    for (int i = 0; i != 1000; i++){
        auto rmps = itensor::randomMPS(Model.sites,512);  
        //auto rmps2 = itensor::randomMPS(Model.sitestan,128);
        std::complex<double> normsq = innerC(rmps,Hexphalfinv,Hexphalfinv,rmps) * innerC(rmps,Hexphalf,Hexphalf,rmps);      

        auto delE = itensor::innerC(rmps,Hcalc3,rmps) / std::real(normsq);
        auto delE2 = itensor::innerC(rmps,Hcalc3inv,rmps) / std::real(normsq);
        E1 += std::real(delE);
        E2 += std::real(delE) * std::real(delE);
        Ec1 += std::real(delE2);
        Ec2 += std::real(delE2) * std::real(delE2);
        Enorm += std::real(normsq);

        std::cout << i << " - " << delE << " - " << delE2 << ", " << normsq << "\n" << std::flush;
    }
    E1 /= 1000.;
    E2 /= 1000.;
    Ec1 /= 1000.;
    Ec2 /= 1000.;
    Enorm /= 1000.;

    
    std::cout << E1 << "--" << std::sqrt(E2-E1*E1) / std::sqrt(1000.) << " -- " << Enorm << "\n";
    std::cout << Ec1 << "--" << std::sqrt(Ec2-Ec1*Ec1) / std::sqrt(1000.) << " -- " << Enorm << "\n";
    std::complex<double> normsqtan = itensor::innerC(Hmpstan,Hmpstan);
    
    /*
    double rnorm = 0;
    for (int i = 0; i != 50; i++){
        auto rmps = itensor::randomMPS(Model.sites,128);
        std::complex<double> normrand = itensor::innerC(rmps,Hexphalf,Hexphalf,rmps);
        std::cout << "RandNorm, " << i << ": " << normrand << "\n" << std::flush;
        rnorm += std::real(normrand);
    }
    rnorm /= 50.;
    std::cout << "Norm: " << rnorm << "\n" << std::flush;
    */
    

    auto tanmps = itensor::applyMPO(H0tan,Hmpstaninv);
    //tanmps = itensor::applyMPO(Hmpotan,tanmps);
    

    std::cout << itensor::innerC(Hmpstaninv,Hcalc4tan,Hmpstaninv) << " " << std::real(normsqtan) << "\n" << std::flush;
    
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

