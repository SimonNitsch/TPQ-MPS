#include <iostream>
#include "itensor/all.h"
#include <array>
#include <vector>
#include <chrono>
#include <complex>
#include <ios>
#include <type_traits>
#include <cmath>
#include <cstdio>
#include "private.h"
#include <string>
#include <algorithm>
#include <cstdlib>
#include "TDVP/tdvp.h"
#include <fstream>
#pragma once

namespace TPQ_MPS{

/*
std::array<itensor::MPO*,2> Create_Heisenberg_Model_1D(int N, double J, double beta, itensor::SpinHalf& sites){
    sites = itensor::SpinHalf(N,{"ConserveQNs=",false});
    auto ampo = itensor::AutoMPO(sites);

    for (int i = 1; i != N; i++){
        ampo += J,"Sz",i,"Sz",i+1;
        ampo += J*0.5,"S+",i,"S-",i+1;
        ampo += J*0.5,"S-",i,"S+",i+1;
    }
    auto H = itensor::toMPO(ampo);
    auto U = itensor::toExpH(ampo,beta);

    std::array<itensor::MPO,2> output;
    output[0] = H;
    output[1] = U;
    return output;
}*/


template<typename T, typename = std::enable_if<std::is_base_of<itensor::SiteSet, T>::value>>

std::array<itensor::MPO,3> Create_Heisenberg_Model_1D(int N, double J, double beta, T& sites, int auxiliaries){
    sites = T(N+2*auxiliaries,{"ConserveQNs=",false});
    auto ampo = itensor::AutoMPO(sites);

    for (int i = 1+auxiliaries; i != N+auxiliaries; i++){
        ampo += J,"Sz",i,"Sz",i+1;
        ampo += J*0.5,"S+",i,"S-",i+1;
        ampo += J*0.5,"S-",i,"S+",i+1;
    }
    auto H = itensor::toMPO(ampo);
    auto U1 = itensor::toExpH(ampo,beta*0.5*(itensor::Cplx_1+itensor::Cplx_i)*0.5);
    auto U2 = itensor::toExpH(ampo,beta*0.5*(itensor::Cplx_1-itensor::Cplx_i)*0.5);
    std::array<itensor::MPO,3> output = {H,U1,U2};
    return output;
}





template<typename T, typename = std::enable_if<std::is_base_of<itensor::SiteSet, T>::value>>
std::array<itensor::MPO,3> Create_Kitaev_Honeycomb_Model_2D(int N, int M, double Kx, double Ky, double Kz, double beta, T& sites, int auxiliaries){
    auto H = Create_Kitaev_Honeycomb_Model_2D(N,M,Kx,Ky,Kz,sites,auxiliaries);
    auto U1 = itensor::expH(H,beta*0.5*(itensor::Cplx_1+itensor::Cplx_i)*0.5);
    auto U2 = itensor::expH(H,beta*0.5*(itensor::Cplx_1-itensor::Cplx_i)*0.5);
    std::array<itensor::MPO,3> output = {H,U1,U2};
    return output;

}

template<typename T, typename = std::enable_if<std::is_base_of<itensor::SiteSet, T>::value>>
itensor::MPO Create_Kitaev_Honeycomb_Model_2D(int N, int M, std::array<double,3>& K, T& sites, int auxiliaries){
    auto ampo = priv::kitaev_honeycomb_2D(N,M,K,sites,auxiliaries);
    auto H = itensor::toMPO(ampo);
    return H;
}


template<typename T, typename = std::enable_if<std::is_base_of<itensor::SiteSet, T>::value>>
itensor::MPO Create_Kitaev_Honeycomb_Model_2D(int N, int M, std::array<double,3>& K, std::array<double,3>& h, T& sites, int auxiliaries){
    auto ampo = priv::kitaev_honeycomb_2D(N,M,K,sites,auxiliaries);
    int length = itensor::length(sites);
    priv::add_magnetic_field(ampo,h,length,auxiliaries);
    auto H = itensor::toMPO(ampo);
    return H;
}


template<typename T, typename = std::enable_if<std::is_base_of<itensor::SiteSet,T>::value>>
itensor::MPO Create_Kitaev_Triangular_Model_2D(int N, int M, std::array<double,3>& K, T& sites, int auxiliaries){
    auto ampo = priv::kitaev_tri_2d(N,M,K,sites,auxiliaries);
    auto H = itensor::toMPO(ampo);
    return H;
}


template<typename T, typename = std::enable_if<std::is_base_of<itensor::SiteSet, T>::value>>
itensor::MPO Create_Kitaev_Triangular_Model_2D(int N, int M, std::array<double,3>& K, std::array<double,3>& h, T& sites, int auxiliaries){
    auto ampo = priv::kitaev_tri_2d(N,M,K,sites,auxiliaries);
    int length = itensor::length(sites);
    priv::add_magnetic_field(ampo,h,length,auxiliaries);
    auto H = itensor::toMPO(ampo);
    return H;
}







template<typename T, typename = std::enable_if<std::is_base_of<itensor::SiteSet, T>::value>>

std::array<itensor::MPO,3> Create_Heisenberg_Model_2D(int N, int M, double J, double beta, T& sites, int auxiliaries){
    sites = T(N*M+2*auxiliaries,{"ConserveQNs=",false});
    auto ampo = itensor::AutoMPO(sites);

    for (int i = 1; i != N; i++){
        for (int j = 0; j != M; j++){
            int n = auxiliaries+j*N+i;
            ampo += J,"Sz",n,"Sz",n+1;
            ampo += J*0.5,"S+",n,"S-",n+1;
            ampo += J*0.5,"S-",n,"S+",n+1;
        }
    }

    for (int i = 1; i <= N; i++){
        for (int j = 0; j != (M-1); j++){
            int n = auxiliaries+j*N+i;
            ampo += J,"Sz",n,"Sz",n+N;
            ampo += J*0.5,"S+",n,"S-",n+N;
            ampo += J*0.5,"S-",n,"S+",n+N;
        }
    }
    
    auto H = itensor::toMPO(ampo);
    auto U1 = itensor::toExpH(ampo,beta*0.5*(itensor::Cplx_1+itensor::Cplx_i)*0.5);
    auto U2 = itensor::toExpH(ampo,beta*0.5*(itensor::Cplx_1-itensor::Cplx_i)*0.5);
    std::array<itensor::MPO,3> output = {H,U1,U2};
    return output;


}




std::vector<std::array<double,2>> Calculate_Energies_WI(int TimeSteps, int Evols, std::array<itensor::MPO,3>& H_U, itensor::SiteSet& sites, int init_rand_sites=32, int data_points=100){
    auto H = H_U[0];
    auto U1 = H_U[1];
    auto U2 = H_U[2];

    std::vector<std::vector<double>> Energies;
    Energies.reserve(Evols);

    int data = std::min(data_points,TimeSteps);
    std::vector<double> E_vec;
    E_vec.reserve(data);
    //auto IState = itensor::InitState(sites);


    for (int i = 0; i != Evols; i++){
        auto t1 = std::chrono::system_clock::now();
        auto psi = itensor::randomMPS(sites,init_rand_sites);

        std::complex<double> E = itensor::innerC(psi,H,psi) / itensor::inner(psi,psi);
        E_vec.push_back(std::real(E));

        int count = 1;

        for (int j = 0; j != TimeSteps; j++){
            //std::cout << "a";
            psi = itensor::applyMPO(U1,psi,{"Method=","DensityMatrix","MaxDim=",256,"Cutoff=",1e-6});
            psi = itensor::applyMPO(U2,psi,{"Method=","DensityMatrix","MaxDim=",256,"Cutoff=",1e-6,"Normalize=",true});

            if (j*data >= count*TimeSteps){
                std::complex<double> E = itensor::innerC(psi,H,psi);
                E_vec.push_back(std::real(E));
                count++;
            }
            
        }
        Energies.push_back(E_vec);
        E_vec.clear();
        auto t2 = std::chrono::system_clock::now();
        auto time = std::chrono::duration<double>(t2-t1);
        std::cout << "Finished Evolution Number " << (i+1) << "/" << Evols << ", Time Needed: " << time.count() << " seconds\n" << std::flush;

    }

    std::vector<std::array<double,2>> Mean_Energies = priv::Mean(Energies);
    return Mean_Energies;

}




std::vector<std::array<double,2>> Calculate_Energies_TDVP(int TimeSteps, std::vector<double> intervals, int Evols, itensor::MPO& H, itensor::SiteSet& sites, int init_rand_sites=32, int Sweeps=5, int data_points=100){
    std::vector<std::vector<double>> Energies;
    Energies.reserve(Evols);

    int data = std::min(data_points,TimeSteps);
    std::vector<double> E_vec;
    E_vec.reserve((data + 1) * intervals.size());

    std::vector<itensor::Cplx> T;
    for (auto & it : intervals){
        itensor::Cplx t = -0.5 * it / static_cast<double>(TimeSteps) * itensor::Cplx_1;
        T.emplace_back(t);
    }
    

    for (int i = 0; i != Evols; i++){
        auto t1 = std::chrono::system_clock::now();
        auto psi = itensor::randomMPS(sites,init_rand_sites);
        auto t = T.begin();

        std::complex<double> E = itensor::innerC(psi,H,psi) / itensor::inner(psi,psi);
        E_vec.push_back(std::real(E));
        priv::tdvp_loop(E_vec,H,psi,*t,Sweeps,TimeSteps,data);
        t++;
        
        for (; t != T.end(); t++){
            E_vec.push_back(0);
            priv::tdvp_loop(E_vec,H,psi,*t,Sweeps,TimeSteps,data);
        }

        
        Energies.push_back(E_vec);
        E_vec.clear();
        auto t2 = std::chrono::system_clock::now();
        auto time = std::chrono::duration<double>(t2-t1);
        std::cout << "Finished Evolution Number " << (i+1) << "/" << Evols << ", Time Needed: " << time.count() << " seconds\n" << std::flush;


    }
    std::vector<std::array<double,2>> Mean_Energies = priv::Mean(Energies);
    return Mean_Energies;
    


}



double DMRG(itensor::MPO& H, itensor::SiteSet& sites, int Sweeps=10, bool tdvp = false){
    auto psi0 = itensor::randomMPS(sites,1);

    if (tdvp){
        itensor::Cplx t = -0.1 * itensor::Cplx_1;
        for (int i = 0; i != 5; i++){
            itensor::tdvp(psi0,H,t,6);
        }
    }

    auto [energy, psi] = itensor::dmrg(H,psi0,Sweeps,{"Quiet",true});
    double e = energy;
    return e;

}








template<std::size_t n, typename T>
void Save_Data_bin(std::string& filename, std::vector<std::array<T,n>>& vec, int data_points=0){
    std::ofstream file(filename,std::ios::binary);
    int length = vec.size();

    if (data_points!=0){
        length = std::min(length,data_points);        
    }

    for (int i = 0; i != length; i++){
        int next_index = (i*vec.size())/length;
        std::array<T,n> point = vec[next_index];
        file.write(reinterpret_cast<const char*>(point.data()),point.size()*sizeof(T));
    }
    file.close();

    std::cout << "Data saved as: " << filename << "\n" << std::flush;

}


template<std::size_t n, typename T>
void Save_Data_txt(std::string& filename, std::vector<std::array<T,n>>& vec){
    std::string ending = ".txt";
    std::string full_file = filename + ending;
    std::ofstream file(full_file);

    for (auto& varr : vec){
        auto v = varr.begin();
        file << *v;
        v++;

        for (; v != varr.end(); v++){
            file << ", " << *v;
        }
        file << "\n";

    }
    file.close();
    std::cout << "Data saved as: " << full_file << "\n" << std::flush;
}


template<typename T>
void Save_Data_txt(std::string& filename, T v){
    std::vector<std::array<T,1>> vec;
    std::array<T,1> v_array = {v};
    vec.push_back(v_array);

    Save_Data_txt(filename, vec);
}



}


