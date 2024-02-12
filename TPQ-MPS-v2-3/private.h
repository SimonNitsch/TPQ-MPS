#include <iostream>
#include "itensor/all.h"
#include <array>
#include <vector>
#include <cmath>
#include "TDVP/tdvp.h"
#include "TDVP/basisextension.h"
#include <type_traits>
#include <exception>
#pragma once

/*
This is the header file I use for backend functions
For this purpose the namespace priv is defined
Kinda like the private methods of a class, only that up until now I felt like splitting the code up like I do is a bit easier to handle
than packing everything in one big class.
Maybe I'll change that in the future though.
*/

namespace priv{

// The TDVP loop for Imaginary Time Evolution
// Probably responsibly for more than 95% of the computation time
void tdvp_loop(std::vector<double>& E_vec, itensor::MPO& H, itensor::MPS& psi, itensor::Cplx t, int Sweeps, int TimeSteps, int data){
    
    int count = 0;
    for (int j = 0; j != TimeSteps; j++){
        double E = itensor::tdvp(psi,H,t,Sweeps,{"DoNormalize",true,
                                                "Quiet",true,
                                                "ErrGoal",1E-7});
            
        if (j*data >= count*TimeSteps){
            E_vec.push_back(E);
            count++;
        }
    }
}



std::vector<std::array<double,2>> Mean(std::vector<std::vector<double>>& M){
    int M0 = M[0].size();
    int M1 = M.size();
    std::vector<std::array<double,2>> vec;
    vec.reserve(M0);
    
    for (int i = 0; i != M0; i++){
        double v = 0;
        double v2 = 0;
        for (int j = 0; j != M1; j++){
            v += M[j][i];
            v2 += M[j][i] * M[j][i];
        }
        double mean = v/static_cast<double>(M1);
        double std = std::sqrt(v2/M1 - mean*mean);
        std::array<double,2> proto_vec = {mean,std};
        vec.emplace_back(proto_vec);
    }
    
    return vec;

}



std::array<int,3> get_neighbour_data_hex(int N, int M, int pos){
    int X = (pos-1) % N;
    int Y = (pos-1) / N;
    std::array<int,3> neighbours{};

    if (X != (N-1)){
        neighbours[0] = pos+1;
    }
    if (Y != 0){
        neighbours[1] = pos-N+1;
    }
    if (pos != 1 && pos != M*N){
        neighbours[2] = pos-1;
    }
    return neighbours;
    
}



std::array<int,3> get_neighbour_data_tri(int N, int M, int pos){
    int X = (pos-1) / N;
    int Y = (pos-1) % N;
    std::array<int,3> neighbours = {0,0,0};

    if (X != (M-1)){
        neighbours[0] = pos + N;
    }
    if (Y == (N-1)){
        neighbours[1] = pos + 1 - N;
    }
    else{
        neighbours[1] = pos + 1;
    }
    if (X != (M-1)){
        if (Y == 0){
            neighbours[2] = pos + 2*N - 1;
        }
        else{
            neighbours[2] = pos + N - 1;
        }
    }
    return neighbours;

}

std::array<int,3> get_neighbour_data(int N, int M, int pos, std::string form){
    std::array<int,3> neighbour;
    if (form=="Honeycomb"){
        neighbour = get_neighbour_data_hex(N,M,pos);
    }
    else if (form=="Triangular"){
        neighbour = get_neighbour_data_tri(N,M,pos);
    }
    else{
        throw std::invalid_argument("Invalid Lattice Shape");
    }
    return neighbour;
}



template<typename T, typename = std::enable_if<std::is_base_of<itensor::SiteSet, T>::value>>
itensor::AutoMPO kitaev_honeycomb_2D(int N, int M, std::array<double,3>& K, T& sites, int aux){
    sites = T(N*M+2*aux, {"ConserveQNs=",false});
    auto ampo = itensor::AutoMPO(sites);

    std::vector<int> full_points;
    full_points.reserve(((N+1)/2)*M);
    for (int m = 0; m != M; m++){
        for (int n = 1; n <= N; n+=2){
            full_points.push_back(n+N*m);
        }   
    }

    for (int& i : full_points){
        std::array<int,3> neighbours = get_neighbour_data(N,M,i,"Honeycomb");
        if (neighbours[0] != 0){
            ampo += K[0],"Sx",i+aux,"Sx",neighbours[0]+aux;
        }
        if (neighbours[1] != 0){
            ampo += K[1],"Sy",i+aux,"Sy",neighbours[1]+aux;
        }
        if (neighbours[2] != 0){
            ampo += K[2],"Sz",i+aux,"Sz",neighbours[2]+aux;
        }
    }

    return ampo;
}






template<typename T, typename = std::enable_if<std::is_base_of<itensor::SiteSet, T>::value>>
itensor::AutoMPO kitaev_tri_2d(int N, int M, std::array<double,3>& K, T& sites, int aux){
    sites = T(N*M+2*aux,{"ConserveQNs=",false});
    auto ampo = itensor::AutoMPO(sites);

    for (int i = 1; i != (N*M+1); i++){
        std::array<int,3> neighbours = get_neighbour_data(N,M,i,"Triangular");
        if (neighbours[0] != 0){
            ampo += K[0],"Sx",i+aux,"Sx",neighbours[0]+aux;
        }
        if (neighbours[1] != 0){
            ampo += K[1],"Sy",i+aux,"Sy",neighbours[1]+aux;
        }
        if (neighbours[2] != 0){
            ampo += K[2],"Sz",i+aux,"Sz",neighbours[2]+aux;
        }
    }

    return ampo; 
    
}



void add_heisenberg_interaction(itensor::AutoMPO& ampo, double J, int N, int M, int aux, std::string form){
    std::vector<int> full_points;

    if (form=="Honeycomb"){
        full_points.reserve(((N+1)/2)*M);
        for (int m = 0; m != M; m++){
            for (int n = 1; n <= N; n+=2){
                full_points.push_back(n+N*m);
            }   
        }
    }   
    else if (form=="Triangular"){
        full_points.reserve(N*M);
        for (int i = 1; i <= N*M; i++){
            full_points.push_back(i);
        }
    }  

    for (int& i : full_points){
        std::array<int,3> neighbours = get_neighbour_data(N,M,i,form);
        if (neighbours[0] != 0){
            ampo += J,"Sz",i+aux,"Sz",neighbours[0]+aux;
            ampo += J*0.5,"S+",i+aux,"S-",neighbours[0]+aux;
            ampo += J*0.5,"S-",i+aux,"S+",neighbours[0]+aux;
        }
        if (neighbours[1] != 0){
            ampo += J,"Sz",i+aux,"Sz",neighbours[1]+aux;
            ampo += J*0.5,"S+",i+aux,"S-",neighbours[1]+aux;
            ampo += J*0.5,"S-",i+aux,"S+",neighbours[1]+aux;
        }
        if (neighbours[2] != 0){
            ampo += J,"Sz",i+aux,"Sz",neighbours[2]+aux;
            ampo += J*0.5,"S+",i+aux,"S-",neighbours[2]+aux;
            ampo += J*0.5,"S-",i+aux,"S+",neighbours[2]+aux;
        }   
    }

}






void add_magnetic_field(itensor::AutoMPO& ampo, std::array<double,3>& h, int length, int aux){
    for (int i = (aux+1); i != (length-aux+1); i++){
        ampo += -1*h[0],"Sx",i;
        ampo += -1*h[1],"Sy",i;
        ampo += -1*h[2],"Sz",i;
    }
}








}
