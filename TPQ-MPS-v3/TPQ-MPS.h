#include <iostream>
#include <array>
#include <vector>
#include <chrono>
#include "itensor/all.h"
#include <ios>
#include <type_traits>
#include <cmath>
#include <cstdio>
#include <string>
#include <algorithm>
#include <cstdlib>
#include "TDVP/tdvp.h"
#include <fstream>
#include "Hamiltonian.h"
#include <exception>
#pragma once


namespace TPQ_MPS{



class Kitaev_Model{

    private:
    itensor::AutoMPO ampo;
    int Lattice_Type; // 1 for Honeycomb, 2 for Triangular
    int Calc_Type; // 1 for DMRG, 2 for TDVP

    itensor::MPO H0;
    Hamiltonian H_Details;
    itensor::SiteSet sites;

    double GSE;
    std::vector<std::array<double,2>> E;
    std::vector<std::array<double,2>> Cv;
    std::vector<std::array<double,2>> S;


    std::array<int,3> get_neighbour_data_hex(int LX, int LY, int pos);
    std::array<int,3> get_neighbour_data_tri(int LX, int LY, int pos);

    std::array<int,3> get_neighbour_data(int LX, int LY, int pos){
        std::array<int,3> n;
        if (Lattice_Type == 1){
            n = get_neighbour_data_hex(LX,LY,pos);
        } else if (Lattice_Type == 2){
            n = get_neighbour_data_tri(LX,LY,pos);
        }
        return n;
    }


    void add_kitaev_interaction(int LX, int LY, std::vector<int>& p_vec, int aux);
    void add_magnetic_interaction(int LX, int LY, int aux);
    void add_heisenberg_interaction(int LX, int LY, std::vector<int>& p_vec, int aux);
    void add_gamma_interaction(int LX, int LY, std::vector<int>& p_vec, int aux);
    void add_gammaq_interaction(int LX, int LY, std::vector<int>& p_vec, int aux);

    std::vector<std::array<double,2>> Mean(std::vector<std::vector<double>>& M);
    void tdvp_loop(std::vector<double>& E_vec, itensor::MPS& psi, itensor::Cplx t, int Sweeps, int TimeSteps);


    template<std::size_t n, typename T>
    void save_data(std::string filename, std::vector<std::array<T,n>>& vec);
    
    template<typename T>
    void save_data(std::string filename, T v);



    public:
    template<typename T, typename = std::enable_if<std::is_base_of<itensor::SiteSet, T>::value>>
    Kitaev_Model(int LX, int LY, Hamiltonian H_Details, T& sites, int aux, std::string shape){
        sites = T((LX*LY+2*aux),{"ConserveQNs=",false});
        auto ampo = itensor::AutoMPO(sites); 
        
        this -> ampo = ampo;
        this -> H_Details = H_Details;
        this -> sites = sites;

        std::vector<int> full_points;
        if (shape == "Honeycomb"){
            Lattice_Type = 1;
            full_points.reserve(((LY+1)/2)*LX);
            for (int m = 0; m != LX; m++){
                for (int n = 1; n <= LY; n+=2){
                    full_points.push_back(n+LX*m);
                }   
            }
        } else if (shape == "Triangular"){
            Lattice_Type = 2;
            full_points.reserve(LX*LY);
            for (int i = 1; i != LX*LY+1; i++){
                full_points.push_back(i);
            }

        } else {
            std::invalid_argument("Invalid Lattice Shape");
        }

        add_kitaev_interaction(LX,LY,full_points,aux);
        add_magnetic_interaction(LX,LY,aux);
        add_heisenberg_interaction(LX,LY,full_points,aux);
        add_gamma_interaction(LX,LY,full_points,aux);
        add_gammaq_interaction(LX,LY,full_points,aux);
        
        H0 = itensor::toMPO(ampo);

    }

/*

    template<typename T, typename = std::enable_if<std::is_base_of<itensor::SiteSet, T>::value>>
    Kitaev_Model* Honeycomb(int LX, int LY, Hamiltonian H_Details, T& sites, int aux){
        sites = T((LX*LY+2*aux),{"ConserveQNs=",false});
        auto ampo = itensor::AutoMPO(sites);
        

        std::vector<int> full_points;
        full_points.reserve(((LY+1)/2)*LX);
        for (int m = 0; m != LX; m++){
            for (int n = 1; n <= LY; n+=2){
                full_points.push_back(n+LX*m);
            }   
        }

        add_kitaev_interaction(LX,LY,full_points,aux);
        add_magnetic_interaction(LX,LY,aux);
        add_heisenberg_interaction(LX,LY,full_points,aux);
        add_gamma_interaction(LX,LY,full_points,aux);
        add_gammaq_interaction(LX,LY,full_points,aux);
        
        auto H0 = itensor::toMPO(ampo);
        H_list[0] = &H0;
        return this;
    }

    template<typename T, typename = std::enable_if<std::is_base_of<itensor::SiteSet, T>::value>>
    Kitaev_Model* Triangular(int LX, int LY, Hamiltonian H_Details, T& sites, int aux){
        sites = T((LX*LY+2*aux),{"ConserveQNs=",false});
        auto ampo = itensor::AutoMPO(sites);
        Kitaev_Model(ampo,H_Details,sites,2);

        std::vector<int> full_points;
        full_points.reserve(LX*LY);
        for (int i = 1; i != LX*LY+1; i++){
            full_points.push_back(i);
        }

        add_kitaev_interaction(LX,LY,full_points,aux);
        add_magnetic_interaction(LX,LY,aux);
        add_heisenberg_interaction(LX,LY,full_points,aux);
        add_gamma_interaction(LX,LY,full_points,aux);
        add_gammaq_interaction(LX,LY,full_points,aux);

        auto H0 = itensor::toMPO(ampo);
        H_list[0] = &H0;        
        return this;
    }
*/


    void Time_Evolution(int TimeSteps, std::vector<double> intervals, int Evols, bool Heat_Capacity=true, int init_rand_sites=32, int Sweeps=5);
    std::array<std::vector<std::array<double,2>>,2> Calculate_Heat_Capacity(int TimeSteps, std::vector<double>& intervals, std::vector<std::vector<double>>& Energies);

    Hamiltonian Get_Constants(){
        return H_Details;
    }

    void Print_Interactions(){
        itensor::PrintData(ampo);
    }

    void Save(std::string x){
        if (Calc_Type == 1){
            save_data(x,GSE);
        }
        else if (Calc_Type == 2){
            save_data(x,E);
        }
    }

    
    void DMRG(int Sweeps=10){
        Calc_Type = 1;
        auto psi0 = itensor::randomMPS(sites,16);
        auto [energy, psi] = itensor::dmrg(H0,psi0,Sweeps,{"Quiet",true});
        GSE = energy;
    }





};






}










