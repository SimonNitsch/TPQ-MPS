#include "TPQ-MPS.h"
#pragma once


namespace TPQ_MPS{


std::array<int,3> Kitaev_Model::get_neighbour_data_hex(int LX, int LY, int pos){
    int Y = (pos-1) % LY;
    int X = (pos-1) / LY;
    std::array<int,3> neighbours{};

    if (Y != (LY-1)){
        neighbours[0] = pos+1;
    }
    if (X != 0){
        neighbours[1] = pos-LY+1;
    }
    if (pos != 1 && pos != LX*LY){
        neighbours[2] = pos-1;
    }
    return neighbours;

}

/*
std::array<int,3> Kitaev_Model::get_neighbour_data_hex(int LX, int LY, int pos){
    int Y = (pos-1) % LY;
    int X = (pos-1) / LY;
    std::array<int,3> neighbours{};

    if (Y != (LY-1)){
        neighbours[0] = pos+1;
    }
    if (X != 0){
        neighbours[1] = pos-LY+1;
    }
    else{
        neighbours[1] = pos+(LX-1)*LY+1;
    }
    if (Y != 0){
        neighbours[2] = pos-1;
    } 
    else{
        neighbours[2] = pos+LY-1;
    }
    return neighbours;

} */

std::array<int,3> Kitaev_Model::get_neighbour_data_tri(int LX, int LY, int pos){
    int X = (pos-1) / LY;
    int Y = (pos-1) % LY;
    std::array<int,3> neighbours;

    if (X != (LX-1)){
        neighbours[0] = pos + LY;
        if (Y == 0){
            neighbours[2] = pos + 2*LY - 1;
        }
        else{
            neighbours[2] = pos + LY - 1;
        }
    }
    if (Y == (LY-1)){
        neighbours[1] = pos + 1 - LY;
    }
    else{
        neighbours[1] = pos + 1;
    }

    return neighbours;
}





void Kitaev_Model::add_kitaev_interaction(int LX, int LY, std::vector<int>& p_vec, int aux){
    double Kx = H_Details.get("Kx");
    double Ky = H_Details.get("Ky");
    double Kz = H_Details.get("Kz");
    
    for (int& i : p_vec){
        std::array<int,3> n = get_neighbour_data(LX,LY,i);
        if (n[0] != 0){
            ampo += Kx,"Sx",i+aux,"Sx",n[0]+aux;
        }
        if (n[1] != 0){
            ampo += Ky,"Sy",i+aux,"Sy",n[1]+aux;
        }
        if (n[2] != 0){
            ampo += Kz,"Sz",i+aux,"Sz",n[2]+aux;
        }
    }

}


void Kitaev_Model::add_magnetic_interaction(int LX, int LY, int aux){
    double hx = H_Details.get("hx");
    double hy = H_Details.get("hy");
    double hz = H_Details.get("hz");

    for (int i = aux+1; i != LX*LY+aux+1; i++){
        ampo += -1*hx,"Sx",i;
        ampo += -1*hy,"Sy",i;
        ampo += -1*hz,"Sz",i;
    }
}



void Kitaev_Model::add_heisenberg_interaction(int LX, int LY, std::vector<int>& p_vec, int aux){
    double J = H_Details.get("J");

    for (int& i : p_vec){
        std::array<int,3> n = get_neighbour_data(LX,LY,i);
        if (n[0] != 0){
            ampo += J,"Sz",i+aux,"Sz",n[0]+aux;
            ampo += J*0.5,"S+",i+aux,"S-",n[0]+aux;
            ampo += J*0.5,"S-",i+aux,"S+",n[0]+aux;
        }
        if (n[1] != 0){
            ampo += J,"Sz",i+aux,"Sz",n[1]+aux;
            ampo += J*0.5,"S+",i+aux,"S-",n[1]+aux;
            ampo += J*0.5,"S-",i+aux,"S+",n[1]+aux;
        }
        if (n[2] != 0){
            ampo += J,"Sz",i+aux,"Sz",n[2]+aux;
            ampo += J*0.5,"S+",i+aux,"S-",n[2]+aux;
            ampo += J*0.5,"S-",i+aux,"S+",n[2]+aux;
        }
    }

}


void Kitaev_Model::add_gamma_interaction(int LX, int LY, std::vector<int>& p_vec, int aux){
    double G = H_Details.get("Gamma");

    for (int& i : p_vec){
        std::array<int,3> n = get_neighbour_data(LX,LY,i);
        if (n[0] != 0){
            ampo += G,"Sy",i+aux,"Sz",n[0]+aux;
            ampo += G,"Sz",i+aux,"Sy",n[0]+aux;
        }
        if (n[1] != 0){
            ampo += G,"Sx",i+aux,"Sz",n[1]+aux;
            ampo += G,"Sz",i+aux,"Sx",n[1]+aux;
        }
        if (n[2] != 0){
            ampo += G,"Sx",i+aux,"Sy",n[2]+aux;
            ampo += G,"Sy",i+aux,"Sx",n[2]+aux;
        }
    }

}



void Kitaev_Model::add_gammaq_interaction(int LX, int LY, std::vector<int>& p_vec, int aux){
    double GQ = H_Details.get("GammaQ");

    for (int& i : p_vec){
        std::array<int,3> n = get_neighbour_data(LX,LY,i);
        if (n[0] != 0){
            ampo += GQ,"Sy",i+aux,"Sx",n[0]+aux;
            ampo += GQ,"Sx",i+aux,"Sy",n[0]+aux;
            ampo += GQ,"Sz",i+aux,"Sx",n[0]+aux;
            ampo += GQ,"Sx",i+aux,"Sz",n[0]+aux;
        }
        if (n[1] != 0){
            ampo += GQ,"Sx",i+aux,"Sy",n[1]+aux;
            ampo += GQ,"Sy",i+aux,"Sx",n[1]+aux;
            ampo += GQ,"Sz",i+aux,"Sy",n[1]+aux;
            ampo += GQ,"Sy",i+aux,"Sz",n[1]+aux;
        }
        if (n[2] != 0){
            ampo += GQ,"Sx",i+aux,"Sz",n[2]+aux;
            ampo += GQ,"Sz",i+aux,"Sx",n[2]+aux;
            ampo += GQ,"Sz",i+aux,"Sy",n[2]+aux;
            ampo += GQ,"Sy",i+aux,"Sz",n[2]+aux;
        }
    }

}



itensor::MPO Kitaev_Model::honeycomb_flux_operator(int LX, int LY, int aux){
    int Py = std::max((LY/2)-1,0);
    int Px = LX-1;
    auto fluxop = itensor::AutoMPO(sites);

    for (int i = 0; i != Px; i++){
        for (int j = 0; j != Py; j++){
            int f = LY*i + 2*j + 2 + aux;
            fluxop += 64/Px/Py,"Sx",f,"Sy",f+1,"Sz",f+2,"Sx",f+LY+1,"Sy",f+LY,"Sz",f+LY-1;
        }
    }
    auto fluxH = itensor::toMPO(fluxop);
    return fluxH;
}


std::array<itensor::MPO,3> Kitaev_Model::magnetization_operators(int LX, int LY, int aux){
    auto mx = itensor::AutoMPO(sites);
    auto my = itensor::AutoMPO(sites);
    auto mz = itensor::AutoMPO(sites);

    for (int i = aux+1; i <= LX*LY+aux; i++){
        mx += "Sx",i;
        my += "Sy",i;
        mz += "Sz",i;
    }
    std::array<itensor::MPO,3> arr = {itensor::toMPO(mx),itensor::toMPO(my),itensor::toMPO(mz)};
    return arr;
}







std::vector<std::array<double,2>> Kitaev_Model::Mean(std::vector<std::vector<double>>& M){
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



void Kitaev_Model::tdvp_loop(std::vector<double>& E_vec, std::vector<double>& C_vec, std::vector<double>& S_vec, std::vector<double>& W_vec, std::array<std::vector<double>,3>& Chi_vec, itensor::MPS& psi, itensor::Cplx& t, int TimeSteps, itensor::Args& args, itensor::Sweeps& Sweeps, double& cb){
    std::complex<double> tcompl2 = t;
    double t_beta = std::real(tcompl2) * 2;

    for (int j = 0; j != TimeSteps; j++){
        double cb_old = cb;
        if (cb_old == 0){
            cb_old = 1;
        }
        cb += t_beta;
        std::cout << cb <<  " " << cb_old << std::flush;

        double E = itensor::tdvp(psi,H0,t,Sweeps,args);

        std::complex<double> w = itensor::innerC(psi,H_flux,psi);
        double c = cb * cb * (E_vec.back() - E) / t_beta;
        double s = S_vec.back() + t_beta * 0.5 * (c/cb + C_vec.back()/cb_old);

        for (int i = 0; i != Chi_vec.size(); i++){
            std::complex<double> m = itensor::innerC(psi,M[i],psi);
            double mr = std::real(m);
            std::complex<double> m2 = itensor::innerC(psi,M[i],M[i],psi);
            double m2r = std::real(m2);
            double chi = cb  * (m2r - mr*mr);
            Chi_vec[i].push_back(chi);
        }        

        E_vec.push_back(E);
        C_vec.push_back(c);
        S_vec.push_back(s);
        W_vec.push_back(std::real(w)); 
    }
}






template<std::size_t n, typename T>
void Kitaev_Model::save_data(std::string filename, std::vector<std::array<T,n>>& vec){
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
void Kitaev_Model::save_data(std::string filename, T v){
    std::vector<std::array<T,1>> vec;
    std::array<T,1> v_array = {v};
    vec.push_back(v_array);

    save_data(filename, vec);
}







void Kitaev_Model::Time_Evolution(int TimeSteps, std::vector<double> intervals, int Evols, int max_sites, int init_rand_sites, std::string TDVP_Type){
    Calc_Type = false;
    
    std::vector<std::vector<double>> Energies;
    Energies.reserve(Evols);
    std::vector<std::vector<double>> Flux;

    std::vector<double> E_vec;
    E_vec.reserve((TimeSteps) * intervals.size() + 1);

    std::vector<std::vector<double>> Capacity;
    std::vector<double> C_vec;
    std::vector<std::vector<double>> Entropy;
    std::vector<double> S_vec;
    Capacity.reserve(Evols);
    Entropy.reserve(Evols);
    C_vec.reserve(TimeSteps * intervals.size() + 1);
    S_vec.reserve(TimeSteps * intervals.size() + 1);

    std::vector<double> W_vec;    
    Flux.reserve(Evols);
    W_vec.reserve((TimeSteps) * intervals.size() + 1);

    std::array<std::vector<std::vector<double>>,3> Suscept;
    for (auto& i : Suscept){
        i.reserve(Evols);
    }
    std::array<std::vector<double>,3> Suscept_vec;
    for (auto& i : Suscept_vec){
        i.reserve(TimeSteps * intervals.size() + 1);
    }

    std::cout << "Intervals: ";
    for (auto& i : intervals){
        std::cout << i << " ";
    }
    std::cout << "\n\n";
    std::cout << "Time Steps per Interval: " << TimeSteps << "\n";
    std::cout << "Initial Random State Bond Dimension: " << init_rand_sites << "\n";


    auto Sweeps = itensor::Sweeps(1);
    itensor::Args tdvp_args;
    if (TDVP_Type == "TwoSite") {
        tdvp_args = itensor::Args({"Silent",true,"ErrGoal",1E-7});
    } else {
        tdvp_args = itensor::Args({"Silent",true,"ErrGoal",1E-7,"NumCenter",1});
    }
    Sweeps.maxdim() = max_sites;
    Sweeps.cutoff() = 1e-10;
    Sweeps.mindim() = init_rand_sites;
    Sweeps.niter() = 30;

    std::vector<itensor::Cplx> T;
    for (auto & it : intervals){
        itensor::Cplx t = -0.5 * it / static_cast<double>(TimeSteps) * itensor::Cplx_1;
        T.emplace_back(t);
    }
    auto t0 = std::chrono::system_clock::now();
    

    for (int i = 0; i != Evols; i++){
        double curr_beta = 0;
        auto t1 = std::chrono::system_clock::now();
        auto psi = itensor::randomMPS(sites,init_rand_sites);

        std::complex<double> E = itensor::innerC(psi,H0,psi) / itensor::inner(psi,psi);
        E_vec.push_back(std::real(E));
        W_vec.push_back(0);
        C_vec.push_back(0);
        S_vec.push_back(0);
        for (auto& i : Suscept_vec){
            i.push_back(0);
        }
        
        for (auto t = T.begin(); t != T.end(); t++){
            tdvp_loop(E_vec,C_vec,S_vec,W_vec,Suscept_vec,psi,*t,TimeSteps,tdvp_args,Sweeps,curr_beta);
        }

        Energies.push_back(E_vec);
        E_vec.clear();
        Capacity.push_back(C_vec);
        C_vec.clear();
        Entropy.push_back(S_vec);
        S_vec.clear();
        Flux.push_back(W_vec);
        W_vec.clear();

        for (int i = 0; i != Suscept.size(); i++){
            Suscept[i].push_back(Suscept_vec[i]);
            Suscept_vec[i].clear();
        }
        

        auto t2 = std::chrono::system_clock::now();
        auto time = std::chrono::duration<double>(t2-t1);
        std::cout << "Finished Evolution Number " << (i+1) << "/" << Evols << ", Time Needed: " << time.count() << " Seconds\n" << std::flush;


    }

    E = Mean(Energies);
    Cv = Mean(Capacity);
    S = Mean(Entropy);
    W = Mean(Flux);
    Chix = Mean(Suscept[0]);
    Chiy = Mean(Suscept[1]);
    Chiz = Mean(Suscept[2]);

    auto t3 = std::chrono::system_clock::now();
    auto time_total = std::chrono::duration<double>(t3-t0);

    auto hours = std::chrono::duration_cast<std::chrono::hours>(time_total);
    auto minutes = std::chrono::duration_cast<std::chrono::minutes>(time_total-hours);
    auto seconds = std::chrono::duration_cast<std::chrono::seconds>(time_total-hours-minutes);
    std::cout << "Finished Imaginary Time Evolution, Time Needed: " << hours.count() << " Hours, " << minutes.count() << " Minutes, " << seconds.count() << " Seconds\n" << std::flush;
}




}





