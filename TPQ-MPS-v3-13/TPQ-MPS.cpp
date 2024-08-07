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


std::array<int,3> Kitaev_Model::get_neighbour_data_hex_periodic(int LX, int LY, int pos){
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

} 


std::array<int,3> Kitaev_Model::get_neighbour_data_hex_rev(int LX, int LY, int pos){
    int Y = (pos-1) % LY;
    int X = (pos-1) / LY;
    std::array<int,3> neighbours{};

    if (X != LX-1){
        neighbours[0] = pos+LY+1;
    }
    else{
        neighbours[0] = pos-(LX-1)*LY+1;
    }

    neighbours[1] = pos+1;

    if (Y != 0){
        neighbours[2] = pos-1;
    }
    else{
        neighbours[2] = pos+LY-1;
    }
    return neighbours;

}


std::array<int,3> Kitaev_Model::get_neighbour_data_hex_rev2(int LX, int LY, int pos){
    int Y = LY - 1 - (pos-1) % LY;
    int X = (pos-1) / LY;
    std::array<int,3> neighbours{};

    neighbours[0] = pos+1;

    if (X != 0){
        neighbours[1] = pos-LY+1;
    }
    else {
        neighbours[1] = pos+(LX-1)*LY+1;
    }  
    if (Y != 0){
        neighbours[2] = pos-1;
    }
    else {
        neighbours[2] = pos+LY-1;
    } 
    return neighbours;

}





std::array<int,3> Kitaev_Model::get_neighbour_data_tri(int LX, int LY, int pos){
    int X = (pos-1) / LY;
    int Y = (pos-1) % LY;
    std::array<int,3> neighbours{};

    if (X != (LX-1)){
        neighbours[0] = pos + LY;
        if (Y == 0){
            neighbours[2] = pos + 2*LY - 1;
        }
        else{
            neighbours[2] = pos + LY - 1;
        }
    }
    else{
        neighbours[0] = pos - (LX-1) * LY;
        if (Y == 0){
            neighbours[2] = pos - (LX-2) * LY - 1;
        }
        else{
            neighbours[2] = pos - (LX-1) * LY - 1;
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



itensor::MPO Kitaev_Model::honeycomb_flux_operator_half(int LX, int LY, int aux){
    int Py = std::max((LY/2)-1,0);
    int Px = LX-1;
    auto fluxop = itensor::AutoMPO(sites);
    double Wfac = 64. / static_cast<double>(Px) / static_cast<double>(Py);

    for (int i = 0; i != Px; i++){
        for (int j = 0; j != Py; j++){
            int f = LY*i + 2*j + 2 + aux;
            fluxop += Wfac,"Sx",f,"Sy",f+1,"Sz",f+2,"Sx",f+LY+1,"Sy",f+LY,"Sz",f+LY-1;
        }
    }
    auto fluxH = itensor::toMPO(fluxop);
    return fluxH;
}




itensor::MPO Kitaev_Model::honeycomb_flux_operator(int LX, int LY, int aux){
    int Py = std::max((LY/2)-1,0);
    int Px = LX-1;
    auto fluxop = itensor::AutoMPO(sites);
    double Wfac = 1. / static_cast<double>(Px) / static_cast<double>(Py);

    for (int i = 0; i != Px; i++){
        for (int j = 0; j != Py; j++){
            int f = LY*i + 2*j + 2 + aux;
            fluxop += -1*Wfac,"Expx",f,"Expy",f+1,"Expz",f+2,"Expx",f+LY+1,"Expy",f+LY,"Expz",f+LY-1;
        }
    }    
    auto fluxH = itensor::toMPO(fluxop);
    return fluxH;
}





std::array<std::array<itensor::MPO,3>,2> Kitaev_Model::magnetization_operators(int LX, int LY, int aux){
    auto mx = itensor::AutoMPO(sites);
    auto my = itensor::AutoMPO(sites);
    auto mz = itensor::AutoMPO(sites);
    auto mx2 = itensor::AutoMPO(sites);
    auto my2 = itensor::AutoMPO(sites);
    auto mz2 = itensor::AutoMPO(sites);

    for (int i = aux+1; i <= LX*LY+aux; i++){
        mx += "Sx",i;
        my += "Sy",i;
        mz += "Sz",i;
        for (int j = aux+1; j <= LX*LY+aux; j++){
            mx += "Sx",i,"Sx",j;
            my += "Sy",i,"Sy",j;
            mz += "Sz",i,"Sz",j;
        }
    }
    std::array<itensor::MPO,3> m1 = {itensor::toMPO(mx),itensor::toMPO(my),itensor::toMPO(mz)};
    std::array<itensor::MPO,3> m2 = {itensor::toMPO(mx2),itensor::toMPO(my2),itensor::toMPO(mz2)};
    std::array<std::array<itensor::MPO,3>,2> m = {m1,m2};

    return m;
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



int Kitaev_Model::tdvp_loop(std::vector<double>& E_vec, std::vector<double>& C_vec, std::vector<double>& S_vec, std::vector<double>& W_vec, std::array<std::vector<double>,3>& M_vec, std::array<std::vector<double>,3>& M_vec2, itensor::MPS& psi, itensor::Cplx& t, int TimeSteps, itensor::Args& args, itensor::Sweeps& sweeps, double& cb){
    std::complex<double> tcompl2 = t;
    double t_beta = std::real(tcompl2) * -2;
    int max_bond = 0;

    for (int j = 0; j != TimeSteps; j++){
        double cb_old = cb;
        if (cb_old == 0){
            cb_old = 1;
        }
        cb += t_beta;
        double E = itensor::tdvp(psi,H0,t,sweeps,args);
        max_bond = std::max(max_bond,itensor::maxLinkDim(psi));

        std::complex<double> w = itensor::innerC(psi,H_flux,psi);
        double c = cb * cb * (E_vec.back() - E) / t_beta;
        double s = S_vec.back() + t_beta * 0.5 * (c/cb + C_vec.back()/cb_old);

        E_vec.push_back(E);
        C_vec.push_back(c);
        S_vec.push_back(s);
        W_vec.push_back(std::real(w)); 

        for (int indi = 0; indi != 3; indi++){
            std::complex<double> m = itensor::innerC(psi,M[indi],psi);
            M_vec[indi].push_back(std::real(m));
            std::complex<double> m2 = itensor::innerC(psi,M2[indi],psi);
            M_vec2[indi].push_back(std::real(m2));
        }
    }
    return max_bond;
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


void Kitaev_Model::chi_int(itensor::MPS& psi, double n, double t, std::array<std::vector<double>,3>& chi_vec, double step, itensor::Sweeps& sweeps, itensor::Args& args){
    int steps = std::ceil(t);
    itensor::Cplx dt = t / static_cast<double>(steps) * itensor::Cplx_1;
    auto appl_args = itensor::Args({"Cutoff=",1e-4,"MaxDim=",32});

    itensor::MPS sp;
    sp = itensor::applyMPO(M[0],psi,appl_args);
    itensor::println(sp);
    std::cout << 3333333 << "\n" << std::flush;
    std::complex<double> ss = itensor::innerC(sp,sp);
    std::cout << ss << "\n" << std::flush;
    sp.normalize();
    std::cout << steps << " " << dt << "\n" << std::flush;


    double E = itensor::tdvp(sp,H0,dt,sweeps,args);
    std::cout << 0 << " " << std::flush;
    


/*
    std::array<itensor::MPS*,3> Spsi;

    for (int i = 0; i != 3; i++){
        *(Spsi[i]) = itensor::applyMPO(M[i],psi);

        for (int j = 0; j != steps; j++){
            double E = itensor::tdvp(*(Spsi[i]),H0,dt,sweeps,args);
            std::cout << j << " " << std::flush;
        }
        double chip = itensor::inner(*(Spsi[i]),*(Spsi[i])) * n;
        double chiold = chi_vec[i].back();
        chi_vec[i].push_back(step/2 * (chip + chiold));
    }
    */

}




void Kitaev_Model::Time_Evolution(int TimeSteps, std::vector<double> intervals, int Evols, int max_sites, int init_rand_sites, std::string TDVP_Type){
    this->Calc_Type = false;
    
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




    std::array<std::vector<std::vector<double>>,3> Magnetization;
    for (auto& i : Magnetization){
        i.reserve(Evols);
    }
    std::array<std::vector<double>,3> Mag_vec;
    for (auto& i : Mag_vec){
        i.reserve(TimeSteps * intervals.size() + 1);
    }



    std::array<std::vector<std::vector<double>>,3> Magnetization2;
    for (auto& i : Magnetization2){
        i.reserve(Evols);
    }
    std::array<std::vector<double>,3> Mag_vec2;
    for (auto& i : Mag_vec2){
        i.reserve(TimeSteps * intervals.size() + 1);
    }



    std::cout << "Intervals: ";
    for (auto& i : intervals){
        std::cout << i << " ";
    }
    std::cout << "\n\n";
    std::cout << "Time Steps per Interval: " << TimeSteps << "\n";
    std::cout << "Initial Random State Bond Dimension: " << init_rand_sites  << "\n";
    std::cout << "Maximum Bond Dimension: " << max_sites << "\n" << std::flush;


    auto Sweeps = itensor::Sweeps(1);
    itensor::Args tdvp_args;
    if (TDVP_Type == "TwoSite") {
        tdvp_args = itensor::Args({"Silent",true,"ErrGo\na << std::flush;l",1E-7});
    } else {
        tdvp_args = itensor::Args({"Silent",true,"ErrGoal",1E-7,"NumCenter",1});
    }
    Sweeps.maxdim() = max_sites;
    Sweeps.cutoff() = 1e-10;
    Sweeps.mindim() = init_rand_sites;
    Sweeps.niter() = 30;
    std::cout << "TDVP Technique used: " << TDVP_Type << "\n";

    std::vector<itensor::Cplx> T;
    for (auto & it : intervals){
        itensor::Cplx t = -0.5 * it / static_cast<double>(TimeSteps) * itensor::Cplx_1;
        T.emplace_back(t);
    }
    auto t0 = std::chrono::system_clock::now();
    int max_bond = 0;
    
    std::cout << "\n\n";
    for (int i = 0; i != Evols; i++){
        double curr_beta = 0;
        auto t1 = std::chrono::system_clock::now();
        auto psi = itensor::randomMPS(sites,init_rand_sites);

        std::complex<double> E = itensor::innerC(psi,H0,psi) / itensor::inner(psi,psi);
        E_vec.push_back(std::real(E));
        W_vec.push_back(0);
        C_vec.push_back(0);
        S_vec.push_back(0);
        for (auto& i : Mag_vec){
            i.push_back(0);
        }
        for (auto& i : Mag_vec2){
            i.push_back(0);
        }
        //double n = 1;

        for (auto t = T.begin(); t != T.end(); t++){
            int curbond = tdvp_loop(E_vec,C_vec,S_vec,W_vec,Mag_vec,Mag_vec2,psi,*t,TimeSteps,tdvp_args,Sweeps,curr_beta);
            max_bond = std::max(max_bond,curbond);
        }
        S_vec = S_vec.back() - S_vec;

        Energies.push_back(E_vec);
        E_vec.clear();
        Capacity.push_back(C_vec);
        C_vec.clear();
        Entropy.push_back(S_vec);
        S_vec.clear();
        Flux.push_back(W_vec);
        W_vec.clear();

        for (int i = 0; i != Magnetization.size(); i++){
            Magnetization[i].push_back(Mag_vec[i]);
            Mag_vec[i].clear();
            Magnetization2[i].push_back(Mag_vec2[i]);
            Mag_vec2[i].clear();
        }
        

        auto t2 = std::chrono::system_clock::now();
        auto time = std::chrono::duration<double>(t2-t1);
        std::cout << "Finished Evolution Number " << (i+1) << "/" << Evols << ", Time Needed: " << time.count() << " Seconds\n" << std::flush;


    }

    E = Mean(Energies);
    Cv = Mean(Capacity);
    S = Mean(Entropy);
    W = Mean(Flux);

    Mx = Mean(Magnetization[0]);
    My = Mean(Magnetization[1]);
    Mz = Mean(Magnetization[2]);

    Mx2 = Mean(Magnetization2[0]);
    My2 = Mean(Magnetization2[1]);
    Mz2 = Mean(Magnetization2[2]);

    auto t3 = std::chrono::system_clock::now();
    auto time_total = std::chrono::duration<double>(t3-t0);

    auto hours = std::chrono::duration_cast<std::chrono::hours>(time_total);
    auto minutes = std::chrono::duration_cast<std::chrono::minutes>(time_total-hours);
    auto seconds = std::chrono::duration_cast<std::chrono::seconds>(time_total-hours-minutes);
    std::cout << "\n\n";
    std::cout << "Finished Imaginary Time Evolution, Time Needed: " << hours.count() << " Hours, " << minutes.count() << " Minutes, " << seconds.count() << " Seconds\n" << std::flush;

    if (max_bond == max_sites){
        std::cout << "Warning: MPS bond dimension has reached the bond dimension maximum. For more accurate results increase max_sites.\n\n" << std::flush;
    }
    else {
        std::cout << "Maximum MPS bond dimension: " << max_bond << "\n\n" << std::flush;
    }
}




}






