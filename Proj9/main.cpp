#include "TDVP-MPS/main.h"

using namespace TDVP_MPS;

int main(){
char env[] = "MKL_NUM_THREADS=4";
putenv(env);

int LX = 2;
int LY = 2;

Hamiltonian H_Details;
H_Details.set("K",1);

int S2 = 1;

int aux = 2;
int aux2 = 2;

std::string shape = "HoneycombPeriodic";

std::vector<int> timesteps = {10,10,10,30,20,40};
std::vector<double> intervals = {3,5,12,80,100,300};


int Evols = 2;

std::string foldername = "Evol";

auto Model = Kitaev_Model(LX,LY,H_Details,S2,aux,aux2,shape);
Model.TPQ_MPS(timesteps,intervals,Evols);
Model.Save(foldername);


}
