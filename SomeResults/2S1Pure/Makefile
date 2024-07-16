
CC=mpicxx -std=c++17 -fPIC

COPT= -I${MKLROOT}/include -I ~/itensor -I ~ -I ~/eigen-3.4.0

LOPT= -lstdc++fs  -L/opt/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_rt -lmkl_core -lpthread -L/opt/intel/oneapi/compiler/2021.2.0/linux/compiler/lib/intel64_lin -liomp5 -L ~/itensor-mt/lib -litensor

tdvp: Makefile tdvp.cpp
	$(CC) tdvp.cpp  $(LOPT) $(COPT) -o $@ 
