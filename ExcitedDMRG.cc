/* May 2022
 * Author: Andreas Burger
 * With generous help of Heitor Peres Casagrande
 *
 * Find excited states using DMRG and the ITensor library
 * Examples: 1d Ising and 1d Heisenberg
 * Spin +-1
 *
 * Ising: 			H = - h Sx - J SzSz | J = -2, ..., 2 | h = 1
 * Heisenberg: 	H = - h (SxSx + SySy) - J SzSz | J = -2, ..., 2 | h = 1
 *
 */

#include "itensor/all.h"
#include <bits/stdc++.h> // strings
#include <fstream> // write to file

using namespace itensor;

// ising: True: Ising. False: Heisenberg
// filename: file to save resulting energies to
// N: spin chain length
void find_excited_states(bool ising=true, int N=32, std::string filename="ising32"){

	// file to save resulting energies to
	std::ofstream energies("./logs/" + filename + "_energies" + ".csv");
	//energies.open ("./logs/" + filename + "_energies" + ".csv", std::ios_base::app);
	energies << "#J,E0,E1,E2" <<std::endl;

	// Initialize the site degrees of freedom.
	auto sites = SpinOne(N,{"ConserveQNs=", false}); //make a chain of N spin 1's

	// Fields
	Real h = 1.0;
	Real J_max = 2.0;
	Real step_size = 0.1;
	for (Real J = -J_max; J < J_max+step_size; J+=step_size){
		// Use the AutoMPO feature to create the model.
		auto ampo = AutoMPO(sites);
		if (ising){
			// Model: Ising (Spin +-1)
			// H = - h Sx - J SzSz
			// J = -2, ..., 2 | h = 1
			for(int i = 1; i < N; ++i) {
				ampo += -J, "Sz", i, "Sz", i + 1;
			}
			for(int i = 1; i <= N; ++i) {
				ampo += -h,"Sx",i;
			}
		} else {
			// Model: Heisenberg (Spin +-1)
			// H = - h (SxSx + SySy) - J SzSz
			// J = -2, ..., 2 | h = 1
			for(int i = 1; i < N; ++i){
				ampo += -J,"Sz",i,"Sz",i+1;
			}
			for(int i = 1; i < N; ++i){
				ampo += -h,"Sx",i,"Sx",i+1;
			}
			for(int i = 1; i < N; ++i){
				ampo += -h,"Sy",i,"Sy",i+1;
			}
		}
		auto H = toMPO(ampo);

		// Set the parameters controlling the accuracy of the DMRG
		// calculation for each DMRG sweep.
		auto sweeps = Sweeps(30);
		sweeps.maxdim() = 10,20,100,100,200;
		sweeps.cutoff() = 1E-10;
		sweeps.niter() = 2;
		sweeps.noise() = 1E-7,1E-8,0.0;
		println(sweeps);

		// Begin the DMRG calculation
		// for the ground state
		auto [en0,psi0] = dmrg(H,randomMPS(sites),sweeps,{"Quiet=",true});

		println("\n----------------------\n");
		println(filename + " at J = ", J);

		// Make a vector of previous wavefunctions;
		// code will penalize future wavefunctions
		// for having any overlap with these
		auto wfs = std::vector<MPS>(1);
		wfs.at(0) = psi0;

		// First excited state
		// Weight option sets the energy penalty for psi1 having any overlap with psi0
		auto [en1,psi1] = dmrg(H,wfs,randomMPS(sites),sweeps,{"Quiet=",true,"Weight=",20.0});

		// wavefunctions of first ground state and first excited state
		auto wfs2 = std::vector<MPS>(2);
		wfs2.at(0) = psi0;
		wfs2.at(1) = psi1;

		// Second excited state
		auto [en2,psi2] = dmrg(H, wfs2, randomMPS(sites), sweeps,{"Quiet=",true,"Weight=",20.0});

		// Print the final energies reported by DMRG
		printfln("\nGround State Energy = %.10f",en0);
		printfln("Excited State 1 Energy = %.10f",en1);
		printfln("Excited State 2 Energy = %.10f",en2);

		// The expected gap of the transverse field Ising
		// model is given by Eg = 2*|h-1|
		// (The DMRG gap will have finite-size corrections.)
		printfln("\nDMRG energy gap ex1-gs0 = %.10f",en1-en0);
		printfln("DMRG energy gap ex2-gs0 = %.10f",en2-en0);

		// The overlap <psi0|psi1> should be very close to zero
		if (ising){
			printfln("\nOverlap <psi0|psi1> = %.2E",inner(psi0,psi1));
			printfln("Overlap <psi0|psi2> = %.2E",inner(psi0,psi2));
			printfln("Overlap <psi1|psi2> = %.2E",inner(psi1,psi2));
		} else {
			printfln("\nOverlap <psi0|psi1> = %.2E",innerC(psi0,psi1));
			printfln("Overlap <psi0|psi2> = %.2E",innerC(psi0,psi2));
			printfln("Overlap <psi1|psi2> = %.2E",innerC(psi1,psi2));
		}

		// Save energies log file
		energies <<J <<"," <<en0 <<"," <<en1 <<"," <<en2 <<std::endl;
	}
	energies.close();
}

int main(){

	// run some examples
	find_excited_states(true, 16, "ising16");
	//find_excited_states(true, 32, "ising32");
	//find_excited_states(false, 16, "heisenberg16");
	//find_excited_states(false, 32, "heisenberg32");
	//find_excited_states(true, 16, "ising16_u1");
	//find_excited_states(false, 16, "heisenberg16_u1");

	return 0;
}