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
#include <chrono> // timings

using namespace itensor;

// model: 'ising' or 'heisenberg'
// filename: file to save resulting energies to
// N: spin chain length
void find_excited_states(std::string model="ising", int N=32, float spin=1.0, int states=3){

	// file to save resulting energies to
	std::string filename = model + std::to_string(N);
	std::string spinname = "spinHalf";
	std::ofstream energies("./logs_" + spinname + "/" + filename + "_energies" + ".csv");
	energies << "#J,E0,E1,E2" <<std::endl;
	std::ofstream timings("./logs_" + spinname + "/" + filename + "_time" + ".csv");
	timings << "#J,E0,E1,E2" <<std::endl;

	// Initialize the site degrees of freedom.
	// TODO change by hand if you want spin 1 / spin 1/2
	auto lattice = SpinHalf(N,{"ConserveQNs=", false}); //make a chain of N spin 0.5's
	//auto lattice = SpinOne(N,{"ConserveQNs=", false}); //make a chain of N spin 1's
	
	// Set the parameters controlling the accuracy of the DMRG
	// calculation for each DMRG sweep.
	auto sweeps = Sweeps(30);
	sweeps.maxdim() = 10,20,100,100,200,300;
	sweeps.cutoff() = 1E-10;
	sweeps.niter() = 2;
	sweeps.noise() = 1E-7,1E-8,0.0;
	//println(sweeps);
		
	// model parameters
	Real h = 1.0;
	Real J_max = 2.0;
	Real step_size = 0.1;

	for (Real J = -J_max; J < J_max+step_size; J+=step_size){
		auto start = std::chrono::high_resolution_clock::now();
		println("\n----------------------\n");
		println(filename + " at J = ", J);

		// Use the AutoMPO feature to create the model.
		auto ampo = AutoMPO(lattice);
		if (model == "ising"){
			// Model: Ising (Spin +-1)
			// H = - h Sx - J SzSz
			// J = -2, ..., 2 | h = 1
			for(int i = 1; i < N; ++i) {
				ampo += -4*J, "Sz", i, "Sz", i + 1;
			}
			for(int i = 1; i <= N; ++i) {
				ampo += -2*h,"Sx",i;
			}
		} else if (model == "heisenberg") {
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

		// containers to save energies to
		std::vector<float> energies_states;
		energies_states.push_back(J);

		// Begin the DMRG calculation
		// for the ground state
		auto [en0, psi0] = dmrg(H,randomMPS(lattice),sweeps,{"Quiet=",true});
		energies_states.push_back(en0);

		// Make a vector of previous wavefunctions;
		// code will penalize future wavefunctions
		// for having any overlap with these
		auto wfs = std::vector<MPS>(1);
		wfs.at(0) = psi0;

		// excited states
		for (int state = 1; state < states; ++state){
			auto [en,psi] = dmrg(H,wfs,randomMPS(lattice),sweeps,{"Quiet=",true,"Weight=",20.0});
			wfs.push_back(psi);
			//wfs.at(state) = psi;
			energies_states.push_back(en);
		}

		auto stop = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
		timings <<J <<"," << duration.count() << std::endl;
		// Save timings to log file
		// printfln("Overlap <psi1|psi2> = %.2E",inner(psi1,psi2));
		// printfln("\nOverlap <psi0|psi1> = %.2E",innerC(psi0,psi1))

		// Save energies to log file
		for (int e = 0; e < states; ++e){
			energies <<energies_states[e] <<",";
		}
		energies <<energies_states[states];
		energies <<std::endl;
	}
	energies.close();
}

int main(){

	// run some examples
	// don't forget to change lattice to spinOne or spinHalf!
	//find_excited_states("ising", 4,  0.5, 3);
	//find_excited_states("ising", 8,  0.5, 3);
	find_excited_states("ising", 16, 0.5, 3);
	//find_excited_states("ising", 32, 0.5, 3);
	//find_excited_states("ising", 64, 0.5, 3);
	//find_excited_states("ising", 128, 0.5, 3);
	
	
	
	
	
	
	
	
	

	return 0;
}
