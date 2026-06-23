//mot_mult.cpp
#include "hopf_algebroid.h"
#include "mot_steenrod.h"
#include"tao_bockstein.h"
#include "matrices_stream.h"

int main(int agrc, char **argv){
	//get the maximal degree and the length of the resolution
	int max_deg = std::atoi(argv[1]);
	int resolution_length = std::atoi(argv[2]);
	string directory(argv[1]);
	directory += "_";
	string director = directory;
	
	//initialize matrix class
	matrix<tauPoly>::moduleOper = &tau_module_oper;
	matrix<motSteenrod>::moduleOper = &motSteenrod_module_oper;
	//the matrix of the coactions
	matrix_mem<motSteenrod> coa;
	//the operator for motivic dual steenrid algebra
	MotSteenrodOp MOP(&coa, max_deg);
	//initialize the list of monomials
	MOP.init_mon_array(directory + "ex2poly_index");
	//the complex of the primitives
	motComplex Complex;
	//load the maps
	Complex.load(resolution_length, director + "mot_gens", director + "mot_res");
	
	//computet the tau-Bockstein
	auto tables = make_table(Complex);
	//output the tables
	std::fstream taufile(director + "tau_bockstein.txt", std::ios::out);
	taufile << output_tables(tables);
	
	//load the new generators with diagonize the boundary maps
	std::ofstream ngf(director + "new_generators");
	std::vector<cycle_data> cycles_table(tables.size());
	make_cycle_tables(tables,cycles_table);
	for(int i=0;i<tables.size();++i){
		ngf << output(cycles_table[i], i);
	}
	
	//compute the h_n multiplication for every n whose degree 2^n fits both
	//within the computed resolution (max_deg) and within xi_1's packed
	//exponent range (xnMaxExpo[1]), beyond which h_n would alias to a
	//lower class instead of failing cleanly
	int max_n = 0;
	while((1 << (max_n+1)) <= max_deg && (unsigned)(1 << (max_n+1)) < xnMaxExpo[1])
		++max_n;
	for(int n=0; n<=max_n; ++n)
		multiplication_table(MOP.hi(n), 1<<n, director + "mot_gens", director + "mot_res", director + "h" + std::to_string(n) + ".txt", resolution_length, MOP, cycles_table);
}