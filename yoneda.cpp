//yoneda.cpp
//General Yoneda products for the C-motivic minimal resolution, via chain-map lifting.
//
//This generalizes mot_mult.cpp (which handles only multiplication by the
//filtration-1 classes h_n) to multiplication by an arbitrary Ext class beta,
//by lifting the chain map that represents beta up through the resolution.
//
//Bootstrap (see README / memory): e2p -> motTab -> mr_ex -> mr_mot must be run
//first, producing <md>_mot_gens, <md>_mot_res, <md>_mot_maps{i}, <md>_mot_gens{i}.

#include "hopf_algebroid.h"
#include "mot_steenrod.h"
#include "matrices_mem.h"
#include "tao_bockstein.h"
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>

//----------------------------------------------------------------------------
// Per-step resolution data, loaded from disk.
//----------------------------------------------------------------------------
struct ResStep {
	FreeMotCoMod F;                 // the cofree comodule F_i (cogenerators = generators of C_i)
	matrix_mem<tauPoly> inj;        // inj_i : X_i -> F_i  (rows indexed by X_i basis)
	matrix_mem<tauPoly> qut;        // qut_i : F_i -> X_{i+1}
	int Xrank;                      // rank of X_i (= inj_i.rank = number of rows)
};

//load one cofree comodule from the per-step generators file <pre>mot_gens{i}
static void load_F(const string &pre, int i, FreeMotCoMod &F){
	std::fstream f(pre + "mot_gens" + std::to_string(i), std::ios::in | std::ios::binary);
	if(!f.is_open()){ std::cerr << "cannot open mot_gens" << i << "\n"; std::exit(1); }
	F.load(f);
}

//load inj_i and qut_i from <pre>mot_maps{i}
static void load_maps(const string &pre, int i, matrix_mem<tauPoly> &inj, matrix_mem<tauPoly> &qut){
	std::fstream f(pre + "mot_maps" + std::to_string(i), std::ios::in | std::ios::binary);
	if(!f.is_open()){ std::cerr << "cannot open mot_maps" << i << "\n"; std::exit(1); }
	inj.clear(); qut.clear();
	inj.load(f);
	qut.load(f);
}

typedef vectors<matrix_index, tauPoly> TVec;

//tau-aware (Hermite-style) echelon decomposition over F2[tau].
//Absorbs rows (each a "row" vector in some column space, paired with a "tag"
//vector recording the linear combination of original generators that produced
//it).  Pivots are kept with MINIMAL tau-power at each leading position, so that
//any element of the F2[tau]-span reduces to zero.
//
//Used two ways:
//  * invert inj_i on its image:  absorb (inj.find(j), e_j); then solve(w) yields
//    the unique u with inj_i(u)=w for w in the image.
//  * section qut_i:             absorb (qut.find(p), e_p); then solve(y) yields a
//    v with qut_i(v)=y (qut surjective => always solvable).
struct TauEchelon {
	//position -> (row, tag); row.dataArray[0].ind == position (leading), minimal tau power
	std::map<matrix_index, std::pair<TVec,TVec>> piv;

	static inline tauPoly lead_pow(const TVec &v){ return v.dataArray[0].coeficient; }

	void absorb(TVec row, TVec tag){
		while(!tau_module_oper.isZero(row)){
			matrix_index p = row.dataArray[0].ind;
			tauPoly a = lead_pow(row);
			auto it = piv.find(p);
			if(it == piv.end()){ piv.emplace(p, std::make_pair(std::move(row), std::move(tag))); return; }
			tauPoly b = lead_pow(it->second.first);
			if(a >= b){
				tauPoly f = a - b;   // tau^{a-b}; char 2 so subtract == add
				row = tau_module_oper.add(row, tau_module_oper.scalor_mult(f, it->second.first));
				tag = tau_module_oper.add(tag, tau_module_oper.scalor_mult(f, it->second.second));
			} else {
				//new row has the smaller tau power: it becomes the pivot, reprocess old
				std::pair<TVec,TVec> old = std::move(it->second);
				it->second = std::make_pair(std::move(row), std::move(tag));
				row = std::move(old.first); tag = std::move(old.second);
			}
		}
		//row reduced to zero: a relation; nothing to record
	}

	//reduce w using the pivots; on success residual==0 and u holds the preimage
	//(combination of tags) with image(u) == w.  Returns false if w is not in span.
	bool solve(TVec w, TVec &u) const {
		u.clear();
		while(!tau_module_oper.isZero(w)){
			matrix_index p = w.dataArray[0].ind;
			tauPoly a = lead_pow(w);
			auto it = piv.find(p);
			if(it == piv.end()) return false;
			tauPoly b = lead_pow(it->second.first);
			if(a < b) return false;          // would need a negative tau power
			tauPoly f = a - b;
			w = tau_module_oper.add(w, tau_module_oper.scalor_mult(f, it->second.first));
			u = tau_module_oper.add(u, tau_module_oper.scalor_mult(f, it->second.second));
		}
		return true;
	}
};

//Build a TauEchelon from the rows of a matrix M (row j tagged by e_j), giving a
//solver for "M^T-style" preimages: solve(w) finds u with sum_j u_j M.find(j) == w.
static void build_echelon(matrix_mem<tauPoly> &M, TauEchelon &ech){
	for(unsigned j=0; j<M.rank; ++j)
		ech.absorb(M.find(j), tau_module_oper.singleton(j));
}

//defined in tao_bockstein.cpp: express a cogenerator-combination in terms of the
//chosen cycle representatives (the tau-Bockstein basis of Ext).
vectors<matrix_index,tauPoly> find_cycle(cycle_data &table, vectors<matrix_index,tauPoly> v);

//---- chain-map lifting helpers -------------------------------------------------

//comodule coaction of a general (tauPoly) element e of cofree comodule F,
//extended base-linearly from the coaction on basis elements.
static vectors<matrix_index,motSteenrod> coaction_of(const TVec &e, FreeMotCoMod &F, MotSteenrodOp &MOP){
	vectors<matrix_index,motSteenrod> res;
	for(auto &tm : e.dataArray){
		auto cap = F.coaction(tm.ind);
		res = motSteenrod_module_oper.add(res,
			motSteenrod_module_oper.scalor_mult(MOP.etaL(tm.coeficient), cap));
	}
	return res;
}

//Apply the comodule map  f = adjoint(psi) : Fsrc -> Ftgt  to an element e, where
//psi : Fsrc -> cogens(Ftgt) is a base map.  By the cofree adjunction
//   f(e) = (id_A (x) psi)( coaction_{Fsrc}(e) )  in Ftgt.
static TVec apply_chainmap(const TVec &e, FreeMotCoMod &Fsrc, matrix_mem<tauPoly> &psi,
                           FreeMotCoMod &Ftgt, MotSteenrodOp &MOP){
	auto coa = coaction_of(e, Fsrc, MOP);          // sum_k A_k [k]
	unsigned ngens = Ftgt.position_of_gens.size();
	TVec res;
	for(auto &tm : coa.dataArray){
		motSteenrod A = tm.coeficient;
		auto psik = psi.find(tm.ind);              // psi([k]) : tauPoly vec over cogens(Ftgt)
		for(auto &ct : psik.dataArray){            // (c, j)
			matrix_index pos = Ftgt.position_of_gens[ct.ind];
			//the cofree summand of cogenerator ct.ind occupies [pos, end); terms placed
			//beyond it would be degree-overflow (product degree > max_deg) and must be
			//dropped, both for correctness and to avoid running off the module.
			matrix_index end = (ct.ind + 1 < ngens) ? Ftgt.position_of_gens[ct.ind+1]
			                                        : (matrix_index)Ftgt.total_rank;
			motSteenrod Ac = motSteenrod_oper.multiply(A, MOP.etaR(ct.coeficient));
			auto contrib = MOP.algebroid2vector(Ac, pos);
			std::function<bool(matrix_index)> in_summand = [end](matrix_index n){ return n < end; };
			res = tau_module_oper.add(res, tau_module_oper.filter(in_summand, contrib));
		}
	}
	return res;
}

//project a cofree-comodule element to its cogenerators (tauPoly vec over cogen index)
static TVec proj_cogens(const TVec &v, FreeMotCoMod &F){
	std::function<matrix_index(matrix_index)> rule = [&F](matrix_index n){
		int h = F.find_index(n);
		return (h == FreeMotCoMod::invalid_pos) ? (matrix_index)-1 : (matrix_index)h; };
	return tau_module_oper.filtered_reindex(rule, v);
}

//format a cogenerator-combination at filtration nm, like mot_mult's output_cycles
static string fmt(int nm, const TVec &v){
	if(v.size() == 0) return "0";
	string res;
	for(auto &tm : v.dataArray){
		if(!res.empty()) res += "+";
		res += "t^" + std::to_string(tm.coeficient) + "{" + std::to_string(nm) + "-" + std::to_string(tm.ind) + "}";
	}
	return res;
}

//----------------------------------------------------------------------------
int main(int argc, char **argv){
	if(argc < 3){
		std::cerr << "usage: " << argv[0] << " <max_deg> <resolution_length>\n";
		return 1;
	}
	int max_deg = std::atoi(argv[1]);
	int resolution_length = std::atoi(argv[2]);
	string pre = string(argv[1]) + "_";

	//module operations
	matrix<tauPoly>::moduleOper = &tau_module_oper;
	matrix<motSteenrod>::moduleOper = &motSteenrod_module_oper;
	curtis_table<tauPoly>::ModOper = &tau_module_oper;

	//motivic dual Steenrod algebra (needed for coactions used later in lifting)
	matrix_mem<motSteenrod> coa;
	MotSteenrodOp MOP(&coa, max_deg);
	MOP.init_mon_array(pre + "ex2poly_index");
	MOP.generate_cofree_coaction(pre + "mot_deltas", pre + "poly_exponents");

	//load the resolution, one step at a time
	std::vector<ResStep> R(resolution_length + 1);
	std::vector<TauEchelon> invInj(resolution_length + 1);  // invInj[i]: invert inj_i on its image
	std::vector<TauEchelon> secQut(resolution_length + 1);  // secQut[i]: section of qut_i
	for(int i=0; i<=resolution_length; ++i){
		load_F(pre, i, R[i].F);
		load_maps(pre, i, R[i].inj, R[i].qut);
		R[i].Xrank = R[i].inj.rank;
		build_echelon(R[i].inj, invInj[i]);
		build_echelon(R[i].qut, secQut[i]);
	}

	//----- solver self-check mode -----
	if(argc == 3){
		bool ok = true;
		for(int i=0; i<=resolution_length; ++i){
			for(unsigned j=0; j<R[i].inj.rank; ++j){
				TVec u;
				bool solved = invInj[i].solve(R[i].inj.find(j), u);
				auto diff = tau_module_oper.add(u, tau_module_oper.minus(tau_module_oper.singleton(j)));
				if(!solved || !tau_module_oper.isZero(diff)){ ok=false;
					std::cerr << "invInj FAIL step " << i << " row " << j << "\n"; }
			}
			if(i < resolution_length)
				for(int y=0; y<R[i+1].Xrank; ++y){
					TVec v; bool solved = secQut[i].solve(tau_module_oper.singleton(y), v);
					auto diff = tau_module_oper.add(R[i].qut.maps_to(v), tau_module_oper.minus(tau_module_oper.singleton(y)));
					if(!solved || !tau_module_oper.isZero(diff)){ ok=false;
						std::cerr << "secQut FAIL step " << i << " y " << y << "\n"; }
				}
		}
		std::cout << (ok ? "SOLVER CHECK: success\n" : "SOLVER CHECK: FAIL\n");
		return ok ? 0 : 1;
	}

	//----- build the tau-Bockstein cycle representatives (the basis of Ext) -----
	//(the resolution/tau-Bockstein routines print progress to stdout; silence it so
	// the product output is clean)
	std::streambuf *saved = std::cout.rdbuf();
	std::ofstream devnull("/dev/null");
	std::cout.rdbuf(devnull.rdbuf());
	motComplex Complex;
	Complex.load(resolution_length, pre + "mot_gens", pre + "mot_res");
	auto tables = make_table(Complex);
	std::vector<cycle_data> cyc(tables.size());
	make_cycle_tables(tables, cyc);
	std::cout.rdbuf(saved);

	//----- product modes -----
	// argc==5 : per-generator table for beta=(bs,bb)
	// argc==7 : single pair  alpha=(as,aa) * beta=(bs,bb)
	int bs, bb;
	if(argc == 5){ bs = std::atoi(argv[3]); bb = std::atoi(argv[4]); }
	else if(argc == 7){ bs = std::atoi(argv[5]); bb = std::atoi(argv[6]); }
	else { std::cerr << "usage: " << argv[0] << " md len            (solver check)\n"
	                 << "       " << argv[0] << " md len bs bb       (table for beta=(bs,bb))\n"
	                 << "       " << argv[0] << " md len as aa bs bb (alpha*beta)\n"; return 1; }

	if(bs < 0 || bs > resolution_length){ std::cerr << "beta filtration out of range\n"; return 1; }

	//lift the chain map for beta=(bs,bb): M[i] : cogens(F_i) -> cogens(F_{i+bs})
	int maxlev = resolution_length - bs;
	std::vector<matrix_mem<tauPoly>> psi(maxlev + 1);
	psi[0].set_rank(R[0].F.total_rank);
	int cog0 = R[0].F.position_of_gens[0];
	for(int p=0; p<R[0].F.total_rank; ++p)
		psi[0].insert(p, (p==cog0) ? tau_module_oper.singleton(bb) : tau_module_oper.zero());
	for(int i=0; i<maxlev; ++i){
		psi[i+1].set_rank(R[i+1].F.total_rank);
		for(int x=0; x<R[i+1].F.total_rank; ++x){
			TVec u;
			if(invInj[i+1].solve(tau_module_oper.singleton(x), u)){
				TVec v; secQut[i].solve(u, v);
				TVec fv  = apply_chainmap(v, R[i].F, psi[i], R[i+bs].F, MOP);
				TVec dfv = R[i+bs+1].inj.maps_to( R[i+bs].qut.maps_to(fv) );
				psi[i+1].insert(x, proj_cogens(dfv, R[i+bs+1].F));
			} else
				psi[i+1].insert(x, tau_module_oper.zero());
		}
	}
	//chain map on cogenerators: M[i].find(a) = proj_cogens(f_i(cogenerator a))
	std::vector<matrix_mem<tauPoly>> M(maxlev + 1);
	for(int i=0; i<=maxlev; ++i){
		M[i].set_rank(R[i].F.generators.rank);
		for(unsigned a=0; a<R[i].F.generators.rank; ++a){
			TVec cog = tau_module_oper.singleton(R[i].F.position_of_gens[a]);
			TVec fa  = apply_chainmap(cog, R[i].F, psi[i], R[i+bs].F, MOP);
			M[i].insert(a, proj_cogens(fa, R[i+bs].F));
		}
	}

	//product of class (s,a) with beta: name M[s](cycle_rep(s,a)) in filtration s+bs
	auto product = [&](int s, int a) -> TVec {
		TVec rep = (cyc[s].count(a)) ? cyc[s].at(a) : tau_module_oper.singleton(a);
		TVec img = M[s].maps_to(rep);
		return find_cycle(cyc[s+bs], img);
	};

	if(argc == 7){
		int as = std::atoi(argv[3]), aa = std::atoi(argv[4]);
		if(as < 0 || as > maxlev){ std::cerr << "alpha filtration out of lifted range (<= " << maxlev << ")\n"; return 1; }
		TVec res = product(as, aa);
		std::cout << "{" << as << "-" << aa << "} * {" << bs << "-" << bb << "}  =  "
		          << fmt(as+bs, res) << "\n";
	} else {
		for(int i=0; i<=maxlev; ++i)
			for(unsigned a=0; a<R[i].F.generators.rank; ++a)
				std::cout << "{" << i << "-" << a << "}\t->\t" << fmt(i+bs, product(i,a)) << "\n";
	}
	return 0;
}
