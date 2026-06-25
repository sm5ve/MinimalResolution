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
#include <stdexcept>

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
static void missing_step_file(const string &path, int i){
	std::cerr << "cannot open " << path << "\n"
	             "  This file is produced by mr_mot. Either mr_mot was never run for this\n"
	             "  <md>, or it was run with a resolution length <= " << i << " (this program\n"
	             "  needs per-step data up through step " << i << "). Bootstrap order:\n"
	             "      e2p <md>  ->  motTab <md>  ->  mr_ex <md> <len_ex>  ->  mr_mot <md> <len>\n"
	             "  with <len> >= " << i << ", and <md> matching the <md> passed to this program.\n";
	std::exit(1);
}

static void load_F(const string &pre, int i, FreeMotCoMod &F){
	string path = pre + "mot_gens" + std::to_string(i);
	std::fstream f(path, std::ios::in | std::ios::binary);
	if(!f.is_open()) missing_step_file(path, i);
	F.load(f);
}

//load inj_i and qut_i from <pre>mot_maps{i}
static void load_maps(const string &pre, int i, matrix_mem<tauPoly> &inj, matrix_mem<tauPoly> &qut){
	string path = pre + "mot_maps" + std::to_string(i);
	std::fstream f(path, std::ios::in | std::ios::binary);
	if(!f.is_open()) missing_step_file(path, i);
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
static void print_usage(const char *prog){
	std::cerr <<
		"Yoneda (composition) products on the C-motivic minimal resolution.\n"
		"\n"
		"Bootstrap (run once per <md>, in this order -- motTab is easy to forget,\n"
		"and mr_ex's length must be strictly greater than the length given to\n"
		"mr_mot/this program):\n"
		"    e2p <md>\n"
		"    motTab <md>\n"
		"    mr_ex <md> <len_ex>\n"
		"    mr_mot <md> <len>\n"
		"\n"
		"Usage:\n"
		"  " << prog << " <md> <len>\n"
		"      Solver self-check. Verifies that the tau-aware solvers invert\n"
		"      iota_i (injection) and qut_i (quotient map) correctly at every\n"
		"      step of the resolution loaded for max_deg=<md>, length=<len>.\n"
		"      Prints \"SOLVER CHECK: success\", or reports the failing step/row\n"
		"      and exits nonzero. No products are computed in this mode.\n"
		"\n"
		"  " << prog << " <md> <len> <bs> <bb>\n"
		"      Table mode. Let beta be the cogenerator {bs-bb} (filtration bs,\n"
		"      generator index bb of F_bs), i.e. Ext class {bs-bb}. Lifts the\n"
		"      chain map representing beta and prints alpha*beta, named in the\n"
		"      tau-Bockstein basis, for every generator alpha={s-a} with\n"
		"      s+bs <= <len> -- i.e. the full multiplication-by-beta table.\n"
		"\n"
		"  " << prog << " <md> <len> <as> <aa> <bs> <bb>\n"
		"      Single product mode. Prints just alpha*beta for the two specific\n"
		"      classes alpha={as-aa}, beta={bs-bb} (same lift as table mode,\n"
		"      but only the one entry is reported).\n"
		"\n"
		"<md> is HALF the maximal topological degree; <len> is the resolution\n"
		"length, and must not exceed the length mr_mot was run with for this <md>.\n";
}

int main(int argc, char **argv){
	if(argc < 3 || string(argv[1]) == "-h" || string(argv[1]) == "--help"){
		print_usage(argv[0]);
		return (argc < 3) ? 1 : 0;
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
	else { print_usage(argv[0]); return 1; }

	if(bs < 0 || bs > resolution_length){ std::cerr << "beta filtration out of range\n"; return 1; }

	//lift the chain map for beta=(bs,bb): M[i] : cogens(F_i) -> cogens(F_{i+bs})
	int maxlev = resolution_length - bs;
	std::vector<matrix_mem<tauPoly>> psi(maxlev + 1);
	psi[0].set_rank(R[0].F.total_rank);
	int cog0 = R[0].F.position_of_gens[0];
	// Use the tau-Bockstein cycle representative of beta={bs-bb} as the seed of the
	// chain map.  singleton(bb) would be wrong whenever the representative involves
	// a tau power or a linear combination of cogenerators.
	TVec beta_rep = (cyc[bs].count(bb)) ? cyc[bs].at(bb) : tau_module_oper.singleton(bb);
	for(int p=0; p<R[0].F.total_rank; ++p)
		psi[0].insert(p, (p==cog0) ? beta_rep : tau_module_oper.zero());

	for(int i=0; i<maxlev; ++i){
		psi[i+1].set_rank(R[i+1].F.total_rank);

		// Position -> cogenerator index map for F_{i+1}, for O(log n) lookup.
		std::map<int,int> pos_to_cog_next;
		for(unsigned a=0; a<R[i+1].F.generators.rank; ++a)
			pos_to_cog_next[(int)R[i+1].F.position_of_gens[a]] = (int)a;

		for(int x=0; x<R[i+1].F.total_rank; ++x){
			TVec u;
			if(invInj[i+1].solve(tau_module_oper.singleton(x), u)){
				TVec v; secQut[i].solve(u, v);
				TVec fv  = apply_chainmap(v, R[i].F, psi[i], R[i+bs].F, MOP);
				TVec dfv = R[i+bs+1].inj.maps_to( R[i+bs].qut.maps_to(fv) );
				psi[i+1].insert(x, proj_cogens(dfv, R[i+bs+1].F));
			} else {
				// Failing position: only cogenerator unit positions can be recovered.
				// Find j_a: the X_{i+1} index j such that inj_{i+1}(e_j) contains x.
				// Chain condition: sum_{y in inj(j_a)} coef_y * psi[i+1](y) = dfv.
				// Solve for psi[i+1](x) = (dfv + sum_{y!=x} coef_y*psi[i+1](y)) / coef_x.
				auto it = pos_to_cog_next.find(x);
				if(it != pos_to_cog_next.end()){
					int j_a = -1;
					TVec injja;
					for(unsigned j=0; j<(unsigned)R[i+1].Xrank && j_a<0; ++j){
						injja = R[i+1].inj.find(j);
						for(auto &tm : injja.dataArray)
							if((int)tm.ind == x){ j_a = (int)j; break; }
					}
					if(j_a >= 0){
						TVec v; secQut[i].solve(tau_module_oper.singleton(j_a), v);
						TVec fv  = apply_chainmap(v, R[i].F, psi[i], R[i+bs].F, MOP);
						TVec dfv = R[i+bs+1].inj.maps_to( R[i+bs].qut.maps_to(fv) );
						TVec rhs = proj_cogens(dfv, R[i+bs+1].F);
						tauPoly coef_x = 0;
						for(auto &tm : injja.dataArray){
							if((int)tm.ind == x){ coef_x = tm.coeficient; continue; }
							auto psi_y = psi[i+1].find(tm.ind);
							rhs = tau_module_oper.add(rhs, tau_module_oper.scalor_mult(tm.coeficient, psi_y));
						}
						// coef_x = tau^k; if k>0 divide rhs by tau^k (multiply by tau^{-k})
						if(coef_x > 0)
							rhs = tau_module_oper.scalor_mult(-coef_x, rhs);
						psi[i+1].insert(x, rhs);
					} else
						psi[i+1].insert(x, tau_module_oper.zero());
				} else
					psi[i+1].insert(x, tau_module_oper.zero());
			}
		}
	}
	// Chain map on cogenerators: M[i].find(a) = proj_cogens(f_i(cogenerator a)).
	// Cycle correction: if cog_a is a strict cycle (d_i(cog_a)=0), then
	// f_i(cog_a) must also be a cycle.  If it isn't, add cogenerators of F_{i+bs}
	// whose boundary equals the unwanted boundary, cancelling it mod 2.
	std::vector<matrix_mem<tauPoly>> M(maxlev + 1);
	for(int i=0; i<=maxlev; ++i){
		M[i].set_rank(R[i].F.generators.rank);
		for(unsigned a=0; a<R[i].F.generators.rank; ++a){
			int pa = (int)R[i].F.position_of_gens[a];
			TVec cog = tau_module_oper.singleton(pa);
			TVec fa  = apply_chainmap(cog, R[i].F, psi[i], R[i+bs].F, MOP);
			// Check whether cog_a is a strict cycle.
			if(i+1 < (int)R.size() && i+bs+1 < (int)R.size()){
				TVec d_cog = R[i+1].inj.maps_to(R[i].qut.maps_to(cog));
				if(tau_module_oper.isZero(d_cog)){
					// cog_a is a strict cycle; verify f_i(cog_a) lies in ker(d_{i+bs}).
					TVec dfa = R[i+bs+1].inj.maps_to(R[i+bs].qut.maps_to(fa));
					if(!tau_module_oper.isZero(dfa)){
						// Collect positions already in fa to avoid cancelling existing terms.
						std::set<int> fa_pks;
						for(auto &tm : fa.dataArray) fa_pks.insert((int)tm.ind);
						for(unsigned k=0; k<R[i+bs].F.generators.rank && !tau_module_oper.isZero(dfa); ++k){
							int pk = (int)R[i+bs].F.position_of_gens[k];
							if(fa_pks.count(pk)) continue;  // skip: already in fa
							TVec dk = R[i+bs+1].inj.maps_to(R[i+bs].qut.maps_to(tau_module_oper.singleton(pk)));
							// dk == dfa  ↔  dk + dfa == 0 over F_2
							if(tau_module_oper.isZero(tau_module_oper.add(dk, dfa))){
								fa  = tau_module_oper.add(fa,  tau_module_oper.singleton(pk));
								dfa = tau_module_oper.add(dfa, dk);  // = 0 now
							}
						}
					}
				}
			}
			M[i].insert(a, proj_cogens(fa, R[i+bs].F));
		}
	}

	//product of class (s,a) with beta: name M[s](cycle_rep(s,a)) in filtration s+bs.
	//Returns an empty vector and sets ok=false when the image falls outside the
	//cycle table (degree-boundary overflow: the product lives beyond max_deg).
	auto product = [&](int s, int a, bool &ok) -> TVec {
		TVec rep = (cyc[s].count(a)) ? cyc[s].at(a) : tau_module_oper.singleton(a);
		TVec img = M[s].maps_to(rep);
		try {
			ok = true;
			return find_cycle(cyc[s+bs], img);
		} catch (const std::out_of_range &) {
			ok = false;
			return TVec{};
		}
	};

	if(argc == 7){
		int as = std::atoi(argv[3]), aa = std::atoi(argv[4]);
		if(as < 0 || as > maxlev){ std::cerr << "alpha filtration out of lifted range (<= " << maxlev << ")\n"; return 1; }
		bool ok;
		TVec res = product(as, aa, ok);
		if(!ok)
			std::cout << "{" << as << "-" << aa << "} * {" << bs << "-" << bb << "}  =  ?\n";
		else
			std::cout << "{" << as << "-" << aa << "} * {" << bs << "-" << bb << "}  =  "
			          << fmt(as+bs, res) << "\n";
	} else {
		for(int i=0; i<=maxlev; ++i)
			for(unsigned a=0; a<R[i].F.generators.rank; ++a){
				bool ok;
				TVec res = product(i, a, ok);
				std::cout << "{" << i << "-" << a << "}\t->\t"
				          << (ok ? fmt(i+bs, res) : "?") << "\n";
			}
	}
	return 0;
}
