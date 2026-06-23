// Standalone capability test for the WIDE_EXPONENT toggle: confirms
// MotSteenrodOp::hi(n) raises the n<=7 ceiling under EXPONENT_WIDTH=64,
// and that the default 32-bit build still hits the assert at n=8 unchanged.
#include "mot_steenrod.h"
#include <iostream>

int main(){
	MotSteenrodOp op(nullptr, 0);
	std::cout << "EXPONENT_WIDTH=" << EXPONENT_WIDTH
	          << " xnMaxExpo[1]=" << xnMaxExpo[1] << "\n";
	for(int n = 7; n <= 10; ++n){
		motSteenrod cls = op.hi(n);
		exponent e = cls.dataArray[0].ind;
		int got = xnVal(e, 1);
		bool ok = (got == (1 << n));
		std::cout << "n=" << n << " -> xi_1^" << got
		          << " (expected " << (1 << n) << ") "
		          << (ok ? "OK" : "MISMATCH") << "\n";
	}
	return 0;
}
