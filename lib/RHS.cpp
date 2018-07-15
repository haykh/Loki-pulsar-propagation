#include <math.h>
#include <vector>

#include "aux/read_write.h"
#include "RHS.h"
#include "constants.h"
#include "process_functions.h"
#include "initialize.h"

void RHS(double R, double *f, double *dydx) {
	double coeff = constants::R_star * Globals::omega / (2.0 * constants::c);

	double LL = Lambda (R);
	double QQ = Q (R);
	double BB = BetaB (R);
	double DD = delta (R);

	dydx[0] = coeff * (-LL / QQ - LL * cos(2 * f[0] - 2 * BB - 2 * DD) * sinh(2 * f[1]));
	dydx[1] = coeff * LL * sin(2 * f[0] - 2 * BB - 2 * DD) * cosh(2 * f[1]);
}

void RHS_float (float R, float *f, float *dydx) {
	float coeff = (float)(constants::R_star * Globals::omega / (2.0 * constants::c));

	float LL = (float)(Lambda((float)(R)));
	float QQ = (float)(Q((float)(R)));
	float BB = (float)(BetaB((float)(R)));
	float DD = (float)(delta((float)(R)));

	dydx[0] = coeff * (-LL / QQ - LL * cos(2 * f[0] - 2 * BB - 2 * DD) * sinh(2 * f[1]));
	dydx[1] = coeff * LL * sin(2 * f[0] - 2 * BB - 2 * DD) * cosh(2 * f[1]);
}
