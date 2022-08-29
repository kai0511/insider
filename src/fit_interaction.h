#ifndef __FIT_INTERACTION__
#define __FIT_INTERACTION__

#include "../inst/include/insider_types.h"

void fit_interaction(const mat& residual, const mat& train_indicator, const mat& interactions, const uvec& interaction_indicator, const mat& column_factor, 
                     const double lambda, const int tuning, const double& tol /* = 1e-10 */, const int n_cores /* = 10 */);
#endif