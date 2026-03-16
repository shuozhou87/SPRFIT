/*
 * spr_optim.h - Optimization and cost functions for SPR fitting
 */

#ifndef SPR_OPTIM_H
#define SPR_OPTIM_H

#include "spr_types.h"

/* Set up the global fitting context */
void optim_setup(SPRData *data, const FitConfig *cfg);

/* Clean up */
void optim_cleanup(void);

/*
 * Run the full fitting pipeline:
 *   1. Multi-restart Nelder-Mead global optimization
 *   2. Refinement pass with all data points
 *   3. Compute statistics (R², RMS, chi²)
 *   4. Steady-state affinity analysis
 */
int optim_fit(SPRData *data, const FitConfig *cfg, FitResult *result);

#endif /* SPR_OPTIM_H */
