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

/* Build parameter layout from model, advanced config, and data */
ParamLayout build_param_layout(ModelType model, const AdvancedConfig *adv,
                               const SPRData *data);

/* Apply RI + drift corrections to a simulation buffer */
void apply_local_corrections(const double *params, const ParamLayout *lay,
                             int cycle_idx, const Cycle *cy, double *buf,
                             double t_assoc_end);

/*
 * Run the full fitting pipeline:
 *   1. Multi-restart Nelder-Mead global optimization
 *   2. Staged optimization for local params (RI, drift)
 *   3. Compute statistics (R², RMS, chi²)
 *   4. Steady-state affinity analysis
 */
int optim_fit(SPRData *data, const FitConfig *cfg, FitResult *result);

#endif /* SPR_OPTIM_H */
