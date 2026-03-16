/*
 * spr_models.h - SPR binding model simulations
 */

#ifndef SPR_MODELS_H
#define SPR_MODELS_H

#include "spr_types.h"

/*
 * Simulate a single MCK cycle with the specified model.
 * params layout depends on model:
 *   LANGMUIR:      [log10(ka), log10(kd), Rmax]
 *   HETEROGENEOUS: [log10(ka1), log10(kd1), Rmax1, log10(ka2), log10(kd2), Rmax2]
 *   TWOSTATE:      [log10(ka1), log10(kd1), log10(ka2), log10(kd2), Rmax]
 */
void simulate_mck_cycle(ModelType model, const double *params,
                        const Cycle *cy, double *result,
                        double t_assoc_end, double dt, double tc);

/*
 * Simulate an SCK trace (single continuous trace with staircase injections).
 * Same params layout as MCK. tc=0 disables mass transport.
 */
void simulate_sck_trace(ModelType model, const double *params,
                        const SPRData *data, double *result,
                        double dt, double tc);

#endif /* SPR_MODELS_H */
