/*
 * spr_models.c - SPR binding model simulations
 *
 * Models:
 *   1:1 Langmuir:      dR/dt = ka*C*(Rmax-R) - kd*R
 *   Heterogeneous:     Two independent 1:1 sites
 *   Two-State:         A+B <=> AB <=> AB*  (conformational change)
 */

#include "spr_models.h"
#include <math.h>

/* ================================================================ */
/*  1:1 Langmuir - RK4                                               */
/* ================================================================ */

static void langmuir_rk4(double ka, double kd, double Rmax,
                         double C_M, double dt, double *R) {
    double k1 = dt * (ka * C_M * (Rmax - *R) - kd * (*R));
    double k2 = dt * (ka * C_M * (Rmax - (*R + k1/2)) - kd * (*R + k1/2));
    double k3 = dt * (ka * C_M * (Rmax - (*R + k2/2)) - kd * (*R + k2/2));
    double k4 = dt * (ka * C_M * (Rmax - (*R + k3))   - kd * (*R + k3));
    *R += (k1 + 2*k2 + 2*k3 + k4) / 6.0;
    if (*R < 0) *R = 0;
}

/* ================================================================ */
/*  1:1 Langmuir with Mass Transport - RK4 (Quasi-Steady-State)      */
/*  Cs equilibrates fast: 0 = tc*(C-Cs) - ka*Cs*(Rmax-R) + kd*R     */
/*    => Cs = (tc*C + kd*R) / (tc + ka*(Rmax-R))                    */
/*  dR/dt = ka*Cs*(Rmax-R) - kd*R                                   */
/* ================================================================ */

static inline double qss_cs(double tc, double C, double ka, double kd,
                            double Rmax, double R) {
    double free = Rmax - R;
    if (free < 0) free = 0;
    double denom = tc + ka * free;
    return (denom > 1e-30) ? (tc * C + kd * R) / denom : C;
}

static void langmuir_tc_rk4(double ka, double kd, double Rmax,
                            double tc, double C_M, double dt, double *R) {
    double cs, free;

    cs = qss_cs(tc, C_M, ka, kd, Rmax, *R);
    free = Rmax - *R; if (free < 0) free = 0;
    double k1 = dt * (ka * cs * free - kd * (*R));

    double r2 = *R + k1/2;
    cs = qss_cs(tc, C_M, ka, kd, Rmax, r2);
    free = Rmax - r2; if (free < 0) free = 0;
    double k2 = dt * (ka * cs * free - kd * r2);

    double r3 = *R + k2/2;
    cs = qss_cs(tc, C_M, ka, kd, Rmax, r3);
    free = Rmax - r3; if (free < 0) free = 0;
    double k3 = dt * (ka * cs * free - kd * r3);

    double r4 = *R + k3;
    cs = qss_cs(tc, C_M, ka, kd, Rmax, r4);
    free = Rmax - r4; if (free < 0) free = 0;
    double k4 = dt * (ka * cs * free - kd * r4);

    *R += (k1 + 2*k2 + 2*k3 + k4) / 6.0;
    if (*R < 0) *R = 0;
}

/* ================================================================ */
/*  Heterogeneous Ligand - Two independent 1:1 sites                 */
/* ================================================================ */

static void heterogeneous_rk4(double ka1, double kd1, double Rmax1,
                              double ka2, double kd2, double Rmax2,
                              double C_M, double dt, double *R1, double *R2) {
    /* Site 1 */
    double k1a = dt * (ka1 * C_M * (Rmax1 - *R1) - kd1 * (*R1));
    double k2a = dt * (ka1 * C_M * (Rmax1 - (*R1 + k1a/2)) - kd1 * (*R1 + k1a/2));
    double k3a = dt * (ka1 * C_M * (Rmax1 - (*R1 + k2a/2)) - kd1 * (*R1 + k2a/2));
    double k4a = dt * (ka1 * C_M * (Rmax1 - (*R1 + k3a))   - kd1 * (*R1 + k3a));
    *R1 += (k1a + 2*k2a + 2*k3a + k4a) / 6.0;
    if (*R1 < 0) *R1 = 0;

    /* Site 2 */
    double k1b = dt * (ka2 * C_M * (Rmax2 - *R2) - kd2 * (*R2));
    double k2b = dt * (ka2 * C_M * (Rmax2 - (*R2 + k1b/2)) - kd2 * (*R2 + k1b/2));
    double k3b = dt * (ka2 * C_M * (Rmax2 - (*R2 + k2b/2)) - kd2 * (*R2 + k2b/2));
    double k4b = dt * (ka2 * C_M * (Rmax2 - (*R2 + k3b))   - kd2 * (*R2 + k3b));
    *R2 += (k1b + 2*k2b + 2*k3b + k4b) / 6.0;
    if (*R2 < 0) *R2 = 0;
}

/* ================================================================ */
/*  Heterogeneous Ligand with Mass Transport - RK4 (QSS)             */
/*  Cs = (tc*C + kd1*R1 + kd2*R2) / (tc + ka1*(Rmax1-R1) + ka2*(Rmax2-R2)) */
/*  dR1/dt = ka1*Cs*(Rmax1-R1) - kd1*R1                             */
/*  dR2/dt = ka2*Cs*(Rmax2-R2) - kd2*R2                             */
/* ================================================================ */

static inline double qss_cs_het(double tc, double C,
                                double ka1, double kd1, double Rmax1, double R1,
                                double ka2, double kd2, double Rmax2, double R2) {
    double f1 = Rmax1 - R1; if (f1 < 0) f1 = 0;
    double f2 = Rmax2 - R2; if (f2 < 0) f2 = 0;
    double denom = tc + ka1*f1 + ka2*f2;
    return (denom > 1e-30) ? (tc*C + kd1*R1 + kd2*R2) / denom : C;
}

static void heterogeneous_tc_rk4(double ka1, double kd1, double Rmax1,
                                 double ka2, double kd2, double Rmax2,
                                 double tc, double C_M, double dt,
                                 double *R1, double *R2) {
    double cs, f1, f2;

    /* k1 */
    cs = qss_cs_het(tc, C_M, ka1, kd1, Rmax1, *R1, ka2, kd2, Rmax2, *R2);
    f1 = Rmax1 - *R1; if (f1 < 0) f1 = 0;
    f2 = Rmax2 - *R2; if (f2 < 0) f2 = 0;
    double k1a = dt*(ka1*cs*f1 - kd1*(*R1));
    double k1b = dt*(ka2*cs*f2 - kd2*(*R2));

    /* k2 */
    double ra2 = *R1+k1a/2, rb2 = *R2+k1b/2;
    cs = qss_cs_het(tc, C_M, ka1, kd1, Rmax1, ra2, ka2, kd2, Rmax2, rb2);
    f1 = Rmax1-ra2; if (f1 < 0) f1 = 0;
    f2 = Rmax2-rb2; if (f2 < 0) f2 = 0;
    double k2a = dt*(ka1*cs*f1 - kd1*ra2);
    double k2b = dt*(ka2*cs*f2 - kd2*rb2);

    /* k3 */
    double ra3 = *R1+k2a/2, rb3 = *R2+k2b/2;
    cs = qss_cs_het(tc, C_M, ka1, kd1, Rmax1, ra3, ka2, kd2, Rmax2, rb3);
    f1 = Rmax1-ra3; if (f1 < 0) f1 = 0;
    f2 = Rmax2-rb3; if (f2 < 0) f2 = 0;
    double k3a = dt*(ka1*cs*f1 - kd1*ra3);
    double k3b = dt*(ka2*cs*f2 - kd2*rb3);

    /* k4 */
    double ra4 = *R1+k3a, rb4 = *R2+k3b;
    cs = qss_cs_het(tc, C_M, ka1, kd1, Rmax1, ra4, ka2, kd2, Rmax2, rb4);
    f1 = Rmax1-ra4; if (f1 < 0) f1 = 0;
    f2 = Rmax2-rb4; if (f2 < 0) f2 = 0;
    double k4a = dt*(ka1*cs*f1 - kd1*ra4);
    double k4b = dt*(ka2*cs*f2 - kd2*rb4);

    *R1 += (k1a + 2*k2a + 2*k3a + k4a) / 6.0;
    *R2 += (k1b + 2*k2b + 2*k3b + k4b) / 6.0;
    if (*R1 < 0) *R1 = 0;
    if (*R2 < 0) *R2 = 0;
}

/* ================================================================ */
/*  Two-State Conformational Change - RK4                            */
/*  A + B <=> AB <=> AB*                                             */
/*  dAB/dt  = ka1*C*(Rmax-AB-ABs) - kd1*AB - ka2*AB + kd2*ABs       */
/*  dABs/dt = ka2*AB - kd2*ABs                                      */
/*  Signal  = AB + ABs                                               */
/* ================================================================ */

static void twostate_rk4(double ka1, double kd1, double ka2, double kd2,
                         double Rmax, double C_M, double dt,
                         double *AB, double *ABs) {
    double free;

    /* k1 */
    free = Rmax - *AB - *ABs;
    double k1a = dt * (ka1*C_M*free - kd1*(*AB) - ka2*(*AB) + kd2*(*ABs));
    double k1b = dt * (ka2*(*AB) - kd2*(*ABs));

    /* k2 */
    double ab2 = *AB + k1a/2, abs2 = *ABs + k1b/2;
    free = Rmax - ab2 - abs2;
    double k2a = dt * (ka1*C_M*free - kd1*ab2 - ka2*ab2 + kd2*abs2);
    double k2b = dt * (ka2*ab2 - kd2*abs2);

    /* k3 */
    double ab3 = *AB + k2a/2, abs3 = *ABs + k2b/2;
    free = Rmax - ab3 - abs3;
    double k3a = dt * (ka1*C_M*free - kd1*ab3 - ka2*ab3 + kd2*abs3);
    double k3b = dt * (ka2*ab3 - kd2*abs3);

    /* k4 */
    double ab4 = *AB + k3a, abs4 = *ABs + k3b;
    free = Rmax - ab4 - abs4;
    double k4a = dt * (ka1*C_M*free - kd1*ab4 - ka2*ab4 + kd2*abs4);
    double k4b = dt * (ka2*ab4 - kd2*abs4);

    *AB  += (k1a + 2*k2a + 2*k3a + k4a) / 6.0;
    *ABs += (k1b + 2*k2b + 2*k3b + k4b) / 6.0;

    if (*AB < 0) *AB = 0;
    if (*ABs < 0) *ABs = 0;
    if (*AB + *ABs > Rmax) {
        double s = Rmax / (*AB + *ABs + 1e-30);
        *AB *= s;
        *ABs *= s;
    }
}

/* ================================================================ */
/*  Two-State with Mass Transport - RK4 (QSS)                       */
/*  Cs = (tc*C + kd1*AB) / (tc + ka1*(Rmax-AB-ABs))                 */
/*  dAB/dt  = ka1*Cs*(Rmax-AB-ABs) - kd1*AB - ka2*AB + kd2*ABs     */
/*  dABs/dt = ka2*AB - kd2*ABs                                      */
/* ================================================================ */

static inline double qss_cs_2s(double tc, double C, double ka1, double kd1,
                               double Rmax, double AB, double ABs) {
    double free = Rmax - AB - ABs;
    if (free < 0) free = 0;
    double denom = tc + ka1 * free;
    return (denom > 1e-30) ? (tc * C + kd1 * AB) / denom : C;
}

static void twostate_tc_rk4(double ka1, double kd1, double ka2, double kd2,
                            double Rmax, double tc, double C_M, double dt,
                            double *AB, double *ABs) {
    double cs, free;

    /* k1 */
    cs = qss_cs_2s(tc, C_M, ka1, kd1, Rmax, *AB, *ABs);
    free = Rmax - *AB - *ABs; if (free < 0) free = 0;
    double k1a = dt*(ka1*cs*free - kd1*(*AB) - ka2*(*AB) + kd2*(*ABs));
    double k1b = dt*(ka2*(*AB) - kd2*(*ABs));

    /* k2 */
    double ab2 = *AB+k1a/2, abs2 = *ABs+k1b/2;
    cs = qss_cs_2s(tc, C_M, ka1, kd1, Rmax, ab2, abs2);
    free = Rmax - ab2 - abs2; if (free < 0) free = 0;
    double k2a = dt*(ka1*cs*free - kd1*ab2 - ka2*ab2 + kd2*abs2);
    double k2b = dt*(ka2*ab2 - kd2*abs2);

    /* k3 */
    double ab3 = *AB+k2a/2, abs3 = *ABs+k2b/2;
    cs = qss_cs_2s(tc, C_M, ka1, kd1, Rmax, ab3, abs3);
    free = Rmax - ab3 - abs3; if (free < 0) free = 0;
    double k3a = dt*(ka1*cs*free - kd1*ab3 - ka2*ab3 + kd2*abs3);
    double k3b = dt*(ka2*ab3 - kd2*abs3);

    /* k4 */
    double ab4 = *AB+k3a, abs4 = *ABs+k3b;
    cs = qss_cs_2s(tc, C_M, ka1, kd1, Rmax, ab4, abs4);
    free = Rmax - ab4 - abs4; if (free < 0) free = 0;
    double k4a = dt*(ka1*cs*free - kd1*ab4 - ka2*ab4 + kd2*abs4);
    double k4b = dt*(ka2*ab4 - kd2*abs4);

    *AB  += (k1a + 2*k2a + 2*k3a + k4a) / 6.0;
    *ABs += (k1b + 2*k2b + 2*k3b + k4b) / 6.0;

    if (*AB < 0) *AB = 0;
    if (*ABs < 0) *ABs = 0;
    if (*AB + *ABs > Rmax) {
        double s = Rmax / (*AB + *ABs + 1e-30);
        *AB *= s;
        *ABs *= s;
    }
}

/* ================================================================ */
/*  Analytical 1:1 Langmuir Solution                                 */
/*  Association: R(t) = Req * (1 - exp(-kobs*t))                     */
/*    where kobs = ka*C + kd, Req = ka*C*Rmax / (ka*C + kd)         */
/*  Dissociation: R(t) = R0 * exp(-kd*(t - t0))                     */
/* ================================================================ */

static void langmuir_analytical(double ka, double kd, double Rmax,
                                double C_M, double t_assoc_end,
                                const double *time, int fit_start,
                                int npoints, double *result) {
    double kobs = ka * C_M + kd;
    double Req = (kobs > 1e-30) ? ka * C_M * Rmax / kobs : 0.0;

    /* R at end of association */
    double R_at_assoc_end = Req * (1.0 - exp(-kobs * t_assoc_end));

    for (int i = fit_start; i < npoints; i++) {
        double t = time[i];
        if (t <= t_assoc_end) {
            result[i] = Req * (1.0 - exp(-kobs * t));
        } else {
            result[i] = R_at_assoc_end * exp(-kd * (t - t_assoc_end));
        }
    }
}

/* ================================================================ */
/*  MCK Cycle Simulation                                             */
/* ================================================================ */

void simulate_mck_cycle(ModelType model, const double *params,
                        const Cycle *cy, double *result,
                        double t_assoc_end, double dt_max, double tc) {
    double C_M = cy->conc_nM * 1e-9;

    /* Fill pre-injection region */
    for (int i = 0; i < cy->fit_start; i++) result[i] = 0.0;

    if (model == MODEL_LANGMUIR) {
        double ka = pow(10.0, params[0]);
        double kd = pow(10.0, params[1]);
        double Rmax = params[2];

        if (tc > 0) {
            /* QSS mass transport model */
            double R = 0.0, t = 0.0;
            for (int i = cy->fit_start; i < cy->npoints; i++) {
                double target = cy->time[i];
                while (t < target - 1e-8) {
                    double ddt = target - t;
                    if (ddt > dt_max) ddt = dt_max;
                    if (t < t_assoc_end && t + ddt > t_assoc_end)
                        ddt = t_assoc_end - t;
                    double C = (t < t_assoc_end) ? C_M : 0.0;
                    langmuir_tc_rk4(ka, kd, Rmax, tc, C, ddt, &R);
                    t += ddt;
                }
                result[i] = R;
            }
        } else {
            /* Use analytical solution — no RK4 needed */
            langmuir_analytical(ka, kd, Rmax, C_M, t_assoc_end,
                               cy->time, cy->fit_start, cy->npoints, result);
        }

    } else if (model == MODEL_HETEROGENEOUS) {
        double ka1 = pow(10.0, params[0]);
        double kd1 = pow(10.0, params[1]);
        double Rmax1 = params[2];
        double ka2 = pow(10.0, params[3]);
        double kd2 = pow(10.0, params[4]);
        double Rmax2 = params[5];

        if (tc > 0) {
            double R1 = 0.0, R2 = 0.0, t = 0.0;
            for (int i = cy->fit_start; i < cy->npoints; i++) {
                double target = cy->time[i];
                while (t < target - 1e-8) {
                    double ddt = target - t;
                    if (ddt > dt_max) ddt = dt_max;
                    if (t < t_assoc_end && t + ddt > t_assoc_end)
                        ddt = t_assoc_end - t;
                    double C = (t < t_assoc_end) ? C_M : 0.0;
                    heterogeneous_tc_rk4(ka1, kd1, Rmax1, ka2, kd2, Rmax2,
                                        tc, C, ddt, &R1, &R2);
                    t += ddt;
                }
                result[i] = R1 + R2;
            }
        } else {
            double R1 = 0.0, R2 = 0.0, t = 0.0;
            for (int i = cy->fit_start; i < cy->npoints; i++) {
                double target = cy->time[i];
                while (t < target - 1e-8) {
                    double ddt = target - t;
                    if (ddt > dt_max) ddt = dt_max;
                    if (t < t_assoc_end && t + ddt > t_assoc_end)
                        ddt = t_assoc_end - t;
                    double C = (t < t_assoc_end) ? C_M : 0.0;
                    heterogeneous_rk4(ka1, kd1, Rmax1, ka2, kd2, Rmax2, C, ddt, &R1, &R2);
                    t += ddt;
                }
                result[i] = R1 + R2;
            }
        }

    } else if (model == MODEL_TWOSTATE) {
        double ka1 = pow(10.0, params[0]);
        double kd1 = pow(10.0, params[1]);
        double ka2 = pow(10.0, params[2]);
        double kd2 = pow(10.0, params[3]);
        double Rmax = params[4];

        if (tc > 0) {
            double AB = 0.0, ABs = 0.0, t = 0.0;
            for (int i = cy->fit_start; i < cy->npoints; i++) {
                double target = cy->time[i];
                while (t < target - 1e-8) {
                    double ddt = target - t;
                    if (ddt > dt_max) ddt = dt_max;
                    if (t < t_assoc_end && t + ddt > t_assoc_end)
                        ddt = t_assoc_end - t;
                    double C = (t < t_assoc_end) ? C_M : 0.0;
                    twostate_tc_rk4(ka1, kd1, ka2, kd2, Rmax, tc, C, ddt,
                                   &AB, &ABs);
                    t += ddt;
                }
                result[i] = AB + ABs;
            }
        } else {
            double AB = 0.0, ABs = 0.0, t = 0.0;
            for (int i = cy->fit_start; i < cy->npoints; i++) {
                double target = cy->time[i];
                while (t < target - 1e-8) {
                    double ddt = target - t;
                    if (ddt > dt_max) ddt = dt_max;
                    if (t < t_assoc_end && t + ddt > t_assoc_end)
                        ddt = t_assoc_end - t;
                    double C = (t < t_assoc_end) ? C_M : 0.0;
                    twostate_rk4(ka1, kd1, ka2, kd2, Rmax, C, ddt, &AB, &ABs);
                    t += ddt;
                }
                result[i] = AB + ABs;
            }
        }
    }
}

/* ================================================================ */
/*  SCK Trace Simulation                                             */
/* ================================================================ */

/* Get concentration at time t from injection schedule */
static double sck_conc_at(const SPRData *data, double t) {
    for (int i = 0; i < data->n_injections; i++) {
        if (t >= data->inj_start[i] && t < data->inj_end[i])
            return data->inj_conc_nM[i] * 1e-9;
    }
    return 0.0;  /* Dissociation or between injections */
}

void simulate_sck_trace(ModelType model, const double *params,
                        const SPRData *data, double *result,
                        double dt_max, double tc) {
    const Cycle *cy = &data->cycles[0];

    for (int i = 0; i < cy->fit_start; i++) result[i] = 0.0;

    if (model == MODEL_LANGMUIR) {
        double ka = pow(10.0, params[0]);
        double kd = pow(10.0, params[1]);
        double Rmax = params[2];

        if (tc > 0) {
            double R = 0.0, t = 0.0;
            for (int i = cy->fit_start; i < cy->npoints; i++) {
                double target = cy->time[i];
                while (t < target - 1e-8) {
                    double dt = target - t;
                    if (dt > dt_max) dt = dt_max;
                    double C = sck_conc_at(data, t + dt/2);
                    langmuir_tc_rk4(ka, kd, Rmax, tc, C, dt, &R);
                    t += dt;
                }
                result[i] = R;
            }
        } else {
            double R = 0.0, t = 0.0;
            for (int i = cy->fit_start; i < cy->npoints; i++) {
                double target = cy->time[i];
                while (t < target - 1e-8) {
                    double dt = target - t;
                    if (dt > dt_max) dt = dt_max;
                    double C = sck_conc_at(data, t + dt/2);
                    langmuir_rk4(ka, kd, Rmax, C, dt, &R);
                    t += dt;
                }
                result[i] = R;
            }
        }

    } else if (model == MODEL_HETEROGENEOUS) {
        double ka1 = pow(10.0, params[0]);
        double kd1 = pow(10.0, params[1]);
        double Rmax1 = params[2];
        double ka2 = pow(10.0, params[3]);
        double kd2 = pow(10.0, params[4]);
        double Rmax2 = params[5];

        if (tc > 0) {
            double R1 = 0.0, R2 = 0.0, t = 0.0;
            for (int i = cy->fit_start; i < cy->npoints; i++) {
                double target = cy->time[i];
                while (t < target - 1e-8) {
                    double dt = target - t;
                    if (dt > dt_max) dt = dt_max;
                    double C = sck_conc_at(data, t + dt/2);
                    heterogeneous_tc_rk4(ka1, kd1, Rmax1, ka2, kd2, Rmax2,
                                        tc, C, dt, &R1, &R2);
                    t += dt;
                }
                result[i] = R1 + R2;
            }
        } else {
            double R1 = 0.0, R2 = 0.0, t = 0.0;
            for (int i = cy->fit_start; i < cy->npoints; i++) {
                double target = cy->time[i];
                while (t < target - 1e-8) {
                    double dt = target - t;
                    if (dt > dt_max) dt = dt_max;
                    double C = sck_conc_at(data, t + dt/2);
                    heterogeneous_rk4(ka1, kd1, Rmax1, ka2, kd2, Rmax2, C, dt, &R1, &R2);
                    t += dt;
                }
                result[i] = R1 + R2;
            }
        }

    } else if (model == MODEL_TWOSTATE) {
        double ka1 = pow(10.0, params[0]);
        double kd1 = pow(10.0, params[1]);
        double ka2 = pow(10.0, params[2]);
        double kd2 = pow(10.0, params[3]);
        double Rmax = params[4];

        if (tc > 0) {
            double AB = 0.0, ABs = 0.0, t = 0.0;
            for (int i = cy->fit_start; i < cy->npoints; i++) {
                double target = cy->time[i];
                while (t < target - 1e-8) {
                    double dt = target - t;
                    if (dt > dt_max) dt = dt_max;
                    double C = sck_conc_at(data, t + dt/2);
                    twostate_tc_rk4(ka1, kd1, ka2, kd2, Rmax, tc, C, dt,
                                   &AB, &ABs);
                    t += dt;
                }
                result[i] = AB + ABs;
            }
        } else {
            double AB = 0.0, ABs = 0.0, t = 0.0;
            for (int i = cy->fit_start; i < cy->npoints; i++) {
                double target = cy->time[i];
                while (t < target - 1e-8) {
                    double dt = target - t;
                    if (dt > dt_max) dt = dt_max;
                    double C = sck_conc_at(data, t + dt/2);
                    twostate_rk4(ka1, kd1, ka2, kd2, Rmax, C, dt, &AB, &ABs);
                    t += dt;
                }
                result[i] = AB + ABs;
            }
        }
    }
}
