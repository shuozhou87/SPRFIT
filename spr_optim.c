/*
 * spr_optim.c - Multi-threaded Nelder-Mead optimizer for SPR fitting
 *
 * Optimizations:
 *   - pthreads: parallel random restarts across CPU cores
 *   - Thread-local buffers for simulation results
 *   - Early stagnation detection in Nelder-Mead
 */

#include "spr_optim.h"
#include "spr_models.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <pthread.h>

/* ================================================================ */
/*  Thread-Safe Fitting Context                                      */
/* ================================================================ */

/* Shared (read-only during fitting) */
static SPRData    *g_data = NULL;
static FitConfig   g_cfg;
static ParamLayout g_layout;

/* Thread-local simulation buffer */
static __thread double *tl_buf = NULL;
static __thread int     tl_skip = 1;

void optim_setup(SPRData *data, const FitConfig *cfg) {
    g_data = data;
    g_cfg  = *cfg;
    /* Main thread buffer */
    tl_buf  = (double*)malloc(sizeof(double) * MAX_POINTS);
    tl_skip = cfg->fit_skip;
}

void optim_cleanup(void) {
    free(tl_buf);
    tl_buf = NULL;
    g_data = NULL;
}

/* ================================================================ */
/*  Parameter Layout Builder                                         */
/* ================================================================ */

ParamLayout build_param_layout(ModelType model, const AdvancedConfig *adv,
                               const SPRData *data) {
    ParamLayout lay;
    memset(&lay, 0, sizeof(lay));

    lay.n_global = MODEL_NPARAMS[model];

    /* Build active cycle mapping */
    lay.n_active = 0;
    for (int c = 0; c < data->ncycles; c++) {
        if (data->excluded[c]) {
            lay.cycle_to_active[c] = -1;
        } else {
            lay.active_map[lay.n_active] = c;
            lay.cycle_to_active[c] = lay.n_active;
            lay.n_active++;
        }
    }

    int idx = lay.n_global;

    /* RI: one param per active cycle (MCK), or one param (SCK) */
    if (adv->enable_ri) {
        lay.off_ri = idx;
        idx += lay.n_active;
    } else {
        lay.off_ri = -1;
    }

    /* Drift: single global param */
    if (adv->enable_drift) {
        lay.off_drift = idx;
        idx += 1;
    } else {
        lay.off_drift = -1;
    }

    /* Mass transport: single global param (log10 space) */
    if (adv->enable_tc) {
        lay.off_tc = idx;
        idx += 1;
    } else {
        lay.off_tc = -1;
    }

    lay.total = idx;
    return lay;
}

/* ================================================================ */
/*  Cost Function (Global SSR)                                       */
/* ================================================================ */

static double cost(const double *params) {
    ModelType model = g_cfg.model;

    /* Validate Rmax parameters are positive */
    if (model == MODEL_LANGMUIR && params[2] <= 0) return 1e30;
    if (model == MODEL_HETEROGENEOUS && (params[2] <= 0 || params[5] <= 0)) return 1e30;
    if (model == MODEL_TWOSTATE && params[4] <= 0) return 1e30;

    double ssr = 0.0;
    int skip = tl_skip;
    double t_assoc = g_cfg.t_assoc_end;
    double tc = (g_layout.off_tc >= 0) ? pow(10.0, params[g_layout.off_tc]) : 0.0;

    if (g_data->mode == MODE_MCK) {
        for (int c = 0; c < g_data->ncycles; c++) {
            if (g_data->excluded[c]) continue;
            const Cycle *cy = &g_data->cycles[c];
            simulate_mck_cycle(model, params, cy, tl_buf, t_assoc, g_cfg.dt, tc);

            /* Get local/global correction params for this cycle */
            int ai = g_layout.cycle_to_active[c];
            double ri = (g_layout.off_ri >= 0) ? params[g_layout.off_ri + ai] : 0.0;
            double drift = (g_layout.off_drift >= 0) ? params[g_layout.off_drift] : 0.0;

            for (int i = cy->fit_start; i < cy->npoints; i += skip) {
                if (cy->skip[i]) continue;
                double sim = tl_buf[i];
                /* Apply RI: step function during association */
                if (ri != 0.0 && cy->time[i] >= 0.0 && cy->time[i] < t_assoc)
                    sim += ri;
                /* Apply drift: linear over all time */
                if (drift != 0.0 && cy->time[i] >= 0.0)
                    sim += drift * cy->time[i];
                double r = cy->resp[i] - sim;
                ssr += r * r;
            }
        }
    } else {
        /* SCK: single trace */
        simulate_sck_trace(model, params, g_data, tl_buf, g_cfg.dt, tc);
        const Cycle *cy = &g_data->cycles[0];

        /* Get RI and drift corrections (cycle 0) */
        double ri = 0.0, drift = 0.0;
        if (g_layout.off_ri >= 0) {
            int ai = g_layout.cycle_to_active[0];
            if (ai >= 0) ri = params[g_layout.off_ri + ai];
        }
        if (g_layout.off_drift >= 0)
            drift = params[g_layout.off_drift];

        for (int i = cy->fit_start; i < cy->npoints; i += skip) {
            if (cy->skip[i]) continue;
            double sim = tl_buf[i];
            if (ri != 0.0 && cy->time[i] >= 0.0 && cy->time[i] < t_assoc)
                sim += ri;
            if (drift != 0.0 && cy->time[i] >= 0.0)
                sim += drift * cy->time[i];
            double r = cy->resp[i] - sim;
            ssr += r * r;
        }
    }

    return ssr;
}

/* ================================================================ */
/*  Nelder-Mead Simplex Optimizer                                    */
/* ================================================================ */

typedef struct { double p[MAX_PARAMS]; double v; } Vertex;

static int cmp_vertex(const void *a, const void *b) {
    double va = ((const Vertex*)a)->v;
    double vb = ((const Vertex*)b)->v;
    return (va > vb) - (va < vb);
}

/* Thread-safe random using seed pointer */
static double rand_r_double(unsigned int *seed) {
    return (double)rand_r(seed) / RAND_MAX;
}

static void nelder_mead_ex(int nparams, double *best, double *best_val,
                        const double lo[], const double hi[],
                        int maxiter, double tol, unsigned int *seed,
                        const double *x0) {
    const int N = nparams;

    /* For high-dimensional problems, allocate on heap */
    Vertex *s = (Vertex*)malloc(sizeof(Vertex) * (N + 1));

    /* Initialize simplex: vertex 0 = x0 if provided, rest random */
    int start = 0;
    if (x0) {
        for (int j = 0; j < N; j++)
            s[0].p[j] = x0[j];
        s[0].v = cost(s[0].p);
        start = 1;
    }
    for (int i = start; i <= N; i++) {
        for (int j = 0; j < N; j++)
            s[i].p[j] = lo[j] + (hi[j] - lo[j]) * rand_r_double(seed);
        s[i].v = cost(s[i].p);
    }

    double prev_best = s[0].v;
    int stagnation = 0;
    const int STAG_LIMIT = 100;  /* Stop if no improvement for 100 iters */

    double *c  = (double*)malloc(sizeof(double) * N);
    double *xr = (double*)malloc(sizeof(double) * N);
    double *xe = (double*)malloc(sizeof(double) * N);
    double *xc = (double*)malloc(sizeof(double) * N);

    for (int iter = 0; iter < maxiter; iter++) {
        qsort(s, N + 1, sizeof(Vertex), cmp_vertex);

        double range = s[N].v - s[0].v;
        if (range < tol * (fabs(s[0].v) + 1e-30)) break;

        /* Early stagnation detection */
        if (s[0].v < prev_best - 1e-10 * (fabs(prev_best) + 1e-30)) {
            prev_best = s[0].v;
            stagnation = 0;
        } else {
            stagnation++;
            if (stagnation >= STAG_LIMIT) break;
        }

        /* Centroid (exclude worst) */
        memset(c, 0, sizeof(double) * N);
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
                c[j] += s[i].p[j] / N;

        /* Reflection */
        for (int j = 0; j < N; j++) xr[j] = c[j] + (c[j] - s[N].p[j]);
        double vr = cost(xr);

        if (vr < s[0].v) {
            for (int j = 0; j < N; j++) xe[j] = c[j] + 2.0 * (xr[j] - c[j]);
            double ve = cost(xe);
            if (ve < vr) { memcpy(s[N].p, xe, sizeof(double)*N); s[N].v = ve; }
            else         { memcpy(s[N].p, xr, sizeof(double)*N); s[N].v = vr; }
        } else if (vr < s[N-1].v) {
            memcpy(s[N].p, xr, sizeof(double)*N); s[N].v = vr;
        } else {
            if (vr < s[N].v) {
                for (int j = 0; j < N; j++) xc[j] = c[j] + 0.5*(xr[j] - c[j]);
            } else {
                for (int j = 0; j < N; j++) xc[j] = c[j] + 0.5*(s[N].p[j] - c[j]);
            }
            double vc = cost(xc);
            if (vc < s[N].v) {
                memcpy(s[N].p, xc, sizeof(double)*N); s[N].v = vc;
            } else {
                for (int i = 1; i <= N; i++) {
                    for (int j = 0; j < N; j++)
                        s[i].p[j] = s[0].p[j] + 0.5*(s[i].p[j] - s[0].p[j]);
                    s[i].v = cost(s[i].p);
                }
            }
        }
    }

    qsort(s, N + 1, sizeof(Vertex), cmp_vertex);
    memcpy(best, s[0].p, sizeof(double) * N);
    *best_val = s[0].v;

    free(s);
    free(c);
    free(xr);
    free(xe);
    free(xc);
}

static void nelder_mead(int nparams, double *best, double *best_val,
                        const double lo[], const double hi[],
                        int maxiter, double tol, unsigned int *seed) {
    nelder_mead_ex(nparams, best, best_val, lo, hi, maxiter, tol, seed, NULL);
}

/* ================================================================ */
/*  Steady-State Affinity (MCK only)                                 */
/* ================================================================ */

static void compute_req(const SPRData *data, double *req, double t_lo, double t_hi) {
    for (int c = 0; c < data->ncycles; c++) {
        const Cycle *cy = &data->cycles[c];
        double sum = 0;
        int count = 0;
        for (int i = cy->fit_start; i < cy->npoints; i++) {
            if (cy->skip[i]) continue;
            if (cy->time[i] >= t_lo && cy->time[i] <= t_hi) {
                sum += cy->resp[i];
                count++;
            }
        }
        req[c] = (count > 0) ? sum / count : 0.0;
    }
}

/* Simple 2-parameter Nelder-Mead for Req = Rmax * C / (KD + C)
 * kinetic_KD_nM: hint from kinetic fit (kd/ka * 1e9), used to seed restart 0 */
static void ss_fit(const SPRData *data, const double *req,
                   double kinetic_KD_nM, FitResult *result) {
    double max_req = 0;
    int n_active = 0;
    for (int c = 0; c < data->ncycles; c++) {
        if (data->excluded[c]) continue;
        if (req[c] > max_req) max_req = req[c];
        n_active++;
    }
    if (n_active < 2 || max_req < 1e-6) {
        result->ss_KD_nM = 0;
        result->ss_Rmax = 0;
        result->ss_R2 = 0;
        return;
    }

    /* Cost function for isotherm */
    double best_ssr = 1e30;
    double best_KD = 1.0, best_Rmax = max_req;

    unsigned int ss_seed = (unsigned int)(time(NULL) ^ 0xA5A5A5A5u);
    int n_restarts = 100;

    for (int restart = 0; restart < n_restarts; restart++) {
        double log_KD, Rm;

        if (restart == 0 && kinetic_KD_nM > 0) {
            /* Seed first restart at the kinetic KD */
            log_KD = log10(kinetic_KD_nM);
            Rm = max_req * 1.2;
        } else if (restart == 1 && kinetic_KD_nM > 0) {
            /* Second restart: offset from kinetic KD */
            log_KD = log10(kinetic_KD_nM) + 0.5;
            Rm = max_req * 2.0;
        } else {
            /* Random starting point */
            log_KD = -2.0 + 10.0 * rand_r_double(&ss_seed);
            Rm = max_req * (0.5 + 19.5 * rand_r_double(&ss_seed));
        }

        /* Simplex vertices */
        double s[3][2] = {
            {log_KD, Rm},
            {log_KD + 0.5, Rm * 1.1},
            {log_KD - 0.5, Rm * 0.9}
        };
        double sv[3];

        for (int k = 0; k < 3; k++) {
            double KD_nM = pow(10.0, s[k][0]);
            double Rmax = s[k][1];
            double ssr = 0;
            if (Rmax <= 0 || KD_nM <= 0) { sv[k] = 1e30; continue; }
            for (int c = 0; c < data->ncycles; c++) {
                if (data->excluded[c]) continue;
                double C = data->cycles[c].conc_nM;
                double pred = Rmax * C / (KD_nM + C);
                double r = req[c] - pred;
                ssr += r * r;
            }
            sv[k] = ssr;
        }

        for (int iter = 0; iter < 5000; iter++) {
            /* Sort */
            for (int i = 0; i < 2; i++)
                for (int j = i+1; j < 3; j++)
                    if (sv[j] < sv[i]) {
                        double tmp[2]; memcpy(tmp, s[i], 16);
                        memcpy(s[i], s[j], 16); memcpy(s[j], tmp, 16);
                        double tv = sv[i]; sv[i] = sv[j]; sv[j] = tv;
                    }

            if (sv[2] - sv[0] < 1e-15 * (fabs(sv[0]) + 1e-30)) break;

            double c0 = (s[0][0] + s[1][0]) / 2;
            double c1 = (s[0][1] + s[1][1]) / 2;

            /* Reflect */
            double xr[2] = {c0 + (c0 - s[2][0]), c1 + (c1 - s[2][1])};
            double KD_nM = pow(10.0, xr[0]);
            double Rmax = xr[1];
            double vr = 0;
            if (Rmax <= 0 || KD_nM <= 0) { vr = 1e30; }
            else {
                for (int cc = 0; cc < data->ncycles; cc++) {
                    if (data->excluded[cc]) continue;
                    double C = data->cycles[cc].conc_nM;
                    double pred = Rmax * C / (KD_nM + C);
                    double r = req[cc] - pred;
                    vr += r * r;
                }
            }

            if (vr < sv[0]) {
                /* Expand */
                double xe[2] = {c0 + 2*(xr[0]-c0), c1 + 2*(xr[1]-c1)};
                KD_nM = pow(10.0, xe[0]); Rmax = xe[1];
                double ve = 0;
                if (Rmax <= 0 || KD_nM <= 0) { ve = 1e30; }
                else {
                    for (int cc = 0; cc < data->ncycles; cc++) {
                        if (data->excluded[cc]) continue;
                        double C = data->cycles[cc].conc_nM;
                        double pred = Rmax * C / (KD_nM + C);
                        double r = req[cc] - pred;
                        ve += r * r;
                    }
                }
                if (ve < vr) { memcpy(s[2], xe, 16); sv[2] = ve; }
                else { memcpy(s[2], xr, 16); sv[2] = vr; }
            } else if (vr < sv[1]) {
                memcpy(s[2], xr, 16); sv[2] = vr;
            } else {
                /* Contract */
                double xc[2] = {c0 + 0.5*(s[2][0]-c0), c1 + 0.5*(s[2][1]-c1)};
                KD_nM = pow(10.0, xc[0]); Rmax = xc[1];
                double vc = 0;
                if (Rmax <= 0 || KD_nM <= 0) { vc = 1e30; }
                else {
                    for (int cc = 0; cc < data->ncycles; cc++) {
                        if (data->excluded[cc]) continue;
                        double C = data->cycles[cc].conc_nM;
                        double pred = Rmax * C / (KD_nM + C);
                        double r = req[cc] - pred;
                        vc += r * r;
                    }
                }
                if (vc < sv[2]) { memcpy(s[2], xc, 16); sv[2] = vc; }
                else {
                    /* Shrink */
                    for (int i = 1; i < 3; i++) {
                        s[i][0] = s[0][0] + 0.5*(s[i][0] - s[0][0]);
                        s[i][1] = s[0][1] + 0.5*(s[i][1] - s[0][1]);
                        KD_nM = pow(10.0, s[i][0]); Rmax = s[i][1];
                        double ssr = 0;
                        if (Rmax <= 0 || KD_nM <= 0) { ssr = 1e30; }
                        else {
                            for (int cc = 0; cc < data->ncycles; cc++) {
                                if (data->excluded[cc]) continue;
                                double C = data->cycles[cc].conc_nM;
                                double pred = Rmax * C / (KD_nM + C);
                                double r = req[cc] - pred;
                                ssr += r * r;
                            }
                        }
                        sv[i] = ssr;
                    }
                }
            }
        }

        if (sv[0] < best_ssr) {
            best_ssr = sv[0];
            best_KD = pow(10.0, s[0][0]);
            best_Rmax = s[0][1];
        }
    }

    result->ss_KD_nM = best_KD;
    result->ss_Rmax = best_Rmax;

    /* Compute R² */
    double ss_ssr = 0, ss_sstot = 0, ss_mean = 0;
    int n = 0;
    for (int c = 0; c < data->ncycles; c++) {
        if (data->excluded[c]) continue;
        ss_mean += req[c];
        n++;
    }
    ss_mean /= (n > 0 ? n : 1);
    for (int c = 0; c < data->ncycles; c++) {
        if (data->excluded[c]) continue;
        double C = data->cycles[c].conc_nM;
        double pred = best_Rmax * C / (best_KD + C);
        double r = req[c] - pred;
        ss_ssr += r * r;
        double d = req[c] - ss_mean;
        ss_sstot += d * d;
    }
    result->ss_R2 = (ss_sstot > 0) ? 1.0 - ss_ssr / ss_sstot : 0.0;
}

/* ================================================================ */
/*  Parallel Restart Worker                                          */
/* ================================================================ */

typedef struct {
    int nparams;
    double lo[MAX_PARAMS], hi[MAX_PARAMS];
    int nrestarts;
    int nm_maxiter;
    double nm_tol;
    int fit_skip;
    unsigned int seed;
    double best_params[MAX_PARAMS];
    double best_val;
    double x0[MAX_PARAMS];  /* seed point for first restart */
    int has_x0;             /* whether to use x0 */
} RestartArg;

static void *restart_worker(void *arg) {
    RestartArg *ta = (RestartArg*)arg;
    tl_buf = (double*)malloc(sizeof(double) * MAX_POINTS);
    tl_skip = ta->fit_skip;
    ta->best_val = 1e30;

    for (int r = 0; r < ta->nrestarts; r++) {
        double best[MAX_PARAMS], best_val;
        /* Seed first restart with x0 if provided */
        const double *seed_pt = (r == 0 && ta->has_x0) ? ta->x0 : NULL;
        nelder_mead_ex(ta->nparams, best, &best_val,
                   ta->lo, ta->hi, ta->nm_maxiter, ta->nm_tol, &ta->seed,
                   seed_pt);
        if (best_val < ta->best_val) {
            memcpy(ta->best_params, best, sizeof(double) * ta->nparams);
            ta->best_val = best_val;
        }
    }

    free(tl_buf);
    tl_buf = NULL;
    return NULL;
}

/* ================================================================ */
/*  Apply RI + Drift to simulation buffer (for output generation)    */
/* ================================================================ */

void apply_local_corrections(const double *params, const ParamLayout *lay,
                             int cycle_idx, const Cycle *cy, double *buf,
                             double t_assoc_end) {
    int ai = lay->cycle_to_active[cycle_idx];
    if (ai < 0) return;  /* excluded cycle */

    double ri = (lay->off_ri >= 0) ? params[lay->off_ri + ai] : 0.0;
    double drift = (lay->off_drift >= 0) ? params[lay->off_drift] : 0.0;

    if (ri == 0.0 && drift == 0.0) return;

    for (int i = cy->fit_start; i < cy->npoints; i++) {
        if (ri != 0.0 && cy->time[i] >= 0.0 && cy->time[i] < t_assoc_end)
            buf[i] += ri;
        if (drift != 0.0 && cy->time[i] >= 0.0)
            buf[i] += drift * cy->time[i];
    }
}

/* ================================================================ */
/*  Main Fitting Pipeline                                            */
/* ================================================================ */

int optim_fit(SPRData *data, const FitConfig *cfg, FitResult *result) {
    memset(result, 0, sizeof(FitResult));
    result->model = cfg->model;
    result->nparams = MODEL_NPARAMS[cfg->model];

    optim_setup(data, cfg);
    srand((unsigned)time(NULL) ^ (unsigned)(size_t)data->filename);

    /* Build parameter layout */
    g_layout = build_param_layout(cfg->model, &cfg->advanced, data);
    result->layout = g_layout;
    int np = g_layout.total;

    int has_local = (g_layout.off_ri >= 0 || g_layout.off_drift >= 0);
    int has_tc = (g_layout.off_tc >= 0);

    /* Estimate Rmax from max observed response */
    double max_resp = 0;
    for (int c = 0; c < data->ncycles; c++) {
        if (data->excluded[c]) continue;
        for (int i = data->cycles[c].fit_start; i < data->cycles[c].npoints; i++)
            if (!data->cycles[c].skip[i] && data->cycles[c].resp[i] > max_resp)
                max_resp = data->cycles[c].resp[i];
    }

    double rmax_hi = (max_resp < 0.5) ? 5.0 : max_resp * 3.0;
    double rmax_lo = (max_resp < 0.5) ? 0.01 : max_resp * 0.1;

    /* Set search bounds based on model */
    int ng = g_layout.n_global;
    double lo[MAX_PARAMS], hi[MAX_PARAMS];

    if (cfg->model == MODEL_LANGMUIR) {
        lo[0] = 2.0;  hi[0] = 8.0;    /* log10(ka) */
        lo[1] = -6.0; hi[1] = 2.0;    /* log10(kd) */
        lo[2] = rmax_lo; hi[2] = rmax_hi; /* Rmax */
    } else if (cfg->model == MODEL_HETEROGENEOUS) {
        lo[0] = 2.0;  hi[0] = 8.0;    /* log10(ka1) */
        lo[1] = -6.0; hi[1] = 2.0;    /* log10(kd1) */
        lo[2] = rmax_lo * 0.1; hi[2] = rmax_hi; /* Rmax1 */
        lo[3] = 1.0;  hi[3] = 7.0;    /* log10(ka2) */
        lo[4] = -4.0; hi[4] = 2.0;    /* log10(kd2) */
        lo[5] = rmax_lo * 0.1; hi[5] = rmax_hi; /* Rmax2 */
    } else if (cfg->model == MODEL_TWOSTATE) {
        lo[0] = 2.0;  hi[0] = 8.0;    /* log10(ka1) */
        lo[1] = -6.0; hi[1] = 2.0;    /* log10(kd1) */
        lo[2] = -5.0; hi[2] = 1.0;    /* log10(ka2) */
        lo[3] = -6.0; hi[3] = 1.0;    /* log10(kd2) */
        lo[4] = rmax_lo; hi[4] = rmax_hi; /* Rmax */
    }

    /* Bounds for local parameters */
    if (g_layout.off_ri >= 0) {
        for (int a = 0; a < g_layout.n_active; a++) {
            lo[g_layout.off_ri + a] = -5.0;    /* RI: -5 to +5 RU */
            hi[g_layout.off_ri + a] = 5.0;
        }
    }
    if (g_layout.off_drift >= 0) {
        lo[g_layout.off_drift] = -0.1;  /* Drift: -0.1 to +0.1 RU/s (single global) */
        hi[g_layout.off_drift] = 0.1;
    }
    if (g_layout.off_tc >= 0) {
        lo[g_layout.off_tc] = 4.0;     /* log10(tc): 1e4 to 1e10 RU·M⁻¹·s⁻¹ */
        hi[g_layout.off_tc] = 10.0;
    }

    double global_best[MAX_PARAMS] = {0};
    double global_best_val = 1e30;

    fprintf(stderr, "Fitting %s model (Rmax range: %.2f - %.2f)",
            MODEL_NAMES[cfg->model], rmax_lo, rmax_hi);
    if (has_local || has_tc)
        fprintf(stderr, " [advanced: %s%s%s, %d params]",
                g_layout.off_ri >= 0 ? "RI " : "",
                g_layout.off_drift >= 0 ? "Drift " : "",
                g_layout.off_tc >= 0 ? "TC" : "",
                np);
    fprintf(stderr, "...\n");

    /* ── Stage 1: Fit global params only (if advanced params present) ── */
    if (has_local || has_tc) {
        /* First pass: global-only search (set local params to 0) */
        int np_global = ng;
        int nrestarts_s1 = cfg->nrestarts;

        int ncpus = 8;
        int nthreads = ncpus;
        if (nthreads > nrestarts_s1) nthreads = nrestarts_s1;

        pthread_t threads[32];
        RestartArg targs[32];

        int restarts_per = nrestarts_s1 / nthreads;
        int remainder = nrestarts_s1 % nthreads;

        /* Temporarily set layout to global-only for stage 1 */
        ParamLayout saved_layout = g_layout;
        g_layout.off_ri = -1;
        g_layout.off_drift = -1;
        g_layout.off_tc = -1;
        g_layout.total = ng;

        for (int t = 0; t < nthreads; t++) {
            targs[t].nparams = np_global;
            memcpy(targs[t].lo, lo, sizeof(double) * np_global);
            memcpy(targs[t].hi, hi, sizeof(double) * np_global);
            targs[t].nrestarts = restarts_per + (t < remainder ? 1 : 0);
            targs[t].nm_maxiter = cfg->nm_maxiter;
            targs[t].nm_tol = cfg->nm_tol;
            targs[t].fit_skip = cfg->fit_skip;
            targs[t].seed = (unsigned int)(time(NULL) ^ (t * 2654435761u));
            targs[t].best_val = 1e30;
            targs[t].has_x0 = 0;
            pthread_create(&threads[t], NULL, restart_worker, &targs[t]);
        }

        for (int t = 0; t < nthreads; t++) {
            pthread_join(threads[t], NULL);
            if (targs[t].best_val < global_best_val) {
                memcpy(global_best, targs[t].best_params, sizeof(double) * np_global);
                global_best_val = targs[t].best_val;
            }
        }

        /* Restore full layout */
        g_layout = saved_layout;

        /* Initialize advanced params after stage 1 */
        for (int j = ng; j < np; j++)
            global_best[j] = 0.0;
        /* tc starts at midpoint of log10 range */
        if (g_layout.off_tc >= 0)
            global_best[g_layout.off_tc] = 7.0;  /* log10(1e7) */

        fprintf(stderr, "  Stage 1 (global only): SSR=%.4f\n", global_best_val);
        /* Recompute SSR with full layout (including advanced params at init)
         * so Stage 2 comparison uses the same cost function */
        {
            double *saved_buf = tl_buf;
            int saved_skip = tl_skip;
            tl_buf = (double*)malloc(sizeof(double) * MAX_POINTS);
            tl_skip = 1;
            global_best_val = cost(global_best);
            free(tl_buf);
            tl_buf = saved_buf;
            tl_skip = saved_skip;
        }

        /* ── Stage 2: Refine all params together (parallel) ── */
        tl_skip = 1;
        {
            double lo2[MAX_PARAMS], hi2[MAX_PARAMS];

            /* Tight bounds around global params */
            for (int j = 0; j < ng; j++) {
                if (j == 2 || j == 5 || (cfg->model == MODEL_TWOSTATE && j == 4)) {
                    lo2[j] = global_best[j] * 0.7;
                    hi2[j] = global_best[j] * 1.3;
                } else {
                    lo2[j] = global_best[j] - 0.5;
                    hi2[j] = global_best[j] + 0.5;
                }
            }
            /* Full range for local params */
            for (int j = ng; j < np; j++) {
                lo2[j] = lo[j];
                hi2[j] = hi[j];
            }

            /* Seed one vertex at stage 1 best + initialized advanced params */
            /* Then parallel restarts explore around it */
            int nrestarts_s2 = 16;
            int nthreads_s2 = 8;
            if (nthreads_s2 > nrestarts_s2) nthreads_s2 = nrestarts_s2;

            int rps2 = nrestarts_s2 / nthreads_s2;
            int rem2 = nrestarts_s2 % nthreads_s2;

            for (int t = 0; t < nthreads_s2; t++) {
                targs[t].nparams = np;
                memcpy(targs[t].lo, lo2, sizeof(double) * np);
                memcpy(targs[t].hi, hi2, sizeof(double) * np);
                targs[t].nrestarts = rps2 + (t < rem2 ? 1 : 0);
                targs[t].nm_maxiter = cfg->nm_maxiter * 2;
                targs[t].nm_tol = cfg->nm_tol;
                targs[t].fit_skip = 1;
                targs[t].seed = (unsigned int)(time(NULL) ^ ((t + 100) * 2654435761u));
                targs[t].best_val = 1e30;
                /* Seed every thread with stage 1 result + initialized advanced params */
                memcpy(targs[t].x0, global_best, sizeof(double) * np);
                targs[t].has_x0 = 1;
                pthread_create(&threads[t], NULL, restart_worker, &targs[t]);
            }

            for (int t = 0; t < nthreads_s2; t++) {
                pthread_join(threads[t], NULL);
                if (targs[t].best_val < global_best_val) {
                    memcpy(global_best, targs[t].best_params, sizeof(double) * np);
                    global_best_val = targs[t].best_val;
                }
            }
        }

        fprintf(stderr, "  Stage 2 (all params): SSR=%.4f\n",
                cost(global_best));

    } else {
        /* ── Standard fitting: no local params ── */
        int ncpus = 8;
        int nthreads = ncpus;
        if (nthreads > cfg->nrestarts) nthreads = cfg->nrestarts;

        pthread_t threads[32];
        RestartArg targs[32];

        int restarts_per = cfg->nrestarts / nthreads;
        int remainder = cfg->nrestarts % nthreads;

        for (int t = 0; t < nthreads; t++) {
            targs[t].nparams = np;
            memcpy(targs[t].lo, lo, sizeof(double) * np);
            memcpy(targs[t].hi, hi, sizeof(double) * np);
            targs[t].nrestarts = restarts_per + (t < remainder ? 1 : 0);
            targs[t].nm_maxiter = cfg->nm_maxiter;
            targs[t].nm_tol = cfg->nm_tol;
            targs[t].fit_skip = cfg->fit_skip;
            targs[t].seed = (unsigned int)(time(NULL) ^ (t * 2654435761u));
            targs[t].best_val = 1e30;
            targs[t].has_x0 = 0;
            pthread_create(&threads[t], NULL, restart_worker, &targs[t]);
        }

        for (int t = 0; t < nthreads; t++) {
            pthread_join(threads[t], NULL);
            if (targs[t].best_val < global_best_val) {
                memcpy(global_best, targs[t].best_params, sizeof(double) * np);
                global_best_val = targs[t].best_val;
            }
        }

        fprintf(stderr, "  Coarse search: SSR=%.4f", global_best_val);
        if (cfg->model == MODEL_LANGMUIR)
            fprintf(stderr, "  ka=%.3e  kd=%.3e  Rmax=%.3f",
                    pow(10, global_best[0]), pow(10, global_best[1]), global_best[2]);
        fprintf(stderr, " (%d threads)\n", nthreads);

        /* Refinement with all data points */
        tl_skip = 1;
        {
            unsigned int ref_seed = (unsigned int)time(NULL);
            double lo2[MAX_PARAMS], hi2[MAX_PARAMS];
            for (int j = 0; j < np; j++) {
                if (j == 2 || j == 5 || (cfg->model == MODEL_TWOSTATE && j == 4)) {
                    lo2[j] = global_best[j] * 0.8;
                    hi2[j] = global_best[j] * 1.2;
                } else {
                    lo2[j] = global_best[j] - 0.5;
                    hi2[j] = global_best[j] + 0.5;
                }
            }
            for (int r = 0; r < 5; r++) {
                double best[MAX_PARAMS], best_val;
                nelder_mead(np, best, &best_val, lo2, hi2, cfg->nm_maxiter, cfg->nm_tol, &ref_seed);
                if (best_val < cost(global_best)) {
                    memcpy(global_best, best, sizeof(double) * np);
                    global_best_val = best_val;
                }
            }
        }
    }

    /* Store results */
    memcpy(result->params, global_best, sizeof(double) * np);
    result->nparams = np;

    if (cfg->model == MODEL_LANGMUIR) {
        result->ka = pow(10.0, global_best[0]);
        result->kd = pow(10.0, global_best[1]);
        result->Rmax = global_best[2];
        result->KD_M = result->kd / result->ka;
    } else if (cfg->model == MODEL_HETEROGENEOUS) {
        result->ka = pow(10.0, global_best[0]);
        result->kd = pow(10.0, global_best[1]);
        result->Rmax = global_best[2];
        result->KD_M = result->kd / result->ka;
        result->ka2 = pow(10.0, global_best[3]);
        result->kd2 = pow(10.0, global_best[4]);
        result->Rmax2 = global_best[5];
        result->KD2_M = result->kd2 / result->ka2;
    } else if (cfg->model == MODEL_TWOSTATE) {
        result->ka = pow(10.0, global_best[0]);
        result->kd = pow(10.0, global_best[1]);
        result->ka2 = pow(10.0, global_best[2]);
        result->kd2 = pow(10.0, global_best[3]);
        result->Rmax = global_best[4];
        result->KD_M = result->kd / result->ka;
    }

    /* Extract local parameters */
    for (int c = 0; c < data->ncycles; c++) {
        int ai = g_layout.cycle_to_active[c];
        result->ri[c] = (ai >= 0 && g_layout.off_ri >= 0) ? global_best[g_layout.off_ri + ai] : 0.0;
    }
    result->drift = (g_layout.off_drift >= 0) ? global_best[g_layout.off_drift] : 0.0;
    result->tc = (g_layout.off_tc >= 0) ? pow(10.0, global_best[g_layout.off_tc]) : 0.0;

    /* Compute statistics */
    double ssr = 0, ss_tot = 0, mean_resp = 0;
    int n_total = 0;

    for (int c = 0; c < data->ncycles; c++) {
        if (data->excluded[c]) continue;
        for (int i = data->cycles[c].fit_start; i < data->cycles[c].npoints; i++) {
            if (data->cycles[c].skip[i]) continue;
            mean_resp += data->cycles[c].resp[i];
            n_total++;
        }
    }
    mean_resp /= (n_total > 0 ? n_total : 1);

    for (int c = 0; c < data->ncycles; c++) {
        if (data->excluded[c]) continue;
        if (data->mode == MODE_MCK) {
            simulate_mck_cycle(cfg->model, global_best, &data->cycles[c], tl_buf,
                              cfg->t_assoc_end, cfg->dt, result->tc);
            apply_local_corrections(global_best, &g_layout, c, &data->cycles[c],
                                   tl_buf, cfg->t_assoc_end);
        } else {
            simulate_sck_trace(cfg->model, global_best, data, tl_buf, cfg->dt, result->tc);
            apply_local_corrections(global_best, &g_layout, c, &data->cycles[c],
                                   tl_buf, cfg->t_assoc_end);
        }

        for (int i = data->cycles[c].fit_start; i < data->cycles[c].npoints; i++) {
            if (data->cycles[c].skip[i]) continue;
            double r = data->cycles[c].resp[i] - tl_buf[i];
            ssr += r * r;
            double d = data->cycles[c].resp[i] - mean_resp;
            ss_tot += d * d;
        }
    }

    result->ssr = ssr;
    result->R2 = (ss_tot > 0) ? 1.0 - ssr / ss_tot : 0.0;
    result->chi2 = (n_total > np) ? ssr / (n_total - np) : ssr;
    result->rms = sqrt(ssr / (n_total > 0 ? n_total : 1));
    result->n_points = n_total;

    fprintf(stderr, "Done: ka=%.4e  kd=%.4e  KD=%.4e M  Rmax=%.3f RU\n",
            result->ka, result->kd, result->KD_M, result->Rmax);
    if (cfg->model == MODEL_HETEROGENEOUS)
        fprintf(stderr, "      ka2=%.4e  kd2=%.4e  KD2=%.4e M  Rmax2=%.3f RU\n",
                result->ka2, result->kd2, result->KD2_M, result->Rmax2);
    if (cfg->model == MODEL_TWOSTATE)
        fprintf(stderr, "      ka2=%.4e  kd2=%.4e\n", result->ka2, result->kd2);
    if (has_local || has_tc) {
        fprintf(stderr, "  Advanced params:");
        if (g_layout.off_ri >= 0) {
            fprintf(stderr, " RI=[");
            for (int a = 0; a < g_layout.n_active; a++)
                fprintf(stderr, "%s%.3f", a > 0 ? "," : "", result->ri[g_layout.active_map[a]]);
            fprintf(stderr, "]");
        }
        if (g_layout.off_drift >= 0) {
            fprintf(stderr, " Drift=%.6f", result->drift);
        }
        if (g_layout.off_tc >= 0) {
            fprintf(stderr, " tc=%.4e", result->tc);
        }
        fprintf(stderr, "\n");
    }

    /* Steady-state affinity (MCK only) */
    if (data->mode == MODE_MCK && data->ncycles > 1) {
        double req[MAX_CYCLES];
        double req_hi = cfg->t_assoc_end - 2.0;
        double req_lo = req_hi - 15.0;
        if (req_lo < 5.0) req_lo = 5.0;
        compute_req(data, req, req_lo, req_hi);
        memcpy(result->ss_req, req, sizeof(double) * data->ncycles);
        double kinetic_KD_nM = result->KD_M * 1e9;
        ss_fit(data, req, kinetic_KD_nM, result);
        fprintf(stderr, "Steady-state: KD=%.2f nM  Rmax=%.3f RU  R2=%.5f\n",
                result->ss_KD_nM, result->ss_Rmax, result->ss_R2);
    }

    /* U-value: parameter uniqueness via condition number of J^T·J
     * U = sqrt(λ_max / λ_min) of the normal matrix.
     * Only uses global model parameters (ka, kd, Rmax etc.).
     * U < 15: parameters well-determined; U > 25: significantly correlated. */
    {
        int ng = MODEL_NPARAMS[cfg->model];  /* global params only */
        double tc_val = result->tc;

        /* Compute J^T·J by numerical Jacobian (finite differences on SSR) */
        /* JtJ[i][j] = sum_k (dr_k/dp_i * dr_k/dp_j) */
        /* We approximate via: JtJ[i][j] ≈ (SSR(+i,+j) - SSR(+i,-j) - SSR(-i,+j) + SSR(-i,-j)) / (4*hi*hj) */
        /* But more efficiently: compute residual vectors and form J^T·J directly */

        /* Collect residual indices and count */
        int n_res = 0;
        for (int c = 0; c < data->ncycles; c++) {
            if (data->excluded[c]) continue;
            const Cycle *cy = &data->cycles[c];
            for (int i = cy->fit_start; i < cy->npoints; i++)
                if (!cy->skip[i]) n_res++;
        }

        /* Allocate Jacobian columns (ng columns, n_res rows each) */
        double **J = (double**)malloc(sizeof(double*) * ng);
        for (int j = 0; j < ng; j++)
            J[j] = (double*)calloc(n_res, sizeof(double));

        double *sim0 = (double*)malloc(sizeof(double) * MAX_POINTS);
        double *sim1 = (double*)malloc(sizeof(double) * MAX_POINTS);
        double params_pert[MAX_PARAMS];
        memcpy(params_pert, global_best, sizeof(double) * np);

        /* For each global parameter, compute column of Jacobian */
        for (int p = 0; p < ng; p++) {
            double h = fabs(global_best[p]) * 1e-5;
            if (h < 1e-10) h = 1e-10;

            memcpy(params_pert, global_best, sizeof(double) * np);
            params_pert[p] = global_best[p] + h;

            int ri = 0;
            for (int c = 0; c < data->ncycles; c++) {
                if (data->excluded[c]) continue;
                const Cycle *cy = &data->cycles[c];

                /* Simulate at optimum and perturbed */
                if (data->mode == MODE_MCK) {
                    simulate_mck_cycle(cfg->model, global_best, cy, sim0,
                                      cfg->t_assoc_end, cfg->dt, tc_val);
                    apply_local_corrections(global_best, &g_layout, c, cy, sim0,
                                           cfg->t_assoc_end);
                    simulate_mck_cycle(cfg->model, params_pert, cy, sim1,
                                      cfg->t_assoc_end, cfg->dt, tc_val);
                    apply_local_corrections(global_best, &g_layout, c, cy, sim1,
                                           cfg->t_assoc_end);
                } else {
                    simulate_sck_trace(cfg->model, global_best, data, sim0,
                                      cfg->dt, tc_val);
                    simulate_sck_trace(cfg->model, params_pert, data, sim1,
                                      cfg->dt, tc_val);
                }

                for (int i = cy->fit_start; i < cy->npoints; i++) {
                    if (cy->skip[i]) continue;
                    /* J[p][ri] = d(residual)/d(param) = -(d(sim)/d(param)) */
                    J[p][ri] = -(sim1[i] - sim0[i]) / h;
                    ri++;
                }
            }
        }

        /* Normalize each Jacobian column to unit length so J^T·J becomes
         * a correlation matrix (removes parameter scale effects) */
        for (int p = 0; p < ng; p++) {
            double norm = 0;
            for (int k = 0; k < n_res; k++)
                norm += J[p][k] * J[p][k];
            norm = sqrt(norm);
            if (norm > 1e-30)
                for (int k = 0; k < n_res; k++)
                    J[p][k] /= norm;
        }

        /* Form J^T·J (ng × ng correlation matrix) */
        double JtJ[6][6] = {{0}};  /* max 6 global params */
        for (int i = 0; i < ng; i++)
            for (int j = i; j < ng; j++) {
                double s = 0;
                for (int k = 0; k < n_res; k++)
                    s += J[i][k] * J[j][k];
                JtJ[i][j] = s;
                JtJ[j][i] = s;
            }

        /* Eigenvalues via Jacobi iteration (for small symmetric matrix) */
        double A[6][6];
        memcpy(A, JtJ, sizeof(A));

        for (int iter = 0; iter < 200; iter++) {
            /* Find largest off-diagonal */
            int pi = 0, pj = 1;
            double maxoff = 0;
            for (int i = 0; i < ng; i++)
                for (int j = i+1; j < ng; j++)
                    if (fabs(A[i][j]) > maxoff) {
                        maxoff = fabs(A[i][j]);
                        pi = i; pj = j;
                    }
            if (maxoff < 1e-20) break;

            /* Jacobi rotation */
            double theta = 0.5 * atan2(2.0 * A[pi][pj],
                                        A[pi][pi] - A[pj][pj]);
            double c_r = cos(theta), s_r = sin(theta);

            double Ai[6], Aj[6];
            for (int k = 0; k < ng; k++) {
                Ai[k] =  c_r * A[pi][k] + s_r * A[pj][k];
                Aj[k] = -s_r * A[pi][k] + c_r * A[pj][k];
            }
            for (int k = 0; k < ng; k++) {
                A[pi][k] = Ai[k]; A[k][pi] = Ai[k];
                A[pj][k] = Aj[k]; A[k][pj] = Aj[k];
            }
            double app =  c_r * Ai[pi] + s_r * Ai[pj];
            double aqq = -s_r * Aj[pi] + c_r * Aj[pj];
            A[pi][pi] = app;
            A[pj][pj] = aqq;
            A[pi][pj] = 0; A[pj][pi] = 0;
        }

        /* Extract eigenvalues (diagonal of A) */
        double lam_min = 1e30, lam_max = 0;
        for (int i = 0; i < ng; i++) {
            double lam = fabs(A[i][i]);
            if (lam > lam_max) lam_max = lam;
            if (lam < lam_min) lam_min = lam;
        }

        result->u_value = (lam_min > 1e-30) ? sqrt(lam_max / lam_min) : 999.0;

        fprintf(stderr, "U-value: %.1f", result->u_value);
        if (result->u_value < 15)
            fprintf(stderr, " (good: parameters well-determined)\n");
        else if (result->u_value < 25)
            fprintf(stderr, " (borderline)\n");
        else
            fprintf(stderr, " (high: parameters may be correlated)\n");

        for (int j = 0; j < ng; j++) free(J[j]);
        free(J);
        free(sim0);
        free(sim1);
    }

    optim_cleanup();
    return 0;
}
