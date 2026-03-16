/*
 * spr_types.h - Common data structures for SPR fitting
 *
 * Supports both MCK (Multi-Cycle Kinetics) and SCK (Single-Cycle Kinetics).
 */

#ifndef SPR_TYPES_H
#define SPR_TYPES_H

#define MAX_CYCLES      30
#define MAX_POINTS      2500
#define MAX_PARAMS      8
#define MAX_INJECTIONS  10
#define DISPLAY_SKIP    5      /* Output every Nth point for display */

/* ── Data Mode ─────────────────────────────────────────────────── */

typedef enum {
    MODE_MCK = 0,   /* Multi-Cycle Kinetics: separate cycles, global fit */
    MODE_SCK = 1    /* Single-Cycle Kinetics: staircase injection, single trace */
} DataMode;

/* ── Binding Models ────────────────────────────────────────────── */

typedef enum {
    MODEL_LANGMUIR      = 0,   /* 1:1 Langmuir: ka, kd, Rmax (3 params) */
    MODEL_HETEROGENEOUS = 1,   /* Two independent sites: ka1,kd1,Rmax1,ka2,kd2,Rmax2 (6 params) */
    MODEL_TWOSTATE      = 2    /* Conformational change: ka1,kd1,ka2,kd2,Rmax (5 params) */
} ModelType;

static const char *MODEL_NAMES[] = {
    "langmuir", "heterogeneous", "twostate"
};

static const int MODEL_NPARAMS[] = { 3, 6, 5 };

/* ── Per-Cycle Data ────────────────────────────────────────────── */

typedef struct {
    double time[MAX_POINTS];
    double resp[MAX_POINTS];
    double inst_resp[MAX_POINTS];   /* Instrument-fitted response */
    int    skip[MAX_POINTS];        /* 1 = exclude from fitting (artifact) */
    int    npoints;
    double conc_nM;
    double baseline;
    int    fit_start;               /* First index where time >= 0 */
} Cycle;

/* ── Dataset Container ─────────────────────────────────────────── */

typedef struct {
    Cycle    cycles[MAX_CYCLES];
    int      ncycles;
    char     filename[512];
    DataMode mode;

    /* SCK-specific: injection schedule */
    double   inj_conc_nM[MAX_INJECTIONS];   /* Concentration per injection */
    double   inj_start[MAX_INJECTIONS];     /* Start time of each injection */
    double   inj_end[MAX_INJECTIONS];       /* End time of each injection */
    int      n_injections;

    /* Cycle exclusion mask (set by CLI --exclude) */
    int      excluded[MAX_CYCLES];          /* 1 = excluded from fitting */
} SPRData;

/* ── Fit Configuration ─────────────────────────────────────────── */

typedef struct {
    ModelType model;
    int       nrestarts;        /* Number of random restarts */
    int       nm_maxiter;       /* Nelder-Mead max iterations */
    double    nm_tol;           /* Convergence tolerance */
    int       fit_skip;         /* Use every Nth point during fitting */
    double    guard_time;       /* Seconds to skip after artifact boundaries */
    double    t_assoc_end;      /* End of association phase (default 60s) */
    int       t_assoc_end_set;  /* 1 = user explicitly set, 0 = auto-detect */
    double    dt;               /* Max RK4 step size */
} FitConfig;

/* Default configuration */
static inline FitConfig fit_config_default(void) {
    FitConfig cfg;
    cfg.model      = MODEL_LANGMUIR;
    cfg.nrestarts  = 40;
    cfg.nm_maxiter = 3000;
    cfg.nm_tol     = 1e-12;
    cfg.fit_skip   = 3;
    cfg.guard_time = 2.0;
    cfg.t_assoc_end = 60.0;
    cfg.t_assoc_end_set = 0;
    cfg.dt         = 0.2;
    return cfg;
}

/* ── Fit Result ────────────────────────────────────────────────── */

typedef struct {
    ModelType model;
    double    params[MAX_PARAMS];   /* Raw fitted parameters */
    int       nparams;

    /* Derived kinetic values */
    double    ka, kd, KD_M;            /* For langmuir / primary site */
    double    Rmax;
    double    ka2, kd2, KD2_M, Rmax2;  /* For heterogeneous / two-state */

    /* Fit quality */
    double    ssr, R2, chi2, rms;
    int       n_points;

    /* Steady-state affinity */
    double    ss_KD_nM, ss_Rmax, ss_R2;
    double    ss_req[MAX_CYCLES];
} FitResult;

#endif /* SPR_TYPES_H */
