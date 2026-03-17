/*
 * spr_io.c - Data I/O for SPR fitting (MCK and SCK formats)
 */

#include "spr_io.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ================================================================ */
/*  Helpers                                                          */
/* ================================================================ */

static double parse_conc(const char *field) {
    const char *p = strstr(field, "Conc ");
    if (!p) return -1.0;
    return atof(p + 5);
}

/* Parse multiple concentrations from SCK header: "Conc 3.9 15.6 62.5 ... nM" */
static int parse_sck_concs(const char *field, double *concs, int max_concs) {
    const char *p = strstr(field, "Conc ");
    if (!p) return 0;
    p += 5;

    int n = 0;
    while (n < max_concs) {
        char *end;
        double v = strtod(p, &end);
        if (end == p) break;
        /* Stop if we hit "nM" or other non-numeric */
        if (*end != ' ' && *end != '\t' && *end != '\0' &&
            *end != '.' && (*end < '0' || *end > '9')) {
            /* Check if it's a number (could be negative exponent etc.) */
            if (*end == 'n' || *end == 'N') { /* "nM" unit marker */
                concs[n++] = v;
                break;
            }
        }
        concs[n++] = v;
        p = end;
        while (*p == ' ' || *p == '\t') p++;
        if (*p == '\0' || *p == '\n' || *p == '\r') break;
    }
    return n;
}

static void preprocess_cycle(Cycle *cy, double guard_time) {
    /* Pass 1: Mark exact-zero points as artifacts */
    for (int j = 0; j < cy->npoints; j++)
        cy->skip[j] = (cy->resp[j] == 0.0) ? 1 : 0;

    /* Pass 2: Guard band after zero-regions */
    for (int j = 1; j < cy->npoints; j++) {
        if (cy->resp[j] != 0.0 && cy->resp[j-1] == 0.0) {
            double guard_end = cy->time[j] + guard_time;
            for (int k = j; k < cy->npoints && cy->time[k] < guard_end; k++)
                cy->skip[k] = 1;
        }
    }

    /* Find fit_start (first point where t >= 0) */
    cy->fit_start = cy->npoints;
    for (int j = 0; j < cy->npoints; j++) {
        if (cy->time[j] >= 0.0) { cy->fit_start = j; break; }
    }

    /* Compute baseline from non-artifact points with t < 0 */
    double sum = 0.0;
    int count = 0;
    for (int j = 0; j < cy->npoints; j++) {
        if (cy->time[j] < 0.0 && !cy->skip[j]) {
            sum += cy->resp[j];
            count++;
        }
    }
    cy->baseline = (count > 0) ? sum / count : 0.0;

    /* Subtract baseline */
    for (int j = 0; j < cy->npoints; j++)
        cy->resp[j] -= cy->baseline;

    /* Log */
    int n_skip = 0;
    for (int j = cy->fit_start; j < cy->npoints; j++)
        if (cy->skip[j]) n_skip++;
    fprintf(stderr, "  Cycle %d (%.2f nM): %d/%d points excluded as artifacts\n",
            0, cy->conc_nM, n_skip, cy->npoints - cy->fit_start);
}

void preprocess_data(SPRData *data, double guard_time) {
    for (int i = 0; i < data->ncycles; i++)
        preprocess_cycle(&data->cycles[i], guard_time);
}

/* ================================================================ */
/*  MCK Format Reader                                                */
/* ================================================================ */

/* Check if a column header ends with "_X" (or "_X " with trailing space).
 * This identifies the time column of a data pair. */
static int ends_with_X(const char *header) {
    int len = (int)strlen(header);
    /* Trim trailing whitespace */
    while (len > 0 && (header[len-1] == ' ' || header[len-1] == '\r'
                        || header[len-1] == '\n')) len--;
    return (len >= 2 && header[len-2] == '_' && header[len-1] == 'X');
}

/* Check if a column header is an instrument-fitted column.
 * Fitted columns have "; Fitted_X" or "; Fitted_Y" as part of their name.
 * These contain Biacore's own fit curves, not raw experimental data.
 * Example data:   "Run 1; Cy 5; Conc 100 nM_X"
 * Example fitted:  "Run 1; Cy 5; Conc 100 nM; Fitted_X" */
static int is_fitted_column(const char *header) {
    return (strstr(header, "; Fitted_") != NULL ||
            strstr(header, ";Fitted_") != NULL) ? 1 : 0;
}

/* Check if a column header marks a cycle as excluded.
 * Example: "Run 1; Cy 15; Conc 100000 nM; Excluded_X" */
static int is_excluded_cycle(const char *header) {
    return (strstr(header, "Excluded") != NULL ||
            strstr(header, "excluded") != NULL) ? 1 : 0;
}

int read_mck_data(const char *filename, SPRData *data, const FitConfig *cfg) {
    FILE *f = fopen(filename, "r");
    if (!f) { fprintf(stderr, "Error: cannot open %s\n", filename); return -1; }

    strncpy(data->filename, filename, 511);
    data->filename[511] = '\0';
    data->mode = MODE_MCK;

    char line[131072];
    if (!fgets(line, sizeof(line), f)) { fclose(f); return -1; }

    /* Split header by tabs */
    char *hcopy = strdup(line);
    char *headers[512];
    int nh = 0;
    {
        char *p = hcopy;
        while (nh < 512) {
            char *tab = strchr(p, '\t');
            if (tab) { *tab = '\0'; headers[nh++] = p; p = tab + 1; }
            else     { headers[nh++] = p; break; }
        }
    }

    /* Scan headers column-by-column to find experimental data columns.
     *
     * Biacore export columns are NOT always in uniform 4-column groups.
     * Excluded cycles may have only 2 columns (X, Y without Fitted pair).
     * We identify data columns by:
     *   1. Header ends with "_X" (time column of a pair)
     *   2. Header does NOT contain "; Fitted_" (not an instrument-fit column)
     *   3. Next column (col+1) is the corresponding _Y response column
     * The Fitted_X/Fitted_Y columns for each cycle (instrument fit) are
     * stored as inst_resp but are NOT treated as separate experimental cycles.
     */
    int data_col[MAX_CYCLES];      /* Column index of the _X column for each cycle */
    int fit_col[MAX_CYCLES];       /* Column index of the Fitted_X column (-1 if none) */
    int cycle_excluded[MAX_CYCLES];
    double cycle_conc[MAX_CYCLES];
    int ncycles = 0;

    for (int c = 0; c < nh && ncycles < MAX_CYCLES; c++) {
        if (!ends_with_X(headers[c])) continue;      /* Not an _X column */
        if (is_fitted_column(headers[c])) continue;   /* Instrument fit, skip */

        /* This is an experimental data _X column */
        double conc = parse_conc(headers[c]);
        if (conc < 0) continue;

        data_col[ncycles] = c;
        fit_col[ncycles] = -1;  /* Look for corresponding Fitted_X */
        cycle_conc[ncycles] = conc;
        cycle_excluded[ncycles] = is_excluded_cycle(headers[c]);

        /* Search following columns for the matching "; Fitted_X" */
        for (int fc = c + 2; fc < nh && fc < c + 4; fc++) {
            if (is_fitted_column(headers[fc]) && ends_with_X(headers[fc])) {
                fit_col[ncycles] = fc;
                break;
            }
        }

        ncycles++;
    }

    /* Track unique concentrations for replicate logging */
    double seen_concs[MAX_CYCLES];
    int    seen_count[MAX_CYCLES];
    int    n_unique_concs = 0;

    for (int i = 0; i < ncycles; i++) {
        data->cycles[i].conc_nM = cycle_conc[i];
        data->cycles[i].npoints = 0;

        if (cycle_excluded[i]) {
            data->excluded[i] = 1;
            fprintf(stderr, "  Cycle %d (%.2f nM): auto-excluded from Biacore header\n",
                    i + 1, cycle_conc[i]);
        }

        /* Track duplicates */
        int rep_num = 1;
        int found = 0;
        for (int k = 0; k < n_unique_concs; k++) {
            if (fabs(seen_concs[k] - cycle_conc[i]) < cycle_conc[i] * 1e-4) {
                seen_count[k]++;
                rep_num = seen_count[k];
                found = 1;
                break;
            }
        }
        if (!found && n_unique_concs < MAX_CYCLES) {
            seen_concs[n_unique_concs] = cycle_conc[i];
            seen_count[n_unique_concs] = 1;
            n_unique_concs++;
        }
        if (rep_num > 1) {
            fprintf(stderr, "  Cycle %d (%.2f nM): replicate #%d\n",
                    i + 1, cycle_conc[i], rep_num);
        }
    }
    data->ncycles = ncycles;

    fprintf(stderr, "  %d experimental cycles (%d unique concentrations) from %d columns\n",
            ncycles, n_unique_concs, nh);

    /* Read data rows */
    while (fgets(line, sizeof(line), f)) {
        int len = (int)strlen(line);
        while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r'))
            line[--len] = '\0';
        if (len == 0) continue;

        char *fields[512];
        int nf = 0;
        char *p = line;
        while (nf < 512) {
            char *tab = strchr(p, '\t');
            if (tab) { *tab = '\0'; fields[nf++] = p; p = tab + 1; }
            else     { fields[nf++] = p; break; }
        }

        for (int i = 0; i < ncycles; i++) {
            int xc = data_col[i];
            int yc = xc + 1;
            if (yc >= nf) continue;
            if (fields[xc][0] == '\0') continue;

            int n = data->cycles[i].npoints;
            if (n >= MAX_POINTS) continue;

            data->cycles[i].time[n] = atof(fields[xc]);
            data->cycles[i].resp[n] = atof(fields[yc]);

            /* Instrument-fitted response (from Fitted_Y column if present) */
            data->cycles[i].inst_resp[n] = 0.0;
            if (fit_col[i] >= 0 && fit_col[i] + 1 < nf &&
                fields[fit_col[i] + 1][0] != '\0') {
                data->cycles[i].inst_resp[n] = atof(fields[fit_col[i] + 1]);
            }

            data->cycles[i].npoints++;
        }
    }
    fclose(f);
    free(hcopy);

    return 0;
}

/* ================================================================ */
/*  MCK Timing Auto-Detection                                        */
/* ================================================================ */

/* Detect association end time from zero-gap in the highest-conc cycle.
 * Biacore exports have a characteristic zero-region between association
 * and dissociation phases. We find the last non-zero time before the
 * gap in the cycle with the highest concentration. */
double detect_mck_assoc_end(const SPRData *data) {
    /* Find cycle with highest concentration */
    int best = 0;
    double max_conc = data->cycles[0].conc_nM;
    for (int i = 1; i < data->ncycles; i++) {
        if (data->cycles[i].conc_nM > max_conc) {
            max_conc = data->cycles[i].conc_nM;
            best = i;
        }
    }

    const Cycle *cy = &data->cycles[best];

    /* Scan for zero-gap: find the first zero-region after t > 0 that is
     * followed by non-zero data (i.e., the gap between assoc and dissoc).
     * The association end is the time of the last non-zero point before the gap.
     * Note: called before preprocessing so resp values are raw from file. */
    for (int j = 1; j < cy->npoints - 1; j++) {
        if (cy->time[j] <= 0.0) continue;
        /* Look for transition: non-zero → zero */
        if (cy->resp[j] != 0.0 && cy->resp[j+1] == 0.0) {
            /* Verify there's non-zero data after the gap */
            for (int k = j + 2; k < cy->npoints; k++) {
                if (cy->resp[k] != 0.0) {
                    double t_end = cy->time[j];
                    fprintf(stderr, "Auto-detected association end: %.1f s "
                            "(from cycle %d, %.1f nM)\n",
                            t_end, best + 1, max_conc);
                    return t_end;
                }
            }
        }
    }

    /* Fallback: use default 60s */
    fprintf(stderr, "Could not auto-detect association end, using default 60.0 s\n");
    return 60.0;
}

/* ================================================================ */
/*  SCK Format Reader                                                */
/* ================================================================ */

int read_sck_data(const char *filename, SPRData *data, const FitConfig *cfg) {
    FILE *f = fopen(filename, "r");
    if (!f) { fprintf(stderr, "Error: cannot open %s\n", filename); return -1; }

    strncpy(data->filename, filename, 511);
    data->filename[511] = '\0';
    data->mode = MODE_SCK;

    char line[131072];
    if (!fgets(line, sizeof(line), f)) { fclose(f); return -1; }

    /* Parse concentrations from header */
    data->n_injections = parse_sck_concs(line, data->inj_conc_nM, MAX_INJECTIONS);
    if (data->n_injections == 0) {
        fprintf(stderr, "Error: no concentrations found in SCK header\n");
        fclose(f);
        return -1;
    }

    fprintf(stderr, "SCK mode: %d injections, concs:", data->n_injections);
    for (int i = 0; i < data->n_injections; i++)
        fprintf(stderr, " %.1f", data->inj_conc_nM[i]);
    fprintf(stderr, " nM\n");

    /* Read the single trace as cycle 0 */
    data->ncycles = 1;
    Cycle *cy = &data->cycles[0];
    cy->npoints = 0;
    cy->conc_nM = 0;  /* Not meaningful for SCK; injection concs in inj_conc_nM */

    while (fgets(line, sizeof(line), f)) {
        int len = (int)strlen(line);
        while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r')) line[--len] = '\0';
        if (len == 0) continue;
        if (cy->npoints >= MAX_POINTS) break;

        double t = 0, r = 0, ft = 0, fr = 0;
        int nf = sscanf(line, "%lf %lf %lf %lf", &t, &r, &ft, &fr);
        if (nf >= 2) {
            cy->time[cy->npoints] = t;
            cy->resp[cy->npoints] = r;
            cy->inst_resp[cy->npoints] = (nf >= 4) ? fr : 0.0;
            cy->skip[cy->npoints] = 0;
            cy->npoints++;
        }
    }
    fclose(f);

    /* ── Auto-detect injection boundaries ──
     *
     * Biacore SCK exports have two artifact patterns at injection transitions:
     *
     * Pattern A (time gaps): Missing data points create time gaps > 1s.
     *   Each injection transition has one gap at start and one at end.
     *
     * Pattern B (response artifacts): No time gaps, but response values
     *   are replaced with the time value at boundaries (resp ≈ time).
     *   These show as huge response jumps (dr > 50 RU in one step).
     *
     * In both cases, boundaries come in pairs for each injection:
     *   boundary[2k]   = start of injection k
     *   boundary[2k+1] = end of injection k
     */

    /* Compute median dt */
    double median_dt = 0.1;
    if (cy->npoints > 10) {
        double *dts = (double*)malloc(sizeof(double) * (cy->npoints - 1));
        for (int j = 0; j < cy->npoints - 1; j++)
            dts[j] = cy->time[j+1] - cy->time[j];
        /* Simple selection sort for median (small array in practice) */
        int ndt = cy->npoints - 1;
        for (int i = 0; i < ndt/2 + 1; i++) {
            int mi = i;
            for (int j = i+1; j < ndt; j++)
                if (dts[j] < dts[mi]) mi = j;
            double tmp = dts[i]; dts[i] = dts[mi]; dts[mi] = tmp;
        }
        median_dt = dts[ndt/2];
        free(dts);
    }

    /* Detect boundaries — try time gaps first, then response artifacts */
    double gap_before[64], gap_after[64]; /* boundary pairs */
    int n_boundaries = 0;

    /* Method A: Time gap detection (dt > 4 * median_dt) */
    for (int j = 0; j < cy->npoints - 1 && n_boundaries < 64; j++) {
        double dt = cy->time[j+1] - cy->time[j];
        if (dt > median_dt * 4.0 && dt > 0.5) {
            gap_before[n_boundaries] = cy->time[j];
            gap_after[n_boundaries] = cy->time[j+1];
            n_boundaries++;
        }
    }

    /* Method B: Response artifact detection (if no time gaps found)
     * Biacore artifacts: resp ≈ time at boundaries, causing huge jumps */
    if (n_boundaries == 0) {
        /* First pass: mark artifact points (resp ≈ time for t > 5) */
        int *is_artifact = (int*)calloc(cy->npoints, sizeof(int));

        /* Mark candidate artifact points: resp ≈ time (Biacore-specific) */
        for (int j = 0; j < cy->npoints; j++) {
            if (cy->time[j] > 5.0 &&
                fabs(cy->resp[j] - cy->time[j]) < 1.0) {
                is_artifact[j] = 1;
            }
        }
        /* Also mark t=0 exact zeros */
        for (int j = 0; j < cy->npoints; j++) {
            if (fabs(cy->time[j]) < 0.01 && cy->resp[j] == 0.0)
                is_artifact[j] = 1;
        }

        /* Filter: only keep runs of >= 8 consecutive candidates.
         * This rejects false positives where response coincidentally
         * equals time (e.g., 750 RU response at t=750s). Real Biacore
         * artifacts always produce long runs (15-25 points). */
        {
            int run_start = -1;
            for (int j = 0; j <= cy->npoints; j++) {
                int art = (j < cy->npoints) ? is_artifact[j] : 0;
                if (art && run_start < 0) {
                    run_start = j;
                } else if (!art && run_start >= 0) {
                    int run_len = j - run_start;
                    if (run_len < 8) {
                        /* Too short — clear these candidates */
                        for (int k = run_start; k < j; k++)
                            is_artifact[k] = 0;
                    }
                    run_start = -1;
                }
            }
        }

        /* Group consecutive artifact points into boundary regions */
        int in_region = 0;
        double region_start = 0;
        for (int j = 0; j < cy->npoints; j++) {
            if (is_artifact[j] && !in_region) {
                in_region = 1;
                region_start = (j > 0) ? cy->time[j-1] : cy->time[j];
            } else if (!is_artifact[j] && in_region) {
                if (n_boundaries < 64) {
                    gap_before[n_boundaries] = region_start;
                    gap_after[n_boundaries] = cy->time[j];
                    n_boundaries++;
                }
                in_region = 0;
            }
        }

        /* Mark artifact points as skip */
        for (int j = 0; j < cy->npoints; j++)
            if (is_artifact[j]) cy->skip[j] = 1;

        /* If the first boundary is well after t=0, the t=0 single-zero
         * marker was filtered out. Prepend a synthetic boundary so that
         * the first injection starts at t≈0. */
        if (n_boundaries > 0 && gap_after[0] > 10.0) {
            /* Find first point at or after t=0 */
            double t0_after = 0.1;
            for (int j = 0; j < cy->npoints; j++) {
                if (cy->time[j] >= 0.0 && !is_artifact[j]) {
                    t0_after = cy->time[j];
                    break;
                }
            }
            /* Shift all boundaries right by 1 */
            for (int b = n_boundaries; b > 0; b--) {
                gap_before[b] = gap_before[b-1];
                gap_after[b] = gap_after[b-1];
            }
            /* Insert synthetic boundary at t=0 */
            gap_before[0] = -0.5;
            gap_after[0] = t0_after;
            n_boundaries++;
        }

        free(is_artifact);
    } else {
        /* For time-gap data, mark points in gap regions as skip
         * (there are no actual data points in the gaps, but mark the
         *  guard band around gap edges) */
    }

    /* Guard band: skip points within 2s of any boundary edge */
    for (int b = 0; b < n_boundaries; b++) {
        for (int j = 0; j < cy->npoints; j++) {
            if (cy->time[j] >= gap_before[b] - 0.5 &&
                cy->time[j] <= gap_after[b] + 0.5)
                cy->skip[j] = 1;
        }
    }

    /* ── Derive injection schedule from boundaries ──
     * Boundaries come in pairs: even = injection start, odd = injection end.
     * If first boundary is very close to t=0, it's the pre-injection marker.
     */
    int n_inj_detected = n_boundaries / 2;

    if (n_inj_detected == data->n_injections) {
        /* Perfect match: boundary[2k] and boundary[2k+1] bracket injection k */
        for (int i = 0; i < data->n_injections; i++) {
            data->inj_start[i] = gap_after[2*i];      /* after the start gap */
            data->inj_end[i]   = gap_before[2*i + 1];  /* before the end gap */
        }
    } else if (n_boundaries > 0) {
        /* Heuristic: try to match boundaries to injections */
        /* If we have 2*N boundaries but N != n_injections, log a warning
         * and use the boundaries we found */
        fprintf(stderr, "  Warning: detected %d boundaries for %d injections\n",
                n_boundaries, data->n_injections);
        int n_use = (n_boundaries / 2 < data->n_injections) ?
                    n_boundaries / 2 : data->n_injections;
        for (int i = 0; i < n_use; i++) {
            data->inj_start[i] = gap_after[2*i];
            data->inj_end[i]   = gap_before[2*i + 1];
        }
        /* If we detected more injections than header says, trim */
        if (n_inj_detected < data->n_injections)
            data->n_injections = n_inj_detected;
    } else {
        /* Fallback: no boundaries detected, use equal-duration schedule */
        fprintf(stderr, "  Warning: no injection boundaries detected, using default schedule\n");
        double total_time = cy->time[cy->npoints-1];
        double cycle_dur = total_time / data->n_injections;
        double assoc_dur = cycle_dur * 0.67;
        for (int i = 0; i < data->n_injections; i++) {
            data->inj_start[i] = i * cycle_dur;
            data->inj_end[i] = i * cycle_dur + assoc_dur;
        }
    }

    /* Log detected injection schedule */
    for (int i = 0; i < data->n_injections; i++) {
        fprintf(stderr, "  Injection %d: %.1f nM, t=[%.1f, %.1f] (%.1fs)\n",
                i+1, data->inj_conc_nM[i],
                data->inj_start[i], data->inj_end[i],
                data->inj_end[i] - data->inj_start[i]);
    }

    /* ── Baseline subtraction ── */
    double sum = 0.0;
    int count = 0;
    for (int j = 0; j < cy->npoints; j++) {
        if (cy->time[j] < 0.0 && !cy->skip[j]) {
            sum += cy->resp[j];
            count++;
        }
    }
    cy->baseline = (count > 0) ? sum / count : 0.0;
    for (int j = 0; j < cy->npoints; j++)
        cy->resp[j] -= cy->baseline;

    /* Find fit_start */
    cy->fit_start = cy->npoints;
    for (int j = 0; j < cy->npoints; j++) {
        if (cy->time[j] >= 0.0) { cy->fit_start = j; break; }
    }

    /* Count skipped points */
    int n_skip = 0;
    for (int j = cy->fit_start; j < cy->npoints; j++)
        if (cy->skip[j]) n_skip++;

    fprintf(stderr, "  SCK trace: %d points, t=[%.1f, %.1f], %d injections, %d artifacts skipped\n",
            cy->npoints, cy->time[0], cy->time[cy->npoints-1],
            data->n_injections, n_skip);

    return 0;
}

/* ================================================================ */
/*  Auto-detect Format                                               */
/* ================================================================ */

int read_spr_data(const char *filename, SPRData *data, const FitConfig *cfg) {
    /* Peek at header to determine format */
    FILE *f = fopen(filename, "r");
    if (!f) { fprintf(stderr, "Error: cannot open %s\n", filename); return -1; }

    char line[131072];
    if (!fgets(line, sizeof(line), f)) { fclose(f); return -1; }
    fclose(f);

    /* SCK headers list multiple concentrations: "Conc 8 40 200 1000 5000 nM".
     * MCK headers have one concentration per column: "Conc 100 nM".
     * Check for multi-concentration header first — this is the definitive
     * SCK indicator, even if the file has many tab-columns (e.g., SCK
     * with replicates exported as separate column groups). */
    double test_concs[MAX_INJECTIONS];
    int n_concs = parse_sck_concs(line, test_concs, MAX_INJECTIONS);

    if (n_concs > 1) {
        fprintf(stderr, "Auto-detected SCK format (multi-concentration header)\n");
        return read_sck_data(filename, data, cfg);
    }

    /* Count tabs to distinguish single-conc formats:
     * MCK: many columns (4 per cycle, so >= 8 tabs for 2+ cycles)
     * SCK: 4 columns (3 tabs) */
    int ntabs = 0;
    for (char *p = line; *p; p++)
        if (*p == '\t') ntabs++;

    if (ntabs <= 4) {
        fprintf(stderr, "Auto-detected SCK format\n");
        return read_sck_data(filename, data, cfg);
    } else {
        fprintf(stderr, "Auto-detected MCK format\n");
        return read_mck_data(filename, data, cfg);
    }
}

/* ================================================================ */
/*  Reference Subtraction                                            */
/* ================================================================ */

int subtract_reference(SPRData *data, const char *ref_filename, const FitConfig *cfg) {
    static SPRData ref;
    memset(&ref, 0, sizeof(ref));

    if (read_mck_data(ref_filename, &ref, cfg) != 0) return -1;

    if (ref.ncycles != data->ncycles) {
        fprintf(stderr, "Error: reference has %d cycles, data has %d\n",
                ref.ncycles, data->ncycles);
        return -1;
    }

    fprintf(stderr, "Subtracting reference: %s\n", ref_filename);

    for (int c = 0; c < data->ncycles; c++) {
        Cycle *dc = &data->cycles[c];
        Cycle *rc = &ref.cycles[c];

        if (dc->npoints != rc->npoints) {
            fprintf(stderr, "  Warning: cycle %d point count mismatch (%d vs %d)\n",
                    c+1, dc->npoints, rc->npoints);
        }

        int n = dc->npoints < rc->npoints ? dc->npoints : rc->npoints;
        for (int j = 0; j < n; j++) {
            dc->resp[j] -= rc->resp[j];
            if (rc->skip[j]) dc->skip[j] = 1;
        }
    }

    return 0;
}

/* ================================================================ */
/*  JSON Output Helpers                                              */
/* ================================================================ */

void json_print_double_array(const double *arr, int n, int skip, const char *fmt) {
    printf("[");
    int first = 1;
    for (int i = 0; i < n; i += skip) {
        if (!first) printf(",");
        printf(fmt, arr[i]);
        first = 0;
    }
    printf("]");
}

void json_print_int_array(const int *arr, int n, int skip) {
    printf("[");
    int first = 1;
    for (int i = 0; i < n; i += skip) {
        if (!first) printf(",");
        printf("%d", arr[i]);
        first = 0;
    }
    printf("]");
}
