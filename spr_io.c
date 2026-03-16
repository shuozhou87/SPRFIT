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
    cy->conc_nM = 0;  /* Not meaningful for SCK */

    while (fgets(line, sizeof(line), f)) {
        int len = (int)strlen(line);
        while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r')) line[--len] = '\0';
        if (len == 0) continue;

        if (cy->npoints >= MAX_POINTS * MAX_CYCLES) break;  /* Safety */

        double t, r, ft, fr;
        if (sscanf(line, "%lf %lf %lf %lf", &t, &r, &ft, &fr) >= 2) {
            /* Store in cycle 0 (we'll use a larger buffer approach) */
            if (cy->npoints < MAX_POINTS) {
                cy->time[cy->npoints] = t;
                cy->resp[cy->npoints] = r;
                cy->inst_resp[cy->npoints] = fr;
                cy->skip[cy->npoints] = 0;
                cy->npoints++;
            }
        }
    }
    fclose(f);

    /* Compute baseline from t < 0 */
    double sum = 0.0;
    int count = 0;
    for (int j = 0; j < cy->npoints; j++) {
        if (cy->time[j] < 0.0) { sum += cy->resp[j]; count++; }
    }
    cy->baseline = (count > 0) ? sum / count : 0.0;
    for (int j = 0; j < cy->npoints; j++)
        cy->resp[j] -= cy->baseline;

    cy->fit_start = cy->npoints;
    for (int j = 0; j < cy->npoints; j++) {
        if (cy->time[j] >= 0.0) { cy->fit_start = j; break; }
    }

    /* Auto-detect injection boundaries from data:
     * SCK injections are typically equal-duration association phases
     * followed by a final long dissociation. We use the default schedule
     * of 120s association + 60s dissociation per injection, with the last
     * injection having a long dissociation. */
    double assoc_dur = 120.0;   /* Default association duration */
    double short_dissoc = 60.0; /* Short dissociation between injections */

    for (int i = 0; i < data->n_injections; i++) {
        double t_start = i * (assoc_dur + short_dissoc);
        data->inj_start[i] = t_start;
        data->inj_end[i] = t_start + assoc_dur;
    }

    fprintf(stderr, "  SCK trace: %d points, t=[%.1f, %.1f], %d injections\n",
            cy->npoints, cy->time[0], cy->time[cy->npoints-1], data->n_injections);

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

    /* Count tabs to distinguish formats:
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
    SPRData ref;
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
