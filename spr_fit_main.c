/*
 * spr_fit_main.c - Unified SPR Fitting Entry Point
 *
 * Supports MCK and SCK data formats, multiple binding models.
 *
 * Usage: ./spr_fit <datafile> [options]
 *   --ref <reffile>           Reference file for double referencing
 *   --model <name>            langmuir | heterogeneous | twostate
 *   --exclude <indices>       Comma-separated cycle indices to exclude (0-based)
 *   --sck                     Force SCK mode
 *   --mck                     Force MCK mode
 *   --restarts <N>            Number of random restarts (default 40)
 *   --assoc-end <T>           Association end time in seconds (default 60)
 *
 * Output: JSON to stdout, log to stderr
 */

#include "spr_types.h"
#include "spr_io.h"
#include "spr_models.h"
#include "spr_optim.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ================================================================ */
/*  CLI Argument Parsing                                             */
/* ================================================================ */

static void print_usage(const char *prog) {
    fprintf(stderr, "Usage: %s <datafile> [options]\n", prog);
    fprintf(stderr, "  --ref <file>        Reference file for subtraction\n");
    fprintf(stderr, "  --model <name>      langmuir | heterogeneous | twostate\n");
    fprintf(stderr, "  --exclude <idx,...>  Cycle indices to exclude (0-based)\n");
    fprintf(stderr, "  --sck / --mck       Force data format\n");
    fprintf(stderr, "  --restarts <N>      Random restarts (default 40)\n");
    fprintf(stderr, "  --assoc-end <T>     Association end time (default 60)\n");
}

static ModelType parse_model(const char *name) {
    if (strcmp(name, "heterogeneous") == 0) return MODEL_HETEROGENEOUS;
    if (strcmp(name, "twostate") == 0) return MODEL_TWOSTATE;
    return MODEL_LANGMUIR;
}

static void parse_exclude(const char *str, SPRData *data) {
    char *copy = strdup(str);
    char *tok = strtok(copy, ",");
    while (tok) {
        int idx = atoi(tok);
        if (idx >= 0 && idx < MAX_CYCLES)
            data->excluded[idx] = 1;
        tok = strtok(NULL, ",");
    }
    free(copy);
}

/* ================================================================ */
/*  Main                                                             */
/* ================================================================ */

int main(int argc, char *argv[]) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }

    const char *datafile = argv[1];
    const char *reffile = NULL;
    int force_mode = -1;  /* -1=auto, 0=MCK, 1=SCK */

    FitConfig cfg = fit_config_default();
    SPRData data;
    memset(&data, 0, sizeof(data));

    /* Parse options */
    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "--ref") == 0 && i+1 < argc) {
            reffile = argv[++i];
        } else if (strcmp(argv[i], "--model") == 0 && i+1 < argc) {
            cfg.model = parse_model(argv[++i]);
        } else if (strcmp(argv[i], "--exclude") == 0 && i+1 < argc) {
            parse_exclude(argv[++i], &data);
        } else if (strcmp(argv[i], "--sck") == 0) {
            force_mode = 1;
        } else if (strcmp(argv[i], "--mck") == 0) {
            force_mode = 0;
        } else if (strcmp(argv[i], "--restarts") == 0 && i+1 < argc) {
            cfg.nrestarts = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--assoc-end") == 0 && i+1 < argc) {
            cfg.t_assoc_end = atof(argv[++i]);
            cfg.t_assoc_end_set = 1;
        } else {
            /* Legacy: second positional arg is reference file */
            if (!reffile) reffile = argv[i];
        }
    }

    /* Read data */
    int ret;
    if (force_mode == 1)
        ret = read_sck_data(datafile, &data, &cfg);
    else if (force_mode == 0)
        ret = read_mck_data(datafile, &data, &cfg);
    else
        ret = read_spr_data(datafile, &data, &cfg);

    if (ret != 0) return 1;

    /* Auto-detect association end for MCK if not explicitly set
     * (must be done BEFORE preprocessing, which modifies resp values) */
    if (data.mode == MODE_MCK && !cfg.t_assoc_end_set) {
        cfg.t_assoc_end = detect_mck_assoc_end(&data);
    }

    /* Preprocess (baseline subtraction, artifact marking) */
    preprocess_data(&data, cfg.guard_time);

    /* Reference subtraction */
    int has_ref = 0;
    if (reffile) {
        if (subtract_reference(&data, reffile, &cfg) != 0) {
            fprintf(stderr, "Warning: reference subtraction failed, proceeding without\n");
        } else {
            has_ref = 1;
        }
    }

    /* Log */
    fprintf(stderr, "File: %s  Mode: %s  Model: %s  Cycles: %d  Concs:",
            datafile,
            data.mode == MODE_MCK ? "MCK" : "SCK",
            MODEL_NAMES[cfg.model],
            data.ncycles);
    for (int i = 0; i < data.ncycles; i++) {
        if (data.excluded[i])
            fprintf(stderr, " [%.2f]", data.cycles[i].conc_nM);
        else
            fprintf(stderr, " %.2f", data.cycles[i].conc_nM);
    }
    fprintf(stderr, " nM\n");

    /* Fit */
    FitResult result;
    if (optim_fit(&data, &cfg, &result) != 0) {
        fprintf(stderr, "Fitting failed\n");
        return 1;
    }

    /* Output JSON (patch reference field) */
    /* We output directly and replace reference flag */
    double *buf = (double*)malloc(sizeof(double) * MAX_POINTS);

    printf("{\n");
    printf("  \"filename\": \"%s\",\n", data.filename);
    printf("  \"mode\": \"%s\",\n", data.mode == MODE_MCK ? "mck" : "sck");
    printf("  \"model\": \"%s\",\n", MODEL_NAMES[result.model]);
    printf("  \"ncycles\": %d,\n", data.ncycles);
    printf("  \"reference\": %s,\n", has_ref ? "true" : "false");
    if (has_ref) printf("  \"ref_file\": \"%s\",\n", reffile);
    printf("  \"t_assoc_end\": %.2f,\n", cfg.t_assoc_end);

    printf("  \"ka\": %.8e,\n", result.ka);
    printf("  \"kd\": %.8e,\n", result.kd);
    printf("  \"KD_M\": %.8e,\n", result.KD_M);
    printf("  \"KD_nM\": %.6f,\n", result.KD_M * 1e9);
    printf("  \"Rmax\": %.6f,\n", result.Rmax);

    if (result.model == MODEL_HETEROGENEOUS) {
        printf("  \"ka2\": %.8e,\n", result.ka2);
        printf("  \"kd2\": %.8e,\n", result.kd2);
        printf("  \"KD2_M\": %.8e,\n", result.KD2_M);
        printf("  \"Rmax2\": %.6f,\n", result.Rmax2);
    }
    if (result.model == MODEL_TWOSTATE) {
        printf("  \"ka2\": %.8e,\n", result.ka2);
        printf("  \"kd2\": %.8e,\n", result.kd2);
    }

    printf("  \"R2\": %.8f,\n", result.R2);
    printf("  \"chi2\": %.8f,\n", result.chi2);
    printf("  \"rms\": %.8f,\n", result.rms);
    printf("  \"ssr\": %.6f,\n", result.ssr);
    printf("  \"n_points\": %d,\n", result.n_points);

    printf("  \"ss_KD_nM\": %.6f,\n", result.ss_KD_nM);
    printf("  \"ss_Rmax\": %.6f,\n", result.ss_Rmax);
    printf("  \"ss_R2\": %.8f,\n", result.ss_R2);
    printf("  \"ss_req\": [");
    for (int i = 0; i < data.ncycles; i++) {
        if (i) printf(", ");
        printf("%.6f", result.ss_req[i]);
    }
    printf("],\n");

    printf("  \"concentrations\": [");
    for (int i = 0; i < data.ncycles; i++) {
        if (i) printf(", ");
        printf("%.4f", data.cycles[i].conc_nM);
    }
    printf("],\n");

    printf("  \"excluded\": [");
    {
        int first = 1;
        for (int i = 0; i < data.ncycles; i++) {
            if (data.excluded[i]) {
                if (!first) printf(", ");
                printf("%d", i);
                first = 0;
            }
        }
    }
    printf("],\n");

    /* Per-cycle data */
    printf("  \"cycles\": [\n");
    for (int c = 0; c < data.ncycles; c++) {
        const Cycle *cy = &data.cycles[c];

        if (data.mode == MODE_MCK)
            simulate_mck_cycle(cfg.model, result.params, cy, buf,
                              cfg.t_assoc_end, cfg.dt);
        else
            simulate_sck_trace(cfg.model, result.params, &data, buf, cfg.dt);

        double cy_ssr = 0;
        int cy_n = 0;
        for (int i = cy->fit_start; i < cy->npoints; i++) {
            if (cy->skip[i]) continue;
            double r = cy->resp[i] - buf[i];
            cy_ssr += r * r;
            cy_n++;
        }

        printf("    {\n");
        printf("      \"conc_nM\": %.4f,\n", cy->conc_nM);
        printf("      \"npoints\": %d,\n", cy->npoints);
        printf("      \"baseline\": %.6f,\n", cy->baseline);
        printf("      \"excluded\": %s,\n", data.excluded[c] ? "true" : "false");
        printf("      \"cycle_ssr\": %.6f,\n", cy_ssr);
        printf("      \"cycle_rms\": %.6f,\n", sqrt(cy_ssr / (cy_n > 0 ? cy_n : 1)));

        printf("      \"time\": ");
        json_print_double_array(cy->time, cy->npoints, DISPLAY_SKIP, "%.3f");
        printf(",\n");

        printf("      \"response\": ");
        json_print_double_array(cy->resp, cy->npoints, DISPLAY_SKIP, "%.4f");
        printf(",\n");

        printf("      \"fitted\": ");
        json_print_double_array(buf, cy->npoints, DISPLAY_SKIP, "%.4f");
        printf(",\n");

        printf("      \"inst_fitted\": ");
        json_print_double_array(cy->inst_resp, cy->npoints, DISPLAY_SKIP, "%.4f");
        printf(",\n");

        printf("      \"skip\": ");
        json_print_int_array(cy->skip, cy->npoints, DISPLAY_SKIP);
        printf("\n");

        printf("    }%s\n", c < data.ncycles - 1 ? "," : "");
    }
    printf("  ]\n");
    printf("}\n");

    free(buf);
    return 0;
}
