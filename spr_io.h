/*
 * spr_io.h - Data I/O for SPR fitting (MCK and SCK formats)
 */

#ifndef SPR_IO_H
#define SPR_IO_H

#include "spr_types.h"

/* Read MCK data file (multi-cycle, tab-separated columns per cycle) */
int read_mck_data(const char *filename, SPRData *data, const FitConfig *cfg);

/* Read SCK data file (single trace, 4 columns: time, resp, fit_time, fit_resp) */
int read_sck_data(const char *filename, SPRData *data, const FitConfig *cfg);

/* Auto-detect format and read */
int read_spr_data(const char *filename, SPRData *data, const FitConfig *cfg);

/* Auto-detect association end time from zero-gap in MCK data (call before preprocess) */
double detect_mck_assoc_end(const SPRData *data);

/* Preprocess all cycles (baseline subtraction, artifact marking) */
void preprocess_data(SPRData *data, double guard_time);

/* Subtract reference channel data (double referencing) */
int subtract_reference(SPRData *data, const char *ref_filename, const FitConfig *cfg);

/* JSON output helpers */
void json_print_double_array(const double *arr, int n, int skip, const char *fmt);
void json_print_int_array(const int *arr, int n, int skip);

#endif /* SPR_IO_H */
