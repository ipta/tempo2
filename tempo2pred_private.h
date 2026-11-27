#ifndef TEMPO2PRED_PRIVATE_H
#define TEMPO2PRED_PRIVATE_H

#include "tempo2pred.h"
#include "TKlongdouble.h"

typedef struct
{
    int nx, ny;
    longdouble *coeff;
} Cheby2D;

struct ChebyModel
{
    char psrname[64];
    char sitename[64];
    longdouble mjd_start, mjd_end;
    longdouble freq_start, freq_end;
    longdouble dispersion_constant;
    Cheby2D cheby;
    Cheby2D frequency_cheby;
};

struct T1Polyco
{
    char psrname[64];
    char date_string[10];
    char utc_string[13];
    longdouble mjd_mid;
    double dm;
    double doppler;
    double log10rms;
    longdouble reference_phase;
    longdouble frequency_psr_0;
    char sitename[5];
    int span;
    int ncoeff;
    double frequency_obs;
    double binary_phase;
    double binary_frequency;
    longdouble coeff[32];
    int extended_format;
};

#endif
