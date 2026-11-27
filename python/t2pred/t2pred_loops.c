#include <tempo2pred.h>
#include "t2pred_loops.h"
#include <math.h>

void T2Predictor_GetPhase_array_ld(const T2Predictor *t2p, long double *mjd, const int nmjd, long double freq, double* out) {
    int i;
    for (i=0; i < nmjd; ++i){
        out[i] = T2Predictor_GetPhase(t2p,mjd[i],freq);
    }
}
void T2Predictor_GetFrequency_array_ld(const T2Predictor *t2p, long double *mjd, const int nmjd, long double freq, double* out) {
    int i;
    for (i=0; i < nmjd; ++i){
        out[i] = T2Predictor_GetFrequency(t2p,mjd[i],freq);
    }
}

double get_FractionalPhase(const T2Predictor *t2p, int intmjd, double fracmjd,long double freq) {
#if defined(T2PREDICTOR_API_VERSION) && T2PREDICTOR_API_VERSION >= 2
    T2PredictorMJD mjd = { intmjd, fracmjd };
    T2PredictorPhase phase;
    T2PredictorStatus status = T2Predictor_GetPhaseSplit(t2p, mjd, freq, &phase);
    if (status != T2PREDICTOR_OK) {
        fprintf(stderr, "get_FractionalPhase: predictor error %d\n", (int)status);
        return NAN;
    }
    return (double)phase.fraction;
#else
    long double mjd = (long double) intmjd + (long double) fracmjd;
    long double phase = T2Predictor_GetPhase(t2p,mjd,freq);
    long double floorphase = floorl(phase);
    return (double)(phase-floorphase);
#endif
}





void T2Predictor_GetPhase_array(const T2Predictor *t2p, double refmjd, double *mjd, const int nmjd, double freq, double* out) {
    int i;
    for (i=0; i < nmjd; ++i){
        out[i] = T2Predictor_GetPhase(t2p,mjd[i],freq);
    }
}
void T2Predictor_GetFrequency_array(const T2Predictor *t2p, double refmjd, double *mjd, const int nmjd, double freq, double* out) {
    int i;
    for (i=0; i < nmjd; ++i){
        out[i] = T2Predictor_GetFrequency(t2p,mjd[i],freq);
    }
}


