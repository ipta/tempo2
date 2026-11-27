#include <cstdio>
#include <cstdlib>
#include <tempo2pred.h>

int main (int argc, char** argv){
    T2Predictor pred;

    char* pred_fname = argv[1];
    double mjd=atof(argv[2]);

    T2Predictor_Init(&pred);
    if (T2Predictor_Read(&pred, pred_fname) != 0) {
        fprintf(stderr, "Unable to read predictor %s\n", pred_fname);
        T2Predictor_Destroy(&pred);
        return 1;
    }

    const double psr_freq = (double)T2Predictor_GetFrequency(&pred,mjd,1400);

    printf("%lf %24.20lf %24.20lf\n",mjd,psr_freq,1.0/psr_freq);

    T2Predictor_Destroy(&pred);
    return 0;
}
