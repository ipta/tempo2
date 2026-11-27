#include "tempo2pred.h"
#include <cassert>

int main(int argc, char **argv)
{
#if !defined(T2PREDICTOR_API_VERSION) || T2PREDICTOR_API_VERSION < 2
#error Predictor split API is unavailable
#endif
    T2Predictor predictor;
    T2Predictor_Init(&predictor);
    if (argc > 1) {
        unsigned count;
        T2PredictorSegmentInfo info;
        T2PredictorPhase phase, reread_phase;
        T2Predictor reread;
        long double frequency;
        assert(T2Predictor_Read(&predictor, argv[1]) == 0);
        assert(T2Predictor_GetSegmentCount(&predictor, &count) == T2PREDICTOR_OK && count);
        assert(T2Predictor_ValidateSegments(&predictor) == T2PREDICTOR_OK);
        assert(T2Predictor_GetSegmentInfo(&predictor, 0, &info) == T2PREDICTOR_OK);
        assert(T2Predictor_GetSegmentPhaseSplit(&predictor, 0, info.midpoint,
                    info.frequency_start, &phase) == T2PREDICTOR_OK);
        assert(T2Predictor_GetSegmentFrequencySplit(&predictor, 0, info.midpoint,
                    info.frequency_start, &frequency) == T2PREDICTOR_OK);
        assert(frequency > 0);
        FILE *file = tmpfile();
        assert(file);
        T2Predictor_FWrite(&predictor, file);
        rewind(file);
        T2Predictor_Init(&reread);
        assert(T2Predictor_FRead(&reread, file) == 0);
        assert(T2Predictor_GetSegmentPhaseSplit(&reread, 0, info.midpoint,
                    info.frequency_start, &reread_phase) == T2PREDICTOR_OK);
        assert(phase.turns == reread_phase.turns && phase.fraction == reread_phase.fraction);
        T2Predictor_Destroy(&reread);
        fclose(file);
    }
    T2Predictor_Destroy(&predictor);
    return 0;
}
