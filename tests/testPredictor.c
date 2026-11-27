#define _POSIX_C_SOURCE 200809L
#include "tempo2pred.h"
#include "../python/t2pred/t2pred_loops.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

static const char cheby_text[] =
    "ChebyModelSet 2 segments\n"
    "ChebyModel BEGIN\n"
    "PSRNAME J1234\nSITENAME test\n"
    "TIME_RANGE 60000 60002\nFREQ_RANGE 1000 2000\n"
    "DISPERSION_CONSTANT 0\nNCOEFF_TIME 2\nNCOEFF_FREQ 1\n"
    "COEFFS 4000000000\nCOEFFS 172800\nChebyModel END\n"
    "ChebyModel BEGIN\n"
    "PSRNAME J1234\nSITENAME test\n"
    "TIME_RANGE 60001 60003\nFREQ_RANGE 1000 2000\n"
    "DISPERSION_CONSTANT 0\nNCOEFF_TIME 2\nNCOEFF_FREQ 1\n"
    "COEFFS -8\nCOEFFS 172800\nChebyModel END\n";

static const char t1_text[] =
    "TEMPO2: POLYCO TEMPO1 emulation\n"
    "J1234\n20200101\n1200\n60001.00000000000001\n0\n0\n0\n"
    "1000000000.125000000123456\n1\nsite\n1440\n2\n1400\n0\n0\n0\n0\n";

static const char overflow_text[] =
    "ChebyModelSet 1 segments\nChebyModel BEGIN\n"
    "PSRNAME J1234\nSITENAME test\nTIME_RANGE 60000 60002\n"
    "FREQ_RANGE 1000 2000\nDISPERSION_CONSTANT 0\n"
    "NCOEFF_TIME 2\nNCOEFF_FREQ 1\n"
    "COEFFS 40000000000000000000\nCOEFFS 0\nChebyModel END\n";

static FILE *text_file(const char *text)
{
    FILE *f = tmpfile();
    assert(f);
    assert(fputs(text, f) >= 0);
    rewind(f);
    return f;
}

static void replace_predictor(T2Predictor *p, const char *text,
        T2PredictorKind kind)
{
    FILE *f = text_file(text);
    assert(T2Predictor_FRead(p, f) == 0);
    assert(fclose(f) == 0);
    assert(T2Predictor_Kind(p) == kind);
}

static void assert_invalid_ranges(const char *time_range,
        const char *frequency_range, const char *dispersion)
{
    char text[512];
    T2Predictor p;
    FILE *f;
    int length = snprintf(text, sizeof(text),
            "ChebyModelSet 1 segments\nChebyModel BEGIN\n"
            "PSRNAME J1234\nSITENAME test\n"
            "TIME_RANGE %s\nFREQ_RANGE %s\nDISPERSION_CONSTANT %s\n"
            "NCOEFF_TIME 2\nNCOEFF_FREQ 1\nCOEFFS 1\nCOEFFS 0\n"
            "ChebyModel END\n",
            time_range, frequency_range, dispersion);
    assert(length > 0 && (size_t)length < sizeof(text));
    T2Predictor_Init(&p);
    f = text_file(text);
    assert(T2Predictor_FRead(&p, f) == 0);
    assert(fclose(f) == 0);
    assert(T2Predictor_ValidateSegments(&p) == T2PREDICTOR_INVALID);
    T2Predictor_Destroy(&p);
}

static void assert_invalid_polyco_numbers(const char *epoch,
        const char *reference_phase, const char *frequency,
        const char *coefficient)
{
    char text[512];
    T2Predictor p;
    FILE *f;
    int length = snprintf(text, sizeof(text),
            "TEMPO2: POLYCO TEMPO1 emulation\n"
            "J1234\n20200101\n1200\n%s\n0\n0\n0\n%s\n%s\n"
            "site\n1440\n2\n1400\n0\n0\n%s\n0\n",
            epoch, reference_phase, frequency, coefficient);
    assert(length > 0 && (size_t)length < sizeof(text));
    T2Predictor_Init(&p);
    f = text_file(text);
    assert(T2Predictor_FRead(&p, f) != 0);
    assert(fclose(f) == 0);
    assert(T2Predictor_Kind(&p) == NonePredType);
    T2Predictor_Destroy(&p);
}

int main(void)
{
    T2Predictor p, copy, reread;
    T2PredictorMJD epoch = {60001, 1e-14L};
    T2PredictorPhase phase, unchanged = {123, 0.25L};
    T2PredictorSegmentInfo info;
    long double freq;
    unsigned n;
    FILE *f;
    char filename[] = "testPredictor-XXXXXX";
    int fd, nplan;
    long double plan[2];
    long double legacy_epoch = 60002.5L, legacy_phase;
    T2PredictorMJD start = {60001, 0}, end = {60001, 1e-5L};

    T2Predictor_Init(&p);
    f = text_file(cheby_text);
    assert(T2Predictor_FRead(&p, f) == 0);
    fclose(f);
    assert(T2Predictor_GetSegmentCount(&p, &n) == T2PREDICTOR_OK && n == 2);
    assert(T2Predictor_GetPhase(&p, 60001.0L, 1400) == 1000000000.0L);
    assert(T2Predictor_GetFrequency(&p, 60001.0L, 1400) == 1.0L);
    assert(T2Predictor_GetSegmentInfo(&p, 0, &info) == T2PREDICTOR_OK);
    assert(info.time_coefficients == 2 && info.frequency_coefficients == 1);
    assert(info.start.day == 60000 && info.end.day == 60002);
    assert(info.midpoint.day == 60001 && info.length_days == 2);
    assert(T2Predictor_GetStartMJDSplit(&p, &epoch) == T2PREDICTOR_OK && epoch.day == 60000);
    assert(T2Predictor_GetEndMJDSplit(&p, &epoch) == T2PREDICTOR_OK && epoch.day == 60003);
    assert(T2Predictor_ValidateSegments(&p) == T2PREDICTOR_OK);
    assert_invalid_ranges("60001 60001", "1000 2000", "0");
    assert_invalid_ranges("60002 60000", "1000 2000", "0");
    assert_invalid_ranges("60000 60002", "1000 1000", "0");
    assert_invalid_ranges("60000 60002", "2000 1000", "0");
    assert_invalid_ranges("nan 60002", "1000 2000", "0");
    assert_invalid_ranges("60000 60002", "1000 inf", "0");
    assert_invalid_ranges("60000 60002", "1000 2000", "nan");
    epoch.day = 60001;
    epoch.fraction = 1e-14L;
    assert(T2Predictor_GetPhaseSplit(&p, epoch, 1400, &phase) == T2PREDICTOR_OK);
    assert(phase.turns == 1000000000);
#if defined(TEST_PREDICTOR_PRECISE)
    assert(fabsl(phase.fraction - 8.64e-10L) < 1e-10L);
    assert(fabs(get_FractionalPhase(&p, 60001, 1e-14, 1400) - 8.64e-10) < 1e-10);
#endif
    assert(T2Predictor_GetFrequencySplit(&p, epoch, 1400, &freq) == T2PREDICTOR_OK);
    assert(fabsl(freq - 1) < 1e-12L);
#if defined(TEST_PREDICTOR_PRECISE)
    epoch.fraction = 1.5e-14L;
    assert(T2Predictor_GetPhaseSplit(&p, epoch, 1400, &phase) == T2PREDICTOR_OK);
    assert(fabsl(phase.fraction - 1.296e-9L) < 1e-10L);
#endif
    epoch.day = 60002;
    epoch.fraction = -1;
    assert(T2Predictor_GetPhaseSplit(&p, epoch, 1400, &phase) == T2PREDICTOR_OK);
    assert(phase.turns == 1000000000 && phase.fraction == 0);
    epoch.day = 60002;
    epoch.fraction = 0;
    assert(T2Predictor_GetSegmentPhaseSplit(&p, 1, epoch, 1400, &phase) == T2PREDICTOR_OK);
    assert(phase.turns == -2 && phase.fraction == 0);
    epoch.fraction = 1e-5L;
    assert(T2Predictor_GetSegmentPhaseSplit(&p, 1, epoch, 1400, &phase) == T2PREDICTOR_OK);
    assert(phase.turns == -2 && fabsl(phase.fraction - 0.864L) < 1e-6L);
    epoch.day = 60001;
    epoch.fraction = 0.5L;
    assert(T2Predictor_GetPhaseSplit(&p, epoch, 1400, &phase) == T2PREDICTOR_OK);
    assert(phase.turns == 1000043200);
    assert(T2Predictor_GetSegmentPhaseSplit(&p, 1, epoch, 1400, &phase) == T2PREDICTOR_OK);
    assert(phase.turns == -43202);
    epoch.day = 59999;
    phase = unchanged;
    assert(T2Predictor_GetPhaseSplit(&p, epoch, 1400, &phase) == T2PREDICTOR_OUT_OF_RANGE);
    assert(ChebyModelSet_OutOfRange && phase.turns == unchanged.turns &&
            phase.fraction == unchanged.fraction);
    assert(T2Predictor_GetSegmentInfo(&p, 2, &info) == T2PREDICTOR_OUT_OF_RANGE);
    epoch.fraction = NAN;
    assert(T2Predictor_GetPhaseSplit(&p, epoch, 1400, &phase) == T2PREDICTOR_INVALID);

    fd = mkstemp(filename);
    assert(fd >= 0);
    f = fdopen(fd, "w");
    assert(f);
    T2Predictor_FWrite(&p, f);
    assert(fclose(f) == 0);
    assert(T2Predictor_GetPlanSplit(filename, start, end, 1, 1400,
                &phase, &nplan, plan, 2) == T2PREDICTOR_OK);
    assert(nplan == 1 && phase.turns == 1000000000);
#if defined(TEST_PREDICTOR_PRECISE)
    assert(fabsl(plan[0] - 1) < 1e-8L);
#else
    assert(fabsl(plan[0] - 1) < 1e-4L);
#endif
    assert(T2Predictor_GetPlan(filename, 60001.0L, 60001.00001L,
                1.0L, 1400, &legacy_phase, &nplan, plan) == 0);
    assert(nplan == 1 && legacy_phase == 1000000000.0L);
    f = fopen(filename, "w");
    assert(f && fputs("ChebyModelSet 0 segments\n", f) >= 0);
    assert(fclose(f) == 0);
    phase = unchanged = (T2PredictorPhase){123, 0.25L};
    nplan = 42;
    assert(T2Predictor_GetPlanSplit(filename, start, end, 1, 1400,
                &phase, &nplan, plan, 2) == T2PREDICTOR_IO_ERROR);
    assert(phase.turns == unchanged.turns && phase.fraction == unchanged.fraction && nplan == 42);
    assert(unlink(filename) == 0);

    T2Predictor_Init(&copy);
    T2Predictor_Copy(&copy, &p);
    assert(T2Predictor_Insert(&copy, &p) == 0);
    assert(T2Predictor_GetSegmentCount(&copy, &n) == T2PREDICTOR_OK && n == 4);
    T2Predictor_Keep(&copy, 1, &legacy_epoch);
    assert(T2Predictor_GetSegmentCount(&copy, &n) == T2PREDICTOR_OK && n == 1);
    T2Predictor_Destroy(&copy);
    T2Predictor_Init(&copy);
    T2Predictor_Copy(&copy, &p);
    epoch.day = 60002;
    epoch.fraction = 0.5L;
    assert(T2Predictor_KeepSplit(&copy, 1, &epoch) == T2PREDICTOR_OK);
    assert(T2Predictor_GetSegmentCount(&copy, &n) == T2PREDICTOR_OK && n == 1);
    assert(T2Predictor_GetSegmentInfo(&copy, 0, &info) == T2PREDICTOR_OK);
    assert(info.start.day == 60001);
    T2Predictor_Init(&reread);
    f = tmpfile();
    assert(f);
    T2Predictor_FWrite(&copy, f);
    rewind(f);
    assert(T2Predictor_FRead(&reread, f) == 0);
    fclose(f);
    assert(T2Predictor_GetPhaseSplit(&copy, epoch, 1400, &phase) == T2PREDICTOR_OK);
    assert(T2Predictor_GetPhaseSplit(&reread, epoch, 1400, &unchanged) == T2PREDICTOR_OK);
    assert(phase.turns == unchanged.turns && phase.fraction == unchanged.fraction);
    T2Predictor_Destroy(&reread);
    T2Predictor_Destroy(&copy);
    T2Predictor_Destroy(&p);

    T2Predictor_Init(&p);
    f = text_file(t1_text);
    assert(T2Predictor_FRead(&p, f) == 0);
    fclose(f);
    assert_invalid_polyco_numbers("60001junk", "1000000000.125", "1",
            "0");
    assert_invalid_polyco_numbers("60001", "1000000000.125junk", "1",
            "0");
    assert_invalid_polyco_numbers("60001", "1000000000.125", "1junk",
            "0");
    assert_invalid_polyco_numbers("60001", "1000000000.125", "1",
            "0junk");
    epoch.day = 60001;
    epoch.fraction = 1e-14L;
    assert(T2Predictor_GetPhaseSplit(&p, epoch, 1400, &phase) == T2PREDICTOR_OK);
    assert(phase.turns == 1000000000);
#if defined(TEST_PREDICTOR_PRECISE)
    assert(fabsl(phase.fraction - 0.125000000123456L) < 1e-13L);
#endif
    T2Predictor_Init(&reread);
    f = tmpfile();
    assert(f);
    T2Predictor_FWrite(&p, f);
    rewind(f);
    assert(T2Predictor_FRead(&reread, f) == 0);
    fclose(f);
    assert(T2Predictor_GetPhaseSplit(&reread, epoch, 1400, &unchanged) == T2PREDICTOR_OK);
    assert(phase.turns == unchanged.turns && phase.fraction == unchanged.fraction);
    T2Predictor_Destroy(&reread);
    T2Predictor_Destroy(&p);

    T2Predictor_Init(&p);
    f = tmpfile();
    assert(f);
    fprintf(f, "%-10s %9s%11s%20.11Lf%21.6f %6.3f%7.3f\n",
            "J1234", "20200101", "1200", 60001.0L, 0.0, 0.0, 0.0);
    fprintf(f, "%20.6Lf%18.12Lf%5s%5d%5d%10.3f%16.14f\n",
            10.125L, 1.0L, "site", 1440, 2, 1400.0, 0.0);
    fprintf(f, "%25.17Le%25.17Le\n", 0.0L, 0.0L);
    rewind(f);
    assert(T2Predictor_FRead(&p, f) == 0);
    fclose(f);
    epoch.day = 60001;
    epoch.fraction = 0;
    assert(T2Predictor_GetPhaseSplit(&p, epoch, 1400, &phase) == T2PREDICTOR_OK);
    assert(phase.turns == 10 && phase.fraction == 0.125L);
    f = tmpfile();
    assert(f);
    T2Predictor_FWrite(&p, f);
    rewind(f);
    {
        char line[128];
        assert(fgets(line, sizeof(line), f));
        assert(strncmp(line, "TEMPO2:", 7) != 0);
    }
    rewind(f);
    T2Predictor_Init(&reread);
    assert(T2Predictor_FRead(&reread, f) == 0);
    fclose(f);
    assert(T2Predictor_GetPhaseSplit(&reread, epoch, 1400, &unchanged) == T2PREDICTOR_OK);
    assert(phase.turns == unchanged.turns && phase.fraction == unchanged.fraction);
    T2Predictor_Destroy(&reread);
    T2Predictor_Destroy(&p);

    T2Predictor_Init(&p);
    f = text_file(overflow_text);
    assert(T2Predictor_FRead(&p, f) == 0);
    fclose(f);
    epoch.day = 60001;
    epoch.fraction = 0;
    phase = unchanged = (T2PredictorPhase){123, 0.25L};
    assert(T2Predictor_GetPhaseSplit(&p, epoch, 1400, &phase) == T2PREDICTOR_OVERFLOW);
    assert(phase.turns == unchanged.turns && phase.fraction == unchanged.fraction);
    T2Predictor_Destroy(&p);

    T2Predictor_Init(&p);
    replace_predictor(&p, cheby_text, Cheby);
    replace_predictor(&p, cheby_text, Cheby);
    replace_predictor(&p, t1_text, T1);
    replace_predictor(&p, t1_text, T1);
    replace_predictor(&p, cheby_text, Cheby);
    assert(T2Predictor_GetPhase(&p, 60001.0L, 1400) == 1000000000.0L);
    replace_predictor(&p, t1_text, T1);
    f = text_file("not a predictor\n");
    assert(T2Predictor_FRead(&p, f) != 0);
    assert(fclose(f) == 0);
    assert(T2Predictor_Kind(&p) == NonePredType);
    replace_predictor(&p, cheby_text, Cheby);
    T2Predictor_Destroy(&p);
    return 0;
}
