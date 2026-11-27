#include "tempo2pred_int.h"

#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static int finite_internal(longdouble value)
{
#ifdef LONGDOUBLE_IS_FLOAT128
    return !isnanq(value) && !isinfq(value);
#else
    return isfinite(value);
#endif
}

static T2PredictorStatus join_mjd(T2PredictorMJD input, longdouble *value)
{
    if (!isfinite(input.fraction))
        return T2PREDICTOR_INVALID;
    *value = (longdouble)input.day + (longdouble)input.fraction;
    return finite_internal(*value) ? T2PREDICTOR_OK : T2PREDICTOR_OVERFLOW;
}

static T2PredictorStatus split_value(longdouble value, int64_t *whole,
        long double *fraction)
{
    const longdouble limit = longdouble(9223372036854775808.0);
    longdouble integer;
    long double remainder;
    int64_t result;

    if (!finite_internal(value))
        return T2PREDICTOR_INVALID;
    integer = floorl(value);
    if (integer < -limit || integer >= limit)
        return T2PREDICTOR_OVERFLOW;
    result = (int64_t)integer;
    remainder = (long double)(value - integer);
    if (remainder >= 1.0L) {
        if (result == INT64_MAX)
            return T2PREDICTOR_OVERFLOW;
        ++result;
        remainder = 0;
    }
    *whole = result;
    *fraction = remainder;
    return T2PREDICTOR_OK;
}

static T2PredictorStatus split_mjd(longdouble value, T2PredictorMJD *out)
{
    return split_value(value, &out->day, &out->fraction);
}

static T2PredictorStatus split_phase(longdouble value, T2PredictorPhase *out)
{
    return split_value(value, &out->turns, &out->fraction);
}

static T2PredictorStatus segment(const T2Predictor *p, unsigned index,
        const ChebyModel **out)
{
    if (!p || !out)
        return T2PREDICTOR_INVALID;
    if (p->kind != Cheby)
        return T2PREDICTOR_WRONG_KIND;
    if (p->modelset.cheby.nsegments < 0 ||
            index >= (unsigned)p->modelset.cheby.nsegments)
        return T2PREDICTOR_OUT_OF_RANGE;
    *out = &p->modelset.cheby.segments[index];
    return T2PREDICTOR_OK;
}

static T2PredictorStatus evaluate(const T2Predictor *p, const ChebyModel *cm,
        longdouble mjd, longdouble freq, int phase, longdouble *out)
{
    longdouble value;
    const T1Polyco *t1;

    if (!p || !out || !finite_internal(mjd) || !finite_internal(freq) || freq <= 0)
        return T2PREDICTOR_INVALID;
    if (cm) {
        value = phase ? ChebyModel_GetPhase(cm, mjd, freq)
                      : ChebyModel_GetFrequency(cm, mjd, freq);
    } else if (p->kind == Cheby) {
        cm = ChebyModelSet_GetNearest(&p->modelset.cheby, mjd);
        if (!cm)
            return T2PREDICTOR_OUT_OF_RANGE;
        value = phase ? ChebyModel_GetPhase(cm, mjd, freq)
                      : ChebyModel_GetFrequency(cm, mjd, freq);
    } else if (p->kind == T1) {
        if (p->modelset.t1.nsegments <= 0 || !p->modelset.t1.segments)
            return T2PREDICTOR_OUT_OF_RANGE;
        t1 = T1PolycoSet_GetNearest(&p->modelset.t1, mjd);
        if (fabsl(mjd - t1->mjd_mid) > (longdouble)t1->span / 2880)
            return T2PREDICTOR_OUT_OF_RANGE;
        value = phase ? T1Polyco_GetPhase(t1, mjd, freq)
                      : T1Polyco_GetFrequency(t1, mjd, freq);
    } else {
        return T2PREDICTOR_WRONG_KIND;
    }
    if (!finite_internal(value))
        return T2PREDICTOR_INVALID;
    *out = value;
    return T2PREDICTOR_OK;
}

static T2PredictorStatus evaluate_split(const T2Predictor *p,
        const ChebyModel *cm, T2PredictorMJD epoch, long double freq,
        int phase, longdouble *out)
{
    longdouble mjd;
    T2PredictorStatus status = join_mjd(epoch, &mjd);
    if (status != T2PREDICTOR_OK)
        return status;
    return evaluate(p, cm, mjd, (longdouble)freq, phase, out);
}

static T2PredictorStatus phase_split(const T2Predictor *p, const ChebyModel *cm,
        T2PredictorMJD epoch, long double freq, T2PredictorPhase *out)
{
    longdouble value;
    T2PredictorPhase result;
    T2PredictorStatus status;
    if (!out)
        return T2PREDICTOR_INVALID;
    status = evaluate_split(p, cm, epoch, freq, 1, &value);
    if (status == T2PREDICTOR_OK)
        status = split_phase(value, &result);
    if (status == T2PREDICTOR_OK)
        *out = result;
    return status;
}

static T2PredictorStatus frequency_split(const T2Predictor *p,
        const ChebyModel *cm, T2PredictorMJD epoch, long double freq,
        long double *out)
{
    longdouble value;
    long double result;
    T2PredictorStatus status;
    if (!out)
        return T2PREDICTOR_INVALID;
    status = evaluate_split(p, cm, epoch, freq, 0, &value);
    if (status != T2PREDICTOR_OK)
        return status;
    result = (long double)value;
    if (!isfinite(result))
        return T2PREDICTOR_OVERFLOW;
    *out = result;
    return T2PREDICTOR_OK;
}

T2PredictorStatus T2Predictor_GetPhaseSplit(const T2Predictor *p,
        T2PredictorMJD epoch, long double freq, T2PredictorPhase *out)
{
    return phase_split(p, NULL, epoch, freq, out);
}

T2PredictorStatus T2Predictor_GetFrequencySplit(const T2Predictor *p,
        T2PredictorMJD epoch, long double freq, long double *out)
{
    return frequency_split(p, NULL, epoch, freq, out);
}

static T2PredictorStatus boundary(const T2Predictor *p, int end,
        T2PredictorMJD *out)
{
    T2PredictorMJD result;
    T2PredictorStatus status;
    longdouble value;
    if (!p || !out)
        return T2PREDICTOR_INVALID;
    if (p->kind == Cheby) {
        int n = p->modelset.cheby.nsegments;
        const ChebyModel *cm;
        if (n <= 0 || !p->modelset.cheby.segments)
            return T2PREDICTOR_OUT_OF_RANGE;
        cm = &p->modelset.cheby.segments[end ? n - 1 : 0];
        value = end ? cm->mjd_end : cm->mjd_start;
    } else if (p->kind == T1) {
        int n = p->modelset.t1.nsegments;
        const T1Polyco *t1;
        if (n <= 0 || !p->modelset.t1.segments)
            return T2PREDICTOR_OUT_OF_RANGE;
        t1 = &p->modelset.t1.segments[end ? n - 1 : 0];
        value = t1->mjd_mid + (end ? 1 : -1) * (longdouble)t1->span / 2880;
    } else {
        return T2PREDICTOR_WRONG_KIND;
    }
    status = split_mjd(value, &result);
    if (status == T2PREDICTOR_OK)
        *out = result;
    return status;
}

T2PredictorStatus T2Predictor_GetStartMJDSplit(const T2Predictor *p,
        T2PredictorMJD *out)
{
    return boundary(p, 0, out);
}

T2PredictorStatus T2Predictor_GetEndMJDSplit(const T2Predictor *p,
        T2PredictorMJD *out)
{
    return boundary(p, 1, out);
}

T2PredictorStatus T2Predictor_GetSegmentPhaseSplit(const T2Predictor *p,
        unsigned index, T2PredictorMJD epoch, long double freq,
        T2PredictorPhase *out)
{
    const ChebyModel *cm;
    T2PredictorStatus status = segment(p, index, &cm);
    return status == T2PREDICTOR_OK ? phase_split(p, cm, epoch, freq, out) : status;
}

T2PredictorStatus T2Predictor_GetSegmentFrequencySplit(const T2Predictor *p,
        unsigned index, T2PredictorMJD epoch, long double freq, long double *out)
{
    const ChebyModel *cm;
    T2PredictorStatus status = segment(p, index, &cm);
    return status == T2PREDICTOR_OK ? frequency_split(p, cm, epoch, freq, out) : status;
}

T2PredictorStatus T2Predictor_GetSegmentCount(const T2Predictor *p, unsigned *out)
{
    if (!p || !out)
        return T2PREDICTOR_INVALID;
    if (p->kind != Cheby)
        return T2PREDICTOR_WRONG_KIND;
    if (p->modelset.cheby.nsegments < 0)
        return T2PREDICTOR_INVALID;
    *out = (unsigned)p->modelset.cheby.nsegments;
    return T2PREDICTOR_OK;
}

T2PredictorStatus T2Predictor_GetSegmentInfo(const T2Predictor *p,
        unsigned index, T2PredictorSegmentInfo *out)
{
    const ChebyModel *cm;
    T2PredictorSegmentInfo result;
    T2PredictorStatus status;
    if (!out)
        return T2PREDICTOR_INVALID;
    status = segment(p, index, &cm);
    if (status != T2PREDICTOR_OK)
        return status;
    status = split_mjd(cm->mjd_start, &result.start);
    if (status == T2PREDICTOR_OK)
        status = split_mjd(cm->mjd_end, &result.end);
    if (status == T2PREDICTOR_OK)
        status = split_mjd(cm->mjd_start + (cm->mjd_end - cm->mjd_start) / 2,
                &result.midpoint);
    if (status != T2PREDICTOR_OK)
        return status;
    result.time_coefficients = cm->cheby.nx;
    result.frequency_coefficients = cm->cheby.ny;
    result.length_days = (long double)(cm->mjd_end - cm->mjd_start);
    result.frequency_start = (long double)cm->freq_start;
    result.frequency_end = (long double)cm->freq_end;
    result.dispersion_constant = (long double)cm->dispersion_constant;
    if (!isfinite(result.length_days) || !isfinite(result.frequency_start) ||
            !isfinite(result.frequency_end) || !isfinite(result.dispersion_constant))
        return T2PREDICTOR_INVALID;
    *out = result;
    return T2PREDICTOR_OK;
}

T2PredictorStatus T2Predictor_ValidateSegments(const T2Predictor *p)
{
    unsigned n, i;
    T2PredictorStatus status = T2Predictor_GetSegmentCount(p, &n);
    if (status != T2PREDICTOR_OK)
        return status;
    for (i = 0; i < n; ++i) {
        const ChebyModel *cm = &p->modelset.cheby.segments[i];
        const Cheby2D *polys[2] = { &cm->cheby, &cm->frequency_cheby };
        int j;
        if (!finite_internal(cm->mjd_start) ||
                !finite_internal(cm->mjd_end) ||
                cm->mjd_start >= cm->mjd_end ||
                !finite_internal(cm->freq_start) ||
                !finite_internal(cm->freq_end) ||
                cm->freq_start >= cm->freq_end ||
                !finite_internal(cm->dispersion_constant))
            return T2PREDICTOR_INVALID;
        for (j = 0; j < 2; ++j) {
            size_t k, count;
            if (polys[j]->nx <= 0 || polys[j]->ny <= 0 || !polys[j]->coeff ||
                    (size_t)polys[j]->nx > SIZE_MAX / (size_t)polys[j]->ny)
                return T2PREDICTOR_INVALID;
            count = (size_t)polys[j]->nx * (size_t)polys[j]->ny;
            for (k = 0; k < count; ++k)
                if (!finite_internal(polys[j]->coeff[k]))
                    return T2PREDICTOR_INVALID;
        }
    }
    return T2PREDICTOR_OK;
}

T2PredictorStatus T2Predictor_KeepSplit(T2Predictor *p, unsigned count,
        const T2PredictorMJD *epochs)
{
    longdouble *values;
    unsigned i;
    T2PredictorStatus status = T2PREDICTOR_OK;
    if (!p || (count && !epochs))
        return T2PREDICTOR_INVALID;
    if (p->kind != Cheby)
        return T2PREDICTOR_WRONG_KIND;
    if (count && SIZE_MAX / count < sizeof(*values))
        return T2PREDICTOR_OVERFLOW;
    values = count ? malloc((size_t)count * sizeof(*values)) : NULL;
    if (count && !values)
        return T2PREDICTOR_NOMEM;
    for (i = 0; i < count; ++i) {
        status = join_mjd(epochs[i], &values[i]);
        if (status != T2PREDICTOR_OK)
            break;
    }
    if (status == T2PREDICTOR_OK &&
            ChebyModelSet_Keep(&p->modelset.cheby, count, values) != 0)
        status = T2PREDICTOR_NOMEM;
    free(values);
    return status;
}

T2PredictorStatus T2Predictor_GetPlanSplit(const char *filename,
        T2PredictorMJD start, T2PredictorMJD end, long double step,
        long double freq, T2PredictorPhase *phase0, int *nsegments,
        long double *frequencies, unsigned capacity)
{
    return T2Predictor_GetPlanSplit_Ext(filename, start, end, step, freq,
            NULL, 0, NULL, 0, phase0, nsegments, frequencies, capacity);
}

T2PredictorStatus T2Predictor_GetPlanSplit_Ext(const char *filename,
        T2PredictorMJD start, T2PredictorMJD end, long double step,
        long double freq, char *psrname, size_t psrname_size,
        char *sitename, size_t sitename_size, T2PredictorPhase *phase0,
        int *nsegments, long double *frequencies, unsigned capacity)
{
    T2Predictor pred;
    T2PredictorPhase first;
    T2PredictorStatus status;
    longdouble first_mjd, end_mjd, last_mjd, mjd, last_phase, phase, count;
    long double *values = NULL;
    FILE *f;
    int i, n;
    const char *psr, *site;

    if (!filename || !phase0 || !nsegments || !isfinite(step) || step <= 0 ||
            !isfinite(freq) || freq <= 0 || (capacity && !frequencies) ||
            (psrname && !psrname_size) || (sitename && !sitename_size))
        return T2PREDICTOR_INVALID;
    status = join_mjd(start, &first_mjd);
    if (status == T2PREDICTOR_OK)
        status = join_mjd(end, &end_mjd);
    if (status != T2PREDICTOR_OK)
        return status;
    if (end_mjd <= first_mjd)
        return T2PREDICTOR_INVALID;
    count = -floorl(-((end_mjd - first_mjd) * 86400 / (longdouble)step));
    if (!finite_internal(count) || count > INT_MAX || count > capacity)
        return T2PREDICTOR_OVERFLOW;
    n = (int)count;
    if (!frequencies || n <= 0)
        return T2PREDICTOR_INVALID;
    if (SIZE_MAX / (size_t)n < sizeof(*values))
        return T2PREDICTOR_OVERFLOW;
    f = fopen(filename, "r");
    if (!f)
        return T2PREDICTOR_IO_ERROR;
    T2Predictor_Init(&pred);
    i = T2Predictor_FRead(&pred, f);
    fclose(f);
    if (i != 0)
        return T2PREDICTOR_IO_ERROR;
    psr = T2Predictor_GetPSRName(&pred);
    site = T2Predictor_GetSiteName(&pred);
    if ((psrname && strlen(psr) >= psrname_size) ||
            (sitename && strlen(site) >= sitename_size)) {
        status = T2PREDICTOR_OVERFLOW;
        goto done;
    }
    values = malloc((size_t)n * sizeof(*values));
    if (!values) {
        status = T2PREDICTOR_NOMEM;
        goto done;
    }
    status = evaluate(&pred, NULL, first_mjd, (longdouble)freq, 1, &last_phase);
    if (status == T2PREDICTOR_OK)
        status = split_phase(last_phase, &first);
    last_mjd = first_mjd;
    for (i = 0; status == T2PREDICTOR_OK && i < n; ++i) {
        mjd = first_mjd + (longdouble)(i + 1) * (longdouble)step / 86400;
        if (mjd > end_mjd)
            mjd = end_mjd;
        if (mjd <= last_mjd) {
            status = T2PREDICTOR_INVALID;
            break;
        }
        status = evaluate(&pred, NULL, mjd, (longdouble)freq, 1, &phase);
        if (status == T2PREDICTOR_OK) {
            values[i] = (long double)((phase - last_phase) / (mjd - last_mjd) / 86400);
            if (!isfinite(values[i]))
                status = T2PREDICTOR_OVERFLOW;
            last_phase = phase;
            last_mjd = mjd;
        }
    }
    if (status == T2PREDICTOR_OK) {
        if (psrname)
            strcpy(psrname, psr);
        if (sitename)
            strcpy(sitename, site);
        *phase0 = first;
        *nsegments = n;
        memcpy(frequencies, values, (size_t)n * sizeof(*values));
    }
done:
    free(values);
    T2Predictor_Destroy(&pred);
    return status;
}
