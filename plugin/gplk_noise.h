/**
 * @file
 * @brief Read-only reconstructed timing-noise quantities for gplk.
 * Values and uncertainties returned here are seconds at each TOA's frequency.
 * Conversion to display units belongs in setPlot/doPlot, not in these helpers.
 * See documentation/GPLK_GIZA_DEVELOPER_GUIDE.md for extension points and caveats.
 * TEMPO2 GPLv3-or-later.
 */
#ifndef TEMPO2_GPLK_NOISE_H
#define TEMPO2_GPLK_NOISE_H
#include "tempo2.h"
#include <cfloat>
#include <cmath>
#include <cstring>

// Keep this range contiguous: setPlot, error scaling and W-key cycling use it.
// The older modes 16/17/19 have different representations; do not alias them.
enum { GPLK_NOISE_RED = 20, GPLK_NOISE_DM, GPLK_NOISE_CHROM,
       GPLK_NOISE_TOTAL, GPLK_NOISE_CLEAN };

static bool gplkFinite(double x) { return x == x && x <= DBL_MAX && x >= -DBL_MAX; }

static int gplkNoiseMode(const char *name)
{
    if (std::strcmp(name, "red") == 0) return GPLK_NOISE_RED;
    if (std::strcmp(name, "dm") == 0) return GPLK_NOISE_DM;
    if (std::strcmp(name, "chromatic") == 0) return GPLK_NOISE_CHROM;
    if (std::strcmp(name, "total") == 0) return GPLK_NOISE_TOTAL;
    if (std::strcmp(name, "clean") == 0) return GPLK_NOISE_CLEAN;
    return 0;
}

static const char *gplkNoiseName(int mode)
{
    switch (mode) {
    case GPLK_NOISE_RED: return "Red-noise delay";
    case GPLK_NOISE_DM: return "DM-noise delay";
    case GPLK_NOISE_CHROM: return "Chromatic-noise delay";
    case GPLK_NOISE_TOTAL: return "Total red + DM + chromatic delay";
    case GPLK_NOISE_CLEAN: return "Postfit minus red + DM + chromatic noise";
    default: return "Postfit residual";
    }
}

static bool gplkNoiseConfigured(const pulsar &p, int mode)
{
    if (mode == GPLK_NOISE_RED) return p.TNRedAmp != 0 && p.TNRedGam != 0;
    if (mode == GPLK_NOISE_DM) return p.TNDMAmp != 0 && p.TNDMGam != 0;
    if (mode == GPLK_NOISE_CHROM)
        return p.TNChromAmp != 0 && p.TNChromGam != 0 && p.TNChromIdx != 0;
    return false;
}

// These fields are already delays in seconds at each TOA's observing frequency.
// Gate reads on the model: absent components are not initialised by the core.
static double gplkNoiseComponent(const pulsar &p, const observation &o, int mode,
                                 bool uncertainty = false)
{
    if (!gplkNoiseConfigured(p, mode)) return 0;
    double value = 0;
    if (mode == GPLK_NOISE_RED) value = uncertainty ? o.TNRedErr : o.TNRedSignal;
    if (mode == GPLK_NOISE_DM) value = uncertainty ? o.TNDMErr : o.TNDMSignal;
    if (mode == GPLK_NOISE_CHROM) value = uncertainty ? o.TNChromErr : o.TNChromSignal;
    return gplkFinite(value) ? (uncertainty ? std::fabs(value) : value) : 0;
}

/**
 * @brief Get one component, their sum, or raw postfit residual minus their sum.
 * @param mode One of the GPLK_NOISE_* modes; callers validate it before dispatch.
 * @return Delay in seconds. No observation, parameter or TNsubtract flag changes.
 */
static double gplkNoiseValue(const pulsar &p, const observation &o, int mode)
{
    if (mode >= GPLK_NOISE_RED && mode <= GPLK_NOISE_CHROM)
        return gplkNoiseComponent(p, o, mode);
    const double total = gplkNoiseComponent(p, o, GPLK_NOISE_RED) +
                         gplkNoiseComponent(p, o, GPLK_NOISE_DM) +
                         gplkNoiseComponent(p, o, GPLK_NOISE_CHROM);
    return mode == GPLK_NOISE_CLEAN ? static_cast<double>(o.residual) - total : total;
}

/** Return a standard error in seconds; Total has none, Clean uses TOA errors only. */
static double gplkNoiseUncertainty(const pulsar &p, const observation &o, int mode)
{
    if (mode == GPLK_NOISE_CLEAN) return o.toaErr * 1e-6;
    // A combined posterior uncertainty needs cross-component covariances,
    // which the per-TOA fields do not contain. Do not invent a quadrature sum.
    return mode == GPLK_NOISE_TOTAL ? 0 : gplkNoiseComponent(p, o, mode, true);
}

/**
 * Require configured components and at least one non-deleted sample with a
 * nonzero finite reconstruction or positive uncertainty. This is availability
 * of stored data, not a statistical test for the presence of timing noise.
 */
static bool gplkNoiseAvailable(const pulsar &p, int mode)
{
    for (int c = GPLK_NOISE_RED; c <= GPLK_NOISE_CHROM; ++c) {
        if (mode <= GPLK_NOISE_CHROM && mode != c) continue;
        if (!gplkNoiseConfigured(p, c)) continue;
        for (int i = 0; i < p.nobs; ++i) {
            if (!p.obsn[i].deleted && (gplkNoiseComponent(p, p.obsn[i], c) != 0 ||
                                      gplkNoiseComponent(p, p.obsn[i], c, true) > 0)) return true;
        }
    }
    return false;
}

#endif
