/**
 * @file
 * @brief Explicit overlay controls, solar geometry and fitted-model drawing.
 * Context is owned by graphicalInterface and passed to the viewer/menu.
 * All curves borrow current plotted x/id arrays; input data are never changed.
 */
#ifndef TEMPO2_GPLK_OVERLAYS_H
#define TEMPO2_GPLK_OVERLAYS_H
#include "gplk_giza.h"
#include "gplk_noise.h"
#include <vector>

struct GplkOverlayControls {
    bool solar, shapiro, clearAllRequested;
    int noiseMode; // zero=off; red/DM/chromatic/total GPLK_NOISE_* modes otherwise
    GplkOverlayControls() : solar(false), shapiro(false), clearAllRequested(false), noiseMode(0) {}
};

static bool gplkOverlayTimeAxis(int axis) { return axis==3 || axis==10 || axis==12; }
static bool gplkOverlayDelayAxis(int axis)
{ return axis==1 || axis==2 || (axis>=GPLK_NOISE_RED && axis<=GPLK_NOISE_CLEAN); }

/** Match setPlot's seconds-to-display conversion, including its phase convention. */
static double gplkOverlayDelay(double seconds, double unit, int phase, double spin)
{ return seconds/unit*(phase==1 ? spin : 1.0); }

/**
 * Geometric topocentric Sun/pulsar separation at one TOA, in degrees.
 * Uses the same ephemeris/observer vectors as shapiro_delay.C, reversed to point
 * observer->Sun. atan2(|a cross b|, a dot b) handles conjunction/opposition.
 * This is not an aberrated apparent position or an interpolated conjunction.
 */
static bool gplkSolarSeparation(const observation &o, double &degrees)
{
    if (!o.delayCorr) return false;
    double sun[3], direction[3], sn=0, pn=0;
    const bool spacecraft=std::strcmp(o.telID,"STL_FBAT")==0;
    for (int k=0;k<3;++k) {
        sun[k]=o.sun_ssb[k]-(spacecraft ? 0 : o.earth_ssb[k])-o.observatory_earth[k];
        direction[k]=o.psrPos[k];
        if (!gplkFinite(sun[k]) || !gplkFinite(direction[k])) return false;
        sn+=sun[k]*sun[k]; pn+=direction[k]*direction[k];
    }
    if (!(sn>0) || !(pn>0) || !gplkFinite(sn) || !gplkFinite(pn)) return false;
    sn=std::sqrt(sn); pn=std::sqrt(pn);
    for (int k=0;k<3;++k) { sun[k]/=sn; direction[k]/=pn; }
    const double cross[]={sun[1]*direction[2]-sun[2]*direction[1],
                          sun[2]*direction[0]-sun[0]*direction[2],
                          sun[0]*direction[1]-sun[1]*direction[0]};
    const double dot=sun[0]*direction[0]+sun[1]*direction[1]+sun[2]*direction[2];
    degrees=std::atan2(std::sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]),dot)*180/M_PI;
    return true;
}

static bool gplkHasSolarGeometry(const pulsar &p)
{
    for (int i=0;i<p.nobs;++i) {
        double angle;
        if (!p.obsn[i].deleted && gplkSolarSeparation(p.obsn[i],angle)) return true;
    }
    return false;
}

static void gplkCycleNoiseOverlay(const pulsar &p, GplkOverlayControls &state)
{
    for (int attempt=0;attempt<5;++attempt) {
        state.noiseMode=state.noiseMode==0 ? GPLK_NOISE_RED :
                        state.noiseMode>=GPLK_NOISE_TOTAL ? 0 : state.noiseMode+1;
        if (state.noiseMode==0 || gplkNoiseAvailable(p,state.noiseMode)) return;
    }
    state.noiseMode=0;
}

/** Solar curve shares x positions, but has its own fixed 0..180 degree right axis. */
static bool gplkDrawSolarOverlay(const pulsar &p, int n, const float *x, const int *id,
                                 float xmin, float xmax, float ymin, float ymax)
{
    std::vector<float> sx, sy;
    for (int i=0;i<n;++i) {
        double angle;
        if (p.obsn[id[i]].deleted || !gplkSolarSeparation(p.obsn[id[i]],angle)) continue;
        sx.push_back(x[i]); sy.push_back(ymin+(ymax-ymin)*angle/180);
    }
    if (sx.size()<2) return false;
    giza_save();
    gplkPrepareOverlayColours();
    gplkDrawOverlayLine(sx.size(),&sx[0],&sy[0],GIZA_LS_DASH_DOT,GPLK_SOLAR);
    giza_set_colour_index(GPLK_SOLAR);
    giza_set_line_width(1.0);
    giza_axis_float("NTS",xmax,ymin,xmax,ymax,0,180,30,3,0.4f,0,0.5f,0.3f,0);
    giza_annotate("R",3.4,0.5,0.5,"Sun-pulsar separation (deg)");
    giza_restore();
    return true;
}

static bool gplkDrawShapiroOverlay(const pulsar &p, int n, const float *x, const int *id,
                                   double unit, int phase, double offset)
{
    std::vector<float> sx, sy;
    for (int i=0;i<n;++i) {
        const observation &o=p.obsn[id[i]];
        if (o.deleted || !gplkFinite(o.shapiroDelaySun)) continue;
        sx.push_back(x[i]);
        sy.push_back(gplkOverlayDelay(o.shapiroDelaySun,unit,phase,p.param[param_f].val[0])-offset);
    }
    if (sx.empty()) return false;
    gplkDrawOverlayLine(sx.size(),&sx[0],&sy[0],GIZA_LS_DASHED);
    return true;
}

/** Fitted values at actual observing frequencies, never random noise samples. */
static bool gplkDrawNoiseOverlay(const pulsar &p, int n, const float *x, const int *id,
                                 int mode, double unit, int phase, bool centre)
{
    if (!gplkNoiseAvailable(p,mode)) return false;
    std::vector<float> nx, ny, ne;
    double mean=0;
    for (int i=0;i<n;++i) {
        const observation &o=p.obsn[id[i]];
        if (o.deleted) continue;
        nx.push_back(x[i]);
        ny.push_back(gplkOverlayDelay(gplkNoiseValue(p,o,mode),unit,phase,p.param[param_f].val[0]));
        ne.push_back(std::fabs(gplkOverlayDelay(gplkNoiseUncertainty(p,o,mode),unit,phase,p.param[param_f].val[0])));
        mean+=ny.back();
    }
    if (nx.empty()) return false;
    if (centre) {
        mean/=ny.size();
        for (unsigned int i=0;i<ny.size();++i) ny[i]-=mean;
    }
    giza_save();
    gplkPrepareOverlayColours();
    // Only achromatic red noise is joined. Connecting DM/chromatic values
    // across observing bands would falsely suggest a single continuous delay.
    if (mode==GPLK_NOISE_RED)
        gplkDrawOverlayLine(nx.size(),&nx[0],&ny[0],GIZA_LS_DOT,GPLK_NOISE_OVERLAY);
    giza_set_colour_index(GPLK_NOISE_OVERLAY);
    giza_set_character_height(gplkPlotCharacterHeight(0.9f));
    giza_set_line_width(1.5);
    if (mode!=GPLK_NOISE_TOTAL) // Combined posterior covariance is not available.
        giza_error_bars_float(6,nx.size(),&nx[0],&ny[0],&ne[0],0.7f);
    gplkDrawPoints(nx.size(),&nx[0],&ny[0],6); // hollow squares, unlike residual circles/bins
    giza_restore();
    return true;
}

#endif
