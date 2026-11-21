/** TEMPO2 spectral estimators, borrowed data only. GPLv3-or-later.
 * Included after GplkSpectrum/GplkSpectrumSample definitions. Native spectra keep
 * their own Fourier grid and yr^3 normalization; display zoom never resamples.
 */
#ifndef TEMPO2_GPLK_SPECTRUM_NATIVE_H
#define TEMPO2_GPLK_SPECTRUM_NATIVE_H
#include "TKspectrum.h"
#include "TKmatrix.h"
#include <memory>

static bool gplkSpectrumIncludes(const pulsar &p,const observation &o) {
    if (o.deleted) return false;
    const parameter &start=p.param[param_start],&finish=p.param[param_finish];
    if (start.aSize>0 && start.paramSet[0] && start.fitFlag[0] &&
        (start.fitFlag[0]==2 ? o.bat : o.sat)<start.val[0]-START_FINISH_DELTA) return false;
    if (finish.aSize>0 && finish.paramSet[0] && finish.fitFlag[0] &&
        (finish.fitFlag[0]==2 ? o.bat : o.sat)>finish.val[0]+START_FINISH_DELTA) return false;
    return true;
}
struct GplkSpectrumTimeOrder {
    bool operator()(const GplkSpectrumSample &a,const GplkSpectrumSample &b) const { return a.day<b.day; }
};
struct GplkSpectrumMatrix {
    double **data;
    explicit GplkSpectrumMatrix(int n):data(malloc_uinv(n)) {}
    ~GplkSpectrumMatrix() { if (data) free_uinv(data); }
private:
    GplkSpectrumMatrix(const GplkSpectrumMatrix &);
    GplkSpectrumMatrix &operator=(const GplkSpectrumMatrix &);
};

static bool gplkComputeNativeSpectrum(const pulsar &p,int mode,GplkSpectrumMethod method,
                                     const std::string &model,double fmin,double fmax,
                                     GplkSpectrum &out,std::string &error)
{
    out=GplkSpectrum();
    if (fmin<0 || fmax<0 || !gplkFinite(fmin) || !gplkFinite(fmax) || (fmax && fmin>=fmax)) {
        error="invalid native frequency bounds (cycles/day)"; return false;
    }
    std::vector<GplkSpectrumSample> data;
    for (int i=0;i<p.nobs;++i) {
        const observation &o=p.obsn[i];
        if (!gplkSpectrumIncludes(p,o)) continue;
        const double y=mode>=GPLK_NOISE_RED ? gplkNoiseValue(p,o,mode) :
                       mode==1 ? double(o.prefitResidual) : double(o.residual);
        const double day=double(o.sat),sigma=o.toaErr*1e-6;
        if (gplkFinite(day) && gplkFinite(y) && gplkFinite(sigma) && sigma>0)
            data.push_back(GplkSpectrumSample(day,y,sigma,i));
    }
    std::stable_sort(data.begin(),data.end(),GplkSpectrumTimeOrder());
    const int n=data.size();
    if (n<5 || !(data.back().day>data.front().day)) {
        error="at least five valid TOAs spanning more than one epoch are needed"; return false;
    }
    std::string covariance=model;
    if (method==GPLK_SPECTRUM_CHOLESKY) {
        if (covariance.empty() || covariance=="NULL") {
            error="choose a covariance model with Model, -spectrum-model or -dcf"; return false;
        }
        if (covariance=="PSRJ") covariance=std::string("covarFunc.dat_")+p.name;
        if (covariance.compare(0,5,"MODEL")!=0) {
            FILE *file=std::fopen(covariance.c_str(),"r");
            if (!file) { error="cannot read covariance model: "+covariance; return false; }
            const int first=std::fgetc(file); std::fclose(file);
            if (first==EOF) { error="covariance model is empty"; return false; }
        }
    }
    const double span=data.back().day-data.front().day;
    std::vector<double> x(n),y(n),sigma(n),frequency(n+2),power(n+2),real(n+2),imag(n+2);
    std::vector<int> ids(n);
    for (int i=0;i<n;++i) {
        // Match cholSpectra's site-time convention. Keep times/values/errors/IDs
        // together when sorting; TKspectrum sorts only x and y internally.
        x[i]=data[i].day-52601.0; y[i]=data[i].value; sigma[i]=data[i].sigma; ids[i]=data[i].observation;
    }
    int bins=0;
    if (method==GPLK_SPECTRUM_TEMPO2) {
        TKspectrum(&x[0],&y[0],&sigma[0],n,0,0,0,0,0,4,1,1,1,
                   &frequency[0],&power[0],&bins,0,0,&real[0],&imag[0]);
    } else {
        GplkSpectrumMatrix lower(n);
        // getCholeskyMatrix writes detL. Use a metadata copy so diagnostics do
        // not change the live fit; model parsing borrows observations read-only.
        std::auto_ptr<pulsar> context(new pulsar(p));
        // Native covariance construction expects extra covariance in the same
        // order as these filtered, sorted rows, not original observation order.
        std::auto_ptr<GplkSpectrumMatrix> extra;
        if (p.ToAextraCovar) {
            extra.reset(new GplkSpectrumMatrix(n));
            for (int i=0;i<n;++i) for (int j=0;j<n;++j)
                extra->data[i][j]=p.ToAextraCovar[ids[i]][ids[j]];
            context->ToAextraCovar=extra->data;
        }
        getCholeskyMatrix(lower.data,covariance.c_str(),context.get(),&x[0],&y[0],&sigma[0],n,0,&ids[0]);
        bins=calcSpectra(lower.data,&x[0],&y[0],n,&frequency[0],&power[0],-1);
    }
    out.psd=true; out.samples=n; out.span=span;
    if (method==GPLK_SPECTRUM_CHOLESKY) out.covarianceModel=covariance;
    for (int i=0;i<bins;++i) {
        if (!gplkFinite(frequency[i]) || !gplkFinite(power[i]) || frequency[i]<=0 || power[i]<0) continue;
        if ((fmin && frequency[i]<fmin) || (fmax && frequency[i]>fmax)) continue;
        out.frequency.push_back(frequency[i]); out.power.push_back(power[i]);
        if (power[i]>out.power[out.peak]) out.peak=out.power.size()-1;
    }
    if (out.frequency.size()<2) { error="fewer than two native Fourier bins in the requested range"; return false; }
    if (!(out.power[out.peak]>0)) { error="the selected series has no positive native PSD samples"; return false; }
    return true;
}
#endif
