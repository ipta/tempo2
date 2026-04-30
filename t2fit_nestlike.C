#include <cmath>
#include <assert.h>
#include <string.h>
#include <cstdio>
#include "t2fit_nestlike.h"
#include "shapelet.h"
#include "enum_str.h"

#include <map>
#include <tuple>
#include <vector>

struct CoeffCache {
    bool computed = false;
    double c[4] = {0.0, 0.0, 0.0, 0.0}; // up to cubic
    double start = 0.0;
    double finish = 0.0;
    double pepoch = 0.0;
    int nobs = 0;
};

double legendre(double x, int n) {
    switch (n) {
        case 0: return 1.0;
        case 1: return x;
        case 2: return 0.5 * (3*x*x - 1);
        case 3: return 0.5 * (5*x*x*x - 3*x);
        default: assert(0); return 0.0; // Only support up to n=3
    }
}

// Compute best-fit Legendre polynomial coefficients for projecting a pure
// sin or cos basis function (with no frequency-dependent prefactor) onto the
// fitted polynomial orders.  Fills entry.c[] via least-squares normal
// equations and marks entry.computed = true.
static void computePolySubtractCoeffs(
        pulsar *psr, int ipsr, param_label poly_param,
        double freq, bool is_cos,
        double alpha, double beta, double pepoch,
        CoeffCache &entry)
{
    int N = psr[ipsr].nobs;
    if (N < 2 || alpha == 0.0) {
        entry.computed = true;
        return;
    }

    int fitOrders[4];
    int nfit = 0;
    if (poly_param==param_f) {
        fitOrders[nfit++] = 0; // Always fit the constant term for red noise subtraction
        for (int n = 1; n < 4; n++) {
            if (psr[ipsr].param[param_f].fitFlag[n-1] == 1)
                fitOrders[nfit++] = n;
        }
    } else {
        for (int n = 0; n < 4; n++) {
            if (psr[ipsr].param[poly_param].fitFlag[n] == 1)
                fitOrders[nfit++] = n;
        }
    } 

    if (nfit == 0) {
        entry.computed = true;
        return;
    }

    std::vector<double> xi(N), f(N);
    for (int i = 0; i < N; i++) {
        double xj = psr[ipsr].obsn[i].bbat - pepoch;
        xi[i] = (xj - beta) / alpha;
        f[i]  = is_cos ? cos(2.0 * M_PI * freq * xj)
                       : sin(2.0 * M_PI * freq * xj);
    }

    // Normal equations: (A^T A) c = A^T f,  A[i][j] = P_{fitOrders[j]}(xi[i])
    double ATA[4][4] = {};
    double ATf[4]    = {};
    for (int i = 0; i < N; i++) {
        double P[4];
        for (int j = 0; j < nfit; j++)
            P[j] = legendre(xi[i], fitOrders[j]);
        for (int j = 0; j < nfit; j++) {
            ATf[j] += P[j] * f[i];
            for (int l = 0; l < nfit; l++)
                ATA[j][l] += P[j] * P[l];
        }
    }

    // Gauss-Jordan elimination with partial pivoting
    double aug[4][5] = {};
    for (int j = 0; j < nfit; j++) {
        for (int l = 0; l < nfit; l++)
            aug[j][l] = ATA[j][l];
        aug[j][nfit] = ATf[j];
    }
    for (int col = 0; col < nfit; col++) {
        int pivot = col;
        for (int row = col+1; row < nfit; row++)
            if (fabs(aug[row][col]) > fabs(aug[pivot][col]))
                pivot = row;
        for (int l = 0; l <= nfit; l++)
            std::swap(aug[col][l], aug[pivot][l]);
        if (fabs(aug[col][col]) < 1e-14) continue; // singular column
        double inv = 1.0 / aug[col][col];
        for (int row = 0; row < nfit; row++) {
            if (row == col) continue;
            double factor = aug[row][col] * inv;
            for (int l = col; l <= nfit; l++)
                aug[row][l] -= factor * aug[col][l];
        }
        for (int l = col+1; l <= nfit; l++)
            aug[col][l] *= inv;
        aug[col][col] = 1.0;
    }

    for (int j = 0; j < nfit; j++)
        entry.c[fitOrders[j]] = aug[j][nfit];
    entry.computed = true;
}

double t2FitFunc_nestlike_red(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k) {
    static std::map<std::tuple<int,int,param_label>, CoeffCache> cache;
    double maxtspan = psr[ipsr].param[param_finish].val[0] - psr[ipsr].param[param_start].val[0];
    double RedFLow = pow(10., psr[ipsr].TNRedFLow);
    double freq; 


    bool is_cos;
    switch (label) {
        case param_red_cos: is_cos = true;  break;
        case param_red_sin: is_cos = false; break;
        default: assert(0); return 0.0;
    }

    if (k >= psr[ipsr].TNRedC){ // This is in the log freq zone!
        int subharm = k - psr[ipsr].TNRedC + 1;
        double freq0 = RedFLow*(1.0/maxtspan);
        freq = freq0 * pow(psr[ipsr].TNRed_log_factor,-subharm);
    } else {
        freq = RedFLow*((double)(k+1.0))/(maxtspan);
    }

    // Compute orthogonalised trig value (no freq-dependent prefactor)
    double trig = is_cos ? cos(2.0*M_PI*freq*x) : sin(2.0*M_PI*freq*x);

    if (psr[ipsr].TNsubtractPoly&2){
        // subtract polynomial from the basis
        double a = psr[ipsr].param[param_start].val[0];
        double b = psr[ipsr].param[param_finish].val[0];
        double alpha = (b - a)/2.0;
        double beta  = (b + a)/2.0;

        auto key = std::make_tuple(ipsr, k, label);        
        const double pepoch = psr[ipsr].param[param_pepoch].val[0];
        const int nobs = psr[ipsr].nobs;


        auto &entry = cache[key];
        if (entry.computed &&
            (entry.start != a || entry.finish != b || entry.pepoch != pepoch || entry.nobs != nobs)) {
            entry.computed = false;
            entry.c[0] = 0.0;
            entry.c[1] = 0.0;
            entry.c[2] = 0.0;
            entry.c[3] = 0.0;
        }

        if (!entry.computed) {
            logmsg("Computing polynomial subtraction coefficients for pulsar %d, k=%d, label=%s", ipsr, k, label_str[label]);
            computePolySubtractCoeffs(psr, ipsr, param_f, freq, is_cos,
                                      alpha, beta, pepoch, entry);
            entry.start  = a;
            entry.finish = b;
            entry.pepoch = pepoch;
            entry.nobs   = nobs;
        }

        

        if (alpha != 0.0) {
            double xi_val = (x - beta)/alpha;
            for (int n = 0; n < 4; n++) {
                trig -= entry.c[n] * legendre(xi_val, n);
            }
        }
    }
    
    return trig;
}

double t2FitFunc_nestlike_red_dm(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k) {
    double maxtspan = psr[ipsr].param[param_finish].val[0] - psr[ipsr].param[param_start].val[0];
    double freq;
    if (k >= psr[ipsr].TNDMC){ // This is in the log freq zone!
        int subharm = k - psr[ipsr].TNDMC + 1;
        double freq0 = (1.0/maxtspan);
        freq = freq0 * pow(psr[ipsr].TNDM_log_factor,-subharm);
    } else {
        freq = ((double)(k+1.0))/(maxtspan);
    }

    bool is_cos;
    switch (label) {
        case param_red_dm_cos: is_cos = true;  break;
        case param_red_dm_sin: is_cos = false; break;
        default: assert(0); return 0.0;
    }

    // Compute orthogonalised trig value (no freq-dependent prefactor)
    double trig = is_cos ? cos(2.0*M_PI*freq*x) : sin(2.0*M_PI*freq*x);

    if (psr[ipsr].TNsubtractPoly&2) {
        static std::map<std::tuple<int,int,param_label>, CoeffCache> cache;
        double a      = psr[ipsr].param[param_start].val[0];
        double b      = psr[ipsr].param[param_finish].val[0];
        double alpha  = (b - a) / 2.0;
        double beta   = (b + a) / 2.0;
        double pepoch = psr[ipsr].param[param_pepoch].val[0];
        int    nobs   = psr[ipsr].nobs;
        auto   key    = std::make_tuple(ipsr, k, label);
        auto  &entry  = cache[key];

        if (entry.computed &&
            (entry.start != a || entry.finish != b || entry.pepoch != pepoch || entry.nobs != nobs)) {
            entry = CoeffCache{};
        }
        if (!entry.computed) {
            logmsg("Computing polynomial subtraction coefficients (dm) for pulsar %d, k=%d, label=%s", ipsr, k, label_str[label]);
            computePolySubtractCoeffs(psr, ipsr, param_dm, freq, is_cos, alpha, beta, pepoch, entry);
            entry.start  = a;
            entry.finish = b;
            entry.pepoch = pepoch;
            entry.nobs   = nobs;
        }

        if (alpha != 0.0) {
            double xi_val = (x - beta) / alpha;
            for (int n = 0; n < 4; n++)
                trig -= entry.c[n] * legendre(xi_val, n);
        }
    }

    double kappa = DM_CONST*1e-12;
    return trig / (kappa * pow((double)psr[ipsr].obsn[ipos].freqSSB, 2));
}

double t2FitFunc_nestlike_red_chrom(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k) {
    double maxtspan = psr[ipsr].param[param_finish].val[0] - psr[ipsr].param[param_start].val[0];
    double freq;
    if (k >= psr[ipsr].TNChromC){ // This is in the log freq zone!
        int subharm = k - psr[ipsr].TNChromC + 1;
        double freq0 = (1.0/maxtspan);
        freq = freq0 * pow(psr[ipsr].TNChrom_log_factor,-subharm);
    } else {
        freq = ((double)(k+1.0))/(maxtspan);
    }

    bool is_cos;
    switch (label) {
        case param_red_chrom_cos: is_cos = true;  break;
        case param_red_chrom_sin: is_cos = false; break;
        default: assert(0); return 0.0;
    }

    // Compute orthogonalised trig value (no freq-dependent prefactor)
    double trig = is_cos ? cos(2.0*M_PI*freq*x) : sin(2.0*M_PI*freq*x);

    if (psr[ipsr].TNsubtractPoly&2) {
        static std::map<std::tuple<int,int,param_label>, CoeffCache> cache;
        double a      = psr[ipsr].param[param_start].val[0];
        double b      = psr[ipsr].param[param_finish].val[0];
        double alpha  = (b - a) / 2.0;
        double beta   = (b + a) / 2.0;
        double pepoch = psr[ipsr].param[param_pepoch].val[0];
        int    nobs   = psr[ipsr].nobs;
        auto   key    = std::make_tuple(ipsr, k, label);
        auto  &entry  = cache[key];

        if (entry.computed &&
            (entry.start != a || entry.finish != b || entry.pepoch != pepoch || entry.nobs != nobs)) {
            entry = CoeffCache{};
        }
        if (!entry.computed) {
            logmsg("Computing polynomial subtraction coefficients (chrom) for pulsar %d, k=%d, label=%s", ipsr, k, label_str[label]);
            computePolySubtractCoeffs(psr, ipsr, param_cm, freq, is_cos, alpha, beta, pepoch, entry);
            entry.start  = a;
            entry.finish = b;
            entry.pepoch = pepoch;
            entry.nobs   = nobs;
        }

        if (alpha != 0.0) {
            double xi_val = (x - beta) / alpha;
            for (int n = 0; n < 4; n++)
                trig -= entry.c[n] * legendre(xi_val, n);
        }
    }

    double index  = psr[ipsr].TNChromIdx;
    double prefac = 1.0;
    return prefac * trig / pow((double)psr[ipsr].obsn[ipos].freqSSB / 1.4e9, index);
}





void t2UpdateFunc_nestlike_red(pulsar *psr, int ipsr ,param_label label,int k, double val, double err) {
    if (k==0 && label==param_red_sin){
        double maxtspan = psr[ipsr].param[param_finish].val[0] - psr[ipsr].param[param_start].val[0];
        double RedFLow = pow(10., psr[ipsr].TNRedFLow);
        double freq = RedFLow*((double)(k+1.0))/(maxtspan);
        //printf("WAVE_FREQ %.18lg\n",freq*365.25);
        //printf("Temponest equivilent fitwaves parameters\n");
        //printf("========================================\n");
        //printf("\n");
        //printf("WAVE_OM %.18lg\n",freq*2.0*M_PI);
        //printf("WAVEEPOCH %.18lg\n",(double)psr[ipsr].param[param_pepoch].val[0]);
        if (writeResiduals&4){
            FILE *fout;
            fout = fopen("tnred.meta","w");
            if (!fout){
                printf("Unable to open file tnred.meta for writing\n");
            }
            fprintf(fout,"RED_OMEGA %lg\n",freq*2.0*M_PI);
            fprintf(fout,"RED_EPOCH %lg\n",(double)psr[ipsr].param[param_pepoch].val[0]);
            fclose(fout);
        }
    }
    /*
    if (label==param_red_sin){
        printf("WAVE_SIN%d %.18lg\n",k+1,-val);
    }
    if (label==param_red_cos){
        printf("WAVE_COS%d %.18lg\n",k+1,-val);
    }
    printf("\n");
    */

    //    logmsg("%d %s %d %lg %lg",ipsr,label_str[label],k,val,err);
    for (int iobs = 0; iobs < psr[ipsr].nobs; ++iobs){
        double x = (double)(psr[ipsr].obsn[iobs].bbat - psr[ipsr].param[param_pepoch].val[0]);
        double y = t2FitFunc_nestlike_red(psr,ipsr,x,iobs,label,k);
        psr[ipsr].obsn[iobs].TNRedSignal  += y *val;
        psr[ipsr].obsn[iobs].TNRedErr     += pow(y*err,2);
        //fprintf(stderr, "are we here????\n");

    }
}


void t2UpdateFunc_nestlike_red_dm(pulsar *psr, int ipsr ,param_label label,int k, double val, double err) {

    if (k==0 && label==param_red_dm_sin){
        if (writeResiduals&4){
            double maxtspan = psr[ipsr].param[param_finish].val[0] - psr[ipsr].param[param_start].val[0];
            double freq = ((double)(k+1.0))/(maxtspan);
            FILE *fout;
            fout = fopen("tnreddm.meta","w");
            if (!fout){
                printf("Unable to open file tnred.meta for writing\n");
            }
            fprintf(fout,"REDDM_OMEGA %lg\n",freq*2.0*M_PI);
            fprintf(fout,"REDDM_EPOCH %lg\n",(double)psr[ipsr].param[param_pepoch].val[0]);
            fprintf(fout,"REDDM_DMEPOCH %lg\n",(double)psr[ipsr].param[param_dmepoch].val[0]);
            fprintf(fout,"REDDM_DM %lg\n",(double)psr[ipsr].param[param_dm].prefit[0]);
            fprintf(fout,"REDDM_DM1 %lg\n",(double)psr[ipsr].param[param_dm].prefit[1]);
            fprintf(fout,"REDDM_DM2 %lg\n",(double)psr[ipsr].param[param_dm].prefit[2]);
            fclose(fout);
        }
    }
    logdbg("%d %s %d %lg %lg",ipsr,label_str[label],k,val,err);
    for (int iobs = 0; iobs < psr[ipsr].nobs; ++iobs){
        double x = (double)(psr[ipsr].obsn[iobs].bbat - psr[ipsr].param[param_pepoch].val[0]);
        double y = t2FitFunc_nestlike_red_dm(psr,ipsr,x,iobs,label,k);
        psr[ipsr].obsn[iobs].TNDMSignal  += y *val;
        psr[ipsr].obsn[iobs].TNDMErr     += pow(y*err,2);

    }
}



void t2UpdateFunc_nestlike_red_chrom(pulsar *psr, int ipsr ,param_label label,int k, double val, double err) {

    if (k==0 && label==param_red_chrom_sin){
        if (writeResiduals&4){
            double maxtspan = psr[ipsr].param[param_finish].val[0] - psr[ipsr].param[param_start].val[0];
            double freq = ((double)(k+1.0))/(maxtspan);
            FILE *fout; 
            fout = fopen("tnreddm.meta","w");
            if (!fout){
                printf("Unable to open file tnred.meta for writing\n");
            }
            fprintf(fout,"REDDM_OMEGA %lg\n",freq*2.0*M_PI);
            fprintf(fout,"REDDM_EPOCH %lg\n",(double)psr[ipsr].param[param_pepoch].val[0]);
            fprintf(fout,"REDDM_DMEPOCH %lg\n",(double)psr[ipsr].param[param_dmepoch].val[0]);
            fprintf(fout,"REDDM_DM %lg\n",(double)psr[ipsr].param[param_dm].prefit[0]);
            fprintf(fout,"REDDM_DM1 %lg\n",(double)psr[ipsr].param[param_dm].prefit[1]);
            fprintf(fout,"REDDM_DM2 %lg\n",(double)psr[ipsr].param[param_dm].prefit[2]);
            fclose(fout);
        }
    }
    logmsg("%d %s %d %lg %lg",ipsr,label_str[label],k,val,err);
    for (int iobs = 0; iobs < psr[ipsr].nobs; ++iobs){
        double x = (double)(psr[ipsr].obsn[iobs].bbat - psr[ipsr].param[param_pepoch].val[0]);
        double y = t2FitFunc_nestlike_red_chrom(psr,ipsr,x,iobs,label,k);
        psr[ipsr].obsn[iobs].TNChromSignal  += y *val;
        psr[ipsr].obsn[iobs].TNChromErr     += pow(y*err,2);

    }
}




double t2FitFunc_nestlike_jitter(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k){
    const double dt=1./SECDAY;
    const double epoch = psr[ipsr].obsn[k].bat;
    const double batmin = epoch-dt;
    const double batmax = epoch+dt;
    const double bat = psr[ipsr].obsn[ipos].bat;
    if (bat > batmin && bat < batmax) return 1;
    else return 0;


}

void t2UpdateFunc_nestlike_jitter(pulsar *psr, int ipsr ,param_label label,int k, double val, double err) {
    logmsg("%d %s %.2f %lg %lg",ipsr,label_str[label],(double)psr[ipsr].obsn[k].bat,val,err);
}





double t2FitFunc_nestlike_band(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k){
    //int totalBandNoiseCoeff=0; // unused

    int iband = 0;
    int ichan = k;
    while (ichan >= psr[ipsr].TNBandNoiseC[iband]){
        ichan -= psr[ipsr].TNBandNoiseC[iband];
        ++iband;
    }

    // we are in channel ichan, band iband.


    double ret = 0;

    double BandLF = psr[ipsr].TNBandNoiseLF[iband];
    double BandHF = psr[ipsr].TNBandNoiseHF[iband];
    //double BandAmp=pow(10.0, psr[ipsr].TNBandNoiseAmp[iband]);
    //double BandSpec=psr[ipsr].TNBandNoiseGam[iband];


    double maxtspan = psr[ipsr].param[param_finish].val[0] - psr[ipsr].param[param_start].val[0];
    double freq = ((double)(ichan+1.0))/(maxtspan);

    // check that we are indeed in a band
    if (psr[ipsr].obsn[ipos].freq > BandLF && psr[ipsr].obsn[ipos].freq < BandHF) {
        switch (label){
            case param_band_red_cos:
                ret = cos(2.0*M_PI*freq*x);
                break;
            case param_band_red_sin:
                ret = sin(2.0*M_PI*freq*x);
                break;
            default:
                assert(0);
                break;
        }
    }
    return ret;
}


double t2FitFunc_nestlike_group(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k){

    //int totalGroupCoeff=0;

    int igroup = 0;
    int ichan = k;
    while (ichan >= psr[ipsr].TNGroupNoiseC[igroup]){
        ichan -= psr[ipsr].TNGroupNoiseC[igroup];
        ++igroup;
    }

    // we are in channel ichan, group igroup.


    double ret=0.0;

    //double GroupAmp=pow(10.0, psr[ipsr].TNGroupNoiseAmp[igroup]);
    //double GroupSpec=psr[ipsr].TNGroupNoiseGam[igroup];
    //int GroupC=psr[ipsr].TNGroupNoiseC[igroup];

    double maxtspan = psr[ipsr].param[param_finish].val[0] - psr[ipsr].param[param_start].val[0];
    double freq = ((double)(ichan+1.0))/(maxtspan);

    bool ingrp=false;

    for (int j=0; j < psr[ipsr].obsn[ipos].nFlags; j++) {
        //Check Group Noise Flag
        if (strcmp(
                    psr[ipsr].obsn[ipos].flagID[j],
                    psr[ipsr].TNGroupNoiseFlagID[igroup]
                  )==0 ) {
            if (strcmp(
                        psr[ipsr].obsn[ipos].flagVal[j],
                        psr[ipsr].TNGroupNoiseFlagVal[igroup]
                      )==0 ) {
                ingrp=true;
                break;
            }
        }
    }

    if(ingrp){ // we are in this group
        switch (label){
            case param_group_red_cos:
                ret = cos(2.0*M_PI*freq*x);
                break;
            case param_group_red_sin:
                ret = sin(2.0*M_PI*freq*x);
                break;
            default:
                assert(0);
                break;
        }

    }
    return ret;

}



void t2UpdateFunc_nestlike_band(pulsar *psr, int ipsr ,param_label label,int k, double val, double err) {
    int iband = 0;
    int ichan = k;
    while (ichan >= psr[ipsr].TNBandNoiseC[iband]){
        ichan -= psr[ipsr].TNBandNoiseC[iband];
        ++iband;
    }
    
    // we are in channel ichan, band iband.
    logmsg("%d %s %d %d %d %lg %lg",ipsr,label_str[label],k,iband,ichan,val,err);
    for (int iobs = 0; iobs < psr[ipsr].nobs; ++iobs){
        double x = (double)(psr[ipsr].obsn[iobs].bbat - psr[ipsr].param[param_pepoch].val[0]);
        double y = t2FitFunc_nestlike_band(psr,ipsr,x,iobs,label,k);
        psr[ipsr].obsn[iobs].TNRedSignal  += y *val;
        psr[ipsr].obsn[iobs].TNRedErr     += pow(y*err,2);

    }
}

void t2UpdateFunc_nestlike_group(pulsar *psr, int ipsr ,param_label label,int k, double val, double err) {
    int igroup = 0;
    int ichan = k;
    while (ichan >= psr[ipsr].TNGroupNoiseC[igroup]){
        ichan -= psr[ipsr].TNGroupNoiseC[igroup];
        ++igroup;
    }

    // we are in channel ichan, group igroup.
    logmsg("%d %s %d %d %d %lg %lg",ipsr,label_str[label],k,igroup,ichan,val,err);
    for (int iobs = 0; iobs < psr[ipsr].nobs; ++iobs){
        double x = (double)(psr[ipsr].obsn[iobs].bbat - psr[ipsr].param[param_pepoch].val[0]);
        double y = t2FitFunc_nestlike_group(psr,ipsr,x,iobs,label,k);
        psr[ipsr].obsn[iobs].TNRedSignal  += y *val;
        psr[ipsr].obsn[iobs].TNRedErr     += pow(y*err,2);

    }
}


double t2FitFunc_nestlike_shape_red(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k) {
    int ishape = k/MAX_TNShapeCoef;
    int icoef = k%MAX_TNShapeCoef;
    double pos = psr[ipsr].TNShapeletEvPos[ishape];
    double width = psr[ipsr].TNShapeletEvWidth[ishape];
    double t = psr[ipsr].obsn[ipos].bat;
    int N=icoef+1;
    double shapecoef[N];
    for (int ic=0; ic < N ;++ic)shapecoef[ic]=0;
    shapecoef[icoef]=1.0;
    double shape = evaluateShapelet(N,pos,width,shapecoef,t);
    return shape;
}

double t2FitFunc_nestlike_shape_dm(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k) {
    int ishape = k/MAX_TNShapeCoef;
    int icoef = k%MAX_TNShapeCoef;
    double pos = psr[ipsr].TNShapeletEvPos[ishape];
    double width = psr[ipsr].TNShapeletEvWidth[ishape];
    double t = psr[ipsr].obsn[ipos].bat;
    int N=icoef+1;
    double shapecoef[N];
    for (int ic=0; ic < N ;++ic)shapecoef[ic]=0;
    shapecoef[icoef]=1.0;
    double shape = evaluateShapelet(N,pos,width,shapecoef,t);
    return shape *pow((double)psr[ipsr].obsn[ipos].freqSSB,-2)*1e12/DM_CONST;
}

void t2UpdateFunc_nestlike_shape(pulsar *psr, int ipsr ,param_label label,int k, double val, double err) {
    int ishape = k/MAX_TNShapeCoef;
    int icoef = k%MAX_TNShapeCoef;
    psr[ipsr].TNShapeletEvCoef[ishape][icoef] += val;
    psr[ipsr].TNShapeletEvCoefErr[ishape][icoef] = err;
}
