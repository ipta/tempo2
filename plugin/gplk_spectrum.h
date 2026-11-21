/**
 * @file
 * @brief Weighted GLS and native TEMPO2 spectral diagnostics. GPLv3-or-later.
 * GLS normalization: (chi2_constant - chi2_sinusoid) / chi2_constant.
 * See Zechmeister & Kuerster (2009), doi:10.1051/0004-6361:200811296.
 * Pure calculation is separated from Giza/device handling below. See also
 * documentation/GPLK_GIZA_DEVELOPER_GUIDE.md for units and interpretation limits.
 */
#ifndef TEMPO2_GPLK_SPECTRUM_H
#define TEMPO2_GPLK_SPECTRUM_H
#include "gplk_noise.h"
#include "gplk_giza.h"
#include <algorithm>
#include <vector>

struct GplkSpectrumSample {
    // day: epoch in days; value and sigma must share units (seconds in gplk).
    double day, value, sigma;
    int observation;
    GplkSpectrumSample(double t, double y, double e,int id=-1) : day(t), value(y), sigma(e),observation(id) {}
};
struct GplkSpectrum {
    // Parallel arrays: cycles/day, standard power [0,1], weighted window power.
    std::vector<double> frequency, power, window;
    // Metadata are valid only after successful computation. peak indexes arrays.
    unsigned int peak;
    double span;
    unsigned int samples;
    bool psd;
    std::string covarianceModel;
    GplkSpectrum():peak(0),span(0),samples(0),psd(false) {}
};
static double gplkSpectrumMin = 0, gplkSpectrumMax = 0;
enum GplkSpectrumMethod { GPLK_SPECTRUM_GLS,GPLK_SPECTRUM_TEMPO2,GPLK_SPECTRUM_CHOLESKY };
static GplkSpectrumMethod gplkSpectrumMethod=GPLK_SPECTRUM_GLS;
static std::string gplkSpectrumModel,gplkSpectrumOutput;
static bool gplkSetSpectrumMethod(const char *name) {
    if (std::strcmp(name,"gls")==0) gplkSpectrumMethod=GPLK_SPECTRUM_GLS;
    else if (std::strcmp(name,"tempo2")==0) gplkSpectrumMethod=GPLK_SPECTRUM_TEMPO2;
    else if (std::strcmp(name,"cholesky")==0) gplkSpectrumMethod=GPLK_SPECTRUM_CHOLESKY;
    else return false;
    return true;
}
static const char *gplkSpectrumMethodName(GplkSpectrumMethod method) {
    return method==GPLK_SPECTRUM_GLS ? "Weighted GLS" :
           method==GPLK_SPECTRUM_TEMPO2 ? "TEMPO2 weighted LS" : "TEMPO2 Cholesky";
}

/**
 * @brief Fit constant + sine + cosine on a uniform frequency grid, without graphics.
 * @param input Samples; invalid values/nonpositive errors are excluded.
 * @param fmin Lower frequency in cycles/day; zero requests the span-based default.
 * @param fmax Upper frequency in cycles/day; zero requests the cadence-based default.
 * @param out Result arrays, cleared on entry; only read metadata after success.
 * @param error Explanation on failure; not cleared/defined as a status on success.
 * @return False for insufficient/degenerate data or invalid frequency bounds.
 */
static bool gplkComputeSpectrum(const std::vector<GplkSpectrumSample> &input,
                                double fmin, double fmax, GplkSpectrum &out,
                                std::string &error)
{
    out=GplkSpectrum();
    std::vector<GplkSpectrumSample> data;
    std::vector<double> times;
    double minimumSigma=DBL_MAX;
    for (unsigned int i=0; i<input.size(); ++i) {
        const GplkSpectrumSample &s=input[i];
        if (!gplkFinite(s.day) || !gplkFinite(s.value) || !gplkFinite(s.sigma) || s.sigma<=0) continue;
        data.push_back(s); times.push_back(s.day);
        minimumSigma=std::min(minimumSigma,s.sigma);
    }
    if (data.size()<5) { error="at least five finite TOAs with positive errors are needed"; return false; }
    std::sort(times.begin(),times.end());
    const double span=times.back()-times.front();
    if (!(span>0)) { error="TOAs must span more than one epoch"; return false; }
    std::vector<double> gaps;
    for (unsigned int i=1; i<times.size(); ++i)
        if (times[i]>times[i-1]) gaps.push_back(times[i]-times[i-1]);
    std::sort(gaps.begin(),gaps.end());
    if (fmin==0) fmin=1.0/span;
    // A convenient cadence-based default bound, not an irregular-sampling Nyquist limit.
    if (fmax==0) fmax=std::max(2*fmin,0.5/gaps[gaps.size()/2]);
    if (!gplkFinite(fmin) || !gplkFinite(fmax) || fmin<=0 || fmax<=fmin) {
        error="frequency bounds must satisfy 0 < fmin < fmax (cycles/day)"; return false;
    }
    // A common weight factor cancels from normalized power. Scaling by the
    // smallest sigma avoids directly forming very large 1/sigma^2 values.
    std::vector<double> weight(data.size()), centred(data.size());
    double sumWeight=0, offset=0;
    for (unsigned int i=0; i<data.size(); ++i) {
        const double ratio=minimumSigma/data[i].sigma;
        weight[i]=ratio*ratio;
        sumWeight+=weight[i]; offset+=weight[i]*(data[i].value-data[0].value);
    }
    // Compute the mean relative to one sample to reduce offset cancellation.
    const double mean=data[0].value+offset/sumWeight;
    double variance=0;
    for (unsigned int i=0; i<data.size(); ++i) {
        centred[i]=data[i].value-mean;
        variance+=weight[i]*centred[i]*centred[i];
    }
    if (!(variance>0) || !gplkFinite(variance)) {
        error="the selected series has no finite variation"; return false;
    }
    const double desired=std::ceil(5*span*(fmax-fmin))+1;
    const int bins=static_cast<int>(std::max(256.0,std::min(8192.0,desired)));
    out.peak=0; out.span=span; out.samples=data.size();
    for (int j=0; j<bins; ++j) {
        const double f=fmin+(fmax-fmin)*j/(bins-1);
        double c=0,s=0,cc=0,ss=0,cs=0,yc=0,ys=0;
        for (unsigned int i=0; i<data.size(); ++i) {
            const double angle=2*M_PI*f*(data[i].day-times.front());
            const double co=std::cos(angle), si=std::sin(angle), w=weight[i];
            c+=w*co; s+=w*si; cc+=w*co*co; ss+=w*si*si;
            cs+=w*co*si; yc+=w*centred[i]*co; ys+=w*centred[i]*si;
        }
        // Window = |sum(w * exp(i*phase))|^2 / sum(w)^2, not a fitted signal.
        out.window.push_back((c*c+s*s)/(sumWeight*sumWeight));
        // Remove weighted sine/cosine means: eliminate the floating intercept
        // from the 3-column normal equations, leaving a 2x2 least-squares solve.
        cc-=c*c/sumWeight; ss-=s*s/sumWeight; cs-=c*s/sumWeight;
        const double determinant=cc*ss-cs*cs;
        double power=0;
        if (cc>0 && ss>0 && determinant>DBL_EPSILON*cc*ss)
            power=(ss*yc*yc+cc*ys*ys-2*cs*yc*ys)/(determinant*variance);
        power=std::max(0.0,std::min(1.0,power));
        out.frequency.push_back(f); out.power.push_back(power);
        if (power>out.power[out.peak]) out.peak=j;
    }
    return true;
}

#include "gplk_spectrum_native.h"
#include <cerrno>

/** Display bounds are independent of spectral samples and survive redraw/resize.
 * x is log10(cycles/year); y is normalized GLS power or log10(PSD / yr^3). */
struct GplkSpectrumView {
    double x0,x1,y0,y1,fullX0,fullX1,fullY0,fullY1;
    bool initialize(const GplkSpectrum &s) {
        if (s.frequency.size()<2) return false;
        fullX0=std::log10(s.frequency.front()*365.25);
        fullX1=std::log10(s.frequency.back()*365.25);
        fullY0=0; fullY1=1.05;
        if (s.psd) {
            double lo=DBL_MAX,hi=-DBL_MAX;
            for (unsigned int i=0;i<s.power.size();++i) if (s.power[i]>0) {
                const double v=std::log10(s.power[i]); lo=std::min(lo,v); hi=std::max(hi,v);
            }
            if (lo==DBL_MAX) return false;
            const double pad=std::max(.2,(hi-lo)*.08);
            fullY0=lo-pad; fullY1=hi+pad;
        }
        reset(); return true;
    }
    void reset() { x0=fullX0; x1=fullX1; y0=fullY0; y1=fullY1; }
    bool rectangle(double ax,double ay,double bx,double by) {
        const double lo=std::max(x0,std::min(ax,bx)),hi=std::min(x1,std::max(ax,bx));
        const double low=std::max(y0,std::min(ay,by)),high=std::min(y1,std::max(ay,by));
        if (!(hi-lo>1e-8 && high-low>1e-8)) return false;
        x0=lo; x1=hi; y0=low; y1=high; return true;
    }
    void scale(double factor) {
        const double mx=(x0+x1)/2,my=(y0+y1)/2,dx=(x1-x0)*factor/2,dy=(y1-y0)*factor/2;
        x0=std::max(fullX0,mx-dx); x1=std::min(fullX1,mx+dx);
        y0=std::max(fullY0,my-dy); y1=std::min(fullY1,my+dy);
    }
};

static bool gplkBuildSpectrum(const pulsar &p,int mode,GplkSpectrumMethod method,
                              const std::string &model,GplkSpectrum &result,std::string &error)
{
    error.clear();
    if (mode>=GPLK_NOISE_RED && !gplkNoiseAvailable(p,mode)) {
        error=std::string("no reconstructed ")+gplkNoiseName(mode); return false;
    }
    if (method!=GPLK_SPECTRUM_GLS)
        return gplkComputeNativeSpectrum(p,mode,method,model,gplkSpectrumMin,gplkSpectrumMax,result,error);
    std::vector<GplkSpectrumSample> data;
    for (int i=0;i<p.nobs;++i) {
        const observation &o=p.obsn[i];
        if (!gplkSpectrumIncludes(p,o)) continue;
        const double value=mode>=GPLK_NOISE_RED ? gplkNoiseValue(p,o,mode) :
                           mode==1 ? double(o.prefitResidual) : double(o.residual);
        data.push_back(GplkSpectrumSample(double(o.bat),value,o.toaErr*1e-6,i));
    }
    return gplkComputeSpectrum(data,gplkSpectrumMin,gplkSpectrumMax,result,error);
}

static bool gplkWriteSpectrum(const GplkSpectrum &s,GplkSpectrumMethod method,int mode,
                              const std::string &path,std::string &error)
{
    // Exclusive creation protects previous scientific results from overwrites.
    FILE *file=std::fopen(path.c_str(),"wx");
    if (!file) { error=std::string("cannot save spectrum: ")+std::strerror(errno); return false; }
    std::fprintf(file,"# method: %s\n# series: %s\n# power units: %s\n",gplkSpectrumMethodName(method),
                  mode==1 ? "prefit" : gplkNoiseName(mode),s.psd ? "yr^3" : "fractional variance explained");
    if (!s.covarianceModel.empty()) {
        std::fprintf(file,"# covariance_model: ");
        for (unsigned int i=0;i<s.covarianceModel.size();++i) {
            const char c=s.covarianceModel[i];
            if (c=='\n') std::fputs("\\n",file);
            else if (c=='\r') std::fputs("\\r",file);
            else if (c=='\t') std::fputs("\\t",file);
            else std::fputc(c,file);
        }
        std::fprintf(file,"\n");
    }
    std::fprintf(file,"# frequency_cycles_per_day\tfrequency_cycles_per_year\tpower%s\n",s.psd ? "" : "\tweighted_window");
    for (unsigned int i=0;i<s.frequency.size();++i) {
        std::fprintf(file,"%.17g\t%.17g\t%.17g",s.frequency[i],s.frequency[i]*365.25,s.power[i]);
        if (!s.psd) std::fprintf(file,"\t%.17g",s.window[i]);
        std::fprintf(file,"\n");
    }
    const bool failed=std::ferror(file)!=0;
    const int closed=std::fclose(file);
    if (failed || closed!=0) { error="error writing spectrum file"; return false; }
    return true;
}

// Full-surface rectangles shared by spectrum toolbar drawing and input.
static GplkRect gplkSpectrumButton(int i,float height) {
    const int column=i%4,row=i/4;
    return GplkRect(column*.25f+.008f,(column+1)*.25f-.008f,
                   1-(32+row*30)/height,1-(8+row*30)/height);
}
static bool gplkSpectrumPrompt(const char *prompt,std::string &value) {
    float vx0,vx1,vy0,vy1,wx0,wx1,wy0,wy1,width,height;
    giza_get_viewport_float(GIZA_UNITS_NORMALIZED,&vx0,&vx1,&vy0,&vy1);
    giza_get_window_float(&wx0,&wx1,&wy0,&wy1);
    giza_get_paper_size_float(GIZA_UNITS_PIXELS,&width,&height);
    std::string input;
    bool accepted=false;
    giza_save();
    while (true) {
        giza_set_viewport_float(0,1,0,1); giza_set_window_float(0,1,0,1);
        giza_begin_buffer();
        giza_set_colour_index(0); giza_set_fill(GIZA_FILL_SOLID);
        giza_rectangle_float(.005f,.995f,0,90/height);
        giza_set_colour_index(1); giza_set_character_height_float(.85f*600/height);
        gplkLiteralText(.015f,66/height,prompt);
        gplkLiteralText(.015f,44/height,"Enter: apply | Esc: cancel");
        const int count=std::max(8,int(width/8)-4);
        const std::string visible=input.size()>unsigned(count) ? "..."+input.substr(input.size()-count) : input;
        gplkLiteralText(.015f,17/height,visible+"|");
        giza_end_buffer();
        float x=0,y=0; char key=0;
        if (!gplkReadEvent(&x,&y,&key) || key==27 || key==0) break;
        if (key==13 || key==10) { accepted=!input.empty(); break; }
        if ((key==8 || key==127) && !input.empty()) input.erase(input.size()-1);
        else if (key>=32 && key<=126 && input.size()<MAX_FILELEN-1) input+=key;
    }
    giza_restore(); giza_set_viewport_float(vx0,vx1,vy0,vy1); giza_set_window_float(wx0,wx1,wy0,wy1);
    if (accepted) value=input;
    return accepted;
}

/** Independent Giza child window. Methods call native numerical APIs only;
 * opening/changing a diagnostic never refits or reorders live observations. */
static bool gplkShowSpectrum(const pulsar &p,int mode,const char *device,bool hasParent=false,
                             const char *fitModel="NULL")
{
    GplkSpectrumMethod method=gplkSpectrumMethod;
    std::string model=gplkSpectrumModel.empty() ? fitModel : gplkSpectrumModel;
    GplkSpectrum result; GplkSpectrumView view;
    std::string status;
    if (!gplkBuildSpectrum(p,mode,method,model,result,status) || !view.initialize(result)) {
        std::fprintf(stderr,"gplk: cannot draw spectrum: %s\n",status.empty() ? "no positive PSD samples" : status.c_str());
        return false;
    }
    if (!gplkSpectrumOutput.empty() && !gplkWriteSpectrum(result,method,mode,gplkSpectrumOutput,status)) {
        std::fprintf(stderr,"gplk: %s\n",status.c_str()); return false;
    }
    int previous=0; if (hasParent) giza_get_device_id(&previous);
    if (!gplkOpenDevice(device,1,1)) return false;
    gplkSetAspect(.75f); gplkSetFont(1);
    const bool interactive=giza_device_has_cursor()!=0;
    const bool parentResponsiveFonts=gplkResponsiveFonts; gplkResponsiveFonts=false;
    bool close=false,zoom=false,sync=true;
    while (!close) {
        giza_begin_buffer();
        if (sync && interactive) giza_change_page();
        sync=false;
        float width,height; giza_get_paper_size_float(GIZA_UNITS_PIXELS,&width,&height);
        const float fontScale=600/std::max(1.0f,height);
        if (gplkTheme==GPLK_CLASSIC) {
            double r,g,b; giza_get_colour_representation(0,&r,&g,&b);
            const bool dark=r+g+b<1.5;
            gplkSetRGB(GPLK_GRID,dark ? 0x2c3b4e : 0xe1e8f0);
            gplkSetRGB(GPLK_ACCENT,dark ? 0x70b8ff : 0x2463a6);
        }
        giza_draw_background();
        giza_set_viewport_float(0,1,0,1); giza_set_window_float(0,1,0,1);
        const char *buttons[]={"GLS [1]","TEMPO2 WLS [2]","Cholesky [3]","Model [m]","Zoom [z]","Reset [u]","Save data [s]","Close [q]"};
        for (int i=0;interactive && i<8;++i) {
            const GplkRect b=gplkSpectrumButton(i,height);
            const bool active=i<3 ? int(method)==i : i==4 && zoom;
            const bool enabled=i!=2 || (!model.empty() && model!="NULL");
            gplkDrawButton(b.left,b.right,b.bottom,b.top,buttons[i],active,enabled,false,fontScale);
        }
        const GplkRect plot(std::min(.18f,85/width),1-std::min(.10f,55/width),
                           std::min(.38f,(interactive ? 140 : 100)/height),
                           1-std::min(.30f,(interactive ? 110 : 70)/height));
        plot.select(); giza_set_window(view.x0,view.x1,view.y0,view.y1);
        giza_set_character_height(.95*fontScale); gplkDrawPanel(view.x0,view.x1,view.y0,view.y1);
        giza_set_colour_index(1); giza_set_line_width(1);
        giza_box("BCLNST",0,0,result.psd ? "BCLNST" : "BNST",0,0);
        char title[256]; std::snprintf(title,sizeof(title),"%s | %s | %s",p.name,
                                      mode==1 ? "Prefit" : gplkNoiseName(mode),gplkSpectrumMethodName(method));
        giza_label("Frequency (cycles/year)",result.psd ? "PSD (yr^3)" : "Fractional variance explained",title);
        std::vector<double> x,y;
        for (unsigned int i=0;i<result.frequency.size();++i) {
            if (result.psd && result.power[i]<=0) continue;
            x.push_back(std::log10(result.frequency[i]*365.25));
            y.push_back(result.psd ? std::log10(result.power[i]) : result.power[i]);
        }
        if (!result.psd) {
            giza_set_colour_index(gplkTheme==GPLK_CLASSIC ? 14 : GPLK_MUTED);
            giza_set_line_style(GIZA_LS_DASHED); giza_line(x.size(),&x[0],&result.window[0]);
        }
        giza_set_colour_index(2); giza_set_line_style(GIZA_LS_SOLID); giza_set_line_width(1.5);
        giza_line(x.size(),&x[0],&y[0]);
        const double peak=result.frequency[result.peak];
        giza_set_colour_index(3);
        giza_single_point(std::log10(peak*365.25),result.psd ? std::log10(result.power[result.peak]) : result.power[result.peak],17);
        char annotation[256];
        std::snprintf(annotation,sizeof(annotation),"Peak %.5g days | %.4g %s | %u TOAs",1/peak,result.power[result.peak],
                      result.psd ? "yr^3" : "power",result.samples);
        giza_set_colour_index(1); giza_set_character_height(.75*fontScale);
        giza_annotate("T",1.1,0,0,annotation);
        giza_set_viewport_float(0,1,0,1); giza_set_window_float(0,1,0,1);
        giza_set_character_height(.75*fontScale);
        // Place footer rows in pixels: Giza annotation offsets scale with each
        // line's font and would otherwise overlap the larger frequency label.
        giza_text_float(plot.left,(interactive ? 76 : 35)/height,
                        result.psd ? "Native Fourier grid | site epochs | PSD in yr^3" :
                        "Weighted floating-mean GLS | barycentric epochs | dashed: sampling window");
        if (interactive)
            giza_text_float(plot.left,55/height,"z: two-corner zoom | +/-: zoom | u/right click: reset | q: close");
        giza_text_float(plot.left,(interactive ? 35 : 15)/height,
                        "Valid non-deleted TOAs in START/FINISH range | display zoom does not resample");
        const std::string message=zoom ? "Zoom: click two opposite corners; Esc cancels" :
                                  !status.empty() ? status :
                                  method==GPLK_SPECTRUM_CHOLESKY ? "Covariance model: "+result.covarianceModel : "";
        gplkLiteralText(.015f,8/height,message.substr(0,std::max(10,int(width/7))));
        plot.select(); giza_set_window(view.x0,view.x1,view.y0,view.y1);
        giza_end_buffer();
        if (!interactive) break;
        float mx=0,my=0; char key=0;
        if (!gplkReadEvent(&mx,&my,&key)) break;
        if (key==0 || key=='R') { sync=true; continue; }
        const float nx=plot.left+(mx-view.x0)/(view.x1-view.x0)*(plot.right-plot.left);
        const float ny=plot.bottom+(my-view.y0)/(view.y1-view.y0)*(plot.top-plot.bottom);
        int button=-1;
        if (key=='A') for (int i=0;i<8;++i) if (gplkSpectrumButton(i,height).contains(nx,ny)) button=i;
        if (key>='1' && key<='3') button=key-'1';
        if (key=='q' || button==7) { close=true; continue; }
        if (key==27) { zoom=false; status.clear(); continue; }
        if (key=='u' || key=='r' || key=='X' || button==5) { view.reset(); zoom=false; status.clear(); continue; }
        if (key=='+' || key=='=') { view.scale(.5); continue; }
        if (key=='-' || key=='_') { view.scale(2); continue; }
        if (key=='z' || button==4) { zoom=!zoom; continue; }
        if (key=='s' || button==6) {
            std::string path;
            if (gplkSpectrumPrompt("Save full spectrum to a new file",path)) {
                if (gplkWriteSpectrum(result,method,mode,path,status)) status="Saved "+path;
            }
            continue;
        }
        if (button==3 || key=='m') {
            std::string candidate;
            if (gplkSpectrumPrompt("Covariance model file",candidate)) {
                if (method==GPLK_SPECTRUM_CHOLESKY) {
                    GplkSpectrum trial; GplkSpectrumView trialView;
                    if (gplkBuildSpectrum(p,mode,method,candidate,trial,status) && trialView.initialize(trial)) {
                        result=trial; view=trialView; model=candidate; status="Covariance model updated";
                    }
                } else { model=candidate; status="Covariance model selected; choose Cholesky"; }
                gplkSpectrumModel=model;
            }
            continue;
        }
        if (button>=0 && button<3) {
            const GplkSpectrumMethod candidate=static_cast<GplkSpectrumMethod>(button);
            if (candidate!=method) {
                GplkSpectrum trial; GplkSpectrumView trialView;
                if (gplkBuildSpectrum(p,mode,candidate,model,trial,status) && trialView.initialize(trial)) {
                    result=trial; view=trialView; method=candidate; gplkSpectrumMethod=method;
                    zoom=false; status=std::string("Method: ")+gplkSpectrumMethodName(method);
                }
            }
            continue;
        }
        if (key=='A' && zoom && plot.contains(nx,ny)) {
            float bx=mx,by=my;
            if (!gplkReadBand(GIZA_BAND_RECTANGLE,0,mx,my,&bx,&by,&key)) break;
            float w,h; giza_get_paper_size_float(GIZA_UNITS_PIXELS,&w,&h);
            if (std::fabs(w-width)>.5f || std::fabs(h-height)>.5f) { sync=true; status="Resize cancelled zoom selection"; }
            else if (key=='A') { view.rectangle(mx,my,bx,by); status.clear(); }
            else if (key=='q') close=true;
            zoom=false;
        }
    }
    giza_close_device(); gplkResponsiveFonts=parentResponsiveFonts;
    if (previous>0) giza_select_device(previous);
    return true;
}

#endif
