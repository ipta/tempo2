// Headless integration checks using the real native Giza library.
#include "../plugin/gplk_giza.h"
#include "../plugin/gplk_spectrum.h"
#include "../plugin/gplk_overlays.h"
#include "../plugin/gplk_ui.h"
#include <memory>
#include <gtest/gtest.h>
#include <src/gtest_main.cc>

static bool near(float a, float b) { return std::fabs(a - b) < 1e-4f; }

#define CHECK(condition) do { if (!(condition)) { \
    std::fprintf(stderr, "testGplkGiza:%d: %s\n", __LINE__, #condition); \
    giza_close_devices(); return 1; } } while (0)

static int checkNativeGiza()
{
    CHECK(gplkOpenDevice("testGplkGiza.svg/svg", 1, 1));
    CHECK(!giza_device_has_cursor());
    gplkSetAspect(0.8f);
    float width, height;
    giza_get_paper_size_float(GIZA_UNITS_INCHES, &width, &height);
    CHECK(width > 0 && std::fabs(height / width - 0.8f) < 0.01f);
    gplkSetAspect(0); // device dimensions must remain valid

    CHECK(gplkSetNamedColour(0, "white"));
    CHECK(gplkSetNamedColour(1, "black"));
    CHECK(gplkSetNamedColour(2, "#123456"));
    float r, g, b;
    giza_get_colour_representation_float(2, &r, &g, &b);
    CHECK(near(r, 18.0f/255) && near(g, 52.0f/255) && near(b, 86.0f/255));
    CHECK(!gplkSetNamedColour(2, "not-a-colour"));
    giza_get_colour_representation_float(2, &r, &g, &b);
    CHECK(near(r, 18.0f/255)); // invalid names must preserve the existing colour

    gplkSetFont(3);
    giza_begin_buffer();
    giza_draw_background();
    giza_set_viewport_float(0.2f, 0.9f, 0.25f, 0.8f);
    giza_set_window_float(-10, 30, -2, 6);
    float x1, x2, y1, y2, sx1, sx2, sy1, sy2;
    gplkGetViewportPixels(&x1, &x2, &y1, &y2);
    gplkGetSurfaceSize(GIZA_UNITS_PIXELS, &sx1, &sx2, &sy1, &sy2);
    CHECK(near(x1 / sx2, 0.2f) && near(x2 / sx2, 0.9f));
    CHECK(near(y1 / sy2, 0.25f) && near(y2 / sy2, 0.8f));
    giza_get_window_float(&x1, &x2, &y1, &y2);
    CHECK(near(x1, -10) && near(x2, 30) && near(y1, -2) && near(y2, 6));
    float chx, chy;
    giza_get_character_size_float(GIZA_UNITS_WORLD, &chx, &chy);
    CHECK(chx > 0 && chy > 0);
    const float literalWidth=gplkLiteralText(0,0,"backend_01",false);
    CHECK(literalWidth>0 && literalWidth<15); // must fit this 40-world-unit window
    giza_box_float("BCNST", 0, 0, "BCNST", 0, 0);
    // The secondary axis needs explicit tick options in native Giza.
    giza_axis_float("NTS", 30, -2, 30, 6, -0.2f, 0.6f, 0, 0, 0.5f, 0, 0.5f, 0, 0);
    giza_label("Epoch", "Residual (\\gms)", "Native Giza gplk");
    const float x[] = {0, 10, 20}, y[] = {0, 3, 1};
    const float lo[] = {-0.2f, 2.8f, 0.8f}, hi[] = {0.2f, 3.2f, 1.2f};
    giza_set_colour_index(2);
    gplkDrawPoints(3, x, y, 16);
    int buffering=0;
    giza_get_buffering(&buffering);
    CHECK(buffering!=0); // point batches must not prematurely present a frame
    giza_error_bars_vert_float(3, x, lo, hi, 1);
    giza_end_buffer();
    giza_close_devices();

    FILE *svg = std::fopen("testGplkGiza.svg", "rb");
    CHECK(svg != NULL);
    std::fseek(svg, 0, SEEK_END);
    const long size = std::ftell(svg);
    std::fclose(svg);
    CHECK(size > 1000);
    GplkSpectrum dummy; std::string saveError;
    CHECK(!gplkWriteSpectrum(dummy,GPLK_SPECTRUM_GLS,2,"testGplkGiza.svg",saveError));
    svg=std::fopen("testGplkGiza.svg","rb"); CHECK(svg!=NULL);
    std::fseek(svg,0,SEEK_END); CHECK(std::ftell(svg)==size); std::fclose(svg);
    std::remove("testGplkGiza.svg");

    CHECK(gplkOpenDevice("/null", 2, 2));
    int column, row;
    giza_get_panel(&column, &row);
    CHECK(column == 1 && row == 1);
    float cx = 7, cy = 9;
    char key = 'x';
    CHECK(!gplkReadBand(GIZA_BAND_NONE, 0, 0, 0, &cx, &cy, &key));
    CHECK(key == 'q' && cx == 7 && cy == 9);
    giza_close_devices();
    return 0;
}

TEST(testGplkGiza, nativeRenderer)
{
    EXPECT_EQ(checkNativeGiza(), 0);
}

TEST(testGplkGiza, responsiveLayoutAndLargeFlagGroups)
{
    GplkLayout small,wide,tall;
    small.update(700,620,35); wide.update(1500,620,35); tall.update(700,1200,35);
    EXPECT_GT(wide.columns,small.columns);
    EXPECT_GT(tall.visibleRows,small.visibleRows);
    EXPECT_GT((wide.plot.right-wide.plot.left)*wide.width,
              (small.plot.right-small.plot.left)*small.width);
    EXPECT_GT((tall.plot.top-tall.plot.bottom)*tall.height,
              (small.plot.top-small.plot.bottom)*small.height);
    EXPECT_LT(small.jumps.top,small.bottom.top);
    EXPECT_LT(small.bottom.top,small.plot.bottom);
    EXPECT_LT(small.plot.top,small.top.bottom);
    EXPECT_LT(small.plot.right,small.colours.left);
    EXPECT_GE(small.colours.bottom,small.bottom.top);

    std::auto_ptr<pulsar> p(new pulsar());
    std::vector<observation> observations(151);
    p->nobs=observations.size(); p->obsn=&observations[0];
    for (int i=0;i<150;++i) {
        observations[i].nFlags=2;
        std::strcpy(observations[i].flagID[0],"-backend");
        std::sprintf(observations[i].flagVal[0],"backend_%d",i);
        std::strcpy(observations[i].flagID[1],"-selectAll");
        std::strcpy(observations[i].flagVal[1],"off");
    }
    GplkColourGroups colours; colours.mode="-backend"; colours.rebuild(*p);
    ASSERT_EQ(colours.values.size(),151u); // includes missing-flag observations
    EXPECT_EQ(colours.values.back(),"(missing)");
    EXPECT_NE(colours.group[0],colours.group[149]);
    EXPECT_EQ(colours.group[150],150);
    const int ids[]={149,0,150};
    colours.groupPositions(3,ids);
    ASSERT_EQ(colours.plotGroups[149].size(),1u);
    EXPECT_EQ(colours.plotGroups[149][0],0); // display positions differ from observation IDs
    EXPECT_EQ(colours.plotGroups[0][0],1);
    EXPECT_EQ(colours.plotGroups[150][0],2);
    EXPECT_EQ(std::count(colours.choices.begin(),colours.choices.end(),"-backend"),1);
    EXPECT_EQ(std::count(colours.choices.begin(),colours.choices.end(),"-selectAll"),0);
    for (unsigned int i=0;i<colours.values.size();++i) {
        EXPECT_GE(GplkColourGroups::colour(i),2);
        EXPECT_LE(GplkColourGroups::colour(i),12);
    }
    p->nJumps=1;
    std::strcpy(p->jumpStr[1],"-backend backend_0");
    GplkUiState ui; ui.rebuild(*p,colours);
    ASSERT_EQ(ui.jumps.size(),1u); // jump buttons exclude color legend entries
    ASSERT_EQ(ui.legend.size(),151u);
    gplkLayout.update(1200,800,ui.jumps.size());
    const GplkRect selector=gplkColourSelector();
    char key='A'; bool grow=false;
    EXPECT_TRUE(gplkHandleUi((selector.left+selector.right)/2,(selector.bottom+selector.top)/2,
                            key,*p,colours,ui,grow));
    EXPECT_TRUE(ui.dropdown);
    EXPECT_EQ(p->fitJump[1],0); // right-panel clicks cannot toggle jumps
    ui.dropdown=false;
    const GplkRect jump=gplkFooterCell(0);
    key='A';
    EXPECT_TRUE(gplkHandleUi((jump.left+jump.right)/2,
                            (jump.bottom+jump.top)/2*gplkLayout.bottom.top,
                            key,*p,colours,ui,grow));
    EXPECT_EQ(p->fitJump[1],1);
    EXPECT_EQ(colours.mode,"-backend"); // bottom clicks cannot change grouping
}

TEST(testGplkGiza, fixedInstrumentCorrectionsAndJumpLabels)
{
    std::auto_ptr<pulsar> p(new pulsar());
    const char *lines[]={"-sys UWL_PDFB4 0 1", "-medusa_58925_jump 1 -2e-7 1",
                         "-medusa_58925_jump 2 -2e-7 1", "-sys 1 -2e-7 1",
                         "-group_jump 1 -2e-7 1", "-medusa_notmjd_jump 1 -2e-7 1",
                         "-pdfb4_59000_jump 1 3e-7 0", "-medusa_5892_jump 1 0 1"};
    p->nJumps=sizeof(lines)/sizeof(lines[0]);
    for (int i=1;i<=p->nJumps;++i) {
        std::strcpy(p->jumpStr[i],lines[i-1]); p->fitJump[i]=1;
    }
    p->jumpVal[2]=-2e-7; p->jumpValErr[2]=1e-8; p->jumpSAT[2]=1;
    EXPECT_EQ(gplkDescribeJump(lines[0]).label,"UWL_PDFB4");
    EXPECT_EQ(gplkDescribeJump(lines[1]).label,"medusa_58925_jump");
    EXPECT_EQ(gplkDescribeJump("-blahblah 1 8989.9898").label,"blahblah");
    EXPECT_EQ(gplkDescribeJump("-blahblah 1 8989.9898 0").label,"blahblah");
    EXPECT_EQ(gplkDescribeJump("-arbitrary_name 1 -2e-7 1").label,"arbitrary_name");
    EXPECT_EQ(gplkDescribeJump("-sys 1 -2e-7 1").label,"1");
    EXPECT_EQ(gplkDescribeJump("-g 1 -2e-7 1").label,"1");
    EXPECT_FALSE(gplkDescribeJump("-blahblah 1 8989.9898").correction);
    EXPECT_EQ(gplkDescribeJump("MJD 58000 59000 0 1").label,"MJD 58000..59000");
    EXPECT_TRUE(gplkDescribeJump(lines[1]).correction);
    EXPECT_TRUE(gplkDescribeJump(lines[6]).correction);
    for (int i=2;i<6;++i) EXPECT_FALSE(gplkDescribeJump(lines[i]).correction);
    EXPECT_FALSE(gplkDescribeJump(lines[7]).correction);
    EXPECT_EQ(gplkKeepCorrectionsFixed(p.get(),1),2);
    EXPECT_EQ(gplkKeepCorrectionsFixed(p.get(),1),0); // repeated fits cannot re-enable corrections
    EXPECT_EQ(p->fitJump[2],0); EXPECT_EQ(p->fitJump[7],0);
    EXPECT_EQ(p->fitJump[1],1); EXPECT_EQ(p->fitJump[4],1);
    EXPECT_EQ(p->jumpVal[2],-2e-7); EXPECT_EQ(p->jumpValErr[2],1e-8);
    EXPECT_EQ(p->jumpSAT[2],1); EXPECT_STREQ(p->jumpStr[2],lines[1]);
    GplkColourGroups colours; GplkUiState ui; ui.rebuild(*p,colours);
    ASSERT_EQ(ui.jumps.size(),unsigned(p->nJumps));
    EXPECT_TRUE(ui.jumps[1].fixed); EXPECT_FALSE(ui.jumps[1].active);
    EXPECT_FALSE(ui.jumps[0].fixed); EXPECT_TRUE(ui.jumps[0].active);
    gplkLayout.update(1200,800,ui.jumps.size());
    const GplkRect correction=gplkFooterCell(1);
    char key='A'; bool grow=false;
    EXPECT_TRUE(gplkHandleUi((correction.left+correction.right)/2,
                            (correction.bottom+correction.top)/2*gplkLayout.bottom.top,
                            key,*p,colours,ui,grow));
    EXPECT_EQ(p->fitJump[2],0); EXPECT_EQ(p->jumpVal[2],-2e-7);
    EXPECT_EQ(p->jumpValErr[2],1e-8);
    EXPECT_TRUE(gplkToggleJumpFit(*p,1)); EXPECT_EQ(p->fitJump[1],0);
    EXPECT_TRUE(gplkToggleJumpFit(*p,1)); EXPECT_EQ(p->fitJump[1],1);
    EXPECT_FALSE(gplkToggleJumpFit(*p,0)); EXPECT_FALSE(gplkToggleJumpFit(*p,p->nJumps+1));
    // A file reload can restore the trailing fit=1; the next fit locks it again.
    p->fitJump[2]=1;
    EXPECT_EQ(gplkKeepCorrectionsFixed(p.get(),1),1);
    EXPECT_EQ(p->jumpVal[2],-2e-7);
}

TEST(testGplkGiza, noiseDelaysAndUncertainties)
{
    std::auto_ptr<pulsar> p(new pulsar());
    observation o[2] = {};
    p->nobs=2; p->obsn=o;
    o[0].TNRedSignal=1e-6; o[0].TNRedErr=0.2e-6;
    o[0].TNDMSignal=2e-6; o[0].TNDMErr=0.3e-6;
    o[0].TNChromSignal=-0.5e-6; o[0].TNChromErr=0.1e-6;
    o[0].residual=4e-6; o[0].toaErr=2;
    // Absent model fields can contain arbitrary values; never read them.
    EXPECT_FALSE(gplkNoiseAvailable(*p,GPLK_NOISE_TOTAL));
    EXPECT_EQ(gplkNoiseValue(*p,o[0],GPLK_NOISE_TOTAL),0);
    p->TNRedAmp=-13; p->TNRedGam=4;
    p->TNDMAmp=-14; p->TNDMGam=3;
    p->TNChromAmp=-14; p->TNChromGam=3; p->TNChromIdx=4;
    EXPECT_TRUE(gplkNoiseAvailable(*p,GPLK_NOISE_RED));
    EXPECT_NEAR(gplkNoiseValue(*p,o[0],GPLK_NOISE_TOTAL),2.5e-6,1e-18);
    EXPECT_NEAR(gplkNoiseValue(*p,o[0],GPLK_NOISE_CLEAN),1.5e-6,1e-18);
    EXPECT_NEAR(gplkNoiseUncertainty(*p,o[0],GPLK_NOISE_DM),0.3e-6,1e-18);
    EXPECT_EQ(gplkNoiseUncertainty(*p,o[0],GPLK_NOISE_TOTAL),0);
    EXPECT_NEAR(gplkNoiseUncertainty(*p,o[0],GPLK_NOISE_CLEAN),2e-6,1e-18);
    o[0].deleted=1;
    EXPECT_FALSE(gplkNoiseAvailable(*p,GPLK_NOISE_TOTAL));
    EXPECT_NEAR(static_cast<double>(o[0].residual),4e-6,1e-18);
}

TEST(testGplkGiza, irregularWeightedSpectrum)
{
    std::vector<GplkSpectrumSample> data;
    const double frequency=0.037;
    for (int i=0; i<100; ++i) {
        double t=i*3.2+0.8*std::sin(i*0.73);
        data.push_back(GplkSpectrumSample(t,5+2*std::sin(2*M_PI*frequency*t)+
                                         0.3*std::cos(2*M_PI*frequency*t),0.5+(i%7)*0.1));
    }
    // The large-error outlier should have negligible influence.
    data.push_back(GplkSpectrumSample(120.1,1000,1e6));
    GplkSpectrum result;
    std::string error;
    ASSERT_TRUE(gplkComputeSpectrum(data,0.01,0.1,result,error)) << error;
    EXPECT_NEAR(result.frequency[result.peak],frequency,0.0004);
    EXPECT_GT(result.power[result.peak],0.98);
    // Independent NumPy weighted least-squares reference (intercept + sine + cosine).
    EXPECT_NEAR(result.power[0],0.007393132645637457,1e-12);
    EXPECT_NEAR(result.power[50],0.0012407580095477044,1e-12);
    EXPECT_NEAR(result.power[130],0.0006657530375474563,1e-12);
    EXPECT_NEAR(result.power[255],0.00016638590606610215,1e-12);
    for (unsigned int i=0; i<result.power.size(); ++i) {
        EXPECT_GE(result.power[i],0); EXPECT_LE(result.power[i],1);
        EXPECT_GE(result.window[i],0); EXPECT_LE(result.window[i],1+1e-12);
    }
    const std::vector<double> powers=result.power;
    for (unsigned int i=0; i<data.size(); ++i) {
        data[i].day+=50000; data[i].value+=100; data[i].sigma*=20;
    }
    ASSERT_TRUE(gplkComputeSpectrum(data,0.01,0.1,result,error));
    for (unsigned int i=0; i<powers.size(); ++i) EXPECT_NEAR(powers[i],result.power[i],1e-10);
    EXPECT_FALSE(gplkComputeSpectrum(data,0.1,0.01,result,error));
    for (unsigned int i=0; i<data.size(); ++i) data[i].value=1;
    EXPECT_FALSE(gplkComputeSpectrum(data,0.01,0.1,result,error));
    for (unsigned int i=0; i<data.size(); ++i) data[i].day=1;
    EXPECT_FALSE(gplkComputeSpectrum(data,0.01,0.1,result,error));
}

TEST(testGplkGiza, spectrumZoomPreservesSamples)
{
    GplkSpectrum s;
    s.frequency.push_back(.01); s.frequency.push_back(.1);
    s.power.push_back(.2); s.power.push_back(.8);
    const std::vector<double> original=s.power;
    GplkSpectrumView view;
    ASSERT_TRUE(view.initialize(s));
    const double lo=view.x0,hi=view.x1;
    EXPECT_TRUE(view.rectangle(hi-.1,.9,lo+.1,.1)); // reverse drag
    EXPECT_NEAR(view.x0,lo+.1,1e-12); EXPECT_NEAR(view.x1,hi-.1,1e-12);
    EXPECT_FALSE(view.rectangle(lo,.5,lo,.5));
    view.scale(.5); EXPECT_LT(view.x1-view.x0,hi-lo);
    view.reset(); EXPECT_EQ(view.x0,lo); EXPECT_EQ(view.x1,hi);
    EXPECT_EQ(s.power,original);
    s.psd=true; s.power[0]=1e-30; s.power[1]=1e-25;
    ASSERT_TRUE(view.initialize(s));
    EXPECT_LT(view.y0,-30); EXPECT_GT(view.y1,-25);
}

TEST(testGplkGiza, nativeSpectraAndCovarianceLeaveFitUnchanged)
{
    std::auto_ptr<pulsar> p(new pulsar());
    std::strcpy(p->name,"NativeSpectrumTest"); p->detL=7;
    const int n=48;
    std::vector<observation> observations(n);
    p->obsn=&observations[0]; p->nobs=n;
    std::vector<double> times(n);
    for (int i=0;i<n;++i) times[i]=3*i+.13*std::sin(i*.7);
    const double span=times.back()-times.front(),fundamental=1/(span*n/(n-1.0));
    for (int i=0;i<n;++i) {
        const double t=times[n-1-i]; // native input sorting must retain error/ID alignment
        observation &o=observations[i]; o.sat=52601+t; o.bat=o.sat;
        o.residual=3e-6*std::cos(2*M_PI*4*fundamental*t)-1e-6*std::sin(2*M_PI*4*fundamental*t);
        o.toaErr=.5+.1*(i%7);
    }
    const longdouble firstResidual=observations[0].residual,firstDay=observations[0].sat;
    GplkSpectrum chol,weighted; std::string error;
    ASSERT_FALSE(gplkComputeNativeSpectrum(*p,2,GPLK_SPECTRUM_CHOLESKY,"/no/such/gplk-model",0,0,chol,error));
    ASSERT_TRUE(gplkComputeNativeSpectrum(*p,2,GPLK_SPECTRUM_CHOLESKY,"MODEL T2",0,0,chol,error)) << error;
    ASSERT_EQ(chol.frequency.size(),unsigned(n/2-1)); EXPECT_TRUE(chol.psd);
    EXPECT_NEAR(chol.frequency[3],4*fundamental,1e-13);
    // A pure cosine/sine at a native Fourier frequency has known coefficients,
    // irrespective of the heteroscedastic weights. TEMPO2's PSD units are yr^3.
    const double expected=span/365.25/2/std::pow(365.25*86400,2)*10e-12;
    EXPECT_NEAR(chol.power[3]/expected,1,1e-9);
    EXPECT_EQ(p->detL,7); EXPECT_EQ(observations[0].residual,firstResidual);
    EXPECT_EQ(observations[0].sat,firstDay);
    GplkSpectrum red;
    ASSERT_TRUE(gplkComputeNativeSpectrum(*p,2,GPLK_SPECTRUM_CHOLESKY,
                                         "MODEL T2\nMODEL T2PowerLaw 4 1e-23 0.5",0,0,red,error));
    EXPECT_NEAR(red.power[3]/expected,1,1e-8);
    EXPECT_NE(red.power[2],chol.power[2]); // a nonzero covariance changes off-frequency weighting
    ASSERT_TRUE(gplkComputeNativeSpectrum(*p,2,GPLK_SPECTRUM_TEMPO2,"",0,0,weighted,error)) << error;
    std::vector<double> x(n),y(n),e(n),fx(n),fy(n),re(n),im(n);
    for (int i=0;i<n;++i) {
        const observation &o=observations[n-1-i]; x[i]=double(o.sat)-52601;
        y[i]=double(o.residual); e[i]=o.toaErr*1e-6;
    }
    int bins=0;
    TKspectrum(&x[0],&y[0],&e[0],n,0,0,0,0,0,4,1,1,1,&fx[0],&fy[0],&bins,0,0,&re[0],&im[0]);
    ASSERT_EQ(weighted.power.size(),unsigned(bins));
    for (int i=0;i<bins;++i) {
        EXPECT_EQ(weighted.frequency[i],fx[i]); EXPECT_EQ(weighted.power[i],fy[i]);
    }
    ASSERT_TRUE(gplkComputeNativeSpectrum(*p,2,GPLK_SPECTRUM_CHOLESKY,"MODEL T2",1.9*fundamental,6.1*fundamental,chol,error));
    EXPECT_EQ(chol.frequency.size(),5u);
    observations[0].deleted=1; observations[1].toaErr=0;
    ASSERT_TRUE(gplkComputeNativeSpectrum(*p,2,GPLK_SPECTRUM_TEMPO2,"",0,0,weighted,error));
    EXPECT_EQ(weighted.samples,unsigned(n-2));
    int isSet=1,fit=1; longdouble start=52601+times[10];
    p->param[param_start].aSize=1; p->param[param_start].paramSet=&isSet;
    p->param[param_start].fitFlag=&fit; p->param[param_start].val=&start;
    ASSERT_TRUE(gplkComputeNativeSpectrum(*p,2,GPLK_SPECTRUM_TEMPO2,"",0,0,weighted,error));
    EXPECT_EQ(weighted.samples,unsigned(n-12));
}

TEST(testGplkGiza, solarGeometryAndOverlayUnits)
{
    observation o={};
    o.delayCorr=1;
    o.psrPos[0]=2; // Non-unit vectors are normalized by the geometry helper.
    o.earth_ssb[0]=1; o.observatory_earth[0]=1;
    o.sun_ssb[0]=10;
    double degrees=-1;
    ASSERT_TRUE(gplkSolarSeparation(o,degrees));
    EXPECT_NEAR(degrees,0,1e-12);
    o.sun_ssb[0]=2; o.sun_ssb[1]=8;
    ASSERT_TRUE(gplkSolarSeparation(o,degrees));
    EXPECT_NEAR(degrees,90,1e-12);
    o.sun_ssb[0]=-8; o.sun_ssb[1]=0;
    ASSERT_TRUE(gplkSolarSeparation(o,degrees));
    EXPECT_NEAR(degrees,180,1e-12);
    // Spacecraft ephemeris uses observatory_earth as its SSB position, matching
    // the STL_FBAT branch in shapiro_delay.C; do not subtract Earth twice.
    std::strcpy(o.telID,"STL_FBAT");
    o.earth_ssb[0]=100; o.sun_ssb[0]=1; o.sun_ssb[1]=8;
    ASSERT_TRUE(gplkSolarSeparation(o,degrees));
    EXPECT_NEAR(degrees,90,1e-12);
    o.psrPos[0]=0;
    EXPECT_FALSE(gplkSolarSeparation(o,degrees));
    o.psrPos[0]=1; o.delayCorr=0;
    EXPECT_FALSE(gplkSolarSeparation(o,degrees));
    EXPECT_NEAR(gplkOverlayDelay(4e-6,1e-6,-1,10),4,1e-12);
    EXPECT_NEAR(gplkOverlayDelay(4e-6,1e-3,-1,10),0.004,1e-12);
    EXPECT_NEAR(gplkOverlayDelay(4e-6,1,1,10),4e-5,1e-15);
}
