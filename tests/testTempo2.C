#include <gtest/gtest.h>
#include <src/gtest_main.cc>
#include "tempo2.h"
#include <string>
#include <unistd.h>
#ifndef DATDIR
#define DATDIR .
#endif

#ifdef LONGDOUBLE_IS_FLOAT128
#include <quadmath.h>
#endif

namespace {

const char *ELL1_BASE =
    "PSRJ J1234+5678\n"
    "F0 1\n"
    "PEPOCH 58000\n"
    "BINARY ELL1\n"
    "A1 1\n"
    "TASC 58000\n"
    "EPS1 0\n"
    "EPS2 0\n";

const char *DDGR_BASE =
    "PSRJ J1234+5678\n"
    "F0 1\n"
    "PEPOCH 58000\n"
    "BINARY DDGR\n"
    "PB 0.5\n"
    "A1 1\n"
    "T0 58000\n"
    "ECC 0.1\n"
    "OM 45\n"
    "M2 0.3\n"
    "MTOT 1.7\n";

void readParText(pulsar *psr, const std::string &text)
{
    FILE *fin = tmpfile();
    ASSERT_TRUE(fin != NULL);
    ASSERT_NE(fputs(text.c_str(),fin),EOF);
    rewind(fin);
    psr->noWarnings = 2;
    readSimpleParfile(fin,psr);
    fclose(fin);
}

double ell1Delay(const std::string &extra, double dt)
{
    pulsar psr;
    MAX_PSR = 1;
    initialise(&psr,0);
    readParText(&psr,std::string(ELL1_BASE)+extra);
    psr.obsn[0].bbat =
            longdouble(58000.0)+static_cast<longdouble>(dt)/SECDAYl;
    const double delay = ELL1model(&psr,0,0,-1,0);
    destroyOne(&psr);
    return delay;
}

} // namespace

TEST(testTempo2h, maxpsrset){
    ASSERT_GT(MAX_PSR,0);
}

TEST(testLongDouble, precision){
    EXPECT_EQ(sizeof(longdouble),16);
#ifdef LONGDOUBLE_IS_FLOAT128
    EXPECT_GT(FLT128_DIG,32);
#else
    EXPECT_GT(LDBL_DIG,17);
#endif
}

TEST(testLongDouble, printAndParse){
    longdouble ld = longdouble(123.0);
    char sb[1024];
    ld_sprintf(sb,"%.1Lf",ld);
    ASSERT_STREQ(sb,"123.0");
    ld = longdouble(50000.12345678912345);
    ld_sprintf(sb,"%.14Lf",ld);

    ASSERT_STREQ("50000.12345678912345",sb);

    longdouble ld2 = parse_longdouble(sb);
    ASSERT_EQ(ld,ld2);

    unsigned uu = 2147483649;
    long long int ll = 2147483649L;
    long long int ll2 = 17179869185L;

    ld_sprintf(sb,"%% %d %.1f %u %.1lf %lld %lld %s %.1Lf %% %c",-1,0.1,uu,0.1,ll,ll2,"t",ld,'x');
    ASSERT_STREQ("% -1 0.1 2147483649 0.1 2147483649 17179869185 t 50000.1 % x",sb);
}

TEST(testELL1BinaryPhase, pbFb2MatchesExplicitZeroFb1)
{
    const double dt = 1.0e6;
    const double sparse = ell1Delay(
            "PB 0.5\n"
            "FB2 1e-18\n",dt);
    const double canonical = ell1Delay(
            "PB 0.5\n"
            "FB1 0\n"
            "FB2 1e-18\n",dt);
    EXPECT_NEAR(sparse,canonical,1e-14);
}

TEST(testELL1BinaryPhase, sparseFb2ChangesDelay)
{
    const double dt = 1.0e6;
    const double withFb2 = ell1Delay(
            "PB 0.5\n"
            "FB2 1e-18\n",dt);
    const double withoutFb2 = ell1Delay(
            "PB 0.5\n"
            "FB2 0\n",dt);
    EXPECT_GT(fabs(withFb2-withoutFb2),1e-6);
}

TEST(testELL1BinaryPhase, pbdotSuppliesMissingFb1)
{
    const double dt = 1.0e6;
    const double pbSeconds = 0.5*SECDAY;
    const double pbdot = 1.0e-8;
    const double fb1 = -pbdot/(pbSeconds*pbSeconds);
    // Nested scopes keep only one stack pulsar live at a time:
    // sizeof(pulsar) is ~5.4 MiB and the default stack is only 8 MiB.
    const longdouble bbat =
            longdouble(58000.0)+static_cast<longdouble>(dt)/SECDAYl;
    double hybridDelay = 0.0;
    double canonicalDelay = 0.0;

    {
        pulsar hybrid;
        MAX_PSR = 1;
        initialise(&hybrid,0);
        readParText(&hybrid,std::string(ELL1_BASE)+
                "PB 0.5\n"
                "PBDOT 1e-8\n"
                "FB2 1e-18\n");
        hybrid.obsn[0].bbat = bbat;
        hybridDelay = ELL1model(&hybrid,0,0,-1,0);
        destroyOne(&hybrid);
    }
    {
        pulsar canonical;
        MAX_PSR = 1;
        initialise(&canonical,0);
        readParText(&canonical,std::string(ELL1_BASE)+
                "PB 0.5\n"
                "FB1 0\n"
                "FB2 1e-18\n");
        canonical.param[param_fb].val[1] = fb1;
        canonical.obsn[0].bbat = bbat;
        canonicalDelay = ELL1model(&canonical,0,0,-1,0);
        destroyOne(&canonical);
    }
    EXPECT_NEAR(hybridDelay,canonicalDelay,1e-14);
}

TEST(testELL1BinaryPhase, fb0PbdotSuppliesMissingFb1)
{
    const double dt = 1.0e6;
    const double fb0 = 2.3148148148148148e-5;
    const double pbdot = 1.0e-8;
    const double fb1 = -pbdot*fb0*fb0;
    const longdouble bbat =
            longdouble(58000.0)+static_cast<longdouble>(dt)/SECDAYl;
    double hybridDelay = 0.0;
    double canonicalDelay = 0.0;

    {
        pulsar hybrid;
        MAX_PSR = 1;
        initialise(&hybrid,0);
        readParText(&hybrid,std::string(ELL1_BASE)+
                "FB0 2.3148148148148148e-5\n"
                "PBDOT 1e-8\n"
                "FB2 1e-18\n");
        hybrid.obsn[0].bbat = bbat;
        hybridDelay = ELL1model(&hybrid,0,0,-1,0);
        destroyOne(&hybrid);
    }
    {
        pulsar canonical;
        MAX_PSR = 1;
        initialise(&canonical,0);
        readParText(&canonical,std::string(ELL1_BASE)+
                "FB0 2.3148148148148148e-5\n"
                "FB1 0\n"
                "FB2 1e-18\n");
        canonical.param[param_fb].val[1] = fb1;
        canonical.obsn[0].bbat = bbat;
        canonicalDelay = ELL1model(&canonical,0,0,-1,0);
        destroyOne(&canonical);
    }
    EXPECT_NEAR(hybridDelay,canonicalDelay,1e-14);
}

TEST(testELL1BinaryPhase, missingInteriorFb2IsZero)
{
    const double dt = 1.0e5;
    const double sparse = ell1Delay(
            "PB 0.5\n"
            "FB1 -1e-12\n"
            "FB3 1e-20\n",dt);
    const double canonical = ell1Delay(
            "PB 0.5\n"
            "FB1 -1e-12\n"
            "FB2 0\n"
            "FB3 1e-20\n",dt);
    EXPECT_NEAR(sparse,canonical,1e-14);
}

TEST(testELL1BinaryPhase, fb0Fb2MatchesEquivalentPbSeries)
{
    const double dt = 1.0e6;
    const double fromFb0 = ell1Delay(
            "FB0 2.3148148148148148e-5\n"
            "FB2 1e-18\n",dt);
    const double fromPb = ell1Delay(
            "PB 0.5\n"
            "FB1 0\n"
            "FB2 1e-18\n",dt);
    EXPECT_NEAR(fromFb0,fromPb,1e-14);
}

TEST(testELL1BinaryPhase, ordinaryPbdotPathIsUnchanged)
{
    pulsar psr;
    MAX_PSR = 1;
    initialise(&psr,0);
    readParText(&psr,std::string(ELL1_BASE)+
            "PB 0.5\n"
            "PBDOT 1e-8\n");

    const double requestedDt = 1.0e6;
    const double pb = 0.5*SECDAY;

    psr.obsn[0].bbat = longdouble(58000.0)+
            static_cast<longdouble>(requestedDt)/SECDAYl;
    const double effectiveDt =
            (static_cast<double>(psr.obsn[0].bbat)-
             static_cast<double>(psr.param[param_tasc].val[0]))*SECDAY;
    const double expectedOrbits =
            effectiveDt/pb-
            0.5e-8*(effectiveDt/pb)*(effectiveDt/pb);
    const double phase = 2.0*M_PI*
            (expectedOrbits-floor(expectedOrbits));
    const double dre = sin(phase);
    const double drep = cos(phase);
    const double drepp = -sin(phase);
    const double an = 2.0*M_PI/pb;
    const double expectedDelay =
            -dre*(1.0-an*drep+(an*drep)*(an*drep)
                    +0.5*an*an*dre*drepp);

    EXPECT_NEAR(ELL1model(&psr,0,0,-1,0),expectedDelay,1e-14);
    destroyOne(&psr);
}

TEST(testBinaryParfile, explicitFb1DisablesPbdot)
{
    pulsar psr;
    MAX_PSR = 1;
    initialise(&psr,0);
    readParText(&psr,std::string(ELL1_BASE)+
            "PB 0.5\n"
            "PBDOT 1e-12 1\n"
            "FB1 -1e-18\n");

    EXPECT_EQ(psr.param[param_fb].paramSet[1],1);
    EXPECT_EQ(psr.param[param_pbdot].paramSet[0],0);
    EXPECT_EQ(psr.param[param_pbdot].fitFlag[0],0);
    destroyOne(&psr);
}

TEST(testBinaryParfile, pbdotRemainsActiveWithoutExplicitFb1)
{
    pulsar psr;
    MAX_PSR = 1;
    initialise(&psr,0);
    readParText(&psr,std::string(ELL1_BASE)+
            "PB 0.5\n"
            "PBDOT 1e-12 1\n"
            "FB2 1e-18\n");

    EXPECT_EQ(psr.param[param_fb].paramSet[1],0);
    EXPECT_EQ(psr.param[param_pbdot].paramSet[0],1);
    EXPECT_EQ(psr.param[param_pbdot].fitFlag[0],1);
    destroyOne(&psr);
}

TEST(testELL1BinaryPhase, explicitFb1HasZeroPbdotDerivative)
{
    pulsar psr;
    MAX_PSR = 1;
    initialise(&psr,0);
    readParText(&psr,std::string(ELL1_BASE)+"PB 0.5\n");

    // Bypass post-parse consistency handling to exercise the defensive
    // derivative branch for direct programmatic construction.
    psr.param[param_pbdot].paramSet[0] = 1;
    psr.param[param_pbdot].val[0] = 1e-12;
    psr.param[param_fb].paramSet[1] = 1;
    psr.param[param_fb].val[1] = -1e-18;
    psr.obsn[0].bbat =
            longdouble(58000.0)+longdouble(1.0e5)/SECDAYl;

    EXPECT_EQ(ELL1model(&psr,0,0,param_pbdot,0),0.0);
    destroyOne(&psr);
}

TEST(testBinaryParfile, ddgrRejectsFb0)
{
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    EXPECT_EXIT(
        {
            ASSERT_EQ(dup2(STDERR_FILENO,STDOUT_FILENO),STDOUT_FILENO);
            pulsar psr;
            MAX_PSR = 1;
            initialise(&psr,0);
            readParText(&psr,std::string(DDGR_BASE)+"FB0 2e-5\n");
        },
        ::testing::ExitedWithCode(1),
        "BINARY DDGR does not support FB parameters");
}

TEST(testBinaryParfile, ddgrRejectsHigherFb)
{
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    EXPECT_EXIT(
        {
            ASSERT_EQ(dup2(STDERR_FILENO,STDOUT_FILENO),STDOUT_FILENO);
            pulsar psr;
            MAX_PSR = 1;
            initialise(&psr,0);
            readParText(&psr,std::string(DDGR_BASE)+"FB1 -1e-18\n");
        },
        ::testing::ExitedWithCode(1),
        "BINARY DDGR does not support FB parameters");
}

TEST(testBinaryParfile, ddgrWithoutFbParses)
{
    pulsar psr;
    MAX_PSR = 1;
    initialise(&psr,0);
    readParText(&psr,DDGR_BASE);
    EXPECT_STREQ(psr.binaryModel,"DDGR");
    EXPECT_EQ(psr.param[param_fb].paramSet[0],0);
    EXPECT_EQ(psr.param[param_fb].paramSet[1],0);
    destroyOne(&psr);
}

TEST(testFormResiduals, basicBATs){
    pulsar _psr;
    pulsar *psr = &_psr;
    MAX_PSR=1;
    char timFile[MAX_PSR][MAX_FILELEN],parFile[MAX_PSR][MAX_FILELEN];
    int npsr=1;
    initialise(psr,0); /* Initialise all */

    strcpy(parFile[0],DATDIR "/test1.par");
    strcpy(timFile[0],DATDIR "/test1.tim");

    readParfile(psr,parFile,timFile,npsr);
    readTimfile(psr,timFile,npsr);
    psr->noWarnings=2;
    formBatsAll(psr,npsr);
    formResiduals(psr,npsr,0);
    ASSERT_LT(static_cast<double>(fabsl(psr[0].obsn[1].residual)),1e-9l);
    ASSERT_LT(static_cast<double>(fabsl(psr[0].obsn[2].residual-longdouble(0.4))),1e-9l);
}




TEST(testFormResiduals, subtractBATs){
    pulsar _psr;
    pulsar *psr = &_psr;
    MAX_PSR=1;
    char timFile[MAX_PSR][MAX_FILELEN],parFile[MAX_PSR][MAX_FILELEN];
    int npsr=1;
    initialise(psr,0); /* Initialise all */

    strcpy(parFile[0],DATDIR "/test1.par");
    strcpy(timFile[0],DATDIR "/test1.tim");

    readParfile(psr,parFile,timFile,npsr);
    readTimfile(psr,timFile,npsr);

    for(int iobs = 0; iobs < psr->nobs; iobs++){
        psr[0].obsn[iobs].bat = psr[0].obsn[iobs].sat;
        psr[0].obsn[iobs].bbat = psr[0].obsn[iobs].sat;
        psr[0].obsn[iobs].delayCorr=0;
    }
    formResiduals(psr,npsr,0);
    for(int iobs = 0; iobs < psr->nobs; iobs++){
        psr[0].obsn[iobs].sat -= psr[0].obsn[iobs].residual/SECDAYl;
        psr[0].obsn[iobs].bat = psr[0].obsn[iobs].sat;
        psr[0].obsn[iobs].bbat = psr[0].obsn[iobs].sat;
    }
    formResiduals(psr,npsr,0);

    for(int iobs = 0; iobs < psr->nobs; iobs++){
        EXPECT_LT(static_cast<double>(fabsl(psr[0].obsn[iobs].residual)),TEST_DELTA) << "Precision lost in formResiduals (s)";
    }

    for(int iobs = 1; iobs < psr->nobs; iobs++){
        psr[0].obsn[iobs].sat += longdouble(4e-9)/SECDAYl;
        psr[0].obsn[iobs].bat = psr[0].obsn[iobs].sat;
        psr[0].obsn[iobs].bbat = psr[0].obsn[iobs].sat;
    }
    formResiduals(psr,npsr,0);

    for(int iobs = 1; iobs < psr->nobs; iobs++){
        EXPECT_LT(static_cast<double>(fabsl(psr[0].obsn[iobs].residual-longdouble(4.0e-9))),TEST_DELTA) << "Precision lost in formResiduals (s)";
    }
}


TEST(testFormResiduals, subtractSATs){
    pulsar _psr;
    pulsar *psr = &_psr;
    MAX_PSR=1;
    char timFile[MAX_PSR][MAX_FILELEN],parFile[MAX_PSR][MAX_FILELEN];
    int npsr=1;
    initialise(psr,0); /* Initialise all */

    strcpy(parFile[0],DATDIR "/test1.par");
    strcpy(timFile[0],DATDIR "/test2.tim");

    readParfile(psr,parFile,timFile,npsr);
    readTimfile(psr,timFile,npsr);
    psr->noWarnings=2;
    formBatsAll(psr,npsr);
    formResiduals(psr,npsr,0);
    for(int iobs = 0; iobs < psr->nobs; iobs++){
        psr[0].obsn[iobs].sat -= psr[0].obsn[iobs].residual/SECDAYl;
    }
    formBatsAll(psr,npsr);
    formResiduals(psr,npsr,0);
    for(int iobs = 0; iobs < psr->nobs; iobs++){
        psr[0].obsn[iobs].sat -= psr[0].obsn[iobs].residual/SECDAYl;
    }
    formBatsAll(psr,npsr);
    formResiduals(psr,npsr,0);
    for(int iobs = 0; iobs < psr->nobs; iobs++){
        psr[0].obsn[iobs].sat -= psr[0].obsn[iobs].residual/SECDAYl;
    }
    formBatsAll(psr,npsr);
    formResiduals(psr,npsr,0);

    for(int iobs = 0; iobs < psr->nobs; iobs++){
        EXPECT_LT(static_cast<double>(fabsl(psr[0].obsn[iobs].residual)),TEST_DELTA) << "Precision lost in formBats or formResiduals";
    }
}

TEST(testFormBats, offsetSATs){
    pulsar _psr;
    pulsar *psr = &_psr;
    MAX_PSR=1;
    char timFile[MAX_PSR][MAX_FILELEN],parFile[MAX_PSR][MAX_FILELEN];
    int npsr=1;
    initialise(psr,0); /* Initialise all */

    strcpy(parFile[0],DATDIR "/test1.par");
    strcpy(timFile[0],DATDIR "/test2.tim");

    readParfile(psr,parFile,timFile,npsr);
    readTimfile(psr,timFile,npsr);
    psr->noWarnings=2;
    formBatsAll(psr,npsr);

    for(int iobs = 0; iobs < psr->nobs; iobs++){
        psr->obsn[iobs].prefitResidual = psr->obsn[iobs].bat;
        psr->obsn[iobs].sat += longdouble(4e-9)/SECDAYl;
    }
    formBatsAll(psr,npsr);

    for(int iobs = 0; iobs < psr->nobs; iobs++){
        EXPECT_LT(static_cast<double>(fabsl(psr->obsn[iobs].bat - psr->obsn[iobs].prefitResidual - longdouble(4e-9)/SECDAYl)),TEST_DELTA) << "Precision lost in formBats";
    }

}
