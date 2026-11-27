//  Copyright (C) 2006,2007,2008,2009, George Hobbs, Russel Edwards

/*
 *    This file is part of TEMPO2. 
 * 
 *    TEMPO2 is free software: you can redistribute it and/or modify 
 *    it under the terms of the GNU General Public License as published by 
 *    the Free Software Foundation, either version 3 of the License, or 
 *    (at your option) any later version. 
 *    TEMPO2 is distributed in the hope that it will be useful, 
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of 
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 *    GNU General Public License for more details. 
 *    You should have received a copy of the GNU General Public License 
 *    along with TEMPO2.  If not, see <http://www.gnu.org/licenses/>. 
 */

/*
 *    If you use TEMPO2 then please acknowledge it by citing 
 *    Hobbs, Edwards & Manchester (2006) MNRAS, Vol 369, Issue 2, 
 *    pp. 655-672 (bibtex: 2006MNRAS.369..655H)
 *    or Edwards, Hobbs & Manchester (2006) MNRAS, VOl 372, Issue 4,
 *    pp. 1549-1574 (bibtex: 2006MNRAS.372.1549E) when discussing the
 *    timing model.
 */

#include "tempo2pred_int.h"
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/* Tempo1 polyco functions destined for tempo2pred library (C linkage) */

longdouble T1Polyco_GetPhase(const T1Polyco *t1p, longdouble mjd, longdouble freq)
{
    longdouble dt, arg, phase, spin_freq;
    int ic;

    /* Evaluate polynomial in phase and frequency */
    dt = (mjd-t1p->mjd_mid)*1440;
    printf("WARNIGN WARNGING: Silly update in polyco read\n");
    printf("Have %.20lg %.20lg\n",(double)mjd,(double)t1p->mjd_mid);
    printf("Diff = %.15lg\n",(double)(mjd-t1p->mjd_mid));

    phase = t1p->reference_phase + dt*60*t1p->frequency_psr_0 + t1p->coeff[0];
    spin_freq = t1p->frequency_psr_0 + t1p->coeff[1]/60;
    arg = dt;
    for (ic=1; ic<t1p->ncoeff; ic++)
    {
        phase += t1p->coeff[ic]*arg;
        if (ic < t1p->ncoeff-1)
            spin_freq += (ic+1)*t1p->coeff[ic+1]*arg/60;
        arg *= dt;
    }
    /* compute dispersion phase delay in observatory frame */
    phase -= spin_freq * t1p->dm/longdouble(2.41e-4) /
        (1.0L+t1p->doppler*longdouble(1.0e-4)) *
        ( 1.0L/(freq*freq) - 1.0L/(t1p->frequency_obs*t1p->frequency_obs));
    return phase;
}

longdouble T1Polyco_GetFrequency(const T1Polyco *t1p, longdouble mjd,
        longdouble freq)
{
    longdouble dt, arg, spin_freq;
    int ic;

    dt = (mjd-t1p->mjd_mid)*1440;
    spin_freq = t1p->frequency_psr_0 + t1p->coeff[1]/60;
    arg = dt;

    for (ic=2; ic<t1p->ncoeff; ic++)
    {
        spin_freq += ic*t1p->coeff[ic]*arg/60;
        arg *= dt;
    }

    return spin_freq;
}

void T1Polyco_Write(const T1Polyco *t1p, FILE *fout)
{
    int i;

    if (!t1p->extended_format) {
        if (t1p->psrname[0]=='J')
            fprintf(fout,"%-10.10s ",t1p->psrname+1);
        else
            fprintf(fout,"%-10.10s ",t1p->psrname);
        fprintf(fout,"%9.9s%11.11s",t1p->date_string,t1p->utc_string);
        ld_fprintf(fout,"%20.11Lf",t1p->mjd_mid);
        fprintf(fout,"%21.6lf %6.3lf%7.3lf\n",t1p->dm,t1p->doppler,t1p->log10rms);
        ld_fprintf(fout,"%20.6Lf%18.12Lf",t1p->reference_phase,t1p->frequency_psr_0);
        fprintf(fout,"%5s%5d%5d%10.3lf%16.14lf\n",t1p->sitename,t1p->span,
                t1p->ncoeff,t1p->frequency_obs,t1p->binary_phase);
        for (i=0;i<t1p->ncoeff;i++) {
            ld_fprintf(fout,"%25.17Le",t1p->coeff[i]);
            if ((i+1)%3==0)
                fprintf(fout,"\n");
        }
        if (t1p->ncoeff%3)
            fprintf(fout,"\n");
        return;
    }

    fprintf(fout, "TEMPO2: POLYCO TEMPO1 emulation\n%s\n%s\n%s\n",
            t1p->psrname, t1p->date_string, t1p->utc_string);
    ld_fprintf(fout, "%.36Lg\n", t1p->mjd_mid);
    fprintf(fout, "%.17g\n%.17g\n%.17g\n", t1p->dm, t1p->doppler, t1p->log10rms);
    ld_fprintf(fout, "%.36Lg\n%.36Lg\n", t1p->reference_phase, t1p->frequency_psr_0);
    fprintf(fout, "%s\n%d\n%d\n%.17g\n%.17g\n%.17g\n", t1p->sitename,
            t1p->span, t1p->ncoeff, t1p->frequency_obs,
            t1p->binary_phase, t1p->binary_frequency);
    for (i=0;i<t1p->ncoeff;i++)
        ld_fprintf(fout, "%.36Lg\n", t1p->coeff[i]);
}

/* helper functions for T1Polyco_Read */
void T1P_grabString(char *str, int istart, int nchar, char *out)
{
    char grabbed[64];

    strncpy(grabbed, str+istart, nchar);
    grabbed[nchar] = '\0';
    sscanf(grabbed, "%s", out);// remove leading/trailing whitespace
}

longdouble T1P_grabLongDouble(char *str, int istart, int nchar)
{
    char grabbed[64];
    char *c;

    T1P_grabString(str, istart, nchar, grabbed);

    /* change D to e for exponential notation*/
    for (c=grabbed; *c; c++)
        if (*c == 'D')
            *c = 'e';

    return parse_longdouble(grabbed);
}

static int parse_polyco_number(const char *text, longdouble *value)
{
    char *end;

#ifdef LONGDOUBLE_IS_FLOAT128
    *value = strtoflt128(text, &end);
#elif defined(LONGDOUBLE_IS_IEEE754)
    *value = strtold(text, &end);
#else
    strtold(text, &end);
    *value = parse_longdouble(text);
#endif
    return end != text && *end == '\0';
}

longdouble T1P_grabInt(char *str, int istart, int nchar)
{
    char grabbed[64];
    int res;

    T1P_grabString(str, istart, nchar, grabbed);
    sscanf(grabbed, "%d", &res);

    return res;
}

int T1Polyco_Read_NewFormat(T1Polyco *t1p, FILE *f)
{
    int ic;
    char junk[1024], mjd_mid[128], reference_phase[128], frequency_psr_0[128], coeff[128];

    if (fscanf(f, "%63s %9s %12s %127s %lf %lf %lf %127s %127s %4s %d %d %lf %lf %lf",
                t1p->psrname, t1p->date_string, t1p->utc_string,
                mjd_mid, &t1p->dm, &t1p->doppler, &t1p->log10rms,
                reference_phase, frequency_psr_0,
                t1p->sitename, &t1p->span, &t1p->ncoeff,
                &t1p->frequency_obs, &t1p->binary_phase, &t1p->binary_frequency)
            != 15)
        return -1;
    if (t1p->ncoeff < 2 || t1p->ncoeff > 32)
        return -1;
    if (!parse_polyco_number(mjd_mid, &t1p->mjd_mid) ||
            !parse_polyco_number(reference_phase, &t1p->reference_phase) ||
            !parse_polyco_number(frequency_psr_0, &t1p->frequency_psr_0))
        return -1;

    for (ic=0; ic < t1p->ncoeff; ic++){
        if (fscanf(f, "%127s", coeff)!=1)
            return -1;
        if (!parse_polyco_number(coeff, &t1p->coeff[ic]))
            return -1;
    }
    fgets(junk, 1024, f); /* eat final newline */
    return 0;
}

int T1Polyco_Read(T1Polyco *t1p, FILE *f)
{
    char line[128];
    int ic;

    if (fgets(line, 128, f)!=line)
        return -1;

    if (!strncasecmp(line, "TEMPO2:", 6))
    {
        t1p->extended_format = 1;
        return T1Polyco_Read_NewFormat(t1p, f);
    }

    t1p->extended_format = 0;

    T1P_grabString(line, 0, 11, t1p->psrname);
    T1P_grabString(line, 11, 9, t1p->date_string);
    T1P_grabString(line, 20, 11, t1p->utc_string);
    t1p->mjd_mid = T1P_grabLongDouble(line, 31, 20);
    t1p->dm = T1P_grabLongDouble(line, 51, 21);
    t1p->doppler = T1P_grabLongDouble(line, 73, 6);
    t1p->log10rms = T1P_grabLongDouble(line, 79, 7);

    if (fgets(line, 128, f)!=line)
        return -1;
    t1p->reference_phase = T1P_grabLongDouble(line, 0, 20);
    t1p->frequency_psr_0 = T1P_grabLongDouble(line, 21, 18);
    T1P_grabString(line, 38, 5, t1p->sitename);
    t1p->span = T1P_grabInt(line, 43, 5);
    t1p->ncoeff = T1P_grabInt(line, 48, 5);
    if (t1p->ncoeff < 2 || t1p->ncoeff > 32)
        return -1;
    t1p->frequency_obs = T1P_grabLongDouble(line, 53, 10);
    t1p->binary_phase = T1P_grabLongDouble(line, 63, 7);
    t1p->binary_frequency = T1P_grabLongDouble(line, 70, 9);

    for (ic=0; ic < t1p->ncoeff; ic+=3)
    {
        if (fgets(line, 128, f)!=line)
            return -1;
        t1p->coeff[ic] = T1P_grabLongDouble(line, 0, 25);
        if (ic+1 < t1p->ncoeff)
            t1p->coeff[ic+1] = T1P_grabLongDouble(line, 25, 25);
        if (ic+2 < t1p->ncoeff)
            t1p->coeff[ic+2] = T1P_grabLongDouble(line, 50, 25);
    }

    return 0;
}

T1Polyco *T1PolycoSet_GetNearest(const T1PolycoSet *t1ps, longdouble mjd)
{
    int inearest=0;
    longdouble best_offset=1e6, offset;
    int iseg;
    for (iseg=0; iseg < t1ps->nsegments ; iseg++)
    {
        offset = fabsl(t1ps->segments[iseg].mjd_mid - mjd);
        if (offset < best_offset)
        {
            inearest = iseg;
            best_offset = offset;
        }
    }

    return &t1ps->segments[inearest];
}

longdouble T1PolycoSet_GetPhase(const T1PolycoSet *t1ps, longdouble mjd, longdouble freq)
{
    return T1Polyco_GetPhase(T1PolycoSet_GetNearest(t1ps, mjd), mjd, freq);
}

longdouble T1PolycoSet_GetFrequency(const T1PolycoSet *t1ps, longdouble mjd, longdouble freq)
{
    return T1Polyco_GetFrequency(T1PolycoSet_GetNearest(t1ps, mjd), mjd, freq);
}

void T1PolycoSet_Write(const T1PolycoSet *t1ps, FILE *f)
{
    int iseg;
    for (iseg=0; iseg < t1ps->nsegments ; iseg++)
        T1Polyco_Write(&t1ps->segments[iseg], f);
}

int T1PolycoSet_Read(T1PolycoSet *t1ps, FILE *f)
{
    int nallocated = 5;
    int iseg=0;

    t1ps->segments = (T1Polyco *)malloc(nallocated*sizeof(T1Polyco));
    t1ps->nsegments = 0;
    if (!t1ps->segments)
        return -1;

    while (T1Polyco_Read(&t1ps->segments[iseg], f)==0)
    {
        iseg++;
        if (iseg >= nallocated)
        {
            T1Polyco *segments;
            nallocated += 5;
            segments = (T1Polyco *)realloc(t1ps->segments,
                    nallocated*sizeof(T1Polyco));
            if (!segments) {
                free(t1ps->segments);
                t1ps->segments = NULL;
                return -1;
            }
            t1ps->segments = segments;
        }
    }
    t1ps->nsegments = iseg;

    if (t1ps->nsegments > 0)
        return 0;
    else
    {
        free(t1ps->segments);
        t1ps->segments = NULL;
        return -1;
    }
}

void T1PolycoSet_Init(T1PolycoSet *t1ps)
{
  t1ps->segments = NULL;
  t1ps->nsegments = 0;
}

void T1PolycoSet_Destroy(T1PolycoSet *t1ps)
{
  if (t1ps->segments)
    free(t1ps->segments);
  T1PolycoSet_Init(t1ps);
}
