//  Copyright (C) 2006,2007,2008,2009, George Hobbs, Russell Edwards

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

//-*-C++-*-
#ifndef TEMPO2PRED_H
#define TEMPO2PRED_H

#include <stdio.h>
#include <stdint.h>

#define T2PREDICTOR_API_VERSION 2

#ifdef __cplusplus
#include <string>
extern "C" {
#endif

    /*******************************************************************/
    /* TEMPO2 Predictive mode library interface for C/C++              */
    /* Initial version: Russell Edwards 14 Feb 2006.                   */
    /*******************************************************************/

    /************************* C TYPES *********************************/
    /* Note: Most of these types are for internal use only. The        */
    /*       external caller need only know about T2Predictor and      */
    /*       its associated functions declared further down            */
    /*******************************************************************/

    /* Segment contents are private and depend on the configured precision. */
    typedef struct ChebyModel ChebyModel;
    typedef struct T1Polyco T1Polyco;

    /* Set of Chebyshev predictive phase models */
    typedef struct
    {
        ChebyModel *segments;
        int nsegments;
    } ChebyModelSet;

    /* Flag set by ChebyModelSet_GetNearest when mjd is out of range */
    extern int ChebyModelSet_OutOfRange;

    /* Set of TEMPO1-style polycos */
    typedef struct
    {
        T1Polyco *segments;
        int nsegments;
    } T1PolycoSet;

    /* Types of models available */
    typedef enum  {NonePredType, Cheby, T1} T2PredictorKind;

    /* Generic set-of-predictive models type */
    typedef struct
    {
        T2PredictorKind kind;
        union
        {
            ChebyModelSet cheby;
            T1PolycoSet t1;
        } modelset;
    } T2Predictor;

    /************************* C FUNCTIONS ******************************/

    /* Initialization, Copying, Building, etc */
    void T2Predictor_Init(T2Predictor *t2p);
    void T2Predictor_Copy(T2Predictor *into_t2p, const T2Predictor *from_t2p);
    int  T2Predictor_Insert(T2Predictor *into_t2p, const T2Predictor *from_t2p);
    void T2Predictor_Keep(T2Predictor*, unsigned nmjd, const long double* mjd);
    void T2Predictor_Destroy(T2Predictor *t2p);

    /* Initialize before reading. FRead replaces existing contents; a failed
       parse leaves an empty predictor. Read leaves it unchanged if open fails. */
    int T2Predictor_Read(T2Predictor *t2p, char *fname);
    int T2Predictor_FRead(T2Predictor *t2p, FILE *f);
    void T2Predictor_Write(const T2Predictor *t2p, char *fname);
    void T2Predictor_FWrite(const T2Predictor *t2p, FILE *f);

    /* Information */
    char * T2Predictor_GetPSRName(T2Predictor *t2p);
    char * T2Predictor_GetSiteName(T2Predictor *t2p);
    long double T2Predictor_GetStartMJD(T2Predictor *t2p);
    long double T2Predictor_GetEndMJD(T2Predictor *t2p);
    long double T2Predictor_GetStartFreq(T2Predictor *t2p); // MHz
    long double T2Predictor_GetEndFreq(T2Predictor *t2p);  // MHz
    T2PredictorKind T2Predictor_Kind(T2Predictor *t2p);

    typedef struct {
        int64_t day;
        long double fraction;
    } T2PredictorMJD;

    typedef struct {
        int64_t turns;
        long double fraction;
    } T2PredictorPhase;

    typedef enum {
        T2PREDICTOR_OK = 0,
        T2PREDICTOR_INVALID = -1,
        T2PREDICTOR_WRONG_KIND = -2,
        T2PREDICTOR_OUT_OF_RANGE = -3,
        T2PREDICTOR_OVERFLOW = -4,
        T2PREDICTOR_NOMEM = -5,
        T2PREDICTOR_IO_ERROR = -6
    } T2PredictorStatus;

    typedef struct {
        T2PredictorMJD start, end, midpoint;
        int time_coefficients, frequency_coefficients;
        long double length_days, frequency_start, frequency_end;
        long double dispersion_constant;
    } T2PredictorSegmentInfo;

    /*
     * Split values represent day + fraction and turns + fraction. Outputs
     * are normalized to 0 <= fraction < 1, including for negative values.
     * Input fractions may be any finite value. Outputs are unchanged on error.
     * Precision is limited by the configured internal longdouble type.
     */
    T2PredictorStatus T2Predictor_GetPhaseSplit(const T2Predictor *,
            T2PredictorMJD, long double freq, T2PredictorPhase *);
    T2PredictorStatus T2Predictor_GetFrequencySplit(const T2Predictor *,
            T2PredictorMJD, long double freq, long double *);
    T2PredictorStatus T2Predictor_GetStartMJDSplit(const T2Predictor *, T2PredictorMJD *);
    T2PredictorStatus T2Predictor_GetEndMJDSplit(const T2Predictor *, T2PredictorMJD *);
    T2PredictorStatus T2Predictor_KeepSplit(T2Predictor *, unsigned,
            const T2PredictorMJD *);

    /* Segment indexes become invalid after read, insert, keep or destroy.
     * Indexed evaluations permit extrapolation, as the original segment
     * functions do. Whole-predictor Chebyshev phase/frequency calls update
     * ChebyModelSet_OutOfRange; indexed evaluations do not change it.
     * Split T1 evaluations require the epoch to be within the nearest
     * polyco's nominal span; the legacy calls allow extrapolation.
     */
    T2PredictorStatus T2Predictor_GetSegmentCount(const T2Predictor *, unsigned *);
    T2PredictorStatus T2Predictor_GetSegmentInfo(const T2Predictor *, unsigned,
            T2PredictorSegmentInfo *);
    T2PredictorStatus T2Predictor_ValidateSegments(const T2Predictor *);
    T2PredictorStatus T2Predictor_GetSegmentPhaseSplit(const T2Predictor *,
            unsigned, T2PredictorMJD, long double freq, T2PredictorPhase *);
    T2PredictorStatus T2Predictor_GetSegmentFrequencySplit(const T2Predictor *,
            unsigned, T2PredictorMJD, long double freq, long double *);

    /* The output frequency array has capacity elements. */
    T2PredictorStatus T2Predictor_GetPlanSplit(const char *, T2PredictorMJD,
            T2PredictorMJD, long double step, long double freq,
            T2PredictorPhase *, int *, long double *, unsigned capacity);
    T2PredictorStatus T2Predictor_GetPlanSplit_Ext(const char *, T2PredictorMJD,
            T2PredictorMJD, long double step, long double freq,
            char *psrname, size_t psrname_size, char *sitename, size_t sitename_size,
            T2PredictorPhase *, int *,
            long double *, unsigned capacity);

    /* Prediction */
    long double T2Predictor_GetPhase(const T2Predictor *t2p, long double mjd,
            long double freq); // freq in MHz
    long double T2Predictor_GetFrequency(const T2Predictor *t2p, long double mjd,
            long double freq); // freq in MHz
    /* This function makes a piecewise approximation. At the moment it just
       interpolates between phase evaluations. In future it will minimise
       the mean offset. Returns 0 on success */
    int T2Predictor_GetPlan(char *filename,
            long double mjd_start,
            long double mjd_end,
            long double step, // seconds
            long double freq, // MHz
            // remaining arguments are returned
            long double *phase0,
            int *nsegments,
            long double *pulse_frequencies);
    int T2Predictor_GetPlan_Ext(char *filename,
            long double mjd_start,
            long double mjd_end,
            long double step, // seconds
            long double freq, // MHz
            // remaining arguments are returned
            char *psrname, char *sitename,
            long double *phase0,
            int *nsegments,
            long double *pulse_frequencies);


#ifdef __cplusplus
}
#endif

#endif
