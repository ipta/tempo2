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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tempo2.h"
#include "t2fit.h"

typedef struct FlagOutlierConfig {
	double sigmaThreshold;
	char outputTim[MAX_FILELEN + 1];
	char flagID[MAX_FLAG_LEN];
	char flagValue[MAX_FLAG_LEN];
	int subtractRed;
	int subtractDM;
	int subtractChrom;
	int subtractGroup;
	int medianCenter;
	int tempDeleteRefit;
	double prefitMaxDeleteFrac;
	int verbose;
} FlagOutlierConfig;

typedef struct ScoredObs {
	size_t idx;
	double score;
} ScoredObs;

typedef struct DetectionStats {
	long totalObs;
	long activeObs;
	long skippedSigma;
	long validForCenter;
	long outliers;
	double zMedian;
} DetectionStats;

static void printUsage(void)
{
	printf("flagOutliers output plugin options:\n");
	printf("  -sigma <value>       Outlier threshold in sigma units (default 10)\n");
	printf("  -out <filename>      Output tim filename (default flagOutliers.tim)\n");
	printf("  -flagID <name>       Flag key to set (default -distrust)\n");
	printf("  -flagVal <value>     Flag value to set (default outlier)\n");
	printf("  -noTNRed             Do not subtract TNRedSignal\n");
	printf("  -noTNDM              Do not subtract TNDMSignal\n");
	printf("  -noTNChrom           Do not subtract TNChromSignal\n");
	printf("  -noTNGroup           Do not subtract TNGroupSignal\n");
	printf("  -noMedianCenter      Do not subtract global median of normalized residuals\n");
	printf("  -tempDeleteRefit     Enable temporary delete+refit pre-stage (default on)\n");
	printf("  -noTempDeleteRefit   Disable temporary delete+refit pre-stage\n");
	printf("  -prefitMaxDeleteFrac Max fraction deleted in pre-stage (default 0.10)\n");
	printf("  -verbose             Print per-outlier details\n");
	printf("  -h / -help           Show this help\n");
}

static int cmpDouble(const void *a, const void *b)
{
	const double da = *(const double *)a;
	const double db = *(const double *)b;
	if (da < db)
		return -1;
	if (da > db)
		return 1;
	return 0;
}

static int cmpScoredObsDesc(const void *a, const void *b)
{
	const ScoredObs *sa = (const ScoredObs *)a;
	const ScoredObs *sb = (const ScoredObs *)b;
	if (sa->score > sb->score)
		return -1;
	if (sa->score < sb->score)
		return 1;
	if (sa->idx < sb->idx)
		return -1;
	if (sa->idx > sb->idx)
		return 1;
	return 0;
}

static double medianOfArray(const double *arr, size_t n)
{
	double *tmp;
	double med;

	if (n == 0)
		return 0.0;

	tmp = (double *)malloc(n * sizeof(double));
	if (tmp == NULL)
		return 0.0;

	memcpy(tmp, arr, n * sizeof(double));
	qsort(tmp, n, sizeof(double), cmpDouble);

	if ((n % 2) == 1)
		med = tmp[n / 2];
	else
		med = 0.5 * (tmp[(n / 2) - 1] + tmp[n / 2]);

	free(tmp);
	return med;
}

static double getCorrectedResidual(const observation *obs, const FlagOutlierConfig *cfg)
{
	double correctedResidual = (double)obs->residual;
	if (cfg->subtractRed)
		correctedResidual -= obs->TNRedSignal;
	if (cfg->subtractDM)
		correctedResidual -= obs->TNDMSignal;
	if (cfg->subtractChrom)
		correctedResidual -= obs->TNChromSignal;
	if (cfg->subtractGroup)
		correctedResidual -= obs->TNGroupSignal;
	return correctedResidual;
}

static double getWhiteSigma(const observation *obs)
{
	/* toaErr is already updated by read/preProcess EFAC/EQUAD logic. */
	return (double)obs->toaErr * 1.0e-6;
}

static void detectOutliers(pulsar *psr, int npsr, const FlagOutlierConfig *cfg,
					   unsigned char *outlierMask, DetectionStats *stats,
					   double *absCenteredZ)
{
	int p;
	size_t maxObs = 0;
	double *rawZ = NULL;
	long validForCenter = 0;
	double zMedian = 0.0;
	size_t idx = 0;

	stats->totalObs = 0;
	stats->activeObs = 0;
	stats->skippedSigma = 0;
	stats->validForCenter = 0;
	stats->outliers = 0;
	stats->zMedian = 0.0;

	for (p = 0; p < npsr; p++)
		maxObs += (size_t)psr[p].nobs;

	if (maxObs > 0) {
		rawZ = (double *)malloc(maxObs * sizeof(double));
		if (rawZ == NULL) {
			printf("ERROR: unable to allocate memory for outlier detection\n");
			return;
		}
	}

	idx = 0;
	for (p = 0; p < npsr; p++) {
		int i;
		for (i = 0; i < psr[p].nobs; i++) {
			observation *obs = &psr[p].obsn[i];
			double correctedResidual;
			double sigma;
			double z;

			outlierMask[idx] = 0;
			if (absCenteredZ != NULL)
				absCenteredZ[idx] = -1.0;
			stats->totalObs++;

			if (obs->deleted != 0) {
				idx++;
				continue;
			}

			correctedResidual = getCorrectedResidual(obs, cfg);
			sigma = getWhiteSigma(obs);
			if (!isfinite(sigma) || sigma <= 0.0) {
				idx++;
				continue;
			}

			z = correctedResidual / sigma;
			if (!isfinite(z)) {
				idx++;
				continue;
			}

			rawZ[validForCenter++] = z;
			idx++;
		}
	}

	if (cfg->medianCenter && validForCenter > 0)
		zMedian = medianOfArray(rawZ, (size_t)validForCenter);

	idx = 0;
	for (p = 0; p < npsr; p++) {
		int i;
		for (i = 0; i < psr[p].nobs; i++) {
			observation *obs = &psr[p].obsn[i];
			double correctedResidual;
			double sigma;
			double z;
			double zUsed;

			if (obs->deleted != 0) {
				idx++;
				continue;
			}

			stats->activeObs++;
			correctedResidual = getCorrectedResidual(obs, cfg);
			sigma = getWhiteSigma(obs);

			if (!isfinite(sigma) || sigma <= 0.0) {
				stats->skippedSigma++;
				idx++;
				continue;
			}

			z = correctedResidual / sigma;
			if (!isfinite(z)) {
				stats->skippedSigma++;
				idx++;
				continue;
			}

			zUsed = z - zMedian;
			if (absCenteredZ != NULL)
				absCenteredZ[idx] = fabs(zUsed);
			if (fabs(zUsed) > cfg->sigmaThreshold) {
				outlierMask[idx] = 1;
				stats->outliers++;
			}

			idx++;
		}
	}

	stats->validForCenter = validForCenter;
	stats->zMedian = zMedian;
	free(rawZ);
}

static int findFlagIndex(const observation *obs, const char *flagID)
{
	int k;
	for (k = 0; k < obs->nFlags; k++) {
		if (strcmp(obs->flagID[k], flagID) == 0)
			return k;
	}
	return -1;
}

static int parseArgs(int argc, char *argv[], FlagOutlierConfig *cfg)
{
	int i;

	cfg->sigmaThreshold = 10.0;
	strcpy(cfg->outputTim, "flagOutliers.tim");
	strcpy(cfg->flagID, "-distrust");
	strcpy(cfg->flagValue, "outlier");
	cfg->subtractRed = 1;
	cfg->subtractDM = 1;
	cfg->subtractChrom = 1;
	cfg->subtractGroup = 1;
	cfg->medianCenter = 1;
	cfg->tempDeleteRefit = 1;
	cfg->prefitMaxDeleteFrac = 0.10;
	cfg->verbose = 0;

	for (i = 0; i < argc; i++) {
		if (strcasecmp(argv[i], "-sigma") == 0) {
			if (i + 1 >= argc) {
				printf("ERROR: -sigma requires a value\n");
				return -1;
			}
			cfg->sigmaThreshold = atof(argv[++i]);
		} else if (strcasecmp(argv[i], "-out") == 0) {
			if (i + 1 >= argc) {
				printf("ERROR: -out requires a filename\n");
				return -1;
			}
			strncpy(cfg->outputTim, argv[++i], MAX_FILELEN);
			cfg->outputTim[MAX_FILELEN] = '\0';
		} else if (strcasecmp(argv[i], "-flagID") == 0) {
			if (i + 1 >= argc) {
				printf("ERROR: -flagID requires a value\n");
				return -1;
			}
			strncpy(cfg->flagID, argv[++i], MAX_FLAG_LEN - 1);
			cfg->flagID[MAX_FLAG_LEN - 1] = '\0';
		} else if (strcasecmp(argv[i], "-flagVal") == 0) {
			if (i + 1 >= argc) {
				printf("ERROR: -flagVal requires a value\n");
				return -1;
			}
			strncpy(cfg->flagValue, argv[++i], MAX_FLAG_LEN - 1);
			cfg->flagValue[MAX_FLAG_LEN - 1] = '\0';
		} else if (strcasecmp(argv[i], "-noTNRed") == 0) {
			cfg->subtractRed = 0;
		} else if (strcasecmp(argv[i], "-noTNDM") == 0) {
			cfg->subtractDM = 0;
		} else if (strcasecmp(argv[i], "-noTNChrom") == 0) {
			cfg->subtractChrom = 0;
		} else if (strcasecmp(argv[i], "-noTNGroup") == 0) {
			cfg->subtractGroup = 0;
		} else if (strcasecmp(argv[i], "-noMedianCenter") == 0) {
			cfg->medianCenter = 0;
		} else if (strcasecmp(argv[i], "-tempDeleteRefit") == 0) {
			cfg->tempDeleteRefit = 1;
		} else if (strcasecmp(argv[i], "-noTempDeleteRefit") == 0) {
			cfg->tempDeleteRefit = 0;
		} else if (strcasecmp(argv[i], "-prefitMaxDeleteFrac") == 0) {
			if (i + 1 >= argc) {
				printf("ERROR: -prefitMaxDeleteFrac requires a value\n");
				return -1;
			}
			cfg->prefitMaxDeleteFrac = atof(argv[++i]);
		} else if (strcasecmp(argv[i], "-verbose") == 0) {
			cfg->verbose = 1;
		} else if (strcasecmp(argv[i], "-h") == 0 || strcasecmp(argv[i], "-help") == 0) {
			printUsage();
			return 1;
		}
	}

	if (!isfinite(cfg->sigmaThreshold) || cfg->sigmaThreshold <= 0.0) {
		printf("ERROR: invalid sigma threshold: %g\n", cfg->sigmaThreshold);
		return -1;
	}

	if (!isfinite(cfg->prefitMaxDeleteFrac) || cfg->prefitMaxDeleteFrac < 0.0 || cfg->prefitMaxDeleteFrac > 1.0) {
		printf("ERROR: invalid prefit max delete fraction: %g\n", cfg->prefitMaxDeleteFrac);
		return -1;
	}

	if (cfg->flagID[0] == '\0' || cfg->flagValue[0] == '\0') {
		printf("ERROR: flag ID and value must not be empty\n");
		return -1;
	}

	return 0;
}

#ifdef __cplusplus
extern "C" {
#endif

int tempoOutput(int argc, char *argv[], pulsar *psr, int npsr)
{
	FlagOutlierConfig cfg;
	DetectionStats finalStats;
	DetectionStats prelimStats;
	long outliers = 0;
	long addedFlags = 0;
	long existingFlags = 0;
	long tempDeleted = 0;
	long refitPerformed = 0;
	long prefitDeleteCap = 0;
	size_t maxObs = 0;
	unsigned char *outlierMask = NULL;
	unsigned char *origDeleted = NULL;
	double *absCenteredZ = NULL;
	double effectivePrefitThreshold = 0.0;
	int p;

	int parseStatus = parseArgs(argc, argv, &cfg);
	if (parseStatus != 0)
		return (parseStatus > 0) ? 0 : -1;

	for (p = 0; p < npsr; p++)
		maxObs += (size_t)psr[p].nobs;

	if (maxObs > 0) {
		outlierMask = (unsigned char *)calloc(maxObs, sizeof(unsigned char));
		origDeleted = (unsigned char *)calloc(maxObs, sizeof(unsigned char));
		absCenteredZ = (double *)malloc(maxObs * sizeof(double));
		if (outlierMask == NULL || origDeleted == NULL || absCenteredZ == NULL) {
			printf("ERROR: unable to allocate memory for outlier workflow\n");
			free(outlierMask);
			free(origDeleted);
			free(absCenteredZ);
			return -1;
		}
	}

	detectOutliers(psr, npsr, &cfg, outlierMask, &prelimStats, absCenteredZ);
	effectivePrefitThreshold = cfg.sigmaThreshold;

	if (cfg.tempDeleteRefit) {
		if (prelimStats.activeObs > 0) {
			prefitDeleteCap = (long)floor(cfg.prefitMaxDeleteFrac * (double)prelimStats.activeObs);
			if (prefitDeleteCap < 0)
				prefitDeleteCap = 0;
		}

		if (prelimStats.outliers > prefitDeleteCap) {
			ScoredObs *scored = NULL;
			long nScored = 0;
			size_t idx = 0;

			if (prelimStats.outliers > 0)
				scored = (ScoredObs *)malloc((size_t)prelimStats.outliers * sizeof(ScoredObs));

			if (scored != NULL || prelimStats.outliers == 0) {
				for (idx = 0; idx < maxObs; idx++) {
					if (!outlierMask[idx])
						continue;
					scored[nScored].idx = idx;
					scored[nScored].score = absCenteredZ[idx];
					nScored++;
				}

				qsort(scored, (size_t)nScored, sizeof(ScoredObs), cmpScoredObsDesc);

				for (idx = 0; idx < maxObs; idx++)
					outlierMask[idx] = 0;

				if (prefitDeleteCap > nScored)
					prefitDeleteCap = nScored;

				for (idx = 0; idx < (size_t)prefitDeleteCap; idx++)
					outlierMask[scored[idx].idx] = 1;

				if (prefitDeleteCap > 0)
					effectivePrefitThreshold = scored[prefitDeleteCap - 1].score;
				else
					effectivePrefitThreshold = HUGE_VAL;

				free(scored);
			} else {
				printf("WARNING: unable to allocate memory for prefit cap selection; using uncapped preliminary mask\n");
			}
		}

		size_t idx = 0;
		for (p = 0; p < npsr; p++) {
			int i;
			for (i = 0; i < psr[p].nobs; i++) {
				origDeleted[idx] = (unsigned char)(psr[p].obsn[i].deleted != 0);
				if (!origDeleted[idx] && outlierMask[idx]) {
					psr[p].obsn[i].deleted = 1;
					tempDeleted++;
				}
				idx++;
			}
		}

		if (tempDeleted > 0) {
			formBatsAll(psr, npsr);
			formResiduals(psr, npsr, 1);
			t2Fit(psr, (unsigned int)npsr, covarFuncFile);
			formBatsAll(psr, npsr);
			formResiduals(psr, npsr, 1);
			refitPerformed = 1;
		}

		/* Restore original deleted state; final flagging should be on full input set. */
		{
			size_t ridx = 0;
			for (p = 0; p < npsr; p++) {
				int i;
				for (i = 0; i < psr[p].nobs; i++) {
					psr[p].obsn[i].deleted = origDeleted[ridx] ? 1 : 0;
					ridx++;
				}
			}
			formBatsAll(psr, npsr);
			formResiduals(psr, npsr, 1);
		}
	}

	/* Final outlier pass after optional refit. */
	detectOutliers(psr, npsr, &cfg, outlierMask, &finalStats, NULL);
	outliers = finalStats.outliers;

	{
		size_t idx = 0;
		for (p = 0; p < npsr; p++) {
			int i;
			for (i = 0; i < psr[p].nobs; i++) {
				observation *obs = &psr[p].obsn[i];

				if (obs->deleted != 0 || !outlierMask[idx]) {
					idx++;
					continue;
				}

				{
					int fidx = findFlagIndex(obs, cfg.flagID);
					if (fidx < 0) {
						if (obs->nFlags < MAX_FLAGS) {
							fidx = obs->nFlags;
							strcpy(obs->flagID[fidx], cfg.flagID);
							strcpy(obs->flagVal[fidx], cfg.flagValue);
							obs->nFlags++;
							addedFlags++;
						} else {
							printf("WARNING: cannot add flag to observation %d for pulsar %s (MAX_FLAGS reached)\n", i, psr[p].name);
						}
					} else {
						existingFlags++;
					}
				}

				if (cfg.verbose) {
					double correctedResidual = getCorrectedResidual(obs, &cfg);
					double sigma = getWhiteSigma(obs);
					double zRaw = correctedResidual / sigma;
					double zUsed = zRaw - finalStats.zMedian;
					printf("Outlier: psr=%s iobs=%d sat=%.15Lf z_raw=%g z_centered=%g z_med=%g resid=%g sigma=%g\n",
						   psr[p].name, i, obs->sat, zRaw, zUsed, finalStats.zMedian, correctedResidual, sigma);
				}

				idx++;
			}
		}
	}

	free(outlierMask);
	free(origDeleted);
	free(absCenteredZ);

	writeTim(cfg.outputTim, psr, "tempo2");

	printf("flagOutliers summary:\n");
	printf("  total observations: %ld\n", finalStats.totalObs);
	printf("  active observations: %ld\n", finalStats.activeObs);
	printf("  skipped (invalid sigma/z): %ld\n", finalStats.skippedSigma);
	printf("  valid for centering: %ld\n", finalStats.validForCenter);
	printf("  z-median subtracted: %.6g\n", finalStats.zMedian);
	printf("  preliminary outliers: %ld\n", prelimStats.outliers);
	printf("  prefit max delete fraction: %.3f\n", cfg.prefitMaxDeleteFrac);
	printf("  prefit delete cap: %ld\n", prefitDeleteCap);
	printf("  effective prefit threshold: %s", isfinite(effectivePrefitThreshold) ? "" : "inf");
	if (isfinite(effectivePrefitThreshold))
		printf("%.6g\n", effectivePrefitThreshold);
	else
		printf("\n");
	printf("  temp deleted for refit: %ld\n", tempDeleted);
	printf("  refit performed: %s\n", refitPerformed ? "yes" : "no");
	printf("  outliers found: %ld\n", outliers);
	printf("  new flags added: %ld\n", addedFlags);
	printf("  already flagged: %ld\n", existingFlags);
	printf("  output tim file: %s\n", cfg.outputTim);

	return 0;
}

#ifdef __cplusplus
}
#endif
