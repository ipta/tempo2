#include <tempo2.h>
#include <stdint.h>

static constexpr uint64_t T2_TNSUBPOLY_FLAG_POSTFIT  = 1ULL;
static constexpr uint64_t T2_TNSUBPOLY_FLAG_NESTLIKE = 2ULL;

static inline int t2ClampSubPolyOrder(int order)
{
	if (order < 0) return 0;
	if (order > 3) return 3;
	return order;
}

static inline int t2DecodeTNSubPolyOrdStored(uint64_t packed, int shift)
{
	return t2ClampSubPolyOrder((int)((packed >> shift) & 0xFFULL));
}


static inline uint64_t t2EncodeTNsubtractPoly(uint64_t flags, int redOrder, int dmOrder, int chromOrder)
{
	uint64_t v = (flags & 0xFFULL);
	v |= ((uint64_t)(t2ClampSubPolyOrder(redOrder)   & 0xFF) << 8);
	v |= ((uint64_t)(t2ClampSubPolyOrder(dmOrder)    & 0xFF) << 16);
	v |= ((uint64_t)(t2ClampSubPolyOrder(chromOrder) & 0xFF) << 24);
	return v;
}

static inline uint64_t t2DecodeTNsubtractPolyFlags(uint64_t packed)
{
	return (packed & 0xFFULL);
}


static inline int t2DecodeTNRedSubPolyOrd(uint64_t packed)
{
	return t2DecodeTNSubPolyOrdStored(packed, 8);
}

static inline int t2DecodeTNDMSubPolyOrd(uint64_t packed)
{
	return t2DecodeTNSubPolyOrdStored(packed, 16);
}

static inline int t2DecodeTNChromSubPolyOrd(uint64_t packed)
{
	return t2DecodeTNSubPolyOrdStored(packed, 24);
}

static inline uint64_t t2SetTNRedSubPolyOrd(uint64_t packed, int order)
{
	packed &= ~(0xFFULL << 8);
	packed |= ((uint64_t)(t2ClampSubPolyOrder(order) & 0xFF) << 8);
	return packed;
}

static inline uint64_t t2SetTNDMSubPolyOrd(uint64_t packed, int order)
{
	packed &= ~(0xFFULL << 16);
	packed |= ((uint64_t)(t2ClampSubPolyOrder(order) & 0xFF) << 16);
	return packed;
}

static inline uint64_t t2SetTNChromSubPolyOrd(uint64_t packed, int order)
{
	packed &= ~(0xFFULL << 24);
	packed |= ((uint64_t)(t2ClampSubPolyOrder(order) & 0xFF) << 24);
	return packed;
}

double t2FitFunc_nestlike_red(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);
void t2UpdateFunc_nestlike_red(pulsar *psr, int ipsr ,param_label label,int k, double val, double err);

double t2FitFunc_nestlike_red_dm(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);
void t2UpdateFunc_nestlike_red_dm(pulsar *psr, int ipsr ,param_label label,int k, double val, double err);


double t2FitFunc_nestlike_red_chrom(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);
void t2UpdateFunc_nestlike_red_chrom(pulsar *psr, int ipsr ,param_label label,int k, double val, double err);


double t2FitFunc_nestlike_jitter(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);



void t2UpdateFunc_nestlike_jitter(pulsar *psr, int ipsr ,param_label label,int k, double val, double err);


double t2FitFunc_nestlike_band(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);

void t2UpdateFunc_nestlike_band(pulsar *psr, int ipsr ,param_label label,int k, double val, double err);


double t2FitFunc_nestlike_group(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);
void t2UpdateFunc_nestlike_group(pulsar *psr, int ipsr ,param_label label,int k, double val, double err);


double t2FitFunc_nestlike_shape_red(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);
double t2FitFunc_nestlike_shape_dm(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k);


void t2UpdateFunc_nestlike_shape(pulsar *psr, int ipsr ,param_label label,int k, double val, double err);
