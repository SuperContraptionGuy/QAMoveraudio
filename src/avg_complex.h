#pragma once

typedef struct
{
	double complex *samples;
	size_t numSamples;
	size_t position;
	size_t count;
} avg_complex_t;

int avgCplxInit(avg_complex_t *avg, size_t length);
int avgCplxReInit(avg_complex_t *avg, size_t length);
void avgCplxDestroy(avg_complex_t *avg);
int avgCplxAppend(avg_complex_t *avg, double complex value);
double complex avgCplxGet(avg_complex_t *avg);
