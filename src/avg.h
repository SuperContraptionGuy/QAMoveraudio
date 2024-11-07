#pragma once

typedef struct
{
	double *samples;
	size_t numSamples;
	size_t position;
	size_t count;
} avg_t;

int avgInit(avg_t *avg, size_t length);
int avgReInit(avg_t *avg, size_t length);
int avgReInitCopy(avg_t *avg, size_t length);
void avgDestroy(avg_t *avg);
int avgAppend(avg_t *avg, double value);
double avgGet(avg_t *avg);
