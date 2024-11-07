#include <stdlib.h>

#include <string.h>
#include <math.h>

#include "avg.h"

static int avgCheckLength(avg_t *avg, size_t length)
{
	// Check average pointer
	if (avg == NULL)
	{
		return -1;
	}

	// Check length is at least 2
	if (length < 2)
	{
		return -2;
	}

	// Check we have valid samples
	if (avg->samples == NULL)
	{
		return -3;
	}

	return 0;
}

static int avgCheck(avg_t *avg)
{
	// Check average pointer
	if (avg == NULL)
	{
		return -1;
	}

	// Check length is at least 2
	if (avg->numSamples < 2)
	{
		return -2;
	}

	// Check we have valid samples
	if (avg->samples == NULL)
	{
		return -3;
	}

	return 0;
}

int avgInit(avg_t *avg, size_t length)
{
	int ret = avgCheckLength(avg, length);

	// Samples must be NULL, otherwise it's already allocated
	if (ret != -3)
	{
		return ret;
	}

	// Allocate new sample buffer
	avg->samples = (double *)malloc(sizeof(double) * length);

	// Check allocation
	if (avg->samples == NULL)
	{
		return -4;
	}

	// Initialize samples to zero
	memset(avg->samples, 0, sizeof(double) * length);

	// Set sample length
	avg->numSamples = length;
	// Move to start
	avg->position = 0;
	avg->count = 0;

	return 0;
}

int avgReInit(avg_t *avg, size_t length)
{
	int ret = avgCheckLength(avg, length);

	// Samples allowed to be NULL
	if ((ret != 0) && (ret != -3))
	{
		return ret;
	}

	// If sample buffer is allocated and the length matches, nothing to do
	if ((avg->samples != NULL) && (length == avg->numSamples))
	{
		return 0;
	}

	// Sample buffer is null, just call Init
	if (avg->samples == NULL)
	{
		return avgInit(avg, length);
	}

	// Reallocate buffer with new length
	double *ptr = (double *)realloc(avg->samples, sizeof(double) * length);

	// Check allocation
	if (ptr == NULL)
	{
		return -3;
	}

	// Set updated sample buffer
	avg->samples = ptr;

	// Initialize samples to zero
	memset(avg->samples, 0, sizeof(double) * length);

	// Set sample length
	avg->numSamples = length;
	// Move to start
	avg->position = 0;
	avg->count = 0;

	return 0;
}

int avgReInitCopy(avg_t *avg, size_t length)
{
	int ret = avgCheckLength(avg, length);

	// Samples allowed to be NULL
	if ((ret != 0) && (ret != -3))
	{
		return ret;
	}

	// If sample buffer is allocated and the length matches, nothing to do
	if ((avg->samples != NULL) && (length == avg->numSamples))
	{
		return 0;
	}

	// Sample buffer is null, just call Init
	if (avg->samples == NULL)
	{
		return avgInit(avg, length);
	}

	// Allocate buffer with new length
	double *ptr = (double *)malloc(sizeof(double) * length);

	// Check allocation
	if (ptr == NULL)
	{
		return -3;
	}

	// Determine number of samples to copy
	size_t samplesToCopy = length;

	// avg->numSamples < length => Growing the buffer
	// avg->numSamples > length => Shrinking the buffer
	if (samplesToCopy > avg->numSamples)
	{
		// We're shrinking the buffer, wrap to avg->numSamples
		samplesToCopy = avg->numSamples;
	}

	// Get index of latest sample
	size_t index = avg->position - 1;

	// Copy samples from the old buffer
	for (size_t i = 0; i < samplesToCopy; i++)
	{
		// Copy to destination forwards...
		ptr[i] = avg->samples[index];

		if (index == 0)
		{
			// Wrap around to the end
			index = avg->numSamples - 1;
		}
		else
		{
			// from source in reverse
			index -= 1;
		}
	}

	// Free the old buffer
	free(avg->samples);
	// Set updated sample buffer
	avg->samples = ptr;

	// Set sample length
	avg->numSamples = length;

	// When shrinking:
	//   avg->numSamples > length -> possible that avg->position > length
	// When growing:
	//   avg->numSamples < length -> guaranteed avg->position < length
	if (avg->position > samplesToCopy)
	{
		// Wrap position to the last copied sample when shrinking
		avg->position = samplesToCopy;
	}

	// Same with position,
	if (avg->count > samplesToCopy)
	{
		// Wrap count to new size when shrinking
		avg->count = samplesToCopy;
	}

	return 0;
}

void avgDestroy(avg_t *avg)
{
	int ret = avgCheck(avg);

	// Sample buffer allowed to be NULL
	if ((ret != 0) && (ret != -3))
	{
		return;
	}

	if (avg->samples != NULL)
	{
		free(avg->samples);
		avg->samples = NULL;
	}

	avg->numSamples = 0;
	avg->position = 0;
	avg->count = 0;
}

int avgAppend(avg_t *avg, double value)
{
	int ret = avgCheck(avg);

	if (ret != 0)
	{
		return ret;
	}

	// Implement circular buffer
	if (avg->position >= avg->numSamples)
	{
		// Wrap position to 0
		avg->position = 0;
	}
	else
	{
		// Otherwise, increment position
		avg->position += 1;
	}

	// Keep track of the number of samples written so far up to numSamples
	//   This is needed when calling avgGet() when n < numSamples to properly
	//   compute sum(samples) / n
	if (avg->count < avg->numSamples)
	{
		avg->count += 1;
	}

	// Update value in sample buffer
	avg->samples[avg->position] = value;
	// Move to the next position
	avg->position += 1;

	return 0;
}

double avgGet(avg_t *avg)
{
	int ret = avgCheck(avg);

	if (ret != 0)
	{
		return NAN;
	}

	if (avg->count == 0)
	{
		return 0.0;
	}

	double accumulator = 0.0;

	// Get sum of collected samples
	for (size_t i = 0; i < avg->count; i++)
	{
		accumulator += avg->samples[i];
	}

	// Compute sum(samples) / n
	double average = accumulator / avg->count;

	return average;
}
