#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>
#include <stdint.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <errno.h>


// Debug flag
#define DEBUG_LEVEL 2
typedef union __attribute__((packed))
{
    int32_t value;
    uint8_t byte[sizeof(int32_t)];
} sample_32_converter_t;    // union to help convert from bytes to integer to double

typedef struct 
{
    double sample;      // double value representation of sample
    int sampleRate;     // rate of samples per second
    int sampleIndex;    // index since begining of record
} sample_double_t;

typedef struct 
{
    double complex sample;      // double value representation of sample
    int sampleRate;     // rate of samples per second
    int sampleIndex;    // index since begining of record
} sample_complex_t;

typedef enum
{
    ALLIGNED = 0,
    MIDPOINT = 1,
    NODEBUG,
} dft_debug_t;

typedef enum
{
    NO_LOCK,
    SYMBOL_LOCK,
    PHASE_LOCK,
} timing_lock_t;

// calculate the discrete fourier transform of an array of real values but only at frequency 'k' (k cycles per windowSize samples)
//  debugFlag is to print the right debug info for different situations
//  offset is the starting position of the window in the buffer
//  windowPhase is the offset in time the window is from the start of the audio samples, probably modulo the symbol period
//  carrierPhase is the phase offset of the carrier relative to the windowPhase
double complex dft(double* buffer,
                   int windowSize,
                   int offset,
                   int windowPhase,
                   double carrierPhase,
                   int k,
                   double* rmsOut,
                   dft_debug_t debugFlag,
                   FILE* debug_fd,
                   int debug_n)
{
#if DEBUG_LEVEL <= 1
    (void)debugFlag;
    (void)debug_fd;
    (void)debug_n;
#endif

    // compute DFT (rn just one freq (k), but when we implement OFDM, then at many harmonic(orthogonal) frequencies)
    // here is a simple quadrature detector
    double complex IQ = 0;
    double RMS = 0;

    int bufferIndex;
    double phase;
    int n;  // the sample number since start of audio recording modulo window size

    for(int i = 0; i < windowSize; i++)
    {
        // recovering the time of each sample relative to the real time modulo symbol period.
        // this is important for having a consistant carrier phase between DFTs
        n = i + windowPhase;
        // This does mean there will be a constant phase misalignment between the original carrier and the IQ demod exponential, so IQ will be rotated.
        // That will be corrected for by the carrierPhase parameter by coasta's loop

        // starts at buffer[offset] and wraps around to the beginning of the buffer
        bufferIndex = (i + offset) % windowSize;

        // phase of the complex exponential
        // phasor offsets the cos and sin waves so that their phase is alligned with the real time of the samples, offset by the carrierPhase
        phase = (double)(n) * k / windowSize + carrierPhase;

        // use a summation over the symbol period to separate orthogonal components (quadrature components) at the single frequency
        double complex wave = cexp(I*2*M_PI*phase); // generate a sample from the complex exponential
        double complex value = buffer[bufferIndex] * wave;  // multiply the complex exponential by the input sample
        IQ += value;    // integrate the result over the window

        // compute RMS amplitude for equalization -- this kinda sucks. if there is a DC bias (ie, some low frequency interference)
        // it also comes through on RMS. Gotta fix that
        RMS += pow(buffer[i], 2);


        // this debug define simplifies the function a bit if debugging is disabled
    #if DEBUG_LEVEL > 1
        switch(debugFlag)
        {
            case MIDPOINT:
            {
                // debug graph outputs
                // debugging the phase allignment
                fprintf(debug_fd, "%i %i %f %i %f %i %f\n", debug_n + i, 5, buffer[bufferIndex] + 4, 6, creal(wave) + 4, 7, cimag(wave) + 4);
                // debugging the integral
                //fprintf(debug_fd, "%i %i %f %i %f %i %f %i %f\n", debug_n + i, 11, creal(value) + 6, 12, cimag(value) + 6, 13, creal(IQ) + 6, 14, cimag(IQ) + 6);
                break;
            }

            case ALLIGNED:
            {
                // debug graph outputs
                // debugging the phase allignment
                fprintf(debug_fd, "%i %i %f %i %f %i %f\n", debug_n + i, 0, buffer[bufferIndex], 1, creal(wave), 2, cimag(wave));
                break;
            }
            case NODEBUG:
                break;
        }
    #endif
    }
    // normalization factor (do I need to divide by k?)
    IQ *= sqrt(1. / windowSize);

    // complete the RMS fomula
    RMS = sqrt(1./windowSize * RMS);

#if DEBUG_LEVEL > 1
    switch(debugFlag)
    {
        case MIDPOINT:
        {
            // debug fft plot
            fprintf(debug_fd, "%i %i %f %i %f\n", debug_n, 8, creal(IQ) + 4, 9, cimag(IQ) + 4);
            fprintf(debug_fd, "%i %i %f %i %f\n", debug_n + windowSize, 8, creal(IQ) + 4, 9, cimag(IQ) + 4);
            break;
        }

        case ALLIGNED:
        {
            // debug fft plot
            fprintf(debug_fd, "%i %i %f %i %f\n", debug_n, 3, creal(IQ), 4, cimag(IQ));
            fprintf(debug_fd, "%i %i %f %i %f\n", debug_n + windowSize, 3, creal(IQ), 4, cimag(IQ));
            break;
        }
        case NODEBUG:
            break;
    }
#endif

    *rmsOut = RMS;
    return IQ;
}

typedef struct
{
    FILE* waveformPlotStdin;
    FILE* fftDebuggerStdin;
    FILE* errorPlotStdin;
    FILE* IQplotStdin;
    FILE* IQvstimeStdin;
    FILE* eyeDiagramRealStdin;
    FILE* eyeDiagramImaginaryStdin;
    FILE* filterDebugStdin;
    FILE* QAMdecoderStdin;
    union
    {
        struct
        {
            unsigned int waveformEnabled            : 1;
            unsigned int fftDebugEnabled            : 1;
            unsigned int errorPlotEnabled           : 1;
            unsigned int IQplotEnabled              : 1;
            unsigned int IQvsTimeEnabled            : 1;
            unsigned int eyeDiagramRealEnabled      : 1;
            unsigned int eyeDiagramImaginaryEnabled : 1;
            unsigned int filterDebugEnabled         : 1;
            unsigned int QAMdecoderEnabled          : 1;
        };
            unsigned long int flags;
    };
} debugPlots_t;

typedef struct
{
    int symbolPeriod;
    double **sampleBuffer;
    int k;  // the OFDM channel number
    double windowPhaseReal;

} OFDM_properties_t;

typedef struct
{
    double carrierFrequency;
    double carrierPhase;            // probably not going to use in the end.
    double complex IQsamplingTransform;         // used to transform the IQ samples to compensate for carrier phase and amplitude mismatches (ln(amp) + i*phase)
    double k;  // number of carrier cycles per sample
    double symbolPeriod;    // number of audio samples per symbol
    double selectedIQsamples[4];    // the closest samples to ideal sample time and mid symbol sample time for gardner algorithm
} QAM_properties_t;

typedef struct
{
    double *buffer;        // a pointer to an array for storing samples
    int length;             // the length of the buffer
    int insertionIndex;     // the index where new data is inserted into the circular buffer
    int phase;              // used for some functions choosing when to do periodic calculations with the buffer
    int sampleRate;
} circular_buffer_double_t;
typedef struct
{
    double complex *buffer;// a pointer to an array for storing samples
    int length;             // the length of the buffer
    int insertionIndex;     // the index where new data is inserted into the circular buffer
    int phase;              // used for some functions choosing when to do periodic calculations with the buffer
    int sampleRate;
    int n;                  // index of sample at insertionIndex
} circular_buffer_complex_t;

typedef enum
{
    RETURNED_SAMPLE,
    AWAITING_SAMPLES,
} buffered_data_return_t;


buffered_data_return_t raisedCosFilter(circular_buffer_complex_t inputSamples, sample_complex_t *outputSample, double cutoffFrequency, debugPlots_t debugPlots)
{
    //int symbolPeriod = 64; // audio samples per symbol
    int k = 1; // cycles per period
    //int filterSides = 10;    // number of symbols to either side of current symbol to filter with raised cos
    //int filterLengthSymbols = 2 * filterSides + 1;    // length of raised cos filter in IQ symbols, ie, how many IQ samples we need to generate the current symbol
    //int filterLength = filterLengthSymbols * symbolPeriod;  // length in audio samples
    //double filterCutoffFrequency = cutoffFrequency * 1.25; // cutoff frequency of the low pass filter
    double filterCutoffFrequency = cutoffFrequency * 1.25; // cutoff frequency of the low pass filter
    // array to store time series of filter data
    static double *filter;
    // array to store timeseries of IQ samples
    //static double complex *IQdata;

    static int initialized = 0;         // initialize the filter kernel first, and only once
    static int bufferPrimed = 0;        // get enough samples to start processing them.
    if(!initialized)
    {
        // initialize raised cos filter data
        //filter = malloc(sizeof(double) * filterLength); // this never gets released. so, might wanna fix that TODO
        filter = malloc(sizeof(double) * inputSamples.length);  // filter kernel = size of input buffer
        for(int i = 0; i < inputSamples.length; i++)
        {
            int filterIndex = i - inputSamples.length / 2;    // should go -filterLength/2 -> 0 -> filterLength/2
            // raised cos filter math
            double b = 0.42;    // filter parameter beta, has to do with frequency falloff and time domain fall off
            double filterValue = sin(2*M_PI * filterCutoffFrequency * filterIndex / inputSamples.sampleRate) / (2*M_PI * filterCutoffFrequency * filterIndex / inputSamples.sampleRate) * 
                (cos(2*M_PI * filterCutoffFrequency * filterIndex / inputSamples.sampleRate * b)) / (1 - pow(4 * b * filterCutoffFrequency * filterIndex / inputSamples.sampleRate, 2));
            if(!isfinite(filterValue))   // in case it's undefined, ie divide by zero case
                filterValue = sin(M_PI / 2 / b) / (M_PI / 2 / b);
            if(filterIndex == 0)
                filterValue = 1;    // the math gives a divide by zero at index 0

            filter[i] = filterValue;
        }

        // ensure initialize only runs once
        initialized = 1;
    }
    // wait for enough samples to come in
    if(!bufferPrimed)
    {
        // basically, we need to wait for future samples before we can filter the first sample. skip the remainder of this function until enough samples are collected.
        if(inputSamples.insertionIndex < inputSamples.length / 2)    // filled the buffer half way, now ready to start processing
            return AWAITING_SAMPLES;    // wait for samples

        bufferPrimed = 1;   // otherwise stop waiting and process
    }

    //int relativeSampleIndex = inputSamples.insertionIndex - inputSamples.length / 2;    // 0 is the 'current' time in the circular input buffer. negatives are into the past, positives are into the future
    //  insertionIndex should the the index of the last inserted sample, inserted by the calling function
    outputSample->sample = 0;
    outputSample->sampleRate = inputSamples.sampleRate;
    outputSample->sampleIndex = inputSamples.n - inputSamples.length / 2;   // the index is shifted by half the filter width
    for(int i = 0; i < inputSamples.length; i++)
    {
        // calculate the convoution for the sample at relative index 0

        // generate relative indexes
        int relativeSampleIndex = i - inputSamples.length / 2;     // 0 is the 'current' time in the circular input buffer. negatives are into the past, positives are into the future
        int relativeFilterIndex = -relativeSampleIndex;     // 0 is the center of the filter
        
        // generate absolute index positions
        int sampleIndex = (relativeSampleIndex + inputSamples.insertionIndex - inputSamples.length / 2) % inputSamples.length;
        if(sampleIndex < 0)
            sampleIndex += inputSamples.length; // wrap negative indexes back around to positive values
        int filterIndex = (relativeFilterIndex + inputSamples.length / 2) % inputSamples.length;
        if(filterIndex < 0)
            filterIndex += inputSamples.length; // wrap negative indexes back around to positive values

        // calculate the multiplication and sumation
        outputSample->sample += inputSamples.buffer[sampleIndex] * filter[filterIndex];
        if(debugPlots.filterDebugEnabled && relativeSampleIndex == 0)
            fprintf(debugPlots.filterDebugStdin, "%i %i %f %i %f %i %f\n", outputSample->sampleIndex, 1, creal(inputSamples.buffer[sampleIndex]), 2, cimag(inputSamples.buffer[sampleIndex]), 3, filter[outputSample->sampleIndex%inputSamples.length]);

    }
    //printf("output: n=%i, %f+%fi\n", inputSamples.n, creal(outputSample->sample), cimag(outputSample->sample));
    //outputSample->sample /= 24.;        // bit arbitrary, gotta figure out normalization factor for raised cos filter

    // print out debug info for the filter
    //fprintf(debugPlots.filterDebugStdin, "%i %i %f %i %f %i %f\n", inputSamples.n, 1, creal(inputSamples.buffer[inputSamples.insertionIndex]), 2, cimag(inputSamples.buffer[inputSamples.insertionIndex]), 3, filter[outputSample->sampleIndex%inputSamples.length]);
    if(debugPlots.filterDebugEnabled)
        fprintf(debugPlots.filterDebugStdin, "%i %i %f %i %f\n", outputSample->sampleIndex, 6, creal(outputSample->sample), 7, cimag(outputSample->sample));

    return RETURNED_SAMPLE;
}

buffered_data_return_t demodulateQAM(sample_double_t sample, QAM_properties_t QAMstate, debugPlots_t debugPlots)
{
    // currently does not handle samples after about half the filter size, instead just using them as future samples, but never getting to them. Should add a command to take care of the remainder I guess and cycle the filter function with zeros until it's done TODO
    static int initialized = 0;
    static circular_buffer_complex_t filterInputBuffer;
    static circular_buffer_complex_t timingSyncBuffer;
    if(!initialized)
    {
        // initialization
        filterInputBuffer.length = sample.sampleRate / QAMstate.carrierFrequency * 4 * (10 * 2 + 1);
        filterInputBuffer.insertionIndex = 0;
        //filterInputBuffer.buffer = malloc(sizeof(double complex) * filterInputBuffer.length);   // allocate some space for the buffer
        filterInputBuffer.buffer = calloc(filterInputBuffer.length, sizeof(double complex));   // allocate some space for the buffer, and zero it
        //printf("got %lu bytes from malloc\n", sizeof(double complex) * filterInputBuffer.length);
        if(filterInputBuffer.buffer == NULL)
            fprintf(stderr, "FilterInputBuffer failed to allocate: %s\n", strerror(errno));

        // initialize timing sync buffer
        timingSyncBuffer.length = QAMstate.symbolPeriod * 2;    // give it some extra size.
        timingSyncBuffer.insertionIndex = 0;
        timingSyncBuffer.buffer = calloc(timingSyncBuffer.length, sizeof(double complex));
        if(timingSyncBuffer.buffer == NULL)
            fprintf(stderr, "TimingSyncBuffer failed to allocate: %s\n", strerror(errno));
        // only run once
        initialized = 1;
    }

    if(debugPlots.QAMdecoderEnabled)
        fprintf(debugPlots.QAMdecoderStdin, "%i %i %f\n", sample.sampleIndex, 0, sample.sample);

    // IQ sampling
    //  multiply by sin and cos
    //      this step may have issues as the carrier frequency comes up to a quarter of the sample rate, since
    //      the multiplication step generates frequencies in the IQsample centered a  twice the original carrier 
    //      frequency. The bandwidth of the signal on that elevated carrier may alias. I should double sample rate before this step,
    //      then reduce the sample rate to fraction of the original after filtering out that high frequency stuff.
    double complex wave = cexp(I*(2*M_PI * QAMstate.carrierFrequency * sample.sampleIndex / sample.sampleRate + QAMstate.carrierPhase));
    double complex IQsample = sample.sample * wave;

    if(debugPlots.filterDebugEnabled)
    {
        fprintf(debugPlots.filterDebugStdin, "%i %i %f\n", sample.sampleIndex, 0, sample.sample);
        fprintf(debugPlots.filterDebugStdin, "%i %i %f %i %f\n", sample.sampleIndex, 4, creal(wave), 5, cimag(wave));
    }

    //  low pass filter
    //      convolution with a raised cos filter would probably do fine
    filterInputBuffer.n = sample.sampleIndex;
    filterInputBuffer.sampleRate = sample.sampleRate;
    filterInputBuffer.phase = 0;
    filterInputBuffer.buffer[filterInputBuffer.insertionIndex] = IQsample;
    sample_complex_t filteredIQsample;

    // filter the IQ samples
    buffered_data_return_t returnValue = raisedCosFilter(filterInputBuffer, &filteredIQsample, QAMstate.carrierFrequency, debugPlots);

    filterInputBuffer.insertionIndex = (filterInputBuffer.insertionIndex + 1) % filterInputBuffer.length;

    if(returnValue != RETURNED_SAMPLE) // don't continue processing unless a sample is returned
        return AWAITING_SAMPLES;

    if(debugPlots.QAMdecoderEnabled)
        fprintf(debugPlots.QAMdecoderStdin, "%i %i %f %i %f\n", filteredIQsample.sampleIndex, 5, creal(filteredIQsample.sample), 6, cimag(filteredIQsample.sample));

    // now do some timing alignemnt.
    //      I'm imagining a clock that counts up by a phase quantity scaled by the number of IQ samples elapsed. accumulates phase
    //      The rate of phase change will be adjusted by a PLL
    //      The zero crossing for the clock will be projected so that the two samples it falls between can be chosen.
    //      I'll keep an array of past IQ samples to choose from.
    //      The next filteredIQsample.sampleIndex where calculations will occur is calculated by the projection
    //          concept after PLL updates during said calculations

    // collect samples into a buffer
    timingSyncBuffer.buffer[timingSyncBuffer.insertionIndex] = filteredIQsample.sample;
    timingSyncBuffer.n = filteredIQsample.sampleIndex;
    timingSyncBuffer.insertionIndex = (timingSyncBuffer.insertionIndex + 1) % timingSyncBuffer.length;


    // Choosing samples for timing lock and symbol detection
    static double symbolSamplerAccumulatedPhase = 0;
    static double symbolSamplerPhaseRate = 0;
    static int symbolSamplerNextIndex = 0;  // next index to trigger calculations, should be just after the ideal sample time.

    if(timingSyncBuffer.n < symbolSamplerNextIndex) // check if it's too early
        return AWAITING_SAMPLES;    // it's too early, wait till the right number of IQ samples has passed.

    symbolSamplerAccumulatedPhase += symbolSamplerPhaseRate;    // apply the phase offset
    symbolSamplerNextIndex = (int)ceil(symbolSamplerAccumulatedPhase); // determine the next index to make calculations

    if(symbolSamplerAccumulatedPhase == 0)
    {
        // initialize phase rate
        symbolSamplerPhaseRate = QAMstate.symbolPeriod; // initialize the phase rate to the idealized value
        symbolSamplerNextIndex = symbolSamplerPhaseRate;    // trigger calculations one period from now
        return AWAITING_SAMPLES;
    }

    // Gardner Algorithm
    //  this can happen just once per IQ symbol, so not every audio sample
    int postIdealIndex =  timingSyncBuffer.insertionIndex - 1;
    postIdealIndex =    postIdealIndex < 0 ? timingSyncBuffer.length + postIdealIndex : postIdealIndex;     // wrap to positive
    int preIdealIndex = postIdealIndex - 1;
    preIdealIndex =     preIdealIndex < 0 ? timingSyncBuffer.length + preIdealIndex : preIdealIndex;        // wrap
    int postMidIndex =  postIdealIndex - QAMstate.symbolPeriod / 2;
    postMidIndex =      postMidIndex < 0 ? timingSyncBuffer.length + postMidIndex : postMidIndex;           // wrap
    int preMidIndex =   postMidIndex - 1;
    preMidIndex =       preMidIndex < 0 ? timingSyncBuffer.length + preMidIndex : preMidIndex;              // wrap

    //  Interpolation between samples
    double complex IQmidpoint = (timingSyncBuffer.buffer[preMidIndex] - timingSyncBuffer.buffer[postMidIndex]) * fmod(symbolSamplerAccumulatedPhase, 1) + timingSyncBuffer.buffer[preMidIndex];
    static double complex IQideal = 0;
    double complex IQlast = IQideal;
    IQideal = (timingSyncBuffer.buffer[preIdealIndex] - timingSyncBuffer.buffer[postIdealIndex]) * fmod(symbolSamplerAccumulatedPhase, 1) + timingSyncBuffer.buffer[preIdealIndex];

    // calculate error signal
    double symbolSamplerPhaseErrorEstimate = creal((IQideal - IQlast) * conj(IQmidpoint));  // this should get us a rough

    // this is shitty code, should be a function
    // PID loop
    static double integral = 0;
    integral += symbolSamplerAccumulatedPhase;
    symbolSamplerPhaseRate += integral * 0.00 + symbolSamplerPhaseErrorEstimate * 0.2;
    symbolSamplerPhaseRate = fmin(fmax(symbolSamplerPhaseRate, (double)QAMstate.symbolPeriod / 1.5), (double)QAMstate.symbolPeriod * 1.5);

    if(debugPlots.QAMdecoderEnabled)
    {
        fprintf(debugPlots.QAMdecoderStdin, "%i %i %f", timingSyncBuffer.n, 7, symbolSamplerPhaseErrorEstimate);
    }

    // update the phase rate for symbol sampler

    // Costas loop for Phase locking
    //  happens every IQ symbol
    
    //outputSample = ;  // output the determined IQ data
    return RETURNED_SAMPLE;
 }


/*
double complex demodulateOFDM(int n, double sample, OFDM_properties_t *OFDMstate, debugPlots_t debugPlots)
{
    // use the windowphase to adjust the buffer index position
    //int bufferIndex = (n + windowPhase)%SYMBOL_PERIOD;
    // Let's not use window phase to adjust index, and just pass the window phase to dft like a sane person
    int bufferIndex = n % OFDMstate->symbolPeriod;

    *OFDMstate->sampleBuffer[bufferIndex] = sample;

    static double equalizationFactor = 1;   // used by the equalization PID loop to control the overall volume
    static double RMSsamples[4] = {0};
    double RMS;


#if DEBUG_LEVEL > 0
    // debug waveform plot
    fprintf(debugPlots.waveformPlotStdin, "%i %f %f\n", n, *OFDMstate->sampleBuffer[bufferIndex], equalizationFactor);
#endif

    // array storing super sampled IQ values
    // all the IQ samples stored in these bins.
    //      0 - just before ideal midpoint
    //      1 - just after ideal midpoint
    //      2 - just before ideal
    //      3 - just after ideal
    static double complex IQsamples[4] = {0};

    // interpolated value between two closest actual IQ samples
    static double complex IQmidpoint = 0;   // between ideal sample times
    static double complex IQ = 0;   // at ideal sample time
    static double complex IQlast = 0;   // previous ideal sample time

#if DEBUG_LEVEL > 0
    // oversampling IQ values for a nice eye diagram to help debug stuff
    double complex oversampledIQ = dft(*OFDMstate->sampleBuffer, OFDMstate->symbolPeriod, bufferIndex + 1, n%OFDMstate->symbolPeriod, 0., OFDMstate->k, &RMS, NODEBUG, NULL, n);
    fprintf(debugPlots.eyeDiagramRealStdin, "%f %i %f\n", fmod((double)n / OFDMstate->symbolPeriod, 4), n / 4 / OFDMstate->symbolPeriod, creal(oversampledIQ));
    fprintf(debugPlots.eyeDiagramImaginaryStdin, "%f %i %f\n", fmod((double)n / OFDMstate->symbolPeriod, 4), n / 4 / OFDMstate->symbolPeriod, cimag(oversampledIQ));
#endif

    // deciding when to take IQ samples
    // take a sample right between the ideal samples which are taken at windowPhase
    // if we are half full on the buffer, take an intermidiate IQ sample, for timing sync later
    // compute DFT (rn just one freq, but when we implement OFDM, then at many harmonic(orthogonal) frequencies)
    if(bufferIndex == ((int)floor(OFDMstate->windowPhaseReal) + OFDMstate->symbolPeriod / 2 - 1) % OFDMstate->symbolPeriod)    // just before real window phase offset midpoint
        IQsamples[0] = dft(*OFDMstate->sampleBuffer,
                           OFDMstate->symbolPeriod,
                           bufferIndex + 1,
                           (n+1)%OFDMstate->symbolPeriod,
                           0.,
                           OFDMstate->k,
                           &RMSsamples[0],
                           MIDPOINT,
                           debugPlots.fftDebuggerStdin,
                           n);

    if(bufferIndex == ((int)ceil(OFDMstate->windowPhaseReal) + OFDMstate->symbolPeriod / 2 - 1) % OFDMstate->symbolPeriod)     // just after real window phase offset midpoint;
        IQsamples[1] = dft(*OFDMstate->sampleBuffer,
                           OFDMstate->symbolPeriod,
                           bufferIndex + 1,
                           (n+1)%OFDMstate->symbolPeriod,
                           0.,
                           OFDMstate->k,
                           &RMSsamples[1],
                           MIDPOINT,
                           debugPlots.fftDebuggerStdin,
                           n);

    if(bufferIndex == ((int)floor(OFDMstate->windowPhaseReal) + OFDMstate->symbolPeriod - 1) % OFDMstate->symbolPeriod)        // just before real window phase offset
        IQsamples[2] = dft(*OFDMstate->sampleBuffer,
                           OFDMstate->symbolPeriod,
                           bufferIndex + 1,
                           (n+1)%OFDMstate->symbolPeriod,
                           0.,
                           OFDMstate->k,
                           &RMSsamples[2],
                           ALLIGNED,
                           debugPlots.fftDebuggerStdin,
                           n);

    if(bufferIndex == ((int)ceil(OFDMstate->windowPhaseReal) + OFDMstate->symbolPeriod - 1) % OFDMstate->symbolPeriod)         // just after real window phase offset
        IQsamples[3] = dft(*OFDMstate->sampleBuffer,
                           OFDMstate->symbolPeriod,
                           bufferIndex + 1,
                           (n+1)%OFDMstate->symbolPeriod,
                           0.,
                           OFDMstate->k,
                           &RMSsamples[3],
                           ALLIGNED,
                           debugPlots.fftDebuggerStdin,
                           n);

    // I added another condition to help debounce. sometimes it takes many samples in a row due to changing window offset
    if((bufferIndex == ((int)ceil(OFDMstate->windowPhaseReal) + OFDMstate->symbolPeriod) % OFDMstate->symbolPeriod) && (n - tookSampleAt > OFDMstate->symbolPeriod / 2))  // just after real window phase offset, do the calculations that must be done once per symbol recieved
    {
        // interpolate IQ and IQmidpoint from actual IQ samples in IQsamples
        // just doing linear interpolation
        //      (slope) * interp_X + initial
        //      (final - initial) * interp_X + initial
        IQmidpoint =    (IQsamples[1] - IQsamples[0]) * fmod(OFDMstate->windowPhaseReal, 1) + IQsamples[0];
        IQlast =        IQ;
        IQ =            (IQsamples[3] - IQsamples[2]) * fmod(OFDMstate->windowPhaseReal, 1) + IQsamples[2];

        tookSampleAt = n;

        // now I'm doing a bunch of stuff that happens every IQ sample. This all happens in the timespan of a single audio sample, which is 1/symbolPeriod of the time between IQ samples that could be used, but whatever

        // averaging filter for the equalizer
        static double rmsaverageWindow[44100 * 2 / OFDMstate->symbolPeriod] = {0};
        int rmsaverageSize = 44100 * 2 /OFDMstate->symbolPeriod;
        int rmsaverageIndex = (n / OFDMstate->symbolPeriod) % rmsaverageSize;

        RMS = 0;
        for(int i = 0; i < 4; i++)
        {
            RMS += RMSsamples[i];
        }
        RMS /= 4;
        rmsaverageWindow[rmsaverageIndex] = RMS;
        double rmsaverage = 0;

        for(int i = 0; i < rmsaverageSize; i++)
        {
            rmsaverage += rmsaverageWindow[i];
        }
        rmsaverage /= rmsaverageSize;
        equalizationFactor += (1. / sqrt(2) - rmsaverage) * 0.002;
        // PID for the equalizer, just proportional. with max
        equalizationFactor = fmax(fmin(equalizationFactor, 10000), -10000);
        fprintf(debugPlots.waveformPlotStdin, "%i %f %f %f\n", n, *OFDMstate->sampleBuffer[bufferIndex], equalizationFactor, rmsaverage);

#if DEBUG_LEVEL > 0
        equalizationFactor = 1;     // equalization factor needs to ignore low frequency signals that give DC offset
        //equalizationFactor = 1./0.007;
#endif

        // try to get a phase lock, symbol time lock, frequency lock, and equalization
        // calculate the error signal
        // Gardner Algorithm: Real Part( derivitive of IQ times conjugate of IQ)
        // boobs-alexis
        // basically, it's trying to estimate the error of zero crossings
        double phaseOffsetEstimate = creal((IQ - IQlast) * conj(IQmidpoint));

#if DEBUG_LEVEL > 1
        // imprint the phase offset estimate on the fftDebugger plot to help understand the relationship between the DFT and symbol sync algo
        fprintf(debugPlots.fftDebuggerStdin, "%i %i %f %i %f\n", n, 10, phaseOffsetEstimate + 2, 15, (double)(n % OFDMstate->symbolPeriod) / OFDMstate->symbolPeriod + 2);
#endif

        // state of the timing lock to modulate PID params
        static timing_lock_t lockstate = NO_LOCK;

        // Process Variable (PV, ie phase estimate) filter
        // rolling average of phase offset estimate
        // this may need to be adjusted based on the state of symbol and phase lock achieved
        static double averageWindow[64] = {0};
        const int averageWindowArrayLength = 64;
        static int averageSize;    // adjusted based on lock state
        int averageIndex = (n / OFDMstate->symbolPeriod) % averageWindowArrayLength;
        switch(lockstate)
        {
            case NO_LOCK:
                // smaller averaging window to achieve
                //  faster error signal response
                //  for more aggressive PID tune
                averageSize = 2;
                break;
            case SYMBOL_LOCK:
            case PHASE_LOCK:
                // longer average window to help average out symbol transitions that are not zero crossings
                //  PID tune must be shittier
                averageSize = 64;
                break;
        }

        // add a sample and take the average of last averageSize samples
        averageWindow[averageIndex] = phaseOffsetEstimate;
        double average = 0;
        // step backwards from current index to averageSize indexed back
        for(int i = averageIndex; i > averageIndex - averageSize; i--)
        {
            if(i < 0)
                average += averageWindow[averageWindowArrayLength + i];
            average += averageWindow[i];
        }
        average /= averageSize;

        // PID loop for symbol timing, ie aligning the fft window to the symbol transitions
        static double P_gain;
        static double I_gain;
        static double D_gain;
        switch(lockstate)
        {
            case NO_LOCK:
                P_gain = 2.5;
                I_gain = 0.003;
                D_gain = 0;
                break;
            case SYMBOL_LOCK:
            case PHASE_LOCK:
                P_gain = 0.1;
                I_gain = 0.001;
                D_gain = 0;
                break;
        }

        double error = 0 - average;
        //error *= -1;
        static double errorIntegral = 0;
        errorIntegral += error;

        static double lastError = 0;
        double errorDerivative = error - lastError;
        lastError = error;

        // this may need to be adjusted based on the state of symbol and phase lock achieved
        double phaseAdjustment = errorDerivative * D_gain + errorIntegral * I_gain + error * P_gain;

        OFDMstate->windowPhaseReal += phaseAdjustment;
        OFDMstate->windowPhaseReal = fmod(OFDMstate->windowPhaseReal + phaseAdjustment, OFDMstate->symbolPeriod);
        OFDMstate->windowPhaseReal = OFDMstate->windowPhaseReal < 0 ? OFDMstate->symbolPeriod + OFDMstate->windowPhaseReal : OFDMstate->windowPhaseReal;
        windowPhase = (int)round(OFDMstate->windowPhaseReal) % OFDMstate->symbolPeriod; // quantize the real window phase
        //windowPhase = windowPhase < 0 ? OFDMstate->symbolPeriod + windowPhase : windowPhase;
#if DEBUG_LEVEL > 0
        // some options to overwrite the window phase given by the PID controller
        //windowPhase = (n * 2 / 2000) % OFDMstate->symbolPeriod;
        //windowPhase = (n * 4 * 2/ (OFDMstate->symbolPeriod * 2000) ) % 4 * OFDMstate->symbolPeriod / 4;
        //windowPhase = OFDMstate->symbolPeriod / 4;
        windowPhase = 0;
        OFDMstate->windowPhaseReal = 0.0001;
        //windowPhase = (int)((windowPhase + phaseAdjustment) < 0 ? (OFDMstate->symbolPeriod - windowPhase + phaseAdjustment) : (windowPhase + phaseAdjustment)) % OFDMstate->symbolPeriod;
#endif

        // extract the frequencies to be decoded
        // add the relevant IQ values to an array
        // plot the IQ and a tail with some kind of persistance to give an animated output
        //fprintf(debugPlots.IQplotStdin, "%f %f\n", I, Q);
        //printf("%f %f\n", I, Q);

        // plot with a new color for each window phase
        fprintf(debugPlots.IQplotStdin, "%f %i %f\n", creal(IQ), windowPhase, cimag(IQ));
        //fprintf(debugPlots.IQplotStdin, "%f %i %f\n", creal(IQmidpoint), windowPhase + OFDMstate->symbolPeriod, cimag(IQmidpoint));
        fprintf(debugPlots.errorPlotStdin, "%i, %f %f %f %f %f\n", n / OFDMstate->symbolPeriod, -phaseOffsetEstimate, error, phaseAdjustment, (double)windowPhase / OFDMstate->symbolPeriod, OFDMstate->windowPhaseReal / OFDMstate->symbolPeriod);
        //fprintf(debugPlots.IQvstimeStdin, "%i, %f, %f, %f, %f, %f, %f\n", n / OFDMstate->symbolPeriod % (2*3), creal(IQ), cimag(IQ), creal(IQlast), cimag(IQlast), creal(IQmidpoint), cimag(IQmidpoint));
        fprintf(debugPlots.IQvstimeStdin, "%f, %f, %f\n", (n / OFDMstate->symbolPeriod) % (2*3) + 0., creal(IQ), cimag(IQ));
        fprintf(debugPlots.IQvstimeStdin, "%f, %f, %f\n", (n / OFDMstate->symbolPeriod) % (2*3) + 0.5, creal(IQmidpoint), cimag(IQmidpoint));
    }
}
*/


int main(void)
{
    // length of each symbol in samples, and so the fft window.
    // This is the period of time all the orthogonal symbols will be integrated over
#define SYMBOL_PERIOD 64
    // the OFDM channel number, how many cycles per symbol
    int k = 4;
    // buffer is the length of the symbol period, so that symbols are orthogonal
    double sampleBuffer[SYMBOL_PERIOD] = {0.0};

    //int guardPeriod = 4./1000 * 44100;      // 4ms guard period for room echos
    //int guardPeriod = 0;
    //int totalPeriod = SYMBOL_PERIOD + guardPeriod;
    sample_32_converter_t sampleConvert;
    sample_double_t sample;

    // the offset of sample buffer or fftwindow in time, to get symbol time sync
    int windowPhase = 0;
    //double windowPhaseReal = (double)rand() / RAND_MAX * SYMBOL_PERIOD; // the floating point window phase, to be quantized into windowPhase
    double windowPhaseReal = 0;

    // to help the if statements that choose when to take an OFDM DFT sample not get stuck in a loop when the window offset changes
    int tookSampleAt = 0;

    // prepare file descriptors for output files and streams
    int retval = 0;

    debugPlots_t debugPlots;    // holds all the file descriptors for plots
    char *iq_plot_buffer = NULL;

    // process arguments
    //  none right now

    // for live plotting, pipe to feedgnuplot
    const char *plot =
        "feedgnuplot "
        "--domain --lines --points "
        "--title \"Raw signal Time domain\" "
        "--xlabel \"Time (microphone sample #\" --ylabel \"value\" "
        "--legend 0 \"Signal\" "
        "--legend 1 \"equalization factor\" "
    ;
    //const char *plot = 
        //"hexdump -C "
        //"tee testoutput.txt"
    //;

    // using it to plot the time domain signal
    debugPlots.waveformPlotStdin = popen(plot, "w");

    if (debugPlots.waveformPlotStdin == NULL)
    {
        fprintf(stderr, "Failed to create waveform plot: %s\n", strerror(errno));
        retval = 1;
        goto exit;
    }

    const char *debugPlot =
        "feedgnuplot "
        "--domain --dataid --lines --points "
        "--title \"FFT debuger graph\" "
        "--xlabel \"Time (microphone sample #\" --ylabel \"value\" "
        "--legend 0 \"input samples\" "
        "--legend 1 \"fft real\" "
        "--legend 2 \"fft imaginary\" "
        "--legend 3 \"I decision\" "
        "--legend 4 \"Q decision\" "
        "--legend 5 \"input samples mid\" "
        "--legend 6 \"fft real mid\" "
        "--legend 7 \"fft imaginary mid\" "
        "--legend 8 \"I decision mid\" "
        "--legend 9 \"Q decision mid\" "
        "--legend 11 \"samp*real\" "
        "--legend 12 \"samp*imag\" "
        "--legend 13 \"integral real\" "
        "--legend 14 \"integral imag\" "
        "--legend 10 \"phase error signal\" "
        "--legend 15 \"window phase\" "
    ;

    // using it to plot the time domain signal
    debugPlots.fftDebuggerStdin = popen(debugPlot, "w");

    if (debugPlots.fftDebuggerStdin == NULL)
    {
        fprintf(stderr, "Failed to create fft debug plot: %s\n", strerror(errno));
        retval = 2;
        goto exit;
    }

    const char *errorPlot =
        "feedgnuplot "
        "--domain --lines --points "
        "--title \"Time Domain error signal\" "
        "--legend 0 \"estimated phase offset\" "
        "--legend 1 \"rolling average error signal\" "
        "--legend 2 \"PI filtered signal\" "
        "--legend 3 \"fft window phase shift\" "
        "--legend 4 \"real target phase offset\" "
    ;

    // using it to plot the time domain signal
    debugPlots.errorPlotStdin = popen(errorPlot, "w");

    if (debugPlots.errorPlotStdin == NULL)
    {
        fprintf(stderr, "Failed to create error plot: %s\n", strerror(errno));
        retval = 3;
        goto exit;
    }

    //FILE* IQplotStdin = popen("feedgnuplot --domain --points --title \"IQ plot\" --unset key", "w");

    // Allocate buffer for call to feedgnuplot for IQ plot
    iq_plot_buffer = malloc(120 + (50 * SYMBOL_PERIOD));

    // Check if allocation failed
    if (iq_plot_buffer == NULL)
    {
        fprintf(stderr, "Failed to allocate buffer for IQ plot: %s\n", strerror(errno));
        retval = 4;
        goto exit;
    }

    int stringLength = 0;

    stringLength += snprintf(iq_plot_buffer, 120,
        "feedgnuplot "
        "--dataid --domain --points --maxcurves %i "
        "--title \"IQ plot\" "
        "--xlabel \"I\" --ylabel \"Q\" ",
        SYMBOL_PERIOD * 2 + 1
    );

    if (stringLength < 0)
    {
        fprintf(stderr, "Printing iq plot buffer failed: %s\n", strerror(errno));
        retval = 5;
        goto exit;
    }
    else if (stringLength == 120)
    {
        fprintf(stderr, "Printing iq plot buffer failed: truncated");
        retval = 5;
        goto exit;
    }

    for(int i = 0; i < SYMBOL_PERIOD; i++)
    {
        stringLength += snprintf(iq_plot_buffer + stringLength, 50,
            "--legend %i \"Phase offset %i samples\" ",
            i,
            i
        );

        if (stringLength < 0)
        {
            fprintf(stderr, "Printing iq plot buffer failed: %s\n", strerror(errno));
            retval = 5;
            goto exit;
        }
        else if (stringLength == 50)
        {
            fprintf(stderr, "Printing iq plot buffer failed: truncated");
            retval = 5;
            goto exit;
        }
    }

    debugPlots.IQplotStdin = popen(iq_plot_buffer, "w");

    // Free buffer, even if popen failed
    free(iq_plot_buffer);
    iq_plot_buffer = NULL;

    if (debugPlots.IQplotStdin == NULL)
    {
        fprintf(stderr, "Failed to create IQ plot: %s\n", strerror(errno));
        retval = 6;
        goto exit;
    }

    // for ploting IQ values over time to hopefully obtain an error function
    char eyeDiagramPlotReal[300] = {0};
    sprintf(
        eyeDiagramPlotReal,
        "feedgnuplot "
        "--domain --dataid --lines --points --maxcurves %i "
        "--title \"Eye Diagram, Real part\" "
        "--xlabel \"Time (Audio sample #)\" --ylabel \"I\" ",
        3000);
    debugPlots.eyeDiagramRealStdin = popen(eyeDiagramPlotReal, "w");

    char eyeDiagramPlotImaginary[300] = {0};
    sprintf(
        eyeDiagramPlotImaginary,
        "feedgnuplot "
        "--domain --dataid --lines --points --maxcurves %i "
        "--title \"Eye Diagram, Imaginary part\" "
        "--xlabel \"Time (Audio sample #)\" --ylabel \"Q\" ",
        3000);
    debugPlots.eyeDiagramImaginaryStdin = popen(eyeDiagramPlotImaginary, "w");

    // using it to plot the time domain signal
    debugPlots.IQvstimeStdin = popen(eyeDiagramPlotReal, "w");

    if (debugPlots.IQvstimeStdin == NULL)
    {
        fprintf(stderr, "Failed to create IQ vs Time plot: %s\n", strerror(errno));
        retval = 7;
        goto exit;
    }

    char *filterDebugPlot = 
        "feedgnuplot "
        "--domain --dataid --lines --points "
        "--title \"Filter debugging plot\" "
        "--xlabel \"Time (Audio sample #)\" --ylabel \"value\" "
        "--legend 0 \"Original Audio samples\" "
        "--legend 1 \"Original I\" "
        "--legend 2 \"Original Q\" "
        "--legend 3 \"Filter array\" "
        "--legend 4 \"IQ sampler internal exponential real part\" "
        "--legend 5 \"IQ sampler internal exponential imaginary part\" "
        "--legend 6 \"Filtered I\" "
        "--legend 7 \"Filtered Q\" "

    ;
    debugPlots.filterDebugStdin = popen(filterDebugPlot, "w");
    if(debugPlots.filterDebugStdin == NULL)
    {
        fprintf(stderr, "Failed to create filter debug plot: %s\n", strerror(errno));
        retval = 8;
        goto exit;
    }

    char *QAMdemodulatePlot = 
        "feedgnuplot "
        "--domain --dataid --lines --points "
        "--title \"QAM demodulation debug plot\" "
        "--xlabel \"Time (sample index)\" --ylabel \"value\" "
        "--legend 0 \"Original audio samples\" "
        "--legend 1 \"IQ sampler internal exponential real part\" "
        "--legend 2 \"IQ sampler internal exponential imaginary part\" "
        "--legend 3 \"I prefilter\" "
        "--legend 4 \"Q prefilter\" "
        "--legend 5 \"I\" "
        "--legend 6 \"Q\" "
        "--legend 7 \"Phase error signal\" "
    ;
    debugPlots.QAMdecoderStdin = popen(QAMdemodulatePlot, "w");
    if(debugPlots.QAMdecoderStdin == NULL)
    {
        fprintf(stderr, "Failed to create QAM decoder debug plot: %s\n", strerror(errno));
        retval = 9;
        goto exit;
    }

    debugPlots.flags = 0;   // reset all the flags
    
    // set some debug flags
    debugPlots.waveformEnabled = 1;
    debugPlots.QAMdecoderEnabled = 1;
    debugPlots.filterDebugEnabled = 1;


    // while there is data to recieve, not end of file -> right now just a fixed number of 2000
    for(int audioSampleIndex = 0; audioSampleIndex < SYMBOL_PERIOD * 600; audioSampleIndex++)
    {
        // recieve data on stdin, signed 32bit integer
        for(size_t i = 0; i < sizeof(sampleConvert.value); i++)
        {
            // get the bytes, little endian, of signed 32bit integer, using the struct to do type punning
            sampleConvert.byte[i] = getchar();
        }


        // convert to a double ranged from -1 to 1
        sample.sample = (double)sampleConvert.value / INT32_MAX;
        sample.sampleRate = 44100;
        sample.sampleIndex = audioSampleIndex;

        if(debugPlots.waveformEnabled)
            fprintf(debugPlots.waveformPlotStdin, "%i %f\n", sample.sampleIndex, sample.sample);

        QAM_properties_t QAMstate;
        QAMstate.carrierFrequency = (double)sample.sampleRate / 64;
        QAMstate.carrierPhase = 0;
        QAMstate.k = 1;
        QAMstate.symbolPeriod = (int)(sample.sampleRate / QAMstate.carrierFrequency);
        demodulateQAM(sample, QAMstate, debugPlots);


        //OFDM_properties_t OFDMdemodulateState = {SYMBOL_PERIOD, &sampleBuffer, k};
        //demodulateOFDM(n, sampleDouble, &OFDMdemodulateState, &debugPlots);

    }


exit:
    if (iq_plot_buffer != NULL)
    {
        free(iq_plot_buffer);
        iq_plot_buffer = NULL;
    }

    // close the streams, allowing the plots to render and open a window, then wait for them to terminate in separate threads.
    if((debugPlots.IQvstimeStdin != NULL) && (fork() == 0))
    {
        pclose(debugPlots.IQvstimeStdin);
        return 0;
    }

    if((debugPlots.IQplotStdin != NULL) && (fork() == 0))
    {
        pclose(debugPlots.IQplotStdin);
        return 0;
    }

    if((debugPlots.errorPlotStdin != NULL) && (fork() == 0))
    {
        pclose(debugPlots.errorPlotStdin);
        return 0;
    }

    if((debugPlots.waveformPlotStdin != NULL) && (fork() == 0))
    {
        pclose(debugPlots.waveformPlotStdin);
        return 0;
    }

    if((debugPlots.fftDebuggerStdin != NULL) && (fork() == 0))
    {
        pclose(debugPlots.fftDebuggerStdin);
        return 0;
    }
    if((debugPlots.eyeDiagramRealStdin != NULL) && (fork() == 0))
    {
        pclose(debugPlots.eyeDiagramRealStdin);
        return 0;
    }
    if((debugPlots.eyeDiagramImaginaryStdin != NULL) && (fork() == 0))
    {
        pclose(debugPlots.eyeDiagramImaginaryStdin);
        return 0;
    }
    if((debugPlots.filterDebugStdin != NULL) && (fork() == 0))
    {
        pclose(debugPlots.filterDebugStdin);
        return 0;
    }

    // after all the data is recieved, generate a plot of IQ
    return retval;
}
