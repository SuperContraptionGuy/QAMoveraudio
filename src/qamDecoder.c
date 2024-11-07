#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>
#include <stdint.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <errno.h>

#include "avg.h"

// length of each symbol in samples, and so the fft window.
// This is the period of time all the orthogonal symbols will be integrated over
#define SYMBOL_PERIOD 64

// Debug flag
#define DEBUG_LEVEL 1

typedef struct __attribute__((packed))
{
    union
    {
        int32_t value;
        uint8_t byte[sizeof(int32_t)];
    };
} sample32_t;

typedef enum
{
    ALLIGNED = 0,
    MIDPOINT = 1,
    NODEBUG,
} dft_debug_flag_t;

typedef struct
{
    dft_debug_flag_t flag;
    FILE *fd;
} dft_debug_t;

typedef enum
{
    NO_LOCK,
    SYMBOL_LOCK,
    PHASE_LOCK,
} timing_lock_t;

typedef struct
{
    double samples[SYMBOL_PERIOD];
    int windowSize;
    int offset;
    int windowPhase;
    double carrierPhase;
    int k;
    double rms;
    int n;
} dft_t;

typedef struct
{
    // interpolated value between two closest actual IQ samples
    double complex midpoint;
    double complex last;
    double complex iq;
    int tookSampleAt;
} iq_sample_t;

// calculate the discrete fourier transform of an array of real values but only at frequency 'k' (k cycles per windowSize samples)
//  debugFlag is to print the right debug info for different situations
//  offset is the starting position of the window in the buffer
//  windowPhase is the offset in time the window is from the start of the audio samples, probably modulo the symbol period
//  carrierPhase is the phase offset of the carrier relative to the windowPhase
double complex dft(double* buffer, int windowSize, int offset, int windowPhase, double carrierPhase, int k, double* rmsOut, dft_debug_flag_t debugFlag, FILE* debug_fd, int debug_n)
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
        double complex wave = cexp(I * 2 * M_PI * phase); // generate a sample from the complex exponential
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
    IQ *= sqrt(1.0 / windowSize);

    // complete the RMS fomula
    RMS = sqrt(1.0 / SYMBOL_PERIOD * RMS);

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

double complex dft2(dft_t *config, dft_debug_t *debug)
{
#if DEBUG_LEVEL <= 1
    (void)debug;
#endif

    // compute DFT (rn just one freq (k), but when we implement OFDM, then at many harmonic(orthogonal) frequencies)
    // here is a simple quadrature detector
    double complex IQ = 0;
    // double RMS = 0;

    int bufferIndex;
    double phase;
    int n;  // the sample number since start of audio recording modulo window size

    for(int i = 0; i < config->windowSize; i++)
    {
        // recovering the time of each sample relative to the real time modulo symbol period.
        // this is important for having a consistant carrier phase between DFTs
        n = i + config->windowPhase;
        // This does mean there will be a constant phase misalignment between the original carrier and the IQ demod exponential, so IQ will be rotated.
        // That will be corrected for by the carrierPhase parameter by coasta's loop

        // starts at buffer[offset] and wraps around to the beginning of the buffer
        bufferIndex = (i + config->offset) % config->windowSize;

        // phase of the complex exponential
        // phasor offsets the cos and sin waves so that their phase is alligned with the real time of the samples, offset by the carrierPhase
        phase = (double)(n) * config->k / config->windowSize + config->carrierPhase;

        // use a summation over the symbol period to separate orthogonal components (quadrature components) at the single frequency
        double complex wave = cexp(I * 2 * M_PI * phase); // generate a sample from the complex exponential
        double complex value = config->samples[bufferIndex] * wave;  // multiply the complex exponential by the input sample
        IQ += value;    // integrate the result over the window

        // compute RMS amplitude for equalization -- this kinda sucks. if there is a DC bias (ie, some low frequency interference)
        // it also comes through on RMS. Gotta fix that
        config->rms += pow(config->samples[i], 2);

        // this debug define simplifies the function a bit if debugging is disabled
    #if DEBUG_LEVEL > 1
        switch(debug->flag)
        {
            case MIDPOINT:
            {
                // debug graph outputs
                // debugging the phase allignment
                fprintf(debug->fd,
                    "%i %i %f %i %f %i %f\n",
                    config->n + i,
                    5,
                    config->samples[bufferIndex] + 4,
                    6,
                    creal(wave) + 4,
                    7, cimag(wave) + 4
                );

                // debugging the integral
                // fprintf(debug->fd,
                //     "%i %i %f %i %f %i %f %i %f\n",
                //     debug->n + i,
                //     11,
                //     creal(value) + 6,
                //     12,
                //     cimag(value) + 6,
                //     13,
                //     creal(IQ) + 6,
                //     14,
                //     cimag(IQ) + 6
                // );

                break;
            }

            case ALLIGNED:
            {
                // debug graph outputs
                // debugging the phase allignment
                fprintf(debug->fd,
                    "%i %i %f %i %f %i %f\n",
                    config->n + i,
                    0,
                    config->samples[bufferIndex],
                    1,
                    creal(wave),
                    2,
                    cimag(wave)
                );

                break;
            }
            case NODEBUG:
                break;
        }
    #endif
    }

    // normalization factor (do I need to divide by k?)
    IQ *= sqrt(1.0 / config->windowSize);

    // complete the RMS fomula
    config->rms = sqrt(1.0 / SYMBOL_PERIOD * config->rms);

#if DEBUG_LEVEL > 1
    switch(debug->flag)
    {
        case MIDPOINT:
        {
            // debug fft plot
            fprintf(debug->fd,
                "%i %i %f %i %f\n",
                config->n,
                8,
                creal(IQ) + 4,
                9, cimag(IQ) + 4
            );

            fprintf(debug->fd,
                "%i %i %f %i %f\n",
                config->n + config->windowSize,
                8,
                creal(IQ) + 4,
                9,
                cimag(IQ) + 4
            );

            break;
        }

        case ALLIGNED:
        {
            // debug fft plot
            fprintf(debug->fd,
                "%i %i %f %i %f\n",
                config->n,
                3,
                creal(IQ),
                4,
                cimag(IQ)
            );

            fprintf(debug->fd,
                "%i %i %f %i %f\n",
                config->n + config->windowSize,
                3,
                creal(IQ),
                4,
                cimag(IQ)
            );

            break;
        }

        case NODEBUG:
        {
            break;
        }
    }
#endif
    return IQ;
}

iq_sample_t *getIQ(dft_t *config, int period, double windowPhaseReal, FILE *debug_fd)
{
    static iq_sample_t sample = {0};

    // all the IQ samples stored in these bins.
    //      0 - just before ideal midpoint
    //      1 - just after ideal midpoint
    //      2 - just before ideal
    //      3 - just after ideal
    static double complex IQsamples[4] = {0};

    dft_debug_t iq_midpoint_debug =
    {
        .flag = MIDPOINT,
        .fd = debug_fd,
    };

    dft_debug_t iq_aligned_debug =
    {
        .flag = ALLIGNED,
        .fd = debug_fd,
    };

    // compute DFT (rn just one freq, but when we implement OFDM, then at many harmonic(orthogonal) frequencies)
    // sample.iq = dft2(&dft_config, &iq_aligned_debug);

    // take a sample right between the ideal samples which are taken at windowPhase
    // if we are half full on the buffer, take an intermidiate IQ sample, for timing sync later
    if(config->windowPhase == ((int)floor(windowPhaseReal) + period / 2 - 1) % period)    // just before real window phase offset midpoint
    {
        IQsamples[0] = dft2(config, &iq_midpoint_debug);
    }

    if(config->windowPhase == ((int)ceil(windowPhaseReal) + period / 2 - 1) % period)     // just after real window phase offset midpoint;
    {
        IQsamples[1] = dft2(config, &iq_midpoint_debug);
    }

    if(config->windowPhase == ((int)floor(windowPhaseReal) + period - 1) % period)        // just before real window phase offset
    {
        IQsamples[2] = dft2(config, &iq_aligned_debug);
    }

    if(config->windowPhase == ((int)ceil(windowPhaseReal) + period - 1) % period)         // just after real window phase offset
    {
        IQsamples[3] = dft2(config, &iq_aligned_debug);
    }

    // just after real window phase offset, do the calculations that must be done once per symbol recieved
    if((config->windowPhase == ((int)ceil(windowPhaseReal) + SYMBOL_PERIOD - 1) % SYMBOL_PERIOD) && (config->n - sample.tookSampleAt > SYMBOL_PERIOD / 2))
    {
        // throwing away the last RMS value from MIDPOINT calculation, which is probably a shame and a waste

        // interpolate IQ and IQmidpoint from actual IQ samples in IQsamples
        // just doing linear interpolation
        //      (slope) * interp_X + initial
        //      (final - initial) * interp_X + initial
        sample.midpoint = (IQsamples[1] - IQsamples[0]) * fmod(windowPhaseReal, 1) + IQsamples[0];
        sample.last = sample.iq;
        sample.iq = (IQsamples[3] - IQsamples[2]) * fmod(windowPhaseReal, 1) + IQsamples[2];

        sample.tookSampleAt = config->n;
    }

    return &sample;
}

int main(void)
{
    int retval = 0;

    FILE* waveformPlotStdin = NULL;
    FILE* fftDebuggerStdin = NULL;
    FILE* errorPlotStdin = NULL;
    char *iq_plot_buffer = NULL;
    FILE* IQplotStdin = NULL;
    FILE* IQvstimeStdin = NULL;
    FILE* eyeDiagramRealStdin = NULL;
    FILE* eyeDiagramImaginaryStdin = NULL;

    // the OFDM channel number, how many cycles per symbol
    int k = 4;
    sample32_t sample;

    // buffer is the length of the symbol period, so that symbols are orthogonal
    double sampleBuffer[SYMBOL_PERIOD] = {0.0};

    // the offset of sample buffer or fftwindow in time, to get symbol time sync
    int windowPhase = 0;

    //double windowPhaseReal = (double)rand() / RAND_MAX * SYMBOL_PERIOD; // the floating point window phase, to be quantized into windowPhase
    double windowPhaseReal = 0;

    // int tookSampleAt = 0;

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

    // using it to plot the time domain signal
    waveformPlotStdin = popen(plot, "w");

    if (waveformPlotStdin == NULL)
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
    fftDebuggerStdin = popen(debugPlot, "w");

    if (fftDebuggerStdin == NULL)
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
    errorPlotStdin = popen(errorPlot, "w");

    if (errorPlotStdin == NULL)
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

    IQplotStdin = popen(iq_plot_buffer, "w");

    // Free buffer, even if popen failed
    free(iq_plot_buffer);
    iq_plot_buffer = NULL;

    if (IQplotStdin == NULL)
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
    eyeDiagramRealStdin = popen(eyeDiagramPlotReal, "w");

    char eyeDiagramPlotImaginary[300] = {0};
    sprintf(
        eyeDiagramPlotImaginary,
        "feedgnuplot "
        "--domain --dataid --lines --points --maxcurves %i "
        "--title \"Eye Diagram, Imaginary part\" "
        "--xlabel \"Time (Audio sample #)\" --ylabel \"Q\" ",
        3000);
    eyeDiagramImaginaryStdin = popen(eyeDiagramPlotImaginary, "w");

    // using it to plot the time domain signal
    IQvstimeStdin = popen(eyeDiagramPlotReal, "w");

    if (IQvstimeStdin == NULL)
    {
        fprintf(stderr, "Failed to create IQ vs Time plot: %s\n", strerror(errno));
        retval = 7;
        goto exit;
    }

    dft_t dft_config =
    {
        .samples = {0.0},
        .windowSize = SYMBOL_PERIOD,
        .offset = 0,
        .windowPhase = 0,
        .carrierPhase = 0.0,
        .k = k,
        .rms = 0.0,
    };

    dft_debug_t dft_debug_none =
    {
        .flag = NODEBUG,
        .fd = NULL,
    };

    avg_t avgPhaseOffsetEstimate = {0};
    avgInit(&avgPhaseOffsetEstimate, 2);

    avg_t avgRMSWindow = {0};
    avgInit(&avgRMSWindow, 5);

    iq_sample_t *iq_sample = NULL;

    // while there is data to recieve, not end of file
    for(int n = 0; n < SYMBOL_PERIOD * 2000; n++)
    {
        // recieve data on stdin, signed 32bit integer

        // Let's not use window phase to adjust index, and just pass the window phase to dft like a sane person
        dft_config.offset = (n % dft_config.windowSize) + 1;
        dft_config.windowPhase = n % dft_config.windowSize;
        dft_config.n = n;

        for(size_t i = 0; i < sizeof(sample.value); i++)
        {
            // get the bytes, little endian, of signed 32bit integer, using the struct to do type punning
            sample.byte[i] = getchar();
        }

        // convert to a double ranged from -1 to 1
        static double equalizationFactor = 1;

        // sampleBuffer[bufferIndex] = (double)sample.value / INT32_MAX * equalizationFactor;
        dft_config.samples[dft_config.windowPhase] = (double)sample.value / INT32_MAX * equalizationFactor;

        fprintf(waveformPlotStdin,
            "%i %f %f\n",
            n,
            sampleBuffer[dft_config.windowPhase],
            equalizationFactor
        );

    #if DEBUG_LEVEL > 0
        // oversampling IQ values for a nice eye diagram to help debug stuff
        double complex oversampledIQ = dft2(&dft_config, &dft_debug_none);

        fprintf(eyeDiagramRealStdin,
            "%f %i %f\n",
            fmod((double)n / SYMBOL_PERIOD, 4.0),
            n / 4 / SYMBOL_PERIOD,
            creal(oversampledIQ)
        );

        fprintf(eyeDiagramImaginaryStdin,
            "%f %i %f\n",
            fmod((double)n / SYMBOL_PERIOD, 4.0),
            n / 4 / SYMBOL_PERIOD,
            cimag(oversampledIQ)
        );
    #endif

        iq_sample = getIQ(&dft_config, SYMBOL_PERIOD, windowPhaseReal, fftDebuggerStdin);

        //if((bufferIndex == (windowPhase + SYMBOL_PERIOD - 1) % SYMBOL_PERIOD) && (n - tookSampleAt > SYMBOL_PERIOD / 2))
        // I added another condition to help debounce. sometimes it takes many samples in a row due to changing window offset
        if((dft_config.windowPhase == ((int)ceil(windowPhaseReal) + SYMBOL_PERIOD - 1) % SYMBOL_PERIOD) && (dft_config.n - iq_sample->tookSampleAt > SYMBOL_PERIOD / 2))  // just after real window phase offset, do the calculations that must be done once per symbol recieved
        {
            // now I'm doing a bunch of stuff that happens every IQ sample. This all happens in the timespan of a single audio sample, which is 1/symbolPeriod of the time between IQ samples that could be used, but whatever

            // averaging filter for the equalizer
            avgAppend(&avgRMSWindow, 1.0 / sqrt(2) - dft_config.rms);
            double rmsaverage = avgGet(&avgRMSWindow);
            //equalizationFactor += (1.0 / sqrt(2) - RMS) * 5.;

            // PID for the equalizer, just proportional. with max
            equalizationFactor = fmax(fmin(equalizationFactor + rmsaverage * 2, 1000), 0);

        #if DEBUG_LEVEL > 0
            //equalizationFactor = 1;
            equalizationFactor = 1.0 / 0.007;
        #endif

            // try to get a phase lock, symbol time lock, frequency lock, and equalization
            // calculate the error signal
            // Gardner Algorithm: Real Part( derivitive of IQ times conjugate of IQ)
            // basically, it's trying to estimate the error of zero crossings
            double phaseOffsetEstimate = creal((iq_sample->iq - iq_sample->last) * conj(iq_sample->midpoint));

        #if DEBUG_LEVEL > 1
            fprintf(fftDebuggerStdin,
                "%i %i %f %i %f\n",
                n,
                10,
                phaseOffsetEstimate + 2,
                15,
                (double)(n % SYMBOL_PERIOD) / SYMBOL_PERIOD + 2
            );
        #endif

            static timing_lock_t lockstate = NO_LOCK;

            // Process Variable (PV, ie phase estimate) filter
            // rolling average of phase offset estimate
            // this may need to be adjusted based on the state of symbol and phase lock achieved
            switch(lockstate)
            {
                case NO_LOCK:
                {
                    // smaller averaging window to achieve
                    //  faster error signal response
                    //  for more aggressive PID tune
                    if (avgPhaseOffsetEstimate.numSamples != 2)
                    {
                        avgReInitCopy(&avgPhaseOffsetEstimate, 2);
                    }

                    break;
                }

                case SYMBOL_LOCK:
                    // Fallthrough
                case PHASE_LOCK:
                {
                    // longer average window to help average out symbol transitions that are not zero crossings
                    //  PID tune must be shittier
                    if (avgPhaseOffsetEstimate.numSamples != 64)
                    {
                        avgReInitCopy(&avgPhaseOffsetEstimate, 64);
                    }

                    break;
                }
            }

            // add a sample and take the average of last averageSize samples
            avgAppend(&avgPhaseOffsetEstimate, phaseOffsetEstimate);
            double phaseOffsetEstimateAverage = avgGet(&avgPhaseOffsetEstimate);

            // PID loop for symbol timing, ie aligning the fft window to the symbol transitions
            static double P_gain;
            static double I_gain;
            static double D_gain;

            switch(lockstate)
            {
                case NO_LOCK:
                {
                    //P_gain = 2;
                    P_gain = 0.1;
                    //I_gain = 0.02;
                    I_gain = 0;
                    D_gain = 0;
                    break;
                }

                case SYMBOL_LOCK:
                    // Fallthrough
                case PHASE_LOCK:
                {
                    P_gain = 0.1;
                    I_gain = 0.001;
                    D_gain = 0;
                    break;
                }
            }

            double error = 0 - phaseOffsetEstimateAverage;
            //error *= -1;
            static double errorIntegral = 0;
            errorIntegral += error;

            static double lastError = 0;
            double errorDerivative = error - lastError;
            lastError = error;

            // this may need to be adjusted based on the state of symbol and phase lock achieved
            double phaseAdjustment = errorDerivative * D_gain + errorIntegral * I_gain + error * P_gain;

            windowPhaseReal += phaseAdjustment;
            windowPhaseReal = fmod(windowPhaseReal + phaseAdjustment, SYMBOL_PERIOD);
            windowPhaseReal = windowPhaseReal < 0 ? SYMBOL_PERIOD + windowPhaseReal : windowPhaseReal;
            windowPhase = (int)round(windowPhaseReal) % SYMBOL_PERIOD; // quantize the real window phase
            //windowPhase = windowPhase < 0 ? SYMBOL_PERIOD + windowPhase : windowPhase;

        #if DEBUG_LEVEL > 0
            // some options to overwrite the window phase given by the PID controller
            //windowPhase = (n * 2 / 2000) % SYMBOL_PERIOD;
            //windowPhase = (n * 4 * 2/ (SYMBOL_PERIOD * 2000) ) % 4 * SYMBOL_PERIOD / 4;
            //windowPhase = SYMBOL_PERIOD / 4;
            //windowPhase = 0;
            //windowPhaseReal = 0.1;
            //windowPhase = (int)((windowPhase + phaseAdjustment) < 0 ? (SYMBOL_PERIOD - windowPhase + phaseAdjustment) : (windowPhase + phaseAdjustment)) % SYMBOL_PERIOD;
        #endif

            // extract the frequencies to be decoded
            // add the relevant IQ values to an array
            // plot the IQ and a tail with some kind of persistance to give an animated output
            //fprintf(IQplotStdin, "%f %f\n", I, Q);
            //printf("%f %f\n", I, Q);

            // plot with a new color for each window phase
            fprintf(IQplotStdin,
                "%f %i %f\n",
                creal(iq_sample->iq),
                windowPhase,
                cimag(iq_sample->iq)
            );

            //fprintf(IQplotStdin,
            //     "%f %i %f\n",
            //     creal(iq_sample->midpoint),
            //     windowPhase + SYMBOL_PERIOD,
            //     cimag(iq_sample->midpoint)
            // );

            fprintf(errorPlotStdin,
                "%i, %f %f %f %f %f\n",
                n / SYMBOL_PERIOD,
                -phaseOffsetEstimate,
                error,
                phaseAdjustment,
                (double)windowPhase / SYMBOL_PERIOD,
                windowPhaseReal / SYMBOL_PERIOD
            );

            //fprintf(IQvstimeStdin,
            //     "%i, %f, %f, %f, %f, %f, %f\n",
            //     n / SYMBOL_PERIOD % (2*3),
            //     creal(iq_sample->iq),
            //     cimag(iq_sample->iq),
            //     creal(iq_sample->last),
            //     cimag(iq_sample->last),
            //     creal(iq_sample->midpoint),
            //     cimag(iq_sample->midpoint)
            // );

            fprintf(IQvstimeStdin,
                "%f, %f, %f\n",
                (n / SYMBOL_PERIOD) % (2 * 3) + 0.0,
                creal(iq_sample->iq),
                cimag(iq_sample->iq)
            );

            fprintf(IQvstimeStdin,
                "%f, %f, %f\n",
                (n / SYMBOL_PERIOD) % (2 * 3) + 0.5,
                creal(iq_sample->midpoint),
                cimag(iq_sample->midpoint)
            );
        }

         // I think I'll use feedgnuplot, so pipe the IQ values into feedgnuplot
    }

exit:
    if (iq_plot_buffer != NULL)
    {
        free(iq_plot_buffer);
        iq_plot_buffer = NULL;
    }

    // close the streams, allowing the plots to render and open a window, then wait for them to terminate in separate threads.
    if((IQvstimeStdin != NULL) && (fork() == 0))
    {
        pclose(IQvstimeStdin);
        return 0;
    }

    if((IQplotStdin != NULL) && (fork() == 0))
    {
        pclose(IQplotStdin);
        return 0;
    }

    if((errorPlotStdin != NULL) && (fork() == 0))
    {
        pclose(errorPlotStdin);
        return 0;
    }

    if((waveformPlotStdin != NULL) && (fork() == 0))
    {
        pclose(waveformPlotStdin);
        return 0;
    }

    if((fftDebuggerStdin != NULL) && (fork() == 0))
    {
        pclose(fftDebuggerStdin);
        return 0;
    }
    if((eyeDiagramRealStdin != NULL) && (fork() == 0))
    {
        pclose(eyeDiagramRealStdin);
        return 0;
    }
    if((eyeDiagramImaginaryStdin != NULL) && (fork() == 0))
    {
        pclose(eyeDiagramImaginaryStdin);
        return 0;
    }

    // after all the data is recieved, generate a plot of IQ
    return retval;
}
