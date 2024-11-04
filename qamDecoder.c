#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>
#include <stdint.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <errno.h>

// length of each symbol in samples, and so the fft window.
// This is the period of time all the orthogonal symbols will be integrated over
#define SYMBOL_PERIOD 64

// Debug flag
#define DEBUG_LEVEL 1

typedef struct
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
} dft_debug_t;

// calculate the discrete fourier transform of an array of real values but only at frequency 'k' (k cycles per windowSize samples)
//  debugFlag is to print the right debug info for different situations
double complex dft(double* buffer, int windowSize, int offset, int k, double* rmsOut, dft_debug_t debugFlag, FILE* debug_fd, int debug_n)
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

    for(int i = 0; i < windowSize; i++)
    {
        // starts at buffer[offset] and wraps around to the beginning of the buffer
        bufferIndex = (i + offset + 1) % windowSize;
        // phase of the complex exponential
        // phasor offsets the cos and sin waves so that they allign with the time sequence of data in the half overwritten buffer
        phase = (double)i * k / windowSize;

        // use a summation over the symbol period to separate orthogonal components (quadrature components) at the single frequency
        double complex wave = cexp(I*M_PI*phase); // generate a sample from the complex exponential
        double complex value = buffer[bufferIndex] * wave;  // multiply the complex exponential by the input sample
        IQ += value;    // integrate the result over the window

        // compute RMS amplitude for equalization
        RMS += pow(buffer[i], 2);


        // this debug define simplifies the function a bit if debugging is disabled
    #if DEBUG_LEVEL == 1
        switch(debugFlag)
        {
            case MIDPOINT:
            {
                // debug graph outputs
                fprintf(debug_fd, "%i %i %f %i %f %i %f\n", debug_n + i, 5, buffer[bufferIndex] + 4, 6, creal(wave) + 4, 7, cimag(wave) + 4);
                fprintf(debug_fd, "%i %i %f %i %f %i %f %i %f\n", debug_n + i, 11, creal(value) + 6, 12, cimag(value) + 6, 13, creal(IQ) + 6, 14, cimag(IQ) + 6);
                break;
            }

            case ALLIGNED:
            {
                // debug graph outputs
                fprintf(debug_fd, "%i %i %f %i %f %i %f\n", debug_n + i, 0, buffer[i], 1, creal(wave), 2, cimag(wave));
                break;
            }
        }
    #endif
    }
    // normalization factor (do I need to divide by k?)
    IQ *= sqrt(1. / windowSize);

    // complete the RMS fomula
    RMS = sqrt(1./SYMBOL_PERIOD * RMS);

#if DEBUG_LEVEL == 1
    switch(debugFlag)
    {
        case MIDPOINT:
        {
            // debug fft plot
            fprintf(debug_fd, "%i %i %f %i %f\n", debug_n, 8, creal(IQ) + 4, 9, cimag(IQ) + 4);
            fprintf(debug_fd, "%f %i %f %i %f\n", debug_n + windowSize - 0.01, 8, creal(IQ) + 4, 9, cimag(IQ) + 4);
            break;
        }

        case ALLIGNED:
        {
            // debug fft plot
            fprintf(debug_fd, "%i %i %f %i %f\n", debug_n, 3, creal(IQ), 4, cimag(IQ));
            fprintf(debug_fd, "%f %i %f %i %f\n", debug_n + windowSize - 0.01, 3, creal(IQ), 4, cimag(IQ));
            break;
        }
    }
#endif

    *rmsOut = RMS;
    return IQ;
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

    // the OFDM channel number, how many cycles per symbol
    int k = 4;
    sample32_t sample;

    // buffer is the length of the symbol period, so that symbols are orthogonal
    double sampleBuffer[SYMBOL_PERIOD] = {0.0};

    // the offset of sample buffer or fftwindow in time, to get symbol time sync
    int windowPhase = 0;

    //double windowPhaseReal = (double)rand() / RAND_MAX * SYMBOL_PERIOD; // the floating point window phase, to be quantized into windowPhase
    double windowPhaseReal = 0;

    int tookSampleAt = 0;

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
        "--legend 10 \"phase error signal\" "
        "--legend 11 \"samp*real\" "
        "--legend 12 \"samp*imag\" "
        "--legend 13 \"integral real\" "
        "--legend 14 \"integral imag\" "
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
    const char *iqPlot =
        "feedgnuplot "
        "--domain --lines --points --unset grid "
        "--title \"IQ vs time Eye Diagram\" "
        "--xlabel \"Time (IQsample #)\" --ylabel \"value\" "
        "--legend 0 \"I\" --legend 1 \"Q\" "
        "--legend 2 \"IQlastReal\" --legend 3 \"IQlastImag\" "
        "--legend 4 \"IQmidReal\" --legend 5 \"IQmidImag\" "
    ;

    // using it to plot the time domain signal
    IQvstimeStdin = popen(iqPlot, "w");

    if (IQvstimeStdin == NULL)
    {
        fprintf(stderr, "Failed to create IQ vs Time plot: %s\n", strerror(errno));
        retval = 7;
        goto exit;
    }

    // while there is data to recieve, not end of file
    for(int n = 0; n < SYMBOL_PERIOD * 2000; n++)
    {
        // recieve data on stdin, signed 32bit integer

        // use the windowphase to adjust the buffer index position
        int bufferIndex = (n + windowPhase)%SYMBOL_PERIOD;

        for(size_t i = 0; i < sizeof(sample.value); i++)
        {
            // get the bytes, little endian, of signed 32bit integer, using the struct to do type punning
            sample.byte[i] = getchar();
        }

        // convert to a double ranged from -1 to 1
        static double equalizationFactor = 1;
        sampleBuffer[bufferIndex] = (double)sample.value / INT32_MAX * equalizationFactor;

        fprintf(waveformPlotStdin, "%i %f %f\n", n, sampleBuffer[bufferIndex], equalizationFactor);

        // if we are half full on the buffer, take an intermidiate IQ sample, for timing sync later
        static double complex IQmidpoint = 0;
        double RMS;

        if((bufferIndex == SYMBOL_PERIOD / 2 - 1) && (n - tookSampleAt > 1))
        {
            IQmidpoint = 0;

            // compute DFT half way between symbols (if time is already synced, otherwise it will help sync the time)
            IQmidpoint = dft(sampleBuffer, SYMBOL_PERIOD, bufferIndex, k, &RMS, MIDPOINT, fftDebuggerStdin, n);

        }

        // if the window buffer is filled, ie, we're on the last index of the buffer
        // I added another condition to help debounce. sometimes it takes many samples in a row due to changing window offset
        static double complex IQ = 0;
        static double complex IQlast = 0;

        if((bufferIndex == SYMBOL_PERIOD - 1) && (n - tookSampleAt > SYMBOL_PERIOD / 2))
        {
            tookSampleAt = n;
            IQlast = IQ;
            IQ = 0;

            // compute DFT (rn just one freq, but when we implement OFDM, then at many harmonic(orthogonal) frequencies)
            IQ = dft(sampleBuffer, SYMBOL_PERIOD, bufferIndex, k, &RMS, ALLIGNED, fftDebuggerStdin, n);
            // throwing away the last RMS value from MIDPOINT calculation, which is probably a shame and a waste
            

            // now I'm doing a bunch of stuff that happens every IQ sample. This all happens in the timespan of a single audio sample, which is 1/symbolPeriod of the time between IQ samples that could be used, but whatever

            // averaging filter for the equalizer
            static double rmsaverageWindow[5] = {0};
            int rmsaverageSize = 5;
            int rmsaverageIndex = (n / SYMBOL_PERIOD) % rmsaverageSize;

            rmsaverageWindow[rmsaverageIndex] = 1./sqrt(2) - RMS;
            double rmsaverage = 0;

            for(int i = 0; i < rmsaverageSize; i++)
            {
                rmsaverage += rmsaverageWindow[i];
            }

            rmsaverage /= rmsaverageSize;
            //equalizationFactor += (1. / sqrt(2) - RMS) * 5.;

            // PID for the equalizer, just proportional. with max
            equalizationFactor = fmax(fmin(equalizationFactor + rmsaverage * 2, 1000), 0);

            // try to get a phase lock, symbol time lock, frequency lock, and equalization
            // calculate the error signal
            // Gardner Algorithm: Real Part( derivitive of IQ times conjugate of IQ)
            // basically, it's trying to estimate the error of zero crossings
            double phaseOffsetEstimate = -creal((IQ - IQlast) * conj(IQmidpoint));

            fprintf(fftDebuggerStdin, "%i %i %f\n", n, 10, phaseOffsetEstimate + 2);

            // Process Variable (PV, ie phase estimate) filter
            // rolling average of phase offset estimate
            static double averageWindow[40] = {0};
            int averageSize = 40;
            int averageIndex = (n / SYMBOL_PERIOD) % averageSize;

            averageWindow[averageIndex] = phaseOffsetEstimate;

            double average = 0;

            for(int i = 0; i < averageSize; i++)
            {
                average += averageWindow[i];
            }

            average /= averageSize;

            // PID loop
            double error = 0 - average;
            static double errorIntegral = 0;
            errorIntegral += error;

            static double lastError = 0;
            double errorDerivative = error - lastError;
            lastError = error;

            //double phaseAdjustment = (errorDerivative * -1.0 + errorIntegral * 0.0 + error * 6.0) * 1;
            //double phaseAdjustment = (errorDerivative * 0.0 + errorIntegral * 0.00 + error * 0.1) * 1;
            double phaseAdjustment = (errorDerivative * 0.1 + errorIntegral * 0.010 + error * 0.5) * 1;

            windowPhaseReal += phaseAdjustment;
            //windowPhaseReal = fmod(windowPhaseReal + phaseAdjustment, SYMBOL_PERIOD);
            //windowPhaseReal = windowPhaseReal < 0 ? SYMBOL_PERIOD + windowPhaseReal : windowPhaseReal;
            //windowPhase = (int)round(windowPhaseReal) % SYMBOL_PERIOD;
            windowPhase = (int)round(windowPhaseReal) % SYMBOL_PERIOD; // quantize the real window phase
            windowPhase = windowPhase < 0 ? SYMBOL_PERIOD + windowPhase : windowPhase;
            //windowPhase = (n * 2 / 2000) % SYMBOL_PERIOD;
            //windowPhase = (n * 4 * 2/ (SYMBOL_PERIOD * 2000) ) % 4 * SYMBOL_PERIOD / 4;
            //windowPhase = SYMBOL_PERIOD / 4;
            windowPhase = 0;
            //windowPhase = (int)((windowPhase + phaseAdjustment) < 0 ? (SYMBOL_PERIOD - windowPhase + phaseAdjustment) : (windowPhase + phaseAdjustment)) % SYMBOL_PERIOD;

            // extract the frequencies to be decoded
            // add the relevant IQ values to an array
            // plot the IQ and a tail with some kind of persistance to give an animated output
            //fprintf(IQplotStdin, "%f %f\n", I, Q);
            //printf("%f %f\n", I, Q);

            // plot with a new color for each window phase
            fprintf(IQplotStdin, "%f %i %f\n", creal(IQ), windowPhase, cimag(IQ));
            //fprintf(IQplotStdin, "%f %i %f\n", creal(IQmidpoint), windowPhase + SYMBOL_PERIOD, cimag(IQmidpoint));
            fprintf(errorPlotStdin, "%i, %f %f %f %f %f\n", n / SYMBOL_PERIOD, -phaseOffsetEstimate, error, phaseAdjustment, (double)windowPhase / SYMBOL_PERIOD, windowPhaseReal / SYMBOL_PERIOD);
            //fprintf(IQvstimeStdin, "%i, %f, %f, %f, %f, %f, %f\n", n / SYMBOL_PERIOD % (2*3), creal(IQ), cimag(IQ), creal(IQlast), cimag(IQlast), creal(IQmidpoint), cimag(IQmidpoint));
            fprintf(IQvstimeStdin, "%f, %f, %f\n", (n / SYMBOL_PERIOD) % (2*3) + 0., creal(IQ), cimag(IQ));
            fprintf(IQvstimeStdin, "%f, %f, %f\n", (n / SYMBOL_PERIOD) % (2*3) + 0.5, creal(IQmidpoint), cimag(IQmidpoint));
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

    // after all the data is recieved, generate a plot of IQ
    return retval;
}
