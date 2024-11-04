#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>
#include <stdint.h>
#include <math.h>
#include <complex.h>

typedef struct
{
    union
    {
        int32_t value;
        uint8_t byte[sizeof(int32_t)];
    };
} sample_32;

// calculate the discrete fourier transform of an array of possibly complex values but only at frequency 'k' (k cycles per windowSize samples)
double complex dft(double* buffer, int windowSize, int offset, int k)
{
    double complex IQ = 0;

    int bufferIndex;
    double phase;
    for(int i = 0; i < windowSize; i++)
    {
        // starts at buffer[offset] and wraps around to the beginning of the buffer
        bufferIndex = (i + offset + 1) % windowSize;
        // phase of the complex exponential
        phase = (double)i * k / windowSize;
        IQ += buffer[bufferIndex] * cexp(I*M_2_PI*phase);
    }
    // normalization factor (do I need to divide by k?)
    IQ *= sqrt(1. / windowSize);

    return IQ;
}

int main(void)
{
    // length of each symbol in samples, and so the fft window.
    // This is the period of time all the orthogonal symbols will be integrated over
    #define symbolPeriod 64

    // the OFDM channel number, how many cycles per symbol
    int k = 4;
    sample_32 sample;

    // buffer is the length of the symbol period, so that symbols are orthogonal
    double sampleBuffer[symbolPeriod] = {0.0};

    for(int i = 0; i < symbolPeriod; i++)
    {
        sampleBuffer[i] = 0;
    }

    // the offset of sample buffer or fftwindow in time, to get symbol time sync
    int windowPhase = 0;

    //double windowPhaseReal = (double)rand() / RAND_MAX * symbolPeriod; // the floating point window phase, to be quantized into windowPhase
    double windowPhaseReal = 0;

    int tookSampleAt = 0;

    // process arguments
    //  none right now

    char buffer[10000] = {0};
    int stringLength = 0;

    // for live plotting, pipe to feedgnuplot
    stringLength += sprintf(buffer + stringLength,"feedgnuplot ");
    stringLength += sprintf(buffer + stringLength,"--domain --lines --points ");
    stringLength += sprintf(buffer + stringLength,"--title \"Raw signal Time domain\" ");
    stringLength += sprintf(buffer + stringLength,"--xlabel \"Time (microphone sample #\" --ylabel \"value\" ");
    stringLength += sprintf(buffer + stringLength,"--legend 0 \"Signal\" ");
    stringLength += sprintf(buffer + stringLength,"--legend 1 \"equalization factor\" ");

    // using it to plot the time domain signal
    FILE* waveformPlotStdin = popen(buffer, "w");


    stringLength = 0;
    stringLength += sprintf(buffer + stringLength, "feedgnuplot ");
    stringLength += sprintf(buffer + stringLength, "--domain --dataid --lines --points ");
    stringLength += sprintf(buffer + stringLength, "--title \"FFT debuger graph\" ");
    stringLength += sprintf(buffer + stringLength, "--xlabel \"Time (microphone sample #\" --ylabel \"value\" ");
    stringLength += sprintf(buffer + stringLength, "--legend 0 \"input samples\" ");
    stringLength += sprintf(buffer + stringLength, "--legend 1 \"fft real\" ");
    stringLength += sprintf(buffer + stringLength, "--legend 2 \"fft imaginary\" ");
    stringLength += sprintf(buffer + stringLength, "--legend 3 \"I decision\" ");
    stringLength += sprintf(buffer + stringLength, "--legend 4 \"Q decision\" ");
    stringLength += sprintf(buffer + stringLength, "--legend 5 \"input samples mid\" ");
    stringLength += sprintf(buffer + stringLength, "--legend 6 \"fft real mid\" ");
    stringLength += sprintf(buffer + stringLength, "--legend 7 \"fft imaginary mid\" ");
    stringLength += sprintf(buffer + stringLength, "--legend 8 \"I decision mid\" ");
    stringLength += sprintf(buffer + stringLength, "--legend 9 \"Q decision mid\" ");
    stringLength += sprintf(buffer + stringLength, "--legend 10 \"phase error signal\" ");
    stringLength += sprintf(buffer + stringLength, "--legend 11 \"samp*real\" ");
    stringLength += sprintf(buffer + stringLength, "--legend 12 \"samp*imag\" ");
    stringLength += sprintf(buffer + stringLength, "--legend 13 \"integral real\" ");
    stringLength += sprintf(buffer + stringLength, "--legend 14 \"integral imag\" ");

    // using it to plot the time domain signal
    FILE* fftDebuggerStdin = popen(buffer, "w");


    stringLength = 0;
    stringLength += sprintf(buffer + stringLength, "feedgnuplot ");
    stringLength += sprintf(buffer + stringLength, "--domain --lines --points ");
    stringLength += sprintf(buffer + stringLength, "--title \"Time Domain error signal\" ");
    stringLength += sprintf(buffer + stringLength, "--legend 0 \"estimated phase offset\" ");
    stringLength += sprintf(buffer + stringLength, "--legend 1 \"rolling average error signal\" ");
    stringLength += sprintf(buffer + stringLength, "--legend 2 \"PI filtered signal\" ");
    stringLength += sprintf(buffer + stringLength, "--legend 3 \"fft window phase shift\" ");
    stringLength += sprintf(buffer + stringLength, "--legend 4 \"real target phase offset\" ");

    // using it to plot the time domain signal
    FILE* errorPlotStdin = popen(buffer, "w");
    //FILE* IQplotStdin = popen("feedgnuplot --domain --points --title \"IQ plot\" --unset key", "w");


    stringLength = 0;
    stringLength += sprintf(buffer + stringLength, "feedgnuplot ");
    stringLength += sprintf(buffer + stringLength, "--dataid --domain --points --maxcurves %i ", symbolPeriod * 2 + 1);
    stringLength += sprintf(buffer + stringLength, "--title \"IQ plot\" ");
    stringLength += sprintf(buffer + stringLength, "--xlabel \"I\" --ylabel \"Q\" ");

    for(int i = 0; i < symbolPeriod; i++)
    {
        stringLength += sprintf(buffer + stringLength,"--legend %i \"Phase offset %i samples\" ", i, i);
    }

    FILE* IQplotStdin = popen(buffer, "w");

    // for ploting IQ values over time to hopefully obtain an error function
    stringLength = 0;
    stringLength += sprintf(buffer + stringLength, "feedgnuplot ");
    stringLength += sprintf(buffer + stringLength, "--domain --lines --points --unset grid ");
    stringLength += sprintf(buffer + stringLength, "--title \"IQ vs time Eye Diagram\" ");
    stringLength += sprintf(buffer + stringLength, "--xlabel \"Time (IQsample #)\" --ylabel \"value\" ");
    stringLength += sprintf(buffer + stringLength, "--legend 0 \"I\" --legend 1 \"Q\" ");
    stringLength += sprintf(buffer + stringLength, "--legend 2 \"IQlastReal\" --legend 3 \"IQlastImag\" ");
    stringLength += sprintf(buffer + stringLength, "--legend 4 \"IQmidReal\" --legend 5 \"IQmidImag\" ");
    FILE* IQvstimeStdin = popen(buffer, "w");    // using it to plot the time domain signal

    // while there is data to recieve, not end of file
    for(int n = 0; n < symbolPeriod * 2000; n++)
    {
        // recieve data on stdin, signed 32bit integer

        // use the windowphase to adjust the buffer index position
        int bufferIndex = (n + windowPhase)%symbolPeriod;

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

        if((bufferIndex == symbolPeriod / 2 - 1) && (n - tookSampleAt > 1))
        {
            IQmidpoint = 0;

            // compute DFT mid symbol (if time is already synced, otherwise it will help sync the time)
            for(int i = 0; i < symbolPeriod; i++)
            {
                // phasor offsets the cos and sin waves so that they allign with the time sequence of data in the half overwritten buffer

                // index that is based on the fact we are starting partway through the buffer
                int index = (i + bufferIndex + 1)%symbolPeriod;

                double phase = (double)(i) * k / symbolPeriod;

                //IQmidpoint += sampleBuffer[i] * cos(2*M_PI*phase) + I * sampleBuffer[i] * sin(2*M_PI*phase);

                double complex wave = cexp(I * 2 * M_PI * phase);
                double complex value = sampleBuffer[index] * wave;

                IQmidpoint += value;

                // debug graph outputs
                fprintf(fftDebuggerStdin, "%i %i %f %i %f %i %f\n", n + i, 5, sampleBuffer[index] + 4, 6, creal(wave) + 4, 7, cimag(wave) + 4);
                fprintf(fftDebuggerStdin, "%i %i %f %i %f %i %f %i %f\n", n + i, 11, creal(value) + 6, 12, cimag(value) + 6, 13, creal(IQmidpoint) + 6, 14, cimag(IQmidpoint) + 6);
            }

            // normalization factor
            IQmidpoint *= sqrt(1. / symbolPeriod) / k;

            // debug fft plot
            fprintf(fftDebuggerStdin, "%i %i %f %i %f\n", n, 8, creal(IQmidpoint) + 4, 9, cimag(IQmidpoint) + 4);
            fprintf(fftDebuggerStdin, "%f %i %f %i %f\n", n + symbolPeriod - 0.01, 8, creal(IQmidpoint) + 4, 9, cimag(IQmidpoint) + 4);
            //printf("midpoint IQ: %f + %fi\n", creal(IQmidpoint), cimag(IQmidpoint));
        }

        // if the window buffer is filled, ie, we're on the last index of the buffer
        // I added another condition to help debounce. sometimes it takes many samples in a row due to changing window offset
        static double complex IQ = 0;
        static double complex IQlast = 0;
        double RMS = 0;

        if((bufferIndex == symbolPeriod - 1) && (n - tookSampleAt > symbolPeriod / 2))
        {
            tookSampleAt = n;
            IQlast = IQ;
            IQ = 0;

            // compute DFT (rn just one freq, but when we implement OFDM, then at many harmonic(orthogonal) frequencies)
            // here is a simple quadrature detector
            for(int i = 0; i < symbolPeriod; i++)
            {
                // use a summation over the symbol period to separate orthogonal components (quadrature components) at the single frequency
                double phase = (double)i * k /symbolPeriod;  // one cycle per symbol period

                //double phase = (double)((i + bufferIndex + 1)%symbolPeriod) * k / symbolPeriod;
                //IQ += sampleBuffer[i] * cos(2*M_PI*phase) + I * sampleBuffer[i] * sin(2*M_PI*phase);

                double complex wave = cexp(I * 2 * M_PI * phase);
                IQ += sampleBuffer[i] * wave;

                // compute RMS amplitude for equalization
                RMS += pow(sampleBuffer[i], 2);

                // debug graph outputs
                fprintf(fftDebuggerStdin, "%i %i %f %i %f %i %f\n", n + i, 0, sampleBuffer[i], 1, creal(wave), 2, cimag(wave));
            }

            // normalization with a unitary normalization factor
            IQ *= sqrt(1. / symbolPeriod) / k;

            // complete the RMS fomula
            RMS = sqrt(1./symbolPeriod * RMS);

            // debug fft plot
            fprintf(fftDebuggerStdin, "%i %i %f %i %f\n", n, 3, creal(IQ), 4, cimag(IQ));
            fprintf(fftDebuggerStdin, "%f %i %f %i %f\n", n + symbolPeriod - 0.01, 3, creal(IQ), 4, cimag(IQ));

            // averaging filter for the equalizer
            static double rmsaverageWindow[5] = {0};
            int rmsaverageSize = 5;
            int rmsaverageIndex = (n / symbolPeriod) % rmsaverageSize;

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
            int averageIndex = (n / symbolPeriod) % averageSize;

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
            //windowPhaseReal = fmod(windowPhaseReal + phaseAdjustment, symbolPeriod);
            //windowPhaseReal = windowPhaseReal < 0 ? symbolPeriod + windowPhaseReal : windowPhaseReal;
            //windowPhase = (int)round(windowPhaseReal) % symbolPeriod;
            windowPhase = (int)round(windowPhaseReal) % symbolPeriod; // quantize the real window phase
            windowPhase = windowPhase < 0 ? symbolPeriod + windowPhase : windowPhase;
            //windowPhase = (n * 2 / 2000) % symbolPeriod;
            //windowPhase = (n * 4 * 2/ (symbolPeriod * 2000) ) % 4 * symbolPeriod / 4;
            //windowPhase = symbolPeriod / 4;
            windowPhase = 0;
            //windowPhase = (int)((windowPhase + phaseAdjustment) < 0 ? (symbolPeriod - windowPhase + phaseAdjustment) : (windowPhase + phaseAdjustment)) % symbolPeriod;

            // extract the frequencies to be decoded
            // add the relevant IQ values to an array
            // plot the IQ and a tail with some kind of persistance to give an animated output
            //fprintf(IQplotStdin, "%f %f\n", I, Q);
            //printf("%f %f\n", I, Q);

            // plot with a new color for each window phase
            fprintf(IQplotStdin, "%f %i %f\n", creal(IQ), windowPhase, cimag(IQ));
            //fprintf(IQplotStdin, "%f %i %f\n", creal(IQmidpoint), windowPhase + symbolPeriod, cimag(IQmidpoint));
            fprintf(errorPlotStdin, "%i, %f %f %f %f %f\n", n / symbolPeriod, -phaseOffsetEstimate, error, phaseAdjustment, (double)windowPhase / symbolPeriod, windowPhaseReal / symbolPeriod);
            //fprintf(IQvstimeStdin, "%i, %f, %f, %f, %f, %f, %f\n", n / symbolPeriod % (2*3), creal(IQ), cimag(IQ), creal(IQlast), cimag(IQlast), creal(IQmidpoint), cimag(IQmidpoint));
            fprintf(IQvstimeStdin, "%f, %f, %f\n", (n / symbolPeriod) % (2*3) + 0., creal(IQ), cimag(IQ));
            fprintf(IQvstimeStdin, "%f, %f, %f\n", (n / symbolPeriod) % (2*3) + 0.5, creal(IQmidpoint), cimag(IQmidpoint));
        }

         // I think I'll use feedgnuplot, so pipe the IQ values into feedgnuplot
    }

    // close the streams, allowing the plots to render and open a window, then wait for them to terminate in separate threads.
    if(!fork())
    {
        pclose(IQvstimeStdin);
        return 0;
    }

    if(!fork())
    {
        pclose(IQplotStdin);
        return 0;
    }

    if(!fork())
    {
        pclose(errorPlotStdin);
        return 0;
    }

    if(!fork())
    {
        pclose(waveformPlotStdin);
        return 0;
    }

    if(!fork())
    {
        pclose(fftDebuggerStdin);
        return 0;
    }

    // after all the data is recieved, generate a plot of IQ
    return 0;
}
