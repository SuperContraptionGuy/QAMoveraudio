#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>
#include <stdint.h>
#include <math.h>
#include <complex.h>



int main(int argc, char** args)
{
    typedef struct {
        union {
            int32_t value;
            uint8_t byte[sizeof(int32_t)];
        };
    } sample_32;
    const int symbolPeriod = 32;       // length of each symbol in samples, and so the fft window. This is the period of time all the orthogonal symbols will be integrated over
    sample_32 sample;
    double sampleBuffer[symbolPeriod];  // buffer is the length of the symbol period, so that symbols are orthogonal
    for(int i = 0; i < symbolPeriod; i++)
    {
        sampleBuffer[i] = 0;
    }
    int windowPhase = 0;    // the offset of sample buffer or fftwindow in time, to get symbol time sync
    double windowPhaseReal = (double)rand() / RAND_MAX * symbolPeriod; // the floating point window phase, to be quantized into windowPhase
    //double windowPhaseReal = 0;

    int tookSampleAt = 0;


    // process arguments
    //  none right now

    char buffer[10000];
    int stringLength = 0;
    // for live plotting, pipe to feedgnuplot
    stringLength = 0;
    stringLength += sprintf(buffer + stringLength,"feedgnuplot ");
    stringLength += sprintf(buffer + stringLength,"--domain --lines --points ");
    stringLength += sprintf(buffer + stringLength,"--title \"Raw signal Time domain\" ");
    stringLength += sprintf(buffer + stringLength,"--xlabel \"Time (microphone sample #\" --ylabel \"value\" ");
    FILE* waveformPlotStdin = popen("feedgnuplot --domain --lines --points --title \"Time Domain\" --unset key --unset grid", "w");    // using it to plot the time domain signal
    stringLength = 0;
    stringLength += sprintf(buffer + stringLength,"feedgnuplot ");
    stringLength += sprintf(buffer + stringLength,"--domain --lines --points ");
    stringLength += sprintf(buffer + stringLength,"--title \"Time Domain error signal\" ");
    stringLength += sprintf(buffer + stringLength,"--legend 0 \"estimated phase offset\" ");
    stringLength += sprintf(buffer + stringLength,"--legend 1 \"rolling average error signal\" ");
    stringLength += sprintf(buffer + stringLength,"--legend 2 \"PI filtered signal\" ");
    stringLength += sprintf(buffer + stringLength,"--legend 3 \"fft window phase shift\" ");
    stringLength += sprintf(buffer + stringLength,"--legend 4 \"real target phase offset\" ");
    FILE* errorPlotStdin = popen(buffer, "w");    // using it to plot the time domain signal
    //FILE* IQplotStdin = popen("feedgnuplot --domain --points --title \"IQ plot\" --unset key", "w");
    stringLength = 0;
    stringLength += sprintf(buffer + stringLength,"feedgnuplot ");
    stringLength += sprintf(buffer + stringLength,"--dataid --domain --points ");
    stringLength += sprintf(buffer + stringLength,"--title \"IQ plot\" ");
    stringLength += sprintf(buffer + stringLength,"--xlabel \"I\" --ylabel \"Q\" ");
    for(int i = 0; i < symbolPeriod; i++)
        stringLength += sprintf(buffer + stringLength,"--legend %i \"Phase offset %i samples\" ", i, i);
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
    for(int n = 0;n < symbolPeriod * 2000;n++)
    {
        // recieve data on stdin, two signed 32bit integers

        int bufferIndex = (n + windowPhase)%symbolPeriod;   // use the windowphase to adjust the buffer index position
        for(int i = 0; i < sizeof(sample.value); i++)
        {
            // get the bytes, little endian, of signed 32bit integer, using the struct to do type punning
            sample.byte[i] = getchar();
        }
        // convert to a double ranged from -1 to 1
        sampleBuffer[bufferIndex] = (double)sample.value / INT32_MAX;
        fprintf(waveformPlotStdin, "%i %f\n", n, sampleBuffer[bufferIndex]);

        // if we are half full on the buffer, take an intermidiate IQ sample, for timing sync later
        static double complex IQmidpoint = 0;
        if(bufferIndex == symbolPeriod / 2 - 1 & n - tookSampleAt > 1)
        {
            IQmidpoint = 0;
            // compute DFT mid symbol (if time is already synced, otherwise it will help sync the time)
            for(int i = 0; i < symbolPeriod; i++)
            {
                // phasor offsets the cos and sin waves so that they allign with the time sequence of data in the half overwritten buffer
                double phase = (double)((i + bufferIndex + 1)%symbolPeriod) / symbolPeriod;
                IQmidpoint += sampleBuffer[i] * cos(2*M_PI*phase) + I * sampleBuffer[i] * sin(2*M_PI*phase);
            }
            // normalization factor
            IQmidpoint *= sqrt(1. / symbolPeriod);
            //printf("midpoint IQ: %f + %fi\n", creal(IQmidpoint), cimag(IQmidpoint));
        }
        // if the window buffer is filled, ie, we're on the last index of the buffer
        static double complex IQ = 0;
        static double complex IQlast = 0;
        if(bufferIndex == symbolPeriod - 1 & n - tookSampleAt > symbolPeriod / 2)
        {
            tookSampleAt = n;
            IQlast = IQ;
            IQ = 0;
            // compute DFT (rn just one freq, but when we implement OFDM, then at many harmonic(orthogonal) frequencies)
            // here is a simple quadrature detector
            for(int i = 0; i < symbolPeriod; i++)
            {
                // use a summation over the symbol period to separate orthogonal components (quadrature components) at the single frequency
                double phase = (double)i/symbolPeriod;  // one cycle per symbol period
                IQ += sampleBuffer[i] * cos(2*M_PI*phase) + I * sampleBuffer[i] * sin(2*M_PI*phase);
            }
            // normalization with a unitary normalization factor
            IQ *= sqrt(1. / symbolPeriod);

            // try to get a phase lock, symbol time lock, frequency lock, and equalization
            // calculate the error signal
            //double error = (pow(IQ, 2) - pow(IQlast, 2)) * pow(IQmidpoint, 2);
            double phaseOffsetEstimate = creal((IQ - IQlast) * conj(IQmidpoint));
            // PI control filter
            static double rollingAverage = 0;
            float weight = 0.1;
            rollingAverage = rollingAverage * (1 - weight) + phaseOffsetEstimate * weight;
            double error = 0 - rollingAverage;
            static double errorIntegral = 0;
            errorIntegral += error;
            double phaseAdjustment = (errorIntegral * 0.0 + error * 8) * 1;
            windowPhaseReal += phaseAdjustment;
            //windowPhaseReal = fmod(windowPhaseReal + phaseAdjustment, symbolPeriod);
            //windowPhaseReal = windowPhaseReal < 0 ? symbolPeriod + windowPhaseReal : windowPhaseReal;
            //windowPhase = (int)round(windowPhaseReal) % symbolPeriod;
            windowPhase = (int)round(windowPhaseReal) % symbolPeriod; // quantize the real window phase
            windowPhase = windowPhase < 0 ? symbolPeriod + windowPhase : windowPhase;
            //windowPhase = n * 2 / 2000;
            //windowPhase = (int)((windowPhase + phaseAdjustment) < 0 ? (symbolPeriod - windowPhase + phaseAdjustment) : (windowPhase + phaseAdjustment)) % symbolPeriod;

            // extract the frequencies to be decoded
            // add the relevant IQ values to an array
            // plot the IQ and a tail with some kind of persistance to give an animated output
            //fprintf(IQplotStdin, "%f %f\n", I, Q);
            //printf("%f %f\n", I, Q);
            // plot with a new color for each window phase
            fprintf(IQplotStdin, "%f %i %f\n", creal(IQ), windowPhase, cimag(IQ));
            //fprintf(IQplotStdin, "%f %i %f\n", creal(IQmidpoint), windowPhase + symbolPeriod, cimag(IQmidpoint));
            fprintf(errorPlotStdin, "%i, %f %f %f %f %f\n", n / symbolPeriod, phaseOffsetEstimate, error, phaseAdjustment, (double)windowPhase / symbolPeriod, windowPhaseReal / symbolPeriod);
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

    
    // after all the data is recieved, generate a plot of IQ 
    return 0;
}
