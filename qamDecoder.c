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
#define DEBUG_LEVEL 2

typedef struct __attribute__((packed))
{
    union
    {
        int32_t value;
        uint8_t byte[sizeof(int32_t)];
    };
} sample32_t;

typedef struct
{
    // array storing super sampled IQ values
    // all the IQ samples stored in these bins.
    //      0 - just before ideal midpoint
    //      1 - just after ideal midpoint
    //      2 - just before ideal
    //      3 - just after ideal
    double complex data[4];

    // interpolated value between two closest actual IQ samples

    // between ideal sample times
    double complex midpoint;
     // previous ideal sample time
    double complex last;
    // at ideal sample time
    double complex iq;
    // to help the if statements that choose when to take an OFDM DFT sample not get stuck in a loop when the window offset changes
    int sampleTakenAt;
} iq_sample_t;

typedef struct
{
    double *buffer;
    int bufferIndex;
    int windowSize;
    int n;
    int k;
    int windowPhase;
    double carrierPhase;
    iq_sample_t sample;
} dft_t;

typedef struct
{
    FILE* waveformPlotStdin;
    FILE* fftDebuggerStdin;
    FILE* errorPlotStdin;
    FILE* IQplotStdin;
    FILE* IQvstimeStdin;
    FILE* eyeDiagramRealStdin;
    FILE* eyeDiagramImaginaryStdin;
} output_files_t;

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
double complex dft(dft_t *config, double* rmsOut, dft_debug_t debugFlag, FILE* debug_fd)
{
#if DEBUG_LEVEL <= 1
    (void)debugFlag;
    (void)debug_fd;
#endif

    // compute DFT (rn just one freq (k), but when we implement OFDM, then at many harmonic(orthogonal) frequencies)
    // here is a simple quadrature detector
    double complex IQ = 0;
    double RMS = 0;

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
        bufferIndex = (i + config->bufferIndex) % config->windowSize;

        // phase of the complex exponential
        // phasor offsets the cos and sin waves so that their phase is alligned with the real time of the samples, offset by the carrierPhase
        phase = (double)(n) * config->k / config->windowSize + config->carrierPhase;

        // use a summation over the symbol period to separate orthogonal components (quadrature components) at the single frequency
        double complex wave = cexp(I * 2 * M_PI * phase); // generate a sample from the complex exponential
        double complex value = config->buffer[bufferIndex] * wave;  // multiply the complex exponential by the input sample
        IQ += value;    // integrate the result over the window

        // compute RMS amplitude for equalization -- this kinda sucks. if there is a DC bias (ie, some low frequency interference)
        // it also comes through on RMS. Gotta fix that
        RMS += pow(config->buffer[i], 2);


        // this debug define simplifies the function a bit if debugging is disabled
    #if DEBUG_LEVEL > 1
        switch(debugFlag)
        {
            case MIDPOINT:
            {
                // debug graph outputs

                // debugging the phase allignment
                fprintf(debug_fd, "%i %i %f %i %f %i %f\n", config->n + i, 5, config->buffer[bufferIndex] + 4, 6, creal(wave) + 4, 7, cimag(wave) + 4);

                // debugging the integral
                //fprintf(debug_fd, "%i %i %f %i %f %i %f %i %f\n", debug_n + i, 11, creal(value) + 6, 12, cimag(value) + 6, 13, creal(IQ) + 6, 14, cimag(IQ) + 6);
                break;
            }

            case ALLIGNED:
            {
                // debug graph outputs

                // debugging the phase allignment
                fprintf(debug_fd, "%i %i %f %i %f %i %f\n", config->n + i, 0, config->buffer[bufferIndex], 1, creal(wave), 2, cimag(wave));
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
    RMS = sqrt(1.0 / config->windowSize * RMS);

#if DEBUG_LEVEL > 1
    switch(debugFlag)
    {
        case MIDPOINT:
        {
            // debug fft plot
            fprintf(debug_fd, "%i %i %f %i %f\n", config->n, 8, creal(IQ) + 4, 9, cimag(IQ) + 4);
            fprintf(debug_fd, "%i %i %f %i %f\n", config->n + config->windowSize, 8, creal(IQ) + 4, 9, cimag(IQ) + 4);
            break;
        }

        case ALLIGNED:
        {
            // debug fft plot
            fprintf(debug_fd, "%i %i %f %i %f\n", config->n, 3, creal(IQ), 4, cimag(IQ));
            fprintf(debug_fd, "%i %i %f %i %f\n", config->n + config->windowSize, 3, creal(IQ), 4, cimag(IQ));
            break;
        }
        case NODEBUG:
            break;
    }
#endif

    *rmsOut = RMS;
    return IQ;
}

int32_t getSample32Stdin(void)
{
    sample32_t sample = {0};

    // recieve data on stdin, signed 32bit integer
    for(size_t i = 0; i < sizeof(sample.value); i++)
    {
        // get the bytes, little endian, of signed 32bit integer, using the struct to do type punning
        sample.byte[i] = getchar();
    }

    return sample.value;
}

int takeSample(dft_t *config, double windowPhaseReal, double *rmsOut, FILE *debug_fd)
{
    // deciding when to take IQ samples
    // take a sample right between the ideal samples which are taken at windowPhase
    // if we are half full on the buffer, take an intermidiate IQ sample, for timing sync later
    // compute DFT (rn just one freq, but when we implement OFDM, then at many harmonic(orthogonal) frequencies)

    double rmsSamples[4] = {0};

    int phaseLower = (int)floor(windowPhaseReal);
    int phaseUpper = (int)ceil(windowPhaseReal);

    // just before real window phase offset midpoint
    int beforePhaseMidpoint = (phaseLower + config->windowSize / 2 - 1) % config->windowSize;
    // just after real window phase offset midpoint;
    int afterPhaseMidpoint = (phaseUpper + config->windowSize / 2 - 1) % config->windowSize;
    // just before real window phase offset
    int beforePhase = (phaseLower + config->windowSize - 1) % config->windowSize;
    // just after real window phase offset
    int afterPhase = (phaseUpper + config->windowSize - 1) % config->windowSize;

    if (config->bufferIndex == beforePhaseMidpoint)
    {
        config->sample.data[0] = dft(config, &rmsSamples[0], MIDPOINT, debug_fd);
    }

    if (config->bufferIndex == afterPhaseMidpoint)
    {
        config->sample.data[1] = dft(config, &rmsSamples[1], MIDPOINT, debug_fd);
    }

    if (config->bufferIndex == beforePhase)
    {
        config->sample.data[2] = dft(config, &rmsSamples[2], ALLIGNED, debug_fd);
    }

    if (config->bufferIndex == afterPhase)
    {
        config->sample.data[3] = dft(config, &rmsSamples[3], ALLIGNED, debug_fd);
    }

    // I added another condition to help debounce. sometimes it takes many samples in a row due to changing window offset
    if (config->bufferIndex == (phaseUpper + config->windowSize) % config->windowSize)
    {
        if (config->n - config->sample.sampleTakenAt > config->windowSize / 2)
        {
            // interpolate IQ and IQmidpoint from actual IQ samples in IQsamples
            // just doing linear interpolation
            //      (slope) * interp_X + initial
            //      (final - initial) * interp_X + initial
            config->sample.midpoint = (config->sample.data[1] - config->sample.data[0]) * fmod(windowPhaseReal, 1.0) + config->sample.data[0];
            config->sample.last = config->sample.iq;
            config->sample.iq = (config->sample.data[3] - config->sample.data[2]) * fmod(windowPhaseReal, 1.0) + config->sample.data[2];

            config->sample.sampleTakenAt = config->n;

            double rms = 0;

            for (size_t i = 0; i < 4; i++)
            {
                rms += rmsSamples[i];
            }

            rms /= 4.0;
            *rmsOut = rms;

            return 0;
        }
    }

    return -1;
}

int openFiles(output_files_t *files)
{
    int retval = 0;

    char *iq_plot_buffer = NULL;

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
    files->waveformPlotStdin = popen(plot, "w");

    if (files->waveformPlotStdin == NULL)
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
    files->fftDebuggerStdin = popen(debugPlot, "w");

    if (files->fftDebuggerStdin == NULL)
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
    files->errorPlotStdin = popen(errorPlot, "w");

    if (files->errorPlotStdin == NULL)
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

    files->IQplotStdin = popen(iq_plot_buffer, "w");

    // Free buffer, even if popen failed
    free(iq_plot_buffer);
    iq_plot_buffer = NULL;

    if (files->IQplotStdin == NULL)
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
    files->eyeDiagramRealStdin = popen(eyeDiagramPlotReal, "w");

    char eyeDiagramPlotImaginary[300] = {0};
    sprintf(
        eyeDiagramPlotImaginary,
        "feedgnuplot "
        "--domain --dataid --lines --points --maxcurves %i "
        "--title \"Eye Diagram, Imaginary part\" "
        "--xlabel \"Time (Audio sample #)\" --ylabel \"Q\" ",
        3000);
    files->eyeDiagramImaginaryStdin = popen(eyeDiagramPlotImaginary, "w");

    // using it to plot the time domain signal
    files->IQvstimeStdin = popen(eyeDiagramPlotReal, "w");

    if (files->IQvstimeStdin == NULL)
    {
        fprintf(stderr, "Failed to create IQ vs Time plot: %s\n", strerror(errno));
        retval = 7;
        goto exit;
    }

exit:
    if (iq_plot_buffer != NULL)
    {
        free(iq_plot_buffer);
        iq_plot_buffer = NULL;
    }

    return retval;
}

void closeFiles(output_files_t *files)
{
    // close the streams, allowing the plots to render and open a window, then wait for them to terminate in separate threads.
    if((files->IQvstimeStdin != NULL) && (fork() == 0))
    {
        pclose(files->IQvstimeStdin);
    }

    if((files->IQplotStdin != NULL) && (fork() == 0))
    {
        pclose(files->IQplotStdin);
    }

    if((files->errorPlotStdin != NULL) && (fork() == 0))
    {
        pclose(files->errorPlotStdin);
    }

    if((files->waveformPlotStdin != NULL) && (fork() == 0))
    {
        pclose(files->waveformPlotStdin);
    }

    if((files->fftDebuggerStdin != NULL) && (fork() == 0))
    {
        pclose(files->fftDebuggerStdin);
    }

    if((files->eyeDiagramRealStdin != NULL) && (fork() == 0))
    {
        pclose(files->eyeDiagramRealStdin);
    }

    if((files->eyeDiagramImaginaryStdin != NULL) && (fork() == 0))
    {
        pclose(files->eyeDiagramImaginaryStdin);
    }
}

void plotFiles(
    dft_t *config,
    output_files_t *files,
    int windowPhase,
    double phaseOffsetEstimate,
    double error,
    double phaseAdjustment,
    double windowPhaseReal
)
{
    // extract the frequencies to be decoded
    // add the relevant IQ values to an array
    // plot the IQ and a tail with some kind of persistance to give an animated output
    //fprintf(files.IQplotStdin, "%f %f\n", I, Q);
    //printf("%f %f\n", I, Q);

    // plot with a new color for each window phase
    fprintf(
        files->IQplotStdin,
        "%f %i %f\n",
        creal(config->sample.iq),
        windowPhase,
        cimag(config->sample.iq)
    );

    // fprintf(
    //     files.IQplotStdin,
    //     "%f %i %f\n",
    //     creal(config->sample.midpoint),
    //     windowPhase + SYMBOL_PERIOD,
    //     cimag(config->sample.midpoint)
    // );

    fprintf(
        files->errorPlotStdin,
        "%i, %f %f %f %f %f\n",
        config->n / SYMBOL_PERIOD,
        -phaseOffsetEstimate,
        error,
        phaseAdjustment,
        (double)windowPhase / SYMBOL_PERIOD,
        windowPhaseReal / SYMBOL_PERIOD
    );

    // fprintf(
    //     files.IQvstimeStdin,
    //     "%i, %f, %f, %f, %f, %f, %f\n",
    //     config->n / SYMBOL_PERIOD % (2 * 3),
    //     creal(config->sample.iq),
    //     cimag(config->sample.iq),
    //     creal(config->sample.last),
    //     cimag(config->sample.last),
    //     creal(config->sample.midpoint),
    //     cimag(config->sample.midpoint)
    // );

    fprintf(
        files->IQvstimeStdin,
        "%f, %f, %f\n",
        (config->n / SYMBOL_PERIOD) % (2*3) + 0.0,
        creal(config->sample.iq),
        cimag(config->sample.iq)
    );

    fprintf(
        files->IQvstimeStdin,
        "%f, %f, %f\n",
        (config->n / SYMBOL_PERIOD) % (2*3) + 0.5,
        creal(config->sample.midpoint),
        cimag(config->sample.midpoint)
    );
}

int main(void)
{
    // the OFDM channel number, how many cycles per symbol
    int k = 4;

    // buffer is the length of the symbol period, so that symbols are orthogonal
    double sampleBuffer[SYMBOL_PERIOD] = {0.0};

    //int guardPeriod = 4./1000 * 44100;      // 4ms guard period for room echos
    //int guardPeriod = 0;
    //int totalPeriod = SYMBOL_PERIOD + guardPeriod;

    // the offset of sample buffer or fftwindow in time, to get symbol time sync
    int windowPhase = 0;

    //double windowPhaseReal = (double)rand() / RAND_MAX * SYMBOL_PERIOD; // the floating point window phase, to be quantized into windowPhase
    double windowPhaseReal = 0;

    // prepare file descriptors for output files and streams
    int retval = 0;

    // process arguments
    //  none right now

    output_files_t files = {0};
    retval = openFiles(&files);

    if (retval != 0)
    {
        goto exit;
    }

    dft_t iq =
    {
        .buffer = sampleBuffer,
        .bufferIndex = 0,
        .windowSize = SYMBOL_PERIOD,
        .n = 0,
        .k = k,
        .windowPhase = 0,
        .carrierPhase = 0.0,
        .sample =
        {
            {0}
        },
    };

    // I just realized this whole loop will not really work for QAM signals, I should be using
    // a continuous sin cos multiplier for IQ sampling, not a DFT.
    // and it's not quite right for OFDM since we don't have a guard period yet, but it's closer to
    // that than qam.
    // while there is data to recieve, not end of file -> right now just a fixed number of 2000
    for(int audioSampleIndex = 0; audioSampleIndex < SYMBOL_PERIOD * 2000; audioSampleIndex++)
    {
        // gonna do something hacky here I think. if we are in a guard period, I'm just not gonna process that sample at all, and not going to index it.
        // No I'm not. This is the point where QAM deviates from OFDM. In QAM, we use raised cosine filter or matched root raised cosine filters
        // to combat inter symbol interference (ISI), but in OFDM, we use guard periods. the two are not interchangable because guard periods
        // mess up the gardner algorithm's time synchronization method in QAM, and raised cosine filters mess up the orthogonality between OFDM channels.
        // raised cos filter also causes ISI in back to back OFDM symbols because it blends them together at the symbol edges. I need a totally different IQ sampling scheme for
        // QAM. Continuous IQ sampling rather than DFT sampling
        iq.n = audioSampleIndex;

        // use the windowphase to adjust the buffer index position
        //int bufferIndex = (n + windowPhase)%SYMBOL_PERIOD;
        // Let's not use window phase to adjust index, and just pass the window phase to dft like a sane person
        int bufferIndex = iq.n % SYMBOL_PERIOD;

        iq.bufferIndex = bufferIndex + 1;
        iq.windowPhase = (iq.n + 1) % SYMBOL_PERIOD;

        // get sample from stdin
        int32_t sample = getSample32Stdin();

        // used by the equalization PID loop to control the overall volume
        static double equalizationFactor = 1.0;
        double RMS = 0.0;

        // convert sample to a double ranged from -1 to 1
        sampleBuffer[bufferIndex] = (double)sample / INT32_MAX * equalizationFactor;

    #if DEBUG_LEVEL > 0
        // debug waveform plot
        fprintf(
            files.waveformPlotStdin,
            "%i %f %f\n",
            iq.n,
            sampleBuffer[bufferIndex],
            equalizationFactor
        );
    #endif

    #if DEBUG_LEVEL > 0
        // oversampling IQ values for a nice eye diagram to help debug stuff
        iq.windowPhase = (iq.n + 0) % SYMBOL_PERIOD;
        double complex oversampledIQ = dft(&iq, &RMS, NODEBUG, NULL);
        iq.windowPhase = (iq.n + 1) % SYMBOL_PERIOD;

        fprintf(
            files.eyeDiagramRealStdin,
            "%f %i %f\n",
            fmod((double)iq.n / SYMBOL_PERIOD, 4.0),
            iq.n / 4 / SYMBOL_PERIOD,
            creal(oversampledIQ)
        );

        fprintf(
            files.eyeDiagramImaginaryStdin,
            "%f %i %f\n",
            fmod((double)iq.n / SYMBOL_PERIOD, 4.0),
            iq.n / 4 / SYMBOL_PERIOD,
            cimag(oversampledIQ)
        );
    #endif

        if (takeSample(&iq, windowPhaseReal, &RMS, files.fftDebuggerStdin) == 0)
        {
            // now I'm doing a bunch of stuff that happens every IQ sample. This all happens in the timespan of a single audio sample, which is 1/symbolPeriod of the time between IQ samples that could be used, but whatever

            // averaging filter for the equalizer
            static double rmsaverageWindow[44100 * 2 / SYMBOL_PERIOD] = {0.0};
            int rmsaverageSize = 44100 * 2 / SYMBOL_PERIOD;
            int rmsaverageIndex = (iq.n / SYMBOL_PERIOD) % rmsaverageSize;

            rmsaverageWindow[rmsaverageIndex] = RMS;
            double rmsaverage = 0.0;

            for(int i = 0; i < rmsaverageSize; i++)
            {
                rmsaverage += rmsaverageWindow[i];
            }

            rmsaverage /= rmsaverageSize;
            equalizationFactor += (1.0 / sqrt(2) - rmsaverage) * 0.002;

            // PID for the equalizer, just proportional. with max
            equalizationFactor = fmax(fmin(equalizationFactor, 10000), -10000);

            fprintf(
                files.waveformPlotStdin,
                "%i %f %f %f\n",
                iq.n,
                sampleBuffer[bufferIndex],
                equalizationFactor,
                rmsaverage
            );

        #if DEBUG_LEVEL > 0
            // equalization factor needs to ignore low frequency signals that give DC offset
            equalizationFactor = 1.0;
            //equalizationFactor = 1./0.007;
        #endif

            // try to get a phase lock, symbol time lock, frequency lock, and equalization
            // calculate the error signal
            // Gardner Algorithm: Real Part( derivitive of IQ times conjugate of IQ)
            // boobs-alexis
            // basically, it's trying to estimate the error of zero crossings
            double phaseOffsetEstimate = creal((iq.sample.iq - iq.sample.last) * conj(iq.sample.midpoint));

        #if DEBUG_LEVEL > 1
            // imprint the phase offset estimate on the fftDebugger plot to help understand the relationship between the DFT and symbol sync algo
            fprintf(
                files.fftDebuggerStdin,
                "%i %i %f %i %f\n",
                iq.n,
                10,
                phaseOffsetEstimate + 2.0,
                15,
                (double)(iq.n % SYMBOL_PERIOD) / SYMBOL_PERIOD + 2.0
            );
        #endif

            // state of the timing lock to modulate PID params
            static timing_lock_t lockstate = NO_LOCK;

            // Process Variable (PV, ie phase estimate) filter
            // rolling average of phase offset estimate
            // this may need to be adjusted based on the state of symbol and phase lock achieved
            static double averageWindow[64] = {0};
            const int averageWindowArrayLength = 64;
            static int averageSize;    // adjusted based on lock state
            int averageIndex = (iq.n / SYMBOL_PERIOD) % averageWindowArrayLength;

            switch(lockstate)
            {
                case NO_LOCK:
                {
                    // smaller averaging window to achieve
                    //  faster error signal response
                    //  for more aggressive PID tune
                    averageSize = 2;
                    break;
                }

                case SYMBOL_LOCK:
                    // Fallthrough
                case PHASE_LOCK:
                {
                    // longer average window to help average out symbol transitions that are not zero crossings
                    //  PID tune must be shittier
                    averageSize = 64;
                    break;
                }
            }

            // add a sample and take the average of last averageSize samples
            averageWindow[averageIndex] = phaseOffsetEstimate;
            double average = 0;

            // step backwards from current index to averageSize indexed back
            for(int i = averageIndex; i > averageIndex - averageSize; i--)
            {
                if(i < 0)
                {
                    average += averageWindow[averageWindowArrayLength + i];
                }

                average += averageWindow[i];
            }

            average /= averageSize;

            // PID loop for symbol timing, ie aligning the fft window to the symbol transitions
            static double P_gain = 0.0;
            static double I_gain = 0.0;
            static double D_gain = 0.0;

            switch(lockstate)
            {
                case NO_LOCK:
                {
                    P_gain = 2.5;
                    I_gain = 0.003;
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

            double error = 0.0 - average;
            //error *= -1;
            static double errorIntegral = 0.0;
            errorIntegral += error;

            static double lastError = 0.0;
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
            windowPhase = 0;
            windowPhaseReal = 0.0001;
            //windowPhase = (int)((windowPhase + phaseAdjustment) < 0 ? (SYMBOL_PERIOD - windowPhase + phaseAdjustment) : (windowPhase + phaseAdjustment)) % SYMBOL_PERIOD;
        #endif

            plotFiles(
                &iq,
                &files,
                windowPhase,
                phaseOffsetEstimate,
                error,
                phaseAdjustment,
                windowPhaseReal
            );
        }

         // I think I'll use feedgnuplot, so pipe the IQ values into feedgnuplot
    }

exit:
    closeFiles(&files);

    // after all the data is recieved, generate a plot of IQ
    return retval;
}
