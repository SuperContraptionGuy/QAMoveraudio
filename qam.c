// compile: gcc -lm qam.c
#define _GNU_SOURCE
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <time.h>

#include <complex.h>

#define WARN_UNUSED __attribute__((warn_unused_result))

#define DEBUG_LEVEL 0
#if DEBUG_LEVEL >= 1
    FILE* hexdumpStdIn = NULL;
    FILE* plotStdIn = NULL;
#endif

typedef enum
{
    RIFF_TYPE_PCM = 1,
} riff_type_t;

typedef struct __attribute__((packed))
{
    union
    {
        struct
        {
            char riff[4];               // 0
            uint32_t size;              // 4
            char format[4];             // 8
            char chunk[4];              // 12
            uint32_t length;            // 16
            uint16_t type;              // 20
            uint16_t channels;          // 22
            uint32_t sampleRate;        // 24
            uint32_t dataRate;          // 28
            uint16_t blockSize;         // 32
            uint16_t bitsPerSample;     // 34
            char data[4];               // 36
            uint32_t chunkSize;         // 40
        };
        uint8_t bytes[44];
    };
} riff_header_t;


// this is a bit lame. should use complex datatype
typedef struct
{
    double InPhase;
    double Quadrature;
} iqsample_t;

typedef struct __attribute__((packed))
{
    union
    {
        int32_t value;
        uint8_t bytes[sizeof(int32_t)];
    };
} int32_to_bytes_t;

typedef struct
{
    int sampleRate;
    int channels;

    // carrier frequenncy?
    int baseband;   // if set to 1, IQ sampels are taken as real so OFDM symbols can be transmitted at baseband. Halves the number of channels. to even channels only

    int guardPeriod;
    int symbolPeriod;   // includes guard period and ODFM period
    int ofdmPeriod;     // just the OFDM symbol itslef

    int initialized;

    FILE* dataInput;    // pipe to read data from that will be sent of the channel

    // array length of channels representing the current entire OFDM symbol
    double complex *currentOFDMSymbol;


    struct
    {
        int long frameStart;
        int fieldIndex;     // index of the field within the frame
        enum
        {
            IDLE,
            ACTIVE,
        } frame;    // a series of fields

        int long fieldStart;
        int symbolIndex;    // index of the symbol within the field
        enum
        {
            // Active field types
            PREAMBLE,   // for frame detection and timing synchronization 
            DATA,       // data field, also contains pilot symbols interspersed

        } field;    // a series of OFDM symbols

        int long symbolStart;
        enum
        {
            // all otehr symbol states
            GUARD_PERIOD,   // cyclic prefix for obsorbing channel effects (ISI)
            OFDM_PERIOD,  // transmission of the actual OFDM symbol information
            

        } symbol;   // an atomic OFDM symbol


    } state;
} OFDM_state_t;


// some functions to generate IQ streams with different properties
iqsample_t impulse(int symbolIndex, int impulseTime)
{
    iqsample_t sample = {0, 0};
    if(symbolIndex == impulseTime)
        sample.InPhase = 1;
    return sample;
}

iqsample_t alternateI(int symbolIndex)
{
    iqsample_t sample = {symbolIndex % 2 * 2 - 1, 0};
    //iqsample_t sample = {symbolIndex % 2, 0};
    return sample;
}

iqsample_t alternateQ(int symbolIndex)
{
    iqsample_t sample = {0, symbolIndex % 2 * 2 - 1};
    //iqsample_t sample = {symbolIndex % 2, 0};
    return sample;
}

iqsample_t randomQAM(int square)
{
    // square is the number of states of I and of Q, total states is square squared
    // I and Q values returned between -1 and 1
    iqsample_t sample = {0, 0};
    sample.InPhase = (double)(rand() % square) / (square - 1) * 2 - 1;
    sample.Quadrature = (double)(rand() % square) / (square - 1) * 2 - 1;
    return sample;
}

iqsample_t randomQAM_withPreamble(int symbolIndex, int square)
{
    iqsample_t sample = {0, 0};
    srand((unsigned int)symbolIndex); // to get a stable random sequence everytime a specific symbol is requested
    int preamble = 150; // length of preamble in symbols
    if (symbolIndex < preamble)
    {
        // add a preamble that's easy to get rough time sync to
        sample = alternateI(symbolIndex);
    } else { // choose a new random QAM IQ value at start of every total period
        // then start sending random data
        sample = randomQAM(square);
    }
    return sample;
}

iqsample_t sequentialIQ(int symbolIndex, int square)
{
    iqsample_t symbol;
    // sequentially hit all the IQ values in order in the constelation defined by power
    symbol.InPhase = (double)(symbolIndex % square) / (square - 1) * 2 - 1;
    symbol.Quadrature = (double)(symbolIndex / square % square) / (square - 1) * 2 - 1;
    return symbol;
}

double raisedCosQAM(int n, int sampleRate)
{
    int carrierPeriod = 32;
    int k = 1; // cycles per period
    //double carrierFrequency = 5000;
    double carrierFrequency = (double)sampleRate / carrierPeriod;
    //int carrierFrequency = sampleRate / carrierPeriod;

    //int symbolPeriod = 64; // audio samples per symbol
    //double symbolPeriod = sampleRate / carrierFrequency * k; // audio samples per symbol
    // array to store time series of filter data
    static double *filter = NULL;
    // array to store timeseries of IQ samples
    static iqsample_t *IQdata = NULL;

    if (n < 0)
    {
        if (filter != NULL)
        {
            free(filter);
            filter = NULL;
        }

        if (IQdata != NULL)
        {
            free(IQdata);
            IQdata = NULL;
        }

        return NAN;
    }

    // This is a sorta bug. the fact that everything is based on these integers is an issue. I think it causes discretization of carrier
    // frequency and symbol periods, which means the output frequency is not what you put in. for example, 5000 Hz input turns out to be
    // about 5500 Hz actual
    int symbolPeriod = (int)(sampleRate / carrierFrequency * k); // audio samples per symbol
    int filterSides = 10;    // number of symbols to either side of current symbol to filter with raised cos
    int filterLengthSymbols = 2 * filterSides + 1;    // length of raised cos filter in IQ symbols, ie, how many IQ samples we need to generate the current symbol
    int filterLength = filterLengthSymbols * symbolPeriod;  // length in audio samples


    // phase offset that cycles through sequentially all phase offsets
    //int phaseOffset = n *  2 / 1200 % symbolPeriod;
    int phaseOffset = 0;
    
    /*
    static int phaseOffset = -1;
    if (phaseOffset == -1)
        phaseOffset = rand() % symbolPeriod;
    */
    /*
    static int phaseOffset = -1;
    if (phaseOffset == -1 || (n % (carrierPeriod * 30)) == 0)
        phaseOffset = rand() % symbolPeriod;
    */

    int originalN = n;  // for plotting
    n -= phaseOffset;
    n = n < 0 ? 0 : n - phaseOffset, 0;

    // concept:
    //  generate the raised cos filter data once
    //  generate IQ samples ahead of time, just in time
    //  generate audio samples based on those two pieces of info

    static int initialized = 0;
    if(!initialized)
    {
        // initialize raised cos filter data
        filter = malloc((unsigned int)(sizeof(double) * filterLength)); // this never gets released. so, might wanna fix that TODO
        for(int i = 0; i < filterLength; i++)
        {
            int filter_symbolIndex = i - filterLength / 2;    // should go -filterLength/2 -> 0 -> filterLength/2
            // raised cos filter math
            double b = 0.42;    // filter parameter beta
            double filterValue = sin(M_PI * filter_symbolIndex / symbolPeriod) / (M_PI * filter_symbolIndex / symbolPeriod) * (cos(M_PI * b * filter_symbolIndex / symbolPeriod)) / (1 - pow(2 * b * filter_symbolIndex / symbolPeriod, 2));
            if(!isfinite(filterValue))   // in case it's undefined, ie divide by zero case
                filterValue = (double)symbolPeriod / 2 / b;
            if(filter_symbolIndex == 0)
                filterValue = 1;    // the math gives a divide by zero at index 0

            filter[i] = filterValue;
        }

        // generate enough IQ samples for first audio sample
        IQdata = malloc((unsigned int)(sizeof(iqsample_t) * filterLengthSymbols));  // TODO this is never freed. These should prob be passed in as a paramter and freed somewhere in the larger scope
        for(int i = 0; i < filterLengthSymbols; i++)
        {
            int IQIndex = i - filterLengthSymbols / 2;    // shoud go from -filterlengthsymbols / 2 -> 0 -> filterlengthsymbols / 2
            iqsample_t sample = {0, 0};
            // fill the circular buffer in the right order so that index 0 turns out to be the first IQ sample
            if(IQIndex < 0)
            {
                IQIndex = filterLengthSymbols + IQIndex;    // make positive
                IQdata[IQIndex % filterLengthSymbols] = sample; // the negative time samples are 0
            } else {
                //IQdata[IQIndex % filterLengthSymbols] = impulse(IQIndex, 20);
                //IQdata[IQIndex % filterLengthSymbols] = alternateI(IQIndex);
                //IQdata[IQIndex % filterLengthSymbols] = alternateQ(IQIndex);
                //IQdata[IQIndex % filterLengthSymbols] = randomQAM(2);
                //IQdata[IQIndex % filterLengthSymbols] = sequentialIQ(IQIndex, 4);
                IQdata[IQIndex % filterLengthSymbols] = randomQAM_withPreamble(IQIndex, 2);
            }
        }

        // ensure initialize only runs once
        initialized = 1;
    }

    // current symbol number
    int symbolIndex = n / symbolPeriod;
    int sampleIndex = n % symbolPeriod; // index of each sample in a symbol, where as n is increasing for the whole signal
    //printf("phaseOffset: %i\tn: %i\tsymbolIndex: %i\tsampleIndex: %i\n", phaseOffset, n, symbolIndex, sampleIndex);
    //int sampleIndex = fmod(n, symbolPeriod); // index of each sample in a symbol, where as n is increasing for the whole signal
    // circular buffer index
    int IQsampleIndex = symbolIndex % filterLengthSymbols;    // IQdata[IQsampleIndex] is current IQ sample, indexes above are future IQdata and below past IQ data both wrapping around until the are filterLengthSymbols / 2 away from current sample index

    // generate the next future IQ sample when entering a new symbol period
    //IQdata[(IQsampleIndex + filterLengthSymbols / 2) % filterLengthSymbols] = alternateI(symbolIndex + filterLengthSymbols / 2);
    if(sampleIndex == 0)
    {
        int IQindex = (IQsampleIndex + filterLengthSymbols / 2) % filterLengthSymbols;
        //IQdata[IQindex] = impulse(symbolIndex + filterLengthSymbols / 2, 20);
        //IQdata[IQindex] = alternateI(symbolIndex + filterLengthSymbols / 2);
        //IQdata[IQindex] = alternateQ(symbolIndex + filterLengthSymbols / 2);
        //IQdata[IQindex] = randomQAM(2);
        //IQdata[IQindex] = sequentialIQ(symbolIndex + filterLengthSymbols / 2, 4);
        IQdata[IQindex] = randomQAM_withPreamble(symbolIndex + filterLengthSymbols / 2, 2);
        IQindex = IQsampleIndex;
    #if DEBUG_LEVEL > 0
        fprintf(plotStdIn, "%i %i %f %i %f\n", originalN, 6, IQdata[IQsampleIndex].InPhase, 7, IQdata[IQsampleIndex].Quadrature);
        fprintf(plotStdIn, "%i %i %i\n", originalN, 8, IQindex);
    #endif
    }

    // add up raised cos contributions from all samples in the IQdata array
    iqsample_t filteredIQsample = {0, 0};
    for(int i = 0; i < filterLengthSymbols; i++)
    //int i = symbolIndex % filterLengthSymbols;
    //i = filterLengthSymbols - i - 1;
    {
        int filterIndex = (filterLengthSymbols - i - 1) * symbolPeriod + sampleIndex;    // pick a filter index
        int IQIndex = (IQsampleIndex + i - filterLengthSymbols / 2) % filterLengthSymbols;
        if(IQIndex < 0)
            IQIndex = filterLengthSymbols + IQIndex;    // make sure index is positive and wraps backwards
        filteredIQsample.InPhase += filter[filterIndex] * IQdata[IQIndex].InPhase;
        filteredIQsample.Quadrature += filter[filterIndex] * IQdata[IQIndex].Quadrature;
        //filteredIQsample.I += filter[filterIndex];
        //filteredIQsample.I = IQdata[IQIndex].I;
        //filteredIQsample.I += IQdata[IQIndex].I / filterLengthSymbols;
    #if DEBUG_LEVEL > 0
        //fprintf(plotStdIn, "%i %i %f\n", originalN + i, 8 + symbolIndex, IQdata[IQIndex].I);
    #endif
    }
    // normalization for raised cos filter, prob not correct
    //filteredIQsample.I /= 2;
    //filteredIQsample.Q /= 2;
    filteredIQsample.InPhase *= M_SQRT1_2;    // normalization factor for QAM constellation with max I and Q of 1, converts to max magnitude of 1
    filteredIQsample.Quadrature *= M_SQRT1_2;

    //return filter[symbolIndex%filterLengthSymbols*symbolPeriod+sampleIndex];
    //return IQdata[(IQsampleIndex + n%filterLengthSymbols - filterLengthSymbols / 2) %filterLengthSymbols].I;
    //return IQdata[IQsampleIndex].I;
    //return filteredIQsample.I;
    //return (double)i / filterLengthSymbols;
    //return filter[n%filterLength];

    //if(n % symbolPeriod < filterLengthSymbols)
        //return IQdata[(IQsampleIndex + sampleIndex) % filterLengthSymbols].I;
    //return -0.5;

    //filteredIQsample.I = 1;
    //filteredIQsample.Q = 0;

    double audioSample =
    (
        (filteredIQsample.InPhase) * cos(2.0 * M_PI * sampleIndex * k / symbolPeriod) +
        (filteredIQsample.Quadrature) * sin(2.0 * M_PI * sampleIndex * k / symbolPeriod)
    ) / 2.0;    // normalization factor for sum of sin+cos of amp 1 not exceeding 1
#if DEBUG_LEVEL > 0
    fprintf(plotStdIn, "%i %i %i %i %i %i %f %i %i\n", originalN, 0, symbolIndex, 1, sampleIndex, 2, audioSample, 3, phaseOffset);
    fprintf(plotStdIn, "%i %i %f %i %f\n", originalN, 4, filteredIQsample.InPhase, 5, filteredIQsample.Quadrature);
#endif

   return audioSample;

}

// generate IQ samples corrisponding to one OFDM symbol for baseband transmission (real samples valued only)
double OFDMsymbolBaseband(int n, OFDM_state_t *OFDMstate)
{
    // take the inverse transform of the current OFDM symbol values at the current offset n
    double sample = 0;
    for(int k = 0; k < OFDMstate->channels; k++)
    {
        // add together effectively a cos and sin scaled by the symbol amplitudes at each frequency (channel)
        sample += creal(OFDMstate->currentOFDMSymbol[k] * cexp(I*2*M_PI * k * n / OFDMstate->ofdmPeriod));
    }
    sample /= OFDMstate->channels;  // normalization constant

    return sample;
}

// runs the state machine for an OFDM based communication transmitter
double OFDM(int long n, OFDM_state_t *OFDMstate)
{
    if(OFDMstate->initialized == 0)
    {
        // initialize the state variables
        OFDMstate->state.frame = IDLE;
        OFDMstate->state.frameStart = n;

        OFDMstate->sampleRate = 44100;
        //OFDMstate->guardPeriod = 0.13 * OFDMstate->sampleRate;  // measured impulse response of room lasts like a bit over a tenth of a second would be nice to have this adjusted on the fly
        OFDMstate->guardPeriod = 1<<12;  // 2^12=1024 closest power of 2 to the impulse response length, slightly shorter
        //OFDMstate->guardPeriod = 128;   // faster testing
        //OFDMstate->guardPeriod = 0.001 * OFDMstate->sampleRate;  // measured impulse response of room lasts like a bit over a tenth of a second would be nice to have this adjusted on the fly
        OFDMstate->ofdmPeriod = OFDMstate->guardPeriod * 4;   // dunno what the best OFDM period is compared to the guard period. I assume longer is better for channel efficiency, but maybe it's worse for noise? don't know
        //OFDMstate->ofdmPeriod = OFDMstate->guardPeriod * 8;   // dunno what the best OFDM period is compared to the guard period. I assume longer is better for channel efficiency, but maybe it's worse for noise? don't know
        OFDMstate->symbolPeriod = OFDMstate->guardPeriod + OFDMstate->ofdmPeriod;
        OFDMstate->channels = OFDMstate->ofdmPeriod / 2;  // half due to using real symbols (ie, not modulating to higher frequency carrier wave but staying in baseband)

        // initialize current symbol array
        OFDMstate->currentOFDMSymbol = calloc(OFDMstate->channels, sizeof(double complex));
        if(OFDMstate->currentOFDMSymbol == NULL)
            fprintf(stderr, "Couldn't allocate current OFDM symbol array");


        OFDMstate->initialized = 1;
    }

    double output = 0;

    switch(OFDMstate->state.frame)
    {
        case IDLE:

            // send nothing
            output = 0;

            // check for exit from IDLE frame
            if(n - OFDMstate->state.frameStart >= OFDMstate->symbolPeriod * 1 - 1) // example of state change based on timing
            {
                OFDMstate->state.frame = ACTIVE;
                OFDMstate->state.frameStart = n + 1;    // starts next index
                OFDMstate->state.fieldIndex = 0;
                OFDMstate->state.field = PREAMBLE;
                OFDMstate->state.fieldStart = n + 1;
                OFDMstate->state.symbol = GUARD_PERIOD;
                OFDMstate->state.symbolStart = n + 1;
                OFDMstate->state.symbolIndex = 0;
            }
            break;

        case ACTIVE:
            switch(OFDMstate->state.field)
            {
                case PREAMBLE:

                    output = 0.70; // testing only

                    switch(OFDMstate->state.symbol)
                    {
                        case GUARD_PERIOD:

                            if(n - OFDMstate->state.symbolStart == 0)
                            {
                                // set new symbol for testing
                                // generate a symetric symbol, with even frequency components only for syncronization
                                for(int k = 0; k < OFDMstate->channels; k++)
                                {
                                    if(OFDMstate->state.symbolIndex == 0)
                                    {
                                        if(k % 2 == 0)
                                        //if(k > 16 && k < 24 && k % 2 == 0)
                                        //if(k == OFDMstate->channels / 16)
                                        {
                                            OFDMstate->currentOFDMSymbol[k] = 
                                                rand() % 2 * 2 - 1 +
                                                I*(rand() % 2 * 2 - 1);
                                            //OFDMstate->currentOFDMSymbol[k] = I;    // testing
                                        } else {
                                            OFDMstate->currentOFDMSymbol[k] = 0;
                                        }
                                    } else if(OFDMstate->state.symbolIndex == 1)
                                    {
                                        OFDMstate->currentOFDMSymbol[k] = 
                                            rand() % 2 * 2 - 1 +
                                            I*(rand() % 2 * 2 - 1);
                                        //OFDMstate->currentOFDMSymbol[k] = 0;    // testing
                                    }
                                }
                            }

                            output = OFDMsymbolBaseband(OFDMstate->ofdmPeriod - OFDMstate->guardPeriod + (n - OFDMstate->state.symbolStart),
                                                        OFDMstate);

                            if(OFDMstate->state.symbolIndex == 0)
                                output *= M_SQRT2;  // scale factor for even only channels
                            //output = output * 0.5;

                            // check for exit form GUARD_PERIOD symbol
                            if(n - OFDMstate->state.symbolStart >= OFDMstate->guardPeriod - 1)
                            {
                                OFDMstate->state.symbol = OFDM_PERIOD;
                                OFDMstate->state.symbolStart = n + 1;
                            }
                            break;
                        case OFDM_PERIOD:

                            output = OFDMsymbolBaseband(n - OFDMstate->state.symbolStart,
                                                        OFDMstate);
                            if(OFDMstate->state.symbolIndex == 0)
                                output *= M_SQRT2;
                            //output = output;

                            // check for exit
                            if(n - OFDMstate->state.symbolStart >= OFDMstate->ofdmPeriod - 1)
                            {
                                OFDMstate->state.symbol = GUARD_PERIOD;
                                OFDMstate->state.symbolStart = n + 1;
                                OFDMstate->state.symbolIndex++;
                            }
                            break;
                    }

                    // check for exit from PREAMBLE field
                    if(n - OFDMstate->state.fieldStart >= OFDMstate->symbolPeriod * 2 - 1)   // for 2 symbol preamble
                    {
                        OFDMstate->state.field = DATA;
                        OFDMstate->state.fieldStart = n + 1;
                        OFDMstate->state.symbolIndex = 0;
                    }
                    break;
                case DATA:

                    output = 1; // testing only

                    switch(OFDMstate->state.symbol)
                    {
                        case GUARD_PERIOD:

                            if(n - OFDMstate->state.symbolStart == 0)
                            {
                                // set new symbol for testing
                                for(int k = 0; k < OFDMstate->channels; k++)
                                {
                                    if(k < OFDMstate->channels)
                                    {
                                        OFDMstate->currentOFDMSymbol[k] = 
                                            rand() % 2 * 2 - 1 +
                                            I*(rand() % 2 * 2 - 1);   // QPSK
                                        //OFDMstate->currentOFDMSymbol[k] = 
                                            //(double)(rand() % 4) / 2 - 1; // 4 level
                                            //rand() % 2 +
                                            //I*(rand() % 2); // on off key
                                        //OFDMstate->currentOFDMSymbol[k] = 
                                            //rand() % 3 - 1; // zero sometimes

                                    } else {
                                        OFDMstate->currentOFDMSymbol[k] = 0;
                                    }
                                }
                            }
                            output = OFDMsymbolBaseband(OFDMstate->ofdmPeriod - OFDMstate->guardPeriod + (n - OFDMstate->state.symbolStart),
                                                        OFDMstate);
                            //output = output * 0.5;

                            // check for exit
                            if(n - OFDMstate->state.symbolStart >= OFDMstate->guardPeriod - 1)
                            {
                                OFDMstate->state.symbol = OFDM_PERIOD;
                                OFDMstate->state.symbolStart = n + 1;
                            }
                            break;
                        case OFDM_PERIOD:

                            output = OFDMsymbolBaseband(n - OFDMstate->state.symbolStart,
                                                        OFDMstate);
                            //output = output;

                            // check for exit
                            if(n - OFDMstate->state.symbolStart >= OFDMstate->ofdmPeriod - 1)
                            {
                                OFDMstate->state.symbol = GUARD_PERIOD;
                                OFDMstate->state.symbolStart = n + 1;
                                OFDMstate->state.symbolIndex++;
                            }
                            break;
                    }

                    // check for exit from DATA field
                    if(n - OFDMstate->state.fieldStart >= OFDMstate->symbolPeriod * 10 - 1) // set number of symbols of data for example
                    {
                        OFDMstate->state.frame = IDLE;
                        OFDMstate->state.frameStart = n + 1;
                    }
                    break;
            }
            break;
    }

    return output;
}

// really this is generating a single OFDM channel without guard periods,
// missing some critical infrustructure
static double singleChannelODFM_noguard(int n, int sampleRate)
{
    int symbolPeriod = 256;
    int guardPeriod = 4./1000 * sampleRate;     // I've found that the echos in a room last for about 3ms, durring that period, the symbol is phase offset and otherwise changed due to the last symbol and the transition between symbols
    //int guardPeriod = 0;    // have to disable for QAM, ie, set to 0, instead use a raised cosine filter for ISI combat
    int totalPeriod = symbolPeriod + guardPeriod;
    int k = 16;      // this is effectively the OFDM channel number, how many cycles per sample period

    // generating offsets in time to test the frame time syncronizer in qamDecoder
    // phase offset that cycles through sequentially all phase offsets
    //int phaseOffset = n *  2 / 2000 % symbolPeriod;
    // no phase offset
    //int phaseOffset = 0;
    // random phase offset
    static int phaseOffset = -1;
    if (phaseOffset == -1)
        phaseOffset = rand() % symbolPeriod;

    // apply phase offset
    n += phaseOffset;

    iqsample_t sample;  // The IQ sample currently being worked on, temporal resolution to audio sample
    double audioSample; // the final audio sample
    int symbolIndex = n / totalPeriod;  // symbol number indexed from 0, temporal resolution to IQ sample

    // random IQ in constelation defined by power
    int power = 2;  // log base2 of number of symbols. number of symbols should also be a perfect square
    int symbols = (int)pow(2, power);
    int square = (int)sqrt(symbols);

    // get IQ samples from the IQ sampling function
    sample = randomQAM_withPreamble(symbolIndex, square);

    // amplitude adjustment
    //double totalAmplitude = 0.01;
    double totalAmplitude = 1;
    //double totalAmplitude = 0.5;

    // random noise injection to IQ value
    double randomness = 0.0;
    double randI = ((double)rand() / RAND_MAX * 2 - 1) * randomness;
    double randQ = ((double)rand() / RAND_MAX * 2 - 1) * randomness;

    // implementing guard period
    // turns out this fucks up the gardner algorithm. ONLY for OFDM, not for QAM. Will use this later probs. disable by setting guard period to 0
    // time index of each symbol to resolution of audio sample
    int symbolStep = n % totalPeriod - guardPeriod;   // should be guardPeriod -> symbolPeriod -> 0 -> symbolPeriod as n increases from 0 -> totalPeriod
    audioSample =
    (
        (sample.InPhase + randI) * cos(2.0 * M_PI * symbolStep * k / symbolPeriod) +
        (sample.Quadrature + randQ) * sin(2.0 * M_PI * symbolStep * k / symbolPeriod)
    ) / 2.0 * sqrt(2.0) * totalAmplitude;

   return audioSample;
}

// this is the point where samples are generated
static double WARN_UNUSED calculateSample(int n, int sampleRate)
{

    //double amplitudeScaler = 0.1;
    double amplitudeScaler = 1;
    /*
    if(n == 2500)
        return 1;
    return 0;
    */
    static OFDM_state_t OFDMstate = {0};
    return OFDM(n, &OFDMstate);
    //return raisedCosQAM(n, sampleRate) * amplitudeScaler;
    //return impulse(n % (int)(44100 * 0.13 * 1.5), 0).I;
    //return singleChannelODFM_noguard(n, sampleRate) * amplitudeScaler;
}

// generates a .wav header of 44 bytes long
// length is the number of samples in the file
static int WARN_UNUSED writeHeader(int length, int fileDescriptor)
{
    riff_header_t header =
    {
        .riff = "RIFF",
        // chunk size plus rest of file size I think (so exluding first 8 bytes of header)
        .size = length * 4 + sizeof(riff_header_t) - 8,
        .format = "WAVE",
        // fmt chunk
        .chunk = "fmt ",
        .length = 16,
        .type = RIFF_TYPE_PCM,
        .channels = 1,
        .sampleRate = 44100,
        .dataRate = 176400,
        .blockSize = 4,
        .bitsPerSample = 32,
        .data = "data",
        .chunkSize = length * 4,
    };

#if DEBUG_LEVEL >= 1
    // dummy test samples
    uint8_t dummy[length * 4];
    memset(dummy, 0, length * 4);

    FILE* hexdumpInput = popen("hexdump -C", "w");

    if (hexdumpInput == NULL)
    {
        fprintf(stderr, "Failed to open hexdump: %s\n", strerror(errno));
        goto exit;
    }

    size_t ret = fwrite(header.bytes, sizeof(riff_header_t), 1, hexdumpInput);

    if (ret != sizeof(riff_header_t))
    {
        fprintf(stderr, "Failed to write to hexdump: %s\n", strerror(errno));
        goto exit;
    }

    ret = fwrite(header.bytes, sizeof(riff_header_t), 1, stdout);

    if (ret != sizeof(riff_header_t))
    {
        fprintf(stderr, "Failed to write to stdout: %s\n", strerror(errno));
        goto exit;
    }

    ret = fwrite(dummy, length * 4, 1, stdout);

    if (ret != sizeof(riff_header_t))
    {
        fprintf(stderr, "Failed to write to stdout: %s\n", strerror(errno));
    }

exit:
    if (hexdumpInput != NULL && fork() == 0)
    {
        pclose(hexdumpInput);
        exit(0);
        return 0;
    }
#endif

    ssize_t rets = write(fileDescriptor, header.bytes, sizeof(riff_header_t));

    if (rets < 0)
    {
        return -1;
    }

    return 0;
}

static int WARN_UNUSED generateSamplesAndOutput(char* filenameInput)
{
    int retval = 0;

    FILE* aplayStdIn = NULL;


    int fileDescriptor = -1;

    // audio sample rate
    int sampleRate = 44100;
    // total number of samples to generate
    long length = sampleRate * 6.5;
    // the number of the current sample
    long n = 0;

    // length of the file write buffer, samples times 4 bytes per sample
    const int bufferLength = 100 * 4;
    // the file write buffer, used to buffer the write calls
    uint8_t buffer[bufferLength];
    // number of bytes ready to be written out of the buffer in case we need to flush the buffer before it's full
    int bufferReadyBytes = 0;

    // Whether to send samples over stdout or to file
    int outputstd = 0;

    // User passed '-' -> use stdout
    if (filenameInput[0] == '-')
    {
        outputstd = 1;
    }

    // set up the file descriptors for the various outputs

    // setup a file descriptor for a pipe to aplay command to play the sound through the speakers
    char aplayCommandString[30] = {0};
    int len = snprintf(aplayCommandString, 30, "aplay -f S32_LE -c1 -r %i", sampleRate);

    if (len < 0)
    {
        fprintf(stderr, "Failed to write aplay string: %s\n", strerror(errno));
        retval = 2;
        goto exit;
    }
    else if (len == 30)
    {
        fprintf(stderr, "Failed to write aplay string: truncated\n");
        retval = 2;
        goto exit;
    }

    //puts(aplayCommandString);
    aplayStdIn = popen(aplayCommandString, "w");

    if (aplayStdIn == NULL)
    {
        fprintf(stderr, "Failed to open aplay: %s\n", strerror(errno));
        retval = 3;
        goto exit;
    }

#if DEBUG_LEVEL >= 1
    if (outputstd == 0)
    {
        // for the hex dump of printed bytes.
        hexdumpStdIn = popen("hexdump -C", "w");

        if (hexdumpStdIn == NULL)
        {
            fprintf(stderr, "Failed to open hexdump: %s\n", strerror(errno));
            retval = 4;
            goto exit;
        }
    }
    char *plotstr = 
        "feedgnuplot "
        "--domain --dataid --lines --points --maxcurves 100000 "
        "--title \"Debug modulator\" "
        "--xlabel \"Time (n)\" --ylabel \"value\" "
        "--legend 0 \"Symbol Index\" "
        "--legend 1 \"Sample Index\" "
        "--legend 2 \"Generated Audio Sample\" "
        "--legend 3 \"Phase Offset\" "
        "--legend 4 \"Filtered I\" "
        "--legend 5 \"Filtered Q\" "
        "--legend 6 \"I\" "
        "--legend 7 \"Q\" "
    ;
    plotStdIn = popen(plotstr, "w");
#endif

    if (outputstd == 0)
    {
        // for the file writing
        char filename[30] = {0};
        len = snprintf(filename, 30, "%s.wav", filenameInput);

        if (len < 0)
        {
            fprintf(stderr, "Failed to get filename: %s\n", strerror(errno));
            retval = 5;
            goto exit;
        }
        else if (len == 30)
        {
            fprintf(stderr, "Failed to get filename: truncated\n");
            retval = 5;
            goto exit;
        }

        //puts(filename);

        fileDescriptor = open(filename, O_WRONLY | O_CREAT, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH);

        if (fileDescriptor < 0)
        {
            fprintf(stderr, "Failed to open file: %s\n", strerror(errno));
            retval = 6;
            goto exit;
        }
    }

    // first generate the header
    if (outputstd == 0)
    {
        int ret = writeHeader(length, fileDescriptor);

        if (ret != 0)
        {
            fprintf(stderr, "Failed to write header: %d\n", ret);
            retval = 7;
            goto exit;
        }
    }

    // calculate all the samples
    // seed the random number generator
    srand(time(NULL));
    while(n < length)
    {
        // calculate a chunk of samples until the buffer is full or max is reached. one sample at a time, 4 bytes at a time
        for(bufferReadyBytes = 0; (bufferReadyBytes < bufferLength) && (n < length); bufferReadyBytes += 4, n++)
        {
            // the sample value used in calculations, to be normalized
            double sampleValue;
            // sample value after put into signed integer range, then split into bytes for file writing and audio output
            int32_to_bytes_t normalizedSampleValue;
            // holds each individual byte as it's written out Little Endian style
            char byte;

            // get the double sample value, should be between -1 and 1
            sampleValue = calculateSample(n, sampleRate);
            // calculate the final signed integer to be output as a sample
            // the magnitude of the max is always one smaller than the magnitude of the min
            normalizedSampleValue.value = sampleValue * INT32_MAX;

            // split up the normalized value into individual bytes
            for(int i = 0; i < 4; i++)
            {
                // get the byte from normalized. pointer points to the adress of the first byte in normalized
                byte = normalizedSampleValue.bytes[i];
                // add byte to the buffer
                buffer[bufferReadyBytes + i] = byte;

                // send to the pipes one byte at a time since they are buffered by the OS
                putc(byte, aplayStdIn);

            #if DEBUG_LEVEL >= 1
                if (outputstd == 0)
                {
                #if DEBUG_LEVEL > 1
                    putc(byte, hexdumpStdIn);
                #endif
                }
                else
                {
                    putchar(byte);
                }
            #else
                if (outputstd != 0)
                {
                    putchar(byte);
                }
            #endif
            }
        }

        // write the buffer to the file bufferReadyBytes number of bytes, usually a whole buffer full at a time, until the end.
        if (outputstd == 0)
        {
            ssize_t ret = write(fileDescriptor, buffer, bufferReadyBytes);

            if (ret < 0)
            {
                fprintf(stderr, "Failed to write output buffer: %s\n", strerror(errno));
                retval = 8;
                goto exit;
            }
        }
    }

exit:
    pclose(aplayStdIn);
    if (outputstd == 0)
    {
        close(fileDescriptor);

    #if DEBUG_LEVEL > 0
        pclose(hexdumpStdIn);
    #endif
    }

#if DEBUG_LEVEL > 0
    if(plotStdIn != NULL && fork() == 0)
    {
        // it's holding onto a reference to stdout, gotta close that off
        //freopen("/dev/null", "w", stdout);
        // nope, that wasn't it
        pclose(plotStdIn);
        exit(0);
        return 0;
    }
#endif

    return retval;
}

static void usage(const char *filename)
{
    fprintf(stderr, "Provide file name as parameter. For example:\n");
    fprintf(stderr, "  Wrap the file name in quotes:\n");
    fprintf(stderr, "  ./%s \"file name\"\n", filename);
    fprintf(stderr, "Will generate a file called \"file name.wav\".\n");
}

int main(int argc, char** args)
{
    // extract filename from arguments
    if (argc < 2)
    {
        // or use stdout
        usage(args[0]);
        return 1;
    }
    else if (argc > 2)
    {
        fprintf(stderr, "Too many arguments.\n");
        usage(args[0]);
        return 1;
    }

    return generateSamplesAndOutput(args[1]);
}
