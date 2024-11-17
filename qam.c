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

#define WARN_UNUSED __attribute__((warn_unused_result))

#define DEBUG_LEVEL 0

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


typedef struct
{
    double I;
    double Q;
} iqsample_t;


// some functions to generate IQ streams with different properties
iqsample_t alternateI(int symbolIndex)
{
    iqsample_t sample = {symbolIndex % 2 * 2 - 1, 0};
    //iqsample_t sample = {symbolIndex % 2, 0};
    return sample;
}

iqsample_t randomQAM(int symbolIndex, int square)
{
    // square is the number of states of I and of Q, total states is square squared
    iqsample_t sample = {0, 0};
    sample.I = (double)(rand() % square) / (square - 1) * 2 - 1;
    sample.Q = (double)(rand() % square) / (square - 1) * 2 - 1;
    return sample;
}

iqsample_t randomQAM_withPreamble(int symbolIndex, int square)
{
    iqsample_t sample = {0, 0};
    srand(symbolIndex); // to get a stable random sequence everytime a specific symbol is requested
    int preamble = 150; // length of preamble in symbols
    if (symbolIndex < preamble)
    {
        // add a preamble that's easy to get rough time sync to
        sample = alternateI(symbolIndex);
    } else { // choose a new random QAM IQ value at start of every total period
        // then start sending random data
        sample = randomQAM(symbolIndex - preamble, square);
    }
    return sample;
}

iqsample_t sequentialIQ(int symbolIndex, int square)
{
    iqsample_t symbol;
    // sequentially hit all the IQ values in order in the constelation defined by power
    symbol.I = (double)(symbolIndex % square) / (square - 1) * 2 - 1;
    symbol.Q = (double)(symbolIndex / square % square) / (square - 1) * 2 - 1;
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
    int symbolPeriod = sampleRate / carrierFrequency * k; // audio samples per symbol
    int filterSides = 10;    // number of symbols to either side of current symbol to filter with raised cos
    int filterLengthSymbols = 2 * filterSides + 1;    // length of raised cos filter in IQ symbols, ie, how many IQ samples we need to generate the current symbol
    int filterLength = filterLengthSymbols * symbolPeriod;  // length in audio samples


    // phase offset that cycles through sequentially all phase offsets
    int phaseOffset = n *  2 / 1200 % symbolPeriod;
    //int phaseOffset = 
    n+=phaseOffset;

    // concept:
    //  generate the raised cos filter data once
    //  generate IQ samples ahead of time, just in time
    //  generate audio samples based on those two pieces of info

    static int initialized = 0;
    if(!initialized)
    {
        // initialize raised cos filter data
        filter = malloc(sizeof(double) * filterLength); // this never gets released. so, might wanna fix that TODO
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
        IQdata = malloc(sizeof(iqsample_t) * filterLengthSymbols);  // TODO this is never freed. These should prob be passed in as a paramter and freed somewhere in the larger scope
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
                IQdata[IQIndex % filterLengthSymbols] = alternateI(IQIndex);
                //IQdata[IQIndex % filterLengthSymbols] = randomQAM(IQIndex, 2);
                //IQdata[IQIndex % filterLengthSymbols] = sequentialIQ(IQIndex, 4);
                //IQdata[IQIndex % filterLengthSymbols] = randomQAM_withPreamble(IQIndex, 2);
            }
        }

        // ensure initialize only runs once
        initialized = 1;
    }

    // current symbol number
    int symbolIndex = n / symbolPeriod;
    int sampleIndex = n % symbolPeriod; // index of each sample in a symbol, where as n is increasing for the whole signal
    //int sampleIndex = fmod(n, symbolPeriod); // index of each sample in a symbol, where as n is increasing for the whole signal
    // circular buffer index
    int IQsampleIndex = symbolIndex % filterLengthSymbols;    // IQdata[IQsampleIndex] is current IQ sample, indexes above are future IQdata and below past IQ data both wrapping around until the are filterLengthSymbols / 2 away from current sample index

    // generate the next future IQ sample when entering a new symbol period
    //IQdata[(IQsampleIndex + filterLengthSymbols / 2) % filterLengthSymbols] = alternateI(symbolIndex + filterLengthSymbols / 2);
    if(sampleIndex == 0)
    {
        IQdata[(IQsampleIndex + filterLengthSymbols / 2) % filterLengthSymbols] = alternateI(symbolIndex + filterLengthSymbols / 2);
        //IQdata[(IQsampleIndex + filterLengthSymbols / 2) % filterLengthSymbols] = randomQAM(symbolIndex + filterLengthSymbols / 2, 2);
        //IQdata[(IQsampleIndex + filterLengthSymbols / 2) % filterLengthSymbols] = sequentialIQ(symbolIndex + filterLengthSymbols / 2, 4);
        //IQdata[(IQsampleIndex + filterLengthSymbols / 2) % filterLengthSymbols] = randomQAM_withPreamble(symbolIndex + filterLengthSymbols / 2, 2);
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
        filteredIQsample.I += filter[filterIndex] * IQdata[IQIndex].I;
        filteredIQsample.Q += filter[filterIndex] * IQdata[IQIndex].Q;
        //filteredIQsample.I += filter[filterIndex];
        //filteredIQsample.I = IQdata[IQIndex].I;
        //filteredIQsample.I += IQdata[IQIndex].I / filterLengthSymbols;
    }
    filteredIQsample.I /= 2;
    filteredIQsample.Q /= 2;

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
        (filteredIQsample.I) * cos(2.0 * M_PI * sampleIndex * k / symbolPeriod) +
        (filteredIQsample.Q) * sin(2.0 * M_PI * sampleIndex * k / symbolPeriod)
    ) / 2.0 * sqrt(2.0);

   return audioSample;

}

// really this is generating a single OFDM channel without guard periods
static double singleChannelODFM_noguard(int n, int sampleRate)
{
    int symbolPeriod = 256;
    //int guardPeriod = 4./1000 * sampleRate;     // I've found that the echos in a room last for about 3ms, durring that period, the symbol is phase offset and otherwise changed due to the last symbol and the transition between symbols
    int guardPeriod = 0;    // have to disable for QAM, ie, set to 0, instead use a raised cosine filter for ISI combat
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
    int symbols = pow(2, power);
    int square = sqrt(symbols);

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
        (sample.I + randI) * cos(2.0 * M_PI * symbolStep * k / symbolPeriod) +
        (sample.Q + randQ) * sin(2.0 * M_PI * symbolStep * k / symbolPeriod)
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
    return raisedCosQAM(n, sampleRate) * amplitudeScaler;
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
    if (hexdumpInput != NULL)
    {
        pclose(hexdumpInput);
    }
#endif

    ssize_t rets = write(fileDescriptor, header.bytes, sizeof(riff_header_t));

    if (rets < 0)
    {
        return -1;
    }

    return 0;
}

typedef struct __attribute__((packed))
{
    union
    {
        int32_t value;
        uint8_t bytes[sizeof(int32_t)];
    };
} int32_to_bytes_t;

static int WARN_UNUSED generateSamplesAndOutput(char* filenameInput)
{
    int retval = 0;

    FILE* aplayStdIn = NULL;

#if DEBUG_LEVEL >= 1
    FILE* hexdumpStdIn = NULL;
#endif

    int fileDescriptor = -1;

    // audio sample rate
    int sampleRate = 44100;
    // total number of samples to generate
    long length = sampleRate * 5;
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
                    putc(byte, hexdumpStdIn);
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
    if (outputstd == 0)
    {
    #if DEBUG_LEVEL >= 1
        pclose(hexdumpStdIn);
    #endif

        close(fileDescriptor);
    }

    pclose(aplayStdIn);

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
