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

static double WARN_UNUSED simpleQAM(int n, double t)
{
    int symbolPeriod = 64;
    int k = 4;      // this is effectively the OFDM channel number, how many cycles per sample period
    //uint8_t count = t * 100;

    // generating offsets in time to test the frame time syncronizer in qamDecoder
    //int phaseOffset = n * 4 * 2 / symbolPeriod / 2000 % 4 * symbolPeriod / 4;
    //int phaseOffset = n * 8 * 2 / symbolPeriod / 2000 % 8 * symbolPeriod / 8;
    //int phaseOffset = n *  2 / 2000 % symbolPeriod;
    //int phaseOffset = 3 * symbolPeriod / 4 + 1;
    int phaseOffset = 0;

    /*
    // random phase offset
    static int phaseOffset = -1;
    if (phaseOffset == -1)
        phaseOffset = rand() % symbolPeriod;
        */

    n += phaseOffset;

    long count = n / symbolPeriod;
    int power = 2;  // log base2 of number of symbols. number of symbols should also be a perfect square
    int symbols = pow(2, power);
    // int square = sqrt(symbols);
    uint8_t mask = symbols - 1;
    //printf("symbols: %i, square: %i, mask: %x\n", symbols, square, mask);
    count = count & mask;

    static double oldI = 0;
    static double oldQ = 0;

    // determine the I and Q
    //double I = count % 2 * 2 - 1;
    //double Q = count % 2 * 2 - 1;
    //double Q = count / 2 * 2 - 1;
    //double I = count % 2;
    //double Q = (count + 1) % 2;
    //double I = 0;
    //double Q = 0;

    // sequentially hit all the IQ values in order in the constelation defined by power
    //double I = (double)(count % square) / (square - 1) * 2 - 1;
    //double Q = (double)(count / square) / (square - 1) * 2 - 1;

    // random IQ in constelation defined by power
    static double I = 0;
    static double Q = 0;
    if (n / symbolPeriod < 150)
    {
        // add a preamble that's easy to get rough time sync to
        I = count % 2 * 2 - 1;
    } else if (n % symbolPeriod == 0) {
        // then start sending random data
        //I = ((double)(rand() % 2) * 2 - 1) / 1;
        //Q = ((double)(rand() % 2) * 2 - 1) / 1;
        I = (double)(rand() % square) / (square - 1) * 2 - 1;
        Q = (double)(rand() % square) / (square - 1) * 2 - 1;
    }

    // variables to enable transition IQ values
    static double decayStartTime = 0;
    static double decayStartI = 0;
    static double decayStartQ = 0;
    //double decayTime = (double)symbolPeriod / 4 / n * t;
    double decayTime = 0;

    // This really only works for descrete IQ steps that take place at intervals longer than decayTime
    if (oldI != I || oldQ != Q)
    {
        // this will switch to the old goal, not the old actual value, could cause discontinuities if the old value is not yet equal to the old goal
        decayStartTime = t;
        decayStartI = oldI;
        decayStartQ = oldQ;
    }

    // representing the target values and old target values, before decay
    oldI = I;
    oldQ = Q;

    // linear interpolation, if the decay time has not elapsed
    if (t - decayStartTime < decayTime)
    {
        I = (I - decayStartI) / decayTime * (t - decayStartTime) + decayStartI;
        Q = (Q - decayStartQ) / decayTime * (t - decayStartTime) + decayStartQ;
    }

    double totalAmplitude = 0.01;
    //double totalAmplitude = 1;
    double randomness = 0.0;
    double randI = ((double)rand() / RAND_MAX * 2 - 1) * randomness;
    double randQ = ((double)rand() / RAND_MAX * 2 - 1) * randomness;

    return ((I + randI) * cos(2.0 * M_PI * n * k / symbolPeriod) + (Q + randQ) * sin(2.0 * M_PI * n * k / symbolPeriod)) / 2.0 * sqrt(2.0) * totalAmplitude;

    //return (I * sin(2 * M_PI * t * 600) + Q * cos(2 * M_PI * t * 600))/2 * sqrt(2);
    //return (I * sin(2 * M_PI * t * 6000) + Q * cos(2 * M_PI * t * 6000))/2 * sqrt(2);
}

// this is the point where samples are generated
static double WARN_UNUSED calculateSample(int n, double t)
{
    return simpleQAM(n, t);
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

static int WARN_UNUSED generateSamplesAndOutput(char* filenameInput)
{
    int retval = 0;

    FILE* aplayStdIn = NULL;

#if DEBUG_LEVEL >= 1
    FILE* hexdumpStdIn = NULL;
#endif

    int fileDescriptor = -1;

    int sampleRate = 44100;

    // total number of samples to generate
    long length = 100000;

    // the number of the current sample
    long n = 0;

    // length of the file write buffer, samples times 4 bytes per sample
    const int bufferLength = 10 * 4;

    // the file write buffer, used to buffer the write calls
    uint8_t buffer[bufferLength];

    // number of bytes ready to be written out of the buffer in case we need to flush the buffer before it's full
    int bufferReadyBytes = 0;

    // the sample value used in calculations, to be normalized
    double value;

    // sample value after put into signed integer range
    int32_t normalized = 0;

    // maximum signed integer value used for normalization
    //int32_t max = INT32_MAX > -(long)INT32_MIN ? -INT32_MIN : INT32_MAX;

    // holds each individual byte as it's written out Little Endian style
    char byte;

    // pointer used for breaking up the normalized value into bytes
    unsigned char* pointer = (unsigned char*)&normalized;

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
    while(n < length)
    {
        // calculate a chunk of samples until the buffer is full or max is reached. one sample at a time, 4 bytes at a time
        for(bufferReadyBytes = 0; (bufferReadyBytes < bufferLength) && (n < length); bufferReadyBytes += 4, n++)
        {
            // get the double sample value, should be between -1 and 1
            value = calculateSample(n, (double)n / sampleRate);

            // calculate the final signed integer to be output as a sample

            // the magnitude of the max is always one smaller than the magnitude of the min
            normalized = value * INT32_MAX;

            // split up the normalized value into individual bytes
            for(int i = 0; i < 4; i++)
            {
                // get the byte from normalized. pointer points to the adress of the first byte in normalized
                byte = *(pointer + i);

                // labeling the bytes to make sure all is well
                //byte = n * 4 + i;

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
