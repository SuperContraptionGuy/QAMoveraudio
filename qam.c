// compile: gcc -lm qam.c
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>
#include <math.h>



double simpleQAM(int n, double t)
{
    int symbolPeriod = 32;
    //uint8_t count = t * 100;
    long count = n / symbolPeriod;
    int power = 2;  // log base2 of number of symbols. number of symbols should also be a perfect square
    int symbols = pow(2, power);
    int square = sqrt(symbols);
    uint8_t mask = symbols - 1;
    //printf("symbols: %i, square: %i, mask: %x\n", symbols, square, mask);
    count = count & mask;

    static double oldI = 0;
    static double oldQ = 0;

    // determine the I and Q
    //double I = (double)(count % square) / (square - 1) * 2 - 1;
    //double Q = (double)(count / square) / (square - 1) * 2 - 1;
    //double I = count % 2 * 2 - 1;
    //double Q = (count + 1) % 2;
    //double Q = 0;

    static double I = 0;
    double Q = 0;
    if(n % symbolPeriod == 0)
    {
        I = (double)(rand() % 2) * 2 - 1;
    }

    // variables to enable transition IQ values
    static double decayStartTime = 0;
    static double decayStartI = 0;
    static double decayStartQ = 0;
    //double decayTime = (double)symbolPeriod / 4 / n * t;
    double decayTime = 0;

    // This really only works for descrete IQ steps that take place at intervals longer than decayTime
    if(oldI != I | oldQ != Q)
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
    if(t - decayStartTime < decayTime)
    {
        I = (I - decayStartI) / decayTime * (t - decayStartTime) + decayStartI;
        Q = (Q - decayStartQ) / decayTime * (t - decayStartTime) + decayStartQ;
    }

    double totalAmplitude = 0.01;
    //double totalAmplitude = 1;
    double randomness = 0;
    double randI = ((double)rand() / RAND_MAX * 2 - 1) * randomness;
    double randQ = ((double)rand() / RAND_MAX * 2 - 1) * randomness;
    return ((I + randI) * cos(2*M_PI*n/symbolPeriod) + (Q + randQ) * sin(2*M_PI*n/symbolPeriod))/2*sqrt(2) * totalAmplitude;
    //return (I * sin(2*M_PI*t*600) + Q * cos(2*M_PI*t*600))/2*sqrt(2);
    //return (I * sin(2*M_PI*t*6000) + Q * cos(2*M_PI*t*6000))/2*sqrt(2);
}





// this is the point where samples are generated
double calculateSample(int n, double t) {

    return simpleQAM(n, t);
}


// generates a .wav header of 44 bytes long
// length is the number of samples in the file
void writeHeader(int length, int fileDescriptor)
{
    char header[44];
    // pointers to parameters of length 2 bytes, 4 bytes, or signed 4 bytes.
    // used to write specific bytes in the header with parameter values
    uint16_t *param2;
    uint32_t *param4;
    int32_t *param4s;

    // RIFF chunk
    sprintf(header, "RIFF");

    // chunk size plus rest of file size I think (so exluding first 8 bytes of header)
    param4 = (uint32_t*)(header+4);
    *param4 = length*4+sizeof(header)-8;

    sprintf(header+8, "WAVE");

    // fmt chunk
    sprintf(header+12, "fmt ");

    // length of fmt chunk: 16
    param4 = (uint32_t*)(header+16);
    *param4 = 16;

    // format code: 1 PCM
    param2 = (uint16_t*)(header+20);
    *param2 = 1;

    // number of channels: 1
    param2 = (uint16_t*)(header+22);
    *param2 = 1;

    //sample rate: 44100
    param4 = (uint32_t*)(header+24);
    *param4 = 44100;

    // Data rate: 176400 bytes/s
    param4 = (uint32_t*)(header+28);
    *param4 = 176400;

    // Data block size: 4 bytes
    param2 = (uint16_t*)(header+32);
    *param2 = 4;

    // bits per sample
    param2 = (uint16_t*)(header+34);
    *param2 = 32;

    // data chunk
    sprintf(header+36, "data");

    // data chunk size
    param4 = (uint32_t*)(header+40);
    *param4 = length*4;

    //FILE* hexdumpInput = popen("hexdump -C", "w");
    for(int i = 0; i < 44; i++)
    {
        //putc(header[i], hexdumpInput);
        //putchar(header[i]);
    }
    //dummy test samples
    for(int i = 0; i < length*4; i++)
    {
        //putchar(0);
    }
    //pclose(hexdumpInput);
    write(fileDescriptor, header, sizeof(header));

}

void generateSamplesAndOutput(char* filenameInput)
{
    int sampleRate = 44100;

    // set up the file descriptors for the various outputs

    // setup a file descriptor for a pipe to aplay command to play the sound through the speakers
    char aplayCommandString[30];
    sprintf(aplayCommandString, "aplay -f S32_LE -c1 -r %i", sampleRate);
    //puts(aplayCommandString);
    FILE* aplayStdIn = popen(aplayCommandString, "w");

    // for the hex dump of printed bytes.

    int outputstd = 0;  // send samples over stdout instead
    if(filenameInput[0] == '-')
        outputstd = 1;
    FILE* hexdumpStdIn = NULL;
    if(!outputstd)
        hexdumpStdIn = popen("hexdump -C", "w");


    // for the file writing
    char filename[30];
    sprintf(filename, "%s.wav", filenameInput);
    //puts(filename);
    int fileDescriptor;
    if(!outputstd)
        fileDescriptor = open(filename, O_WRONLY | O_CREAT, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH);


    long length = 100000;           // total number of samples to generate
    long n = 0;                 // the number of the current sample
    const int bufferLength = 10 * 4;    // length of the file write buffer, samples times 4 bytes per sample
    uint8_t buffer[bufferLength];   // the file write buffer, used to buffer the write calls
    int bufferReadyBytes = 0;   // number of bytes ready to be written out of the buffer in case we need to flush the buffer before it's full
    double value;           // the sample value used in calculations, to be normalized
    int32_t normalized = 0;    // sample value after put into signed integer range
    //int32_t max = INT32_MAX > -(long)INT32_MIN ? -INT32_MIN : INT32_MAX;    // maximum signed integer value used for normalization
    char byte;                  // holds each individual byte as it's written out Little Endian style
    unsigned char* pointer = (unsigned char*) &normalized;  // pointer used for breaking up the normalized value into bytes
    

    // first generate the header
    if(!outputstd)
        writeHeader(length, fileDescriptor);


    // calculate all the samples
    while(n < length)
    {
;
        // calculate a chunk of samples until the buffer is full or max is reached. one sample at a time, 4 bytes at a time
        for(bufferReadyBytes = 0; bufferReadyBytes < bufferLength & n < length; bufferReadyBytes+=4, n++)
        {
            // get the double sample value, should be between -1 and 1
            value = calculateSample(n, (double)n / sampleRate);

            // calculate the final signed integer to be output as a sample
            normalized = value * INT32_MAX; // the magnitude of the max is always one smaller than the magnitude of the min
            // split up the normalized value into individual bytes
            for(int i = 0; i < 4; i++)
            {
                // get the byte from normalized. pointer points to the adress of the first byte in normalized
                byte = *(pointer+i);
                //byte = n * 4 + i;     // labeling the bytes to make sure all is well

                // add byte to the buffer
                buffer[bufferReadyBytes+i] = byte;
            
                // send to the pipes one byte at a time since they are buffered by the OS
                putc(byte, aplayStdIn);
                if(!outputstd)
                {
                    //putc(byte, hexdumpStdIn);
                } else {
                    putchar(byte);
                }
            }
        }
        // write the buffer to the file bufferReadyBytes number of bytes, usually a whole buffer full at a time, until the end.
        if(!outputstd)
            write(fileDescriptor, buffer, bufferReadyBytes);
    }

    if(!outputstd)
    {
        pclose(hexdumpStdIn);
        close(fileDescriptor);
    }
    pclose(aplayStdIn);
}

int main(int argc, char** args) {

    // extract filename from arguments
    if(argc == 1)
    {
        // or use stdout
        puts("Provide file name as parameter. For example:");
        puts("  ./a.out \"file name\"");
        puts("will generate a file called \"file name.wav\".");
        return 2;
    } else if(argc > 2) {
        puts("Too many arguments.");
        puts("Accepts one argument, a file name.");
        puts("wrap the file name in quotes, like this:");
        puts("  ./a.out \"file name\"");
        return 2;
    }

    generateSamplesAndOutput(args[1]);

    return 0;
}
