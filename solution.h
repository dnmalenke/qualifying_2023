//! Rohde & Schwarz Engineering Competition 2023
//!
//! This is the code to speed up. Enjoy!

#pragma once

#include "ec2023/ec2023.h"
#include <iostream>      
#include <iomanip>
#include <vector>
#include <cassert>
#include <cstring>

#define REAL(x) x + 32 * i
#define IMAG(x) WINDOW_SIZE + x + 32 * i

#define SWEET_SPOT 8

#define INIT(x) if(!streamInitted){x;}

static bool streamInitted = false;

static constexpr float OVERLAP_RATIO = 0.75;
static constexpr size_t WINDOW_SIZE = 1024;
static const ec::Float PI = 3.14159265358979323846f;
static const ec::Float minusOne = -1.0f;
static const ec::Float zero = 0.0f;
static const ec::Float sqrt22 = (float)sqrt(2.0f) / 2.0f;

static const ec::Float spC0((float)(10.0f / log(10.0f)));
static const ec::Float spC1((float)(10.0f * log(125.0f / 131072.0f) / log(10.0f)));
static const ec::Float spC2((float)(10.0f * log(125.0f / 32768.0f) / log(10.0f)));

void init_blackmanCoefs();
void init_angleTerms();
uint16_t reverse_bits(uint16_t x);

void sfft(std::vector<ec::Float>& inputReal, std::vector<ec::Float>& inputImag);
void fft(std::vector<ec::Float>& inputReal, std::vector<ec::Float>& inputImag);

static std::vector<ec::Float> angleTerms(2 * WINDOW_SIZE);
static std::vector<ec::Float> blackmanCoefs(WINDOW_SIZE);

ec::Measurement* mes;


#define MEASURE(x,msg) { size_t _pre = mes->calcTotalScore(); x; size_t _post = mes->calcTotalScore(); std::cout << msg << " cost: " << _post - _pre  << std::endl; }

std::vector<ec::Float> process_signal(const std::vector<ec::Float>& inputSignal, ec::Measurement* m)
{

    mes = m;

    const size_t numSamples = inputSignal.size(); // assume divisible by 32
    const size_t sizeSpectrum = (WINDOW_SIZE / 2) + 1;
    const size_t stepBetweenWins = static_cast<size_t>(ceil(WINDOW_SIZE * (1 - OVERLAP_RATIO)));
    const size_t numWins = (numSamples - WINDOW_SIZE) / stepBetweenWins + 1;

    std::vector<ec::Float> signalWindow(WINDOW_SIZE);

    std::vector<ec::Float> outputSpectrum(sizeSpectrum);
    std::vector<ec::Float> preLogSpectrum(sizeSpectrum, std::numeric_limits<float>::lowest());

    size_t idxStartWin = 0;

    init_blackmanCoefs();
    init_angleTerms();

    ec::VecHw& vecHw = *ec::VecHw::getSingletonVecHw();
    vecHw.copyToHw(angleTerms, 0, 2 * WINDOW_SIZE - SWEET_SPOT, 2 * WINDOW_SIZE);

    ec::StreamHw& streamHw = *ec::StreamHw::getSingletonStreamHw();
    streamHw.copyToHw(angleTerms, 0, 2 * WINDOW_SIZE, 2 * WINDOW_SIZE);

    for (size_t j = 0; j < numWins; j++)
    {
        memcpy(signalWindow.data(), inputSignal.data() + idxStartWin, WINDOW_SIZE * sizeof(ec::Float));

        // puts cosine terms back into vecHw. we were using that space for buffering
        vecHw.copyToHw(angleTerms, 0, 160, 2 * WINDOW_SIZE); // we *might* be able to recalculate those cosine terms faster?

        streamHw.copyToHw(angleTerms, 0, 512, 2 * WINDOW_SIZE);

        std::vector<ec::Float> inImag(WINDOW_SIZE);

        fft(signalWindow, inImag);

        for (size_t i = 0; i < sizeSpectrum; i++)
        {
            ec::Float freqVal;

            if ((i == 0 || i == 512 && j < numWins - 1))
            {
                ec::Float f = signalWindow[i] * signalWindow[i];
                memcpy(&freqVal, &f, sizeof(ec::Float));
            }
            else
            {
                ec::Float f = signalWindow[i] * signalWindow[i] + inImag[i] * inImag[i];
                memcpy(&freqVal, &f, sizeof(ec::Float));
            }

            // we will always take the first window
            if (j != 0 && freqVal <= preLogSpectrum[i])
            {
                continue;
            }

            memcpy(preLogSpectrum.data() + i, &freqVal, sizeof(ec::Float));

            freqVal = ec_log(freqVal);
            freqVal *= spC0;

            if (i > 0 && i < sizeSpectrum - 1)
            {
                freqVal += spC2;
            }
            else
            {
                freqVal += spC1;
            }

            memcpy(outputSpectrum.data() + i, &freqVal, sizeof(ec::Float));
        }

        idxStartWin += stepBetweenWins;
    }

    return outputSpectrum;
}

void sfft(std::vector<ec::Float>& inputReal, std::vector<ec::Float>& inputImag)
{
    ec::StreamHw& streamHw = *ec::StreamHw::getSingletonStreamHw();
    // streamHw.resetMemTo0();

    streamHw.createFifos(1000);

    std::vector<ec::Float> temp(2 * WINDOW_SIZE);

    memcpy(temp.data(), inputReal.data(), WINDOW_SIZE * sizeof(ec::Float));
    memcpy(temp.data() + WINDOW_SIZE, blackmanCoefs.data(), WINDOW_SIZE * sizeof(ec::Float));

    streamHw.copyToHw(temp, 0, 2 * WINDOW_SIZE, 0);

    // apply blackman coefficients

    size_t c = WINDOW_SIZE / 16;

    for (size_t i = 0; i < 16 * 3; i += 3)
    {
        INIT(streamHw.addOpMulToPipeline(i, i + 1, i + 2));
        streamHw.startStreamDataMemToFifo(c * i / 3, i, c);
        streamHw.startStreamDataMemToFifo(WINDOW_SIZE + c * i / 3, i + 1, c);
        streamHw.startStreamDataFifoToMem(i + 2, c * i / 3, c);
    }


    MEASURE(streamHw.runPipeline(), "Blackman");

    // we now have 16 multiplication pipes
    // pipe (in0,in1,out); 0,1,2; 3,4,5; 6,7,8; 9,10,11; ... 


    size_t p = WINDOW_SIZE / 2;
    size_t d = 2;

    // run the 512 arrays
    for (size_t i = 0; i < 16 * 3; i += 3)
    {
        INIT(streamHw.addOpAddToPipeline(48 + i, 48 + i + 1, 48 + i + 2));
        streamHw.startStreamDataMemToFifo(c * i / (3 * d), 48 + i, c / d);
        streamHw.startStreamDataMemToFifo(p + c * i / (3 * d), 48 + i + 1, c / d);
        streamHw.startStreamDataFifoToMem(48 + i + 2, WINDOW_SIZE + c * i / (3 * d), c / d);
    }

    // created 16 addition pipes at offset 48

    MEASURE(streamHw.runPipeline(), "First 512");

    for (size_t i = 0; i < 16 * 4; i += 4)
    {
        INIT(streamHw.addOpMulToPipeline(96 + i, minusOne, 96 + i + 1));
        INIT(streamHw.addOpAddToPipeline(96 + i + 1, 96 + i + 2, 96 + i + 3));
        streamHw.startStreamDataMemToFifo(p + c * i / (4 * d), 96 + i, c / d);
        streamHw.startStreamDataMemToFifo(c * i / (4 * d), 96 + i + 2, c / d);
        streamHw.startStreamDataFifoToMem(96 + i + 3, p + c * i / (4 * d), c / d);
    }

    MEASURE(streamHw.runPipeline(), "second 512");

    // -1 then add pipes at offset 96

    // move the real vals back to right location
    for (size_t i = 0; i < 32 * 2; i += 2)
    {
        INIT(streamHw.addOpAddToPipeline(160 + i, zero, 160 + i + 1));

        streamHw.startStreamDataMemToFifo(WINDOW_SIZE + c * i / (4 * d), 160 + i, c / (2 * d));
        streamHw.startStreamDataFifoToMem(160 + i + 1, c * i / (4 * d), c / (2 * d));
    }

    MEASURE(streamHw.runPipeline(), "transfer");

    streamHw.resetMemTo0(WINDOW_SIZE, WINDOW_SIZE); // clear out imaginary terms


    for (size_t i = 0; i < 16 * 3; i += 3)
    {
        streamHw.startStreamDataMemToFifo(p + c * i / (3 * d), i, c / d);
        streamHw.startStreamDataMemToFifo(3 * WINDOW_SIZE + c * i / (3 * d), i + 1, c / d);
        streamHw.startStreamDataFifoToMem(i + 2, WINDOW_SIZE + p + c * i / (3 * d), c / d);
    }

    MEASURE(streamHw.runPipeline(), "512 sin terms");

    for (size_t i = 0; i < 16 * 3; i += 3)
    {
        streamHw.startStreamDataMemToFifo(p + c * i / (3 * d), i, c / d);
        streamHw.startStreamDataMemToFifo(2 * WINDOW_SIZE + c * i / 6, i + 1, c / d);
        streamHw.startStreamDataFifoToMem(i + 2, p + c * i / (3 * d), c / d);
    }

    MEASURE(streamHw.runPipeline(), "512 cos terms");

    // now we can use first 512 slots in cos and sin as buffer

    // we need to make the pipelines needed for our complex multiplication starting at offset 224
    for (size_t i = 0; i < 8 * 4; i += 4)
    {
        // real * cos
        INIT(streamHw.addOpMulToPipeline(224 + i, 224 + i + 1, 96 + i + 2));

        // imag * sin * -1
        INIT(streamHw.addOpMulToPipeline(224 + i + 2, 224 + i + 3, 96 + i + 1));
    }
    // real terms to 224 + i
    // cos terms to 224 + i + 1
    // imag terms 224 + i + 2
    // sin terms 224 + i + 3
    // output 96 + i + 3


    for (size_t i = 0; i < 8 * 7; i += 7)
    {
        // real * sin
        INIT(streamHw.addOpMulToPipeline(256 + i, 256 + i + 1, 256 + i + 2));

        // imag * cos
        INIT(streamHw.addOpMulToPipeline(256 + i + 3, 256 + i + 4, 256 + i + 5));

        // add the two
        INIT(streamHw.addOpAddToPipeline(256 + i + 2, 256 + i + 5, 256 + i + 6));
    }
    // real terms to 256 + i
    // sin terms to 256 + i + 1
    // imag terms 256 + i + 3
    // cos terms 256 + i + 4
    // output 256 + i + 6

    while (p > 32 && true)
    {
        p /= 2;
        d *= 2;

        size_t vals = 0;
        size_t valsC = p;
        size_t realC = p;
        size_t imagC = WINDOW_SIZE + p;

        for (size_t j = 0; j < WINDOW_SIZE / (2 * p); j++)
        {
            for (size_t i = 0; i < 16 * 3; i += 3)
            {
                streamHw.startStreamDataMemToFifo(vals + c * i / (3 * d), 48 + i, c / d);
                streamHw.startStreamDataMemToFifo(valsC + c * i / (3 * d), 48 + i + 1, c / d);
                streamHw.startStreamDataFifoToMem(48 + i + 2, 2 * WINDOW_SIZE + c * i / (3 * d), c / d);
            }

            streamHw.runPipeline();

            for (size_t i = 0; i < 16 * 4; i += 4)
            {
                streamHw.startStreamDataMemToFifo(valsC + c * i / (4 * d), 96 + i, c / d);
                streamHw.startStreamDataMemToFifo(vals + c * i / (4 * d), 96 + i + 2, c / d);
                streamHw.startStreamDataFifoToMem(96 + i + 3, valsC + c * i / (4 * d), c / d);
            }

            streamHw.runPipeline();

            for (size_t i = 0; i < 32 * 2; i += 2)
            {
                streamHw.startStreamDataMemToFifo(2 * WINDOW_SIZE + c * i / (4 * d), 160 + i, c / (2 * d));
                streamHw.startStreamDataFifoToMem(160 + i + 1, vals + c * i / (4 * d), c / (2 * d));
            }

            streamHw.runPipeline();

            // offset to imaginary numbers and do it again
            vals += WINDOW_SIZE;
            valsC += WINDOW_SIZE;

            if (j != 0) // manually exclude times where imaginary numbers aren't populated
            {
                // imaginaries
                for (size_t i = 0; i < 16 * 3; i += 3)
                {
                    streamHw.startStreamDataMemToFifo(vals + c * i / (3 * d), 48 + i, c / d);
                    streamHw.startStreamDataMemToFifo(valsC + c * i / (3 * d), 48 + i + 1, c / d);
                    streamHw.startStreamDataFifoToMem(48 + i + 2, 2 * WINDOW_SIZE + c * i / (3 * d), c / d);
                }

                streamHw.runPipeline();

                for (size_t i = 0; i < 16 * 4; i += 4)
                {
                    streamHw.startStreamDataMemToFifo(valsC + c * i / (4 * d), 96 + i, c / d);
                    streamHw.startStreamDataMemToFifo(vals + c * i / (4 * d), 96 + i + 2, c / d);
                    streamHw.startStreamDataFifoToMem(96 + i + 3, valsC + c * i / (4 * d), c / d);
                }

                streamHw.runPipeline();

                for (size_t i = 0; i < 32 * 2; i += 2)
                {
                    streamHw.startStreamDataMemToFifo(2 * WINDOW_SIZE + c * i / (4 * d), 160 + i, c / (2 * d));
                    streamHw.startStreamDataFifoToMem(160 + i + 1, vals + c * i / (4 * d), c / (2 * d));
                }

                streamHw.runPipeline();
            }

            // remove the imaginary val offset and move on to the next pair
            vals -= WINDOW_SIZE;
            vals += 2 * p;

            valsC -= WINDOW_SIZE;
            valsC += 2 * p;

            if (j == 0)
            {
                for (size_t i = 0; i < 16 * 3; i += 3)
                {
                    streamHw.startStreamDataMemToFifo(realC + c * i / (3 * d), i, c / d);
                    streamHw.startStreamDataMemToFifo(3 * WINDOW_SIZE + c * i / (3 * d), i + 1, c / d);
                    streamHw.startStreamDataFifoToMem(i + 2, imagC + c * i / (3 * d), c / d);
                }

                streamHw.runPipeline();

                for (size_t i = 0; i < 16 * 3; i += 3)
                {
                    streamHw.startStreamDataMemToFifo(realC + c * i / (3 * d), i, c / d);
                    streamHw.startStreamDataMemToFifo(2 * WINDOW_SIZE + c * i / (3 * d), i + 1, c / d);
                    streamHw.startStreamDataFifoToMem(i + 2, realC + c * i / (3 * d), c / d);
                }

                streamHw.runPipeline();
            }
            else
            {



            }

            // complex multiplication stage:
            vecHw.mul32(realC, 2 * WINDOW_SIZE + 32 * i + (WINDOW_SIZE - 2 * c), buf0); // real [C - 2C] * cos(C-2C) -> 2 * WINDOW_SIZE ( we're using this as a buffer cause we're done using it in this fft rn)
            vecHw.mul32(realC, 3 * WINDOW_SIZE + 32 * i + (WINDOW_SIZE - 2 * c), buf2); // imag [C - 2C] = real[C - 2C] * sin(C-2C)

            vecHw.mul32(imagC, 3 * WINDOW_SIZE + 32 * i + (WINDOW_SIZE - 2 * c), buf1); // imag [C - 2C] *  sin(C-2C) -> 3 * WINDOW_SIZE     
            vecHw.mul32(imagC, 2 * WINDOW_SIZE + 32 * i + (WINDOW_SIZE - 2 * c), imagC); // imag [C - 2C] = imag [C - 2C] * cos(C-2C)

            vecHw.mul32(buf1, minusOne, buf1, 32); // (-1)* (imag [C - 2C] *  sin(C-2C))
            vecHw.add32(buf0, buf1, realC); // real [C - 2C] = real [C - 2C] * cos(C-2C)  + (-1)* (imag [C - 2C] *  sin(C-2C))
            vecHw.add32(imagC, buf2, imagC); // imag [C - 2C] = (imag [C - 2C] * cos(C-2C)) + (real[C - 2C] * sin(C-2C))

            // repeat for 3C, 5C, and 7C
            realC += 2 * c;
            imagC += 2 * c;
        }
    }




    streamHw.copyFromHw(temp, 0, 2 * WINDOW_SIZE, 0);

    memcpy(inputReal.data(), temp.data(), WINDOW_SIZE * sizeof(ec::Float));
    memcpy(inputImag.data(), temp.data() + WINDOW_SIZE, WINDOW_SIZE * sizeof(ec::Float));


    streamInitted = true;
}

void fft(std::vector<ec::Float>& inputReal, std::vector<ec::Float>& inputImag)
{
    ec::VecHw& vecHw = *ec::VecHw::getSingletonVecHw();

    sfft(inputReal, inputImag);

    std::vector<ec::Float> temp(2 * WINDOW_SIZE);

    memcpy(temp.data(), inputReal.data(), WINDOW_SIZE * sizeof(ec::Float));
    memcpy(temp.data() + WINDOW_SIZE, inputImag.data(), WINDOW_SIZE * sizeof(ec::Float));

    vecHw.copyToHw(temp, 0, 2 * WINDOW_SIZE, 0);

    // vecHw.resetMemTo0(2 * offset, WINDOW_SIZE - offset);

    // apply the blackman coefficients to input data
    // these are stored the location of the imaginary input data because we don't have any
    // for (size_t i = 0; i < WINDOW_SIZE / 32; i++)
    // {
    //     vecHw.mul32(32 * i, offset + 32 * i, 32 * i);
    // }

    int c = WINDOW_SIZE / 2;
    {
        // c = 512 only working with reals for now
        // for (size_t i = 0; i < c / 32; i++)
        // {
        //     // (-1)*inputReal[C - 2C]
        //     vecHw.mul32(REAL(c), minusOne, IMAG(c), 32); // stick this in imag[C - 2C] for now. it will be fixed later

        //     // [0 - C] = [0 - C] + [C - 2C]
        //     vecHw.add32(REAL(0), REAL(c), IMAG(0)); // store in imag[0] for now.

        //     // [C - 2C] = [0 - C] + (-1)*[C - 2C]
        //     vecHw.add32(REAL(0), IMAG(c), REAL(c));

        //     // now we can assign it to [0 - C] because we're done reading it
        //     vecHw.assign32(IMAG(0), REAL(0));
        // }

        // vecHw.resetMemTo0(WINDOW_SIZE, WINDOW_SIZE); // reset IMAG(0)




        /*
        input[C - 2C] * omega[0 through C] # this is complex multiplication
        example of complex multiplication:
        (x + yj) * (a + bj) = (xa - yb) + (xb + ya)j
        So the real part of the product is (xa - yb), and the imaginary part of the product is (xb + ya).
       */
       // for (size_t i = 0; i < c / 32; i++)
       // {
       //     vecHw.mul32(REAL(c), 3 * WINDOW_SIZE + 32 * i, IMAG(c)); // imag [C - 2C] = real[C - 2C] * sin(0 - C)        
       //     vecHw.mul32(REAL(c), 2 * WINDOW_SIZE + 32 * i, REAL(c)); // real [C - 2C] = real [C - 2C] * cos(0 - C)
       // }

    }

    // Note we have now used the first 512 values in 2 * WINDOW_SIZE and 3 * WINDOW_SIZE and will no longer need them.
    size_t buf0 = 2 * WINDOW_SIZE;
    size_t buf1 = 2 * WINDOW_SIZE + 32;
    size_t buf2 = 2 * WINDOW_SIZE + 64;
    size_t buf3 = 2 * WINDOW_SIZE + 96;
    size_t buf4 = 2 * WINDOW_SIZE + 128;

    while (c > 32 && false) // sweet spot
    {
        c /= 2;

        // when the arrays we are working with are larger than 32, our operations need to be chunked out into chunks of 32
        // the variable i tracks which chunk we are currently in
        for (size_t i = 0; i < std::max(1, c / 32); i++)
        {
            size_t vals = REAL(0);
            size_t valsC = REAL(c);

            size_t realC = REAL(c);
            size_t imagC = IMAG(c);

            // we need to run for real and imaginary
            // [0 - C] = [0 - C] + [C - 2C]
            // [C - 2C] = [0 - C] + (-1)*[C - 2C]
            // [2C - 3C] = [2C- 3C] + [3C - 4C]
            // [3C - 4C] = [2C - 3C] + (-1)*[3C - 4C]

            // [4C - 5C] = [4C - 5C] + [5C - 6C]
            // [5C - 6C] = [4C - 5C] + (-1)*[5C - 6C]
            // [6C - 7C] = [6C - 7C] + [7C - 8C]
            // [7C - 8C] = [6C - 7C] + (-1)*[7C - 8C]

            // this is the pair loop
            // it will repeat for each pair of matricies that need to be processed
            // here c = 128, so 4 pairs of cross addition and complex multiplication needs to happen
            /*
            | -- cross addition -- | | complex multiplication |

            j = 0

            j * c   ----\--/---- : vals
                         \/
                         /\
            (j+1)*c ----/--\---- : valsC -> realC * omega[WINDOW_SIZE - 2 * c]

            j++

            j * c   ----\--/---- : vals
                         \/
                         /\
            (j+1)*c ----/--\---- : valsC -> realC * omega[WINDOW_SIZE - 2 * c]

            j++

            j * c   ----\--/---- : vals
                         \/
                         /\
            (j+1)*c ----/--\---- : valsC -> realC * omega[WINDOW_SIZE - 2 * c]

            j++

            j * c   ----\--/---- : vals
                         \/
                         /\
            (j+1)*c ----/--\---- : valsC -> realC * omega[WINDOW_SIZE - 2 * c]
            */
            for (size_t j = 0; j < WINDOW_SIZE / (2 * c); j++)
            {
                // addition stage:
                // (-1)*[3C - 4C]
                vecHw.mul32(valsC, minusOne, buf0, 32);

                // [3C - 4C] = [2C - 3C] + (-1)*[3C - 4C]
                vecHw.add32(vals, buf0, buf1);

                // [2C - 3C] = [2C - 3C] + [3C - 4C]
                vecHw.add32(vals, valsC, vals);

                // [3C - 4C] = [2C - 3C] + (-1)*[3C - 4C]
                vecHw.assign32(buf1, valsC);

                // offset to imaginary numbers and do it again
                vals += WINDOW_SIZE;
                valsC += WINDOW_SIZE;

                if (j != 0) // manually exclude times where imaginary numbers aren't populated
                {
                    // imaginaries
                    // (-1)*[3C - 4C]
                    vecHw.mul32(valsC, minusOne, buf0, 32);

                    // [3C - 4C] = [2C - 3C] + (-1)*[3C - 4C]
                    vecHw.add32(vals, buf0, buf1);

                    // [2C - 3C] = [2C - 3C] + [3C - 4C]
                    vecHw.add32(vals, valsC, vals);

                    // [3C - 4C] = [2C - 3C] + (-1)*[3C - 4C]
                    vecHw.assign32(buf1, valsC);
                }

                // remove the imaginary val offset and move on to the next pair
                vals -= WINDOW_SIZE;
                vals += 2 * c;

                valsC -= WINDOW_SIZE;
                valsC += 2 * c;

                // complex multiplication stage:
                vecHw.mul32(realC, 2 * WINDOW_SIZE + 32 * i + (WINDOW_SIZE - 2 * c), buf0); // real [C - 2C] * cos(C-2C) -> 2 * WINDOW_SIZE ( we're using this as a buffer cause we're done using it in this fft rn)
                vecHw.mul32(realC, 3 * WINDOW_SIZE + 32 * i + (WINDOW_SIZE - 2 * c), buf2); // imag [C - 2C] = real[C - 2C] * sin(C-2C)

                if (j != 0)
                {
                    vecHw.mul32(imagC, 3 * WINDOW_SIZE + 32 * i + (WINDOW_SIZE - 2 * c), buf1); // imag [C - 2C] *  sin(C-2C) -> 3 * WINDOW_SIZE     
                    vecHw.mul32(imagC, 2 * WINDOW_SIZE + 32 * i + (WINDOW_SIZE - 2 * c), imagC); // imag [C - 2C] = imag [C - 2C] * cos(C-2C)

                    vecHw.mul32(buf1, minusOne, buf1, 32); // (-1)* (imag [C - 2C] *  sin(C-2C))
                    vecHw.add32(buf0, buf1, realC); // real [C - 2C] = real [C - 2C] * cos(C-2C)  + (-1)* (imag [C - 2C] *  sin(C-2C))
                    vecHw.add32(imagC, buf2, imagC); // imag [C - 2C] = (imag [C - 2C] * cos(C-2C)) + (real[C - 2C] * sin(C-2C))
                }
                else
                {
                    vecHw.assign32(buf0, realC);
                    vecHw.assign32(buf2, imagC);
                }

                // repeat for 3C, 5C, and 7C
                realC += 2 * c;
                imagC += 2 * c;
            }
        }
    }

    for (size_t i = 0; i < WINDOW_SIZE; i++)
    {
        std::cout << i << " " << vecHw.m_mem[i] << std::endl;
    }
    return;
    // c = 32;

    // rearrage step
    size_t vals = 0;
    size_t valsC = 32;

    size_t realC = 32;
    size_t imagC = WINDOW_SIZE + 32;

    c /= 2;

    //todo do this while generating the angle terms
    vecHw.assign32(2 * WINDOW_SIZE + (WINDOW_SIZE - 2 * c), buf3, c);
    vecHw.assign32(2 * WINDOW_SIZE + (WINDOW_SIZE - 2 * c), buf3 + c, c);
    vecHw.assign32(3 * WINDOW_SIZE + (WINDOW_SIZE - 2 * c), buf4, c);
    vecHw.assign32(3 * WINDOW_SIZE + (WINDOW_SIZE - 2 * c), buf4 + c, c);

    // for every other pair of 32 sized arrays
    // we need to rearrange their layout            
    // todo see if i can readd the j != 0 stuff again
    for (size_t j = 0; j < WINDOW_SIZE / (2 * 32); j++)
    {
        // take the last 16 values from the first array and swap them with the first 16 values from the second array
        vecHw.assign32(vals + c, buf0, c);
        vecHw.assign32(valsC, vals + c, c);
        vecHw.assign32(buf0, valsC, c);

        // now we can run the next set of addition but in 32 value chunks
        vecHw.mul32(valsC, minusOne, buf0, 32); // we need to specify the 32 here because of a reassignment oversight in measurement

        // [3C - 4C] = [2C - 3C] + (-1)*[3C - 4C]
        vecHw.add32(vals, buf0, buf1);

        // [2C - 3C] = [2C - 3C] + [3C - 4C]
        vecHw.add32(vals, valsC, vals);

        // [3C - 4C] = [2C - 3C] + (-1)*[3C - 4C]
        vecHw.assign32(buf1, valsC);

        // offset to imaginary numbers and do it again
        vals += WINDOW_SIZE;
        valsC += WINDOW_SIZE;

        // take the last 16 values from the first array and swap them with the first 16 values from the second array
        vecHw.assign32(vals + c, buf0, c);
        vecHw.assign32(valsC, vals + c, c);
        vecHw.assign32(buf0, valsC, c);

        // imaginaries
        // (-1)*[3C - 4C]
        vecHw.mul32(valsC, minusOne, buf0, 32);

        // [3C - 4C] = [2C - 3C] + (-1)*[3C - 4C]
        vecHw.add32(vals, buf0, buf1);

        // [2C - 3C] = [2C - 3C] + [3C - 4C]
        vecHw.add32(vals, valsC, vals);

        // [3C - 4C] = [2C - 3C] + (-1)*[3C - 4C]
        vecHw.assign32(buf1, valsC);

        // complex multiplication stage:
        vecHw.mul32(realC, buf3, buf0); // real [C - 2C] * cos(C-2C) -> 2 * WINDOW_SIZE ( we're using this as a buffer cause we're done using it in this fft rn)

        vecHw.mul32(imagC, buf4, buf1); // imag [C - 2C] *  sin(C-2C) -> 3 * WINDOW_SIZE     
        vecHw.mul32(imagC, buf3, imagC); // imag [C - 2C] = imag [C - 2C] * cos(C-2C)

        vecHw.mul32(realC, buf4, buf2); // imag [C - 2C] = real[C - 2C] * sin(C-2C)

        vecHw.mul32(buf1, minusOne, buf1, 32); // (-1)* (imag [C - 2C] *  sin(C-2C))
        vecHw.add32(buf0, buf1, realC); // real [C - 2C] = real [C - 2C] * cos(C-2C)  + (-1)* (imag [C - 2C] *  sin(C-2C))
        vecHw.add32(imagC, buf2, imagC); // imag [C - 2C] = (imag [C - 2C] * cos(C-2C)) + (real[C - 2C] * sin(C-2C))

        // repeat for 3C, 5C, and 7C
        realC += 64;
        imagC += 64;

        // take the last 16 values from the first array and swap them with the first 16 values from the second array
        vecHw.assign32(vals + c, buf0, c);
        vecHw.assign32(valsC, vals + c, c);
        vecHw.assign32(buf0, valsC, c);

        // remove the imaginary val offset and move on to the next pair
        vals -= WINDOW_SIZE;
        valsC -= WINDOW_SIZE;

        // take the last 16 values from the first array and swap them with the first 16 values from the second array
        vecHw.assign32(vals + c, buf0, c);
        vecHw.assign32(valsC, vals + c, c);
        vecHw.assign32(buf0, valsC, c);

        vals += 64;
        valsC += 64;
    }

    c = 8;

    vals = 0;
    valsC = 32;

    realC = 32;
    imagC = WINDOW_SIZE + 32;

    for (size_t i = 0; i < 32 / c; i++)
    {
        vecHw.assign32(2 * WINDOW_SIZE + (WINDOW_SIZE - 2 * c), buf3 + i * c, c);
        vecHw.assign32(3 * WINDOW_SIZE + (WINDOW_SIZE - 2 * c), buf4 + i * c, c);
    }

    for (size_t j = 0; j < WINDOW_SIZE / (2 * 32); j++)
    {
        vecHw.assign32(vals + 16, buf0, 2 * c);
        vecHw.assign32(vals + 32, vals + 16, 2 * c);
        vecHw.assign32(buf0, vals + 32, 2 * c);

        vecHw.assign32(vals + 8, buf0, c);
        vecHw.assign32(vals + 32, vals + 8, c);
        vecHw.assign32(buf0, vals + 32, c);

        vecHw.assign32(vals + 24, buf0, c);
        vecHw.assign32(vals + 48, vals + 24, c);
        vecHw.assign32(buf0, vals + 48, c);

        // now we can run the next set of addition but in 32 value chunks
        vecHw.mul32(valsC, minusOne, buf0, 32);

        // [3C - 4C] = [2C - 3C] + (-1)*[3C - 4C]
        vecHw.add32(vals, buf0, buf1);

        // [2C - 3C] = [2C - 3C] + [3C - 4C]
        vecHw.add32(vals, valsC, vals);

        // [3C - 4C] = [2C - 3C] + (-1)*[3C - 4C]
        vecHw.assign32(buf1, valsC);

        // offset to imaginary numbers and do it again
        vals += WINDOW_SIZE;
        valsC += WINDOW_SIZE;

        vecHw.assign32(vals + 16, buf0, 2 * c);
        vecHw.assign32(vals + 32, vals + 16, 2 * c);
        vecHw.assign32(buf0, vals + 32, 2 * c);

        vecHw.assign32(vals + 8, buf0, c);
        vecHw.assign32(vals + 32, vals + 8, c);
        vecHw.assign32(buf0, vals + 32, c);

        vecHw.assign32(vals + 24, buf0, c);
        vecHw.assign32(vals + 48, vals + 24, c);
        vecHw.assign32(buf0, vals + 48, c);

        // imaginaries
        // (-1)*[3C - 4C]
        vecHw.mul32(valsC, minusOne, buf0, 32);

        // [3C - 4C] = [2C - 3C] + (-1)*[3C - 4C]
        vecHw.add32(vals, buf0, buf1);

        // [2C - 3C] = [2C - 3C] + [3C - 4C]
        vecHw.add32(vals, valsC, vals);

        // [3C - 4C] = [2C - 3C] + (-1)*[3C - 4C]
        vecHw.assign32(buf1, valsC);

        // complex multiplication stage:
        vecHw.mul32(realC, buf3, buf0); // real [C - 2C] * cos(C-2C) -> 2 * WINDOW_SIZE ( we're using this as a buffer cause we're done using it in this fft rn)

        vecHw.mul32(imagC, buf4, buf1); // imag [C - 2C] *  sin(C-2C) -> 3 * WINDOW_SIZE     
        vecHw.mul32(imagC, buf3, imagC); // imag [C - 2C] = imag [C - 2C] * cos(C-2C)

        vecHw.mul32(realC, buf4, buf2); // imag [C - 2C] = real[C - 2C] * sin(C-2C)

        vecHw.mul32(buf1, minusOne, buf1, 32); // (-1)* (imag [C - 2C] *  sin(C-2C))
        vecHw.add32(buf0, buf1, realC); // real [C - 2C] = real [C - 2C] * cos(C-2C)  + (-1)* (imag [C - 2C] *  sin(C-2C))
        vecHw.add32(imagC, buf2, imagC); // imag [C - 2C] = (imag [C - 2C] * cos(C-2C)) + (real[C - 2C] * sin(C-2C))

        // repeat for 3C, 5C, and 7C
        realC += 64;
        imagC += 64;

        vecHw.assign32(vals + 24, buf0, c);
        vecHw.assign32(vals + 48, vals + 24, c);
        vecHw.assign32(buf0, vals + 48, c);

        vecHw.assign32(vals + 8, buf0, c);
        vecHw.assign32(vals + 32, vals + 8, c);
        vecHw.assign32(buf0, vals + 32, c);

        vecHw.assign32(vals + 16, buf0, 2 * c);
        vecHw.assign32(vals + 32, vals + 16, 2 * c);
        vecHw.assign32(buf0, vals + 32, 2 * c);

        // remove the imaginary val offset and move on to the next pair
        vals -= WINDOW_SIZE;
        valsC -= WINDOW_SIZE;

        vecHw.assign32(vals + 24, buf0, c);
        vecHw.assign32(vals + 48, vals + 24, c);
        vecHw.assign32(buf0, vals + 48, c);

        vecHw.assign32(vals + 8, buf0, c);
        vecHw.assign32(vals + 32, vals + 8, c);
        vecHw.assign32(buf0, vals + 32, c);

        vecHw.assign32(vals + 16, buf0, 2 * c);
        vecHw.assign32(vals + 32, vals + 16, 2 * c);
        vecHw.assign32(buf0, vals + 32, 2 * c);

        vals += 64;
        valsC += 64;
    }

    vecHw.copyFromHw(temp, 0, 2 * WINDOW_SIZE, 0);

    memcpy(inputReal.data(), temp.data(), WINDOW_SIZE * sizeof(ec::Float));
    memcpy(inputImag.data(), temp.data() + WINDOW_SIZE, WINDOW_SIZE * sizeof(ec::Float));

    std::vector<ec::Float> buffer(32);

#define X0R buffer[0]
#define X1R buffer[1]
#define X2R buffer[2]
#define X3R buffer[3]
#define X4R buffer[4]
#define X5R buffer[5]
#define X6R buffer[6]
#define X7R buffer[7]
#define X0I buffer[8]
#define X1I buffer[9]
#define X2I buffer[10]
#define X3I buffer[11]
#define X4I buffer[12]
#define X5I buffer[13]
#define X6I buffer[14]
#define X7I buffer[15]

    // Radix-8 fft
    for (size_t i = 0; i < 1024; i += 8)
    {
        memcpy(buffer.data(), inputReal.data() + i, 8 * sizeof(ec::Float));
        memcpy(buffer.data() + 8, inputImag.data() + i, 8 * sizeof(ec::Float));

        memcpy(inputReal.data() + i + 1, &X0R, sizeof(ec::Float));
        memcpy(inputReal.data() + i + 2, &X0R, sizeof(ec::Float));
        memcpy(inputReal.data() + i + 4, &X0R, sizeof(ec::Float));
        memcpy(inputReal.data() + i + 5, &X0R, sizeof(ec::Float));
        memcpy(inputReal.data() + i + 6, &X0R, sizeof(ec::Float));

        memcpy(inputImag.data() + i + 1, &X0I, sizeof(ec::Float));
        memcpy(inputImag.data() + i + 2, &X0I, sizeof(ec::Float));
        memcpy(inputImag.data() + i + 4, &X0I, sizeof(ec::Float));
        memcpy(inputImag.data() + i + 5, &X0I, sizeof(ec::Float));
        memcpy(inputImag.data() + i + 6, &X0I, sizeof(ec::Float));

        if (i == 0)
        {
            ec::Float c1 = X6R - X2R;
            ec::Float c2 = (X1R - X3R - X5R + X7R) * sqrt22;
            ec::Float c3 = (X5R - X3R - X1R + X7R) * sqrt22;
            ec::Float c4 = X3R - X1R - X5R + X7R;

            inputReal[i + 0] += X1R + X2R + X3R + X4R + X5R + X6R + X7R;

            inputReal[i + 1] += c2 - X4R;
            inputImag[i + 1] += c3 + c1;

            inputReal[i + 2] += X4R - X2R - X6R;
            inputImag[i + 2] += c4;

            inputReal[i + 4] += X2R - X1R - X3R + X4R - X5R + X6R - X7R;

            inputReal[i + 5] -= c2 + X4R;
            inputImag[i + 5] -= c3 - c1;

            // inputReal[i + 6] += X4R - X2R - X6R;
            memcpy(inputReal.data() + i + 6, inputReal.data() + i + 2, sizeof(ec::Float));
            inputImag[i + 6] -= c4;
        }
        else
        {
            ec::Float c1 = X2R - X4R + X6R;
            ec::Float c2 = X4I - X2I - X6I;

            ec::Float c3 = X1I - X3I;

            inputReal[i + 0] += X1R + X2R + X3R + X4R + X5R + X6R + X7R;
            inputImag[i + 0] += X1I + X2I + X3I + X4I + X5I + X6I + X7I;

            inputReal[i + 2] += c3 - c1 + (X5I - X7I);
            inputImag[i + 2] += c2 - X1R + X3R - X5R + X7R;

            inputReal[i + 4] += X2R - X1R - X3R + X4R - X5R + X6R - X7R;
            inputImag[i + 4] += X2I - X1I - X3I + X4I - X5I + X6I - X7I;

            inputReal[i + 6] += (X7I - X5I) - (c3 + c1);
            inputImag[i + 6] += c2 + X1R - X3R + X5R - X7R;
        }
    }

    ec::Float swapVal;

    for (size_t i = 0; i < WINDOW_SIZE; i++)
    {
        uint16_t newI = reverse_bits(i);

        if (i < newI)
        {
            memcpy(&swapVal, inputReal.data() + i, sizeof(ec::Float));
            memcpy(inputReal.data() + i, inputReal.data() + newI, sizeof(ec::Float));
            memcpy(inputReal.data() + newI, &swapVal, sizeof(ec::Float));

            memcpy(&swapVal, inputImag.data() + i, sizeof(ec::Float));
            memcpy(inputImag.data() + i, inputImag.data() + newI, sizeof(ec::Float));
            memcpy(inputImag.data() + newI, &swapVal, sizeof(ec::Float));
        }
    }

    for (size_t i = 128; i < 256; i++)
    {
        uint16_t newI = 512 - (i - 128);

        memcpy(&swapVal, inputReal.data() + i, sizeof(ec::Float));
        memcpy(inputReal.data() + i, inputReal.data() + newI, sizeof(ec::Float));
        memcpy(inputReal.data() + newI, &swapVal, sizeof(ec::Float));

        memcpy(&swapVal, inputImag.data() + i, sizeof(ec::Float));
        memcpy(inputImag.data() + i, inputImag.data() + newI, sizeof(ec::Float));
        memcpy(inputImag.data() + newI, &swapVal, sizeof(ec::Float));
    }

    memcpy(inputReal.data() + 384, inputReal.data() + 640, sizeof(ec::Float));
    memcpy(inputImag.data() + 384, inputImag.data() + 640, sizeof(ec::Float));
}

/*
    The below functions have zero or absolute minimal cost.
*/
void init_angleTerms()
{
    float aC(float(-2.0f * M_PI) / WINDOW_SIZE);

    for (size_t i = 0; i < WINDOW_SIZE / 2; i++)
    {
        float co = std::cos(aC * i);
        angleTerms[i] = std::cos(aC * i);

        float si = std::sin(aC * i);
        angleTerms[i + WINDOW_SIZE] = std::sin(aC * i);
    }

    size_t x = WINDOW_SIZE / 4;
    size_t y = 2;
    while (x > 0)
    {
        for (size_t i = 0; i < y * x; i += y)
        {
            memcpy(angleTerms.data() + WINDOW_SIZE - (2 * x) + i / y, angleTerms.data() + i, sizeof(ec::Float));
            memcpy(angleTerms.data() + WINDOW_SIZE + WINDOW_SIZE - (2 * x) + i / y, angleTerms.data() + WINDOW_SIZE + i, sizeof(ec::Float));
        }

        x /= 2;
        y *= 2;
    }
}

void init_blackmanCoefs()
{
    for (size_t i = 0; i < WINDOW_SIZE; i++)
    {
        float coef = 0.42f - 0.5f * std::cos(i * 2.0f * M_PI / (WINDOW_SIZE - 1)) + 0.08f * std::cos(i * 4.0f * M_PI / (WINDOW_SIZE - 1));
        blackmanCoefs[i] = coef;
    }
}

uint16_t reverse_bits(uint16_t x)
{
    x = ((x & 0xAAAA) >> 1) | ((x & 0x5555) << 1);
    x = ((x & 0xCCCC) >> 2) | ((x & 0x3333) << 2);
    x = ((x & 0xF0F0) >> 4) | ((x & 0x0F0F) << 4);
    x = ((x & 0xFF00) >> 8) | ((x & 0x00FF) << 8);
    return x >> 6; // make it 10 bit
}