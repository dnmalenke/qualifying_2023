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

static constexpr float OVERLAP_RATIO = 0.75;
static constexpr size_t WINDOW_SIZE = 1024;
static const ec::Float PI = 3.14159265358979323846f;
static const ec::Float minusOne = -1.0f;
static const ec::Float minusTwo = -2.0f;

static const ec::Float spC0((float)(10.0f / log(10.0f)));
static const ec::Float spC1((float)(10.0f * log(125.0f / 131072.0f) / log(10.0f)));
static const ec::Float spC2((float)(10.0f * log(125.0f / 32768.0f) / log(10.0f)));

void init_blackmanCoefs(std::vector<ec::Float>& input);
void init_angleTerms();

void fft(std::vector<ec::Float>& inputReal, std::vector<ec::Float>& inputImag, size_t count);

static std::vector<ec::Float> angleTerms(2 * WINDOW_SIZE);

std::vector<ec::Float> process_signal(const std::vector<ec::Float>& inputSignal)
{
    const size_t numSamples = inputSignal.size(); // assume divisible by 32
    const size_t sizeSpectrum = (WINDOW_SIZE / 2) + 1;
    const size_t stepBetweenWins = static_cast<size_t>(ceil(WINDOW_SIZE * (1 - OVERLAP_RATIO)));
    const size_t numWins = (numSamples - WINDOW_SIZE) / stepBetweenWins + 1;

    std::vector<ec::Float> signalWindow(WINDOW_SIZE);
    std::vector<ec::Float> blackmanCoefs(WINDOW_SIZE);

    std::vector<ec::Float> outputSpectrum(sizeSpectrum, std::numeric_limits<float>::lowest());
    std::vector<ec::Float> preLogSpectrum(sizeSpectrum, std::numeric_limits<float>::lowest());

    size_t idxStartWin = 0;

    init_blackmanCoefs(blackmanCoefs);
    init_angleTerms();


    for (size_t j = 0; j < numWins; j++)
    {
        for (size_t i = 0; i < WINDOW_SIZE; i++)
        {
            signalWindow[i] = inputSignal[i + idxStartWin] * blackmanCoefs[i];
        }

        std::vector<ec::Float> inImag(WINDOW_SIZE);

        ec::VecHw& vecHw = *ec::VecHw::getSingletonVecHw();
        vecHw.resetMemTo0();

        vecHw.copyToHw(angleTerms, 0, 2 * WINDOW_SIZE, 2 * WINDOW_SIZE);

        fft(signalWindow, inImag, WINDOW_SIZE);

        for (size_t i = 0; i < sizeSpectrum; i++)
        {
            ec::Float freqVal = signalWindow[i] * signalWindow[i] + inImag[i] * inImag[i];

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

#define REAL(x) x + 32 * i
#define IMAG(x) WINDOW_SIZE + x + 32 * i

uint16_t reverse_bits(uint16_t x) {
    x = ((x & 0xAAAA) >> 1) | ((x & 0x5555) << 1);
    x = ((x & 0xCCCC) >> 2) | ((x & 0x3333) << 2);
    x = ((x & 0xF0F0) >> 4) | ((x & 0x0F0F) << 4);
    x = ((x & 0xFF00) >> 8) | ((x & 0x00FF) << 8);
    return x >> 6;
}

void fft(std::vector<ec::Float>& inputReal, std::vector<ec::Float>& inputImag, size_t count)
{
    ec::VecHw& vecHw = *ec::VecHw::getSingletonVecHw();

    vecHw.copyToHw(inputReal, 0, WINDOW_SIZE, 0);

    int c = WINDOW_SIZE / 2;

    // c = 512 only working with reals for now
    for (size_t i = 0; i < c / 32; i++)
    {
        // (-1)*inputReal[C - 2C]
        vecHw.mul32(REAL(c), minusOne, IMAG(c)); // stick this in imag[C - 2C] for now. it will be fixed later

        // [0 - C] = [0 - C] + [C - 2C]
        vecHw.add32(REAL(0), REAL(c), IMAG(0)); // store in imag[0] for now.

        // [C - 2C] = [0 - C] + (-1)*[C - 2C]
        vecHw.add32(REAL(0), IMAG(c), REAL(c));

        // now we can assign it to [0 - C] because we're done reading it
        vecHw.assign32(IMAG(0), REAL(0));

        vecHw.mul32(IMAG(0), ec::Float(0.0f), IMAG(0)); // just in case. probably can remove later
    }

    /*
    input[C - 2C] * omega[0 through C] # this is complex multiplication
    example of complex multiplication:
    (x + yj) * (a + bj) = (xa - yb) + (xb + ya)j
    So the real part of the product is (xa - yb), and the imaginary part of the product is (xb + ya).
   */
    for (size_t i = 0; i < c / 32; i++)
    {
        vecHw.mul32(REAL(c), 3 * WINDOW_SIZE + 32 * i, IMAG(c)); // imag [C - 2C] = real[C - 2C] * sin(0 - C)
        vecHw.mul32(REAL(c), 2 * WINDOW_SIZE + 32 * i, REAL(c)); // real [C - 2C] = real [C - 2C] * cos(0 - C)
    }

    // Note we have now used the first 512 values in 2 * WINDOW_SIZE and 3 * WINDOW_SIZE and will no longer need them.

    while (c > 1)
    {
        c /= 2;

        for (size_t i = 0; i < std::max(1, c / 32); i++)
        {
            size_t vals = REAL(0);
            size_t valsC = REAL(c);
            size_t buf0 = 2 * WINDOW_SIZE + 32 * i;
            size_t buf1 = 2 * WINDOW_SIZE + c + 32 * i;

            // we need to run for real and imaginary
            // [0 - C] = [0 - C] + [C - 2C]
            // [C - 2C] = [0 - C] + (-1)*[C - 2C]
            // [2C - 3C] = [2C- 3C] + [3C - 4C]
            // [3C - 4C] = [2C - 3C] + (-1)*[3C - 4C]

            // [4C - 5C] = [4C - 5C] + [5C - 6C]
            // [5C - 6C] = [4C - 5C] + (-1)*[5C - 6C]
            // [6C - 7C] = [6C - 7C] + [7C - 8C]
            // [7C - 8C] = [6C - 7C] + (-1)*[7C - 8C]
            for (size_t j = 0; j < WINDOW_SIZE / (2 * c); j++)
            {
                // (-1)*[3C - 4C]
                vecHw.mul32(valsC, minusOne, buf0, std::min(c, 32));

                // [3C - 4C] = [2C - 3C] + (-1)*[3C - 4C]
                vecHw.add32(vals, buf0, buf1, std::min(c, 32));

                // [2C - 3C] = [2C - 3C] + [3C - 4C]
                vecHw.add32(vals, valsC, vals, std::min(c, 32));

                // [3C - 4C] = [2C - 3C] + (-1)*[3C - 4C]
                vecHw.assign32(buf1, valsC, std::min(c, 32));

                vals += WINDOW_SIZE; // imag
                valsC += WINDOW_SIZE; // imag

                // imaginaries
                // (-1)*[3C - 4C]
                vecHw.mul32(valsC, minusOne, buf0, std::min(c, 32));

                // [3C - 4C] = [2C - 3C] + (-1)*[3C - 4C]
                vecHw.add32(vals, buf0, buf1, std::min(c, 32));

                // [2C - 3C] = [2C - 3C] + [3C - 4C]
                vecHw.add32(vals, valsC, vals, std::min(c, 32));

                // [3C - 4C] = [2C - 3C] + (-1)*[3C - 4C]
                vecHw.assign32(buf1, valsC, std::min(c, 32));

                vals -= WINDOW_SIZE;
                vals += 2 * c;

                valsC -= WINDOW_SIZE;
                valsC += 2 * c;
            }
        }

        // input[C - 2C] * omega[0 through C] 
        // input[3C - 4C] * omega[0 through C] 
        // input[5C - 6C] * omega[0 through C] 
        // input[7C - 8C] * omega[0 through C] 
        for (size_t i = 0; i < std::max(1, c / 32); i++)
        {
            size_t realC = REAL(c);
            size_t imagC = IMAG(c);
            size_t buf0 = 2 * WINDOW_SIZE + 32 * i;
            size_t buf1 = 2 * WINDOW_SIZE + c + 32 * i;
            size_t buf2 = 3 * WINDOW_SIZE + 32 * i;

            // example of complex multiplication:
            // (x + yj) * (a + bj) = (xa - yb) + (xb + ya)j
            // So the real part of the product is (xa - yb), and the imaginary part of the product is (xb + ya).
            for (size_t j = 0; j < WINDOW_SIZE / (2 * c); j++)
            {
                vecHw.mul32(realC, 2 * WINDOW_SIZE + 32 * i + (WINDOW_SIZE - 2 * c), buf0, std::min(c, 32)); // real [C - 2C] * cos(C-2C) -> 2 * WINDOW_SIZE ( we're using this as a buffer cause we're done using it in this fft rn)
                vecHw.mul32(imagC, 3 * WINDOW_SIZE + 32 * i + (WINDOW_SIZE - 2 * c), buf1, std::min(c, 32)); // imag [C - 2C] *  sin(C-2C) -> 3 * WINDOW_SIZE 

                vecHw.mul32(buf1, minusOne, buf1, std::min(c, 32)); // (-1)* (imag [C - 2C] *  sin(C-2C)) THIS CAN BE OPTIMIZED BY PRE-NEGATING THE SIN TERMS

                vecHw.mul32(imagC, 2 * WINDOW_SIZE + 32 * i + (WINDOW_SIZE - 2 * c), imagC, std::min(c, 32)); // imag [C - 2C] = imag [C - 2C] * cos(C-2C)

                vecHw.mul32(realC, 3 * WINDOW_SIZE + 32 * i + (WINDOW_SIZE - 2 * c), buf2, std::min(c, 32)); // imag [C - 2C] = real[C - 2C] * sin(C-2C)

                vecHw.add32(buf0, buf1, realC, std::min(c, 32)); // real [C - 2C] = real [C - 2C] * cos(C-2C)  + (-1)* (imag [C - 2C] *  sin(C-2C))

                vecHw.add32(imagC, buf2, imagC, std::min(c, 32)); // imag [C - 2C] = (imag [C - 2C] * cos(C-2C)) + (real[C - 2C] * sin(C-2C))

                // repeat for 3C, 5C, and 7C
                realC += 2 * c;
                imagC += 2 * c;
            }
        }
    }

    vecHw.copyFromHw(inputReal, 0, WINDOW_SIZE, 0);
    vecHw.copyFromHw(inputImag, WINDOW_SIZE, WINDOW_SIZE, 0);

    // fix order of arrays
    ec::Float temp(0.0f);
    for (size_t i = 0; i < WINDOW_SIZE / 2; i++)
    {
        uint16_t newI = reverse_bits(i);

        if (i < newI)
        {
            memcpy(&temp, inputReal.data() + i, sizeof(ec::Float));
            memcpy(inputReal.data() + i, inputReal.data() + newI, sizeof(ec::Float));
            memcpy(inputReal.data() + newI, &temp, sizeof(ec::Float));

            memcpy(&temp, inputImag.data() + i, sizeof(ec::Float));
            memcpy(inputImag.data() + i, inputImag.data() + newI, sizeof(ec::Float));
            memcpy(inputImag.data() + newI, &temp, sizeof(ec::Float));
        }
    }

    // vecHw.copyFromHw(inputReal, 0, WINDOW_SIZE, 0);
    // vecHw.copyFromHw(inputImag, WINDOW_SIZE, WINDOW_SIZE, 0);
    // for (size_t i = 0; i < WINDOW_SIZE; i++)
    // {
    //     std::cout << "i: " << i << " " << inputReal[i].toFloat() << " " << inputImag[i].toFloat() << std::endl;
    // }


    // std::cout << inputReal[0].toFloat() << " " << inputReal[1].toFloat() << " " << inputReal[2].toFloat() << " " << inputReal[3].toFloat() << std::endl;
    // std::cout << inputImag[0].toFloat() << " " << inputImag[1].toFloat() << " " << inputImag[2].toFloat() << " " << inputImag[3].toFloat() << std::endl;
    // std::cout << inputReal[256].toFloat() << " " << inputReal[257].toFloat() << " " << inputReal[258].toFloat() << " " << inputReal[259].toFloat() << std::endl;
    // exit(1);
}

void init_angleTerms()
{
    size_t idx = 0;

    ec::Float aC(float(-2.0f * M_PI) / WINDOW_SIZE);

    for (size_t i = 0; i < WINDOW_SIZE / 2; i++)
    {
        angleTerms[idx] = ec_cos(aC * i);
        angleTerms[idx + WINDOW_SIZE] = ec_sin(aC * i);
        idx++;
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

/*
    We are computing this operation:

    ec::Float blackmanWinCoef = 0.42f - 0.5f * ec_cos(ec::Float(i) * 2.0f * PI / (WINDOW_SIZE - 1));
    blackmanWinCoef = blackmanWinCoef + 0.08f * ec_cos(ec::Float(i) * 4.0f * PI / (WINDOW_SIZE - 1));

    across the entire input vector.
    this will in-place modify the input vector.

*/
void init_blackmanCoefs(std::vector<ec::Float>& input)
{
    ec::VecHw& vecHw = *ec::VecHw::getSingletonVecHw();
    vecHw.resetMemTo0();

    // we will be assuming input will be size WINDOW_SIZE
    assert(input.size() == WINDOW_SIZE);
    // but we should double check that

    const size_t halfSize = WINDOW_SIZE / 2;


    // micro optimization, we can skip 0 because vecHw memory is reset to 0
    for (size_t i = 1; i < halfSize; ++i)
    {
        input[i] = i;
    }

    vecHw.copyToHw(input, 1, halfSize - 1, 1);

    // TODO see if we can pipeline this to be faster
    // multiplaction by constants
    for (size_t i = 0; i < WINDOW_SIZE / 64; ++i)
    {
        vecHw.mul32(32 * i, ec::Float(float(4.0f * M_PI) / (WINDOW_SIZE - 1)), WINDOW_SIZE + 32 * i);
        vecHw.mul32(32 * i, ec::Float(float(2.0f * M_PI) / (WINDOW_SIZE - 1)), 32 * i);
    }

    // in-place cosine operations
    for (size_t i = 0; i < WINDOW_SIZE / 8; ++i)
    {
        vecHw.cos4(4 * i, 4 * i);
        vecHw.cos4(WINDOW_SIZE + 4 * i, WINDOW_SIZE + 4 * i);
    }

    /*
        outputs of vector calculated cosine values:
        0 - WINDOW_SIZE-1:               ec_cos(ec::Float(i) * 2.0f * PI / (WINDOW_SIZE - 1)); (first cosine)
        WINDOW_SIZE - 2*WINDOW_SIZE - 1: ec_cos(ec::Float(i) * 4.0f * PI / (WINDOW_SIZE - 1)); (second cosine)
    */
    std::vector<ec::Float> cosOut(2 * WINDOW_SIZE);

    vecHw.copyFromHw(cosOut, 0, halfSize, 0);
    vecHw.copyFromHw(cosOut, WINDOW_SIZE, halfSize, WINDOW_SIZE);

    for (size_t i = 0; i < halfSize; i++)
    {
        memcpy(cosOut.data() + halfSize + i, cosOut.data() + halfSize - i - 1, sizeof(ec::Float));
        memcpy(cosOut.data() + WINDOW_SIZE + halfSize + i, cosOut.data() + WINDOW_SIZE + halfSize - i - 1, sizeof(ec::Float));
    }

    ec::StreamHw& streamHw = *ec::StreamHw::getSingletonStreamHw();
    streamHw.resetStreamHw();

    streamHw.copyToHw(cosOut, 0, 2 * WINDOW_SIZE, 0);

    streamHw.createFifos(6);

    // -0.5 * first cosine
    streamHw.addOpMulToPipeline(0, ec::Float(-0.5f), 1);

    // 0.42 + first cosine
    streamHw.addOpAddToPipeline(1, ec::Float(0.42f), 2);

    // 0.08 * second cosine
    streamHw.addOpMulToPipeline(3, ec::Float(0.08f), 4);

    // first cosine + second cosine
    streamHw.addOpAddToPipeline(2, 4, 5);

    // stream first cosine into fifo 0
    streamHw.startStreamDataMemToFifo(0, 0, WINDOW_SIZE);

    // stream second cosine into fifo 3
    streamHw.startStreamDataMemToFifo(WINDOW_SIZE, 3, WINDOW_SIZE);

    // stream output from fifo 5 to memory
    streamHw.startStreamDataFifoToMem(5, 3 * WINDOW_SIZE, WINDOW_SIZE);

    streamHw.runPipeline();

    streamHw.copyFromHw(input, 3 * WINDOW_SIZE, WINDOW_SIZE, 0);
}
