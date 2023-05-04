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


    ec::VecHw& vecHw = *ec::VecHw::getSingletonVecHw();
    vecHw.resetMemTo0();

    vecHw.copyToHw(angleTerms, 0, 2 * WINDOW_SIZE - 64, 0);

    for (size_t j = 0; j < numWins; j++)
    {
        for (size_t i = 0; i < WINDOW_SIZE; i++)
        {
            signalWindow[i] = inputSignal[i + idxStartWin] * blackmanCoefs[i];
        }

        std::vector<ec::Float> inImag(WINDOW_SIZE);

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

void fft(std::vector<ec::Float>& inputReal, std::vector<ec::Float>& inputImag, size_t count)
{
    size_t halfCount = count / 2;
    size_t quarterCount = halfCount / 2;

    std::vector<ec::Float> even(halfCount);
    std::vector<ec::Float> evenI(halfCount);

    std::vector<ec::Float> odd(halfCount);
    std::vector<ec::Float> oddI(halfCount);

    if (count > 4)
    {
        for (size_t i = 0; i < halfCount; i++)
        {
            memcpy(even.data() + i, inputReal.data() + i * 2, sizeof(ec::Float));
            memcpy(evenI.data() + i, inputImag.data() + i * 2, sizeof(ec::Float));

            memcpy(odd.data() + i, inputReal.data() + i * 2 + 1, sizeof(ec::Float));
            memcpy(oddI.data() + i, inputImag.data() + i * 2 + 1, sizeof(ec::Float));
        }

        fft(even, evenI, halfCount);
        fft(odd, oddI, halfCount);
    }
    else
    {
        even[0] = inputReal[0] + inputReal[2];
        evenI[0] = inputImag[0] + inputImag[2];

        even[1] = inputReal[0] - inputReal[2];
        evenI[1] = inputImag[0] - inputImag[2];

        odd[0] = inputReal[1] + inputReal[3];
        oddI[0] = inputImag[1] + inputImag[3];

        odd[1] = inputReal[1] - inputReal[3];
        oddI[1] = inputImag[1] - inputImag[3];
    }

    // negative returns if we run it with too small of data because of copy overhead
    if (halfCount >= 128)
    {
        std::vector<ec::Float> big(count);
        memcpy(big.data(), odd.data(), sizeof(ec::Float) * halfCount);
        memcpy(big.data() + halfCount, oddI.data(), sizeof(ec::Float) * halfCount);

        ec::VecHw& vecHw = *ec::VecHw::getSingletonVecHw();

        vecHw.copyToHw(big, 0, count, 2 * WINDOW_SIZE);

        for (size_t i = 0; i < halfCount / 32; i++)
        {
            vecHw.mul32(WINDOW_SIZE - count + 32 * i, 2 * WINDOW_SIZE + 32 * i, 3 * WINDOW_SIZE + 32 * i); // co * odd
            vecHw.mul32(2 * WINDOW_SIZE - count + 32 * i, (2 * WINDOW_SIZE + halfCount) + 32 * i, (3 * WINDOW_SIZE + 512) + 32 * i); // si * oddI
            vecHw.mul32((3 * WINDOW_SIZE + 512) + 32 * i, minusOne, (3 * WINDOW_SIZE + 512) + 32 * i); // -1 * (si * oddI)

            vecHw.add32(3 * WINDOW_SIZE + 32 * i, (3 * WINDOW_SIZE + 512) + 32 * i, 3 * WINDOW_SIZE + 32 * i); // c1 is at 3 * WINDOW_SIZE

            vecHw.mul32(WINDOW_SIZE - count + 32 * i, 2 * WINDOW_SIZE + halfCount + 32 * i, 2 * WINDOW_SIZE + 512 + 32 * i); // co * oddI... we are done with oddI so we can overwrite
            vecHw.mul32(2 * WINDOW_SIZE - count + 32 * i, 2 * WINDOW_SIZE + 32 * i, 2 * WINDOW_SIZE + 32 * i); // si * odd... we are done with odd so we can overwrite

            vecHw.add32(2 * WINDOW_SIZE + 32 * i, 2 * WINDOW_SIZE + 512 + 32 * i, 3 * WINDOW_SIZE + 512 + 32 * i); // c2 is at 3 * WINDOW_SIZE + 512
        }

        memcpy(big.data(), even.data(), sizeof(ec::Float) * halfCount);
        memcpy(big.data() + halfCount, evenI.data(), sizeof(ec::Float) * halfCount);

        vecHw.copyToHw(big, 0, count, 2 * WINDOW_SIZE);

        for (size_t i = 0; i < halfCount / 32; i++)
        {
            vecHw.add32(2 * WINDOW_SIZE + 32 * i, 3 * WINDOW_SIZE + 32 * i, 2 * WINDOW_SIZE + 32 * i); // even + c1
            vecHw.add32(2 * WINDOW_SIZE + halfCount + 32 * i, 3 * WINDOW_SIZE + 512 + 32 * i, 2 * WINDOW_SIZE + halfCount + 32 * i); // evenI + c2
        }

        vecHw.copyFromHw(big, 2 * WINDOW_SIZE + 1, count - 1, 0);
        memcpy(inputReal.data() + 1, big.data(), sizeof(ec::Float) * (halfCount - 1));
        memcpy(inputImag.data() + 1, big.data() + halfCount, sizeof(ec::Float) * (halfCount - 1));

        for (size_t i = 0; i < halfCount / 32; i++)
        {
            vecHw.mul32(3 * WINDOW_SIZE + 32 * i, minusTwo, 3 * WINDOW_SIZE + 32 * i); // -2 * c1
            vecHw.mul32(3 * WINDOW_SIZE + 512 + 32 * i, minusTwo, 3 * WINDOW_SIZE + 512 + 32 * i); // -2 * c2

            vecHw.add32(2 * WINDOW_SIZE + 32 * i, 3 * WINDOW_SIZE + 32 * i, 2 * WINDOW_SIZE + 32 * i); // even + c1 - 2*c1
            vecHw.add32(2 * WINDOW_SIZE + halfCount + 32 * i, 3 * WINDOW_SIZE + 512 + 32 * i, 2 * WINDOW_SIZE + halfCount + 32 * i); // evenI + c2 - 2*c2
        }

        vecHw.copyFromHw(big, 2 * WINDOW_SIZE + 1, count - 1, 0);
        memcpy(inputReal.data() + halfCount + 1, big.data(), sizeof(ec::Float) * (halfCount - 1));
        memcpy(inputImag.data() + halfCount + 1, big.data() + halfCount, sizeof(ec::Float) * (halfCount - 1));
    }
    else
    {
        std::vector<ec::Float> v1Cache(quarterCount);
        std::vector<ec::Float> v2Cache(quarterCount);
        std::vector<ec::Float> c2Cache(quarterCount);

        for (size_t k = 1; k < halfCount; k++)
        {
            ec::Float c1;
            ec::Float c2;

            if (k == quarterCount)
            {
                memcpy(&inputReal[k], &even[k], sizeof(ec::Float));
                memcpy(&inputReal[halfCount + k], &even[k], sizeof(ec::Float));

                c2 = angleTerms[2 * WINDOW_SIZE - count + k] * odd[k];

                inputImag[k] = evenI[k] + c2;
                inputImag[halfCount + k] = evenI[k] - c2;
            }
            else
            {
                if (k > quarterCount)
                {
                    c1 = (v2Cache[halfCount - k] - v1Cache[halfCount - k]);

                    inputImag[k] = evenI[k] + c2Cache[halfCount - k];
                    inputImag[halfCount + k] = evenI[k] - c2Cache[halfCount - k];
                }
                else
                {
                    v1Cache[k] = angleTerms[WINDOW_SIZE - count + k] * odd[k];
                    v2Cache[k] = angleTerms[2 * WINDOW_SIZE - count + k] * oddI[k];
                    c1 = (v1Cache[k] - v2Cache[k]);

                    c2Cache[k] = (angleTerms[WINDOW_SIZE - count + k] * oddI[k] + angleTerms[2 * WINDOW_SIZE - count + k] * odd[k]);

                    inputImag[k] = evenI[k] + c2Cache[k];
                    inputImag[halfCount + k] = evenI[k] - c2Cache[k];
                }

                inputReal[k] = even[k] + c1;
                inputReal[halfCount + k] = even[k] - c1;
            }
        }
    }

    inputReal[0] = even[0] + odd[0];
    inputImag[0] = evenI[0] + oddI[0];

    inputReal[halfCount] = even[0] - odd[0];
    inputImag[halfCount] = evenI[0] - oddI[0];
}

void init_angleTerms()
{
    size_t count = WINDOW_SIZE;

    size_t idx = 0;

    while (count >= 2)
    {
        ec::Float aC(float(-2.0f * M_PI) / count);

        for (size_t i = 0; i < count / 2; i++)
        {
            if (i != 0)
            {
                if (!(idx >= 977 && idx <= 991 || idx >= 1001 && idx <= 1007 || idx >= 1013 && idx <= 1015 || idx == 1019))
                {
                    angleTerms[idx] = ec_cos(aC * i);
                    angleTerms[idx + WINDOW_SIZE] = ec_sin(aC * i);
                }
            }

            idx++;
        }

        count >>= 1;
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
