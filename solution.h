//! Rohde & Schwarz Engineering Competition 2023
//!
//! This is the code to speed up. Enjoy!

#pragma once

#include "ec2023/ec2023.h"
#include <iostream>      
#include <iomanip>
#include <vector>
#include <cassert>

static constexpr float OVERLAP_RATIO = 0.75;
static constexpr size_t WINDOW_SIZE = 1024;
static const ec::Float PI = 3.14159265358979323846f;

void init_blackmanCoefs(std::vector<ec::Float>& input);
void init_angleTerms();

void fft(std::vector<ec::Float>& inputReal, std::vector<ec::Float>& inputImag, size_t count);

static std::vector<ec::Float> cosAngleTerms(WINDOW_SIZE);
static std::vector<ec::Float> sinAngleTerms(WINDOW_SIZE);

std::vector<ec::Float> process_signal(const std::vector<ec::Float>& inputSignal)
{
    const size_t numSamples = inputSignal.size(); // assume divisible by 32
    const size_t sizeSpectrum = (WINDOW_SIZE / 2) + 1;
    const size_t stepBetweenWins = static_cast<size_t>(ceil(WINDOW_SIZE * (1 - OVERLAP_RATIO)));
    const size_t numWins = (numSamples - WINDOW_SIZE) / stepBetweenWins + 1;

    std::vector<ec::Float> signalWindow(WINDOW_SIZE);
    std::vector<ec::Float> blackmanCoefs(WINDOW_SIZE);

    std::vector<ec::Float> outputSpectrum(sizeSpectrum, std::numeric_limits<float>::lowest());

    size_t idxStartWin = 0;

    init_blackmanCoefs(blackmanCoefs);
    init_angleTerms();

    const ec::Float spC0((float)(10.0f / log(10.0f)));
    const ec::Float spC1((float)(10.0f * log(125.0f / 131072.0f) / log(10.0f)));
    const ec::Float spC2((float)(10.0f * log(125.0f / 32768.0f) / log(10.0f)));

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

            outputSpectrum[i] = ec_max(outputSpectrum[i], freqVal);
        }

        idxStartWin += stepBetweenWins;
    }

    return outputSpectrum;
}

void fft(std::vector<ec::Float>& inputReal, std::vector<ec::Float>& inputImag, size_t count)
{
    std::vector<ec::Float> even(count / 2);
    std::vector<ec::Float> odd(count / 2);
    std::vector<ec::Float> evenI(count / 2);
    std::vector<ec::Float> oddI(count / 2);

    if (count > 4)
    {
        for (size_t i = 0; i < count / 2; i++)
        {
            even[i] = inputReal[i * 2];
            evenI[i] = inputImag[i * 2];

            odd[i] = inputReal[i * 2 + 1];
            oddI[i] = inputImag[i * 2 + 1];
        }

        fft(even, evenI, count / 2);
        fft(odd, oddI, count / 2);
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

    inputReal[0] = even[0] + odd[0];
    inputImag[0] = evenI[0] + oddI[0];

    inputReal[count / 2] = even[0] - odd[0];
    inputImag[count / 2] = evenI[0] - oddI[0];

    if (count / 2 >= 256 && false)
    {
        ec::VecHw& vecHw = *ec::VecHw::getSingletonVecHw();
        vecHw.resetMemTo0();
    }
    else
    {
        for (size_t k = 1; k < count / 2; k++)
        {
            ec::Float co = cosAngleTerms[WINDOW_SIZE - count + k];
            ec::Float si = sinAngleTerms[WINDOW_SIZE - count + k];

            ec::Float c1 = (co * odd[k] - si * oddI[k]);
            ec::Float c2 = (co * oddI[k] + si * odd[k]);

            inputReal[k] = even[k] + c1;
            inputImag[k] = evenI[k] + c2;

            inputReal[count / 2 + k] = even[k] - c1;
            inputImag[count / 2 + k] = evenI[k] - c2;
        }
    }
}

void init_angleTerms()
{
    size_t count = WINDOW_SIZE;

    size_t idx = 0;

    while (count >= 2)
    {
        ec::Float aC = (-2.0f * PI) * (1.0f / count);

        for (size_t i = 0; i < count / 2; i++)
        {
            if (i != 0)
            {
                cosAngleTerms[idx] = ec_cos(aC * i);
                sinAngleTerms[idx] = ec_sin(aC * i);
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

    const ec::Float c1 = 2.0f * PI / (WINDOW_SIZE - 1);
    const ec::Float c2 = 2 * c1;

    // micro optimization, we can skip 0 because vecHw memory is reset to 0
    for (size_t i = 1; i < WINDOW_SIZE; ++i)
    {
        input[i] = i;
    }

    vecHw.copyToHw(input, 1, WINDOW_SIZE - 1, 1);

    // TODO see if we can pipeline this to be faster
    // multiplaction by constants
    for (size_t i = 0; i < WINDOW_SIZE / 32; ++i)
    {
        vecHw.mul32(32 * i, c2, WINDOW_SIZE + 32 * i);
        vecHw.mul32(32 * i, c1, 32 * i);
    }

    // in-place cosine operations
    //TODO: COS WILL REPEAT EVENTUALLY
    for (size_t i = 0; i < WINDOW_SIZE / 4; ++i)
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

    vecHw.copyFromHw(cosOut, 0, 2 * WINDOW_SIZE, 0);

    ec::StreamHw& streamHw = *ec::StreamHw::getSingletonStreamHw();
    streamHw.resetStreamHw();

    streamHw.copyToHw(cosOut, 0, 2 * WINDOW_SIZE, 0);

    streamHw.createFifos(6);

    ec::Float s1(-0.5f);
    ec::Float s2(0.42f);
    ec::Float s3(0.08f);

    // -0.5 * first cosine
    streamHw.addOpMulToPipeline(0, s1, 1);

    // 0.42 + first cosine
    streamHw.addOpAddToPipeline(1, s2, 2);

    // 0.08 * second cosine
    streamHw.addOpMulToPipeline(3, s3, 4);

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
