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
void init_angleTerms(std::vector<ec::Float>& cosInput, std::vector<ec::Float>& sinInput);
void compute_fourier_transform(const std::vector<ec::Float>& input, std::vector<ec::Float>& cosTerms, std::vector<ec::Float>& sinTerms, std::vector<ec::Float>& outputReal, std::vector<ec::Float>& outputImag);

std::vector<ec::Float> process_signal(const std::vector<ec::Float>& inputSignal)
{
    const size_t numSamples = inputSignal.size(); // assume divisible by 32
    const size_t sizeSpectrum = (WINDOW_SIZE / 2) + 1;
    const size_t stepBetweenWins = static_cast<size_t>(ceil(WINDOW_SIZE * (1 - OVERLAP_RATIO)));
    const size_t numWins = (numSamples - WINDOW_SIZE) / stepBetweenWins + 1;

    std::vector<ec::Float> signalWindow(WINDOW_SIZE);
    std::vector<ec::Float> signalFreqReal(WINDOW_SIZE);
    std::vector<ec::Float> signalFreqImag(WINDOW_SIZE);
    std::vector<ec::Float> spectrumWindow(sizeSpectrum);
    std::vector<ec::Float> blackmanCoefs(WINDOW_SIZE);

    // TODO Shrink these to exclude the 1's and 0's and hardcode
    std::vector<ec::Float> cosAngleTerms(WINDOW_SIZE * WINDOW_SIZE);
    std::vector<ec::Float> sinAngleTerms(WINDOW_SIZE * WINDOW_SIZE);
    std::vector<ec::Float> outputSpectrum(sizeSpectrum, std::numeric_limits<float>::lowest());

    size_t idxStartWin = 0;

    init_blackmanCoefs(blackmanCoefs);
    init_angleTerms(cosAngleTerms, sinAngleTerms);

    const ec::Float spC0((float)(10.0f / log(10.0f)));
    const ec::Float spC1((float)(10.0f * log(125.0f / 131072.0f) / log(10.0f)));
    const ec::Float spC2((float)(10.0f * log(125.0f / 32768.0f) / log(10.0f)));

    for (size_t j = 0; j < numWins; j++)
    {
        for (size_t i = 0; i < WINDOW_SIZE; i++)
        {
            signalWindow[i] = inputSignal[i + idxStartWin] * blackmanCoefs[i];
        }

        compute_fourier_transform(signalWindow, cosAngleTerms, sinAngleTerms, signalFreqReal, signalFreqImag);

        for (size_t i = 0; i < sizeSpectrum; i++)
        {
            ec::Float freqVal = signalFreqReal[i] * signalFreqReal[i] + signalFreqImag[i] * signalFreqImag[i];
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

// TODO pipeline and vectorize this
void init_angleTerms(std::vector<ec::Float>& cosInput, std::vector<ec::Float>& sinInput)
{
    const ec::Float angleConst = (-2.0f * PI) * (1.0f / ec::Float(WINDOW_SIZE));

    ec::VecHw& vecHw = *ec::VecHw::getSingletonVecHw();
    vecHw.resetMemTo0();

    for (size_t i = 1; i < WINDOW_SIZE; ++i)
    {
        cosInput[i] = i;
    }

    vecHw.copyToHw(cosInput, 1, WINDOW_SIZE - 1, 1);

    for (size_t i = 0; i < WINDOW_SIZE / 32; ++i)
    {
        vecHw.mul32(32 * i, angleConst, 32 * i);
    }

    for (size_t j = 1; j < WINDOW_SIZE; j++)
    {
        ec::Float jC((float)j);

        for (size_t i = 0; i < WINDOW_SIZE / 32; ++i)
        {
            vecHw.mul32(32 * i, jC, 32 * i + WINDOW_SIZE);
        }

        for (size_t i = 0; i < WINDOW_SIZE / 4; ++i)
        {
            vecHw.sin4(WINDOW_SIZE + 4 * i, 2 * WINDOW_SIZE + 4 * i);
            vecHw.cos4(WINDOW_SIZE + 4 * i, WINDOW_SIZE + 4 * i);
        }

        vecHw.copyFromHw(cosInput, WINDOW_SIZE, WINDOW_SIZE, WINDOW_SIZE * j);
        vecHw.copyFromHw(sinInput, 2 * WINDOW_SIZE, WINDOW_SIZE, WINDOW_SIZE * j);
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

void compute_fourier_transform(const std::vector<ec::Float>& input, std::vector<ec::Float>& cosTerms, std::vector<ec::Float>& sinTerms, std::vector<ec::Float>& outputReal, std::vector<ec::Float>& outputImag)
{
    // we will be assuming input will be size WINDOW_SIZE
    assert(input.size() == WINDOW_SIZE);
    // but we should double check that

    ec::StreamHw& streamHw = *ec::StreamHw::getSingletonStreamHw();
    streamHw.resetStreamHw();

    streamHw.createFifos(11);

    streamHw.addOpMulToPipeline(0, input[0], 1);
    streamHw.addOpMulToPipeline(2, input[1], 3);
    streamHw.addOpAddToPipeline(1, 3, 4);

    streamHw.addOpMulToPipeline(5, input[2], 6);
    streamHw.addOpAddToPipeline(4, 6, 7);

    streamHw.addOpMulToPipeline(8, input[3], 9);
    streamHw.addOpAddToPipeline(7, 9, 10);

    streamHw.copyToHw(cosTerms, 0, 4 * WINDOW_SIZE, 0);

    streamHw.startStreamDataMemToFifo(0, 0, WINDOW_SIZE);
    streamHw.startStreamDataMemToFifo(WINDOW_SIZE, 2, WINDOW_SIZE);
    streamHw.startStreamDataMemToFifo(2 * WINDOW_SIZE, 5, WINDOW_SIZE);
    streamHw.startStreamDataMemToFifo(3 * WINDOW_SIZE, 8, WINDOW_SIZE);

    streamHw.startStreamDataFifoToMem(10, 0, WINDOW_SIZE);

    streamHw.runPipeline();

    streamHw.copyFromHw(outputReal, 0, WINDOW_SIZE, 0);

    streamHw.copyToHw(sinTerms, 0, 4 * WINDOW_SIZE, 0);

    streamHw.startStreamDataMemToFifo(0, 0, WINDOW_SIZE);
    streamHw.startStreamDataMemToFifo(WINDOW_SIZE, 2, WINDOW_SIZE);
    streamHw.startStreamDataMemToFifo(2 * WINDOW_SIZE, 5, WINDOW_SIZE);
    streamHw.startStreamDataMemToFifo(3 * WINDOW_SIZE, 8, WINDOW_SIZE);

    streamHw.startStreamDataFifoToMem(10, 0, WINDOW_SIZE);

    streamHw.runPipeline();

    streamHw.copyFromHw(outputImag, 0, WINDOW_SIZE, 0);

    for (size_t i = 4; i < WINDOW_SIZE; i += 3)
    {
        streamHw.resetStreamHw();

        streamHw.createFifos(10);

        streamHw.addOpMulToPipeline(0, input[i], 1);
        streamHw.addOpMulToPipeline(2, input[i + 1], 3);
        streamHw.addOpAddToPipeline(1, 3, 4);

        streamHw.addOpMulToPipeline(5, input[i + 2], 6);
        streamHw.addOpAddToPipeline(4, 6, 7);

        streamHw.addOpAddToPipeline(7, 8, 9);

        streamHw.copyToHw(cosTerms, i * WINDOW_SIZE, 3 * WINDOW_SIZE, 0);
        streamHw.copyToHw(outputReal, 0, WINDOW_SIZE, 3 * WINDOW_SIZE);

        streamHw.startStreamDataMemToFifo(0, 0, WINDOW_SIZE);
        streamHw.startStreamDataMemToFifo(WINDOW_SIZE, 2, WINDOW_SIZE);
        streamHw.startStreamDataMemToFifo(2 * WINDOW_SIZE, 5, WINDOW_SIZE);
        streamHw.startStreamDataMemToFifo(3 * WINDOW_SIZE, 8, WINDOW_SIZE);

        streamHw.startStreamDataFifoToMem(9, 0, WINDOW_SIZE);

        streamHw.runPipeline();

        streamHw.copyFromHw(outputReal, 0, WINDOW_SIZE, 0);

        // imaginary terms
        streamHw.copyToHw(sinTerms, i * WINDOW_SIZE, 3 * WINDOW_SIZE, 0);
        streamHw.copyToHw(outputImag, 0, WINDOW_SIZE, 3 * WINDOW_SIZE);

        streamHw.startStreamDataMemToFifo(0, 0, WINDOW_SIZE);
        streamHw.startStreamDataMemToFifo(WINDOW_SIZE, 2, WINDOW_SIZE);
        streamHw.startStreamDataMemToFifo(2 * WINDOW_SIZE, 5, WINDOW_SIZE);
        streamHw.startStreamDataMemToFifo(3 * WINDOW_SIZE, 8, WINDOW_SIZE);

        streamHw.startStreamDataFifoToMem(9, 0, WINDOW_SIZE);

        streamHw.runPipeline();

        streamHw.copyFromHw(outputImag, 0, WINDOW_SIZE, 0);
    }

    return;
}

