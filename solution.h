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
void compute_fourier_transform(const std::vector<ec::Float>& input, std::vector<ec::Float>& outputReal, std::vector<ec::Float>& outputImag);

std::vector<ec::Float> process_signal(const std::vector<ec::Float>& inputSignal)
{
    const size_t numSamples = inputSignal.size();
    const size_t sizeSpectrum = (WINDOW_SIZE / 2) + 1;
    const size_t stepBetweenWins = static_cast<size_t>(ceil(WINDOW_SIZE * (1 - OVERLAP_RATIO)));
    const size_t numWins = (numSamples - WINDOW_SIZE) / stepBetweenWins + 1;


    std::vector<ec::Float> signalWindow(WINDOW_SIZE);
    std::vector<ec::Float> signalFreqReal(WINDOW_SIZE);
    std::vector<ec::Float> signalFreqImag(WINDOW_SIZE);
    std::vector<ec::Float> spectrumWindow(sizeSpectrum);
    std::vector<ec::Float> blackmanCoefs(WINDOW_SIZE);
    std::vector<ec::Float> outputSpectrum(sizeSpectrum, std::numeric_limits<float>::lowest());

    size_t idxStartWin = 0;

    init_blackmanCoefs(blackmanCoefs);

    for (size_t j = 0; j < numWins; j++)
    {
        for (size_t i = 0; i < WINDOW_SIZE; i++)
        {
            signalWindow[i] = inputSignal[i + idxStartWin] * blackmanCoefs[i];
        }

        compute_fourier_transform(signalWindow, signalFreqReal, signalFreqImag);

        for (size_t i = 0; i < sizeSpectrum; i++)
        {
            ec::Float freqVal = signalFreqReal[i] * signalFreqReal[i] + signalFreqImag[i] * signalFreqImag[i];
            freqVal = ec_sqrt(freqVal);
            freqVal = freqVal / ec::Float(WINDOW_SIZE);

            if (i > 0 && i < sizeSpectrum - 1) freqVal = freqVal * 2.0f;

            freqVal = freqVal * freqVal;

            freqVal = 10.0f * ec_log10(1000.0f * freqVal);

            outputSpectrum[i] = ec_max(outputSpectrum[i], freqVal);
        }

        idxStartWin += stepBetweenWins;

    }

    return outputSpectrum;
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

    std::vector<ec::Float> indices(WINDOW_SIZE);

    const ec::Float c1 = 2.0f * PI / (WINDOW_SIZE - 1);
    const ec::Float c2 = 2 * c1;

    for (size_t i = 0; i < WINDOW_SIZE; i++)
    {
        indices[i] = i;
    }
    
    vecHw.copyToHw(indices, 0, WINDOW_SIZE, 0);

    // multiplaction by constants
    for (size_t i = 0; i < WINDOW_SIZE / 32; i++)
    {
        vecHw.mul32(32 * i, c2, WINDOW_SIZE + 32 * i);
        vecHw.mul32(32 * i, c1, 32 * i);
    }

    // in-place cosine operations
    for (size_t i = 0; i < WINDOW_SIZE / 4; i++)
    {
        vecHw.cos4(4 * i,  4 * i);
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

    // stream second cosing into fifo 3
    streamHw.startStreamDataMemToFifo(WINDOW_SIZE, 3, WINDOW_SIZE);

    // stream output from fifo 5 to memory
    streamHw.startStreamDataFifoToMem(5, 3 * WINDOW_SIZE, WINDOW_SIZE);

    streamHw.runPipeline();

    streamHw.copyFromHw(input, 3 * WINDOW_SIZE, WINDOW_SIZE, 0);
}

void compute_fourier_transform(const std::vector<ec::Float>& input, std::vector<ec::Float>& outputReal, std::vector<ec::Float>& outputImag)
{
    size_t inputSize = input.size();

    outputReal.clear();
    outputReal.resize(inputSize, 0.0f);
    outputImag.clear();
    outputImag.resize(inputSize, 0.0f);

    for (size_t i = 0; i < inputSize; ++i)
    {
        for (size_t j = 0; j < inputSize; ++j)
        {
            const ec::Float angleTerm = (-2.0f * PI) * ec::Float(i) * j * (1.0f / ec::Float(inputSize));

            outputReal[i] += input[j] * ec_cos(angleTerm);
            outputImag[i] += input[j] * ec_sin(angleTerm);
        }
    }

    return;
}

