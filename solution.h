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
void init_angleTerms(std::vector<std::vector<ec::Float>>& cosInput, std::vector<std::vector<ec::Float>>& sinInput);
void compute_fourier_transform(const std::vector<ec::Float>& input, std::vector<std::vector<ec::Float>>& cosTerms, std::vector<std::vector<ec::Float>>& sinTerms, std::vector<ec::Float>& outputReal, std::vector<ec::Float>& outputImag);

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
    std::vector<std::vector<ec::Float>> cosAngleTerms(WINDOW_SIZE, std::vector<ec::Float>(WINDOW_SIZE));
    std::vector<std::vector<ec::Float>> sinAngleTerms(WINDOW_SIZE, std::vector<ec::Float>(WINDOW_SIZE));
    std::vector<ec::Float> outputSpectrum(sizeSpectrum, std::numeric_limits<float>::lowest());

    size_t idxStartWin = 0;

    init_blackmanCoefs(blackmanCoefs);
    init_angleTerms(cosAngleTerms, sinAngleTerms);

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

// TODO pipeline and vectorize this
void init_angleTerms(std::vector<std::vector<ec::Float>>& cosInput, std::vector<std::vector<ec::Float>>& sinInput)
{
    const ec::Float angleConst = (-2.0f * PI) * (1.0f / ec::Float(WINDOW_SIZE));

    // ec::StreamHw& streamHw = *ec::StreamHw::getSingletonStreamHw();
    // streamHw.resetStreamHw();


    for (size_t i = 0; i < WINDOW_SIZE; ++i)
    {
        if (i == 0)
        {
            for (size_t j = 0; j < WINDOW_SIZE; ++j)
            {
                cosInput[i][j] = 1.0f;
            }

            continue;
        }

        // i != 0 from now on
        for (size_t j = 0; j < WINDOW_SIZE; ++j)
        {
            if (j == 0)
            {
                cosInput[i][j] = 1.0f;
            }
            else
            {
                cosInput[i][j] = ec_cos(angleConst * i * j);
            }
        }

        // sin is just an offset of cos. let's take advantage of that
        if (WINDOW_SIZE % (4 * i) == 0)
        {
            size_t offset = WINDOW_SIZE / (4 * i);

            for (size_t j = 0; j < WINDOW_SIZE - offset; j++)
            {
                sinInput[i][j] = cosInput[i][j + offset];
            }

            for (size_t j = 0; j < offset; j++)
            {
                sinInput[i][j + WINDOW_SIZE - offset] = cosInput[i][j];
            }
        }
        else
        {
            for (size_t j = 1; j < WINDOW_SIZE; j++)
            {
                sinInput[i][j] = ec_sin(angleConst * i * j);
            }
        }
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

void compute_fourier_transform(const std::vector<ec::Float>& input, std::vector<std::vector<ec::Float>>& cosTerms, std::vector<std::vector<ec::Float>>& sinTerms, std::vector<ec::Float>& outputReal, std::vector<ec::Float>& outputImag)
{
    // we will be assuming input will be size WINDOW_SIZE
    assert(input.size() == WINDOW_SIZE);
    // but we should double check that

    outputReal.clear();
    outputReal.resize(WINDOW_SIZE, 0.0f);
    outputImag.clear();
    outputImag.resize(WINDOW_SIZE, 0.0f);

    const ec::Float angleConst = (-2.0f * PI) * (1.0f / ec::Float(WINDOW_SIZE));

    // outputReal[i] += input[j] * cosTerms[i][j];
    // outputImag[i] += input[j] * sinTerms[i][j];

    ec::StreamHw& streamHw = *ec::StreamHw::getSingletonStreamHw();
    streamHw.resetStreamHw();

    streamHw.createFifos(2 * WINDOW_SIZE);

    size_t in_idx = 0;
    size_t i = 0;


    streamHw.addOpMulToPipeline(0, input[0], 1);
    streamHw.addOpMulToPipeline(2, input[1], 3);
    streamHw.addOpAddToPipeline(1, 3, 4);

    streamHw.addOpMulToPipeline(5, input[2], 6);
    streamHw.addOpAddToPipeline(4, 6, 7);

    streamHw.addOpMulToPipeline(8, input[3], 9);
    streamHw.addOpAddToPipeline(7, 9, 10);

    streamHw.copyToHw(cosTerms[0], 0, WINDOW_SIZE, 0);
    streamHw.copyToHw(cosTerms[1], 0, WINDOW_SIZE, WINDOW_SIZE);
    streamHw.copyToHw(cosTerms[2], 0, WINDOW_SIZE, 2 * WINDOW_SIZE);
    streamHw.copyToHw(cosTerms[3], 0, WINDOW_SIZE, 3 * WINDOW_SIZE);

    streamHw.startStreamDataMemToFifo(0, 0, WINDOW_SIZE);
    streamHw.startStreamDataMemToFifo(WINDOW_SIZE, 2, WINDOW_SIZE);
    streamHw.startStreamDataMemToFifo(2 * WINDOW_SIZE, 5, WINDOW_SIZE);
    streamHw.startStreamDataMemToFifo(3 * WINDOW_SIZE, 8, WINDOW_SIZE);

    

    // streamHw.startStreamDataMemToFifo(0, i, WINDOW_SIZE);
    // streamHw.addOpMulToPipeline(0, input[in_idx], i + 1);
    // i++;
    // in_idx++;

    // for (; in_idx < 30; i += 4)
    // {
    //     streamHw.addOpMulToPipeline(i + 2, input[in_idx], i + 3);
    //     streamHw.startStreamDataMemToFifo(0, i + 2, WINDOW_SIZE);
    //     in_idx++;

    //     streamHw.addOpAddToPipeline(i, i + 3, i + 4);
    // }


    // // costerms to fifo

    // // stream fifo back to memory
    streamHw.startStreamDataFifoToMem(10, 0, WINDOW_SIZE);

    // copy input to memory
    //streamHw.copyToHw(input, 0, WINDOW_SIZE, WINDOW_SIZE);

    streamHw.runPipeline();

    streamHw.copyFromHw(outputReal, 0, WINDOW_SIZE, 0);

    std::cout << "A: " << outputReal[0].toFloat() << " " << outputReal[1].toFloat() << " " << outputReal[2].toFloat() << " " << outputReal[3].toFloat() << std::endl;


    // in_idx = 0;
    // i = 0;

    // streamHw.startStreamDataMemToFifo(0, i, WINDOW_SIZE);
    // //streamHw.addOpMulToPipeline(0, input[in_idx], i + 1);
    // i++;
    // in_idx++;

    // for (; in_idx < 30; i += 4)
    // {
    //     // streamHw.addOpMulToPipeline(i + 2, input[in_idx], i + 3);
    //     streamHw.startStreamDataMemToFifo(0, i + 2, WINDOW_SIZE);
    //     in_idx++;

    //     //  streamHw.addOpAddToPipeline(i, i + 3, i + 4);
    // }

    // streamHw.copyToHw(cosTerms[2], 0, WINDOW_SIZE, 0);

    // streamHw.startStreamDataFifoToMem(i, WINDOW_SIZE, WINDOW_SIZE);

    // // copy input to memory
    // //streamHw.copyToHw(input, 0, WINDOW_SIZE, WINDOW_SIZE);

    // streamHw.runPipeline();

    // streamHw.copyFromHw(outputReal, WINDOW_SIZE, WINDOW_SIZE, 0);

    // std::cout << "B: " << outputReal[0].toFloat() << " " << outputReal[1].toFloat() << " " << outputReal[2].toFloat() << " " << outputReal[3].toFloat() << std::endl;
    // std::vector<ec::Float> test{0.5f, 0.6f, 0.7f, 0.8f};

    // std::vector<ec::Float> in{0.5f, 0.6f, 0.7f, 0.8f};


    // streamHw.addOpMulToPipeline(0, in[0], 1);
    // streamHw.addOpMulToPipeline(2, in[1], 3);

    // streamHw.addOpAddToPipeline(1, 3, 4);

    // streamHw.addOpMulToPipeline(5, in[2], 6);
    // streamHw.addOpMulToPipeline(7, in[3], 8);

    // streamHw.addOpAddToPipeline(6, 8, 9);

    // streamHw.addOpAddToPipeline(4, 9, 10);

    // streamHw.copyToHw(test, 0, 4, 0);

    // streamHw.startStreamDataMemToFifo(0, 0, 4);
    // streamHw.startStreamDataMemToFifo(0, 2, 4);
    // streamHw.startStreamDataMemToFifo(0, 5, 4);
    // streamHw.startStreamDataMemToFifo(0, 7, 4);

    // streamHw.startStreamDataFifoToMem(10, 5, 4);

    // std::vector<ec::Float> o(4);

    // streamHw.runPipeline();

    // streamHw.copyFromHw(o, 5, 4, 0);

    // std::cout << o[0].toFloat() << " " << o[1].toFloat() << " " << o[2].toFloat() << " " << o[3].toFloat() << std::endl;





    // streamHw.resetStreamHw();

    // streamHw.createFifos(2 * WINDOW_SIZE);


    // streamHw.addOpAddToPipeline(0, 1, 2);

    // i = 2;
    // for (; in_idx < WINDOW_SIZE; i += 5)
    // {
    //     streamHw.addOpMulToPipeline(i, input[in_idx], i + 1);
    //     in_idx++;

    //     streamHw.addOpMulToPipeline(i + 2, input[in_idx], i + 3);
    //     in_idx++;

    //     streamHw.addOpAddToPipeline(i + 1, i + 3, i + 4);
    // }



    // for (size_t i = 0; i < WINDOW_SIZE; ++i)
    // {
    //     // copy current output back into mem @ 2*window size
    //     streamHw.copyToHw(outputReal, 0, WINDOW_SIZE, 2 * WINDOW_SIZE);

    //     // copy next set of cosine terms into mem 0
    //     streamHw.copyToHw(cosTerms[i], 0, WINDOW_SIZE, 0);

    //     streamHw.runPipeline();

    //     // pull the output from the streamhw so we can push it back
    //     streamHw.copyFromHw(outputReal, 3 * WINDOW_SIZE, WINDOW_SIZE, 0);
    // }

    //std::cout << outputReal[0].toFloat() << " " << outputReal[1].toFloat() << " " << outputReal[2].toFloat() << " " << outputReal[3].toFloat() << std::endl;

    outputReal.clear();
    outputReal.resize(WINDOW_SIZE, 0.0f);

    for (size_t i = 0; i < WINDOW_SIZE; i++)
    {
        for (size_t j = 0; j < 4; j++)
        {
            outputReal[i] += input[j] * cosTerms[i][j];
            outputImag[i] += input[j] * sinTerms[i][j];
        }
    }

    std::cout << "C: " << outputReal[0].toFloat() << " " << outputReal[1].toFloat() << " " << outputReal[2].toFloat() << " " << outputReal[3].toFloat() << std::endl;

    // for (size_t i = 0; i < WINDOW_SIZE; i++)
    // {
    //     for (size_t j = 0; j < WINDOW_SIZE; j++)
    //     {
    //         // outputReal[i] += input[j] * cosTerms[i][j];
    //          //outputImag[i] += input[j] * sinTerms[i][j];
    //     }
    // }

    return;
}

