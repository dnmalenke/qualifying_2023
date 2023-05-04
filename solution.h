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
static const ec::Float minusOne = -1.0f;
static const ec::Float minusTwo = -2.0f;

void init_blackmanCoefs(std::vector<ec::Float>& input);
void init_angleTerms();

void fft(std::vector<ec::Float>& inputs, size_t count);

static std::vector<ec::Float> angleTerms(2 * WINDOW_SIZE);

std::vector<ec::Float> process_signal(const std::vector<ec::Float>& inputSignal)
{
    const size_t numSamples = inputSignal.size(); // assume divisible by 32
    const size_t sizeSpectrum = (WINDOW_SIZE / 2) + 1;
    const size_t stepBetweenWins = static_cast<size_t>(ceil(WINDOW_SIZE * (1 - OVERLAP_RATIO)));
    const size_t numWins = (numSamples - WINDOW_SIZE) / stepBetweenWins + 1;

    std::vector<ec::Float> blackmanCoefs(WINDOW_SIZE);
    std::vector<ec::Float> outputSpectrum(sizeSpectrum, std::numeric_limits<float>::lowest());

    size_t idxStartWin = 0;

    init_blackmanCoefs(blackmanCoefs);
    init_angleTerms();

    const ec::Float spC0((float)(10.0f / log(10.0f)));
    const ec::Float spC1((float)(10.0f * log(125.0f / 131072.0f) / log(10.0f)));
    const ec::Float spC2((float)(10.0f * log(125.0f / 32768.0f) / log(10.0f)));

    // ec::VecHw& vecHw = *ec::VecHw::getSingletonVecHw();
    // vecHw.resetMemTo0();

    // vecHw.copyToHw(angleTerms, 0, 2 * WINDOW_SIZE, 0);

    ec::StreamHw& streamHw = *ec::StreamHw::getSingletonStreamHw();
    streamHw.resetStreamHw();

    streamHw.createFifos(41);

    streamHw.copyToHw(angleTerms, 0, 2 * WINDOW_SIZE, 0);

    streamHw.addOpMulToPipeline(0, 1, 2); // co * odd 
    streamHw.addOpMulToPipeline(3, 4, 5); // si * odd
    streamHw.addOpMulToPipeline(6, 7, 8); // si * oddI
    streamHw.addOpMulToPipeline(8, minusOne, 9); // -1 * (si * oddI)

    streamHw.addOpMulToPipeline(10, 11, 12); // co * oddI

    streamHw.addOpAddToPipeline(9, 2, 13); // -1 * (si * oddI) + (co * odd) -> c1
    streamHw.addOpAddToPipeline(12, 5, 14); // (co * oddI) + (si * odd) -> c2

    streamHw.addOpAddToPipeline(15, 13, 16); // even + c1
    streamHw.addOpAddToPipeline(18, 14, 19); // even[halfCount] + c2

    // generate another set of c1 and c2 but NEGATIVE this time!!
    streamHw.addOpMulToPipeline(20, 21, 22); // co * odd 
    streamHw.addOpMulToPipeline(23, 24, 25); // si * odd
    streamHw.addOpMulToPipeline(26, 27, 28); // si * oddI
    streamHw.addOpMulToPipeline(28, minusOne, 29); // -1 * (si * oddI)

    streamHw.addOpMulToPipeline(30, 31, 32); // co * oddI

    streamHw.addOpAddToPipeline(29, 22, 33); // -1 * (si * oddI) + (co * odd) -> c1
    streamHw.addOpAddToPipeline(32, 25, 34); // (co * oddI) + (si * odd) -> c2
    streamHw.addOpMulToPipeline(33, minusOne, 35); // -1 * c1
    streamHw.addOpMulToPipeline(34, minusOne, 36); // -1 * c2

    streamHw.addOpAddToPipeline(37, 35, 38); // even - c1
    streamHw.addOpAddToPipeline(39, 36, 40); // even[halfCount] -c2

    for (size_t j = 0; j < numWins; j++)
    {
        std::vector<ec::Float> signalWindow(2 * WINDOW_SIZE);

        for (size_t i = 0; i < WINDOW_SIZE; i++)
        {
            signalWindow[i] = inputSignal[i + idxStartWin] * blackmanCoefs[i];
            signalWindow[i + WINDOW_SIZE] = ec::Float(0.0f);
        }

        fft(signalWindow, WINDOW_SIZE);

        for (size_t i = 0; i < sizeSpectrum; i++)
        {
            ec::Float freqVal = signalWindow[i] * signalWindow[i] + signalWindow[WINDOW_SIZE + i] * signalWindow[WINDOW_SIZE + i];
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

void fft(std::vector<ec::Float>& inputs, size_t count)
{
    size_t halfCount = count / 2;

    std::vector<ec::Float> even(count);
    std::vector<ec::Float> odd(count);

    if (count > 4)
    {
        for (size_t i = 0; i < halfCount; i++)
        {
            even[i] = inputs[i * 2];
            even[halfCount + i] = inputs[count + i * 2];

            odd[i] = inputs[i * 2 + 1];
            odd[halfCount + i] = inputs[count + i * 2 + 1];
        }

        fft(even, halfCount);
        fft(odd, halfCount);
    }
    else
    {
        even[0] = inputs[0] + inputs[2];
        even[halfCount] = inputs[count + 0] + inputs[count + 2];

        even[1] = inputs[0] - inputs[2];
        even[halfCount + 1] = inputs[count + 0] - inputs[count + 2];

        odd[0] = inputs[1] + inputs[3];
        odd[halfCount] = inputs[count + 1] + inputs[count + 3];

        odd[1] = inputs[1] - inputs[3];
        odd[halfCount + 1] = inputs[count + 1] - inputs[count + 3];

        //TODO try return here
    }

    // negative returns if we run it with too small of data because of copy overhead
    if (count == 1024)
    {
        std::vector<ec::Float> test2(2 * WINDOW_SIZE);
        ec::StreamHw& streamHw = *ec::StreamHw::getSingletonStreamHw();
        streamHw.resetStreamHw();

        streamHw.createFifos(41);

        streamHw.copyToHw(angleTerms, 0, 512, 0);
        streamHw.copyToHw(angleTerms, 1024, 512, 1024);
        streamHw.copyToHw(odd, 0, count, 2 * WINDOW_SIZE);
        streamHw.copyToHw(even, 0, count, 3 * WINDOW_SIZE);

        streamHw.addOpMulToPipeline(0, 1, 2); // co * odd 
        streamHw.addOpMulToPipeline(3, 4, 5); // si * odd
        streamHw.addOpMulToPipeline(6, 7, 8); // si * oddI
        streamHw.addOpMulToPipeline(8, minusOne, 9); // -1 * (si * oddI)

        streamHw.addOpMulToPipeline(10, 11, 12); // co * oddI

        streamHw.addOpAddToPipeline(9, 2, 13); // -1 * (si * oddI) + (co * odd) -> c1
        streamHw.addOpAddToPipeline(12, 5, 14); // (co * oddI) + (si * odd) -> c2

        streamHw.startStreamDataMemToFifo(WINDOW_SIZE - count, 0, halfCount); // cos to pipe 0
        streamHw.startStreamDataMemToFifo(WINDOW_SIZE - count, 10, halfCount); // cos to pipe 10
        streamHw.startStreamDataMemToFifo(2 * WINDOW_SIZE - count, 3, halfCount); // sin to pipe 3
        streamHw.startStreamDataMemToFifo(2 * WINDOW_SIZE - count, 6, halfCount); // sin to pipe 6

        // odd to pipe 1 and 4
        streamHw.startStreamDataMemToFifo(2 * WINDOW_SIZE, 1, halfCount);
        streamHw.startStreamDataMemToFifo(2 * WINDOW_SIZE, 4, halfCount);

        streamHw.startStreamDataMemToFifo(2 * WINDOW_SIZE + halfCount, 7, halfCount); // oddI to pipe 7
        streamHw.startStreamDataMemToFifo(2 * WINDOW_SIZE + halfCount, 11, halfCount); // oddI to pipe 11

        streamHw.startStreamDataFifoToMem(13, 2 * WINDOW_SIZE, halfCount);
        streamHw.startStreamDataFifoToMem(14, 2 * WINDOW_SIZE + halfCount, halfCount);

        streamHw.runPipeline();

        streamHw.startStreamDataMemToFifo(2 * WINDOW_SIZE, 13, halfCount); // c1 -> 13
        streamHw.startStreamDataMemToFifo(2 * WINDOW_SIZE + halfCount, 14, halfCount); // c2 14

        streamHw.addOpAddToPipeline(15, 13, 16); // even + c1
        streamHw.addOpAddToPipeline(18, 14, 19); // even[halfCount] + c2

        streamHw.startStreamDataMemToFifo(3 * WINDOW_SIZE, 15, halfCount); // even -> 15
        streamHw.startStreamDataMemToFifo(3 * WINDOW_SIZE + halfCount, 18, halfCount); // even[halfCount] 18

        streamHw.startStreamDataFifoToMem(16, 0, halfCount); // inputs[k] = even[k] + c1;
        streamHw.startStreamDataFifoToMem(19, count, halfCount); // inputs[count + k] = even[halfCount + k] + c2;

        streamHw.startStreamDataMemToFifo(2 * WINDOW_SIZE, 33, halfCount); // c1 -> 33
        streamHw.startStreamDataMemToFifo(2 * WINDOW_SIZE + halfCount, 34, halfCount); // c2 34

        streamHw.addOpMulToPipeline(33, minusOne, 35); // -1 * c1
        streamHw.addOpMulToPipeline(34, minusOne, 36); // -1 * c2

        streamHw.addOpAddToPipeline(37, 35, 38); // even - c1
        streamHw.addOpAddToPipeline(39, 36, 40); // even[halfCount] -c2

        streamHw.startStreamDataFifoToMem(38, halfCount, halfCount); // inputs[halfCount + k] = even[k] - c1;
        streamHw.startStreamDataFifoToMem(40, count + halfCount, halfCount); // inputs[count + halfCount + k] = even[halfCount + k] - c2;

        streamHw.runPipeline();

        streamHw.copyFromHw(inputs, 0, 2 * count, 0);

        streamHw.resetStreamHw();

        streamHw.createFifos(41);

        streamHw.copyToHw(angleTerms, 512, 2 * WINDOW_SIZE - 512, 512);

        streamHw.addOpMulToPipeline(0, 1, 2); // co * odd 
        streamHw.addOpMulToPipeline(3, 4, 5); // si * odd
        streamHw.addOpMulToPipeline(6, 7, 8); // si * oddI
        streamHw.addOpMulToPipeline(8, minusOne, 9); // -1 * (si * oddI)

        streamHw.addOpMulToPipeline(10, 11, 12); // co * oddI

        streamHw.addOpAddToPipeline(9, 2, 13); // -1 * (si * oddI) + (co * odd) -> c1
        streamHw.addOpAddToPipeline(12, 5, 14); // (co * oddI) + (si * odd) -> c2

        streamHw.addOpAddToPipeline(15, 13, 16); // even + c1
        streamHw.addOpAddToPipeline(18, 14, 19); // even[halfCount] + c2

        // generate another set of c1 and c2 but NEGATIVE this time!!
        streamHw.addOpMulToPipeline(20, 21, 22); // co * odd 
        streamHw.addOpMulToPipeline(23, 24, 25); // si * odd
        streamHw.addOpMulToPipeline(26, 27, 28); // si * oddI
        streamHw.addOpMulToPipeline(28, minusOne, 29); // -1 * (si * oddI)

        streamHw.addOpMulToPipeline(30, 31, 32); // co * oddI

        streamHw.addOpAddToPipeline(29, 22, 33); // -1 * (si * oddI) + (co * odd) -> c1
        streamHw.addOpAddToPipeline(32, 25, 34); // (co * oddI) + (si * odd) -> c2
        streamHw.addOpMulToPipeline(33, minusOne, 35); // -1 * c1
        streamHw.addOpMulToPipeline(34, minusOne, 36); // -1 * c2

        streamHw.addOpAddToPipeline(37, 35, 38); // even - c1
        streamHw.addOpAddToPipeline(39, 36, 40); // even[halfCount] -c2
    }
    else if (halfCount >= 64)
    {
        // std::vector<ec::Float> test(2 * WINDOW_SIZE);
        // std::vector<ec::Float> test2(2 * WINDOW_SIZE);

        ec::StreamHw& streamHw = *ec::StreamHw::getSingletonStreamHw();

        streamHw.copyToHw(odd, 0, count, 2 * WINDOW_SIZE);
        streamHw.copyToHw(even, 0, count, 3 * WINDOW_SIZE);

        streamHw.startStreamDataMemToFifo(WINDOW_SIZE - count, 0, halfCount); // cos to pipe 0
        streamHw.startStreamDataMemToFifo(WINDOW_SIZE - count, 10, halfCount); // cos to pipe 10
        streamHw.startStreamDataMemToFifo(2 * WINDOW_SIZE - count, 3, halfCount); // sin to pipe 3
        streamHw.startStreamDataMemToFifo(2 * WINDOW_SIZE - count, 6, halfCount); // sin to pipe 6

        // odd to pipe 1 and 4
        streamHw.startStreamDataMemToFifo(2 * WINDOW_SIZE, 1, halfCount);
        streamHw.startStreamDataMemToFifo(2 * WINDOW_SIZE, 4, halfCount);

        streamHw.startStreamDataMemToFifo(2 * WINDOW_SIZE + halfCount, 7, halfCount); // oddI to pipe 7
        streamHw.startStreamDataMemToFifo(2 * WINDOW_SIZE + halfCount, 11, halfCount); // oddI to pipe 11

        streamHw.startStreamDataMemToFifo(3 * WINDOW_SIZE, 15, halfCount); // even -> 15
        streamHw.startStreamDataMemToFifo(3 * WINDOW_SIZE + halfCount, 18, halfCount); // even[halfCount] 18

        streamHw.startStreamDataFifoToMem(16, 2 * WINDOW_SIZE, halfCount); // inputs[k] = even[k] + c1;
        streamHw.startStreamDataFifoToMem(19, 2 * WINDOW_SIZE + count, halfCount); // inputs[count + k] = even[halfCount + k] + c2;

        streamHw.startStreamDataMemToFifo(WINDOW_SIZE - count, 20, halfCount); // cos to pipe 20
        streamHw.startStreamDataMemToFifo(WINDOW_SIZE - count, 30, halfCount); // cos to pipe 30
        streamHw.startStreamDataMemToFifo(2 * WINDOW_SIZE - count, 23, halfCount); // sin to pipe 23
        streamHw.startStreamDataMemToFifo(2 * WINDOW_SIZE - count, 26, halfCount); // sin to pipe 26

        // odd to pipe 1 and 4
        streamHw.startStreamDataMemToFifo(2 * WINDOW_SIZE, 21, halfCount);
        streamHw.startStreamDataMemToFifo(2 * WINDOW_SIZE, 24, halfCount);

        streamHw.startStreamDataMemToFifo(2 * WINDOW_SIZE + halfCount, 27, halfCount); // oddI to pipe 27
        streamHw.startStreamDataMemToFifo(2 * WINDOW_SIZE + halfCount, 31, halfCount); // oddI to pipe 31

        streamHw.startStreamDataMemToFifo(3 * WINDOW_SIZE, 37, halfCount); // even -> 37
        streamHw.startStreamDataMemToFifo(3 * WINDOW_SIZE + halfCount, 39, halfCount); // even[halfCount] 39

        streamHw.startStreamDataFifoToMem(38, 2 * WINDOW_SIZE + halfCount, halfCount); // inputs[halfCount + k] = even[k] - c1;
        streamHw.startStreamDataFifoToMem(40, 2 * WINDOW_SIZE + count + halfCount, halfCount); // inputs[count + halfCount + k] = even[halfCount + k] - c2;

        streamHw.runPipeline();

        streamHw.copyFromHw(inputs, 2 * WINDOW_SIZE, 2 * count, 0);
    }
    else
    {
        for (size_t k = 1; k < halfCount; k++)
        {
            ec::Float co = angleTerms[WINDOW_SIZE - count + k];
            ec::Float si = angleTerms[2 * WINDOW_SIZE - count + k];

            ec::Float c1 = (co * odd[k] - si * odd[halfCount + k]);
            ec::Float c2 = (co * odd[halfCount + k] + si * odd[k]);

            inputs[k] = even[k] + c1;
            inputs[count + k] = even[halfCount + k] + c2;

            inputs[halfCount + k] = even[k] - c1;
            inputs[count + halfCount + k] = even[halfCount + k] - c2;
        }
    }

    inputs[0] = even[0] + odd[0];
    inputs[count + 0] = even[halfCount] + odd[halfCount];

    inputs[halfCount] = even[0] - odd[0];
    inputs[count + halfCount] = even[halfCount] - odd[halfCount];
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
                angleTerms[idx] = ec_cos(aC * i);
                angleTerms[idx + WINDOW_SIZE] = ec_sin(aC * i);
            }
            // else
            // {
            //     angleTerms[idx] = ec::Float(1.0f);
            //     angleTerms[idx + WINDOW_SIZE] = ec::Float(0.0f);
            // }

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
