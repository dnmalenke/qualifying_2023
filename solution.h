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

void fft(std::vector<ec::Float>& inputReal);

static std::vector<ec::Float> angleTerms(2 * WINDOW_SIZE);
static std::vector<ec::Float> blackmanCoefs(WINDOW_SIZE);

std::vector<ec::Float> process_signal(const std::vector<ec::Float>& inputSignal)
{
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

    ec::StreamHw& streamHw = *ec::StreamHw::getSingletonStreamHw();
    streamHw.copyToHw(angleTerms, 0, 2 * WINDOW_SIZE, 2 * WINDOW_SIZE);

    for (size_t j = 0; j < numWins; j++)
    {
        memcpy(signalWindow.data(), inputSignal.data() + idxStartWin, WINDOW_SIZE * sizeof(ec::Float));

        // puts cosine terms back into streamHw. we were using that space for buffering
        streamHw.copyToHw(angleTerms, 0, 256, 2 * WINDOW_SIZE);

        fft(signalWindow);

        for (size_t i = 0; i < sizeSpectrum; i++)
        {
            ec::Float freqVal;
            memcpy(&freqVal, signalWindow.data() + i, sizeof(ec::Float));

            // we will always take the first window
            if (j != 0 && freqVal <= preLogSpectrum[i])
            {
                continue;
            }

            memcpy(preLogSpectrum.data() + i, &freqVal, sizeof(ec::Float));
        }

        idxStartWin += stepBetweenWins;
    }

    for (size_t i = 0; i < sizeSpectrum; i++)
    {
        ec::Float freqVal;

        memcpy(&freqVal, preLogSpectrum.data() + i, sizeof(ec::Float));

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

    return outputSpectrum;
}

void fft(std::vector<ec::Float>& inputReal)
{
    ec::StreamHw& streamHw = *ec::StreamHw::getSingletonStreamHw();

    INIT(streamHw.createFifos(1000));

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

    streamHw.runPipeline();

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

    streamHw.runPipeline();

    for (size_t i = 0; i < 16 * 4; i += 4)
    {
        INIT(streamHw.addOpMulToPipeline(96 + i, minusOne, 96 + i + 1));
        INIT(streamHw.addOpAddToPipeline(96 + i + 1, 96 + i + 2, 96 + i + 3));
        streamHw.startStreamDataMemToFifo(p + c * i / (4 * d), 96 + i, c / d);
        streamHw.startStreamDataMemToFifo(c * i / (4 * d), 96 + i + 2, c / d);
        streamHw.startStreamDataFifoToMem(96 + i + 3, p + c * i / (4 * d), c / d);
    }

    streamHw.runPipeline();

    // -1 then add pipes at offset 96

    // move the real vals back to right location
    for (size_t i = 0; i < 32 * 2; i += 2)
    {
        INIT(streamHw.addOpAddToPipeline(160 + i, zero, 160 + i + 1));

        streamHw.startStreamDataMemToFifo(WINDOW_SIZE + c * i / (4 * d), 160 + i, c / (2 * d));
        streamHw.startStreamDataFifoToMem(160 + i + 1, c * i / (4 * d), c / (2 * d));
    }

    streamHw.runPipeline();

    streamHw.resetMemTo0(WINDOW_SIZE, WINDOW_SIZE); // clear out imaginary terms


    for (size_t i = 0; i < 16 * 3; i += 3)
    {
        streamHw.startStreamDataMemToFifo(p + c * i / (3 * d), i, c / d);
        streamHw.startStreamDataMemToFifo(3 * WINDOW_SIZE + c * i / (3 * d), i + 1, c / d);
        streamHw.startStreamDataFifoToMem(i + 2, WINDOW_SIZE + p + c * i / (3 * d), c / d);
    }

    streamHw.runPipeline();

    for (size_t i = 0; i < 16 * 3; i += 3)
    {
        streamHw.startStreamDataMemToFifo(p + c * i / (3 * d), i, c / d);
        streamHw.startStreamDataMemToFifo(2 * WINDOW_SIZE + c * i / 6, i + 1, c / d);
        streamHw.startStreamDataFifoToMem(i + 2, p + c * i / (3 * d), c / d);
    }

    streamHw.runPipeline();

    // now we can use first 256 slots in cos as buffer

    // we need to make the pipelines needed for our complex multiplication starting at offset 224
    for (size_t i = 0; i < 8 * 4; i += 4)
    {
        // real * cos
        INIT(streamHw.addOpMulToPipeline(224 + i, 224 + i + 1, 96 + i + 2));

        // imag * sin * -1
        INIT(streamHw.addOpMulToPipeline(224 + i + 2, 224 + i + 3, 96 + i));
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

    for (size_t i = 0; i < 32 * 2; i += 2)
    {
        INIT(streamHw.addOpMulToPipeline(312 + i, minusOne, 312 + i + 1));
    }
    // add -1 pipes so we can do the -j

    while (p > 32)
    {
        p /= 2;
        d *= 2;

        size_t vals = WINDOW_SIZE - (2 * p);
        size_t valsC = WINDOW_SIZE - p;
        size_t realC = WINDOW_SIZE - p;
        size_t imagC = 2 * WINDOW_SIZE - p;

        for (int j = WINDOW_SIZE / (2 * p) - 1; j >= 0; j--)
        {
            for (size_t i = 0; i < 16 * 4; i += 4)
            {
                streamHw.startStreamDataMemToFifo(valsC + c * i / (4 * d), 96 + i, c / d);
                streamHw.startStreamDataMemToFifo(vals + c * i / (4 * d), 96 + i + 2, c / d);
                streamHw.startStreamDataFifoToMem(96 + i + 3, 2 * WINDOW_SIZE + c * i / (4 * d), c / d);
            }

            streamHw.runPipeline();


            for (size_t i = 0; i < 16 * 3; i += 3)
            {
                streamHw.startStreamDataMemToFifo(vals + c * i / (3 * d), 48 + i, c / d);
                streamHw.startStreamDataMemToFifo(valsC + c * i / (3 * d), 48 + i + 1, c / d);
                streamHw.startStreamDataFifoToMem(48 + i + 2, vals + c * i / (3 * d), c / d);
            }

            streamHw.runPipeline();

            for (size_t i = 0; i < 32 * 2; i += 2)
            {
                streamHw.startStreamDataMemToFifo(2 * WINDOW_SIZE + c * i / (4 * d), 160 + i, c / (2 * d));
                streamHw.startStreamDataFifoToMem(160 + i + 1, valsC + c * i / (4 * d), c / (2 * d));
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
            vals -= 2 * p;

            valsC -= WINDOW_SIZE;
            valsC -= 2 * p;

            if (j == 0)
            {
                for (size_t i = 0; i < 16 * 3; i += 3)
                {
                    streamHw.startStreamDataMemToFifo(realC + c * i / (3 * d), i, c / d);
                    streamHw.startStreamDataMemToFifo(3 * WINDOW_SIZE + (WINDOW_SIZE - 2 * p) + c * i / (3 * d), i + 1, c / d);
                    streamHw.startStreamDataFifoToMem(i + 2, imagC + c * i / (3 * d), c / d);
                }

                streamHw.runPipeline();

                for (size_t i = 0; i < 16 * 3; i += 3)
                {
                    streamHw.startStreamDataMemToFifo(realC + c * i / (3 * d), i, c / d);
                    streamHw.startStreamDataMemToFifo(2 * WINDOW_SIZE + (WINDOW_SIZE - 2 * p) + c * i / (3 * d), i + 1, c / d);
                    streamHw.startStreamDataFifoToMem(i + 2, realC + c * i / (3 * d), c / d);
                }

                streamHw.runPipeline();
            }
            else if (j == 1 && false)
            {
                for (size_t i = 0; i < 8 * 4; i += 4)
                {
                    streamHw.startStreamDataMemToFifo(realC + c * i / (4 * (d / 2)), 224 + i, 2 * c / d);
                    streamHw.startStreamDataMemToFifo(2 * WINDOW_SIZE + (WINDOW_SIZE - 4 * p) + c * i / (4 * (d / 2)), 224 + i + 1, 2 * c / d);

                    // imag * sin * -1
                    streamHw.startStreamDataMemToFifo(imagC + c * i / (4 * (d / 2)), 224 + i + 2, 2 * c / d);
                    streamHw.startStreamDataMemToFifo(3 * WINDOW_SIZE + (WINDOW_SIZE - 4 * p) + c * i / (4 * (d / 2)), 224 + i + 3, 2 * c / d);

                    streamHw.startStreamDataFifoToMem(96 + i + 3, 2 * WINDOW_SIZE + c * i / (4 * (d / 2)), 2 * c / d);
                }
                // real terms to 224 + i
                // cos terms to 224 + i + 1
                // imag terms 224 + i + 2
                // sin terms 224 + i + 3
                // output 96 + i + 3

                streamHw.runPipeline();


                for (size_t i = 0; i < 8 * 7; i += 7)
                {
                    // real * sin

                    streamHw.startStreamDataMemToFifo(realC + c * i / (7 * (d / 2)), 256 + i, 2 * c / d);
                    streamHw.startStreamDataMemToFifo(3 * WINDOW_SIZE + (WINDOW_SIZE - 4 * p) + c * i / (7 * (d / 2)), 256 + i + 1, 2 * c / d);

                    // imag * cos

                    streamHw.startStreamDataMemToFifo(imagC + c * i / (7 * (d / 2)), 256 + i + 3, 2 * c / d);
                    streamHw.startStreamDataMemToFifo(2 * WINDOW_SIZE + (WINDOW_SIZE - 4 * p) + c * i / (7 * (d / 2)), 256 + i + 4, 2 * c / d);

                    // add the two

                    streamHw.startStreamDataFifoToMem(256 + i + 6, imagC + c * i / (7 * (d / 2)), 2 * c / d);
                }
                // real terms to 256 + i
                // sin terms to 256 + i + 1
                // imag terms 256 + i + 3
                // cos terms 256 + i + 4
                // output 256 + i + 6

                streamHw.runPipeline();


                for (size_t i = 0; i < 32 * 2; i += 2)
                {
                    streamHw.startStreamDataMemToFifo(2 * WINDOW_SIZE + c * i / (4 * d), 160 + i, c / (2 * d));
                    streamHw.startStreamDataFifoToMem(160 + i + 1, realC + c * i / (4 * d), c / (2 * d));
                }

                streamHw.runPipeline();
            }
            else
            {
                for (size_t i = 0; i < 8 * 4; i += 4)
                {
                    streamHw.startStreamDataMemToFifo(realC + c * i / (4 * (d / 2)), 224 + i, 2 * c / d);
                    streamHw.startStreamDataMemToFifo(2 * WINDOW_SIZE + (WINDOW_SIZE - 2 * p) + c * i / (4 * (d / 2)), 224 + i + 1, 2 * c / d);

                    // imag * sin * -1
                    streamHw.startStreamDataMemToFifo(imagC + c * i / (4 * (d / 2)), 224 + i + 2, 2 * c / d);
                    streamHw.startStreamDataMemToFifo(3 * WINDOW_SIZE + (WINDOW_SIZE - 2 * p) + c * i / (4 * (d / 2)), 224 + i + 3, 2 * c / d);

                    streamHw.startStreamDataFifoToMem(96 + i + 3, 2 * WINDOW_SIZE + c * i / (4 * (d / 2)), 2 * c / d);
                }
                // real terms to 224 + i
                // cos terms to 224 + i + 1
                // imag terms 224 + i + 2
                // sin terms 224 + i + 3
                // output 96 + i + 3

                streamHw.runPipeline();


                for (size_t i = 0; i < 8 * 7; i += 7)
                {
                    // real * sin

                    streamHw.startStreamDataMemToFifo(realC + c * i / (7 * (d / 2)), 256 + i, 2 * c / d);
                    streamHw.startStreamDataMemToFifo(3 * WINDOW_SIZE + (WINDOW_SIZE - 2 * p) + c * i / (7 * (d / 2)), 256 + i + 1, 2 * c / d);

                    // imag * cos

                    streamHw.startStreamDataMemToFifo(imagC + c * i / (7 * (d / 2)), 256 + i + 3, 2 * c / d);
                    streamHw.startStreamDataMemToFifo(2 * WINDOW_SIZE + (WINDOW_SIZE - 2 * p) + c * i / (7 * (d / 2)), 256 + i + 4, 2 * c / d);

                    // add the two

                    streamHw.startStreamDataFifoToMem(256 + i + 6, imagC + c * i / (7 * (d / 2)), 2 * c / d);
                }
                // real terms to 256 + i
                // sin terms to 256 + i + 1
                // imag terms 256 + i + 3
                // cos terms 256 + i + 4
                // output 256 + i + 6

                streamHw.runPipeline();


                for (size_t i = 0; i < 32 * 2; i += 2)
                {
                    streamHw.startStreamDataMemToFifo(2 * WINDOW_SIZE + c * i / (4 * d), 160 + i, c / (2 * d));
                    streamHw.startStreamDataFifoToMem(160 + i + 1, realC + c * i / (4 * d), c / (2 * d));
                }

                streamHw.runPipeline();
            }

            // repeat for 3C, 5C, and 7C
            realC -= 2 * p;
            imagC -= 2 * p;
        }
    }

    size_t x = 1;

    while (p > 4)
    {
        p /= 2;
        x *= 2;

        size_t vals = 0;
        size_t  valsC = p;

        size_t realC = p;
        size_t imagC = WINDOW_SIZE + p;

        for (size_t j = 0; j < WINDOW_SIZE / (2 * p); j += x)
        {
            for (size_t i = 0, vOff = 0, iOff = 0; i < 16 * 3; i += 3)
            {
                if (i != 0 && i % (16 * 3 / x) == 0)
                {
                    vOff += 2 * p;
                    iOff = i;
                }

                size_t sIdx = vOff + c * (i - iOff) / (3 * d);

                streamHw.startStreamDataMemToFifo(vals + sIdx, 48 + i, c / d);
                streamHw.startStreamDataMemToFifo(valsC + sIdx, 48 + i + 1, c / d);
                streamHw.startStreamDataFifoToMem(48 + i + 2, 2 * WINDOW_SIZE + c * i / (3 * d), c / d);
            }

            streamHw.runPipeline();

            for (size_t i = 0, vOff = 0, iOff = 0; i < 16 * 4; i += 4)
            {
                if (i != 0 && i % (16 * 4 / x) == 0)
                {
                    vOff += 2 * p;
                    iOff = i;
                }

                size_t sIdx = vOff + c * (i - iOff) / (4 * d);

                streamHw.startStreamDataMemToFifo(valsC + sIdx, 96 + i, c / d);
                streamHw.startStreamDataMemToFifo(vals + sIdx, 96 + i + 2, c / d);
                streamHw.startStreamDataFifoToMem(96 + i + 3, valsC + sIdx, c / d);
            }

            streamHw.runPipeline();

            for (size_t i = 0, vOff = 0, iOff = 0; i < 32 * 2; i += 2)
            {
                if (i != 0 && i % (32 * 2 / x) == 0)
                {
                    vOff += 2 * p;
                    iOff = i;
                }

                size_t sIdx = vOff + c * (i - iOff) / (4 * d);

                streamHw.startStreamDataMemToFifo(2 * WINDOW_SIZE + c * i / (4 * d), 160 + i, c / (2 * d));
                streamHw.startStreamDataFifoToMem(160 + i + 1, vals + sIdx, c / (2 * d));
            }

            streamHw.runPipeline();

            // offset to imaginary numbers and do it again
            vals += WINDOW_SIZE;
            valsC += WINDOW_SIZE;


            // imaginaries
            for (size_t i = 0, vOff = 0, iOff = 0; i < 16 * 3; i += 3)
            {
                if (i != 0 && i % (16 * 3 / x) == 0)
                {
                    vOff += 2 * p;
                    iOff = i;
                }

                size_t sIdx = vOff + c * (i - iOff) / (3 * d);

                streamHw.startStreamDataMemToFifo(vals + sIdx, 48 + i, c / d);
                streamHw.startStreamDataMemToFifo(valsC + sIdx, 48 + i + 1, c / d);
                streamHw.startStreamDataFifoToMem(48 + i + 2, 2 * WINDOW_SIZE + c * i / (3 * d), c / d);
            }

            streamHw.runPipeline();

            for (size_t i = 0, vOff = 0, iOff = 0; i < 16 * 4; i += 4)
            {
                if (i != 0 && i % (16 * 4 / x) == 0)
                {
                    vOff += 2 * p;
                    iOff = i;
                }

                size_t sIdx = vOff + c * (i - iOff) / (4 * d);

                streamHw.startStreamDataMemToFifo(valsC + sIdx, 96 + i, c / d);
                streamHw.startStreamDataMemToFifo(vals + sIdx, 96 + i + 2, c / d);
                streamHw.startStreamDataFifoToMem(96 + i + 3, valsC + sIdx, c / d);
            }

            streamHw.runPipeline();

            for (size_t i = 0, vOff = 0, iOff = 0; i < 32 * 2; i += 2)
            {
                if (i != 0 && i % (32 * 2 / x) == 0)
                {
                    vOff += 2 * p;
                    iOff = i;
                }

                size_t sIdx = vOff + c * (i - iOff) / (4 * d);

                streamHw.startStreamDataMemToFifo(2 * WINDOW_SIZE + c * i / (4 * d), 160 + i, c / (2 * d));
                streamHw.startStreamDataFifoToMem(160 + i + 1, vals + sIdx, c / (2 * d));
            }

            streamHw.runPipeline();

            // remove the imaginary val offset and move on to the next pair
            vals -= WINDOW_SIZE;
            vals += 64;

            valsC -= WINDOW_SIZE;
            valsC += 64;


            // complex multiplication
            for (size_t i = 0, vOff = 0, iOff = 0; i < 8 * 4; i += 4)
            {
                if (i != 0 && i % (8 * 4 / x) == 0)
                {
                    vOff += 2 * p;
                    iOff = i;
                }

                size_t sIdx = c * (i - iOff) / (4 * (d / 2));

                streamHw.startStreamDataMemToFifo(vOff + realC + sIdx, 224 + i, 2 * c / d);
                streamHw.startStreamDataMemToFifo(2 * WINDOW_SIZE + (WINDOW_SIZE - 2 * p) + sIdx, 224 + i + 1, 2 * c / d);

                // imag * sin * -1
                streamHw.startStreamDataMemToFifo(vOff + imagC + sIdx, 224 + i + 2, 2 * c / d);
                streamHw.startStreamDataMemToFifo(3 * WINDOW_SIZE + (WINDOW_SIZE - 2 * p) + sIdx, 224 + i + 3, 2 * c / d);
                streamHw.startStreamDataFifoToMem(96 + i + 3, 2 * WINDOW_SIZE + c * i / (4 * (d / 2)), 2 * c / d);
            }
            // real terms to 224 + i
            // cos terms to 224 + i + 1
            // imag terms 224 + i + 2
            // sin terms 224 + i + 3
            // output 96 + i + 3

            streamHw.runPipeline();


            for (size_t i = 0, vOff = 0, iOff = 0; i < 8 * 7; i += 7)
            {
                if (i != 0 && i % (8 * 7 / x) == 0)
                {
                    vOff += 2 * p;
                    iOff = i;
                }

                size_t sIdx = c * (i - iOff) / (7 * (d / 2));

                streamHw.startStreamDataMemToFifo(vOff + realC + sIdx, 256 + i, 2 * c / d);
                streamHw.startStreamDataMemToFifo(3 * WINDOW_SIZE + (WINDOW_SIZE - 2 * p) + sIdx, 256 + i + 1, 2 * c / d);

                streamHw.startStreamDataMemToFifo(vOff + imagC + sIdx, 256 + i + 3, 2 * c / d);
                streamHw.startStreamDataMemToFifo(2 * WINDOW_SIZE + (WINDOW_SIZE - 2 * p) + sIdx, 256 + i + 4, 2 * c / d);

                streamHw.startStreamDataFifoToMem(256 + i + 6, vOff + imagC + sIdx, 2 * c / d);
            }
            // real terms to 256 + i
            // sin terms to 256 + i + 1
            // imag terms 256 + i + 3
            // cos terms 256 + i + 4
            // output 256 + i + 6

            streamHw.runPipeline();


            for (size_t i = 0, vOff = 0, iOff = 0; i < 32 * 2; i += 2)
            {
                if (i != 0 && i % (32 * 2 / x) == 0)
                {
                    vOff += 2 * p;
                    iOff = i;
                }

                size_t sIdx = vOff + c * (i - iOff) / (4 * d);

                streamHw.startStreamDataMemToFifo(2 * WINDOW_SIZE + c * i / (4 * d), 160 + i, c / (2 * d));
                streamHw.startStreamDataFifoToMem(160 + i + 1, realC + sIdx, c / (2 * d));
            }

            streamHw.runPipeline();

            // repeat for 3C, 5C, and 7C
            realC += 64;
            imagC += 64;
        }
    }

    p /= 2; // 2
    x *= 2; // 16
    // d = 32
    // c = 64

    size_t vals = 0;
    size_t valsC = p;
    size_t vvOff = 0;

    size_t realC = p;
    size_t imagC = WINDOW_SIZE + p;
    size_t sc = 0;

    // run the radix 4 for both real and imaginary
    for (size_t j = 0, i = 0; j < WINDOW_SIZE; j++)
    {
        if ((valsC + vvOff) % 4 != 2 && (valsC + vvOff) % 4 != 3 || (valsC + vvOff) < 6)
        {
            streamHw.startStreamDataMemToFifo(valsC + vvOff, 96 + 4 * i, 1);
            streamHw.startStreamDataMemToFifo(vals + vvOff, 96 + 4 * i + 2, 1);
            streamHw.startStreamDataFifoToMem(96 + 4 * i + 3, valsC + vvOff, 1);
        }

        streamHw.startStreamDataMemToFifo(vals + vvOff, 48 + 3 * i, 1);
        streamHw.startStreamDataMemToFifo(valsC + vvOff, 48 + 3 * i + 1, 1);
        streamHw.startStreamDataFifoToMem(48 + 3 * i + 2, vals + vvOff, 1);

        if (vals >= WINDOW_SIZE)
        {
            vals -= WINDOW_SIZE;
            valsC -= WINDOW_SIZE;

            if (vvOff % 2 == 0)
            {
                vvOff += 1;
            }
            else
            {
                vvOff += 2 * p - 1;
            }
        }
        else
        {
            vals += WINDOW_SIZE;
            valsC += WINDOW_SIZE;
        }

        i++;
        sc += 4;

        if (sc == 32)
        {
            streamHw.runPipeline();
            sc = 0;
            i = 0;
        }
    }

    streamHw.runPipeline();

    // swap real and imaginary at element 3
    streamHw.startStreamDataMemToFifo(3, 160, 1);
    streamHw.startStreamDataMemToFifo(WINDOW_SIZE + 3, 160 + 2, 1);
    streamHw.startStreamDataFifoToMem(160 + 3, 3, 1);
    streamHw.startStreamDataFifoToMem(160 + 1, WINDOW_SIZE + 3, 1);
    streamHw.runPipeline();

    sc = 0;

    for (size_t j = 0, i = 0; j < WINDOW_SIZE; j += 2)
    {
        size_t newI = reverse_bits(j);

        if (newI > 768 || (newI > 256 && newI < 512))
        {
            continue;
        }

        streamHw.startStreamDataMemToFifo(j, 48 + 3 * i, 1);
        streamHw.startStreamDataMemToFifo(j + 1, 48 + 3 * i + 1, 1);
        streamHw.startStreamDataFifoToMem(48 + 3 * i + 2, j, 1);

        streamHw.startStreamDataMemToFifo(j, 96 + 4 * i + 2, 1);
        streamHw.startStreamDataMemToFifo(j + 1, 96 + 4 * i, 1);
        streamHw.startStreamDataFifoToMem(96 + 4 * i + 3, j + 1, 1);

        i++;
        sc += 4;

        if (sc == 32)
        {
            streamHw.runPipeline();

            sc = 0;
            i = 0;
        }
    }

    for (size_t j = WINDOW_SIZE, i = 0; j < 2 * WINDOW_SIZE; j += 2)
    {
        size_t newI = reverse_bits(j - WINDOW_SIZE);

        if (newI > 768 || (newI > 256 && newI < 512))
        {
            continue;
        }

        streamHw.startStreamDataMemToFifo(j, 48 + 3 * i, 1);
        streamHw.startStreamDataMemToFifo(j + 1, 48 + 3 * i + 1, 1);
        streamHw.startStreamDataFifoToMem(48 + 3 * i + 2, j, 1);

        streamHw.startStreamDataMemToFifo(j, 96 + 4 * i + 2, 1);
        streamHw.startStreamDataMemToFifo(j + 1, 96 + 4 * i, 1);
        streamHw.startStreamDataFifoToMem(96 + 4 * i + 3, j + 1, 1);

        i++;
        sc += 4;

        if (sc == 32)
        {
            streamHw.runPipeline();

            sc = 0;
            i = 0;
        }
    }

    for (size_t j = 0, i = 0; j < WINDOW_SIZE; j++)
    {
        size_t newI = reverse_bits(j);

        if (newI > 768 || (newI > 256 && newI < 512))
        {
            continue;
        }

        streamHw.startStreamDataMemToFifo(j, 256 + 7 * i, 1);
        streamHw.startStreamDataMemToFifo(j, 256 + 7 * i + 1, 1);

        streamHw.startStreamDataMemToFifo(WINDOW_SIZE + j, 256 + 7 * i + 3, 1);
        streamHw.startStreamDataMemToFifo(WINDOW_SIZE + j, 256 + 7 * i + 4, 1);

        streamHw.startStreamDataFifoToMem(256 + 7 * i + 6, j, 1);

        i++;
        sc += 4;
        if (sc == 32)
        {
            streamHw.runPipeline();

            sc = 0;
            i = 0;
        }
    }

    streamHw.runPipeline();

    streamHw.copyFromHw(inputReal, 0, WINDOW_SIZE, 0);

    streamInitted = true;

    ec::Float swapVal;

    for (size_t i = 0; i < WINDOW_SIZE; i++)
    {
        uint16_t newI = reverse_bits(i);

        if (i < newI)
        {
            memcpy(&swapVal, inputReal.data() + i, sizeof(ec::Float));
            memcpy(inputReal.data() + i, inputReal.data() + newI, sizeof(ec::Float));
            memcpy(inputReal.data() + newI, &swapVal, sizeof(ec::Float));
        }
    }

    for (size_t i = 512; i < 768; i++)
    {
        size_t newI = 512 - (i - 512);

        memcpy(&swapVal, inputReal.data() + i, sizeof(ec::Float));
        memcpy(inputReal.data() + i, inputReal.data() + newI, sizeof(ec::Float));
        memcpy(inputReal.data() + newI, &swapVal, sizeof(ec::Float));
    }
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

    // for (size_t i = WINDOW_SIZE / 4; i < WINDOW_SIZE / 2; i++)
    // {
    //     angleTerms[i] = 0.0f;

    //     angleTerms[i + WINDOW_SIZE] = 0.0f;
    // }


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
    for (size_t i = 0; i < 512; i++)
    {
        float coef = 0.42f - 0.5f * std::cos(i * 2.0f * M_PI / (WINDOW_SIZE - 1)) + 0.08f * std::cos(i * 4.0f * M_PI / (WINDOW_SIZE - 1));
        blackmanCoefs[i] = coef;
    }

    for (size_t i = 0; i < 512; i++)
    {
        memcpy(blackmanCoefs.data() + 512 + i, blackmanCoefs.data() + (511 - i), sizeof(ec::Float));
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