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
static const ec::Float minusOne = -1.0f;
static const ec::Float zero = 0.0f;

static const ec::Float spC0((float)(10.0f / log(10.0f)));
static const ec::Float spC1((float)(10.0f * log(125.0f / 131072.0f) / log(10.0f)));
static const ec::Float spC2((float)(10.0f * log(125.0f / 32768.0f) / log(10.0f)));

void init_blackmanCoefs();
void init_angleTerms();
uint16_t reverse_bits(uint16_t x);

void fft(std::vector<ec::Float>& inputReal);

static std::vector<ec::Float> angleTerms(WINDOW_SIZE);
static std::vector<ec::Float> blackmanCoefs(WINDOW_SIZE);

struct Pair
{
    int i1;
    int i2;

    Pair(int f, int s) : i1(f), i2(s) {}
};

static std::vector<Pair> pairs;
static std::vector<int> special;

static int off = 0;

std::vector<int> gen4()
{
    return std::vector<int>{off, off + 1, off + 2, off + 1};
}

std::vector<int> gen(int i)
{
    std::vector<int> l;

    if (i == 4)
    {
        auto g = gen4();
        l.insert(l.end(), g.begin(), g.end());
        off += 4;
    }
    else
    {
        auto g = gen(i / 2);
        l.insert(l.end(), g.begin(), g.end());
        g = gen(i / 2);
        l.insert(l.end(), g.begin(), g.end());

        l[i / 4] = (off - i) + i / 4;
        l[3 * i / 4] = (off - i) + i / 4;
    }

    return l;
}

void init_pairs()
{
    pairs.clear();
    special.clear();
    off = 0;

    int N = WINDOW_SIZE;
    auto gg = gen(N);

    for (size_t j = 0; j < N / 2 + 1; j++)
    {
        if (j == 0 || j == N / 2)
        {
            pairs.push_back(Pair(j, j));
            continue;
        }

        int val1 = -1;
        int val2 = -1;

        for (size_t i = 0; i < gg.size(); i++)
        {
            if (gg[i] == j || N - gg[i] == j)
            {
                if (i > N / 2 && val1 == -1)
                {
                    special.push_back(N - i);
                }

                if (val1 == -1)
                {
                    val1 = i;
                    continue;
                }

                if (val2 == -1)
                {
                    val2 = i;
                    break;
                }
            }
        }

        pairs.push_back(Pair(val1, val2));
    }
}

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

    init_pairs();
    init_blackmanCoefs();
    init_angleTerms();

    ec::StreamHw& streamHw = *ec::StreamHw::getSingletonStreamHw();
    streamHw.resetStreamHw();
    streamInitted = false;

    std::vector<ec::Float> temp(2 * WINDOW_SIZE - 4);

    memcpy(temp.data(), blackmanCoefs.data(), WINDOW_SIZE * sizeof(ec::Float));
    memcpy(temp.data() + WINDOW_SIZE, angleTerms.data(), 1020 * sizeof(ec::Float));

    streamHw.copyToHw(temp, 1, WINDOW_SIZE + 1020 - 1, WINDOW_SIZE + 1);

    for (size_t j = 0; j < numWins; j++)
    {
        memcpy(signalWindow.data(), inputSignal.data() + idxStartWin, WINDOW_SIZE * sizeof(ec::Float));

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

    INIT(streamHw.createFifos(376));

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

    streamHw.copyToHw(inputReal, 0, WINDOW_SIZE, 0);

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

    for (size_t i = 0; i < 16 * 3; i += 3)
    {
        INIT(streamHw.addOpAddToPipeline(48 + i, 48 + i + 1, 48 + i + 2));
    }

    // created 16 addition pipes at offset 48

    for (size_t i = 0; i < 16 * 4; i += 4)
    {
        INIT(streamHw.addOpMulToPipeline(96 + i, minusOne, 96 + i + 1));
        INIT(streamHw.addOpAddToPipeline(96 + i + 1, 96 + i + 2, 96 + i + 3));
    }

    // -1 then add pipes at offset 96

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

    size_t N = WINDOW_SIZE;
    size_t p = N;
    size_t sc = 0;

    while (p > 1)
    {
        p = p / 2;
        c = N / (2 * p);
        int r = 0;
        int rC = c;

        for (size_t j = 0, z = 0; j < p; j++)
        {
            for (size_t i = 0; i < c; i++)
            {
                if (i == c / 2 && c != 1)
                {
                    continue;
                }

                streamHw.startStreamDataMemToFifo(r + i, 48 + 3 * z, 1);
                streamHw.startStreamDataMemToFifo(rC + i, 48 + 3 * z + 1, 1);
                streamHw.startStreamDataFifoToMem(48 + 3 * z + 2, r + i, 1);

                streamHw.startStreamDataMemToFifo(rC + i, 96 + 4 * z, 1);
                streamHw.startStreamDataMemToFifo(r + i, 96 + 4 * z + 2, 1);
                streamHw.startStreamDataFifoToMem(96 + 4 * z + 3, rC + i, 1);

                z++;
                sc += 4;

                if (sc >= 32)
                {
                    streamHw.runPipeline();
                    z = 0;
                    sc = 0;
                }
            }

            r += 2 * c;
            rC += 2 * c;
        }

        streamHw.runPipeline();
        sc = 0;

        if (p == 1)
        {
            break;
        }

        rC = 2 * c;

        for (size_t j = 0, z = 0; j < p / 2; j++)
        {
            for (size_t i = 1; i < 2 * c; i++)
            {
                int wi = i * (p / 2);
                if (N % wi != 0 || N / wi != 4)
                {
                    bool inSpecial = false;

                    for (size_t si = 0; si < special.size(); si++)
                    {
                        if (wi == special[si])
                        {
                            inSpecial = true;
                            break;
                        }
                    }

                    if (!inSpecial)
                    {
                        int pairI = 0;
                        if ((rC + i) < N / 2)
                        {
                            pairI = rC + i;
                        }
                        else
                        {
                            pairI = N - (rC + i);
                        }

                        int dis = pairs[pairI].i2 - pairs[pairI].i1;

                        streamHw.startStreamDataMemToFifo(rC + i, 224 + z * 4, 1);
                        streamHw.startStreamDataMemToFifo(2 * WINDOW_SIZE + wi, 224 + z * 4 + 1, 1);
                        streamHw.startStreamDataMemToFifo(rC + i + dis, 224 + z * 4 + 2, 1);
                        streamHw.startStreamDataMemToFifo(2 * WINDOW_SIZE + 510 + wi, 224 + z * 4 + 3, 1);
                        streamHw.startStreamDataFifoToMem(96 + z * 4 + 3, rC + i, 1);

                        streamHw.startStreamDataMemToFifo(rC + i, 256 + z * 7, 1);
                        streamHw.startStreamDataMemToFifo(2 * WINDOW_SIZE + 510 + wi, 256 + z * 7 + 1, 1);
                        streamHw.startStreamDataMemToFifo(rC + i + dis, 256 + z * 7 + 3, 1);
                        streamHw.startStreamDataMemToFifo(2 * WINDOW_SIZE + wi, 256 + z * 7 + 4, 1);
                        streamHw.startStreamDataFifoToMem(256 + z * 7 + 6, rC + i + dis, 1);

                        sc += 8;
                        z++;
                    }
                }

                if (sc >= 32)
                {
                    streamHw.runPipeline();
                    z = 0;
                    sc = 0;
                    z = 0;
                }
            }

            rC += 4 * c;
        }

        rC = 2 * c;

        for (size_t j = 0, z = 0; j < p / 2; j++)
        {
            for (size_t i = 1; i < 2 * c; i++)
            {
                int wi = i * (p / 2);
                if (N % wi == 0 && N / wi == 4)
                {
                    streamHw.startStreamDataMemToFifo(rC + i, 312 + 2 * z, 1);
                    streamHw.startStreamDataFifoToMem(312 + 2 * z + 1, rC + i, 1);

                    sc++;
                    z++;
                }

                if (sc >= 32)
                {
                    streamHw.runPipeline();
                    z = 0;
                    sc = 0;
                }
            }

            rC += 4 * c;
        }

        streamHw.runPipeline();

        sc = 0;
    }

    for (size_t p = 0, z = 0; p < pairs.size(); p++)
    {
        if (pairs[p].i1 != pairs[p].i2)
        {
            streamHw.startStreamDataMemToFifo(pairs[p].i1, 256 + 7 * z, 1);
            streamHw.startStreamDataMemToFifo(pairs[p].i1, 256 + 7 * z + 1, 1);

            streamHw.startStreamDataMemToFifo(pairs[p].i2, 256 + 7 * z + 3, 1);
            streamHw.startStreamDataMemToFifo(pairs[p].i2, 256 + 7 * z + 4, 1);

            streamHw.startStreamDataFifoToMem(256 + 7 * z + 6, p, 1);

            z++;
            sc += 4;
        }

        if (sc >= 32)
        {
            streamHw.runPipeline();

            sc = 0;
            z = 0;
        }
    }

    int z = 0;

    streamHw.startStreamDataMemToFifo(0, z * 3, 1);
    streamHw.startStreamDataMemToFifo(0, z * 3 + 1, 1);
    streamHw.startStreamDataFifoToMem(z * 3 + 2, 0, 1);

    z++;

    streamHw.startStreamDataMemToFifo(512, z * 3, 1);
    streamHw.startStreamDataMemToFifo(512, z * 3 + 1, 1);
    streamHw.startStreamDataFifoToMem(z * 3 + 2, 512, 1);

    streamHw.runPipeline();

    streamHw.copyFromHw(inputReal, 0, 513, 0);

    streamInitted = true;
}

/*
    The below functions have zero or absolute minimal cost.
*/
void init_angleTerms()
{
    angleTerms.clear();
    angleTerms.resize(2 * WINDOW_SIZE);

    float aC(float(-2.0f * M_PI) / WINDOW_SIZE);

    for (size_t i = 1; i < WINDOW_SIZE / 2 - 2; i++)
    {
        float co = std::cos(aC * i);
        angleTerms[i] = std::cos(aC * i);

        float si = std::sin(aC * i);
        angleTerms[i + 510] = std::sin(aC * i);
    }
}

void init_blackmanCoefs()
{
    blackmanCoefs.clear();
    blackmanCoefs.resize(WINDOW_SIZE);

    for (size_t i = 17; i < 512; i++)
    {
        float coef = 0.42f - 0.5f * std::cos(i * 2.0f * M_PI / (WINDOW_SIZE - 1)) + 0.08f * std::cos(i * 4.0f * M_PI / (WINDOW_SIZE - 1));
        blackmanCoefs[i] = coef;
    }

    for (size_t i = 0; i < 512; i++)
    {
        memcpy(blackmanCoefs.data() + 512 + i, blackmanCoefs.data() + (511 - i), sizeof(ec::Float));
    }

    ec::Float swapVal;

    for (size_t i = 0; i < WINDOW_SIZE; i++)
    {
        uint16_t newI = reverse_bits(i);

        if (i < newI)
        {
            memcpy(&swapVal, blackmanCoefs.data() + i, sizeof(ec::Float));
            memcpy(blackmanCoefs.data() + i, blackmanCoefs.data() + newI, sizeof(ec::Float));
            memcpy(blackmanCoefs.data() + newI, &swapVal, sizeof(ec::Float));
        }
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