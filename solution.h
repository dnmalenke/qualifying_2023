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

static constexpr float OVERLAP_RATIO = 0.69;
static constexpr size_t WINDOW_SIZE = 1024;
static const ec::Float PI = 3.14159265358979323846f;
static const ec::Float minusOne = -1.0f;

static const ec::Float spC0((float)(10.0f / log(10.0f)));
static const ec::Float spC1((float)(10.0f * log(125.0f / 131072.0f) / log(10.0f)));
static const ec::Float spC2((float)(10.0f * log(125.0f / 32768.0f) / log(10.0f)));

void init_blackmanCoefs();
void init_angleTerms();
uint16_t reverse_bits(uint16_t x);

void fft(std::vector<ec::Float>& inputReal, std::vector<ec::Float>& inputImag, size_t count);

static std::vector<ec::Float> angleTerms(2 * WINDOW_SIZE);
static std::vector<ec::Float> blackmanCoefs(WINDOW_SIZE);

std::vector<ec::Float> process_signal(const std::vector<ec::Float>& inputSignal)
{
    const size_t numSamples = inputSignal.size(); // assume divisible by 32
    const size_t sizeSpectrum = (WINDOW_SIZE / 2) + 1;
    const size_t stepBetweenWins = static_cast<size_t>(ceil(WINDOW_SIZE * (1 - OVERLAP_RATIO)));
    const size_t numWins = (numSamples - WINDOW_SIZE) / stepBetweenWins + 1;

    std::vector<ec::Float> signalWindow(WINDOW_SIZE);

    std::vector<ec::Float> outputSpectrum(sizeSpectrum, std::numeric_limits<float>::lowest());
    std::vector<ec::Float> preLogSpectrum(sizeSpectrum, std::numeric_limits<float>::lowest());

    size_t idxStartWin = 0;

    init_blackmanCoefs();
    init_angleTerms();

    ec::VecHw& vecHw = *ec::VecHw::getSingletonVecHw();
    vecHw.copyToHw(angleTerms, 0, 2 * WINDOW_SIZE - SWEET_SPOT, 2 * WINDOW_SIZE);

    for (size_t j = 0; j < numWins; j++)
    {
        memcpy(signalWindow.data(), inputSignal.data() + idxStartWin, WINDOW_SIZE * sizeof(ec::Float));

        // puts cosine terms back into vecHw. we were using that space for buffering
        // vecHw.copyToHw(angleTerms, 0, 16, 2 * WINDOW_SIZE); // cost 292. we *might* be able to recalculate those cosine terms faster?

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
    ec::VecHw& vecHw = *ec::VecHw::getSingletonVecHw();

    std::vector<ec::Float> temp(2 * WINDOW_SIZE);

    memcpy(temp.data(), inputReal.data(), 93 * sizeof(ec::Float));
    memcpy(temp.data() + 93, blackmanCoefs.data(), 27 * sizeof(ec::Float));

    vecHw.copyToHw(temp, 0, 120, 0);

    // apply the blackman coefficients to input data
    // these are stored the location of the imaginary input data because we don't have any
    for (size_t i = 0; i < 1; i++)
    {
        vecHw.mul32(32 * i, 93 + 32 * i, 32 * i);
    }

    vecHw.resetMemTo0(512, 512); // reset IMAG(0)

    int c = WINDOW_SIZE / 2;

    // c = 512 only working with reals for now
    for (size_t i = 0; i < c / 32 - 13; i++)
    {
        // (-1)*inputReal[C - 2C]
        vecHw.mul32(REAL(c), minusOne, IMAG(c), 32); // stick this in imag[C - 2C] for now. it will be fixed later

        // [0 - C] = [0 - C] + [C - 2C]
        vecHw.add32(REAL(0), REAL(c), IMAG(0)); // store in imag[0] for now.

        // [C - 2C] = [0 - C] + (-1)*[C - 2C]
        vecHw.add32(REAL(0), IMAG(c), REAL(c));

        // now we can assign it to [0 - C] because we're done reading it
        vecHw.assign32(IMAG(0), REAL(0));
    }

    vecHw.resetMemTo0(WINDOW_SIZE, WINDOW_SIZE); // reset IMAG(0)

    /*
    input[C - 2C] * omega[0 through C] # this is complex multiplication
    example of complex multiplication:
    (x + yj) * (a + bj) = (xa - yb) + (xb + ya)j
    So the real part of the product is (xa - yb), and the imaginary part of the product is (xb + ya).
   */
    for (size_t i = 0; i < c / 32 - 13; i++)
    {
        vecHw.mul32(REAL(c), 3 * WINDOW_SIZE + 32 * i, IMAG(c)); // imag [C - 2C] = real[C - 2C] * sin(0 - C)
        vecHw.mul32(REAL(c), 2 * WINDOW_SIZE + 32 * i, REAL(c)); // real [C - 2C] = real [C - 2C] * cos(0 - C)
    }

    // Note we have now used the first 512 values in 2 * WINDOW_SIZE and 3 * WINDOW_SIZE and will no longer need them.
    size_t buf0 = 2 * WINDOW_SIZE;
    size_t buf1 = 2 * WINDOW_SIZE + 32;
    size_t buf2 = 2 * WINDOW_SIZE + 64;
    size_t buf3 = 2 * WINDOW_SIZE + 96;
    size_t buf4 = 2 * WINDOW_SIZE + 128;

    while (c > SWEET_SPOT) // sweet spot
    {
        c /= 2;

        // when the arrays we are working with are larger than 32, our operations need to be chunked out into chunks of 32
        // the variable i tracks which chunk we are currently in
        for (size_t i = 0; i < std::max(1, c / 32) - 1; i++)
        {
            size_t vals = REAL(0);
            size_t valsC = REAL(c);

            for (size_t j = 1; j < WINDOW_SIZE / (2 * c); j++)
            {
                // [2C - 3C] = [2C - 3C] + [3C - 4C]
                vecHw.add32(vals, valsC, vals, std::min(c, 32));

                vals += 2 * c;
                valsC += 2 * c;
            }
        }


        // rearrage step
        if (c == 32)
        {
            c = 8;
            size_t vals = 0;
            size_t valsC = 32;

            size_t realC = 32;
            size_t imagC = WINDOW_SIZE + 32;

            for (size_t i = 0; i < 32 / c - 1; i++)
            {
                // vecHw.assign32(2 * WINDOW_SIZE + (WINDOW_SIZE - 2 * c), buf3 + i * c, c);
                vecHw.assign32(3 * WINDOW_SIZE + (WINDOW_SIZE - 2 * c), buf4 + i * c, c);
            }

            for (size_t j = 0; j < 46; j++)
            {
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

                if (j != 0)
                {
                    // imaginaries
                    // [3C - 4C] = [2C - 3C] + (-1)*[3C - 4C]
                    vecHw.assign32(buf1, valsC);
                }

                // complex multiplication stage:
                if (j != 0)
                {
                    vecHw.mul32(realC, buf4, buf2); // imag [C - 2C] = real[C - 2C] * sin(C-2C)
                    vecHw.add32(imagC, buf2, imagC); // imag [C - 2C] = (imag [C - 2C] * cos(C-2C)) + (real[C - 2C] * sin(C-2C))
                }

                // repeat for 3C, 5C, and 7C
                realC += 32;
                imagC += 32;

                // remove the imaginary val offset and move on to the next pair
                vals -= WINDOW_SIZE;
                valsC -= WINDOW_SIZE;

                vals += 32;
                valsC += 32;
            }

            c /= 2;
            break;
        }
    }

    vecHw.copyFromHw(temp, 0, WINDOW_SIZE - (512 - 479) + 1, 0);

    memcpy(inputReal.data(), temp.data(), 513 * sizeof(ec::Float));
    memcpy(inputImag.data(), temp.data() + 513, 479 * sizeof(ec::Float));

    // returning here increases score cause then the next window overlaps and takes care of it but that requires extra logs

    std::vector<ec::Float> buffer(32);

    while (c > 1)
    {
        c /= 2;

        size_t vals = 0;
        size_t valsC = c;

        // we need to run for real and imaginary
        // [0 - C] = [0 - C] + [C - 2C]
        // [C - 2C] = [0 - C] + (-1)*[C - 2C]
        // [2C - 3C] = [2C- 3C] + [3C - 4C]
        // [3C - 4C] = [2C - 3C] + (-1)*[3C - 4C]

        // [4C - 5C] = [4C - 5C] + [5C - 6C]
        // [5C - 6C] = [4C - 5C] + (-1)*[5C - 6C]
        // [6C - 7C] = [6C - 7C] + [7C - 8C]
        // [7C - 8C] = [6C - 7C] + (-1)*[7C - 8C]
        for (size_t j = 0; j < WINDOW_SIZE / (2 * c) - 256; j++)
        {
            memcpy(inputReal.data() + valsC, buffer.data(), c * sizeof(ec::Float));

            for (size_t k = 0; k < c; k++)
            {
                if (j == 0)
                {
                    memcpy(buffer.data() + k, inputImag.data() + vals + k, sizeof(ec::Float));
                }
                else
                {
                    if (k == 0 && j == 1)
                    {
                        memcpy(inputImag.data() + vals + k, inputImag.data() + valsC + k, sizeof(ec::Float));
                    }
                }
            }

            memcpy(inputImag.data() + valsC, buffer.data(), c * sizeof(ec::Float));

            vals += 2 * c;
            valsC += 2 * c;
        }

        if (c == 1)
        {
            break;
        }

        size_t realC = c;
        size_t imagC = c;

        // example of complex multiplication:
        // (x + yj) * (a + bj) = (xa - yb) + (xb + ya)j
        // So the real part of the product is (xa - yb), and the imaginary part of the product is (xb + ya).
        for (size_t j = 0; j < WINDOW_SIZE / (2 * c); j++)
        {
            for (size_t k = 0; k < c; k++)
            {
                if (k == 0)
                {
                    if (c != 1)
                    {
                        memcpy(buffer.data() + k, inputReal.data() + realC + k, sizeof(ec::Float));
                    }
                }
            }

            if (c != 1)
            {
                memcpy(inputReal.data() + realC, buffer.data(), c * sizeof(ec::Float));
            }

            // repeat for 3C, 5C, and 7C
            realC += 2 * c;
            imagC += 2 * c;
        }
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