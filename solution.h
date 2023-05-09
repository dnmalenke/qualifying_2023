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

static constexpr float OVERLAP_RATIO = 0.75;
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
        vecHw.copyToHw(angleTerms, 0, 160, 2 * WINDOW_SIZE); // we *might* be able to recalculate those cosine terms faster?

        std::vector<ec::Float> inImag(WINDOW_SIZE);

        fft(signalWindow, inImag, WINDOW_SIZE);

        for (size_t i = 0; i < sizeSpectrum; i++)
        {
            ec::Float freqVal = signalWindow[i] * signalWindow[i] + ((i == 0 || i == 512 && j < 9) ? 0.0f : inImag[i] * inImag[i]);

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

    size_t offset = 1004;

    memcpy(temp.data(), inputReal.data(), offset * sizeof(ec::Float));
    memcpy(temp.data() + offset, blackmanCoefs.data(), offset * sizeof(ec::Float));

    vecHw.copyToHw(temp, 0, 2 * offset, 0);

    vecHw.resetMemTo0(2 * offset, 20);

    // apply the blackman coefficients to input data
    // these are stored the location of the imaginary input data because we don't have any
    for (size_t i = 0; i < WINDOW_SIZE / 32; i++)
    {
        vecHw.mul32(32 * i, offset + 32 * i, 32 * i);
    }

    int c = WINDOW_SIZE / 2;

    // c = 512 only working with reals for now
    for (size_t i = 0; i < c / 32; i++)
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
    for (size_t i = 0; i < c / 32; i++)
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
        for (size_t i = 0; i < std::max(1, c / 32); i++)
        {
            size_t vals = REAL(0);
            size_t valsC = REAL(c);

            size_t realC = REAL(c);
            size_t imagC = IMAG(c);

            // we need to run for real and imaginary
            // [0 - C] = [0 - C] + [C - 2C]
            // [C - 2C] = [0 - C] + (-1)*[C - 2C]
            // [2C - 3C] = [2C- 3C] + [3C - 4C]
            // [3C - 4C] = [2C - 3C] + (-1)*[3C - 4C]

            // [4C - 5C] = [4C - 5C] + [5C - 6C]
            // [5C - 6C] = [4C - 5C] + (-1)*[5C - 6C]
            // [6C - 7C] = [6C - 7C] + [7C - 8C]
            // [7C - 8C] = [6C - 7C] + (-1)*[7C - 8C]

            // this is the pair loop
            // it will repeat for each pair of matricies that need to be processed
            // here c = 128, so 4 pairs of cross addition and complex multiplication needs to happen
            /*
            | -- cross addition -- | | complex multiplication |

            j = 0

            j * c   ----\--/---- : vals
                         \/
                         /\
            (j+1)*c ----/--\---- : valsC -> realC * omega[WINDOW_SIZE - 2 * c]

            j++

            j * c   ----\--/---- : vals
                         \/
                         /\
            (j+1)*c ----/--\---- : valsC -> realC * omega[WINDOW_SIZE - 2 * c]

            j++

            j * c   ----\--/---- : vals
                         \/
                         /\
            (j+1)*c ----/--\---- : valsC -> realC * omega[WINDOW_SIZE - 2 * c]

            j++

            j * c   ----\--/---- : vals
                         \/
                         /\
            (j+1)*c ----/--\---- : valsC -> realC * omega[WINDOW_SIZE - 2 * c]
            */
            for (size_t j = 0; j < WINDOW_SIZE / (2 * c); j++)
            {
                // addition stage:
                // (-1)*[3C - 4C]
                vecHw.mul32(valsC, minusOne, buf0, std::min(c, 32));

                // [3C - 4C] = [2C - 3C] + (-1)*[3C - 4C]
                vecHw.add32(vals, buf0, buf1, std::min(c, 32));

                // [2C - 3C] = [2C - 3C] + [3C - 4C]
                vecHw.add32(vals, valsC, vals, std::min(c, 32));

                // [3C - 4C] = [2C - 3C] + (-1)*[3C - 4C]
                vecHw.assign32(buf1, valsC, std::min(c, 32));

                // offset to imaginary numbers and do it again
                vals += WINDOW_SIZE;
                valsC += WINDOW_SIZE;

                if (j != 0) // manually exclude times where imaginary numbers aren't populated
                {
                    // imaginaries
                    // (-1)*[3C - 4C]
                    vecHw.mul32(valsC, minusOne, buf0, std::min(c, 32));

                    // [3C - 4C] = [2C - 3C] + (-1)*[3C - 4C]
                    vecHw.add32(vals, buf0, buf1, std::min(c, 32));

                    // [2C - 3C] = [2C - 3C] + [3C - 4C]
                    vecHw.add32(vals, valsC, vals, std::min(c, 32));

                    // [3C - 4C] = [2C - 3C] + (-1)*[3C - 4C]
                    vecHw.assign32(buf1, valsC, std::min(c, 32));
                }

                // remove the imaginary val offset and move on to the next pair
                vals -= WINDOW_SIZE;
                vals += 2 * c;

                valsC -= WINDOW_SIZE;
                valsC += 2 * c;

                // complex multiplication stage:
                vecHw.mul32(realC, 2 * WINDOW_SIZE + 32 * i + (WINDOW_SIZE - 2 * c), buf0, std::min(c, 32)); // real [C - 2C] * cos(C-2C) -> 2 * WINDOW_SIZE ( we're using this as a buffer cause we're done using it in this fft rn)

                if (j != 0)
                {
                    vecHw.mul32(imagC, 3 * WINDOW_SIZE + 32 * i + (WINDOW_SIZE - 2 * c), buf1, std::min(c, 32)); // imag [C - 2C] *  sin(C-2C) -> 3 * WINDOW_SIZE     
                    vecHw.mul32(imagC, 2 * WINDOW_SIZE + 32 * i + (WINDOW_SIZE - 2 * c), imagC, std::min(c, 32)); // imag [C - 2C] = imag [C - 2C] * cos(C-2C)
                }

                vecHw.mul32(realC, 3 * WINDOW_SIZE + 32 * i + (WINDOW_SIZE - 2 * c), buf2, std::min(c, 32)); // imag [C - 2C] = real[C - 2C] * sin(C-2C)

                if (j == 0)
                {
                    vecHw.assign32(buf0, realC, std::min(c, 32));
                    vecHw.assign32(buf2, imagC, std::min(c, 32));
                }
                else
                {
                    vecHw.mul32(buf1, minusOne, buf1, std::min(c, 32)); // (-1)* (imag [C - 2C] *  sin(C-2C))
                    vecHw.add32(buf0, buf1, realC, std::min(c, 32)); // real [C - 2C] = real [C - 2C] * cos(C-2C)  + (-1)* (imag [C - 2C] *  sin(C-2C))
                    vecHw.add32(imagC, buf2, imagC, std::min(c, 32)); // imag [C - 2C] = (imag [C - 2C] * cos(C-2C)) + (real[C - 2C] * sin(C-2C))
                }

                // repeat for 3C, 5C, and 7C
                realC += 2 * c;
                imagC += 2 * c;
            }
        }


        // rearrage step
        if (c == 32)
        {
            size_t i = 0;
            size_t vals = REAL(0);
            size_t valsC = REAL(c);

            size_t realC = REAL(c);
            size_t imagC = IMAG(c);

            //todo do this while generating the angle terms
            vecHw.assign32(2 * WINDOW_SIZE + (WINDOW_SIZE - 2 * 16), buf3, 16);
            vecHw.assign32(2 * WINDOW_SIZE + (WINDOW_SIZE - 2 * 16), buf3 + 16, 16);
            vecHw.assign32(3 * WINDOW_SIZE + (WINDOW_SIZE - 2 * 16), buf4, 16);
            vecHw.assign32(3 * WINDOW_SIZE + (WINDOW_SIZE - 2 * 16), buf4 + 16, 16);

            // for every other pair of 32 sized arrays
            // we need to rearrange their layout            
            // todo see if i can readd the j != 0 stuff again
            for (size_t j = 0; j < WINDOW_SIZE / (2 * 32); j++)
            {
                // take the last 16 values from the first array and swap them with the first 16 values from the second array
                vecHw.assign32(vals + 16, buf0, 16);
                vecHw.assign32(valsC, vals + 16, 16);
                vecHw.assign32(buf0, valsC, 16);

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

                // take the last 16 values from the first array and swap them with the first 16 values from the second array
                vecHw.assign32(vals + 16, buf0, 16);
                vecHw.assign32(valsC, vals + 16, 16);
                vecHw.assign32(buf0, valsC, 16);

                // imaginaries
                // (-1)*[3C - 4C]
                vecHw.mul32(valsC, minusOne, buf0, 32);

                // [3C - 4C] = [2C - 3C] + (-1)*[3C - 4C]
                vecHw.add32(vals, buf0, buf1);

                // [2C - 3C] = [2C - 3C] + [3C - 4C]
                vecHw.add32(vals, valsC, vals);

                // [3C - 4C] = [2C - 3C] + (-1)*[3C - 4C]
                vecHw.assign32(buf1, valsC);

                // complex multiplication stage:
                vecHw.mul32(realC, buf3, buf0); // real [C - 2C] * cos(C-2C) -> 2 * WINDOW_SIZE ( we're using this as a buffer cause we're done using it in this fft rn)

                vecHw.mul32(imagC, buf4, buf1); // imag [C - 2C] *  sin(C-2C) -> 3 * WINDOW_SIZE     
                vecHw.mul32(imagC, buf3, imagC); // imag [C - 2C] = imag [C - 2C] * cos(C-2C)

                vecHw.mul32(realC, buf4, buf2); // imag [C - 2C] = real[C - 2C] * sin(C-2C)

                vecHw.mul32(buf1, minusOne, buf1, 32); // (-1)* (imag [C - 2C] *  sin(C-2C))
                vecHw.add32(buf0, buf1, realC); // real [C - 2C] = real [C - 2C] * cos(C-2C)  + (-1)* (imag [C - 2C] *  sin(C-2C))
                vecHw.add32(imagC, buf2, imagC); // imag [C - 2C] = (imag [C - 2C] * cos(C-2C)) + (real[C - 2C] * sin(C-2C))

                // repeat for 3C, 5C, and 7C
                realC += 2 * c;
                imagC += 2 * c;

                // take the last 16 values from the first array and swap them with the first 16 values from the second array
                vecHw.assign32(vals + 16, buf0, 16);
                vecHw.assign32(valsC, vals + 16, 16);
                vecHw.assign32(buf0, valsC, 16);

                // remove the imaginary val offset and move on to the next pair
                vals -= WINDOW_SIZE;
                valsC -= WINDOW_SIZE;

                // take the last 16 values from the first array and swap them with the first 16 values from the second array
                vecHw.assign32(vals + 16, buf0, 16);
                vecHw.assign32(valsC, vals + 16, 16);
                vecHw.assign32(buf0, valsC, 16);

                vals += 2 * c;
                valsC += 2 * c;
            }

            c = 8;

            vals = 0;
            valsC = 32;

            realC = 32;
            imagC = WINDOW_SIZE + 32;

            for (size_t i = 0; i < 32 / c; i++)
            {
                vecHw.assign32(2 * WINDOW_SIZE + (WINDOW_SIZE - 2 * c), buf3 + i * c, c);
                vecHw.assign32(3 * WINDOW_SIZE + (WINDOW_SIZE - 2 * c), buf4 + i * c, c);
            }

            for (size_t j = 0; j < WINDOW_SIZE / (2 * 32); j++)
            {
                vecHw.assign32(vals + 16, buf0, 2 * c);
                vecHw.assign32(vals + 32, vals + 16, 2 * c);
                vecHw.assign32(buf0, vals + 32, 2 * c);

                vecHw.assign32(vals + 8, buf0, c);
                vecHw.assign32(vals + 32, vals + 8, c);
                vecHw.assign32(buf0, vals + 32, c);

                vecHw.assign32(vals + 24, buf0, c);
                vecHw.assign32(vals + 48, vals + 24, c);
                vecHw.assign32(buf0, vals + 48, c);

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

                vecHw.assign32(vals + 16, buf0, 2 * c);
                vecHw.assign32(vals + 32, vals + 16, 2 * c);
                vecHw.assign32(buf0, vals + 32, 2 * c);

                vecHw.assign32(vals + 8, buf0, c);
                vecHw.assign32(vals + 32, vals + 8, c);
                vecHw.assign32(buf0, vals + 32, c);

                vecHw.assign32(vals + 24, buf0, c);
                vecHw.assign32(vals + 48, vals + 24, c);
                vecHw.assign32(buf0, vals + 48, c);

                // imaginaries
                // (-1)*[3C - 4C]
                vecHw.mul32(valsC, minusOne, buf0, 32);

                // [3C - 4C] = [2C - 3C] + (-1)*[3C - 4C]
                vecHw.add32(vals, buf0, buf1);

                // [2C - 3C] = [2C - 3C] + [3C - 4C]
                vecHw.add32(vals, valsC, vals);

                // [3C - 4C] = [2C - 3C] + (-1)*[3C - 4C]
                vecHw.assign32(buf1, valsC);

                // complex multiplication stage:
                vecHw.mul32(realC, buf3, buf0); // real [C - 2C] * cos(C-2C) -> 2 * WINDOW_SIZE ( we're using this as a buffer cause we're done using it in this fft rn)

                vecHw.mul32(imagC, buf4, buf1); // imag [C - 2C] *  sin(C-2C) -> 3 * WINDOW_SIZE     
                vecHw.mul32(imagC, buf3, imagC); // imag [C - 2C] = imag [C - 2C] * cos(C-2C)

                vecHw.mul32(realC, buf4, buf2); // imag [C - 2C] = real[C - 2C] * sin(C-2C)

                vecHw.mul32(buf1, minusOne, buf1, 32); // (-1)* (imag [C - 2C] *  sin(C-2C))
                vecHw.add32(buf0, buf1, realC); // real [C - 2C] = real [C - 2C] * cos(C-2C)  + (-1)* (imag [C - 2C] *  sin(C-2C))
                vecHw.add32(imagC, buf2, imagC); // imag [C - 2C] = (imag [C - 2C] * cos(C-2C)) + (real[C - 2C] * sin(C-2C))

                // repeat for 3C, 5C, and 7C
                realC += 64;
                imagC += 64;

                vecHw.assign32(vals + 24, buf0, c);
                vecHw.assign32(vals + 48, vals + 24, c);
                vecHw.assign32(buf0, vals + 48, c);

                vecHw.assign32(vals + 8, buf0, c);
                vecHw.assign32(vals + 32, vals + 8, c);
                vecHw.assign32(buf0, vals + 32, c);

                vecHw.assign32(vals + 16, buf0, 2 * c);
                vecHw.assign32(vals + 32, vals + 16, 2 * c);
                vecHw.assign32(buf0, vals + 32, 2 * c);

                // remove the imaginary val offset and move on to the next pair
                vals -= WINDOW_SIZE;
                valsC -= WINDOW_SIZE;

                vecHw.assign32(vals + 24, buf0, c);
                vecHw.assign32(vals + 48, vals + 24, c);
                vecHw.assign32(buf0, vals + 48, c);

                vecHw.assign32(vals + 8, buf0, c);
                vecHw.assign32(vals + 32, vals + 8, c);
                vecHw.assign32(buf0, vals + 32, c);

                vecHw.assign32(vals + 16, buf0, 2 * c);
                vecHw.assign32(vals + 32, vals + 16, 2 * c);
                vecHw.assign32(buf0, vals + 32, 2 * c);

                vals += 64;
                valsC += 64;
            }
        }
    }

    vecHw.copyFromHw(temp, 0, 2 * WINDOW_SIZE, 0);

    memcpy(inputReal.data(), temp.data(), WINDOW_SIZE * sizeof(ec::Float));
    memcpy(inputImag.data(), temp.data() + WINDOW_SIZE, WINDOW_SIZE * sizeof(ec::Float));

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
        for (size_t j = 0; j < WINDOW_SIZE / (2 * c); j++)
        {
            for (size_t k = 0; k < c; k++)
            {
                if (c != 1 || reverse_bits(valsC) <= 512)
                {
                    buffer[k] = inputReal[vals + k] - inputReal[valsC + k];
                }

                inputReal[vals + k] += inputReal[valsC + k];
            }

            memcpy(inputReal.data() + valsC, buffer.data(), c * sizeof(ec::Float));

            for (size_t k = 0; k < c; k++)
            {
                if (j == 0)
                {
                    memcpy(buffer.data() + k, inputImag.data() + vals + k, sizeof(ec::Float));
                }
                else
                {
                    if (c != 1 || reverse_bits(valsC) <= 512)
                    {
                        buffer[k] = inputImag[vals + k] - inputImag[valsC + k];
                    }

                    if (k == 0 && j == 1)
                    {
                        memcpy(inputImag.data() + vals + k, inputImag.data() + valsC + k, sizeof(ec::Float));
                    }
                    else
                    {
                        inputImag[vals + k] += inputImag[valsC + k];
                    }
                }
            }

            memcpy(inputImag.data() + valsC, buffer.data(), c * sizeof(ec::Float));

            // vals -= WINDOW_SIZE;
            vals += 2 * c;

            // valsC -= WINDOW_SIZE;
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
                else
                {
                    buffer[k] = inputReal[realC + k] * angleTerms[(WINDOW_SIZE - 2 * c) + k] - (inputImag[imagC + k] * angleTerms[WINDOW_SIZE + (WINDOW_SIZE - 2 * c) + k]);
                }

                if (k != 0)
                {
                    inputImag[imagC + k] = inputImag[imagC + k] * angleTerms[(WINDOW_SIZE - 2 * c) + k] + (inputReal[realC + k] * angleTerms[WINDOW_SIZE + (WINDOW_SIZE - 2 * c) + k]);
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

    // fix order of arrays
    ec::Float swapVal;
    for (size_t i = 0; i < WINDOW_SIZE / 2; i++)
    {
        uint16_t newI = reverse_bits(i);

        if (i < newI)
        {
            memcpy(&swapVal, inputReal.data() + i, sizeof(ec::Float));
            memcpy(inputReal.data() + i, inputReal.data() + newI, sizeof(ec::Float));
            memcpy(inputReal.data() + newI, &swapVal, sizeof(ec::Float));

            memcpy(&swapVal, inputImag.data() + i, sizeof(ec::Float));
            memcpy(inputImag.data() + i, inputImag.data() + newI, sizeof(ec::Float));
            memcpy(inputImag.data() + newI, &swapVal, sizeof(ec::Float));
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