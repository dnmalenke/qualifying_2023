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
void fft2(std::vector<ec::Float>& inputReal, std::vector<ec::Float>& inputImag, size_t count);

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


    for (size_t j = 0; j < numWins; j++)
    {
        for (size_t i = 0; i < WINDOW_SIZE; i++)
        {
            signalWindow[i] = inputSignal[i + idxStartWin] * blackmanCoefs[i];
        }

        std::vector<ec::Float> inImag(WINDOW_SIZE);

        if (false)
        {
            ec::VecHw& vecHw = *ec::VecHw::getSingletonVecHw();
            vecHw.resetMemTo0();

            vecHw.copyToHw(angleTerms, 0, 2 * WINDOW_SIZE - 64, 0);
            fft2(signalWindow, inImag, WINDOW_SIZE);

        }
        else
        {
            ec::VecHw& vecHw = *ec::VecHw::getSingletonVecHw();
            vecHw.resetMemTo0();

            vecHw.copyToHw(angleTerms, 0, 2 * WINDOW_SIZE, 2 * WINDOW_SIZE);

            fft(signalWindow, inImag, WINDOW_SIZE);
        }


        for (size_t i = 0; i < sizeSpectrum; i++)
        {
            ec::Float freqVal = signalWindow[i] * signalWindow[i] + inImag[i] * inImag[i];
            // std::cout << freqVal.toFloat() << std::endl;
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

#define REAL(x) x + 32 * i
#define IMAG(x) WINDOW_SIZE + x + 32 * i

uint16_t reverse_bits(uint16_t x) {
    x = ((x & 0xAAAA) >> 1) | ((x & 0x5555) << 1);
    x = ((x & 0xCCCC) >> 2) | ((x & 0x3333) << 2);
    x = ((x & 0xF0F0) >> 4) | ((x & 0x0F0F) << 4);
    x = ((x & 0xFF00) >> 8) | ((x & 0x00FF) << 8);
    return x >> 6;
}

void fft(std::vector<ec::Float>& inputReal, std::vector<ec::Float>& inputImag, size_t count)
{
    ec::VecHw& vecHw = *ec::VecHw::getSingletonVecHw();
    // vecHw.resetMemTo0();

    for (ec::Float i(1.0f); i <= 1024; i++)
    {
        inputReal[(int)i.toFloat() - 1] = i;
    }

    vecHw.copyToHw(inputReal, 0, WINDOW_SIZE, 0);

    int c = WINDOW_SIZE / 2;

    // c = 512 only working with reals for now
    for (size_t i = 0; i < c / 32; i++)
    {
        // (-1)*inputReal[C - 2C]
        vecHw.mul32(REAL(c), minusOne, IMAG(c)); // stick this in imag[C - 2C] for now. it will be fixed later

        // [0 - C] = [0 - C] + [C - 2C]
        vecHw.add32(REAL(0), REAL(c), IMAG(0)); // store in imag[0] for now.

        // [C - 2C] = [0 - C] + (-1)*[C - 2C]
        vecHw.add32(REAL(0), IMAG(c), REAL(c));

        // now we can assign it to [0 - C] because we're done reading it
        vecHw.assign32(IMAG(0), REAL(0));

        vecHw.mul32(IMAG(0), ec::Float(0.0f), IMAG(0)); // just in case. probably can remove later
    }

    /*
    input[C - 2C] * omega[0 through C] # this is complex multiplication
    example of complex multiplication:
    (x + yj) * (a + bj) = (xa - yb) + (xb + ya)j
    So the real part of the product is (xa - yb), and the imaginary part of the product is (xb + ya).
   */
    for (size_t i = 0; i < c / 32; i++)
    {
        vecHw.mul32(REAL(c), 3 * WINDOW_SIZE + 32 * i + (WINDOW_SIZE - 2 * c), IMAG(c)); // imag [C - 2C] = real[C - 2C] * sin(0 - C)
        vecHw.mul32(REAL(c), 2 * WINDOW_SIZE + 32 * i + (WINDOW_SIZE - 2 * c), REAL(c)); // real [C - 2C] = real [C - 2C] * cos(0 - C)
    }

    // Note we have now used the first 512 values in 2 * WINDOW_SIZE and 3 * WINDOW_SIZE and will no longer need them.

    // c = c / 2; // c = 256

    // for (size_t i = 0; i < c / 32; i++)
    // {
    //     size_t vals = REAL(0);
    //     size_t valsC = REAL(c);
    //     size_t buf0 = 2 * WINDOW_SIZE + 32 * i;
    //     size_t buf1 = 2 * WINDOW_SIZE + c + 32 * i;

    //     // (-1)*[C - 2C]
    //     vecHw.mul32(valsC, minusOne, buf0);

    //     // [0 - C] = [0 - C] + [C - 2C]
    //     vecHw.add32(vals, valsC, buf1);

    //     // [C - 2C] = [0 - C] + (-1)*[C - 2C]
    //     vecHw.add32(vals, buf0, valsC);

    //     // [0 - C] = [0 - C] + [C - 2C]
    //     vecHw.assign32(buf1, vals);

    //     vals += 2 * c; // real[2C - 3C]
    //     valsC += 2 * c; // real[3C - 4C]

    //     // (-1)*[3C - 4C]
    //     vecHw.mul32(valsC, minusOne, buf0);

    //     // [3C - 4C] = [2C - 3C] + (-1)*[3C - 4C]
    //     vecHw.add32(vals, buf0, buf1);

    //     // [2C - 3C] = [2C - 3C] + [3C - 4C]
    //     vecHw.add32(vals, valsC, vals);

    //     // [3C - 4C] = [2C - 3C] + (-1)*[3C - 4C]
    //     vecHw.assign32(buf1, valsC);

    //     vals += WINDOW_SIZE; // imag[2C - 3C]
    //     valsC += WINDOW_SIZE; // imag[3C - 4C]

    //     // we now have imaginaries so we need to work with them.
    //      // (-1)*[3C - 4C]
    //     vecHw.mul32(valsC, minusOne, buf0);

    //     // [3C - 4C] = [2C - 3C] + (-1)*[3C - 4C]
    //     vecHw.add32(vals, buf0, buf1);

    //     // [2C - 3C] = [2C - 3C] + [3C - 4C]
    //     vecHw.add32(vals, valsC, vals);

    //     // [3C - 4C] = [2C - 3C] + (-1)*[3C - 4C]
    //     vecHw.assign32(buf1, valsC);
    // }


    // // input[C - 2C] * omega[0 through C] 
    // // input[3C - 4C] * omega[0 through C] 
    // for (size_t i = 0; i < c / 32; i++)
    // {
    //     size_t realC = REAL(c);
    //     size_t imagC = IMAG(c);
    //     size_t buf0 = 2 * WINDOW_SIZE + 32 * i;
    //     size_t buf1 = 2 * WINDOW_SIZE + c + 32 * i;
    //     size_t buf2 = 3 * WINDOW_SIZE + 32 * i;

    //     for (size_t j = 0; j < 2; j++)
    //     {
    //         vecHw.mul32(realC, 2 * WINDOW_SIZE + 32 * i + (WINDOW_SIZE - 2 * c), buf0); // real [C - 2C] * cos(C-2C) -> 2 * WINDOW_SIZE ( we're using this as a buffer cause we're done using it in this fft rn)
    //         vecHw.mul32(imagC, 3 * WINDOW_SIZE + 32 * i + (WINDOW_SIZE - 2 * c), buf1); // imag [C - 2C] *  sin(C-2C) -> 3 * WINDOW_SIZE 
    //         vecHw.mul32(buf1, minusOne, buf1); // (-1)* (imag [C - 2C] *  sin(C-2C)) THIS CAN BE OPTIMIZED BY PRE-NEGATING THE SIN TERMS

    //         vecHw.mul32(imagC, 2 * WINDOW_SIZE + 32 * i + (WINDOW_SIZE - 2 * c), imagC); // imag [C - 2C] = imag [C - 2C] * cos(C-2C)

    //         vecHw.mul32(realC, 3 * WINDOW_SIZE + 32 * i + (WINDOW_SIZE - 2 * c), buf2); // imag [C - 2C] = real[C - 2C] * sin(C-2C)

    //         vecHw.add32(buf0, buf1, realC); // real [C - 2C] = real [C - 2C] * cos(C-2C)  + (-1)* (imag [C - 2C] *  sin(C-2C))

    //         vecHw.add32(imagC, buf2, imagC); // imag [C - 2C] = (imag [C - 2C] * cos(C-2C)) + (real[C - 2C] * sin(C-2C))

    //         // repeat for 3C
    //         realC += 2 * c;
    //         imagC += 2 * c;
    //     }
    // }

    // vecHw.copyFromHw(inputReal, 0, WINDOW_SIZE, 0);
    // vecHw.copyFromHw(inputImag, WINDOW_SIZE, WINDOW_SIZE, 0);
    // std::cout << inputReal[512].toFloat() << " " << inputReal[513].toFloat() << " " << inputReal[514].toFloat() << " " << inputReal[515].toFloat() << std::endl;
    // std::cout << inputReal[512 + 128].toFloat() << " " << inputReal[513 + 128].toFloat() << " " << inputReal[514 + 128].toFloat() << " " << inputReal[515 + 128].toFloat() << std::endl;
    // std::cout << inputImag[512 + 256].toFloat() << " " << inputImag[513 + 256].toFloat() << " " << inputImag[514 + 256].toFloat() << " " << inputImag[515 + 256].toFloat() << std::endl;
    while (c > 1)
    {
        c /= 2;

        for (size_t i = 0; i < std::max(1, c / 32); i++)
        {
            size_t vals = REAL(0);
            size_t valsC = REAL(c);
            size_t buf0 = 2 * WINDOW_SIZE + 32 * i;
            size_t buf1 = 2 * WINDOW_SIZE + c + 32 * i;

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
                // (-1)*[3C - 4C]
                vecHw.mul32(valsC, minusOne, buf0, std::min(c, 32));

                // [3C - 4C] = [2C - 3C] + (-1)*[3C - 4C]
                vecHw.add32(vals, buf0, buf1, std::min(c, 32));

                // [2C - 3C] = [2C - 3C] + [3C - 4C]
                vecHw.add32(vals, valsC, vals, std::min(c, 32));

                // [3C - 4C] = [2C - 3C] + (-1)*[3C - 4C]
                vecHw.assign32(buf1, valsC, std::min(c, 32));

                vals += WINDOW_SIZE; // imag
                valsC += WINDOW_SIZE; // imag

                // imaginaries
                // (-1)*[3C - 4C]
                vecHw.mul32(valsC, minusOne, buf0, std::min(c, 32));

                // [3C - 4C] = [2C - 3C] + (-1)*[3C - 4C]
                vecHw.add32(vals, buf0, buf1, std::min(c, 32));

                // [2C - 3C] = [2C - 3C] + [3C - 4C]
                vecHw.add32(vals, valsC, vals, std::min(c, 32));

                // [3C - 4C] = [2C - 3C] + (-1)*[3C - 4C]
                vecHw.assign32(buf1, valsC, std::min(c, 32));

                vals -= WINDOW_SIZE;
                vals += 2 * c;

                valsC -= WINDOW_SIZE;
                valsC += 2 * c;
            }
        }

        // input[C - 2C] * omega[0 through C] 
        // input[3C - 4C] * omega[0 through C] 
        // input[5C - 6C] * omega[0 through C] 
        // input[7C - 8C] * omega[0 through C] 
        for (size_t i = 0; i < std::max(1, c / 32); i++)
        {
            size_t realC = REAL(c);
            size_t imagC = IMAG(c);
            size_t buf0 = 2 * WINDOW_SIZE + 32 * i;
            size_t buf1 = 2 * WINDOW_SIZE + c + 32 * i;
            size_t buf2 = 3 * WINDOW_SIZE + 32 * i;

            for (size_t j = 0; j < WINDOW_SIZE / (2 * c); j++)
            {
                vecHw.mul32(realC, 2 * WINDOW_SIZE + 32 * i + (WINDOW_SIZE - 2 * c), buf0, std::min(c, 32)); // real [C - 2C] * cos(C-2C) -> 2 * WINDOW_SIZE ( we're using this as a buffer cause we're done using it in this fft rn)
                vecHw.mul32(imagC, 3 * WINDOW_SIZE + 32 * i + (WINDOW_SIZE - 2 * c), buf1, std::min(c, 32)); // imag [C - 2C] *  sin(C-2C) -> 3 * WINDOW_SIZE 
                vecHw.mul32(buf1, minusOne, buf1); // (-1)* (imag [C - 2C] *  sin(C-2C)) THIS CAN BE OPTIMIZED BY PRE-NEGATING THE SIN TERMS

                vecHw.mul32(imagC, 2 * WINDOW_SIZE + 32 * i + (WINDOW_SIZE - 2 * c), imagC, std::min(c, 32)); // imag [C - 2C] = imag [C - 2C] * cos(C-2C)

                vecHw.mul32(realC, 3 * WINDOW_SIZE + 32 * i + (WINDOW_SIZE - 2 * c), buf2, std::min(c, 32)); // imag [C - 2C] = real[C - 2C] * sin(C-2C)

                vecHw.add32(buf0, buf1, realC, std::min(c, 32)); // real [C - 2C] = real [C - 2C] * cos(C-2C)  + (-1)* (imag [C - 2C] *  sin(C-2C))

                vecHw.add32(imagC, buf2, imagC, std::min(c, 32)); // imag [C - 2C] = (imag [C - 2C] * cos(C-2C)) + (real[C - 2C] * sin(C-2C))

                // repeat for 3C, 5C, and 7C
                realC += 2 * c;
                imagC += 2 * c;
            }
        }
    }


    vecHw.copyFromHw(inputReal, 0, WINDOW_SIZE, 0);
    vecHw.copyFromHw(inputImag, WINDOW_SIZE, WINDOW_SIZE, 0);

    // fix order of arrays
    ec::Float temp(0.0f);
    for (size_t i = 0; i < WINDOW_SIZE / 2; i++)
    {
        uint16_t newI = reverse_bits(i);
        memcpy(&temp, inputReal.data() + i, sizeof(ec::Float));
        memcpy(inputReal.data() + i, inputReal.data() + newI, sizeof(ec::Float));
        memcpy(inputReal.data() + newI, &temp, sizeof(ec::Float));

        memcpy(&temp, inputImag.data() + i, sizeof(ec::Float));
        memcpy(inputImag.data() + i, inputImag.data() + newI, sizeof(ec::Float));
        memcpy(inputImag.data() + newI, &temp, sizeof(ec::Float));
    }

    //  vecHw.copyFromHw(inputReal, 0, WINDOW_SIZE, 0);
    // vecHw.copyFromHw(inputImag, WINDOW_SIZE, WINDOW_SIZE, 0);

    std::cout << inputReal[0].toFloat() << " " << inputReal[1].toFloat() << " " << inputReal[2].toFloat() << " " << inputReal[3].toFloat() << std::endl;
    std::cout << inputImag[0].toFloat() << " " << inputImag[1].toFloat() << " " << inputImag[2].toFloat() << " " << inputImag[3].toFloat() << std::endl;
    // std::cout << inputReal[256].toFloat() << " " << inputReal[257].toFloat() << " " << inputReal[258].toFloat() << " " << inputReal[259].toFloat() << std::endl;
    exit(1);
}

// This algorithm breaks everything up recursively into many function calls
// if we can somehow make it iterative and have access to all of it's working arrays
// it might be possible to combine those working arrays into larger arrays that can be
// fed through a really fast pipeline
void fft2(std::vector<ec::Float>& inputReal, std::vector<ec::Float>& inputImag, size_t count)
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

        fft2(even, evenI, halfCount);
        fft2(odd, oddI, halfCount);
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

        // we need to use the symmetry of v1, v2 and c2 to reduce operations.
        // this is about 25% of the computation power
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
        // this is about 50% of the computation time
        // a lot of it is going towards assignments and multiplication and idk how to make it faster cause we're doing a bunch of little operations.

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
    size_t idx = 0;

    ec::Float aC(float(-2.0f * M_PI) / WINDOW_SIZE);

    for (size_t i = 0; i < WINDOW_SIZE / 2; i++)
    {
        angleTerms[idx] = ec_cos(aC * i);
        angleTerms[idx + WINDOW_SIZE] = ec_sin(aC * i);
        idx++;
    }

    size_t x = WINDOW_SIZE / 4;
    size_t y = 2;
    while (x > 0)
    {
        for (size_t i = 0; i < x; i += y)
        {
            memcpy(angleTerms.data() + WINDOW_SIZE - (2 * x) + i / y, angleTerms.data() + i, sizeof(ec::Float));
            memcpy(angleTerms.data() + WINDOW_SIZE + WINDOW_SIZE - (2 * x) + i / y, angleTerms.data() + WINDOW_SIZE + i, sizeof(ec::Float));
        }

        x /= 2;
        y *= 2;
    }






    // size_t count = WINDOW_SIZE;


    // while (count >= 2)
    // {


    //     for (size_t i = 0; i < count / 2; i++)
    //     {
    //         if (i != 0)
    //         {
    //             // if (!(idx >= 977 && idx <= 991 || idx >= 1001 && idx <= 1007 || idx >= 1013 && idx <= 1015 || idx == 1019))
    //             // {
    //             angleTerms[idx] = ec_cos(aC * i);
    //             angleTerms[idx + WINDOW_SIZE] = ec_sin(aC * i);
    //             // }
    //         }
    //          idx++;

    //     }

    //     count >>= 1;
    // }
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
