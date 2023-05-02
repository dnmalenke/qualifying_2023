//! Rohde & Schwarz Engineering Competition 2023
//!
//! (c) Rohde & Schwarz GmbH & Co. KG, Munich
//!
//! This code is needed to run your solution locally. Any changes will be lost once you upload your solution!

#include "solution.h"
#include <iostream>

bool verify_if_correct(const std::vector<ec::Float>& maxSpectrum, const std::vector<ec::Float>& refSpectrum);

int main()
{
  try
  {
    auto inputSignal = ec::SignalReader::readSignal("./input_signal.txt");
    std::vector<ec::Float> maxSpectrum;
    { // Scope of score measurement - start

      ec::Measurement scoreMeasurement;

      maxSpectrum = process_signal(inputSignal); 
    } // Scope of score measurement - end

    std::vector<ec::Float> refSpectrum = ec::SignalReader::readSignal("./reference_spectrum.txt");
    if (!verify_if_correct(maxSpectrum, refSpectrum))
    {
      std::cout << "Verification of output failed. The computed spectrum is different from the reference spectrum.\n";
      return -1;
    }
  }
  catch (const std::exception& e)
  {
    std::cout << "Uncaught exception: " << e.what() << "\n";
    return -1;
  }
  return 0;
}

// Check if the reversed output signal matches the original input signal.
// The output spectrum is valid if the difference to reference spectrum is smaller than maxDelta.
bool verify_if_correct(const std::vector<ec::Float>& maxSpectrum, const std::vector<ec::Float>& refSpectrum)
{
  const ec::Float maxDelta(0.003f);
  uint32_t numErrs = 0u;

  if (refSpectrum.size() != maxSpectrum.size())
  {
    return false;
  }

  for (size_t I = 0; I < refSpectrum.size(); ++I)
  {
    ec::Float diff = ec_abs(maxSpectrum[I] - refSpectrum[I]);

    if (diff > maxDelta)
    {
      numErrs++;
    }

  }

  if (numErrs > 0)
      return false;

  return true;
}
