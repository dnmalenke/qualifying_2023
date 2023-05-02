//! Rohde & Schwarz Engineering Competition 2023
//!
//! (c) Rohde & Schwarz GmbH & Co. KG, Munich
//!
//! This code is needed to run your solution locally. Any changes will be lost once you upload your solution!

#pragma once

#include "ec2023/ec_float.h"

#include <array>
#include <fstream>
#include <string>
#include <vector>

namespace ec
{

  class SignalReader final
  {
  public:
    static std::vector<Float> readSignal(const std::string& filePath);
  };



  inline std::vector<ec::Float> SignalReader::readSignal(const std::string& filePath)
  {
    std::vector<Float> ret;

    std::ifstream istrm(filePath, std::ios::binary);
    if (!istrm.is_open())
    {
      throw std::runtime_error("Failed to open: " + filePath);
    }

    std::array<char, 16> input = { 0 };

    while (istrm.getline(&input[0], 16, ','))
    {
      std::string date(&input[0], input.size());
      ret.emplace_back(std::stof(date));
    }
    
    return ret;
  }
}
