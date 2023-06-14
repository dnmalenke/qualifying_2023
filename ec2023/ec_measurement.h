//! Rohde & Schwarz Engineering Competition 2023
//!
//! (c) Rohde & Schwarz GmbH & Co. KG, Munich
//!
//! This code is needed to run your solution locally. Any changes will be lost once you upload your solution!

#pragma once

#include <iostream>
#include <fstream>
#include <mutex>
#include <unordered_map>
#include <string>


namespace ec
{

  class Measurement final
  {
  public:

    enum class MeasType
    {
      assignment,
      addition,
      subtraction,
      multiplication,
      division,
      compare,
      sin,
      asin,
      tan,
      atan,
      sqrt,
      abs,
      max,
      exp,
      round,
      log,
      log10,
      pow,

      // ec::VecHw operations
      assign32_VecHw,
      mul32_VecHw,
      add32_VecHw,
      acc32_VecHw,
      cos4_VecHw,
      sin4_VecHw,
      sqrt4_VecHw,
      initCopy_VecHw,
      copy4B_VecHw,
      use_VecHw,

      // ec::StreamHw operations
      add_StreamHw,
      mul_StreamHw,
      initCopy_StreamHw,
      copy4B_StreamHw,
      use_StreamHw

    };

    struct MeasPoint
    {
     

      std::string m_name;
      size_t m_countSpeed = 0;
      size_t m_speedWeight = 1;
      size_t m_countEnergy = 0;
      size_t m_energyWeight = 1;
    };

    Measurement();
    ~Measurement();
    Measurement(const Measurement&) = delete;
    Measurement& operator=(const Measurement&) = delete;


    static void inc(MeasType measType) noexcept;
    static void incSpeed(Measurement::MeasType measType, size_t count) noexcept;
    static void incEnergy(Measurement::MeasType measType, size_t count) noexcept;
    static MeasPoint get(MeasType measType) noexcept;
    static size_t getCountSpeed(MeasType measType) noexcept;
    static bool isFromVecHw(Measurement::MeasType measType) noexcept;
    static bool isFromStreamHw(Measurement::MeasType measType) noexcept;
    static void addUsageOfVecHw() noexcept;
    static void addUsageOfStreamHw() noexcept;

    size_t calcTotalEnergyScore()const noexcept;
    size_t calcTotalSpeedScore() const noexcept;
  private:
    size_t calcTotalSpeedOperations()const noexcept;
    size_t calcTotalEnergyOperations() const noexcept;
    void writeScoreFile()const noexcept;


    static std::unordered_map<MeasType, MeasPoint> s_measPoints;
    inline static std::mutex s_mutex;
    static float s_percentAnalyzed;
    static size_t s_bonusScore;
    static size_t s_pseudoScore;
  };

  inline std::unordered_map<Measurement::MeasType, Measurement::MeasPoint> Measurement::s_measPoints;
  inline float Measurement::s_percentAnalyzed = 0.0f;
  inline size_t Measurement::s_bonusScore = 0;
  inline size_t Measurement::s_pseudoScore = 0;


  inline size_t Measurement::calcTotalSpeedScore() const noexcept
  {
    size_t totalScore = 0;

    for (const auto& [key, value]: s_measPoints)
    {
      totalScore += (value.m_countSpeed * value.m_speedWeight);
    }

    return totalScore;
  }

  inline size_t Measurement::calcTotalEnergyScore() const noexcept
  {
    size_t totalScore = 0;

    for (const auto& [key, value] : s_measPoints)
    {
      totalScore += (value.m_countEnergy * value.m_energyWeight);
    }
    return totalScore;
  }


  inline size_t Measurement::calcTotalSpeedOperations() const noexcept
  {
    size_t totalScore = 0;

    for (const auto& [key, value] : s_measPoints)
    {
      totalScore += value.m_countSpeed;
    }
    return totalScore;
  }

  inline size_t Measurement::calcTotalEnergyOperations() const noexcept
  {
    size_t totalScore = 0;

    for (const auto& [key, value] : s_measPoints)
    {
      totalScore += value.m_countEnergy;
    }
    return totalScore;
  }


  inline void Measurement::writeScoreFile() const noexcept
  {
    std::lock_guard lock(s_mutex);
    std::ofstream ostrm("./scores.json", std::ofstream::out);

    if(!ostrm.is_open())
    {
        std::cout << "Unable to open scores.json for writing.\n";
        return;
        //throw std::runtime_error("Unable to open file: scores.json");
    }

    ostrm << "{\n";
    ostrm << "  \"energy score\" : "     << calcTotalEnergyScore() << ",\n";
    ostrm << "  \"speed score\" : "      << calcTotalSpeedScore()  << ",\n";
    ostrm << "  \"percent analyzed\" : " << s_percentAnalyzed      << ",\n";
    ostrm << "  \"bonus score\" : -"     << s_bonusScore           << ",\n";
    ostrm << "  \"final score\" : "      << (calcTotalEnergyScore() + calcTotalSpeedScore())/2 - s_bonusScore << ",\n";
    ostrm << "  \"total speed relevant operations\" : " << calcTotalSpeedOperations()  << ",\n";
    ostrm << "  \"total energy relevant operations\" : " << calcTotalEnergyOperations() << ",\n";
    
    for (const auto& [key, value] : s_measPoints)
    {
      ostrm << "  \"" << value.m_name << "\" : " << value.m_countSpeed << ",\n";
    }
    ostrm << "  \"pseudo score\" : "     << s_pseudoScore          << "\n";
    ostrm << "}\n";
  }



  inline Measurement::Measurement()
  {
    std::lock_guard lock(s_mutex);

    s_measPoints.clear();

    s_measPoints.insert({MeasType::assignment,      {"assignments",     0, 1}});
    s_measPoints.insert({ MeasType::addition,       {"additions",       0, 2}});
    s_measPoints.insert({ MeasType::subtraction,    {"subtractions",    0, 2}});
    s_measPoints.insert({ MeasType::multiplication, {"multiplications", 0, 5}});
    s_measPoints.insert({ MeasType::division,       {"divisions",       0, 10}});
    s_measPoints.insert({ MeasType::compare,        {"comparisons",     0, 1}});
    s_measPoints.insert({ MeasType::sin,            {"sin/cos",         0, 15}});
    s_measPoints.insert({ MeasType::asin,           {"asin/acos",       0, 18}});
    s_measPoints.insert({ MeasType::tan,            {"tan",             0, 25}});
    s_measPoints.insert({ MeasType::atan,           {"atan",            0, 22}});
    s_measPoints.insert({ MeasType::exp,            {"exp",             0, 14}});
    s_measPoints.insert({ MeasType::sqrt,           {"square roots",    0, 15}});
    s_measPoints.insert({ MeasType::abs,            {"abs",             0, 1}});
    s_measPoints.insert({ MeasType::max,            {"max",             0, 1}});
    s_measPoints.insert({ MeasType::log,            {"log",             0, 26}});
    s_measPoints.insert({ MeasType::round,          {"round",           0, 8}});
    s_measPoints.insert({ MeasType::log10,          {"log10",           0, 27}});
    s_measPoints.insert({ MeasType::pow,            {"pow",             0, 62}});

    s_measPoints.insert({MeasType::assign32_VecHw,  {"assign32_VecHw",  0, 4}});
    s_measPoints.insert({MeasType::mul32_VecHw,     {"mul32_VecHw",     0, 16}});
    s_measPoints.insert({MeasType::add32_VecHw,     {"add32_VecHw",     0, 8}});
    s_measPoints.insert({MeasType::acc32_VecHw,     {"acc32_VecHw",     0, 8}});      
    s_measPoints.insert({MeasType::cos4_VecHw,      {"cos4_VecHw",      0, 16}});
    s_measPoints.insert({MeasType::sin4_VecHw,      {"sin4_VecHw",      0, 16}});
    s_measPoints.insert({MeasType::sqrt4_VecHw,     {"sqrt4_VecHw",     0, 16}});
    s_measPoints.insert({ MeasType::copy4B_VecHw,   {"copy4B_VecHw",    0, 2}});
    s_measPoints.insert({ MeasType::initCopy_VecHw, {"initCopy_VecHw",  0, 100}});
    s_measPoints.insert({MeasType::use_VecHw,        {"use_VecHw",      0, 0}});

    s_measPoints.insert({MeasType::add_StreamHw,    {"add_StreamHw",      0, 2}});
    s_measPoints.insert({MeasType::mul_StreamHw,    {"mul_StreamHw",      0, 5}});
    s_measPoints.insert({MeasType::copy4B_StreamHw, {"copy4B_StreamHw",   0, 2}});
    s_measPoints.insert({MeasType::initCopy_StreamHw, {"initCopy_StreamHw", 0, 100}});
    s_measPoints.insert({MeasType::use_StreamHw,      {"use_StreamHw",    0, 0}});

    s_measPoints[MeasType::assignment].m_energyWeight = 1;
    s_measPoints[MeasType::addition].m_energyWeight = 1;
    s_measPoints[MeasType::subtraction].m_energyWeight = 1;
    s_measPoints[MeasType::multiplication].m_energyWeight = 1;
    s_measPoints[MeasType::division].m_energyWeight = 4;
    s_measPoints[MeasType::compare].m_energyWeight = 1;
    s_measPoints[MeasType::sin].m_energyWeight = 4;
    s_measPoints[MeasType::asin].m_energyWeight = 4;
    s_measPoints[MeasType::tan].m_energyWeight = 4;
    s_measPoints[MeasType::atan].m_energyWeight = 4;
    s_measPoints[MeasType::exp].m_energyWeight = 4;
    s_measPoints[MeasType::sqrt].m_energyWeight = 4;
    s_measPoints[MeasType::abs].m_energyWeight = 1;
    s_measPoints[MeasType::max].m_energyWeight = 1;
    s_measPoints[MeasType::log].m_energyWeight = 4;
    s_measPoints[MeasType::round].m_energyWeight = 1;
    s_measPoints[MeasType::log10].m_energyWeight = 4;
    s_measPoints[MeasType::pow].m_energyWeight = 4;

    s_measPoints[MeasType::assign32_VecHw].m_energyWeight = 32;
    s_measPoints[MeasType::mul32_VecHw].m_energyWeight = 32;
    s_measPoints[MeasType::add32_VecHw].m_energyWeight = 32;
    s_measPoints[MeasType::acc32_VecHw].m_energyWeight = 32;
    s_measPoints[MeasType::sin4_VecHw].m_energyWeight = 4;
    s_measPoints[MeasType::cos4_VecHw].m_energyWeight = 4;
    s_measPoints[MeasType::sqrt4_VecHw].m_energyWeight = 4;
    s_measPoints[MeasType::copy4B_VecHw].m_energyWeight = 1;
    s_measPoints[MeasType::initCopy_VecHw].m_energyWeight = 100;
    s_measPoints[MeasType::use_VecHw].m_energyWeight = 100000;

    s_measPoints[MeasType::add_StreamHw].m_energyWeight = 1;
    s_measPoints[MeasType::mul_StreamHw].m_energyWeight = 2;
    s_measPoints[MeasType::copy4B_StreamHw].m_energyWeight = 1;
    s_measPoints[MeasType::initCopy_StreamHw].m_energyWeight = 100;
    s_measPoints[MeasType::use_StreamHw].m_energyWeight = 100000;
  }

  inline Measurement::~Measurement()
  {
    writeScoreFile();
    s_measPoints.clear();
  }

  inline bool Measurement::isFromVecHw(Measurement::MeasType measType) noexcept
  {
    if (s_measPoints[measType].m_name.find(std::string("_VecHw")) != std::string::npos)
      return true;

    return false;
  }

  inline bool Measurement::isFromStreamHw(Measurement::MeasType measType) noexcept
  {
    if (s_measPoints[measType].m_name.find(std::string("_StreamHw")) != std::string::npos)
      return true;

    return false;
  }

  inline void Measurement::addUsageOfVecHw() noexcept
  {
    if (s_measPoints[MeasType::use_VecHw].m_countSpeed == 0)
    {
      s_measPoints[MeasType::use_VecHw].m_countSpeed = 1;
      s_measPoints[MeasType::use_VecHw].m_countEnergy = 1;
    }
  }

  inline void Measurement::addUsageOfStreamHw() noexcept
  {
    if (s_measPoints[MeasType::use_StreamHw].m_countSpeed == 0)
    {
      s_measPoints[MeasType::use_StreamHw].m_countSpeed = 1;
      s_measPoints[MeasType::use_StreamHw].m_countEnergy = 1;
    }

  }

  inline void Measurement::inc(Measurement::MeasType measType) noexcept
  {
    std::lock_guard lock(s_mutex);

    if (isFromVecHw(measType)) addUsageOfVecHw();
    if (isFromStreamHw(measType)) addUsageOfStreamHw();

    s_measPoints[measType].m_countSpeed += 1;
    s_measPoints[measType].m_countEnergy += 1;

  }

  inline void Measurement::incSpeed(Measurement::MeasType measType, size_t count) noexcept
  {
    std::lock_guard lock(s_mutex);

    if (isFromVecHw(measType)) addUsageOfVecHw();
    if (isFromStreamHw(measType)) addUsageOfStreamHw();

    s_measPoints[measType].m_countSpeed += count;
  }

  inline void Measurement::incEnergy(Measurement::MeasType measType, size_t count) noexcept
  {
    std::lock_guard lock(s_mutex);

    if (isFromVecHw(measType)) addUsageOfVecHw();
    if (isFromStreamHw(measType)) addUsageOfStreamHw();

    s_measPoints[measType].m_countEnergy += count;
  }

  inline ec::Measurement::MeasPoint Measurement::get(MeasType measType) noexcept
  {
    std::lock_guard lock(s_mutex);
    return s_measPoints[measType];
  }

  inline size_t Measurement::getCountSpeed(MeasType measType) noexcept
  {
    std::lock_guard lock(s_mutex);
    return s_measPoints[measType].m_countSpeed;
  }

} // namespace ec
