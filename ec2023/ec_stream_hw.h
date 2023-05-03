//! Rohde & Schwarz Engineering Competition 2023
//!
//! (c) Rohde & Schwarz GmbH & Co. KG, Munich
//!
//! This code is needed to run your solution locally. Any changes will be lost once you upload your solution!

#pragma once
#include <vector>
#include <queue>

static constexpr size_t STREAM_HW_MEM_SIZE = 4 * 1024; // [4B-values]
static constexpr size_t FIFO_MAX_SIZE = 32; // [4B-values]
static constexpr size_t MAX_NUM_MEM_TO_FIFO_CONNECTIONS = 32;
static constexpr size_t MAX_NUM_OPERATIONS = 1024;

namespace ec
{

  class StreamHw final
  {
  private:
    std::vector<float> m_mem; 
    enum class StreamOp_t { mulVecVec, mulVecScal, addVecVec, addVecScal };

    static inline StreamHw* pSingletonStreamHw = nullptr;

    inline StreamHw()
      : m_mem(VEC_HW_MEM_SIZE, 0.0f)
    {
    }

    StreamHw(const StreamHw& other) = delete;

    class StreamOperation
    {
    public:
      StreamHw::StreamOp_t m_operation;
      std::vector<std::queue<float>>& m_fifos;
      std::vector<uint32_t> m_inputFifosIDs;
      std::vector<uint32_t> m_outputFifosIDs;
      size_t m_precountedStepSize = 0;
      ec::Float m_scalarOperand = 0.0f;

      StreamOperation(StreamHw::StreamOp_t operation,
        std::vector<std::queue<float>>& fifos,
        std::vector<uint32_t> inputFifosIDs,
        std::vector<uint32_t> outputFifosIDs)
        : m_operation(operation)
        , m_fifos(fifos)
        , m_inputFifosIDs(inputFifosIDs)
        , m_outputFifosIDs(outputFifosIDs)
      {
        switch (operation)
        {
        case StreamOp_t::mulVecVec:
        case StreamOp_t::addVecVec:
          if (inputFifosIDs.size() != 2)
            throw std::runtime_error("Wrong number of input Fifos for the chosen operation.");
          if (outputFifosIDs.size() != 1)
            throw std::runtime_error("Wrong number of output Fifos for the chosen operation.");
          break;
        case StreamOp_t::mulVecScal:
        case StreamOp_t::addVecScal:
          // cannot happen here, handled in the other constructor.
          break;
        }
      }

      StreamOperation(StreamHw::StreamOp_t operation, std::vector<std::queue<float>>& fifos,
        std::vector<uint32_t> inputFifosIDs, ec::Float scalarOperand,
        std::vector<uint32_t> outputFifosIDs)
        : m_operation(operation)
        , m_fifos(fifos)
        , m_inputFifosIDs(inputFifosIDs)
        , m_outputFifosIDs(outputFifosIDs)
        , m_scalarOperand(scalarOperand)
      {
        switch (operation)
        {
        case StreamOp_t::mulVecScal: //  i.e. vA x scalar:
        case StreamOp_t::addVecScal: //  i.e. vA + scalar:
          if (inputFifosIDs.size() != 1)
            throw std::runtime_error("Wrong number of input Fifos for the chosen operation.");
          if (outputFifosIDs.size() != 1)
            throw std::runtime_error("Wrong number of output Fifos for the chosen operation.");
          break;
        case StreamOp_t::mulVecVec:
        case StreamOp_t::addVecVec:
          // cannot happen here, handled in the other constructor.
          break;

        }

      }

      inline void countWhatIsToProcessInNextStep() noexcept
      {
        size_t numInputVals = countInputToProcess();
        size_t numOutputVals = countOutputToGenerate();
        m_precountedStepSize = std::min(numInputVals, numOutputVals);
      }

      inline size_t countInputToProcess() const noexcept
      {
        size_t numValsReady = std::numeric_limits<size_t>::max();

        for (const auto& fifoID : m_inputFifosIDs)
        {
          numValsReady = std::min<size_t>(numValsReady, m_fifos[fifoID].size());
        }

        return numValsReady;
      }

      inline size_t countOutputToGenerate() const
      {
        size_t numValsReady = std::numeric_limits<size_t>::max();

        for (const auto& fifoID : m_outputFifosIDs)
        {
          numValsReady = std::min<size_t>(numValsReady, FIFO_MAX_SIZE - m_fifos[fifoID].size());
        }

        return numValsReady;
      }

      inline Measurement::MeasType getMeasurementInfoForOp() const noexcept
      {
        switch (m_operation)
        {
        case StreamOp_t::mulVecVec:
        case StreamOp_t::mulVecScal:
          return Measurement::MeasType::mul_StreamHw;
        case StreamOp_t::addVecVec:
        case StreamOp_t::addVecScal:
          return Measurement::MeasType::add_StreamHw;
        }

        return Measurement::MeasType::add_StreamHw;
      }

      inline size_t executePrecountedStep()
      {
        const size_t numValsReady = m_precountedStepSize;

        if (numValsReady == 0)
          return 0;

        if (countInputToProcess() < m_precountedStepSize)
          throw std::runtime_error("Errorneous precounted operation size!");

        if (countOutputToGenerate() < m_precountedStepSize)
          throw std::runtime_error("Errorneous precounted operation size!");

        switch (m_operation)
        {
        case StreamOp_t::mulVecVec:
          for (size_t J = 0; J < numValsReady; J++)
          {
            m_fifos[m_outputFifosIDs[0]].push(m_fifos[m_inputFifosIDs[0]].front() * m_fifos[m_inputFifosIDs[1]].front());

            m_fifos[m_inputFifosIDs[0]].pop();
            if (m_inputFifosIDs[0] != m_inputFifosIDs[1]) // is not vA x vA
              m_fifos[m_inputFifosIDs[1]].pop();
          }
          break;

        case StreamOp_t::mulVecScal: //  i.e. vA x scalar
          for (size_t J = 0; J < numValsReady; J++)
          {
            m_fifos[m_outputFifosIDs[0]].push(m_fifos[m_inputFifosIDs[0]].front() *
              m_scalarOperand.m_value);

            m_fifos[m_inputFifosIDs[0]].pop();
          }
          break;

        case StreamOp_t::addVecVec:
          for (size_t J = 0; J < numValsReady; J++)
          {
            m_fifos[m_outputFifosIDs[0]].push(m_fifos[m_inputFifosIDs[0]].front() + m_fifos[m_inputFifosIDs[1]].front());

            m_fifos[m_inputFifosIDs[0]].pop();
            if (m_inputFifosIDs[0] != m_inputFifosIDs[1]) // is not vA + vA
              m_fifos[m_inputFifosIDs[1]].pop();
          }
          break;
        case StreamOp_t::addVecScal:
          for (size_t J = 0; J < numValsReady; J++)
          {
            m_fifos[m_outputFifosIDs[0]].push(m_fifos[m_inputFifosIDs[0]].front() + m_scalarOperand.m_value);

            m_fifos[m_inputFifosIDs[0]].pop();
          }
          break;

        }

        // reset precounted step size after execution
        m_precountedStepSize = 0;

        return numValsReady;
      }
    };

    class FifoMemConnection
    {
    public:
      size_t m_idxMemStart;
      size_t m_numEl;
      size_t m_idxMemCurrent;
      uint32_t m_idFifo;

      constexpr FifoMemConnection(size_t idxMemStart, size_t numEl, uint32_t idFifo) noexcept
        : m_idxMemStart(idxMemStart), m_numEl(numEl), m_idxMemCurrent(idxMemStart), m_idFifo(idFifo)
      {
      }
    };

    uint32_t m_numFifos = 0; // IDs of the Fifos = 0, 1, ..., m_numFifos-1;
    std::vector<std::queue<float>> m_fifos;
    std::vector<StreamOperation> m_operations;
    std::vector<FifoMemConnection> m_fifoToMemConnections;
    std::vector<FifoMemConnection> m_memToFifoConnections;

    inline size_t pushDataToFifo(FifoMemConnection& connection)
    {
      // fill the input fifo just to 1/2, so that fifos between two operations can be read and filled by 1/2 in a single pipeline step
      size_t numFreeInFifo = (FIFO_MAX_SIZE / 2) - m_fifos[connection.m_idFifo].size();
      size_t numAvailableInMem = (connection.m_idxMemStart + connection.m_numEl) - connection.m_idxMemCurrent;
      size_t numToPush = std::min(numFreeInFifo, numAvailableInMem);

      if (numToPush == 0)
        return 0;

      for (size_t I = 0; I < numToPush; I++)
      {
        m_fifos[connection.m_idFifo].push(m_mem[connection.m_idxMemCurrent + I]);
      }

      connection.m_idxMemCurrent += numToPush;

      return numToPush;
    }

    inline size_t popDataFromFifo(FifoMemConnection& connection) noexcept
    {
      const size_t numElsInFifo = m_fifos[connection.m_idFifo].size();
      const size_t numAvailableInMem = (connection.m_idxMemStart + connection.m_numEl) - connection.m_idxMemCurrent;
      const size_t numToPop = std::min(numElsInFifo, numAvailableInMem);

      if (numToPop == 0)
        return 0;

      for (size_t I = 0; I < numToPop; I++)
      {
        m_mem[connection.m_idxMemCurrent + I] = m_fifos[connection.m_idFifo].front();
        m_fifos[connection.m_idFifo].pop();
      }

      connection.m_idxMemCurrent += numToPop;

      return numToPop;
    }


  public:

    static inline StreamHw* getSingletonStreamHw()
    {
      if (pSingletonStreamHw == nullptr)
      {
        pSingletonStreamHw = new StreamHw;
      }

      return pSingletonStreamHw;
    }

    inline void startStreamDataMemToFifo(size_t idxMemStart, uint32_t idFifo, size_t numEl)
    {
      if (idxMemStart + numEl > m_mem.size())
        throw std::runtime_error("Size of memory of StreamHw to copy from does not fit the elements to be copied.");
      if (idFifo > m_numFifos)
        throw std::runtime_error("Unknown Fifo ID. Crosscheck your settings of numFifos with the IDs of inputFifos");
      if (m_memToFifoConnections.size() >= MAX_NUM_MEM_TO_FIFO_CONNECTIONS)
        throw std::runtime_error("Maximum number of Memory to Fifo connections was reached. See MAX_NUM_MEM_TO_FIFO_CONNECTIONS limitation of StreamHw.");

      const FifoMemConnection fifoToMem(idxMemStart, numEl, idFifo);
      m_memToFifoConnections.push_back(fifoToMem);
    }

    inline void startStreamDataFifoToMem(uint32_t idFifo, size_t idxMemStart, size_t numEl)
    {
      if (idxMemStart + numEl > m_mem.size())
        throw std::runtime_error("Size of memory of StreamHw to copy from does not fit the elements to be copied.");
      if (idFifo > m_numFifos)
        throw std::runtime_error("Unknown Fifo ID. Crosscheck your settings of numFifos with the IDs of inputFifos");
      if (m_fifoToMemConnections.size() >= MAX_NUM_MEM_TO_FIFO_CONNECTIONS)
        throw std::runtime_error("Maximum number of Fifo to Memory connections was reached. See MAX_NUM_MEM_TO_FIFO_CONNECTIONS limitation of StreamHw.");

      FifoMemConnection fifoToMem(idxMemStart, numEl, idFifo);
      m_fifoToMemConnections.push_back(fifoToMem);
    }

    inline void runPipeline()
    {
      // Pipeline is first finished, when no operation can do any progress.

      bool wasSomeProgress = true;

      while (wasSomeProgress)
      {
        wasSomeProgress = false;
        size_t numOpsMaxPipelineStep = 0;
        Measurement::MeasType opTypeMaxPipelineStep = Measurement::MeasType::addition;
        size_t scoreOfOpMaxPipelineStep = 0;

        // Try to stream data from mem to fifos
        for (auto& connection : m_memToFifoConnections)
        {
          size_t numElsPushed = pushDataToFifo(connection);
          wasSomeProgress = numElsPushed || wasSomeProgress;
        }

        // For every operation, count and store how much it can process in this step
        for (auto& operation : m_operations)
        {
          operation.countWhatIsToProcessInNextStep();
        }

        for (auto& operation : m_operations)
        {
          const size_t numOpsDone = operation.executePrecountedStep();
          wasSomeProgress = numOpsDone || wasSomeProgress;

          const size_t scorePerOp = Measurement::get(operation.getMeasurementInfoForOp()).m_runtimeWeight;

          // search which operation is the bottleneck of the pipeline step
          if (numOpsDone * scorePerOp > numOpsMaxPipelineStep * scoreOfOpMaxPipelineStep)
          {
            numOpsMaxPipelineStep = numOpsDone;
            opTypeMaxPipelineStep = operation.getMeasurementInfoForOp();
            scoreOfOpMaxPipelineStep = scorePerOp;
          }
        }

        // Try to stream data from fifos to mem
        for (auto& connection : m_fifoToMemConnections)
        {
          const size_t numElsPoped = popDataFromFifo(connection);
          wasSomeProgress = numElsPoped || wasSomeProgress;
        }

        Measurement::inc(opTypeMaxPipelineStep, numOpsMaxPipelineStep);
      }

      m_fifoToMemConnections.clear();
      m_memToFifoConnections.clear();
    }

    inline void resetStreamHw()
    {
      resetMemTo0(0, STREAM_HW_MEM_SIZE);
      m_numFifos = 0;
      m_fifos.clear();
      m_operations.clear();
    }

    inline void resetMemTo0(size_t idx_start, size_t num)
    {
      if (idx_start + num > m_mem.size())
        throw std::runtime_error("Size of memory of StreamHw to reset does not fit the size of the memory.");

      std::fill(m_mem.begin() + static_cast<long int>(idx_start), m_mem.begin() + static_cast<long int>(idx_start + num), 0.0f);
    }


    inline void resetMemTo0()
    {
      resetMemTo0(0, m_mem.size());
    }

    inline void copyToHw(const std::vector<ec::Float>& data, size_t idx_from, size_t num, size_t idx_to)
    {
      Measurement::inc(Measurement::MeasType::initCopy_StreamHw);

      if (idx_from + num > data.size())
        throw std::runtime_error("Size of vector to copy does not fit the number of elements to be copied.");

      if (idx_to + num > m_mem.size())
        throw std::runtime_error("Size of memory of StreamHw to copy to does not fit the number of elements to be copied.");

      for (size_t I = 0; I < num; I++)
      {
        Measurement::inc(Measurement::MeasType::copy4B_VecHw);
        m_mem[idx_to + I] = data[idx_from + I].m_value;
      }

    }

    inline void copyFromHw(std::vector<ec::Float>& data, size_t idx_from, size_t num, size_t idx_to)
    {
      Measurement::inc(Measurement::MeasType::initCopy_VecHw);

      if (idx_from + num > m_mem.size())
        throw std::runtime_error("Size of memory of StreamHw to copy from does not fit the elements to be copied.");

      if (idx_to + num > data.size())
        throw std::runtime_error("Size of vector to copy to does not fit the elements to be copied.");

      for (size_t I = 0; I < num; I++)
      {
        Measurement::inc(Measurement::MeasType::copy4B_StreamHw);
        data[idx_to + I] = ec::Float(m_mem[idx_from + I]); // Here is added 1x ec::Float assignment to the score additionaly.
      }
    }

    inline void createFifos(uint32_t numFifos)
    {
      m_numFifos = numFifos;
      m_fifos = std::vector<std::queue<float>>(numFifos);
    }

    inline void addOpMulToPipeline(uint32_t inputFifoLeftId, uint32_t inputFifoRightId, uint32_t outputFifoId)
    {
      if (m_operations.size() >= MAX_NUM_OPERATIONS)
        throw std::runtime_error("Maximum number of operations at StreamHw was reached. See MAX_NUM_OPERATIONS limitation of StreamHw.");

      std::vector<uint32_t> inputFifosIds;
      inputFifosIds.push_back(inputFifoLeftId);
      inputFifosIds.push_back(inputFifoRightId);

      std::vector<uint32_t> outputFifosIds;
      outputFifosIds.push_back(outputFifoId);

      StreamOperation op(StreamOp_t::mulVecVec, m_fifos, inputFifosIds, outputFifosIds);
      m_operations.push_back(op);
    }

    inline void addOpMulToPipeline(uint32_t inputFifoLeftId, ec::Float scalarRight, uint32_t outputFifoId)
    {
      if (m_operations.size() >= MAX_NUM_OPERATIONS)
        throw std::runtime_error("Maximum number of operations at StreamHw was reached. See MAX_NUM_OPERATIONS limitation of StreamHw.");

      std::vector<uint32_t> inputFifosIds;
      inputFifosIds.push_back(inputFifoLeftId);

      std::vector<uint32_t> outputFifosIds;
      outputFifosIds.push_back(outputFifoId);

      StreamOperation op(StreamOp_t::mulVecScal, m_fifos, inputFifosIds, scalarRight, outputFifosIds);
      m_operations.push_back(op);
    }

    inline void addOpAddToPipeline(uint32_t inputFifoLeftId, uint32_t inputFifoRightId, uint32_t outputFifoId)
    {
      if (m_operations.size() >= MAX_NUM_OPERATIONS)
        throw std::runtime_error("Maximum number of operations at StreamHw was reached. See MAX_NUM_OPERATIONS limitation of StreamHw.");

      std::vector<uint32_t> inputFifosIds;
      inputFifosIds.push_back(inputFifoLeftId);
      inputFifosIds.push_back(inputFifoRightId);

      std::vector<uint32_t> outputFifosIds;
      outputFifosIds.push_back(outputFifoId);

      StreamOperation op(StreamOp_t::addVecVec, m_fifos, inputFifosIds, outputFifosIds);
      m_operations.push_back(op);
    }

    inline void addOpAddToPipeline(uint32_t inputFifoLeftId, ec::Float scalarRight, uint32_t outputFifoId)
    {
      if (m_operations.size() >= MAX_NUM_OPERATIONS)
        throw std::runtime_error("Maximum number of operations at StreamHw was reached. See MAX_NUM_OPERATIONS limitation of StreamHw.");

      std::vector<uint32_t> inputFifosIds;
      inputFifosIds.push_back(inputFifoLeftId);

      std::vector<uint32_t> outputFifosIds;
      outputFifosIds.push_back(outputFifoId);

      StreamOperation op(StreamOp_t::addVecScal, m_fifos, inputFifosIds, scalarRight, outputFifosIds);
      m_operations.push_back(op);
    }


  };

} // namespace ec
