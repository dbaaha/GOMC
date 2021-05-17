/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#pragma once

#include "OutputAbstracts.h"
#include "MoveSettings.h"
#include "Coordinates.h"
#include "MoveBase.h"
#include <iostream>
#include "GOMC_Config.h"

#include <boost/archive/tmpdir.hpp>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/assume_abstract.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>

class CheckpointOutput : public OutputableBase
{
public:
  CheckpointOutput(System & sys, StaticVals const& statV);

  ~CheckpointOutput()
  {
    if(outputFile)
      fclose(outputFile);
  }

  virtual void DoOutput(const ulong step);
  virtual void DoOutputRestart(const ulong step);  
  virtual void Init(pdb_setup::Atoms const& atoms,
                    config_setup::Output const& output);
  virtual void Sample(const ulong step) {}


private:
  MoveSettings & moveSetRef;
  MoleculeLookup & molLookupRef;
  BoxDimensions & boxDimRef;
  Molecules const & molRef;
  PRNG & prngRef;
  Coordinates & coordCurrRef;
#if GOMC_LIB_MPI
  PRNG & prngPTRef;
#endif

  bool enableParallelTemperingBool;
  std::string filename;
  FILE* outputFile;
  ulong stepsPerCheckpoint;

  /* Variables to set before call to boost serialize */
  char gomc_version[5];

  uint64_t step;

  const int N = 624;


  uint32_t* saveArray;
  uint32_t location;
  uint32_t left;
  uint32_t seed;

  std::vector<uint> originalMoleculeIndicesVec;

  int8_t enableParallelTempering;

#if GOMC_LIB_MPI
  uint32_t* saveArrayPT;
  uint32_t locationPT;
  uint32_t leftPT;
  uint32_t seedPT;
#endif
  /* Variables to set before call to boost serialize */


  void openOutputFile();
  void setGOMCVersion();
  void printGOMCVersion();
  void setParallelTemperingBoolean();
  void printParallelTemperingBoolean();
  void setStepNumber(ulong stepArg);
  void printStepNumber(ulong step);
  void setRandomNumbers();
  void printRandomNumbers();
  void setMoleculesData();
  void printMoleculesData();
  void setSortedMoleculeIndices();
  void printSortedMoleculeIndices();

#if GOMC_LIB_MPI
  void setRandomNumbersParallelTempering();
  void printRandomNumbersParallelTempering();
#endif
  void printMoveSettingsData();
  void setMoveSettingsData();


  void printVector3DDouble(const std::vector< std::vector< std::vector<double> > > &data);
  void printVector3DUint(const std::vector< std::vector< std::vector<uint> > > &data);
  void printVector2DUint(const std::vector< std::vector< uint > > &data);
  void printVector1DDouble(const std::vector< double > &data);
  void printVector1DUint(const std::vector< uint > &data);
  void printArray1DUint(const uint * data, uint count);

  void write_double_binary(double data);
  void write_uint8_binary(int8_t data);
  void write_uint32_binary(uint32_t data);
  void write_uint64_binary(uint64_t data);

  friend class boost::serialization::access;
  // When the class Archive corresponds to an output archive, the
  // & operator is defined similar to <<.  Likewise, when the class Archive
  // is a type of input archive the & operator is defined similar to >>.
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    // GOMC Version
    ar & boost::serialization::make_array<char>(gomc_version, 5);  
    // Step
    ar & step;
    // PRNG Vars
    ar & boost::serialization::make_array<uint32_t>(saveArray, N + 1);  
    ar & location;
    ar & left;
    ar & seed;
    // Move Settings Vectors
    ar & moveSetRef.scale;
    ar & moveSetRef.acceptPercent;
    ar & moveSetRef.accepted;
    ar & moveSetRef.tries;
    ar & moveSetRef.tempAccepted;
    ar & moveSetRef.tempTries;
    ar & moveSetRef.mp_tries;
    ar & moveSetRef.mp_accepted;
    ar & moveSetRef.mp_t_max;
    ar & moveSetRef.mp_r_max;
    // Sorted Molecule Indices
    ar & originalMoleculeIndicesVec;
    ar & boost::serialization::make_array<uint>(molLookupRef.permutedMoleculeIndices, molLookupRef.molLookupCount);  
    // PT boolean
    ar & enableParallelTempering;
    #if GOMC_LIB_MPI
      // PRNG PT Vars
      ar & boost::serialization::make_array<uint32_t>(saveArrayPT, N + 1);  
      ar & locationPT;
      ar & leftPT;
      ar & seedPT;
    #endif
  }

  /* Variables to set before call to boost serialize */
  char gomc_version[5];

  uint64_t step;

  const int N = 624;

  uint32_t* saveArray;
  uint32_t location;
  uint32_t left;
  uint32_t seed;

  std::vector<uint32_t> originalMoleculeIndicesVec;

  int8_t enableParallelTempering;

#if GOMC_LIB_MPI
  uint32_t* saveArrayPT;
  uint32_t locationPT;
  uint32_t leftPT;
  uint32_t seedPT;
#endif
  }

};
