#ifndef _SpectrumPairSet_H_
#define _SpectrumPairSet_H_

// Module Includes
#include "SpectrumPair.h"

// System Includes
#include <fstream>
#include <map>
#include <string>
#include <vector>

namespace specnets
{
  /**
   * TODO: add description
   */
  class SpectrumPairSet
  {
  public:
    SpectrumPairSet(void);
    SpectrumPairSet(unsigned int nsize);


    int loadFromBinaryFile(const std::string & filename);
    bool saveToBinaryFile(const std::string & filename);

    unsigned int size(void) const;
    void resize(unsigned int newSize, SpectrumPair newPairExemplar =
        SpectrumPair());
    SpectrumPair const & operator[](unsigned int index) const;
    SpectrumPair & operator[](unsigned int index);
    void push_back(const SpectrumPair & newPair);

    void sort_pairs();
    
    void sort_pairs_by_index();

    bool getModificationFrequencies(float resolution, 
                                    std::map<float, float> & modFreqs);

  protected:
    std::vector<SpectrumPair> thePairs;
  };

} //namespace specnets

#endif // _SpectrumPairSet_H_
