// Header Include
#include "Logger.h"
#include "SpectrumPairSet.h"

#include <algorithm>
#include <math.h>
#include <iostream>

using namespace std;

namespace specnets
{


bool SpectrumPairComparator (SpectrumPair sp1, SpectrumPair sp2)
{
	float val1 = sp1.score1 + sp1.score2;
	float val2 = sp2.score1 + sp2.score2;
	return (val1 >= val2);
}

bool SpectrumPairComparatorByIndex (SpectrumPair sp1, SpectrumPair sp2)
{       
        float val1 = sp1.spec1;
        float val2 = sp2.spec1;
        return (val1 < val2);
}

	//---------------------------------------------------------------------------
  SpectrumPairSet::SpectrumPairSet(void)
  {
    resize(0);
  }

  //---------------------------------------------------------------------------
  SpectrumPairSet::SpectrumPairSet(unsigned int nsize)
  {
    thePairs.resize(nsize);
  }

  //---------------------------------------------------------------------------
  unsigned int SpectrumPairSet::size(void) const
  {
    return thePairs.size();
  }

  //---------------------------------------------------------------------------
  void SpectrumPairSet::resize(unsigned int newSize,
                               SpectrumPair newPairExemplar)
  {
    thePairs.resize(newSize, newPairExemplar);
  }

  //---------------------------------------------------------------------------
  SpectrumPair const & SpectrumPairSet::operator[](unsigned int index) const
  {
    return thePairs[index];
  }

  //---------------------------------------------------------------------------
  SpectrumPair & SpectrumPairSet::operator[](unsigned int index)
  {
    return thePairs[index];
  }

  //---------------------------------------------------------------------------
  void SpectrumPairSet::sort_pairs()
  {
	  // using function as comp
	    sort (thePairs.begin(), thePairs.end(), SpectrumPairComparator);
  }
  //---------------------------------------------------------------------------
  void SpectrumPairSet::sort_pairs_by_index()
  {
          // using function as comp
            sort (thePairs.begin(), thePairs.end(), SpectrumPairComparatorByIndex);
  }
  //---------------------------------------------------------------------------
  int SpectrumPairSet::loadFromBinaryFile(const std::string & filename)
  {
    FILE *fp;
    unsigned int numEntries, vCount, resIdx, vIdx;
    fp = fopen(filename.c_str(), "rb");
    if (fp == 0)
    {
      ERROR_MSG("Can not open: " << filename);
      return -1;
    }

    fread(&numEntries, sizeof(unsigned int), 1, fp); // Number of entries in the file
    resize(numEntries);
    if (numEntries == 0)
    {
      fclose(fp);
      ERROR_MSG("No entries in file: " << filename);
      return 0;
    }
    vCount = operator[](0).loadSz();
    float *data = (float *)new float[vCount];
    int read_result;
    vector<float> dataV(vCount);
    for (unsigned int resIdx = 0; resIdx < numEntries; resIdx++)
    {
      read_result = fread(data, sizeof(float), vCount, fp);
      if (read_result != vCount)
      {
        resize(0);
        ERROR_MSG("Can not read " << filename);
        return -2;
      }
      for (vIdx = 0; vIdx < 2; vIdx++)
      {
        dataV[vIdx] = data[vIdx] - 1; // Spectrum indices are 1-based
      }
      for (vIdx = 2; vIdx < vCount; vIdx++)
      {
        dataV[vIdx] = data[vIdx];
      }
      operator[](resIdx).load(dataV);
      operator[](resIdx).spec2rev = operator[](resIdx).shift2 > -0.01;
    }
    delete[] data;
    fclose(fp);
    return numEntries;
  }

  //---------------------------------------------------------------------------
  bool SpectrumPairSet::saveToBinaryFile(const std::string & filename)
  {
    FILE *fp;
    vector<float> dataV;
    float data[6];

    fp = fopen(filename.c_str(), "wb");
    if (fp == 0)
    {
      ERROR_MSG("Can not open: " << filename);
      return false;
    }
    unsigned int numEntries = thePairs.size();
    fwrite(&numEntries, sizeof(unsigned int), 1, fp); // Number of entries in the file
    for (unsigned int i = 0; i < numEntries; i++)
    {
      thePairs[i].serialize(dataV);
      for (unsigned int j = 0; j < 2; j++)
      {
        data[j] = dataV[j] + 1; // Spectrum indices are 1-based
      }
      for (unsigned int j = 2; j < dataV.size(); j++)
      {
        data[j] = dataV[j];
      }
      fwrite(data, sizeof(float), dataV.size(), fp);
    }
    fclose(fp);
    return true;
  }

  //---------------------------------------------------------------------------
  void SpectrumPairSet::push_back(const SpectrumPair & newPair)
  {
    thePairs.push_back(newPair);
  }

  //-----------------------------------------------------------------------------
  bool SpectrumPairSet::getModificationFrequencies(float resolution,
                                                   map<float, float> & modFreqs)
  {
    size_t fpsize = thePairs.size();
    DEBUG_VAR(fpsize);

    float totalCounts = 0.0;
    for (size_t i = 0; i < fpsize; i++)
    {
      SpectrumPair sppair = thePairs[i];
      //DEBUG_VAR(sppair.shift1);
      //DEBUG_VAR(sppair.shift2);
      float shift1 = sppair.shift1;
      float shift2 = sppair.shift2;
      // Round to the desired precision
      shift1 = fabs((float)((int)(shift1 / resolution)) * resolution);
      shift2 = fabs((float)((int)(shift2 / resolution)) * resolution);
      //DEBUG_VAR(shift1);
      //DEBUG_VAR(shift2);
      if (shift1 != 0) {
        modFreqs[shift1] += 1.0;
        totalCounts += 1.0;
      }

      if (shift2 != 0) {
        modFreqs[shift2] += 1.0;
        totalCounts += 1.0;
      }
    } // for (size_t i = 0; i < fpsize; i++)

#if 0
    map<int, int> modHist;
    // Count the number of masses that have 'count' number of modifications
    map<float, float>::iterator itr = modFreqs.begin();
    map<float, float>::iterator itrEnd = modFreqs.end();
    for (; itr != itrEnd; itr++) {
      //DEBUG_MSG("modCount: " << itr->first << "  " << itr->second);
      if (itr->first <= 100.0) {
        modHist[itr->second] += 1;
      }
    }

    map<int, int>::iterator itrH = modHist.begin();
    map<int, int>::iterator itrHEnd = modHist.end();

    //for (; itrH != itrHEnd; itrH++) {
    //  DEBUG_MSG("modHist: " << itrH->first << "  " << itrH->second);
    //}
#endif

    map<float, float>::iterator itr = modFreqs.begin();
    map<float, float>::iterator itrEnd = modFreqs.end();
    for (; itr != itrEnd; itr++) {
      //DEBUG_MSG("modCount: " << itr->first << "  " << itr->second);
      if (itr->second <= 1 || itr->first > 90.0) {
        itr->second = 0;
      } else {
        //DEBUG_MSG("modFreqs: " << itr->first << "  " << itr->second);
        itr->second /= totalCounts;
        //DEBUG_MSG("modFreqs: " << itr->first << "  " << itr->second);
      }
    }
    return true;
  }


} // namespace specnets
