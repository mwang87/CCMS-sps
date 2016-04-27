#include "FdrPeptide.h"

namespace specnets
{
  //---------------------------------------------------------------------------

  bool compareFDR(psmPtr i, psmPtr j)
  {
    return i->m_score > j->m_score;
  }
  // -------------------------------------------------------------------------
  bool FdrPeptide::concatenatedTargetDecoy(PeptideSpectrumMatchSet &inputPeptides,
                                           PeptideSpectrumMatchSet &outputPeptides,
                                           double scalingFactor,
                                           bool (*fdrCompareFunction)(psmPtr, psmPtr))
  {
    //filter to top hit per scan
    PeptideSpectrumMatchSet temp;
    if (!filterToTopHit(inputPeptides, temp))
    {
      return false;
    }
    if (!calculatePValues(temp, outputPeptides, scalingFactor, fdrCompareFunction))
    {
      return false;
    }
    return true;
  }
  // -------------------------------------------------------------------------
  bool FdrPeptide::separateTargetDecoy(PeptideSpectrumMatchSet &inputPeptides,
                                       PeptideSpectrumMatchSet &outputPeptides,
                                       double scalingFactor,
                                       bool (*fdrCompareFunction)(psmPtr, psmPtr))
  {
    //separate target and decoy hits
    PeptideSpectrumMatchSet target;
    PeptideSpectrumMatchSet decoy;
    for (unsigned int i = 0; i < inputPeptides.size(); i++)
    {
      psmPtr psm = inputPeptides[i];

      if (psm->m_isDecoy)
      {
        decoy.push_back(psm);
      }
      else
      {
        target.push_back(psm);
      }
    }
    PeptideSpectrumMatchSet targetFiltered;
    PeptideSpectrumMatchSet decoyFiltered;

    if (!filterToTopHit(target, targetFiltered))
    {
      return false;
    }

    if (!filterToTopHit(decoy, decoyFiltered))
    {
      return false;
    }

    PeptideSpectrumMatchSet concatenatedFiltered;
    concatenateTargetDecoy(targetFiltered,decoyFiltered,concatenatedFiltered);

    if (!calculatePValues(concatenatedFiltered, outputPeptides, scalingFactor, fdrCompareFunction))
    {
      return false;
    }
    return true;
  }
  // -------------------------------------------------------------------------
  bool FdrPeptide::calculatePValues(PeptideSpectrumMatchSet &inputPeptides,
                                    PeptideSpectrumMatchSet &outputPeptides,
                                    double scalingFactor,
                                    bool (*fdrCompareFunction)(psmPtr, psmPtr))
  {
    outputPeptides = inputPeptides;


    if (fdrCompareFunction == 0) {
      sort(outputPeptides.m_psmSet.begin(), outputPeptides.m_psmSet.end(), compareFDR);
    } else {
      sort(outputPeptides.m_psmSet.begin(), outputPeptides.m_psmSet.end(), fdrCompareFunction);
    }

    unsigned int correctHits = 0;
    unsigned int incorrectHits = 0;

    for (unsigned int i = 0; i < outputPeptides.size(); i++)
    {
      psmPtr psm = outputPeptides[i];
      if (psm->m_isDecoy)
      {
        incorrectHits++;
      }
      else
      {
        correctHits++;
      }
      if (correctHits > 0)
      {
        double FDR = ((double) incorrectHits * (double) scalingFactor) / (double) correctHits;
        if (FDR > 1)
        {
          outputPeptides[i]->m_fdr = 1;
        }
        else
        {
          outputPeptides[i]->m_fdr = FDR;
        }
      }
      else
      {
        outputPeptides[i]->m_fdr = 1;
      }
    }
    
    float min_FDR = 1.0;
    for(int i = outputPeptides.size() - 1; i>= 0; i--){
        DEBUG_MSG(i<<"\t"<<outputPeptides.size());
        min_FDR = min(outputPeptides[i]->m_fdr, min_FDR);
        outputPeptides[i]->m_fdr = min_FDR;
    }


    if (incorrectHits == 0)
    {
      ERROR_MSG("No decoys set on inputPeptides!");
      return false;
    }
    else
    {
      return true;
    }
  }
  // helper function for calculatePValuesSeparateSearch,
  // adds psm with maximum m_score to map.
  // -------------------------------------------------------------------------
  bool addToMaxMap(PeptideSpectrumMatchSet &inputSet, std::tr1::unordered_map<
        string, psmPtr> &maxPsmMap)
    {
      std::tr1::unordered_map<string, psmPtr>::const_iterator it;

      for (unsigned int i = 0; i < inputSet.size(); i++)
      {
        psmPtr currMatch = inputSet[i];
        PeptideSpectrumMatch psm;

        stringstream ss;
        ss << currMatch->m_scanNum << '_' << currMatch->m_spectrumFile;
        string key = ss.str();

        it = maxPsmMap.find(key);

        if (it != maxPsmMap.end())
        {
          if (maxPsmMap[key]->m_score < currMatch->m_score)
          {
            maxPsmMap[key] = currMatch;
          }
        }
        else
        {
          maxPsmMap[key] = currMatch;
        }
      }
      return true;
    }
  // -------------------------------------------------------------------------
  bool FdrPeptide::mergeTargetDecoy(PeptideSpectrumMatchSet &targetPeptides,
                                    PeptideSpectrumMatchSet &decoyPeptides,
                                    PeptideSpectrumMatchSet &outputPeptides)
  {
    std::tr1::unordered_map<string, psmPtr> maxPsmMap; //map between spectrum pointer and current max psm

    //cycle through target peptides first

    if (!addToMaxMap(targetPeptides, maxPsmMap))
    {
      ERROR_MSG("Unable to filter to top hit!");
      return false;
    }

    if (!addToMaxMap(decoyPeptides, maxPsmMap))
    {
      ERROR_MSG("Unable to filter to top hit!");
      return false;
    }

    std::tr1::unordered_map<string, psmPtr>::const_iterator it;

    for (it = maxPsmMap.begin(); it != maxPsmMap.end(); it++)
    {
      psmPtr currPtr = it->second;
      outputPeptides.push_back(currPtr);
    }
    return true;
  }
  // -------------------------------------------------------------------------
  bool FdrPeptide::filterToTopHit(PeptideSpectrumMatchSet &inputPeptides,
                                  PeptideSpectrumMatchSet &outputPeptides)
  {
    std::tr1::unordered_map<string, psmPtr> maxPsmMap; //map between spectrum pointer and current max psm

    if (!addToMaxMap(inputPeptides, maxPsmMap))
    {
      ERROR_MSG("Unable to filter to top hit!");
      return false;
    }

    std::tr1::unordered_map<string, psmPtr>::const_iterator it;

    for (it = maxPsmMap.begin(); it != maxPsmMap.end(); it++)
    {
      psmPtr currPtr = it->second;
      outputPeptides.push_back(currPtr);
    }
    return true;
  }
  // -------------------------------------------------------------------------
  void FdrPeptide::concatenateTargetDecoy(PeptideSpectrumMatchSet &targetPeptides,
                                          PeptideSpectrumMatchSet &decoyPeptides,
                                          PeptideSpectrumMatchSet &outputPeptides)
  {
    for (unsigned int i = 0; i < targetPeptides.size(); i++)
    {
      psmPtr psm = targetPeptides[i];
      outputPeptides.push_back(psm);
    }

    for (unsigned int i = 0; i < decoyPeptides.size(); i++)
    {
      psmPtr psm = decoyPeptides[i];
      outputPeptides.push_back(psm);
    }
  }
  // -------------------------------------------------------------------------
  bool FdrPeptide::filterByPValue(PeptideSpectrumMatchSet &inputPeptides,
                                  double cutoff)
  {
    int count = 0;
    PeptideSpectrumMatchSet filteredSet;
    int max_index = 0;
    for (unsigned int i = 0; i < inputPeptides.size(); i++)
    {
      psmPtr currMatch = inputPeptides[i];

      if (currMatch->m_fdr == -1)
      {
        ERROR_MSG("PValue not set on input!");
        return false;
      }

      if (currMatch->m_fdr <= cutoff)
      {
        max_index = i;
        //filteredSet.m_psmSet.push_back(currMatch);
      }
    }

    for (unsigned int i = 0; i <= max_index; i++)
    {
      psmPtr currMatch = inputPeptides[i];
      filteredSet.m_psmSet.push_back(currMatch);
    }

    inputPeptides = filteredSet;
    return true;
  }
}
