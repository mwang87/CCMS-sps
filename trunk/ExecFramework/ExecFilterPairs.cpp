// Header Include
#include "ExecFilterPairs.h"

// Module Includes
#include "AlignmentUtils.h"
#include "Logger.h"
#include "Filters.h"
#include "FileUtils.h"

// External Includes
#include "utils.h"  // for Save_binArray only
// System Includes
#include "stdlib.h"

using namespace specnets;
using namespace std;

namespace specnets
{

  // -------------------------------------------------------------------------
  ExecFilterPairs::ExecFilterPairs(void) :
      ExecBase(), m_inputSpectra(0x0), m_inputSpectraMS2(0x0),
          m_inputSpectraMS2_normalized(0x0), ownInput(true),
          m_filteredPairs(0x0), m_ratios(0x0), m_means(0x0), m_varTerms(0x0),
          m_alignStats(0x0), m_specStats(0x0), ownOutput(true)
  {
    m_name = "ExecFilterPairs";
    m_type = "ExecFilterPairs";

    m_inputSpectra = new SpecSet;
    m_inputSpectraMS2 = new SpecSet;
    m_inputSpectraMS2_normalized = new SpecSet;

    m_filteredPairs = new SpectrumPairSet;
    m_ratios = new std::vector<TwoValues<float> >;
    m_means = new std::vector<TwoValues<float> >;
    m_varTerms = new std::vector<float>;
    m_alignStats = new list<vector<float> >;
    m_specStats = new std::vector<vector<float> >;
  }

  // -------------------------------------------------------------------------
  ExecFilterPairs::ExecFilterPairs(const ParameterList & inputParams) :
      ExecBase(inputParams), m_inputSpectra(0x0), m_inputSpectraMS2(0x0),
          m_inputSpectraMS2_normalized(0x0), ownInput(true),
          m_filteredPairs(0x0), m_ratios(0x0), m_means(0x0), m_varTerms(0x0),
          m_alignStats(0x0), m_specStats(0x0), ownOutput(true)
  {
    m_name = "ExecFilterPairs";
    m_type = "ExecFilterPairs";

    m_inputSpectra = new SpecSet;
    m_inputSpectraMS2 = new SpecSet;
    m_inputSpectraMS2_normalized = new SpecSet;

    m_filteredPairs = new SpectrumPairSet;
    m_ratios = new std::vector<TwoValues<float> >;
    m_means = new std::vector<TwoValues<float> >;
    m_varTerms = new std::vector<float>;
    m_alignStats = new list<vector<float> >;
    m_specStats = new std::vector<vector<float> >;
  }

  ExecFilterPairs::ExecFilterPairs(const ParameterList & inputParams,
                                   SpecSet * inputSpectra,
                                   SpecSet * inputSpectraMS2,
                                   SpecSet * inputSpectraMS2_normed) :
      ExecBase(inputParams), m_inputSpectra(inputSpectra),
          m_inputSpectraMS2(inputSpectraMS2),
          m_inputSpectraMS2_normalized(inputSpectraMS2_normed), ownInput(false),
          m_filteredPairs(0x0), m_ratios(0x0), m_means(0x0), m_varTerms(0x0),
          m_alignStats(0x0), m_specStats(0x0), ownOutput(true)
  {
    m_name = "ExecFilterPairs";
    m_type = "ExecFilterPairs";

    m_filteredPairs = new SpectrumPairSet;
    m_ratios = new std::vector<TwoValues<float> >;
    m_means = new std::vector<TwoValues<float> >;
    m_varTerms = new std::vector<float>;
    m_alignStats = new list<vector<float> >;
    m_specStats = new std::vector<vector<float> >;
  }

  // -------------------------------------------------------------------------
  ExecFilterPairs::ExecFilterPairs(const ParameterList & inputParams,
                                   SpecSet * inputSpectra,
                                   SpecSet * inputSpectraMS2,
                                   SpecSet * inputSpectraMS2_normed,
                                   SpectrumPairSet * filteredPairs,
                                   vector<TwoValues<float> > * ratios,
                                   vector<TwoValues<float> > * means,
                                   vector<float> * varTerms,
                                   list<vector<float> > * alignStats,
                                   vector<vector<float> > * specStats,
                                   vector<unsigned int> * idxKept,
                                   vector<TwoValues<float> > * pvalues) :
      ExecBase(inputParams), m_inputSpectra(inputSpectra),
          m_inputSpectraMS2(inputSpectraMS2),
          m_inputSpectraMS2_normalized(inputSpectraMS2_normed), ownInput(false),
          m_filteredPairs(filteredPairs), m_ratios(ratios), m_means(means),
          m_varTerms(varTerms), m_alignStats(alignStats),
          m_specStats(specStats), ownOutput(false)
  {
    m_name = "ExecFilterPairs";
    m_type = "ExecFilterPairs";

  }

  // -------------------------------------------------------------------------
  ExecFilterPairs::~ExecFilterPairs(void)
  {
    if (ownInput)
    {
      delete m_inputSpectra;
      delete m_inputSpectraMS2;
      delete m_inputSpectraMS2_normalized;
    }

    if (ownOutput)
    {
      delete m_filteredPairs;
      delete m_ratios;
      delete m_means;
      delete m_varTerms;
      delete m_alignStats;
      delete m_specStats;
    }
  }

  // -------------------------------------------------------------------------
  ExecBase * ExecFilterPairs::clone(const ParameterList & inputParams) const
  {
    return new ExecFilterPairs(inputParams);
  }

  // -------------------------------------------------------------------------
  bool ExecFilterPairs::invoke(void)
  {

    DEBUG_VAR(m_params.getValueDouble("MAX_PVALUE"));

    int startBaseIdx;
    if (m_params.exists("IDX_START"))
    {
      startBaseIdx = max(0, m_params.getValueInt("IDX_START"));
    }
    else
    {
      startBaseIdx = 0;
    }
    DEBUG_TRACE;
    int endBaseIdx;
    if (m_params.exists("IDX_END"))
    {
      endBaseIdx = max(0, m_params.getValueInt("IDX_END"));
    }
    else
    {
      endBaseIdx = -1;
    }

    DEBUG_TRACE;
    float minOverlap = 0;
    int aaDiff = m_params.getValueInt("AA_DIFF_COUNT");
    double minShift = m_params.getValueDouble("MIN_SHIFT");
    double maxShift = m_params.getValueDouble("MAX_SHIFT");

    DEBUG_TRACE;
    bool alignPA = m_params.getValueBool("PARTIAL_OVERLAPS", false);
    if (alignPA)
    {
      minOverlap = m_params.getValueDouble("MIN_OVERLAP_AREA", 0);
    }

    DEBUG_TRACE;
    short minNumMatchedPeaks = (short)m_params.getValueInt("MIN_MATCHED_PEAKS",
                                                           6);
    float minRatio = m_params.getValueDouble("MIN_RATIO", 0);
    float peakTol = m_params.getValueDouble("TOLERANCE_PEAK", 0.5);
    float pmTol = m_params.getValueDouble("TOLERANCE_PM", 1);
    float resolution = m_params.getValueDouble("RESOLUTION", 0.1);
    int useMinDist = m_params.getValueInt("USE_MIN_DIST_57", 1);
    bool specTypeMSMS = m_params.getValueBool("SPEC_TYPE_MSMS", false);
    float symmetryOffset = specTypeMSMS ? 2 * AAJumps::massHion : 0;
    float ionOffset = specTypeMSMS ? AAJumps::massHion : 0;

    DEBUG_TRACE;
    vector<vector<TTag> > tags(0);
    int tagsMatchFlank = m_params.getValueInt("TAGS_MATCH_FLANK", 1);
    int tagsMatchCount = m_params.getValueInt("TAGS_MATCH_COUNT", 1);

    DEBUG_TRACE;
    if (m_params.exists("INPUT_TAGS"))
    {
      for (int i = 0; i <= 3; i++)
      {
        stringstream s;
        s << i;
        LoadTags(string(m_params.getValue("INPUT_TAGS")) + s.str() + ".txt",
                 tags,
                 tagsMatchFlank);
      }
    }

    DEBUG_TRACE;
    if (endBaseIdx < 0)
    {
      endBaseIdx = m_inputSpectra->size() - 1;
    }

    DEBUG_VAR(m_inputSpectra->size());
    DEBUG_VAR(startBaseIdx);
    DEBUG_VAR(endBaseIdx);

    DEBUG_VAR(aaDiff);
    DEBUG_VAR(minShift);
    DEBUG_VAR(maxShift);
    DEBUG_VAR(pmTol);
    DEBUG_VAR(peakTol);
    DEBUG_VAR(minRatio);
    DEBUG_VAR(minNumMatchedPeaks);
    DEBUG_VAR(m_inputSpectraMS2_normalized);
    if (m_inputSpectraMS2_normalized != 0x0)
      DEBUG_VAR(m_inputSpectraMS2_normalized->size());

    if (startBaseIdx < 0 || startBaseIdx >= m_inputSpectra->size())
    {
      DEBUG_MSG("Invalid start index [" << startBaseIdx << "] (" << m_inputSpectra->size() << " spectra)");
      return false;
    }

    if (endBaseIdx >= m_inputSpectra->size())
    {
      WARN_MSG("IDX_END (" << endBaseIdx << ") was lowered to the index of the last spectrum in the dataset (" << m_inputSpectra->size() - 1 << ")");
    }
    endBaseIdx = min(m_inputSpectra->size() - 1, (unsigned int)endBaseIdx);

    DEBUG_MSG("Computing alignments between spectra " << startBaseIdx << ":" << endBaseIdx << " and all others in the dataset.");

    vector<int> baseSpecIdx(endBaseIdx - startBaseIdx + 1);
    for (int i = 0; i < endBaseIdx - startBaseIdx + 1; i++)
    {
      baseSpecIdx[i] = startBaseIdx + i;
    }

    if (!alignPA)
    {
      DEBUG_TRACE;

      DEBUG_VAR(m_params.getValue("PAIRS_MATCH_MODE",""));
      if(m_params.getValue("PAIRS_MATCH_MODE","") == "cosine")
      {
        float minCosine = (float) m_params.getValueDouble("PAIRS_MIN_COSINE", 0.0);
        DEBUG_VAR(minCosine);

        if (!m_inputSpectraMS2_normalized || m_inputSpectraMS2_normalized->size()==0)
        {
          DEBUG_MSG("No MS2 spectra for cosine calculations (INPUT_SPECS_MS2=" << m_params.getValue("INPUT_SPECS_MS2")<<")");
          return false;
        }

        //(*m_inputSpectraMS2)[0].output_ms2(cout);

        getPairCosines(*m_inputSpectraMS2_normalized,
            baseSpecIdx,
            maxShift,// maximum parent mass difference
            pmTol,
            peakTol,
            minCosine,
            minRatio,
            minNumMatchedPeaks,
            *m_filteredPairs,
            *m_ratios,
            *m_means,
            *m_varTerms,
            resolution,
            symmetryOffset);
      }
      else
      {
        getPairAlignsASP(*m_inputSpectra,
            baseSpecIdx,
            aaDiff,
            minShift,
            maxShift,
            pmTol,
            peakTol,
            minRatio,
            minNumMatchedPeaks,
            *m_filteredPairs,
            *m_ratios,
            *m_means,
            *m_varTerms,
            tags,
            tagsMatchFlank,
            tagsMatchCount,
            resolution,
            symmetryOffset);
      }
    }
    else
    {
      DEBUG_TRACE;
      m_alignStats->clear();
      m_specStats->resize(0);

      list<TwoValues<int> > numMatchedPeaks;
      AAJumps jumps(0); //jumps.alljumps(1000.0,resolution);

      DEBUG_TRACE;
      getPairAlignsPA2(*m_inputSpectra,
          startBaseIdx,
          endBaseIdx,
          peakTol,
          pmTol,
          minRatio,
          minOverlap,
          minNumMatchedPeaks,
          jumps,
          minShift,
          *m_filteredPairs,
          *m_ratios,
          numMatchedPeaks,
          *m_means,
          *m_varTerms,
          *m_alignStats,
          *m_specStats);
    }

    DEBUG_VAR(m_filteredPairs->size());

    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecFilterPairs::loadInputData(void)
  {

    DEBUG_TRACE;
    if (m_params.exists("INPUT_SPECS"))
    {
      if (!m_inputSpectra->LoadSpecSet_pkl(m_params.getValue("INPUT_SPECS").c_str()))
      {
        DEBUG_MSG("Could not load " << m_params.getValue("INPUT_SPECS"));
        return false;
      }
    }
    else if (m_params.exists("INPUT_SPECS_PKLBIN"))
    {
      if (!m_inputSpectra->loadPklBin(m_params.getValue("INPUT_SPECS_PKLBIN").c_str()))
      {
        DEBUG_MSG("Could not load " << m_params.getValue("INPUT_SPECS_PKLBIN"));
        return false;
      }
    }

    DEBUG_TRACE;
    if (m_params.exists("INPUT_SPECS_MS2"))
    {
      if (!m_inputSpectraMS2->loadPklBin(m_params.getValue("INPUT_SPECS_MS2").c_str()))
      {
        DEBUG_MSG("Could not load " << m_params.getValue("INPUT_SPECS_MS2"));
        return false;
      }
    }

    DEBUG_TRACE;
    if (m_inputSpectra->size() == 0 && m_inputSpectraMS2->size()==0)
    {
      DEBUG_MSG("Input spectra size is 0!");
      return false;
    }

    DEBUG_VAR(m_inputSpectra->size());

    // Load amino acid masses
    AAJumps jumps(1);
    if (m_params.exists("AMINO_ACID_MASSES"))
    {
      if (!jumps.loadJumps(m_params.getValue("AMINO_ACID_MASSES").c_str(), true))
      {
        DEBUG_MSG("Could not load " << m_params.getValue("AMINO_ACID_MASSES"));
      }
      DEBUG_MSG("Loaded amino acids from " << m_params.getValue("AMINO_ACID_MASSES"));
    }

    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecFilterPairs::saveOutputData(void)
  {
    if (m_params.exists("OUTPUT_ALIGNS"))
    {
      DEBUG_VAR(m_filteredPairs->size());
      m_filteredPairs->saveToBinaryFile(m_params.getValue("OUTPUT_ALIGNS"));
    }
    if (m_params.exists("OUTPUT_RATIOS"))
    {
      DEBUG_VAR(m_ratios->size());
      DEBUG_MSG("Outputting ratios...");
      Save_binArray(m_params.getValue("OUTPUT_RATIOS").c_str(), *m_ratios);
      DEBUG_MSG("done...");
    }
    if (m_params.exists("OUTPUT_MEANS"))
    {
      DEBUG_VAR(m_means->size());
      DEBUG_MSG("Outputting means...");
      Save_binArray(m_params.getValue("OUTPUT_MEANS").c_str(), *m_means);
      DEBUG_MSG("done...");
    }
    if (m_params.exists("OUTPUT_VARIANCE"))
    {
      DEBUG_VAR(m_varTerms->size());
      DEBUG_MSG("Outputting variance terms...");
      Save_binArray(m_params.getValue("OUTPUT_VARIANCE").c_str(), *m_varTerms);
      DEBUG_MSG("done...");
    }

    string statsFilename = m_params.getValue("OUTPUT_STATS");
    if (m_params.exists("OUTPUT_STATS"))
    {
      string specstatsFilename(statsFilename);
      specstatsFilename += "_specStats.bin";
      Save_binListArray<float, vector<float>, vector<float>::iterator>(specstatsFilename.c_str(),
                                                                       *m_specStats);

      string alignstatsFilename(statsFilename);
      alignstatsFilename += "_alignStats.bin";
      Save_binListArray<float, vector<float>, vector<float>::iterator>(alignstatsFilename.c_str(),
                                                                       *m_alignStats);
    }

    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecFilterPairs::saveInputData(std::vector<std::string> & filenames)
  {
    string dataDir = m_params.getValue("GRID_DATA_DIR");
    if (dataDir.empty())
    {
      dataDir = ".";
    }
    string baseDirectory = dataDir + "/";
    string baseFilename = baseDirectory + getName();

    string spectraFilename = baseDirectory + "specs_scored.pklbin";
    if (!fileExists(spectraFilename))
    {
      m_inputSpectra->savePklBin(spectraFilename.c_str());
      DEBUG_MSG("Saving " << spectraFilename);
    }
    else
    {
      DEBUG_MSG("Not Saving " << spectraFilename << " (already exists)");
    }
    m_params.setValue("INPUT_SPECS_PKLBIN", spectraFilename);

    //SpecSet m_inputSpectra; // the input spectra
    string paramFilename = baseFilename + ".params";
    m_params.writeToFile(paramFilename);

    filenames.push_back(paramFilename); // Parameter file MUST be first in vector
    filenames.push_back(spectraFilename);

    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecFilterPairs::loadOutputData(void)
  {
    if (m_params.exists("OUTPUT_ALIGNS"))
    {
      DEBUG_MSG("Loading alignments...");
      DEBUG_MSG(m_params.getValue("OUTPUT_ALIGNS"));
      if (m_filteredPairs->loadFromBinaryFile(m_params.getValue("OUTPUT_ALIGNS"))
          < 0)
      {
        ERROR_MSG("Could not load: " << m_params.getValue("OUTPUT_ALIGNS"));
        return false;
      }
      DEBUG_MSG("Loaded " << m_filteredPairs->size() << " alignments");
    }

    if (m_params.exists("OUTPUT_RATIOS"))
    {
      DEBUG_MSG("Loading ratios...");
      vector<vector<float> > ratios;
      if (Load_binArray(m_params.getValue("OUTPUT_RATIOS").c_str(), ratios) < 0)
      {
        ERROR_MSG("Could not load: " << m_params.getValue("OUTPUT_RATIOS"));
        return false;
      }
      m_ratios->resize(ratios.size());
      for (int i = 0; i < ratios.size(); i++)
      {
        (*m_ratios)[i][0] = ratios[i][0];
        (*m_ratios)[i][1] = ratios[i][1];
        //DEBUG_MSG("Ratio " << i << " " << (*m_ratios)[i][0] << " " << (*m_ratios)[i][1])
      }
      m_ratios->resize(ratios.size());
      DEBUG_MSG("Loaded " << m_ratios->size() << " ratios");
    }

    if (m_params.exists("OUTPUT_MEANS"))
    {
      DEBUG_MSG("Loading means...");
      vector<vector<float> > means;
      if (Load_binArray(m_params.getValue("OUTPUT_MEANS").c_str(), means) < 0)
      {
        ERROR_MSG("Could not load: " << m_params.getValue("OUTPUT_MEANS"));
        return false;
      }
      m_means->resize(means.size());
      for (int i = 0; i < means.size(); i++)
      {
        (*m_means)[i][0] = means[i][0];
        (*m_means)[i][1] = means[i][1];
        //DEBUG_MSG("Mean " << i << " " << (*m_means)[i][0] << " " << (*m_means)[i][1])
      }
      m_means->resize(means.size());
      DEBUG_MSG("Loaded " << m_means->size() << " means");
    }

    if (m_params.exists("OUTPUT_VARIANCE"))
    {
      DEBUG_MSG("Loading means...");
      vector<vector<float> > varTerms;
      if (Load_binArray(m_params.getValue("OUTPUT_VARIANCE").c_str(), varTerms)
          < 0)
      {
        ERROR_MSG("Could not load: " << m_params.getValue("OUTPUT_VARIANCE"));
        return false;
      }
      m_varTerms->resize(varTerms.size());
      for (int i = 0; i < varTerms.size(); i++)
      {
        (*m_varTerms)[i] = varTerms[i][0];
        //DEBUG_MSG("Variance " << i << " " << (*m_varTerms)[i])
      }
      m_varTerms->resize(varTerms.size());
      DEBUG_MSG("Loaded " << m_varTerms->size() << " varTerms");
    }

    string statsFilename = m_params.getValue("OUTPUT_STATS");
    if (m_params.exists("OUTPUT_STATS"))
    {
      string specstatsFilename(statsFilename);
      specstatsFilename += "_specStats.bin";
      Load_binListArray<float, vector<float>, vector<float>::iterator>(specstatsFilename.c_str(),
                                                                       *m_specStats);

      string alignstatsFilename(statsFilename);
      alignstatsFilename += "_alignStats.bin";
      Load_binListArray<float, vector<float>, list<vector<float> >::iterator>(alignstatsFilename.c_str(),
                                                                              *m_alignStats);
    }

    return true;
  }

  // ------------------------------------------------------------------------

// -------------------------------------------------------------------------

  vector<ExecBase*> const & ExecFilterPairs::split(int numSplit)
  {
    DEBUG_VAR(numSplit);

    if (numSplit < 2)
    {
      DEBUG_MSG("Number split [" << numSplit << "] must be at least 2");
      return m_subModules;
    }

    int spectraSize = m_inputSpectra->size();
    DEBUG_VAR(spectraSize);
    if (spectraSize == 0)
    {
      DEBUG_MSG("Must have at least one spectrum");
      return m_subModules;
    }

    int startBaseIdx;
    if (m_params.exists("IDX_START"))
    {
      startBaseIdx = max(0, m_params.getValueInt("IDX_START"));
    }
    else
    {
      startBaseIdx = 0;
    }
    DEBUG_VAR(startBaseIdx);
    int endBaseIdx;
    if (m_params.exists("IDX_END"))
    {
      endBaseIdx = max(0, m_params.getValueInt("IDX_END"));
      endBaseIdx = min(endBaseIdx, (int)m_inputSpectra->size() - 1);
    }
    else
    {
      endBaseIdx = m_inputSpectra->size() - 1;
    }
    DEBUG_VAR(endBaseIdx);
    if (startBaseIdx == endBaseIdx)
    {
      DEBUG_MSG("Must have more than one row");
      return m_subModules;
    }

    DEBUG_TRACE;

    int startIndex = 0;
    double runningCount = 0.0;
    int nSplit = 0;
    int remainingSplits = numSplit;

    // Keeps track of the computation cost per spectrum
    vector<double> costPerSpectrum(endBaseIdx - startBaseIdx + 1, 0);
    double comparisonsPerCpu = 0;
    for (int i = startBaseIdx; i <= endBaseIdx; i++)
    {
      double totPeaks = 0;
      // every spectrum is aligned to all others
      for (int j = i + 1; j < m_inputSpectra->size(); j++)
      {
        totPeaks += (*m_inputSpectra)[j].size();
      }
      // cost for every pair is defined as numPeaks_spec1 * numPeaks_spec2
      costPerSpectrum[i - startBaseIdx] = ((double)(*m_inputSpectra)[i].size())
          * totPeaks;
      comparisonsPerCpu += costPerSpectrum[i - startBaseIdx];
    }
    comparisonsPerCpu = comparisonsPerCpu / ((double)numSplit);

    DEBUG_VAR(comparisonsPerCpu);
    DEBUG_VAR(m_params.getValueBool("PARTIAL_OVERLAPS", false));
    // Note: "i" is startIdx for current CPU - it should never be equal to endBaseIdx
    for (int i = startBaseIdx; i < endBaseIdx; i++)
    {
      //DEBUG_VAR(i);
      if (nSplit == numSplit)
      {
        ERROR_MSG(i);
        ERROR_MSG(runningCount);
        ERROR_MSG(nSplit);
        ERROR_MSG("Exceeded number of available jobs while still allocating alignments for spectrum "<<i);
        exit(-1);
      }

      double comparisonsThisSpectra = costPerSpectrum[i - startBaseIdx]; //spectraSize - i - 1;
//      DEBUG_VAR(comparisonsThisSpectra);
      runningCount += comparisonsThisSpectra;
//      DEBUG_VAR(runningCount);
//      DEBUG_VAR(comparisonsPerCpu);

      // When there is one comparison left.. its the last comparison needed
      //  Keep in mind that the last spectrum doesn't get compared to anyone else
      //  Everyone has already been compared to it
      if (runningCount >= comparisonsPerCpu || i == endBaseIdx - 1)
      {
        // We have enough comparisons for the CPU
        // Copy the parameters
        ParameterList childParams(m_params);
        // But set the start and end indices
        char buf[128];
        sprintf(buf, "%d", startIndex);
        childParams.setValue("IDX_START", buf);
        sprintf(buf, "%d", i);
        childParams.setValue("IDX_END", buf);

        //DEBUG_MSG("Start [" << startIndex << "] End [" << i << "] Split [" << numSplit << "]");

        //ExecBase * theClone = new ExecFilterPairs(childParams,
        //                                          m_inputSpectra);

        ExecBase * theClone = new ExecFilterPairs(childParams,
                                                  m_inputSpectra,
                                                  m_inputSpectraMS2,
                                                  m_inputSpectraMS2_normalized);

        theClone->setName(makeName(m_name, nSplit));

        // Have to set up the output files also so the params will be correct on reload
        string dataDir = m_params.getValue("GRID_DATA_DIR");
        if (dataDir.empty())
        {
          dataDir = ".";
        }
        string baseName = dataDir + "/" + theClone->getName();
        //DEBUG_VAR(baseName);

        string jobsFilename = baseName + "_means.bin";
        theClone->m_params.setValue("OUTPUT_ALIGNS",
                                    baseName + "_aligns" + ".bin");
        theClone->m_params.setValue("OUTPUT_MEANS",
                                    baseName + "_means_" + ".bin");
        theClone->m_params.setValue("OUTPUT_VARIANCE",
                                    baseName + "_vars_" + ".bin");
        theClone->m_params.setValue("OUTPUT_RATIOS",
                                    baseName + "_ratios_" + ".bin");

        std::string suffix("");
        char bufSplit[128];
        sprintf(bufSplit, "%d", nSplit + 1);
        theClone->m_params.setValue("NUM_SPLIT", bufSplit);

        m_subModules.push_back(theClone);

        // Reset start index and running counts
        //DEBUG_VAR(runningCount);

        // Recompute the loads since we used up "extra" on this node
        remainingSplits--;
//        DEBUG_VAR(remainingSplits);
        if (remainingSplits == 0)
        {
          comparisonsPerCpu = 0;
        }
        else
        {
          comparisonsPerCpu = 0;
          for (int j = i + 1; j <= endBaseIdx; j++)
          {
            comparisonsPerCpu += costPerSpectrum[j - startBaseIdx];
          }
          comparisonsPerCpu = comparisonsPerCpu / ((double)remainingSplits);
          DEBUG_VAR(comparisonsPerCpu);
        }

        runningCount = 0;
        startIndex = i + 1;
        nSplit++;

      } // if (runningCount >= comparisonsPerCpu || i == spectraSize - 1)

    } // for (int i = 0; i <= spectraSize - 1; i++)

    DEBUG_VAR(m_subModules.size());
    return m_subModules;
  }

  // -------------------------------------------------------------------------
  bool ExecFilterPairs::merge(void)
  {

    DEBUG_TRACE;

    unsigned int startIdx = (unsigned int)m_params.getValueInt("IDX_START", 0);
    unsigned int endIdx = (unsigned int)m_params.getValueInt("IDX_END", 0);
    bool alignPA = m_params.getValueBool("PARTIAL_OVERLAPS", false);
    float pmTol = m_params.getValueDouble("TOLERANCE_PM", 1);
    float minRatio = m_params.getValueDouble("MIN_RATIO", 0);
    float minPValue = m_params.getValueDouble("MAX_PVALUE", 0.05);

    // Find out how many total pairs from all the children
    DEBUG_VAR(m_subModules.size());
    int totalPairs = 0;
    for (int i = 0; i < m_subModules.size(); i++)
    {
      ExecFilterPairs * fpe = (ExecFilterPairs *)m_subModules[i];
      DEBUG_VAR(fpe);
      DEBUG_VAR(fpe->m_filteredPairs->size());
      totalPairs += fpe->m_filteredPairs->size();
    }

    DEBUG_VAR(totalPairs);
    if (totalPairs <= 0)
    {
      DEBUG_MSG("Can not merge/filter with " << totalPairs << " input pairs!");
      return false;
    }

    // Resize our array to hold all the pairs
    DEBUG_VAR(m_filteredPairs);
    m_filteredPairs->resize(totalPairs);

    DEBUG_VAR(m_filteredPairs->size());

    // Copy all the result pairs from the children into our array
    int pairCount = 0;
    for (int i = 0; i < m_subModules.size(); i++)
    {
      //DEBUG_VAR(i);
      ExecFilterPairs * fpe = (ExecFilterPairs *)m_subModules[i];
      SpectrumPairSet * sps = fpe->m_filteredPairs;
      for (int k = 0; k < sps->size(); k++)
      {
        //DEBUG_VAR(k);
        (*m_filteredPairs)[pairCount++] = sps->operator[](k);
      }
    }
    DEBUG_MSG("Merged " << totalPairs << " pairs");

    // Read first means/vars
    ExecFilterPairs * fpe = (ExecFilterPairs *)m_subModules[0];
    unsigned int szParams = fpe->m_means->size();
    m_means->resize(szParams);
    m_varTerms->resize(szParams);
    for (int pivot = 0; pivot < szParams; pivot++)
    {
      (*m_means)[pivot][0] = (*fpe->m_means)[pivot].values[0];
      (*m_means)[pivot][1] = (*fpe->m_means)[pivot].values[1];
      (*m_varTerms)[pivot] = (*fpe->m_varTerms)[pivot];
      //DEBUG_MSG("Mean / Var " << pivot << " " << (*m_means)[pivot][0] << " " << (*m_means)[pivot][1] << " " << (*m_varTerms)[pivot]);
    }

    for (int i = 1; i < m_subModules.size(); i++)
    {
      float ratio1, ratio2;

      ExecFilterPairs * fpe = (ExecFilterPairs *)m_subModules[i];
      for (int pivot = 0; pivot < m_inputSpectra->size(); pivot++)
      {
        // Avoid computations when there are no new observations
        if ((*fpe->m_means)[pivot][1] > 0.1)
        {
          ratio1 = (*m_means)[pivot][1]
          / ((*m_means)[pivot][1] + (*fpe->m_means)[pivot].values[1]);

          ratio2 = 1 - ratio1;

          (*m_means)[pivot][0] = (*m_means)[pivot][0] * ratio1
          + (*fpe->m_means)[pivot].values[0] * ratio2;

          (*m_varTerms)[pivot] = (*m_varTerms)[pivot] * ratio1
          + (*fpe->m_varTerms)[pivot] * ratio2;

          (*m_means)[pivot][1] += (*fpe->m_means)[pivot][1];

          //DEBUG_MSG("Mean / Var " << pivot << " " << (*m_means)[pivot][0] << " " << (*m_means)[pivot][1] << " " << (*m_varTerms)[pivot]);
        }
      }
    }

    DEBUG_MSG("Done estimating means/stddevs...");

    // Find out how many total ratios from all the children
    int totalRatios = 0;
    for (int i = 0; i < m_subModules.size(); i++)
    {
      ExecFilterPairs * fpe = (ExecFilterPairs *)m_subModules[i];
      totalRatios += fpe->m_ratios->size();
    }
    DEBUG_VAR(totalRatios);

    // Resize our array to hold all the pairs
    m_ratios->resize(totalRatios);

    // Copy all the ratios from the children into our array
    int ratioCount = 0;
    for (int i = 0; i < m_subModules.size(); i++)
    {
      ExecFilterPairs * fpe = (ExecFilterPairs *)m_subModules[i];
      vector<TwoValues<float> > * pratios = fpe->m_ratios;
      for (int k = 0; k < fpe->m_ratios->size(); k++)
      {
        (*m_ratios)[ratioCount++] = pratios->operator[](k);
      }
    }

    DEBUG_MSG("Merged " << totalRatios << " ratios");

    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecFilterPairs::validateParams(std::string & error)
  {
    m_isValid = false;

    VALIDATE_PARAM_EXIST("AA_DIFF_COUNT");
    VALIDATE_PARAM_EXIST("MIN_SHIFT");
    VALIDATE_PARAM_EXIST("MAX_SHIFT");
    VALIDATE_PARAM_EXIST("MIN_RATIO");
    VALIDATE_PARAM_EXIST("TOLERANCE_PEAK");
    VALIDATE_PARAM_EXIST("TOLERANCE_PM");
    VALIDATE_PARAM_EXIST("RESOLUTION");
    VALIDATE_PARAM_EXIST("USE_MIN_DIST_57");
    VALIDATE_PARAM_EXIST("SPEC_TYPE_MSMS");
    VALIDATE_PARAM_EXIST("TAGS_MATCH_FLANK");
    VALIDATE_PARAM_EXIST("TAGS_MATCH_COUNT");
    VALIDATE_PARAM_EXIST("MAX_PVALUE");

    int startBaseIdx;
    if (m_params.exists("IDX_START"))
    {
      startBaseIdx = max(0, m_params.getValueInt("IDX_START"));
    }
    else
    {
      startBaseIdx = 0;
    }
    int endBaseIdx;
    if (m_params.exists("IDX_END"))
    {
      endBaseIdx = max(0, m_params.getValueInt("IDX_END"));
    }
    else
    {
      endBaseIdx = -1;
    }

    m_isValid = true;
    return true;
  }

} // namespace specnets

