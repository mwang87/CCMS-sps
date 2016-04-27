// Header Include
#include "ExecFilterAligns.h"

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
  ExecFilterAligns::ExecFilterAligns(void)
  {
    m_name = "ExecFilterAligns";
    m_type = "ExecFilterAligns";
  }

  // -------------------------------------------------------------------------
  ExecFilterAligns::ExecFilterAligns(const ParameterList & inputParams) :
    ExecBase(inputParams), m_filteredPairs(0x0), m_ratios(0x0), m_means(0x0), 
        m_varTerms(0x0), ownInput(true), 
        m_idxKept(0x0), m_pvalues(0x0), ownOutput(true)
  {
    m_name = "ExecFilterAligns";
    m_type = "ExecFilterAligns";
  }

  // -------------------------------------------------------------------------
  ExecFilterAligns::ExecFilterAligns(const ParameterList & inputParams,
                                   SpectrumPairSet * filteredPairs,
                                   vector<TwoValues<float> > * ratios,
                                   vector<TwoValues<float> > * means,
                                   vector<float> * varTerms,
                                   vector<unsigned int> * idxKept,
                                   vector<TwoValues<float> > * pvalues) :
    ExecBase(inputParams), m_filteredPairs(filteredPairs), m_ratios(ratios), 
        m_means(means), m_varTerms(varTerms), ownInput(true), 
        m_idxKept(idxKept), m_pvalues(pvalues), ownOutput(true)
  {
    m_name = "ExecFilterAligns";
    m_type = "ExecFilterAligns";
  }

  // -------------------------------------------------------------------------
  ExecFilterAligns::~ExecFilterAligns(void)
  {
  }

  // -------------------------------------------------------------------------
  ExecBase * ExecFilterAligns::clone(const ParameterList & inputParams) const
  {
    return new ExecFilterAligns(inputParams);
  }

  // -------------------------------------------------------------------------
  bool ExecFilterAligns::invoke(void)
  {
    DEBUG_VAR(m_filteredPairs->size());

    if (m_idxKept == 0x0)
    {
      ownOutput = true;
      m_idxKept = new vector<unsigned int> ();
      m_pvalues = new vector<TwoValues<float> > ();
    }

    float minPValue = m_params.getValueDouble("MAX_PVALUE");
    float minRatio = m_params.getValueDouble("MIN_RATIO", 0);
    float pmTol = m_params.getValueDouble("TOLERANCE_PM", 1);
    bool filterTrigs = m_params.getValueBool("FILTER_TRIGS", true);

    DEBUG_VAR(minPValue);
    DEBUG_VAR(minRatio);
    DEBUG_VAR(pmTol);
    DEBUG_VAR(filterTrigs);

    // Convert varTerms to standard deviations
    for (int pivot = 0; pivot < m_means->size(); pivot++)
    {
      (*m_varTerms)[pivot] = sqrt((*m_varTerms)[pivot] - (*m_means)[pivot][0]
          * (*m_means)[pivot][0]);
    }

    DEBUG_VAR(m_means->size());
    DEBUG_VAR(m_varTerms->size());
    DEBUG_VAR(m_ratios->size());

    unsigned int numKept = filterAligns(*m_filteredPairs,
                                        *m_idxKept,
                                        *m_pvalues,
                                        *m_means,
                                        *m_varTerms,
                                        *m_ratios,
                                        minPValue,
                                        minRatio,
                                        pmTol,
                                        filterTrigs);

    DEBUG_VAR(numKept);
    DEBUG_VAR(m_filteredPairs->size());
    DEBUG_VAR(m_pvalues->size());
    DEBUG_VAR(m_idxKept->size());

    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecFilterAligns::loadInputData(void)
  {
    if (m_filteredPairs == 0)
    {
      ownInput = true;
      m_filteredPairs = new SpectrumPairSet();
      m_ratios = new vector<TwoValues<float> > ();
      m_means = new vector<TwoValues<float> > ();
      m_varTerms = new vector<float> ();
    } 
  
    if (m_params.exists("INPUT_ALIGNS"))
    {
      if (!m_filteredPairs->loadFromBinaryFile(m_params.getValue("INPUT_ALIGNS")))
      {
        ERROR_MSG("Could not load: " << m_params.getValue("INPUT_ALIGNS"));
        return false;
      }
    }

    if (m_params.exists("INPUT_RATIOS"))
    {
      if (!Load_binArray(m_params.getValue("INPUT_RATIOS").c_str(), *m_ratios))
      {
        ERROR_MSG("Could not load: " << m_params.getValue("INPUT_RATIOS"));
        return false;
      }
      DEBUG_VAR(m_ratios->size());
    }

    if (m_params.exists("INPUT_MEANS"))
    {
      if (!Load_binArray(m_params.getValue("INPUT_MEANS").c_str(), *m_means))
      {
        ERROR_MSG("Could not load: " << m_params.getValue("INPUT_MEANS"));
        return false;
      }
      DEBUG_VAR(m_means->size());
    }

    if (m_params.exists("INPUT_VARIANCE"))
    {
      if (!Load_binArray(m_params.getValue("INPUT_VARIANCE").c_str(), *m_varTerms))
      {
        ERROR_MSG("Could not load: " << m_params.getValue("INPUT_VARIANCE"));
        return false;
      }
      DEBUG_VAR(m_varTerms->size());
    }
    
    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecFilterAligns::saveOutputData(void)
  {
    if (m_params.exists("OUTPUT_ALIGNS"))
    {
      if (!m_filteredPairs->saveToBinaryFile(m_params.getValue("OUTPUT_ALIGNS")))
      {
        ERROR_MSG("Could not save: " << m_params.getValue("OUTPUT_ALIGNS"));
        return false;
      }
    }
    if (m_params.exists("OUTPUT_PVALUES"))
    {
      if (!Save_binArray(m_params.getValue("OUTPUT_PVALUES").c_str(), *m_pvalues))
      {
        ERROR_MSG("Could not save: " << m_params.getValue("OUTPUT_PVALUES"));
        return false;
      }
    }
    if (m_params.exists("OUTPUT_INDICES"))
    {
      if (!Save_binArray(m_params.getValue("OUTPUT_INDICES").c_str(), *m_idxKept))
      {
        ERROR_MSG("Could not save: " << m_params.getValue("OUTPUT_INDICES"));
        return false;
      }
    }

    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecFilterAligns::saveInputData(std::vector<std::string> & filenames)
  {
    string alignsFilename = m_params.getValue("INPUT_ALIGNS");
    if (!fileExists(alignsFilename))
    {
      m_filteredPairs->saveToBinaryFile(alignsFilename);
      DEBUG_MSG("Saving " << alignsFilename);
    }
    else
    {
      DEBUG_MSG("Not Saving " << alignsFilename << " (already exists)");
    }
    m_params.setValue("INPUT_ALIGNS", alignsFilename);

    string baseFilename = getName();
    string paramFilename = baseFilename + ".params";
    m_params.writeToFile(paramFilename);

    filenames.push_back(paramFilename); // Parameter file MUST be first in vector
    filenames.push_back(alignsFilename);

    return false;
  }

  // -------------------------------------------------------------------------
  bool ExecFilterAligns::loadOutputData(void)
  {
    if (m_idxKept == 0x0)
    {
      ownOutput = true;
      m_idxKept = new vector<unsigned int> ();
      m_pvalues = new vector<TwoValues<float> > ();
    }
    
    if (m_params.exists("OUTPUT_ALIGNS"))
    {
      if (!m_filteredPairs->loadFromBinaryFile(m_params.getValue("OUTPUT_ALIGNS")))
      {
        ERROR_MSG("Could not load: " << m_params.getValue("OUTPUT_ALIGNS"));
        return false;
      }
    }
    if (m_params.exists("OUTPUT_PVALUES"))
    {
      if (!Load_binArray(m_params.getValue("OUTPUT_PVALUES").c_str(), *m_pvalues))
      {
        ERROR_MSG("Could not load: " << m_params.getValue("OUTPUT_PVALUES"));
        return false;
      }
    }
    if (m_params.exists("OUTPUT_INDICES"))
    {
      if (!Load_binArray(m_params.getValue("OUTPUT_INDICES").c_str(), *m_idxKept))
      {
        ERROR_MSG("Could not load: " << m_params.getValue("OUTPUT_INDICES"));
        return false;
      }
    }

    return true;
  }

  // -------------------------------------------------------------------------
  vector<ExecBase*> const & ExecFilterAligns::split(int numSplit)
  {
    return m_subModules;
  }

  // -------------------------------------------------------------------------
  bool ExecFilterAligns::merge(void)
  {
    return false;
  }

  // -------------------------------------------------------------------------
  bool ExecFilterAligns::validateParams(std::string & error)
  {
    m_isValid = false;

    VALIDATE_PARAM_EXIST("MAX_PVALUE");
    VALIDATE_PARAM_EXIST("MIN_RATIO");
    VALIDATE_PARAM_EXIST("TOLERANCE_PM");
    VALIDATE_PARAM_EXIST("FILTER_TRIGS");
    
    m_isValid = true;
    return true;
  }

} // namespace specnets

