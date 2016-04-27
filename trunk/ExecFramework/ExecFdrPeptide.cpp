// Header Include
#include "ExecFdrPeptide.h"

using namespace std;

namespace specnets
{
  // -------------------------------------------------------------------------
  ExecFdrPeptide::ExecFdrPeptide(void) :
    m_peptideResults(0x0), m_targetResults(0x0), m_decoyResults(0x0),
        m_fdrPeptides(0x0), ownInput(true), ownOutput(true)
  {
    m_name = "ExecFdrPeptide";
    m_type = "ExecFdrPeptide";
  }

  // -------------------------------------------------------------------------
  ExecFdrPeptide::ExecFdrPeptide(const ParameterList & inputParams) :
    ExecBase(inputParams), m_peptideResults(0x0), m_targetResults(0x0),
        m_decoyResults(0x0), m_fdrPeptides(0x0), ownInput(true),
        ownOutput(true)
  {
    m_name = "ExecFdrPeptide";
    m_type = "ExecFdrPeptide";
  }
  // -------------------------------------------------------------------------
  ExecFdrPeptide::ExecFdrPeptide(const ParameterList & inputParams,
                                 PeptideSpectrumMatchSet * peptideResults,
                                 PeptideSpectrumMatchSet * fdrPeptides) :
    ExecBase(inputParams), m_peptideResults(peptideResults),
        m_targetResults(0x0), m_decoyResults(0x0), m_fdrPeptides(fdrPeptides),
        ownInput(false), ownOutput(false)
  {
    m_name = "ExecFdrPeptide";
    m_type = "ExecFdrPeptide";
  }
  // -------------------------------------------------------------------------
  ExecFdrPeptide::ExecFdrPeptide(const ParameterList & inputParams,
                                 PeptideSpectrumMatchSet * targetResults,
                                 PeptideSpectrumMatchSet * decoyResults,
                                 PeptideSpectrumMatchSet * fdrPeptides) :
    ExecBase(inputParams), m_targetResults(targetResults),
        m_peptideResults(0x0), m_decoyResults(decoyResults),
        m_fdrPeptides(fdrPeptides), ownInput(false), ownOutput(false)
  {
    m_name = "ExecFdrPeptide";
    m_type = "ExecFdrPeptide";
  }
  // -------------------------------------------------------------------------
  ExecFdrPeptide::~ExecFdrPeptide(void)
  {
    if (ownInput)
    {
      if (m_peptideResults)
      {
        delete m_peptideResults;
        m_peptideResults = 0x0;
      }
      if (m_targetResults)
      {
        delete m_targetResults;
        m_targetResults = 0x0;
      }
      if (m_decoyResults)
      {
        delete m_decoyResults;
        m_decoyResults = 0x0;
      }
    }
    if (ownOutput)
    {
      if (m_fdrPeptides)
      {
        delete m_fdrPeptides;
        m_fdrPeptides = 0x0;
      }
    }
  }

  // -------------------------------------------------------------------------
  ExecBase * ExecFdrPeptide::clone(const ParameterList & inputParams) const
  {
    return new ExecFdrPeptide(inputParams);
  }
  // -------------------------------------------------------------------------
  bool ExecFdrPeptide::invoke(void)
  {
    DEBUG_TRACE;
    if (m_fdrPeptides == 0x0)
    {
      ownOutput = true;
      m_fdrPeptides = new PeptideSpectrumMatchSet;
    }

    double scalingFactor = 1;

    if (m_params.exists("TRUE_FALSE_DATABASE_RATIO"))
    {
      scalingFactor = m_params.getValueDouble("TRUE_FALSE_DATABASE_RATIO");
    }

    if (m_peptideResults != 0x0)
    {
      if (m_params.exists("TDA_TYPE"))
      {
        if (m_params.getValue("TDA_TYPE").compare("concatenated") == 0)
        {
          if (!FdrPeptide::concatenatedTargetDecoy(*m_peptideResults,
                                                   *m_fdrPeptides,
                                                   scalingFactor))
          {
            return false;
          }
        }
        else if (m_params.getValue("TDA_TYPE").compare("separate") == 0)
        {
          if (!FdrPeptide::separateTargetDecoy(*m_peptideResults,
                                               *m_fdrPeptides,
                                               scalingFactor))
          {
            return false;
          }
        }
        else
        {
          ERROR_MSG("Unknown TDA_TYPE! Valid options are 'concatenated' and 'separate'.");
          return false;
        }
      }
      else
      {
        ERROR_MSG("Unknown TDA_TYPE! Valid options are 'concatenated' and 'separate'.");
      }
    }
    else
    {
      PeptideSpectrumMatchSet mergedPeptides;
      FdrPeptide::concatenateTargetDecoy(*m_targetResults,
                             *m_decoyResults,
                             mergedPeptides);

      if (m_params.exists("TDA_TYPE"))
      {
        if (m_params.getValue("TDA_TYPE").compare("concatenated") == 0)
        {
          if (!FdrPeptide::concatenatedTargetDecoy(mergedPeptides,
                                                   *m_fdrPeptides))
          {
            return false;
          }
        }
        else if (m_params.getValue("TDA_TYPE").compare("separate") == 0)
        {
          if (!FdrPeptide::separateTargetDecoy(mergedPeptides, *m_fdrPeptides))
          {
            return false;
          }
        }
        else
        {
          ERROR_MSG("Unknown TDA_TYPE! Valid options are 'concatenated' and 'separate'.");
          return false;
        }
      }
      else
      {
        ERROR_MSG("Unknown TDA_TYPE! Valid options are 'concatenated' and 'separate'.");
      }
    }

    if (m_params.exists("PEPTIDE_FDR_CUTOFF"))
    {
      double cutoff = m_params.getValueDouble("PEPTIDE_FDR_CUTOFF");
      if (!FdrPeptide::filterByPValue(*m_fdrPeptides, cutoff))
      {
        return false;
      }
    }
    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecFdrPeptide::loadInputData(void)
  {
    //Load in statistics if they haven't been passed in by
    //another function
    if (m_peptideResults == 0x0)
    {
      ownInput = true;
      m_peptideResults = new PeptideSpectrumMatchSet;
    }

    //Load in inspect / specnets results
    if (m_params.exists("INPUT_RESULTS"))
    {
      DEBUG_MSG("Input results: " << m_params.getValue("INPUT_RESULTS"));

      if (m_params.getValue("INPUT_RESULTS_TYPE").compare("inspect") == 0)
      {
        if (!m_peptideResults->loadInspectResultsFile(m_params.getValue("INPUT_RESULTS").c_str(),
                                                      m_params.getValueBool("SCAN_ZERO_INDEX")))
        {
          DEBUG_MSG("Could not load" << m_params.getValue("INPUT_RESULTS"));
          return false;
        }
      }
      else
      {
        if (!m_peptideResults->loadSpecnetsResultsFile(m_params.getValue("INPUT_RESULTS").c_str(),
                                                       m_params.getValueBool("SCAN_ZERO_INDEX")))
        {
          DEBUG_MSG("Could not load" << m_params.getValue("INPUT_RESULTS"));
          return false;
        }
      }
    }

    DEBUG_VAR(m_peptideResults->size());

    return true;
    DEBUG_TRACE;
  }

  // -------------------------------------------------------------------------
  bool ExecFdrPeptide::saveOutputData(void)
  {
    if (m_params.exists("OUTPUT_FDR_RESULTS"))
    {
      return m_fdrPeptides->saveToFile(m_params.getValue("OUTPUT_FDR_RESULTS").c_str(), true);
    }
    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecFdrPeptide::saveInputData(std::vector<std::string> & filenames)
  {

  }

  // -------------------------------------------------------------------------
  bool ExecFdrPeptide::loadOutputData(void)
  {
  }

  // -------------------------------------------------------------------------
  vector<ExecBase*> const & ExecFdrPeptide::split(int numSplit)
  {
    //Mr. Stubbs
  }

  // -------------------------------------------------------------------------
  bool ExecFdrPeptide::merge(void)
  {
    //Mr. Stubbs
  }

  // -------------------------------------------------------------------------
  bool ExecFdrPeptide::validateParams(std::string & error)
  {
    m_isValid = false;

    m_isValid = true;
    return true;
  }

}
