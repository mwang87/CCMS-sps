#include "ExecSpecProtAlignTgtDecoy.h"

#include "alignment_modmut.h"
#include "AlignmentPenaltyBased.h"
#include "FdrPeptide.h"
#include "FileUtils.h"
#include "Logger.h"
#include "ParallelPbsExecution.h"
#include "ParallelSgeExecution.h"
#include "ParallelThreadedExecution.h"

// SpecNets Includes
#include "tags.h"
#include <limits.h>


#define DEBUG_SPECPROTALIGN 1

using namespace std;
using namespace specnets;

namespace specnets
{
  bool comparePValue(psmPtr i, psmPtr j)
  {
    return i->m_pValue > j->m_pValue;
  }

  // -------------------------------------------------------------------------
  ExecSpecProtAlignTgtDecoy::ExecSpecProtAlignTgtDecoy(void) :
        m_inputSpectra(0x0), m_prmSpectra(0x0), m_db(0x0), m_dbDecoy(0x0), 
        m_psmSetTag(0x0), m_psmSetTagDecoy(0x0), ownInput(true),
        m_matchedSpectraAll(0x0), 
        m_psmSet(0x0), m_psmSetDecoy(0x0), m_psmSetFdr(0x0), ownOutput(true)
  {
    m_name = "ExecSpecProtAlignTgtDecoy";
    m_type = "ExecSpecProtAlignTgtDecoy";
  }

  // -------------------------------------------------------------------------

  ExecSpecProtAlignTgtDecoy::ExecSpecProtAlignTgtDecoy(const ParameterList & inputParams) :
    ExecBase(inputParams), m_inputSpectra(0x0), m_prmSpectra(0x0), 
        m_db(0x0), m_dbDecoy(0x0), m_psmSetTag(0x0), m_psmSetTagDecoy(0x0), ownInput(true),
        m_matchedSpectraAll(0x0), 
        m_psmSet(0x0), m_psmSetDecoy(0x0), m_psmSetFdr(0x0), ownOutput(true)
  {
    m_name = "ExecSpecProtAlignTgtDecoy";
    m_type = "ExecSpecProtAlignTgtDecoy";
  }

  // -------------------------------------------------------------------------

  ExecSpecProtAlignTgtDecoy::ExecSpecProtAlignTgtDecoy(const ParameterList & inputParams,
                                       SpecSet * inputSpectra,
                                       SpecSet * prmSpectra,
                                       DB_fasta * db,
                                       DB_fasta * dbDecoy,
                                       PeptideSpectrumMatchSet * psmSetTag,
                                       PeptideSpectrumMatchSet * psmSetTagDecoy,
                                       PenaltyMatrix * penaltyMatrixBlosum,
                                       PenaltyMatrix * penaltyMatrixMods,                                             
                                       SpecSet * matchedSpectraAll,
                                       SpecSet * matchedSpectra,
                                       PeptideSpectrumMatchSet * psmSet,
                                       PeptideSpectrumMatchSet * psmSetDecoy,
                                       PeptideSpectrumMatchSet * psmSetFdr) :
    ExecBase(inputParams), m_inputSpectra(inputSpectra), m_prmSpectra(prmSpectra),
        m_db(db), m_dbDecoy(dbDecoy), m_psmSetTag(psmSetTag), m_psmSetTagDecoy(psmSetTagDecoy),
        m_penaltyMatrixBlosum(penaltyMatrixBlosum), 
        m_penaltyMatrixMods(penaltyMatrixMods), ownInput(false),
        m_matchedSpectraAll(matchedSpectraAll), m_matchedSpectra(matchedSpectra),
        m_psmSet(psmSet), m_psmSetDecoy(psmSetDecoy), m_psmSetFdr(psmSetFdr), ownOutput(false)
  {
    m_name = "ExecSpecProtAlignTgtDecoy";
    m_type = "ExecSpecProtAlignTgtDecoy";
  }
                			

  // -------------------------------------------------------------------------

  ExecSpecProtAlignTgtDecoy::~ExecSpecProtAlignTgtDecoy(void)
  {
    if (ownInput) {
      if (m_inputSpectra)
        delete m_inputSpectra;
      if (m_prmSpectra)
        delete m_prmSpectra;
      if (m_db)
        delete m_db;
      if (m_dbDecoy)
        delete m_dbDecoy;
      if (m_psmSetTag)
        delete m_psmSetTag;
      if (m_psmSetTagDecoy)
        delete m_psmSetTagDecoy;
      if (m_penaltyMatrixBlosum)
        delete m_penaltyMatrixBlosum;
      if (m_penaltyMatrixMods)
        delete m_penaltyMatrixMods;
    }
    if (ownOutput) {
      if (m_matchedSpectraAll)
        delete m_matchedSpectraAll;
      if (m_matchedSpectra)
        delete m_matchedSpectra;
      if (m_psmSet)
        delete m_psmSet;
      if (m_psmSetDecoy)
        delete m_psmSetDecoy;
      if (m_psmSetFdr)
        delete m_psmSetFdr;
    }
  }

  // -------------------------------------------------------------------------

  ExecBase * ExecSpecProtAlignTgtDecoy::clone(const ParameterList & inputParams) const
  {
    return new ExecSpecProtAlignTgtDecoy(inputParams);
  }

  // -------------------------------------------------------------------------

  bool ExecSpecProtAlignTgtDecoy::invoke(void)
  {
    if (!m_inputSpectra || m_inputSpectra->size() == 0) {
      ERROR_MSG("ERROR: empty set of input spectra");
      return false;
    }

    if (!m_db or m_db->size() == 0) {
      ERROR_MSG("ERROR: empty database");
      return false;
    }

    if (!m_matchedSpectraAll or !m_matchedSpectra) {
      ERROR_MSG("ERROR: empty returned data pointer");
      return false;
    }

    float peakTol = (float) m_params.getValueDouble("TOLERANCE_PEAK");
    if (DEBUG_SPECPROTALIGN) DEBUG_VAR(peakTol);
    bool gridExecutionFlag = (bool)m_params.getValueInt("GRID_EXECUTION_FLAG", 0);
    if (DEBUG_SPECPROTALIGN) DEBUG_VAR(gridExecutionFlag);
    bool resume = (bool)m_params.getValueInt("GRID_RESUME_FLAG", 0);
    if (DEBUG_SPECPROTALIGN) DEBUG_VAR(resume);

    m_params.setValue("GRID_DATA_DIR", m_params.getValue("GRID_DATA_DIR_TARGET"));
    m_params.setValue("OUTPUT_MATCHED_PEAKS_IDX_ALL", m_params.getValue("OUTPUT_MATCHED_PEAKS_IDX_TGT"));
    m_params.setValue("OUTPUT_MATCHED_SPECS_ALL", m_params.getValue("OUTPUT_MATCHED_SPECS_TGT"));
    m_params.setValue("OUTPUT_PSM_ALL", m_params.getValue("OUTPUT_PSM_TGT"));

    bool returnStatus;
    // -- TARGET PROTEIN ALIGNMENT
    ExecSpecProtAlign moduleSpecProtAlign(m_params,
                                          m_inputSpectra,
                                          m_prmSpectra,
                                          m_db,
                                          m_penaltyMatrixBlosum,
                                          m_penaltyMatrixMods,
                                          m_matchedSpectraAll);

    if (!m_params.exists("GRID_NUMNODES") or m_params.getValueInt("GRID_NUMNODES") <= 0) {
      returnStatus = moduleSpecProtAlign.invoke();
    } else {
      if (DEBUG_SPECPROTALIGN) DEBUG_TRACE;
      int numNodes = m_params.getValueInt("GRID_NUMNODES");

      string gridType = m_params.getValue("GRID_TYPE");
      if (gridType == "pbs") {
        ParallelPbsExecution exec(&moduleSpecProtAlign,
                                    gridExecutionFlag,
                                    !gridExecutionFlag,
                                  resume);
        returnStatus = exec.invoke(numNodes);
      } else if (gridType == "sge") {
        ParallelSgeExecution exec(&moduleSpecProtAlign,
                                    gridExecutionFlag,
                                    !gridExecutionFlag,
                                  resume);
        returnStatus = exec.invoke(numNodes);
      }
      else if (gridType == "threaded")
      {
	ParallelThreadedExecution exec(&moduleSpecProtAlign);
	returnStatus = exec.invoke(numNodes);
      }
    }

    returnStatus = moduleSpecProtAlign.saveOutputData();
    if (DEBUG_SPECPROTALIGN) DEBUG_VAR(returnStatus);

    float fdrThreshold = m_params.getValueFloat("FDR_THRESHOLD", -1.0);
    if (DEBUG_SPECPROTALIGN) DEBUG_VAR(fdrThreshold);

    // -- DECOY PROTEIN ALIGNMENT
    if (fdrThreshold >= 0.0) {

      m_params.setValue("GRID_DATA_DIR", m_params.getValue("GRID_DATA_DIR_DECOY"));
      m_params.setValue("OUTPUT_MATCHED_PEAKS_IDX_ALL", m_params.getValue("OUTPUT_MATCHED_PEAKS_IDX_DEC"));
      m_params.setValue("OUTPUT_MATCHED_SPECS_ALL", m_params.getValue("OUTPUT_MATCHED_SPECS_DEC"));
      m_params.setValue("OUTPUT_PSM_ALL", m_params.getValue("OUTPUT_PSM_DEC"));

      m_inputSpectra->clearPsms(); // Clear all the target PSMs
      m_psmSetTagDecoy->addSpectra(m_inputSpectra); // Associate decoy PSMs

      ExecSpecProtAlign moduleSpecProtAlign2(m_params,
                                            m_inputSpectra,
                                            m_prmSpectra,
                                            m_dbDecoy,
                                            m_penaltyMatrixBlosum,
                                            m_penaltyMatrixMods,
                                            &m_matchedSpectraAllDecoy);

      if (!m_params.exists("GRID_NUMNODES") or m_params.getValueInt("GRID_NUMNODES") <= 0) {
        returnStatus = moduleSpecProtAlign2.invoke();
      } else {
        if (DEBUG_SPECPROTALIGN) DEBUG_TRACE;
        int numNodes = m_params.getValueInt("GRID_NUMNODES");

        string gridType = m_params.getValue("GRID_TYPE");
        if (gridType == "pbs") {
          ParallelPbsExecution exec(&moduleSpecProtAlign2,
                                    gridExecutionFlag,
                                    !gridExecutionFlag,
                                    resume);
          returnStatus = exec.invoke(numNodes);
        } else if (gridType == "sge") {
          ParallelSgeExecution exec(&moduleSpecProtAlign2,
                                    gridExecutionFlag,
                                    !gridExecutionFlag,
                                    resume);
          returnStatus = exec.invoke(numNodes);
        }
	else if (gridType == "threaded")
	{
	  ParallelThreadedExecution exec(&moduleSpecProtAlign2);
	  returnStatus = exec.invoke(numNodes);
	}
      }

      if ((!m_params.exists("GRID_NUMNODES") || m_params.getValueInt("GRID_NUMNODES") > 0)
          && (!resume && !gridExecutionFlag)) {
        return true;
      }
      
      returnStatus = moduleSpecProtAlign2.saveOutputData();
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(returnStatus);

      m_psmSet->getPSMSet(m_matchedSpectraAll);
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(m_psmSet->size());
      m_psmSetDecoy->getPSMSet(&m_matchedSpectraAllDecoy);
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(m_psmSetDecoy->size());

      // Make sure the isDecoy fields are set on the decoy PSMs
      for (int iDecoy = 0; iDecoy < m_psmSetDecoy->size(); iDecoy++) {
        (*m_psmSetDecoy)[iDecoy]->m_isDecoy = true;
      }

      float scalingFactor = (float)m_psmSet->size() / (float)m_psmSetDecoy->size();
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(scalingFactor);

      PeptideSpectrumMatchSet psmSetCombined;
      if (!FdrPeptide::mergeTargetDecoy(*m_psmSet, *m_psmSetDecoy, psmSetCombined)) {
        ERROR_MSG("FDR concatenation failed.");
      }
      if (DEBUG_SPECPROTALIGN) psmSetCombined.saveToFile("debug_psm_combined.txt");
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(psmSetCombined.size());

      if (!FdrPeptide::concatenatedTargetDecoy(psmSetCombined, *m_psmSetFdr, scalingFactor, &comparePValue)) {
        ERROR_MSG("FDR sorting failed.");
      }
      if (DEBUG_SPECPROTALIGN) m_psmSetFdr->saveToFile("debug_psm_sorted.txt");
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(m_psmSetFdr->size());

      double fdrCutoff  = 1.0;
      if (m_params.getValueFloat("FDR_THRESHOLD", -1.0) >= 0.0) {
        fdrCutoff = m_params.getValueFloat("FDR_THRESHOLD", 1.0);  // Purposely 1.0 not -1.0 (should exist anyway)
      }
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(fdrCutoff);

      if (!FdrPeptide::filterByPValue(*m_psmSetFdr, fdrCutoff)) {
        ERROR_MSG("Filter by FDR value failed.");
      }


#if 0  // LARS LARS LARS LARS LARS LARS LARS LARS LARS LARS LARS
      bool penaltyAlign = (bool) m_params.exists("PENALTY_ALIGNMENT") ?
    		((int) m_params.getValueInt("PENALTY_ALIGNMENT") ? 1 : 0) : 0;
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(penaltyAlign);
      if (penaltyAlign) {
        int maxSpectrumGapDSaltons = m_params.exists("MAX_ALIGN_SPECTRUM_GAP_DALTONS") ?
            max(0, m_params.getValueInt("MAX_ALIGN_SPECTRUM_GAP_DALTONS")) : 1000;
        if (DEBUG_SPECPROTALIGN) DEBUG_VAR(maxSpectrumGapDSaltons);
        float penaltyAlpha = m_params.exists("PENALTY_ALIGNMENT_ALPHA") ?
              m_params.getValueFloat("PENALTY_ALIGNMENT_ALPHA") : 1.0;
        if (DEBUG_SPECPROTALIGN) DEBUG_VAR(penaltyAlpha);
        float penaltyBeta = m_params.exists("PENALTY_ALIGNMENT_BETA") ?
              m_params.getValueFloat("PENALTY_ALIGNMENT_BETA") : 1.0;
        if (DEBUG_SPECPROTALIGN) DEBUG_VAR(penaltyBeta);
        float maxModMass = m_params.exists("MAX_MOD_MASS") ? 
              (float) max(0.0, m_params.getValueDouble("MAX_MOD_MASS")) : 372.2;
        float minModMass = m_params.exists("MIN_MOD_MASS") ?
              (float) max(-372.2, m_params.getValueDouble("MIN_MOD_MASS")) : -372.2;
        AlignmentPenaltyBased apb(maxSpectrumGapDSaltons,
                              m_penaltyMatrixMods,
                              m_penaltyMatrixBlosum,
                              penaltyAlpha,
                              penaltyBeta,
                              maxModMass,
                              minModMass);

        for (int i = 0; i < m_psmSetFdr->size(); i++) {
          string stringAnnotation;
          apb.computeAllGapAnnotations((*m_psmSetFdr)[i]->m_annotation, stringAnnotation);
          if (DEBUG_SPECPROTALIGN) DEBUG_VAR(stringAnnotation);
          (*m_psmSetFdr)[i]->m_annotation = stringAnnotation;
        }
      }
#endif

      if (DEBUG_SPECPROTALIGN) m_psmSetFdr->saveToFile("debug_psm_cutoff.txt");
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(m_psmSetFdr->size());
    }

    if ((!m_params.exists("GRID_NUMNODES") || m_params.getValueInt("GRID_NUMNODES") > 0)
        && (!resume && !gridExecutionFlag)) {
      return true;
    }
      

    // Resize to make enough room
    m_matchedSpectraIndices.resize(m_matchedSpectraAll->size());
    m_matchedSpectra->resize(m_matchedSpectraAll->size());

    int keepIdx = 0;
    for (int i = 0; i < m_matchedSpectraAll->size(); i++) {
      //if (DEBUG_SPECPROTALIGN) DEBUG_VAR(i);
      //if (DEBUG_SPECPROTALIGN) DEBUG_VAR((*m_matchedSpectraAll)[i].size());
      //if (DEBUG_SPECPROTALIGN) DEBUG_VAR((*m_matchedSpectraAll)[i].psmList.size());
      if ((*m_matchedSpectraAll)[i].psmList.size() != 0) {
        m_matchedSpectraIndices[keepIdx] = i;
        (*m_matchedSpectra)[keepIdx++] = (*m_matchedSpectraAll)[i];
      }
    }

    if (DEBUG_SPECPROTALIGN) DEBUG_VAR(keepIdx);
    m_matchedSpectra->resize(keepIdx);
    m_matchedSpectraIndices.resize(keepIdx);

    if (DEBUG_SPECPROTALIGN) DEBUG_VAR(m_matchedSpectra->size());
    if (DEBUG_SPECPROTALIGN) DEBUG_VAR(m_matchedSpectraIndices.size());

    return true;
  }

  // -------------------------------------------------------------------------

  bool ExecSpecProtAlignTgtDecoy::loadInputData(void)
  {
    if (ownInput) {
      if (!m_inputSpectra)
        m_inputSpectra = new SpecSet;
      if (!m_prmSpectra)
        m_prmSpectra = new SpecSet;
      if (!m_db)
        m_db = new DB_fasta;
      if (!m_dbDecoy)
        m_dbDecoy = new DB_fasta;
      if (!m_psmSetTag)
        m_psmSetTag = new PeptideSpectrumMatchSet;
      if (!m_psmSetTagDecoy)
        m_psmSetTagDecoy = new PeptideSpectrumMatchSet;
    }
    m_inputSpectra->resize(0);

    if (ownOutput) {
      if (!m_matchedSpectraAll)
        m_matchedSpectraAll = new SpecSet;
      if (!m_psmSet)
        m_psmSet = new PeptideSpectrumMatchSet;
      if (!m_psmSetDecoy)
        m_psmSetDecoy = new PeptideSpectrumMatchSet;
      if (!m_psmSetFdr)
        m_psmSetFdr = new PeptideSpectrumMatchSet;
    }
    m_matchedSpectraAll->resize(0);

    if (m_params.exists("INPUT_SPECS_PKLBIN") and m_params.exists("INPUT_PSM")) {
      if (m_inputSpectra->loadPklBin(m_params.getValue("INPUT_SPECS_PKLBIN").c_str(),
                                     m_params.getValue("INPUT_PSM").c_str(),
                                     m_params.getValue("INPUT_MATCHED_PEAKS_IDX").c_str()) < 0) {
        ERROR_MSG("Error reading input spectra files.");
        return false;
      }
    }

    if (m_params.exists("INPUT_PRM_PKLBIN"))
    {
      if (m_prmSpectra->loadPklBin(m_params.getValue("INPUT_PRM_PKLBIN").c_str())
          < 0)
      {
        ERROR_MSG("Error reading input spectra files.");
        return false;
      }
    }

    if (!m_params.exists("INPUT_FASTA")) {
      ERROR_MSG("Parameters are incomplete. INPUT_FASTA is missing.");
      return false;
    }
    else if (m_db->Load(m_params.getValue("INPUT_FASTA").c_str()) <= 0) {
      ERROR_MSG("Error reading database sequences from "
          << m_params.getValue("INPUT_FASTA"));
      return false;
    }

    if (!m_params.exists("INPUT_FASTA_DECOY")) {
      if (m_dbDecoy->Load(m_params.getValue("INPUT_FASTA_DECOY").c_str()) <= 0) {
        ERROR_MSG("Error reading decoy database sequences from "
            << m_params.getValue("INPUT_FASTA_DECOY"));
        return false;
      }
    }

    if (m_params.exists("INPUT_PSM")) {
      m_psmSetTag->loadFromFile(m_params.getValue("INPUT_PSM").c_str());
    }

    if (m_params.exists("INPUT_PSM_DECOY")) {
      m_psmSetTagDecoy->loadFromFile(m_params.getValue("INPUT_PSM_DECOY").c_str());
    }

    //---------------------------------------------------------------------------
    // Load penalty matrices for new alignment
    //---------------------------------------------------------------------------

    bool penaltyAlign = (bool) m_params.exists("PENALTY_ALIGNMENT") ?
    		((int) m_params.getValueInt("PENALTY_ALIGNMENT") ? 1 : 0) : 0;

    if (penaltyAlign) {

      // Load amino acid masses
      AAJumps jumps(1);
      if (!m_params.exists("AMINO_ACID_MASSES")) {
        jumps.loadJumps(m_params.getValue("AMINO_ACID_MASSES").c_str(), true);
      }

      float resolution = m_params.getValueDouble("ALIGNMENT_RESOLUTION", 0.1);

      if (!m_penaltyMatrixBlosum)
        m_penaltyMatrixBlosum = new PenaltyMatrix(jumps, resolution);
      if (!m_penaltyMatrixMods)
        m_penaltyMatrixMods = new PenaltyMatrix(jumps, resolution);

      if (!m_params.exists("BLOSUM_PENALTY_FILE")) {
        ERROR_MSG("Parameters are incomplete. BLOSUM_FILE is missing and penalty alignment was specified.");
        return false;
      } else {
        string blosumFileName = m_params.getValue("BLOSUM_PENALTY_FILE"); 
        if (!m_penaltyMatrixBlosum->load(blosumFileName) ) {
          ERROR_MSG("Error loading blosum penalties from " << blosumFileName);
          return false;
        }
      }

      if (!m_params.exists("KNOWN_MODS_FILE")) {
        ERROR_MSG("Parameters are incomplete. KNOWN_MODS_FILE is missing and penalty alignment was specified.");
        return false;
      } else {
        string knowmModsFileName = m_params.getValue("KNOWN_MODS_FILE"); 
        m_penaltyMatrixMods->loadKnownModifications(knowmModsFileName);
      }

      if (!m_params.exists("MODS_PENALTY_FILE")) {
        ERROR_MSG("Parameters are incomplete. MODS_PENALTY_FILE is missing and penalty alignment was specified.");
        return false;
      } else {
        string modFileName = m_params.getValue("MODS_PENALTY_FILE"); 
        if (!m_penaltyMatrixMods->load(modFileName) ) {
          ERROR_MSG("Error loading mod penalties from " << modFileName);
          return false;
        }
      }

    } // if (penaltyAlign)

    return true;
  }

  // -------------------------------------------------------------------------

  bool ExecSpecProtAlignTgtDecoy::saveOutputData(void)
  {
    bool penaltyAlign = (bool) m_params.exists("PENALTY_ALIGNMENT") ?
  	  	((int) m_params.getValueInt("PENALTY_ALIGNMENT") ? 1 : 0) : 0;
    DEBUG_VAR(penaltyAlign);

    if (m_params.exists("OUTPUT_MATCHED_PEAKS_IDX")) {
      int keepIdx = 0;
      for (int i = 0; i < m_matchedSpectraAll->size(); i++) {
        if ((*m_matchedSpectraAll)[i].psmList.size() != 0) {
          (*m_matchedSpectra)[keepIdx++] = (*m_matchedSpectraAll)[i];
        }
      }

      DEBUG_MSG("Saving matched peaks...");
      SpecSet tempMatchedPeaks;
      tempMatchedPeaks.resize(m_matchedSpectra->size());
      for (int i = 0; i < m_matchedSpectra->size(); i++) {
        if ((*m_matchedSpectra)[i].psmList.size() != 0) {
          list<psmPtr>::iterator litr = (*m_matchedSpectra)[i].psmList.begin();
          int peakListSize = (*litr)->m_matchedPeaks.size();
          //DEBUG_VAR(peakListSize);
          tempMatchedPeaks[i].resize(peakListSize);
          for (int j = 0; j < peakListSize; j++) {
            tempMatchedPeaks[i][j].set((*litr)->m_matchedPeaks[j][0],
                                       (*litr)->m_matchedPeaks[j][1]);
            //DEBUG_MSG((*litr)->m_matchedPeaks[j][0] << "  " << (*litr)->m_matchedPeaks[j][1]);
          }
        }
      }
      tempMatchedPeaks.savePklBin(m_params.getValue("OUTPUT_MATCHED_PEAKS_IDX").c_str());
    }

    if (DEBUG_SPECPROTALIGN) DEBUG_VAR(m_matchedSpectra->size());
    if (m_matchedSpectra and m_params.exists("OUTPUT_MATCHED_SPECS")) {
      if (m_params.exists("OUTPUT_PSM") && m_params.exists("OUTPUT_MATCHED_PEAKS_IDX_ALL")) {
        DEBUG_MSG("Saving matched specs (3 files)...");
        m_matchedSpectra->savePklBin(m_params.getValue("OUTPUT_MATCHED_SPECS").c_str(),
                                     m_params.getValue("OUTPUT_PSM").c_str(),
                                     m_params.getValue("OUTPUT_MATCHED_PEAKS_IDX_ALL").c_str());

      } else if (m_params.exists("OUTPUT_PSM")) {
        DEBUG_MSG("Saving matched specs (2 files)...");
        m_matchedSpectra->savePklBin(m_params.getValue("OUTPUT_MATCHED_SPECS").c_str(),
                                     m_params.getValue("OUTPUT_PSM").c_str());

      } else {
        DEBUG_MSG("Saving matched specs (1 file)...");
        m_matchedSpectra->savePklBin(m_params.getValue("OUTPUT_MATCHED_SPECS").c_str());
      }
    }

    //=================================================
    // LARS - TEMPORARY MEASURE SO REPORTS KEEP WORKING
    //=================================================
    if (DEBUG_SPECPROTALIGN) DEBUG_VAR(m_matchedSpectraAll->size());
    if (m_params.exists("OUTPUT_MATCHED_PROTS_ALL")) {
      DEBUG_MSG("Saving matched prots all...");
      m_matchedSpectraAll->saveMatchedProts(m_params.getValue("OUTPUT_MATCHED_PROTS_ALL").c_str());
    }
    if (m_params.exists("OUTPUT_MATCHED_PROTS")) {
      DEBUG_MSG("Saving matched prots...");
      m_matchedSpectra->saveMatchedProts(m_params.getValue("OUTPUT_MATCHED_PROTS").c_str());
    }
    if (m_params.exists("OUTPUT_MATCHED_SPECS_IDX")) {
      DEBUG_MSG("Saving matched specs idx...");
      Save_binArray(m_params.getValue("OUTPUT_MATCHED_SPECS_IDX").c_str(),
                    m_matchedSpectraIndices);
    }
    //=================================================

    if (m_params.exists("OUTPUT_PSM_MOD_FREQS")) {
      m_psmSetFdr->saveModMatrix(m_params.getValue("OUTPUT_PSM_MOD_FREQS").c_str(), false);
    }

    if (m_params.exists("OUTPUT_PSM_FDR")) {
      m_psmSetFdr->saveToFile(m_params.getValue("OUTPUT_PSM_FDR").c_str());
    }

    return true;
  }

  // -------------------------------------------------------------------------

  bool ExecSpecProtAlignTgtDecoy::saveInputData(std::vector<std::string> & filenames)
  {
    return false;
  }

  // -------------------------------------------------------------------------

  bool ExecSpecProtAlignTgtDecoy::loadOutputData(void)
  {
    return false;
  }

  // -------------------------------------------------------------------------

  vector<ExecBase *> const & ExecSpecProtAlignTgtDecoy::split(int numSplit)
  {
    m_subModules.resize(0);
    return m_subModules;
  }

  // -------------------------------------------------------------------------

  bool ExecSpecProtAlignTgtDecoy::merge(void)
  {
    return false;
  }


  // -------------------------------------------------------------------------

  bool ExecSpecProtAlignTgtDecoy::validateParams(std::string & error)
  {
    m_isValid = false;

    VALIDATE_PARAM_EXIST("TOLERANCE_PEAK");
    VALIDATE_PARAM_EXIST("TOLERANCE_PM");

    m_isValid = true;
    return true;
  }

// -------------------------------------------------------------------------


}

