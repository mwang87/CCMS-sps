#include "ExecSpecProtAlign.h"

#include "alignment_modmut.h"
#include "AlignmentPenaltyBased.h"
#include "FdrPeptide.h"

// Module Includes
#include "Logger.h"
#include "FileUtils.h"

// SpecNets Includes
#include "tags.h"
#include <limits.h>
#include <time.h>


static bool DEBUG_SPECPROTALIGN = false;
static bool DEBUG_SPECPROTALIGN_SPECTRA = false;
static bool DEBUG_SPECPROTALIGN_TIME = false;
static bool DEBUG_SPECPROTALIGN_SPRINKLE = false;
static bool DEBUG_SPECPROTALIGN_ANNO = false;
static bool DEBUG_SPECPROTALIGN_MERGE = false;

using namespace std;
using namespace specnets;

namespace specnets
{
  class LessThanPredicate
  {
  public:
    LessThanPredicate(float value) :
      m_value(value)
    {
    }
    bool operator()(const psmPtr & value)
    {
      return value->m_score < m_value;
    }
  private:
    float m_value;
  };

  bool PsmDbIndexAnnoSort(psmPtr left, psmPtr right)
  {
    if (left->m_dbIndex == right->m_dbIndex) {
      return left->m_annotation < right->m_annotation;
    } else {
      return left->m_dbIndex < right->m_dbIndex;
    }
  }

  bool PsmDbIndexAnnoUnique(psmPtr left, psmPtr right)
  {
    return (left->m_dbIndex == right->m_dbIndex) && (left->m_annotation == right->m_annotation);
  }

  void DebugPsm(psmPtr & psm)
  {
    DEBUG_MSG(psm->m_scanNum << "  " << 
              psm->m_matchOrientation << "  " << 
              psm->m_annotation << "  " << 
              psm->m_score << "  " << 
              psm->m_dbIndex << "  " << 
              psm->m_protein);
    return;
  } 

  void ClearEndpoints(Spectrum & cSpec)
  {
    if (DEBUG_SPECPROTALIGN_SPECTRA) DEBUG_VAR(cSpec.size());
    if (DEBUG_SPECPROTALIGN_SPECTRA) DEBUG_MSG("Original spectrum");
    if (DEBUG_SPECPROTALIGN_SPECTRA) cSpec.output(cerr);

    // Clean up any pre-existing endpoints so we don't add too many
    list<int> peaksToRemove;
    for (int iPeak = 0; iPeak < cSpec.size(); iPeak++) {
      if (cSpec[iPeak][0] < 57.0 || cSpec[iPeak][0] > cSpec.parentMass - 57.0) {
        peaksToRemove.push_back(iPeak);
      }
    }
    if (peaksToRemove.size() != 0) {
      cSpec.removePeaks(peaksToRemove);
    }
    if (DEBUG_SPECPROTALIGN_SPECTRA) DEBUG_VAR(cSpec.size());
    return;
  }

  // -------------------------------------------------------------------------
  ExecSpecProtAlign::ExecSpecProtAlign(void) :
    m_inputSpectra(0x0), m_prmSpectra(0x0), m_db(0x0), 
    m_penaltyMatrixBlosum(0x0), m_penaltyMatrixMods(0x0),
    ownInput(true), m_matchedSpectraAll(0x0),
        ownOutput(true)
  {
    m_name = "ExecSpecProtAlign";
    m_type = "ExecSpecProtAlign";
  }

  // -------------------------------------------------------------------------
  ExecSpecProtAlign::ExecSpecProtAlign(const ParameterList & inputParams,
                                       SpecSet * inputSpectra,
                                       SpecSet * prmSpectra,
                                       DB_fasta * db,
                                       PenaltyMatrix * penaltyMatrixBlosum,
                                       PenaltyMatrix * penaltyMatrixMods) :
    ExecBase(inputParams), m_inputSpectra(inputSpectra), m_prmSpectra(prmSpectra), m_db(db),
        m_penaltyMatrixBlosum(penaltyMatrixBlosum),
        m_penaltyMatrixMods(penaltyMatrixMods), ownInput(true),
        m_matchedSpectraAll(0x0), ownOutput(true)
  {
    m_name = "ExecSpecProtAlign";
    m_type = "ExecSpecProtAlign";
  }

  // -------------------------------------------------------------------------

  ExecSpecProtAlign::ExecSpecProtAlign(const ParameterList & inputParams) :
    ExecBase(inputParams), m_inputSpectra(0x0), m_prmSpectra(0x0), m_db(0x0), 
       m_penaltyMatrixBlosum(0x0), m_penaltyMatrixMods(0x0), ownInput(true),
       m_matchedSpectraAll(0x0), ownOutput(true)
  {
    m_name = "ExecSpecProtAlign";
    m_type = "ExecSpecProtAlign";
  }

  // -------------------------------------------------------------------------

  ExecSpecProtAlign::ExecSpecProtAlign(const ParameterList & inputParams,
                                       SpecSet * inputSpectra,
                                       SpecSet * prmSpectra,
                                       DB_fasta * db,
                                       PenaltyMatrix * penaltyMatrixBlosum,
                                       PenaltyMatrix * penaltyMatrixMods,
                                       SpecSet * matchedSpectraAll) :
    ExecBase(inputParams), m_inputSpectra(inputSpectra), m_prmSpectra(prmSpectra), m_db(db),
        m_penaltyMatrixBlosum(penaltyMatrixBlosum),
        m_penaltyMatrixMods(penaltyMatrixMods), ownInput(false),
        m_matchedSpectraAll(matchedSpectraAll), ownOutput(false)
  {
    m_name = "ExecSpecProtAlign";
    m_type = "ExecSpecProtAlign";
  }

  // -------------------------------------------------------------------------

  ExecSpecProtAlign::~ExecSpecProtAlign(void)
  {
#if 0
    if (ownInput)
    {
      if (m_inputSpectra) {
        delete m_inputSpectra;
      }
    }
    if (ownInput)
    {
      if (m_db) {
        delete m_db;
      }
    }
    if (ownOutput)
    {
      if (m_matchedSpectraAll) {
        delete m_matchedSpectraAll;
      }
    }
#endif
  }

  // -------------------------------------------------------------------------

  ExecBase * ExecSpecProtAlign::clone(const ParameterList & inputParams) const
  {
    return new ExecSpecProtAlign(inputParams);
  }

  // -------------------------------------------------------------------------

  bool ExecSpecProtAlign::invoke(void)
  {
    if (!m_inputSpectra || m_inputSpectra->size() == 0) {
      ERROR_MSG("ERROR: empty set of input spectra");
      return false;
    }

    // More readable code if we use a reference instead of pointer
    SpecSet & inputSpectra = *m_inputSpectra;

#if 0
    bool usePrmSprinkling = m_params.exists("SPECPROTALIGN_PRM_CONTRIBUTION_MULTIPLIER");
    DEBUG_VAR(m_params.exists("SPECPROTALIGN_PRM_CONTRIBUTION_MULTIPLIER"));
    if (usePrmSprinkling && m_prmSpectra != 0x0 && m_prmSpectra->size() != 0) {
      for (int i = 0; i < inputSpectra[i].size(); i++) {
        inputSpectra[i].resize(0);
        inputSpectra[i].mergePeakList((*m_prmSpectra)[i], 0x0);
      }
    }

    DEBUG_VAR(inputSpectra.size());
    if (m_prmSpectra != 0x0 && m_prmSpectra->size() != 0) {
      DEBUG_MSG("Good PRM");
      inputSpectra = *m_prmSpectra;
      DEBUG_VAR(inputSpectra.size());
      for (int i = 0; i < inputSpectra[i].size(); i++) {
        DEBUG_VAR((*m_inputSpectra)[i].psmList.size());
        // Copy PSMs to the PRM spectrum
        list<psmPtr>::iterator itr = (*m_inputSpectra)[i].psmList.begin();
        list<psmPtr>::iterator itr_end = (*m_inputSpectra)[i].psmList.end();
        for (; itr != itr_end; itr++) {
          inputSpectra[i].psmList.push_back(*itr);
        }
        DEBUG_VAR(inputSpectra[i].psmList.size());
      }
    }
#endif

    if (!m_db or m_db->size() == 0) {
      ERROR_MSG("ERROR: empty database");
      return false;
    }

    if (!m_matchedSpectraAll) {
      m_matchedSpectraAll = new SpecSet;
      ownOutput = true;
    }

    // -----------------------------------------------------
    // Get all the parameters for alignment
    // -----------------------------------------------------
    float pmTol = (float)m_params.getValueDouble("TOLERANCE_PM");
    DEBUG_VAR(pmTol);
    float peakTol = (float)m_params.getValueDouble("TOLERANCE_PEAK");
    DEBUG_VAR(peakTol);

    unsigned int startIdx = m_params.exists("IDX_START")
        ? max(0, m_params.getValueInt("IDX_START")) : 0;
    DEBUG_VAR(startIdx);
    unsigned int endIdx = m_params.exists("IDX_END") ? min(inputSpectra.size()
        - 1, (unsigned int)m_params.getValueInt("IDX_END"))
        : inputSpectra.size() - 1;
    DEBUG_VAR(endIdx);

    unsigned int scanFirst = m_params.getValueInt("SCAN_FIRST", 0);
    DEBUG_VAR(scanFirst);
    unsigned int scanLast = m_params.getValueInt("SCAN_LAST", 
					inputSpectra[inputSpectra.size()-1].scan);
    DEBUG_VAR(scanLast);

    unsigned int startIdx_db = m_params.exists("IDX_START_DB")
        ? max(0, m_params.getValueInt("IDX_START_DB")) : 0;
    unsigned int endIdx_db = m_params.exists("IDX_END_DB") ? min(m_db->size()
        - 1, (unsigned int)m_params.getValueInt("IDX_END_DB")) : m_db->size()
        - 1;

    DEBUG_VAR(m_db->size());
    if (startIdx_db > m_db->size())
      return false;

    bool enforceEndpeaks = m_params.exists("ENFORCE_ENDPEAKS")
        ? m_params.getValueInt("ENFORCE_ENDPEAKS") > 0 : false;
    DEBUG_VAR(enforceEndpeaks);

    float thresholdPercent = m_params.exists("ALIGNMENT_SCORE_THRESHOLD")
        ? m_params.getValueDouble("ALIGNMENT_SCORE_THRESHOLD") : 0.90;
    DEBUG_VAR(thresholdPercent);

    float maxModMass = m_params.exists("MAX_MOD_MASS")
        ? (float)max(0.0, m_params.getValueDouble("MAX_MOD_MASS")) : 372.2;
    float minModMass = m_params.exists("MIN_MOD_MASS")
        ? (float)max(-372.2, m_params.getValueDouble("MIN_MOD_MASS")) : -372.2;
    int showMatchedPRMs = m_params.exists("SHOW_MATCHED_PRMS")
        ? m_params.getValueInt("SHOW_MATCHED_PRMS") : 0;
    int maxDbGapAas = m_params.exists("MAX_ALIGN_DB_GAP_AAS")
        ? max(0, m_params.getValueInt("MAX_ALIGN_DB_GAP_AAS")) : 6;
    DEBUG_VAR(maxDbGapAas);
    int maxSpectrumGapDaltons =
        m_params.exists("MAX_ALIGN_SPECTRUM_GAP_DALTONS")
            ? max(0, m_params.getValueInt("MAX_ALIGN_SPECTRUM_GAP_DALTONS"))
            : 1000;
    DEBUG_VAR(maxSpectrumGapDaltons);

    int maxNumMods = m_params.exists("MAX_NUM_MODS")
        ? max(0, m_params.getValueInt("MAX_NUM_MODS")) : 1;
    int minNumMatchPeaks = m_params.exists("MIN_MATCHED_PEAKS_DB")
        ? max(0, m_params.getValueInt("MIN_MATCHED_PEAKS_DB")) : 7;
    int matchesPerMod = m_params.exists("MATCHES_PER_MOD")
        ? max(0, m_params.getValueInt("MATCHES_PER_MOD")) : 2;
    bool tagParsimony = (bool)m_params.exists("MAX_PARSIMONY")
        ? ((int)m_params.getValueInt("MAX_PARSIMONY") ? 1 : 0) : 0;
    DEBUG_VAR(tagParsimony);

    bool penaltyAlign = (bool)m_params.exists("PENALTY_ALIGNMENT")
        ? ((int)m_params.getValueInt("PENALTY_ALIGNMENT") ? 1 : 0) : 0;
    DEBUG_VAR(penaltyAlign);
    float penaltyAlpha = m_params.exists("PENALTY_ALIGNMENT_ALPHA")
        ? m_params.getValueFloat("PENALTY_ALIGNMENT_ALPHA") : 1.0;
    DEBUG_VAR(penaltyAlpha);
    float penaltyBeta = m_params.exists("PENALTY_ALIGNMENT_BETA")
        ? m_params.getValueFloat("PENALTY_ALIGNMENT_BETA") : 1.0;
    DEBUG_VAR(penaltyBeta);
    bool alignall = (bool)m_params.exists("ALIGN_ALL")
        ? ((int)m_params.getValueInt("ALIGN_ALL") ? 1 : 0) : 0;
    DEBUG_VAR(alignall);

    bool resetEndpoints = (bool)m_params.exists("RESET_ENDPOINTS")
        ? ((int)m_params.getValueInt("RESET_ENDPOINTS") ? 1 : 0) : 0;
    DEBUG_VAR(resetEndpoints);

    bool usePrmSprinkling = m_params.exists("SPECPROTALIGN_PRM_CONTRIBUTION_MULTIPLIER");
    DEBUG_VAR(m_params.exists("SPECPROTALIGN_PRM_CONTRIBUTION_MULTIPLIER"));
    float prmMultiplier = m_params.exists("SPECPROTALIGN_PRM_CONTRIBUTION_MULTIPLIER")
        ? m_params.getValueFloat("SPECPROTALIGN_PRM_CONTRIBUTION_MULTIPLIER") : 0.2;
    DEBUG_VAR(prmMultiplier);

    bool roundAnno = m_params.exists("SPECPROTALIGN_ROUND_ANNOTATION_MAX");
    DEBUG_VAR(roundAnno);
    float roundAnnoMax = m_params.exists("SPECPROTALIGN_ROUND_ANNOTATION_MAX")
        ? m_params.getValueFloat("SPECPROTALIGN_ROUND_ANNOTATION_MAX") : 1.0;
    DEBUG_VAR(roundAnnoMax);

    bool specTypeMSMS = m_params.getValueBool("SPEC_TYPE_MSMS", false);
    float ionOffset = specTypeMSMS ? AAJumps::massHion : 0;

    DEBUG_SPECPROTALIGN = m_params.getValueBool("DEBUG_SPECPROTALIGN");
    DEBUG_VAR(DEBUG_SPECPROTALIGN);
    DEBUG_SPECPROTALIGN_SPECTRA = m_params.getValueBool("DEBUG_SPECPROTALIGN_SPECTRA");
    DEBUG_VAR(DEBUG_SPECPROTALIGN_SPECTRA);
    DEBUG_SPECPROTALIGN_TIME = m_params.getValueBool("DEBUG_SPECPROTALIGN_TIME");
    DEBUG_VAR(DEBUG_SPECPROTALIGN_TIME);
    DEBUG_SPECPROTALIGN_SPRINKLE = m_params.getValueBool("DEBUG_SPECPROTALIGN_SPRINKLE");
    DEBUG_VAR(DEBUG_SPECPROTALIGN_SPRINKLE);
    DEBUG_SPECPROTALIGN_ANNO = m_params.getValueBool("DEBUG_SPECPROTALIGN_ANNO");
    DEBUG_VAR(DEBUG_SPECPROTALIGN_ANNO);
    DEBUG_VAR(m_params.getValue("DEBUG_SPECPROTALIGN_SINGLESPECTRUM"));
    

    int modIdx, specIdx, protIdx;

    Spectrum tmpSpec;
    tmpSpec.reserve(1024);
    Spectrum cSpec;
    cSpec.reserve(1024);
    Spectrum cSpecRev;
    cSpecRev.reserve(1024);
    AMM_match *bestMatch;

    // -----------------------------------------------------
    // Construct the penalty based alignment class
    //   this will initialize all the scoring vectors
    // -----------------------------------------------------
    AlignmentPenaltyBased * apb = 0x0;
    if (penaltyAlign) {
      apb = new AlignmentPenaltyBased(maxSpectrumGapDaltons,
                              m_penaltyMatrixMods,
                              m_penaltyMatrixBlosum,
                              penaltyAlpha,
                              penaltyBeta,
                              maxModMass,
                              minModMass);
    }
    
    DEBUG_VAR(inputSpectra.size());
    for (specIdx = 0; specIdx < inputSpectra.size(); specIdx++) {
      DEBUG_VAR(specIdx);

      int scanNum = inputSpectra[specIdx].scan;
      DEBUG_VAR(scanNum);
      if (scanNum < scanFirst || scanNum > scanLast) {
        continue;
      }

      if ( m_params.exists("DEBUG_SPECPROTALIGN_SINGLESPECTRUM")) {
        if (m_params.getValueInt("DEBUG_SPECPROTALIGN_SINGLESPECTRUM") != scanNum) {
          continue;
        }
      }

      if (DEBUG_SPECPROTALIGN_SPECTRA) DEBUG_VAR(inputSpectra[specIdx].parentMass);
      if (DEBUG_SPECPROTALIGN_SPECTRA) DEBUG_VAR(inputSpectra[specIdx].size());

      // Apparently it is possible to have a contig that is empty
      if (inputSpectra[specIdx].size() == 0) {
        m_matchedSpectraAll->push_back(inputSpectra[specIdx]);
        continue;
      }

#if 0
      if (usePrmSprinkling) {
        SprinkleSpectrumWithPrms(inputSpectra[specIdx], cSpec, prmMultiplier);
        DEBUG_VAR(cSpec.size());
      } else {
        cSpec = inputSpectra[specIdx];
      }
#else
      cSpec = inputSpectra[specIdx];
#endif
      cSpec.psmList.clear();

      // -----------------------------------------------------
      // "Massage" the spectra for better alignment
      // -----------------------------------------------------
      // Smash down the intensities of the end points so they don't "over-contribute"
      for (int iPeak = 0; iPeak < cSpec.size(); iPeak++) {
        if (cSpec[iPeak][0] < 57.0 || cSpec[iPeak][0] > cSpec.parentMass - 57.0) {
          cSpec[iPeak][1] = 0.1;
        }
      }

      if (resetEndpoints) {
        // Clean up any pre-existing endpoints so we don't add too many
        list<int> peaksToRemove;
        for (int iPeak = 0; iPeak < cSpec.size(); iPeak++) {
          if (cSpec[iPeak][0] < 57.0 || cSpec[iPeak][0] > cSpec.parentMass - 57.0) {
            peaksToRemove.push_back(iPeak);
          }
        }
        if (peaksToRemove.size() != 0) {
          cSpec.removePeaks(peaksToRemove);
        }
      }

      if (DEBUG_SPECPROTALIGN_SPECTRA) cSpec.output(cerr);

      // Get the reverse spectrum
      cSpec.reverse(0, &cSpecRev);
      cSpecRev.psmList.clear();

      if (DEBUG_SPECPROTALIGN_SPECTRA) cSpecRev.output(cerr);

      if (resetEndpoints && enforceEndpeaks) {
        cSpec.addZPMpeaks(peakTol, ionOffset, false);
        cSpecRev.addZPMpeaks(peakTol, ionOffset, false);
      }

      // Round all the spectra to "center them at 0.0"
      cSpec.roundPeaks();
      cSpecRev.roundPeaks();

      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(inputSpectra[specIdx].psmList.size())

      // -----------------------------------------------------
      // Coalesce the starting positions on all proteins
      // -----------------------------------------------------
      map<int, set<float> > matchedProtMapForward;
      map<int, set<float> > matchedProtMapReverse;
      if (alignall) {
        DEBUG_MSG("Inserting dummy start mass for every DB index")
        // If ALIGN_ALL is specified insert dummy start mass for every db idx
        // This will cause alignment to every protein in database
        for (int idx = 0; idx < m_db->size(); idx++) {
          matchedProtMapForward[idx].insert(0.0);
          matchedProtMapReverse[idx].insert(0.0);
        }
      } else if (inputSpectra[specIdx].psmList.size() != 0) {
        list<psmPtr>::iterator itr = inputSpectra[specIdx].psmList.begin();
        list<psmPtr>::iterator itrEnd = inputSpectra[specIdx].psmList.end();
        for (; itr != itrEnd; itr++) {
          int idx = (*itr)->m_dbIndex;
          // Only use the forward orientation tags(if tag parsimony was not done)
          if (tagParsimony || (*itr)->m_matchOrientation == 0) {
            if (matchedProtMapForward.find(idx) == matchedProtMapForward.end()) {
              set<float> newList;
              matchedProtMapForward[idx] = newList;
            }
            matchedProtMapForward[idx].insert((*itr)->m_startMass);
            //if (DEBUG_SPECPROTALIGN) DEBUG_MSG(idx << "  " << (*itr)->m_startMass)
          }
          if (tagParsimony || (*itr)->m_matchOrientation == 1) {
            if (matchedProtMapReverse.find(idx) == matchedProtMapReverse.end()) {
              set<float> newList;
              matchedProtMapReverse[idx] = newList;
            }
            matchedProtMapReverse[idx].insert((*itr)->m_startMass);
            //if (DEBUG_SPECPROTALIGN) DEBUG_MSG(idx << "  " << (*itr)->m_startMass)
          }
        } // for (; itr != itrEnd; itr++)
      }
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(matchedProtMapForward.size());
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(matchedProtMapReverse.size());

      DEBUG_MSG("Matching as b...");
      clock_t startTime = clock();
      // -----------------------------------------------------
      // Forward matching
      // -----------------------------------------------------
      map<int, set<float> >::iterator itrMap = matchedProtMapForward.begin();
      map<int, set<float> >::iterator itrMapEnd = matchedProtMapForward.end();
      for (; itrMap != itrMapEnd; itrMap++) {
        int protIdx = itrMap->first;
        //if (DEBUG_SPECPROTALIGN) DEBUG_VAR(protIdx);
        // If parsimony was done in TagSearch then don't use start positions
        if (tagParsimony) {
          itrMap->second.clear();
        }
        Spectrum dbSpec;
        dbSpec.reserve(1024);
        dbSpec = m_db->getMassesSpec(protIdx);
        // Round all the db spectra to "center them at 0.0"
        dbSpec.roundPeaks();

        if (!penaltyAlign) {
          scoreOverlapAMM(cSpec,
                          dbSpec,
                          protIdx,
                          0,
                          itrMap->second,
                          maxNumMods,
                          minNumMatchPeaks,
                          pmTol,
                          peakTol,
                          57,
                          maxModMass,
                          minModMass,
                          enforceEndpeaks);
        } else {
	    string dbString = m_db->getSequence(protIdx);
            apb->computeAlignment(cSpec,
                               dbSpec,
                               m_db->getSequence(protIdx),
                               protIdx,
                               0,
                               itrMap->second,
                               minNumMatchPeaks,
                               maxDbGapAas,
                               pmTol,
                               peakTol,
                               enforceEndpeaks);
        }
      } // for (; itrMap !=  itrMapEnd; itrMap++) {

      float elapsed = ((float)(clock() - startTime)) / CLOCKS_PER_SEC;
      if (DEBUG_SPECPROTALIGN_TIME) DEBUG_MSG(elapsed << " seconds elapsed");
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(cSpec.psmList.size());
      if (DEBUG_SPECPROTALIGN_TIME) DEBUG_MSG((cSpec.psmList.size() == 0 ? 0 : elapsed / (float)cSpec.psmList.size()) << " per protein");

      list<psmPtr>::iterator litr = cSpec.psmList.begin();
      list<psmPtr>::iterator litrEnd = cSpec.psmList.end();
      for (; litr != litrEnd; litr++) {
        (*litr)->m_scanNum = cSpec.scan;
        (*litr)->m_protein = m_db->getID((*litr)->m_dbIndex);
        m_allPsms.push_back(*litr);
      }

      DEBUG_MSG("Matching as y...");
      startTime = clock();
      // -----------------------------------------------------
      // Reverse matching
      // -----------------------------------------------------
      itrMap = matchedProtMapReverse.begin();
      itrMapEnd = matchedProtMapReverse.end();
      for (; itrMap != itrMapEnd; itrMap++) {
        int protIdx = itrMap->first;
        //if (DEBUG_SPECPROTALIGN) DEBUG_VAR(protIdx);
        // If parsimony was done in TagSearch then don't use start positions
        if (tagParsimony) {
          itrMap->second.clear();
        }
        Spectrum dbSpec;
        dbSpec.reserve(1024);
        dbSpec = m_db->getMassesSpec(protIdx);
        // Round all the db spectra to "center them at 0.0"
        dbSpec.roundPeaks();
        if (!penaltyAlign) {
          scoreOverlapAMM(cSpecRev,
                          dbSpec,
                          protIdx,
                          1,
                          itrMap->second,
                          maxNumMods,
                          minNumMatchPeaks,
                          pmTol,
                          peakTol,
                          57,
                          maxModMass,
                          minModMass,
                          enforceEndpeaks);
        } else {
          apb->computeAlignment(cSpecRev,
                               dbSpec,
                               m_db->getSequence(protIdx),
                               protIdx,
                               1,
                               itrMap->second,
                               minNumMatchPeaks,
                               maxDbGapAas,
                               pmTol,
                               peakTol,
                               enforceEndpeaks);
        }
      } // for (; itr != itrEnd; itr++) {

      elapsed = ((float)(clock() - startTime)) / CLOCKS_PER_SEC;
      if (DEBUG_SPECPROTALIGN_TIME) DEBUG_MSG(elapsed << " seconds elapsed");
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(cSpecRev.psmList.size());
      if (DEBUG_SPECPROTALIGN_TIME) DEBUG_MSG((cSpecRev.psmList.size() == 0 ? 0 : elapsed / (float)cSpecRev.psmList.size()) << " per protein");

      litr = cSpecRev.psmList.begin();
      litrEnd = cSpecRev.psmList.end();
      for (; litr != litrEnd; litr++) {
        (*litr)->m_scanNum = cSpecRev.scan;
        (*litr)->m_protein = m_db->getID((*litr)->m_dbIndex);
        m_allPsms.push_back(*litr);
      }

      DEBUG_MSG("Matching Complete");

      sort(m_allPsms.m_psmSet.begin(), m_allPsms.m_psmSet.end(), PsmDbIndexAnnoSort);

      // -----------------------------------------------------
      // Determine the top PSMs to pass on
      // -----------------------------------------------------
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(cSpec.psmList.size());
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(cSpecRev.psmList.size());

      // Find the top score
      float bestScore = -(float)INT_MAX;
      litr = cSpec.psmList.begin();
      litrEnd = cSpec.psmList.end();
      for (; litr != litrEnd; litr++) {
        if (DEBUG_SPECPROTALIGN) DebugPsm(*litr);
        if ((*litr)->m_score > bestScore) {
          bestScore = (*litr)->m_score;
        }
      }
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(bestScore);

      bool bestScoreForward = true;
      litr = cSpecRev.psmList.begin();
      litrEnd = cSpecRev.psmList.end();
      for (; litr != litrEnd; litr++) {
        if (DEBUG_SPECPROTALIGN) DebugPsm(*litr);
        if ((*litr)->m_score > bestScore) {
          bestScore = (*litr)->m_score;
          bestScoreForward = false;
        }
      }
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(bestScore);

      float thresholdScore = bestScore > 0 ? bestScore * thresholdPercent
          : bestScore / thresholdPercent;
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(thresholdScore);
      cSpec.psmList.remove_if(LessThanPredicate(thresholdScore));
      cSpecRev.psmList.remove_if(LessThanPredicate(thresholdScore));
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(cSpec.psmList.size());
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(cSpecRev.psmList.size());
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(bestScoreForward);
      if (bestScoreForward) {
        // "Un-round" all the spectra to put them back like they were before
        cSpec.roundPeaks(1.0 / AA_ROUNDING);
        if (DEBUG_SPECPROTALIGN) DEBUG_VAR(cSpec.psmList.size());
        // Add all the reverse PSMs
        for (list<psmPtr>::iterator iter = cSpecRev.psmList.begin(); iter
            != cSpecRev.psmList.end(); iter++) {
          cSpec.psmList.push_back(*iter);
        }
        if (DEBUG_SPECPROTALIGN) DEBUG_VAR(cSpec.psmList.size());
        m_matchedSpectraAll->push_back(cSpec);
      }
      else
      {
        // "Un-round" all the spectra to put them back like they were before
        cSpecRev.roundPeaks(1.0 / AA_ROUNDING);
        if (DEBUG_SPECPROTALIGN) DEBUG_VAR(cSpecRev.psmList.size());
        // Add all the forward PSMs
        for (list<psmPtr>::iterator iter = cSpec.psmList.begin(); iter
            != cSpec.psmList.end(); iter++) {
          cSpecRev.psmList.push_back(*iter);
        }
        if (DEBUG_SPECPROTALIGN) DEBUG_VAR(cSpecRev.psmList.size());
        m_matchedSpectraAll->push_back(cSpecRev);
      }

    } // for (specIdx = startIdx; specIdx <= endIdx; specIdx++)

    if (DEBUG_SPECPROTALIGN) DEBUG_VAR(m_matchedSpectraAll->size());

    // -----------------------------------------------------
    // We could have multiple hits on one protein
    // So we perform a sort() and unique() to eliminate them
    // -----------------------------------------------------
    for (int i = 0; i < m_matchedSpectraAll->size(); i++) {
      if (DEBUG_SPECPROTALIGN)DEBUG_VAR(i);
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR((*m_matchedSpectraAll)[i].psmList.size());
      (*m_matchedSpectraAll)[i].psmList.sort(PsmDbIndexAnnoSort);
      (*m_matchedSpectraAll)[i].psmList.unique(PsmDbIndexAnnoUnique);
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR((*m_matchedSpectraAll)[i].psmList.size());
    }

    // -----------------------------------------------------
    // Set backwards pointers/scans for PSM--->Spectrum link
    // -----------------------------------------------------
    list<psmPtr>::iterator iter;
    for (int i = 0; i < m_matchedSpectraAll->size(); i++) {
      (*m_matchedSpectraAll)[i].scan = startIdx + i + 1;
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR((*m_matchedSpectraAll)[i].scan);
      for (iter = (*m_matchedSpectraAll)[i].psmList.begin(); iter
          != (*m_matchedSpectraAll)[i].psmList.end(); iter++) {
        (*iter)->m_spectrum = &(*m_matchedSpectraAll)[i];
        //if (DEBUG_SPECPROTALIGN) DEBUG_VAR((*iter)->m_matchedPeaks.size());

        string dbSeqString = m_db->getSequence((*iter)->m_dbIndex);
        if (tagParsimony) {
          string annotation;
          (*iter)->getAnnotationFromMatchedPeaks(m_db->getMassesSpec((*iter)->m_dbIndex),
                                                 dbSeqString,
                                                 annotation);
          (*iter)->m_annotation = (*iter)->m_origAnnotation;
          (*iter)->m_annotation = annotation;
          if (DEBUG_SPECPROTALIGN) DEBUG_VAR((*iter)->m_annotation);
        }

        // For now copy annotation to original annotation(we'll fix this later)
        (*iter)->m_origAnnotation = (*iter)->m_annotation;
        //        (*iter)->m_annotation = stringAnnotation;
      }
    }

    // -----------------------------------------------------
    // Setup mods for Spectral Probability Computation
    // -----------------------------------------------------
    vector<pair<unsigned int, bool> > ntermMods;
    ntermMods.resize(2);
    ntermMods[0].first = 42;
    ntermMods[0].second = false; // N-term acetylation
    ntermMods[0].first = 43;
    ntermMods[0].second = false;
    ntermMods[1].first = 111;
    ntermMods[1].second = true; // N-term PyroGlu

    // Turn known mods from penalty matrix into known mods for spectral probabilities
    vector<unsigned int> mods;
    if (m_penaltyMatrixMods != 0) {
      const map<string, set<float> > & knownMods =
          m_penaltyMatrixMods->getKnownMods();
      map<string, set<float> >::const_iterator itr = knownMods.begin();
      map<string, set<float> >::const_iterator itrEnd = knownMods.end();
      for (; itr != itrEnd; itr++) {
        string stringAA = itr->first;
        //DEBUG_VAR(stringAA);
        float aaMass = m_penaltyMatrixMods->getMass(itr->first);
        //DEBUG_VAR(aaMass);
        if (aaMass == 0)
          continue; // Sanity check
        const set<float> & setMods = itr->second;
        set<float>::const_iterator itrSet = setMods.begin();
        set<float>::const_iterator itrSetEnd = setMods.end();
        for (; itrSet != itrSetEnd; itrSet++) {
          float mod = *itrSet;
          //DEBUG_VAR(mod);
          float totalMass = aaMass + mod;
          //DEBUG_VAR(totalMass);
          mods.push_back(totalMass);
        }
      }
      //DEBUG_VAR(mods.size());
    }

    AAJumps aminoacids(1); // Used to calculate spectral probabilities
    bool isReversed;

    // -----------------------------------------------------
    // Compute Spectral Probability scores
    // -----------------------------------------------------
    clock_t startTime = clock();

    if (DEBUG_SPECPROTALIGN_SPECTRA) DEBUG_VAR(m_matchedSpectraAll->size());
    for (int i = 0; i < m_matchedSpectraAll->size(); i++) {
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(i);

      Spectrum & origSpec = (*m_matchedSpectraAll)[i];

#if DEBUG_SPECPROTALIGN
      origSpec.computeSpectralProbabilities(ntermMods, mods, peakTol, false);
      list<psmPtr>::iterator itr = origSpec.psmList.begin();
      list<psmPtr>::iterator itrEnd = origSpec.psmList.end();
      for (; itr != itrEnd; itr++) {
        DEBUG_MSG("OR  " << (*itr)->m_annotation.c_str() << "  " << (*itr)->m_score << "  " << (*itr)->m_pValue);
      }
      //==========================wc==============
#endif

      if (!penaltyAlign) {

        origSpec.computeSpectralProbabilities(ntermMods, mods, peakTol, false);

      } else if (origSpec.psmList.size() != 0 && usePrmSprinkling && m_prmSpectra != 0x0 && m_prmSpectra->size() != 0) {

        Spectrum tempUnsprinkled = origSpec;

        Spectrum tempSprinkled = (*m_prmSpectra)[i];
        SprinkleSpectrumWithPrms(origSpec, tempSprinkled, prmMultiplier);
        m_matchedSpectraAllSprinkled.push_back(tempSprinkled);

        if (tempUnsprinkled.psmList.size() != 0 && roundAnno) getBestRoundedAnnos(tempUnsprinkled, peakTol, roundAnnoMax);
        computeAllGapAnnotations(*apb, tempUnsprinkled);
        tempUnsprinkled.computeSpectralProbabilities(ntermMods, mods, peakTol, false);
	
        if (tempSprinkled.psmList.size() != 0 && roundAnno) getBestRoundedAnnos(tempSprinkled, peakTol, roundAnnoMax);
        computeAllGapAnnotations(*apb, tempSprinkled);
        tempSprinkled.computeSpectralProbabilities(ntermMods, mods, peakTol, false);

        if (DEBUG_SPECPROTALIGN) DEBUG_VAR(tempUnsprinkled.psmList.size());
        if (DEBUG_SPECPROTALIGN) DEBUG_VAR(tempSprinkled.psmList.size());

        // Copy all the pValues back to the original spectrum
        list<psmPtr>::iterator itr = origSpec.psmList.begin();
        list<psmPtr>::iterator itrEnd = origSpec.psmList.end();
        list<psmPtr>::iterator itr2 = tempUnsprinkled.psmList.begin();
        list<psmPtr>::iterator itr3 = tempSprinkled.psmList.begin();
        for (; itr != itrEnd; itr++, itr2++, itr3++) {
          if (DEBUG_SPECPROTALIGN) DEBUG_MSG("Unsprinkled pvalue = " << (*itr2)->m_pValue);
          if (DEBUG_SPECPROTALIGN) DEBUG_MSG("Sprinkled pvalue = " << (*itr3)->m_pValue);
          if ((*itr2)->m_pValue > (*itr3)->m_pValue) {
            (*itr)->m_pValue = (*itr2)->m_pValue;
            (*itr)->m_annotation = (*itr2)->m_annotation;
            if (DEBUG_SPECPROTALIGN) DEBUG_MSG("UN  " << (*itr)->m_annotation.c_str() << "  " << (*itr)->m_score << "  " << (*itr)->m_pValue);
          } else {
            (*itr)->m_pValue = (*itr3)->m_pValue;
            (*itr)->m_annotation = (*itr3)->m_annotation;
            if (DEBUG_SPECPROTALIGN) DEBUG_MSG("SP  " << (*itr)->m_annotation.c_str() << "  " << (*itr)->m_score << "  " << (*itr)->m_pValue);
          }
          vector<float> modifications;
          vector<unsigned int> positions;
          (*itr)->getModificationsAndPositions(modifications, positions);
          (*itr)->m_numMods = modifications.size();
        }

      } else {

        m_matchedSpectraAllSprinkled.push_back(origSpec);

        if (DEBUG_SPECPROTALIGN) {
          list<psmPtr>::iterator itr = origSpec.psmList.begin();
          list<psmPtr>::iterator itrEnd = origSpec.psmList.end();
          for (; itr != itrEnd; itr++) {
           float origScore = MatchSpecToPeptide(origSpec,
                                          (*itr)->m_annotation.c_str(),
                                          peakTol,
                                          0,
                                          false,
                                          &isReversed,
                                          &aminoacids);
            if (DEBUG_SPECPROTALIGN) DEBUG_MSG("OR  " << (*itr)->m_annotation.c_str() << "  " << origScore << "  " << (*itr)->m_pValue);
          }
        }

        if (origSpec.psmList.size() != 0 && roundAnno) getBestRoundedAnnos(origSpec, peakTol, roundAnnoMax);
        computeAllGapAnnotations(*apb, origSpec);
        origSpec.computeSpectralProbabilities(ntermMods, mods, peakTol, false);

        if (DEBUG_SPECPROTALIGN) {
          list<psmPtr>::iterator itr = origSpec.psmList.begin();
          list<psmPtr>::iterator itrEnd = origSpec.psmList.end();
          for (; itr != itrEnd; itr++) {
            DEBUG_MSG("FI  " << (*itr)->m_annotation.c_str() << "  " << (*itr)->m_score << "  " << (*itr)->m_pValue);
          }
        }
      }

    }

    float elapsed = ((float)(clock() - startTime)) / CLOCKS_PER_SEC;
    if (DEBUG_SPECPROTALIGN_TIME) DEBUG_MSG(elapsed << " seconds elapsed");

    return true;
  }

  // -------------------------------------------------------------------------

  bool ExecSpecProtAlign::loadInputData(void)
  {
    if (ownInput) {
      if (!m_inputSpectra)
        m_inputSpectra = new SpecSet;
      if (!m_prmSpectra)
        m_prmSpectra = new SpecSet;
      if (!m_db)
        m_db = new DB_fasta;
    }
    m_inputSpectra->resize(0);

    if (ownOutput) {
      if (!m_matchedSpectraAll)
        m_matchedSpectraAll = new SpecSet;
    }
    m_matchedSpectraAll->resize(0);

    // Load amino acid masses
    AAJumps jumps(1);
    if (m_params.exists("AMINO_ACID_MASSES")) {
      jumps.loadJumps(m_params.getValue("AMINO_ACID_MASSES").c_str(), true);
    }

    if (m_params.exists("INPUT_SPECS_PKLBIN") and m_params.exists("INPUT_PSM")) {
      if (m_inputSpectra->loadPklBin(m_params.getValue("INPUT_SPECS_PKLBIN").c_str(),
                                     m_params.getValue("INPUT_PSM").c_str()) < 0) {
        ERROR_MSG("Error reading input spectra files.");
        return false;
      }
    }

    if (m_params.exists("INPUT_PRM_PKLBIN")) {
      if (m_prmSpectra->loadPklBin(m_params.getValue("INPUT_PRM_PKLBIN").c_str()) < 0) {
        ERROR_MSG("Error reading input spectra files.");
        return false;
      }
    }

    if (!m_params.exists("INPUT_FASTA")) {
      ERROR_MSG("Parameters are incomplete. INPUT_FASTA is missing.");
      return false;
    } else if (m_db->Load(m_params.getValue("INPUT_FASTA").c_str()) <= 0) {
      ERROR_MSG("Error reading database sequences from "
          << m_params.getValue("INPUT_FASTA"));
      return false;
    }

    //---------------------------------------------------------------------------
    // Load penalty matrices for new alignment
    //---------------------------------------------------------------------------

    bool penaltyAlign = (bool)m_params.exists("PENALTY_ALIGNMENT")
        ? ((int)m_params.getValueInt("PENALTY_ALIGNMENT") ? 1 : 0) : 0;

    if (penaltyAlign) {
      float resolution = m_params.getValueDouble("ALIGNMENT_RESOLUTION", 0.1);
      float unknownPenalty = m_params.getValueDouble("PENALTY_ALIGNMENT_UNKNOWN_PENALTY", 1.0);
      float unknownMultiplier = m_params.getValueDouble("PENALTY_ALIGNMENT_UNKNOWN_MULTIPLIER", 2.0);

      if (!m_penaltyMatrixBlosum) {
        m_penaltyMatrixBlosum = new PenaltyMatrix(jumps, resolution, unknownPenalty, unknownMultiplier);
      }
      if (!m_penaltyMatrixMods) {
        m_penaltyMatrixMods = new PenaltyMatrix(jumps, resolution, unknownPenalty, unknownMultiplier);
      }

      if (!m_params.exists("BLOSUM_PENALTY_FILE"))
      {
        ERROR_MSG("Parameters are incomplete. BLOSUM_FILE is missing and penalty alignment was specified.");
        return false;
      } else {
        string blosumFileName = m_params.getValue("BLOSUM_PENALTY_FILE");
        if (!m_penaltyMatrixBlosum->load(blosumFileName)) {
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
        if (!m_penaltyMatrixMods->load(modFileName)) {
          ERROR_MSG("Error loading mod penalties from " << modFileName);
          return false;
        }
      }

    } // if (penaltyAlign)

    return true;
  }

  // -------------------------------------------------------------------------

  bool ExecSpecProtAlign::saveOutputData(void)
  {
    bool penaltyAlign = (bool)m_params.exists("PENALTY_ALIGNMENT")
        ? ((int)m_params.getValueInt("PENALTY_ALIGNMENT") ? 1 : 0) : 0;
    DEBUG_VAR(penaltyAlign);

    if (m_params.exists("OUTPUT_MATCHED_PEAKS_IDX_ALL")) {

      AAJumps aaJumps(1);
      if (m_params.exists("AMINO_ACID_MASSES")) {
        DEBUG_MSG("Loading amino acid masses...");
        if (!aaJumps.loadJumps(m_params.getValue("AMINO_ACID_MASSES").c_str(),
                               false))
        {
          ERROR_MSG("Error reading input amino acid mass file " << m_params.getValue("AMINO_ACID_MASSES"));
          return false;
        }
      }

      bool tagParsimony = (bool)m_params.exists("MAX_PARSIMONY")
          ? ((int)m_params.getValueInt("MAX_PARSIMONY") ? 1 : 0) : 0;
      DEBUG_VAR(tagParsimony);

      DEBUG_MSG("Saving matched peaks all...");
      SpecSet tempMatchedPeaks;
      tempMatchedPeaks.resize(m_matchedSpectraAll->size());
      for (int i = 0; i < m_matchedSpectraAll->size(); i++) {
        //DEBUG_VAR(i);
        //DEBUG_VAR((*m_matchedSpectraAll)[i].psmList.size());
        if ((*m_matchedSpectraAll)[i].psmList.size() != 0) {
          list<psmPtr>::iterator litr =
              (*m_matchedSpectraAll)[i].psmList.begin();

          if (tagParsimony) {
            //DEBUG_VAR((*litr)->m_annotation);
            //if (penaltyAlign) DEBUG_VAR((*litr)->m_origAnnotation);
            Spectrum locDbSpec = m_db->getMassesSpec((*litr)->m_dbIndex);
            locDbSpec.roundPeaks();

            (*litr)->getMatchedPeaksFromAnnotation(locDbSpec, aaJumps);
            //DEBUG_VAR((*litr)->m_annotation);
          }

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
      tempMatchedPeaks.savePklBin(m_params.getValue("OUTPUT_MATCHED_PEAKS_IDX_ALL").c_str());
    }

    DEBUG_VAR(m_params.getValue("OUTPUT_MATCHED_PEAKS_IDX_ALL").c_str());
    DEBUG_VAR(m_matchedSpectraAll->size());
    if (m_matchedSpectraAll and m_params.exists("OUTPUT_MATCHED_SPECS_ALL")) {
      if (m_params.exists("OUTPUT_PSM_ALL")
          && m_params.exists("OUTPUT_MATCHED_PEAKS_IDX_ALL")) {
        DEBUG_MSG("Saving matched specs all (3 files)...");
        m_matchedSpectraAll->savePklBin(m_params.getValue("OUTPUT_MATCHED_SPECS_ALL").c_str(),
                                        m_params.getValue("OUTPUT_PSM_ALL").c_str(),
                                        m_params.getValue("OUTPUT_MATCHED_PEAKS_IDX_ALL").c_str());

      } else if (m_params.exists("OUTPUT_PSM_ALL")) {
        DEBUG_MSG("Saving matched specs all (2 files)...");
        DEBUG_VAR(m_params.getValue("OUTPUT_MATCHED_SPECS_ALL").c_str());
        DEBUG_VAR(m_params.getValue("OUTPUT_PSM_ALL").c_str());
        m_matchedSpectraAll->savePklBin(m_params.getValue("OUTPUT_MATCHED_SPECS_ALL").c_str(),
                                        m_params.getValue("OUTPUT_PSM_ALL").c_str());

      } else {
        DEBUG_MSG("Saving matched specs all (1 file)...");
        m_matchedSpectraAll->savePklBin(m_params.getValue("OUTPUT_MATCHED_SPECS_ALL").c_str());
      }
    }

    if (m_params.exists("OUTPUT_PSM_ALL")) {
      m_allPsms.saveToFile(string(m_params.getValue("OUTPUT_PSM_ALL") + ".debug").c_str());
    }
    DEBUG_VAR(m_params.getValue("OUTPUT_SPECTRA_SPRINKLED"));
    if (m_matchedSpectraAll && m_params.exists("OUTPUT_SPECTRA_SPRINKLED")) {
      m_matchedSpectraAllSprinkled.savePklBin(m_params.getValue("OUTPUT_SPECTRA_SPRINKLED").c_str());
      DEBUG_MSG("Saving matched specs all (1 file)...");
    }
    DEBUG_TRACE;
    return true;
  }

  // -------------------------------------------------------------------------

  bool ExecSpecProtAlign::saveInputData(std::vector<std::string> & filenames)
  {
    string dataDir = m_params.getValue("GRID_DATA_DIR");
    if (dataDir.empty()) {
      dataDir = ".";
    }
    string baseDirectory = dataDir + "/";
    string baseFilename = baseDirectory + getName();

    string specFilename = baseFilename + "_stars.pklbin";
    string psmFilename = baseFilename + "_tags.txt";
    DEBUG_MSG("Saving " << specFilename);
    m_inputSpectra->savePklBin(specFilename.c_str(), psmFilename.c_str());

    string prmFilename = baseDirectory + "prm_spectra.pklbin";
    if (m_prmSpectra != 0x0) {
      if (!fileExists(prmFilename))
      {
        DEBUG_MSG("Saving " << prmFilename );
        m_prmSpectra->savePklBin(prmFilename.c_str());
      } else {
       DEBUG_MSG("Not Saving " << prmFilename  << " (already exists)");
      }
    }

    DEBUG_VAR(m_db);
    string dbFilename = baseDirectory + "db.fasta";
    if (!fileExists(dbFilename)) {
      DEBUG_MSG("Saving " << dbFilename );
      m_db->Save(dbFilename.c_str());
    } else {
      DEBUG_MSG("Not Saving " << dbFilename << " (already exists)");
    }

    string aaFilename = baseDirectory + "aminoacids.txt";
    if (!fileExists(aaFilename)) {
      DEBUG_MSG("Saving " << aaFilename);
      m_penaltyMatrixMods->saveAminoAcids(aaFilename);
    } else {
      DEBUG_MSG("Not Saving " << aaFilename << " (already exists)");
    }

    string knownmodsFilename = baseDirectory + "knownmods.txt";
    if (!fileExists(knownmodsFilename)) {
      DEBUG_MSG("Saving " << knownmodsFilename);
      m_penaltyMatrixMods->saveKnownMods(knownmodsFilename);
    } else {
      DEBUG_MSG("Not Saving " << knownmodsFilename << " (already exists)");
    }

    string modPenaltiesFilename = baseDirectory + "modpenalties.txt";
    if (!fileExists(modPenaltiesFilename)) {
      DEBUG_MSG("Saving " << modPenaltiesFilename);
      m_penaltyMatrixMods->saveMatrix(modPenaltiesFilename);
    } else {
      DEBUG_MSG("Not Saving " << modPenaltiesFilename << " (already exists)");
    }

    string blossumFilename = baseDirectory + "blosumpenalties.txt";
    if (!fileExists(blossumFilename)) {
      DEBUG_MSG("Saving " << blossumFilename);
      m_penaltyMatrixBlosum->saveMatrix(blossumFilename);
    } else {
      DEBUG_MSG("Not Saving " << blossumFilename << " (already exists)");
    }

    m_params.setValue("INPUT_SPECS_PKLBIN", specFilename);
    if (m_prmSpectra != 0x0) {
      m_params.setValue("INPUT_PRM_PKLBIN", prmFilename );
    }
    m_params.setValue("INPUT_PSM", psmFilename);
    m_params.setValue("INPUT_FASTA", dbFilename);

    m_params.setValue("AMINO_ACID_MASSES", aaFilename);
    m_params.setValue("KNOWN_MODS_FILE", knownmodsFilename);
    m_params.setValue("MODS_PENALTY_FILE", modPenaltiesFilename);
    m_params.setValue("BLOSUM_PENALTY_FILE", blossumFilename);

    string paramFilename = baseFilename + ".params";
    m_params.writeToFile(paramFilename);

    filenames.push_back(paramFilename); // Parameter file MUST be first in vector
    filenames.push_back(specFilename);
    filenames.push_back(psmFilename);
    filenames.push_back(dbFilename);
    filenames.push_back(modPenaltiesFilename);
    filenames.push_back(aaFilename);
    filenames.push_back(knownmodsFilename);
    filenames.push_back(blossumFilename);

    return true;
  }
  

  // -------------------------------------------------------------------------

  bool ExecSpecProtAlign::loadOutputData(void)
  {
    if (m_matchedSpectraAll == 0x0) {
      ownOutput = true;
      m_matchedSpectraAll = new SpecSet;
    }

    if (m_params.exists("OUTPUT_MATCHED_SPECS_ALL")) {
      m_matchedSpectraAll->loadPklBin(m_params.getValue("OUTPUT_MATCHED_SPECS_ALL").c_str(),
                                      m_params.getValue("OUTPUT_PSM_ALL").c_str(),
                                      m_params.getValue("OUTPUT_MATCHED_PEAKS_IDX_ALL").c_str());
      DEBUG_VAR(m_matchedSpectraAll->size());				      
    }

    DEBUG_VAR(m_params.getValue("OUTPUT_SPECTRA_SPRINKLED"));
    if (m_params.exists("OUTPUT_SPECTRA_SPRINKLED")) {
      m_matchedSpectraAllSprinkled.loadPklBin(m_params.getValue("OUTPUT_SPECTRA_SPRINKLED").c_str());
      DEBUG_VAR(m_matchedSpectraAllSprinkled.size());				      
    }

    DEBUG_TRACE;
    return true;
  }

  // -------------------------------------------------------------------------

  vector<ExecBase *> const & ExecSpecProtAlign::split(int numSplit)
  {
    DEBUG_VAR(numSplit);

    if (numSplit < 2) {
      DEBUG_MSG("Number split [" << numSplit << "] must be at least 2");
      return m_subModules;
    }

    DEBUG_VAR(m_inputSpectra);
    int spectraSize = m_inputSpectra->size();
    DEBUG_VAR(spectraSize);
    if (spectraSize == 0) {
      DEBUG_MSG("Must have at least one spectra");
      return m_subModules;
    }

    int startBaseIdx;
    if (m_params.exists("IDX_START")) {
      startBaseIdx = max(0, m_params.getValueInt("IDX_START"));
    } else {
      startBaseIdx = 0;
    }
    DEBUG_VAR(startBaseIdx);
    int endBaseIdx;
    if (m_params.exists("IDX_END")) {
      endBaseIdx = max(0, m_params.getValueInt("IDX_END"));
      endBaseIdx = min(endBaseIdx, (int)spectraSize - 1);
    } else {
      endBaseIdx = spectraSize - 1;
    }
    DEBUG_VAR(endBaseIdx);

    int numSpectraPerSplit = (endBaseIdx - startBaseIdx + 1) / numSplit;
    DEBUG_VAR(numSpectraPerSplit);
    int extraSpectra = (endBaseIdx - startBaseIdx + 1) % numSplit;
    DEBUG_VAR(extraSpectra);

    int indexStart = startBaseIdx;
    int indexEnd = startBaseIdx + numSpectraPerSplit - 1;
    if (extraSpectra > 0) {
      indexEnd++;
      extraSpectra--;
    }

    for (int i = 0; i < numSplit; i++) {
      if (startBaseIdx >= spectraSize) {
        break;
      }

      //DEBUG_VAR(i);

      // Copy the parameters
      ParameterList childParams(m_params);
      // Set the start and end indices
      char buf[128];
      sprintf(buf, "%d", indexStart);
      childParams.setValue("IDX_START", buf);
      sprintf(buf, "%d", indexEnd);
      childParams.setValue("IDX_END", buf);

      SpecSet * starSpecSet = new SpecSet;
      for (int iSpec = indexStart; iSpec <= indexEnd; iSpec++) {
        starSpecSet->push_back((*m_inputSpectra)[iSpec]);
      }

      sprintf(buf, "%d", (*starSpecSet)[0].scan );
      childParams.setValue("SCAN_FIRST", buf);
      sprintf(buf, "%d", (*starSpecSet)[starSpecSet->size() - 1].scan);
      childParams.setValue("SCAN_LAST", buf);

      // Make a clone of this module
      ExecBase * theClone = new ExecSpecProtAlign(childParams,
                                                  starSpecSet,
                                                  m_prmSpectra,
                                                  m_db,
                                                  m_penaltyMatrixBlosum,
                                                  m_penaltyMatrixMods);

      // Give it a new name based on the split
      theClone->setName(makeName(m_name, i));

      // Have to set up the output files also so the params will be correct on reload
      string dataDir = m_params.getValue("GRID_DATA_DIR");
      if (dataDir.empty()) {
        dataDir = ".";
      }
      string baseName = dataDir + "/" + theClone->getName();
      //DEBUG_VAR(baseName);

      theClone->m_params.setValue("OUTPUT_MATCHED_SPECS_ALL", baseName
          + "_contigs_all.pklbin");
      theClone->m_params.setValue("OUTPUT_PSM_ALL", baseName + "_psm_all.txt");
      theClone->m_params.setValue("OUTPUT_MATCHED_PEAKS_IDX_ALL", baseName
          + "_midx_all.pklbin");
      //theClone->m_params.removeParam("OUTPUT_MATCHED_PEAKS_IDX_ALL");
      theClone->m_params.setValue("OUTPUT_SPECTRA_SPRINKLED", baseName
            + "_sprinkled.pklbin");

      std::string suffix("");
      char bufSplit[128];
      sprintf(bufSplit, "%d", i + 1);
      theClone->m_params.setValue("NUM_SPLIT", bufSplit);

      m_subModules.push_back(theClone);

      indexStart = indexEnd + 1;
      indexEnd = indexStart + numSpectraPerSplit - 1;
      if (extraSpectra > 0) {
        indexEnd++;
        extraSpectra--;
      }

    } // for (int i = 0; i < numSplit; i++)

    DEBUG_VAR(m_subModules.size());

    return m_subModules;
  }

  // -------------------------------------------------------------------------

  bool ExecSpecProtAlign::merge(void)
  {
    if (m_matchedSpectraAll == 0x0) {
      ownOutput = true;
      m_matchedSpectraAll = new SpecSet;
    }

    DEBUG_VAR(m_subModules.size());
    int iSpec = 0;
    for (int i = 0; i < m_subModules.size(); i++) {
      ExecSpecProtAlign * espa = (ExecSpecProtAlign*)m_subModules[i];
      if (DEBUG_SPECPROTALIGN_MERGE) DEBUG_VAR(espa);
      if (DEBUG_SPECPROTALIGN_MERGE) DEBUG_VAR(espa->m_matchedSpectraAll->size());
      for (int j = 0; j < espa->m_matchedSpectraAll->size(); j++) {
        if (DEBUG_SPECPROTALIGN_MERGE) DEBUG_VAR(iSpec);
        m_matchedSpectraAll->push_back(espa->m_matchedSpectraAll->operator[](j));
        if (DEBUG_SPECPROTALIGN_MERGE) DEBUG_VAR((*m_matchedSpectraAll)[iSpec].size());
        if (DEBUG_SPECPROTALIGN_MERGE) DEBUG_VAR((*m_matchedSpectraAll)[iSpec].psmList.size());
        iSpec++;
      }
    }

    // Set backwards pointers/scans for PSM--->Spectrum link
    list<psmPtr>::iterator iter;
    for (int i = 0; i < m_matchedSpectraAll->size(); i++) {
      if (DEBUG_SPECPROTALIGN_MERGE) DEBUG_VAR((*m_matchedSpectraAll)[i].scan);
      //(*m_matchedSpectraAll)[i].scan = i + 1;
      for (iter = (*m_matchedSpectraAll)[i].psmList.begin(); iter
          != (*m_matchedSpectraAll)[i].psmList.end(); iter++) {
        (*iter)->m_spectrum = &(*m_matchedSpectraAll)[i];
        (*iter)->m_scanNum = (*m_matchedSpectraAll)[i].scan;
        (*iter)->m_spectrumFile.clear();
        if (DEBUG_SPECPROTALIGN_MERGE) DEBUG_VAR((*iter)->m_origAnnotation);
        if (DEBUG_SPECPROTALIGN_MERGE) DEBUG_VAR((*iter)->m_annotation);
      }
    }

    DEBUG_VAR(m_subModules.size());
    for (int i = 0; i < m_subModules.size(); i++) {
      ExecSpecProtAlign * espa = (ExecSpecProtAlign*)m_subModules[i];
      if (DEBUG_SPECPROTALIGN_MERGE) DEBUG_VAR(espa->m_matchedSpectraAllSprinkled.size());
      for (int j = 0; j < espa->m_matchedSpectraAllSprinkled.size(); j++) {
        m_matchedSpectraAllSprinkled.push_back(espa->m_matchedSpectraAllSprinkled[j]);
      }
    }

    return true;
  }

  // -------------------------------------------------------------------------

  bool ExecSpecProtAlign::validateParams(std::string & error)
  {
    m_isValid = false;

    VALIDATE_PARAM_EXIST("TOLERANCE_PEAK");
    VALIDATE_PARAM_EXIST("TOLERANCE_PM");

    m_isValid = true;
    return true;
  }

  // -------------------------------------------------------------------------

  void ExecSpecProtAlign::SprinkleSpectrumWithPrms(Spectrum & inSpec, Spectrum & outSpec, float prmMultiplier)
  {
    if (DEBUG_SPECPROTALIGN_SPRINKLE) {
      list<psmPtr>::iterator itr = inSpec.psmList.begin();
      list<psmPtr>::iterator itr_end = inSpec.psmList.end();
      for (; itr != itr_end; itr++) {
        psmPtr psm = *itr;
        if (DEBUG_SPECPROTALIGN_SPRINKLE) DEBUG_MSG(psm->m_scanNum << "\t" << psm->m_annotation << "\t" << psm->m_score << "\t" << psm->m_pValue);
      }
    }

    int scan = inSpec.scan;
    int scanIndex = scan - 1;
    if (m_prmSpectra != 0x0 && m_prmSpectra->size() != 0) {
      outSpec = inSpec;

      float avgIntensityStar = 0.0;
      for (int iPeak = 0; iPeak < outSpec.size(); iPeak++) {
        avgIntensityStar += outSpec[iPeak][1];
        if (DEBUG_SPECPROTALIGN_SPRINKLE) DEBUG_MSG(outSpec[iPeak][0] << "  " << outSpec[iPeak][1]);
      }
      avgIntensityStar /= (float)outSpec.size();
      if (DEBUG_SPECPROTALIGN_SPRINKLE) DEBUG_VAR(avgIntensityStar);

      if (DEBUG_SPECPROTALIGN_SPRINKLE) DEBUG_VAR(scanIndex);
      if (DEBUG_SPECPROTALIGN_SPRINKLE) DEBUG_VAR(m_prmSpectra->size());
      Spectrum tempPrmSpectrum = (*m_prmSpectra)[scanIndex];

      float avgIntensityPrm = 0.0;
      if (DEBUG_SPECPROTALIGN_SPRINKLE) DEBUG_VAR(tempPrmSpectrum.parentMass);
      tempPrmSpectrum.filterLowMassPeaks(57.0);        
      tempPrmSpectrum.filterHighMassPeaks(tempPrmSpectrum.parentMass - 57.0);
      for (int iPeak = 0; iPeak < tempPrmSpectrum.size(); iPeak++) {
        if (DEBUG_SPECPROTALIGN_SPRINKLE) DEBUG_MSG(tempPrmSpectrum[iPeak][0] << "  " << tempPrmSpectrum[iPeak][1]);
        avgIntensityPrm += tempPrmSpectrum[iPeak][1];
      }
      avgIntensityPrm /= (float)tempPrmSpectrum.size();
      if (DEBUG_SPECPROTALIGN_SPRINKLE) DEBUG_VAR(avgIntensityPrm);

      // Remove any peaks that already exist in the star spectrum
      list<int> peaksToRemove;
      for (int iPeak = 0; iPeak < outSpec.size(); iPeak++) {
        int iIndex = tempPrmSpectrum.findPeaks(outSpec[iPeak][0], outSpec.getTolerance(iPeak));
        if (iIndex != -1) {
          peaksToRemove.push_back(iIndex);
        }
      }
      if (DEBUG_SPECPROTALIGN_SPRINKLE) DEBUG_VAR(peaksToRemove.size());
      if (peaksToRemove.size() != 0) {
        tempPrmSpectrum.removePeaks(peaksToRemove);
      }

      // Smash down the intensities of the PRM spectra by desired factor
      for (int iPeak = 0; iPeak < tempPrmSpectrum.size(); iPeak++) {
        tempPrmSpectrum[iPeak][1] *= avgIntensityStar / avgIntensityPrm  * prmMultiplier;
        if (DEBUG_SPECPROTALIGN_SPRINKLE) DEBUG_MSG(tempPrmSpectrum[iPeak][0] << "  " << tempPrmSpectrum[iPeak][1]);
      }

      if (DEBUG_SPECPROTALIGN_SPRINKLE) DEBUG_VAR(outSpec.size());
      outSpec.mergeClosestPeaks(tempPrmSpectrum, 1);

      if (DEBUG_SPECPROTALIGN_SPRINKLE){
        DEBUG_VAR(outSpec.size());
        for (int iPeak = 0; iPeak < outSpec.size(); iPeak++) {
          DEBUG_MSG(outSpec[iPeak][0] << "  " << outSpec[iPeak][1]);
        }
      }
    }
  }

  // -------------------------------------------------------------------------
  void ExecSpecProtAlign::makeAnnotationFromData(string & cleanAnnotation,
                                                 vector<float> & modifications,
                                                 vector<unsigned int> & positions,
                                                 vector<unsigned int> & lengths,
                                                 string & completeAnnotation)
  {
      if (DEBUG_SPECPROTALIGN_ANNO) DEBUG_VAR(cleanAnnotation );

      char modBuf[100];
      int iMod = 0;
      int iLength = 0;

      if (modifications.size() > 0 && positions[0] == 0) {
        completeAnnotation += "[";
        sprintf(modBuf, "%.3f", modifications[iMod]);
        completeAnnotation += modBuf;
        completeAnnotation += "]";
        iMod++;
      }

      for (int i = 0; i < cleanAnnotation.length(); i++) {
        //if (DEBUG_SPECPROTALIGN_ANNO) DEBUG_VAR(i);
        //if (DEBUG_SPECPROTALIGN_ANNO) DEBUG_VAR(iMod);
        if (iMod < modifications.size() && i == positions[iMod] - lengths[iMod]) {
          if (lengths[iMod] == 0) {
            completeAnnotation += "[";
            sprintf(modBuf, "%.3f", modifications[iMod]);
            completeAnnotation += modBuf;
            completeAnnotation += "]";
            iMod++;
            continue;
          } else {
            completeAnnotation += "(";
          }
          iLength=0;
        }
        //if (DEBUG_SPECPROTALIGN_ANNO) DEBUG_VAR(completeAnnotation);
        completeAnnotation += cleanAnnotation[i];
        //if (DEBUG_SPECPROTALIGN_ANNO) DEBUG_VAR(completeAnnotation);
        iLength++;
        //if (DEBUG_SPECPROTALIGN_ANNO) DEBUG_VAR(iLength);
        if (iMod < modifications.size() && i == positions[iMod] - 1) {
          sprintf(modBuf, "%.3f", modifications[iMod]);
          completeAnnotation += ",";
          completeAnnotation += modBuf;
          completeAnnotation += ")";
          iMod++;
        }
        //if (DEBUG_SPECPROTALIGN_ANNO) DEBUG_VAR(completeAnnotation);
        //if (DEBUG_SPECPROTALIGN_ANNO) DEBUG_VAR(iMod);
      }

      if (iMod < modifications.size()) {
        sprintf(modBuf, "%.3f", modifications[iMod]);
        completeAnnotation += ",";
        completeAnnotation += modBuf;
        completeAnnotation += ")";
        iMod++;
      }
 
      //if (DEBUG_SPECPROTALIGN_ANNO) DEBUG_VAR(completeAnnotation);

    return;
  }

  // -------------------------------------------------------------------------
  void ExecSpecProtAlign::computeAllGapAnnotations(AlignmentPenaltyBased & apb, Spectrum & inSpec)
  {
    list<psmPtr>::iterator itr = inSpec.psmList.begin();
    list<psmPtr>::iterator itrEnd = inSpec.psmList.end();
    for (; itr != itrEnd; itr++) {
      string fullAnnotation;
      apb.computeAllGapAnnotations((*itr)->m_annotation, fullAnnotation);
      (*itr)->m_annotation = fullAnnotation;
    }
    return;
  }

  // -------------------------------------------------------------------------
  void ExecSpecProtAlign::getBestRoundedAnnos(Spectrum & inSpec, float peakTol, float roundAnnoMax)
  {
    AAJumps aminoacids(1); // Used to calculate spectral probabilities

    if (DEBUG_SPECPROTALIGN_ANNO) DEBUG_VAR(inSpec.psmList.size());
    // Loop over all PSMs
    list<psmPtr>::iterator itr = inSpec.psmList.begin();
    list<psmPtr>::iterator itrEnd = inSpec.psmList.end();
    for (; itr != itrEnd; itr++) {

      if (DEBUG_SPECPROTALIGN_ANNO) DEBUG_VAR((*itr)->m_annotation);
      string cleanAnnotation; 
      static string aminoAcids("ABCDEFGHIJKLMNOPQRSTUVWXYZ");
      for (int iChar = 0; iChar < (*itr)->m_annotation.length(); iChar++) {
        if (aminoAcids.find_first_of((*itr)->m_annotation[iChar]) != string::npos) {
          cleanAnnotation += (*itr)->m_annotation[iChar];
        }
      }
      if (DEBUG_SPECPROTALIGN_ANNO) DEBUG_VAR(cleanAnnotation );

      vector<float> modifications;
      vector<unsigned int> positions;
      vector<unsigned int> lengths;
      (*itr)->getModificationsAndPositions(modifications, positions, lengths);

      if (DEBUG_SPECPROTALIGN_ANNO){
        for (int i = 0; i < modifications.size(); i++) {
          DEBUG_MSG(i << "  " << modifications[i] << "  " << positions[i] << "  " << lengths[i]);
        }
      }

      vector<float> lowerModifications;
      vector<float> upperModifications;
      vector<unsigned int> lowerPositions;
      vector<unsigned int> lowerLengths;
      vector<unsigned int> upperPositions;
      vector<unsigned int> upperLengths;
      for (int i = 0; i < modifications.size(); i++) {
        if (abs(modifications[i]) < roundAnnoMax) {
          if (modifications[i] < 0) {
            upperModifications.push_back(floor(modifications[i]));
          } else {
            upperModifications.push_back(ceil(modifications[i]));
          }
          upperPositions.push_back(positions[i]);
          upperLengths.push_back(lengths[i]);
        } else {
          lowerModifications.push_back(modifications[i]);
          lowerPositions.push_back(positions[i]);
          lowerLengths.push_back(lengths[i]);
          upperModifications.push_back(modifications[i]);
          upperPositions.push_back(positions[i]);
          upperLengths.push_back(lengths[i]);
        }
      }

      string lowerAnnotation;
      makeAnnotationFromData(cleanAnnotation, 
                             lowerModifications, 
                             lowerPositions, 
                             lowerLengths, 
                             lowerAnnotation);
      if (DEBUG_SPECPROTALIGN_ANNO) DEBUG_VAR(lowerAnnotation);

      string upperAnnotation;
      makeAnnotationFromData(cleanAnnotation, 
                             upperModifications, 
                             upperPositions, 
                             upperLengths, 
                             upperAnnotation);
      if (DEBUG_SPECPROTALIGN_ANNO) DEBUG_VAR(upperAnnotation);

      bool isReversed;
      float lowerMatchScore = MatchSpecToPeptide(inSpec,
                                          lowerAnnotation.c_str(),
                                          peakTol,
                                          0,
                                          false,
                                          &isReversed,
                                          &aminoacids);
      if (DEBUG_SPECPROTALIGN_ANNO) DEBUG_VAR(lowerMatchScore);

      float upperMatchScore = MatchSpecToPeptide(inSpec,
                                          upperAnnotation.c_str(),
                                          peakTol,
                                          0,
                                          false,
                                          &isReversed,
                                          &aminoacids);
      if (DEBUG_SPECPROTALIGN_ANNO) DEBUG_VAR(upperMatchScore);

      if (upperMatchScore > lowerMatchScore) {
        (*itr)->m_annotation = upperAnnotation;
      } else {
        (*itr)->m_annotation = lowerAnnotation;
      }
      if (DEBUG_SPECPROTALIGN_ANNO) DEBUG_VAR((*itr)->m_annotation);

    } // for (; itr != itrEnd; itr++) {

  }


}
