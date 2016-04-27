//

// Module Includes
#include "AlignmentPenaltyBased.h"
#include "CommandLineParser.h"
#include "Constants.h"
#include "FdrPeptide.h"
#include "Logger.h"
#include "ExecMsCluster.h"
#include "ExecGenoMS.h"
#include "ExecAlignment.h"
#include "ExecDeconvoluteMS2.h"
#include "ExecFilterAligns.h"
#include "ExecFilterContigPairs.h"
#include "ExecFilterPairs.h"
#include "ExecFilterStarPairs.h"
#include "ExecReportSpsplot.h"
#include "ExecReportSPSStats.h"
//#include "ExecReportProteinCoverage.h"
#include "ExecAssembly.h"
#include "ExecMergeConvert.h"
#include "ExecMetaAssembly.h"
#include "ExecPrmScoring.h"
#include "ExecProtProtAlign.h"
#include "ExecSpecProtAlignTgtDecoy.h"
#include "ExecTagSearch.h"
#include "ExecHomologyAssembly.h"
#include "ExecMainSpecnets.h"
#include "FileUtils.h"
#include "ParallelThreadedExecution.h"
#include "ParallelPbsExecution.h"
#include "ParallelSgeExecution.h"
#include "ParameterList.h"
#include "PeptideSpectrumMatchSet.h"

// Specnets Includes
#include "utils.h"  // stringSplit
#include "clusters.h"
#include "ClusterData.h"
#include "abruijn.h"
#include "copyright.h"
#include "SpecSet.h"
#include "Specific.h"
#include "StatusFile.h"

// System Includes
#include <fstream>
#include <iostream>
#include <string.h>
#include <string>
#include <unistd.h>
#include <algorithm>
#include <stdlib.h>

using namespace specnets;
using namespace std;

#define DEBUG_SPEC_ALIGN 0

const string DEFAULT_BIN_FILES_FILENAME = "spectra/bin_files.txt";
const string DEFAULT_INPUT_FILES_LIST = "spectra/pklbin_files.txt";
const string DEFAULT_DECONV_FILES_LIST = "spectra/deconv_files.txt";
const string DEFAULT_USER_INPUT_FILES_LIST = "spectra/input_index.txt";
const string DEFAULT_INPUT_FILE_BASE = "spectra/specs_ms_";

const string DEFAULT_INPUT_MAPPING = "spectra/input_mapping.bin";
const string DEFAULT_DECONV_MODEL = "model_isoenv.bin";

const string CLUSTER_MSCLUST = "MSCluster";
const string CLUSTER_PRMS = "PrmClust";

enum Stage
{
  STAGE_BEGIN = 0,
  STAGE_MS2DECONV = 1,
  STAGE_MSCLUSTER = 2,
  STAGE_SCORING = 3,
  STAGE_GENOMS = 4,
  STAGE_FILTERPAIRS = 5,
  STAGE_FILTERALIGNS = 6,
  STAGE_ALIGNMENT = 7,
  STAGE_FILTERSTARPAIRS = 8,
  STAGE_ASSEMBLY = 9,
  STAGE_METAASSEMBLY = 10,
  STAGE_TAGSEARCH = 11,
  STAGE_SPECNETS = 12,
  STAGE_CONTIGPROTALIGN = 13,
  STAGE_SPECPROTALIGN = 14,
  STAGE_PROTPROTALIGN = 15,
  STAGE_HOMOLOGYASSEMBLY = 16,
  STAGE_MERGE = 17,
  STAGE_STATPROTSEQS = 18,
  STAGE_REPORT = 19
};

// System Defines
#define TEST_VALID {          \
  DEBUG_VAR(isValid);         \
  if (!isValid) {             \
    ERROR_MSG(errorString);   \
    return false;             \
  }                           \
}

#define TEST_RETURN(__module, __var) {                                                    \
  if(__var.size() == 0) {                                                               \
    ERROR_MSG("invoking " << __module << ": empty return set \"" << #__var << "\".");     \
  /*  return false; */                                                                     \
  }                                                                                   \
}

#define TEST_RETURN_STATUS(__module) {                            \
  DEBUG_VAR( returnStatus);                                     \
  if (!returnStatus) {                                          \
    ERROR_MSG("invoking " << __module << " exited in error.");    \
    return false;                                              \
  }                                                             \
}

#define TEST_SAVE_OUPUT_DATA(__module) {                                      \
  DEBUG_VAR(returnStatus);                                                  \
  if (!returnStatus) {                                                      \
    ERROR_MSG("Saving output data for " << __module << " exited in error.");  \
    return false;                                                           \
  }                                                                         \
}

//-----------------------------------------------------------------------------

char *curLibPath = NULL;

//-----------------------------------------------------------------------------
void makeSpectrumTagsFromContig(SpecSet & matchedContigs,
                                vector<
                                    vector<sps::tuple<unsigned int, float, bool> > > & outputAssembledShifts,
                                SpecSet & starSpectra,
                                SpecSet & starAlignSpectra,
                                PeptideSpectrumMatchSet & psmSetContigTag)
{
  for (int i = 0; i < matchedContigs.size(); i++)
  {
    if (DEBUG_SPEC_ALIGN)
      DEBUG_VAR(i);
    if (DEBUG_SPEC_ALIGN)
      DEBUG_VAR(matchedContigs[i].psmList.size());
    if (matchedContigs[i].psmList.size() == 0)
    {
      continue;
    }
    list<psmPtr>::iterator itr = matchedContigs[i].psmList.begin();
    list<psmPtr>::iterator itr_end = matchedContigs[i].psmList.end();
    for (; itr != itr_end; itr++)
    {
      psmPtr origPsm = *itr;
      if (DEBUG_SPEC_ALIGN)
        DEBUG_VAR(origPsm->m_scanNum);
      int contigIndex = origPsm->m_scanNum - 1;
      if (DEBUG_SPEC_ALIGN)
        DEBUG_VAR(contigIndex);
      if (DEBUG_SPEC_ALIGN)
        DEBUG_VAR(outputAssembledShifts[contigIndex].size());
      for (int j = 0; j < outputAssembledShifts[contigIndex].size(); j++)
      {
        int spectrumIndex = outputAssembledShifts[contigIndex][j].m0;
        int spectrumMassShift = (int)outputAssembledShifts[contigIndex][j].m1;

        if (DEBUG_SPEC_ALIGN)
          DEBUG_MSG("  " << contigIndex << "  " << spectrumIndex << "  " << spectrumMassShift);

        psmPtr ptrSpectrumPsm(new PeptideSpectrumMatch);
        *ptrSpectrumPsm = *origPsm;

        if (DEBUG_SPEC_ALIGN)
          DEBUG_VAR(ptrSpectrumPsm->m_startMass);
        if (origPsm->m_matchOrientation == 1)
        {
          // If the contig was reversed then we need to "reverse" the starting location for the spectrum within the contig
          if (DEBUG_SPEC_ALIGN)
            DEBUG_VAR(matchedContigs[i].parentMass);
          if (DEBUG_SPEC_ALIGN)
            DEBUG_VAR(starAlignSpectra[spectrumIndex].parentMass);
          ptrSpectrumPsm->m_startMass += matchedContigs[i].parentMass
              - spectrumMassShift - starAlignSpectra[spectrumIndex].parentMass;
          //ptrSpectrumPsm->m_startMass -= -(starAlignSpectra[spectrumIndex].parentMass - matchedContigs[i].parentMass - spectrumMassShift);
        }
        else
        {
          ptrSpectrumPsm->m_startMass += spectrumMassShift;
        }
        if (DEBUG_SPEC_ALIGN)
          DEBUG_VAR(ptrSpectrumPsm->m_startMass);

#if 1
        // Is the spectrum reversed wrt the contig? If so we have to reverse it!
        if (outputAssembledShifts[contigIndex][j].m2 == 1)
        {
          if (DEBUG_SPEC_ALIGN)
            DEBUG_MSG("REVERSING THE REVERSAL FLAG");
          ptrSpectrumPsm->m_matchOrientation = !origPsm->m_matchOrientation;
        }
#endif

        if (DEBUG_SPEC_ALIGN)
          DEBUG_VAR(ptrSpectrumPsm->m_scanNum);
        ptrSpectrumPsm->m_scanNum = starSpectra[spectrumIndex].scan;
        if (DEBUG_SPEC_ALIGN)
          DEBUG_VAR(ptrSpectrumPsm->m_scanNum);

        // Set the spectrum pointer to this spectrum
        ptrSpectrumPsm->m_spectrum = &starSpectra[spectrumIndex];

        psmSetContigTag.push_back(ptrSpectrumPsm);
        starAlignSpectra[spectrumIndex].psmList.push_back(ptrSpectrumPsm);
        if (DEBUG_SPEC_ALIGN)
          DEBUG_TRACE;
      } // for (int j = 0; j < outputAssembledShifts[scanIndex].size(); j++) {
    } // for ( ; itr != itr_end; itr++) {
  } // for (int i = 0; i < matchedContigs.size(); i++) {

  return;
}

//-----------------------------------------------------------------------------
string getCurrentDirectory(void)
{
  string currentWorkindDir;
  currentWorkindDir = getcwd(NULL, 1024);
  return currentWorkindDir;

  //int dummy;
  //dummy = spsSystem("pwd > pwd.tmp");
  //ifstream ifs("pwd.tmp", ios::binary);
  //char buf[1024];
  //ifs.getline(buf, 1024);
  //return buf;
}

//-----------------------------------------------------------------------------
string getProjPath(const ParameterList & pl, const string & addPath)
{
  string projDir = pl.getValue("PROJECT_DIR", "");
  bool relDir = pl.getValueBool("RELATIVE_DIR", true);

  return getPath(projDir, addPath, relDir);
}

//-----------------------------------------------------------------------------
string getCurrentTimeString(void)
{
  time_t rawtime;
  struct tm * timeinfo;
  time(&rawtime);
  timeinfo = localtime(&rawtime);

  return asctime(timeinfo);
}

void addDefaultParameterValues(ParameterList &p)
{
  // Basic parameters
  p.addIfDoesntExist("TOLERANCE_PEAK", "0.4");
  p.addIfDoesntExist("TOLERANCE_PM", "1.5");
  p.addIfDoesntExist("RESOLUTION", "0.1");

  // Preprocessing parameters
  p.addIfDoesntExist("CLUSTER_MIN_SIZE", "1");
  p.addIfDoesntExist("CLUSTER_MODEL", "LTQ_TRYP");
  p.addIfDoesntExist("MSCLUSTER_MIX_PROB ", "0.05");
  p.addIfDoesntExist("INSTRUMENT_TYPE", "IT");
  p.addIfDoesntExist("MIN_SPECTRUM_QUALITY", "0");
  p.addIfDoesntExist("CORRECT_PM", "no");
  p.addIfDoesntExist("GUESS_CHARGE", "no");
  p.addIfDoesntExist("PEPNOVO_OUTDIR", "spectra");

  // Alignment parameters
  p.addIfDoesntExist("AA_DIFF_COUNT", "2");
  p.addIfDoesntExist("MIN_SHIFT", "0");
  p.addIfDoesntExist("MIN_MOD_MASS", "-100");
  p.addIfDoesntExist("MAX_MOD_MASS", "100");
  p.addIfDoesntExist("MAX_NUM_MODS", "1");
  p.addIfDoesntExist("MIN_RATIO", "0.35");

  p.addIfDoesntExist("MAX_PVALUE", "0.045");
  p.addIfDoesntExist("MIN_MATCHED_PEAKS", "6"); // Minimum number of matched peaks to consider a spectral alignment
  p.addIfDoesntExist("MAX_AA_JUMP", "2");
  p.addIfDoesntExist("MIN_OVERLAP_AREA", "0.45");
  p.addIfDoesntExist("PENALTY_PTM", "-200"); // Set to "0" for SpecNets
  p.addIfDoesntExist("PENALTY_SAME_VERTEX", "-1000000");
  p.addIfDoesntExist("FILTER_TRIGS", "no"); // Set to "no" for SpecNets
  p.addIfDoesntExist("PARTIAL_OVERLAPS", "1"); // Set to "0" for SpecNets
  p.addIfDoesntExist("TAGS_FILTER", "");
  p.addIfDoesntExist("TAGS_MATCH_FLANK", "1");
  p.addIfDoesntExist("TAGS_MATCH_COUNT", "2");

  // Comparative Shotgun Protein Sequencing (CSPS) parameters
  //p.addIfDoesntExist("CLUSTALW_EXE_DIR",      "");
  // p.addIfDoesntExist("CLUSTALW_MINSCORE", "250");
  p.addIfDoesntExist("CLUSTALW_MINSCORE", "1000000"); // Disabled by default - activate explicitly (i.e., set to ~250) for cSPS projects
  //  p.addIfDoesntExist("FORCE_REFERENCE",       "-1");

  // De novo sequencing parameters
  p.addIfDoesntExist("SPSPATH_MIN_NUM_PEAKS", "5");
  p.addIfDoesntExist("ADD_ENDPOINTS", "0");
  p.addIfDoesntExist("SPSPATH_MIN_NUM_SPECS", "2");
  p.addIfDoesntExist("PARALLEL_PATHS", "0");
  p.addIfDoesntExist("SPS_MIN_EDGES_TO_COMPONENT", "1");
  p.addIfDoesntExist("MIN_METACONTIG_SCORE", "3.35");
  p.addIfDoesntExist("MIN_METACONTIG_SIZE", "0");

  // tagsearch/matchma parameters
  p.addIfDoesntExist("TAG_LEN", "6");
  p.addIfDoesntExist("DOUBLE_AA_JUMPS", "1");
  p.addIfDoesntExist("MATCH_TAG_FLANKING_MASSES", "0"); // Set to 2 for SpecNets
  p.addIfDoesntExist("MAX_NUM_TAGS", "0");
  p.addIfDoesntExist("MAX_NUM_MODS", "2");
  p.addIfDoesntExist("MIN_MATCHED_PEAKS_DB", "7"); // Minimum number of matched peaks between spectrum/database to accept PSM
  p.addIfDoesntExist("TAG_MATCH_TOP_SCORING_ONLY", "1");

  // Networks parameters (pathproj)
  //  p.addIfDoesntExist("MIN_PERC_EXPINT",   "0.01");
  //  p.addIfDoesntExist("MIN_PERC_TP",       "0.01");

  // Grid parameters
  p.addIfDoesntExist("GRID_TYPE", "sge");
  p.addIfDoesntExist("GRID_NUMNODES", "-1");
  p.addIfDoesntExist("GRID_NUMCPUS", "1");
  p.addIfDoesntExist("GRID_EXE_DIR", "");
  p.addIfDoesntExist("GRID_SGE_EXE_DIR", "");
  p.addIfDoesntExist("GRID_PARAMS", "-l h_vmem=3G");

  string currDir = getCurrentDirectory();
  p.addIfDoesntExist("PROJECT_DIR", currDir);
  p.addIfDoesntExist("RELATIVE_DIR", "1");

  // Reporting parameters
  p.addIfDoesntExist("REPORT_DIR", "report");
  p.addIfDoesntExist("REPORT_DYNAMIC", "1");
  p.addIfDoesntExist("SPS_PROJECTS", "./sps_projects.txt");

}

/**
 * Loads raw MS/MS spectra from user, converts them to pklbin, then saves pklbin
 *   spectra along with lists of filenames to the project spectra directory. ExecMergeConvert
 *   is used to do any pre-processing of input spectra (set tolerances, set activation, etc.)
 * @param ip list of parameters from user
 * @param inputFileNames input file names delimeted by ";"
 * @param loadedSpectra output collection of SpecSets representing the converted spectra
 * @return true if successful, false if not
 */
bool loadInitialData(ParameterList & ip,
                     string inputFileNames,
                     vector<string>& loadedSpectraFiles)
{
  DEBUG_TRACE;
  //---------------------------------
  // Load spectra
  //---------------------------------

  string pklbinPath;

  // List of generated data files
  vector<string> pklbinFileList;
  //
  vector<string> filenames;
  ostringstream cmd;
  stringSplit(inputFileNames, filenames, ";");

  vector<string> set_filenames;

  if (ip.exists("SET_FILENAMES"))
  {
    stringSplit(ip.getValue("SET_FILENAMES"), set_filenames, ";");
    if (set_filenames.size() != filenames.size())
    {
      ERROR_MSG("SET_FILENAMES length (" << set_filenames.size() << ") must match INPUT_SPECS_MS length (" << filenames.size() << ")");
      return false;
    }
  }
  vector<pair<int, int> > loadedIndices;

  loadedSpectraFiles.resize(filenames.size());

  DEBUG_TRACE;
  unsigned int baseIdx, totalSpectrumCount = 0;

  for (unsigned int i = 0; i < filenames.size(); i++)
  {
    DEBUG_VAR (filenames[i]);

    ParameterList loadParams(ip);
    loadParams.setValue("INPUT_SPECTRA", filenames[i]);
    if (set_filenames.size() > 0)
    {
      loadParams.setValue("SET_FILENAMES", set_filenames[i]);
    }
    SpecSet tempSpecs;

    ExecMergeConvert* loader = new ExecMergeConvert(loadParams, &tempSpecs);

    // let ExecMergeConvert take care of loading spectra
    if (!loader->loadInputData())
    {
      delete loader;
      return false;
    }

    // let ExecMergeConvert pre-process spectra (set tolerances, activation, etc)
    if (!loader->invoke())
    {
      delete loader;
      return false;
    }

    cmd.str("");
    // generate sequential name, 1 based index
    cmd << DEFAULT_INPUT_FILE_BASE << i + 1 << ".pklbin";
    pklbinPath = getProjPath(ip, cmd.str());
    DEBUG_VAR(pklbinPath);
    pklbinFileList.push_back(pklbinPath);

    loadedSpectraFiles[i] = pklbinPath;

    delete loader;

    if (!ExecMergeConvert::saveSpecset(pklbinPath, &tempSpecs))
    {
      return false;
    }

    DEBUG_TRACE;

    if(tempSpecs.size() == 0)
    {
      ERROR_MSG("Input spectra set is empty for file " << filenames[i]);
      return false;
    }

    totalSpectrumCount += tempSpecs.size();

  }

  // Check if there are no spectra, in which case it fails
  if(totalSpectrumCount == 0)
  {
    ERROR_MSG("Input spectra set is empty");
    return false;
  }

  if (!ip.exists("SET_FILENAMES"))
  {
    set_filenames = filenames;
  }

  // save input file list
  string aux;
  aux = getProjPath(ip, DEFAULT_USER_INPUT_FILES_LIST);
  if (!writeFileIndex(aux.c_str(), set_filenames))
  {
    return false;
  }

  DEBUG_TRACE;

  // save pklbin file list
  aux = getProjPath(ip, DEFAULT_INPUT_FILES_LIST);
  if (!writeFileIndex(aux.c_str(), pklbinFileList))
  {
    return false;
  }

  DEBUG_TRACE;
  return true;
}

/**
 * Does the merging of separate sets of spectra into one. This sets scan #s in the merged set so
 *  they are unique and saves a mapping to keep track of them.
 * @param ip input list of user parameters
 * @param loadedSpectra input set of SpecSets, represents current working set. This will have size 0 after this function returns
 * @param ms2spectra output set of merged spectra w/ modified scan #s
 * @return true if successful, false if not
 */
bool mergeSeparateSpecs(ParameterList& ip,
                        vector<SpecSet>& loadedSpectra,
                        SpecSet& ms2spectra,
                        vector<vector<int> >& specMapping)
{
  DEBUG_VAR(loadedSpectra.size());

  if (loadedSpectra.size() == 1)
  {
    ms2spectra.clear();
    ms2spectra.swap(loadedSpectra[0]);
    loadedSpectra.resize(0);

    for (int i = 0; i < ms2spectra.size(); i++)
    {
      vector<int> aux;
      aux.push_back(0);
      aux.push_back(i);
      specMapping.push_back(aux);
    }
  }
  else
  {
    // count # of spectra
    unsigned int totalSpectrumCount = 0;
    for (unsigned int i = 0; i < loadedSpectra.size(); i++)
    {
      totalSpectrumCount += loadedSpectra[i].size();
    }
    totalSpectrumCount = 0;

    for (unsigned int i = 0; i < loadedSpectra.size(); i++)
    {

      ms2spectra.swapAppendSpecSet(loadedSpectra[i], false);

      for (unsigned int j = totalSpectrumCount; j < ms2spectra.size(); j++)
      {
        // make sure scan #s are unique
        ms2spectra[j].scan = j + 1;

        vector<int> aux;
        aux.push_back(i);
        aux.push_back(j - totalSpectrumCount);
        specMapping.push_back(aux);
      }
      totalSpectrumCount = ms2spectra.size();
      DEBUG_VAR(totalSpectrumCount);
    }
  }

  // save input mapping between input spectra files and combined spectra files
  string aux3 = getProjPath(ip, DEFAULT_INPUT_MAPPING);
  if (!Save_binArray(aux3.c_str(), specMapping))
  {
    ERROR_MSG("Failed to save to \'" << DEFAULT_INPUT_MAPPING << "\'");
    return false;
  }
  loadedSpectra.resize(0);
  return true;
}

/**
 * Performs the MS2 deconvolution step. Since this occurs before clustering, we need to keep
 *    the spectra in separate files.
 * @param ip list of user parameters
 * @param inputSpecs set of input SpecSets
 * @param outputSpecs set of output SpecSets
 * @return true if successful, false if not
 */
bool performMS2Deconv(ParameterList & ip,
                      vector<string>& inputSpecs,
                      vector<string>& outputSpecs)
{
  ParameterList deconvParams;
  deconvParams.addIfExists(ip, "MAX_KLDiv");
  deconvParams.addIfExists(ip, "OUTPUT_SPECTRA_PATH");

  // Need the input iso envelope. It should be in the execution directory
  string isoModelPath = getPath(ip.getValue("EXE_DIR"), "resources", false);
  isoModelPath = getPath(isoModelPath, DEFAULT_DECONV_MODEL, false);

  outputSpecs.resize(inputSpecs.size());
  for (unsigned int i = 0; i < inputSpecs.size(); i++)
  {
    FilenameManager mngr(inputSpecs[i]);
    string outFile = getPath(mngr.path, mngr.filename, false);
    outFile += "_z.pklbin";
    outputSpecs[i] = outFile;
    deconvParams.setValue("OUTPUT_SPECTRA", outputSpecs[i]);
    deconvParams.setValue("INPUT_SPECTRA", inputSpecs[i]);
    deconvParams.setValue("INPUT_ISO_ENV", isoModelPath);
    IsoEnvelope isoModel;

    SpecSet inputSpecs;
    SpecSet outputSpecs;
    ExecDeconvoluteMS2 module(deconvParams,
                              &inputSpecs,
                              &isoModel,
                              &outputSpecs);

    if (!module.loadInputData())
    {
      return false;
    }

    if (!module.invoke())
    {
      return false;
    }

    if (!module.saveOutputData())
    {
      return false;
    }

    DEBUG_VAR(outputSpecs.size());
  }

  string deconvFilesPath = getProjPath(ip, DEFAULT_DECONV_FILES_LIST);

  if (!writeFileIndex(deconvFilesPath, outputSpecs))
  {
    return false;
  }

  return true;
}

//-----------------------------------------------------------------------------
bool performMsCluster(ParameterList & ip,
                      vector<string>& inputFilesList,
                      SpecSet &clustSpectra)
{
  DEBUG_TRACE;

  ParameterList msclusterParams;
  msclusterParams.addIfExists(ip, "CLUSTER_MIN_SIZE");
  if (ip.getValue("CLUSTER_TOOL", CLUSTER_MSCLUST) != CLUSTER_MSCLUST)
  {
    msclusterParams.setValue("CLUSTER_MIN_SIZE", "0");
  }
  msclusterParams.addIfExists(ip, "EXE_DIR");
  msclusterParams.addIfExists(ip, "TOLERANCE_PEAK");
  msclusterParams.addIfExists(ip, "TOLERANCE_PM");
  msclusterParams.addIfExists(ip, "TOLERANCE_PM_PPM");
  msclusterParams.addIfExists(ip, "CLUSTER_MODEL");
  msclusterParams.addIfExists(ip, "MIN_SPECTRUM_QUALITY");
  msclusterParams.addIfExists(ip, "RESET_MSCLUSTER_PMS");
  msclusterParams.addIfExists(ip, "MSCLUSTER_MIX_PROB");
  msclusterParams.addIfExists(ip, "CLUST_RANK_FILTER");

  msclusterParams.setValue("PROJECT_DIR", getProjPath(ip, SPECTRA_DIR));

  string outputPaths = getProjPath(ip, SPECS_MS_PATH);
  outputPaths += ";";
  outputPaths += getProjPath(ip, SPECS_MS_MGF_PATH);
  msclusterParams.setValue("OUTPUT_SPECTRA", outputPaths);

  msclusterParams.setValue("OUTPUT_CLUSTERS",
      getProjPath(ip, SPECS_MS_CLUST_PATH));

  DEBUG_TRACE;
  stringstream aux;
  msclusterParams.print(aux);
  DEBUG_MSG(aux.str());
  ClusterSet clusters;

  ExecMsCluster moduleMsCluster(msclusterParams,
      &inputFilesList,
      &clusters,
      &clustSpectra);

  if (!moduleMsCluster.invoke())
  {
    return false;
  }

  if (!moduleMsCluster.saveOutputData())
  {
    return false;
  }

  DEBUG_VAR(clustSpectra.size());

  return true;
}

//-----------------------------------------------------------------------------
bool performScoring(ParameterList & ip,
                    SpecSet & inputSpectra,
                    SpecSet & outputSpectra)
{
  DEBUG_TRACE;

  ParameterList scoringParams;
  scoringParams.addIfExists(ip, "EXE_DIR");
  scoringParams.addIfExists(ip, "OUTPUT_SPECTRA_PATH");
  scoringParams.addIfExists(ip, "PEPNOVO_EXE_DIR");
  scoringParams.addIfExists(ip, "PEPNOVO_MODEL_DIR");
  scoringParams.addIfExists(ip, "TOLERANCE_PEAK");
  scoringParams.addIfExists(ip, "TOLERANCE_PM");
  scoringParams.addIfExists(ip, "TOLERANCE_PM_PPM");
  scoringParams.addIfExists(ip, "CLUSTER_MODEL");
  scoringParams.addIfExists(ip, "MIN_SPECTRUM_QUALITY");
  scoringParams.addIfExists(ip, "GUESS_CHARGE");
  scoringParams.addIfExists(ip, "CORRECT_PM");
  scoringParams.addIfExists(ip, "PEPNOVO_MODEL");
  scoringParams.addIfExists(ip, "INSTRUMENT_TYPE");
  scoringParams.addIfExists(ip, "MIN_SILAC_COSINE");

  scoringParams.setValue("INPUT_CLUSTERS",
      getProjPath(ip, SPECS_MS_CLUST_PATH));

  if (ip.getValueInt("MERGE_SAME_PREC", 0) > 0
      || (ip.getValueInt("CLUSTER_MIN_SIZE", 0) >= 1
          && ip.getValue("CLUSTER_TOOL", "") == CLUSTER_PRMS))
  {
    scoringParams.addIfExists(ip, "MERGE_SAME_PREC");
    scoringParams.addIfExists(ip, "CLUSTER_MIN_SIZE");
    scoringParams.setValue("OUTPUT_CLUSTERS",
        getProjPath(ip, SPECS_SCORED_CLUST_PATH));
  }

  scoringParams.addIfExists(ip, "PRM_RANK_FILTER");
  scoringParams.addIfExists(ip, "PEPNOVO_MODEL");
  scoringParams.addIfExists(ip, "SKIP_PEPNOVO_INVOKE");
  scoringParams.addIfExists(ip, "PEPNOVO_PTMS");
  scoringParams.addIfExists(ip, "PEPNOVO_OUTDIR");
  scoringParams.addIfExists(ip, "NUM_CONSECUTIVE");
  scoringParams.addIfExists(ip, "PRM_CLUSTER_RATIO");
  scoringParams.addIfExists(ip, "BOOST_SILAC_PRMS");
  scoringParams.addIfExists(ip, "SILAC_SCAN_RANGE");
  scoringParams.addIfExists(ip, "FILTER_NONSILAC_PRMS");

  if (scoringParams.getValueInt("BOOST_SILAC_PRMS", 0) > 0)
  {
    scoringParams.setValue("OUTPUT_SILAC_PAIRS", getProjPath(ip,
            ALIGNS_SILAC_PATH));
  }

  scoringParams.addIfExists(ip, "PRM_SPEC_IDS");
  scoringParams.addIfExists(ip, "PRM_SPEC_ID_FORMAT");

  scoringParams.addIfExists(ip, "PROJECT_DIR");
  /*scoringParams.setValue("INPUT_SPECTRA",
   getProjPath(ip, "spectra/specs_ms.pklbin"));
   */
  // Free memory for PepNovo
  inputSpectra.clear();

  string ms2SpecPath = getProjPath(ip, SPECS_MS_PATH);
  scoringParams.setValue("PEPNOVO_OUTDIR", getProjPath(ip, SPECTRA_DIR));
  scoringParams.setValue("INPUT_SPECTRA", ms2SpecPath);
  scoringParams.setValue("OUTPUT_SPECTRA", getProjPath(ip, SPECS_SCORED_PATH));
  scoringParams.setValue("PEPNOVO_INPUT_MGF",
      getProjPath(ip, "spectra/specs_ms_pepnovo.mgf"));
  scoringParams.setValue("PEPNOVO_OUTPUT_PRMS",
      getProjPath(ip, "spectra/specs_scored.prms"));
  scoringParams.setValue("ENFORCE_DA_TOL", "1");

  //  specProtAlignParams.writeToFile(getProjPath(ip, "debug_specprotalign.params"));
  scoringParams.writeToFile(getProjPath(ip, "debug_scoring.params"));

  ClusterSet outputClusters;

  ExecPrmScoring moduleScoring(scoringParams, &outputSpectra, &outputClusters);

  string errorString;
  bool isValid = moduleScoring.validateParams(errorString);
  TEST_VALID;

  DEBUG_TRACE;

  //moduleScoring.loadInputData();

  if (!moduleScoring.loadInputData())
  {
    return false;
  }

  DEBUG_TRACE;

  bool returnStatus = moduleScoring.invoke();

  // Test for return status
  TEST_RETURN_STATUS("moduleScoring");

  DEBUG_TRACE;

  if (!moduleScoring.saveOutputData())
  {
    return false;
  }

  if (!inputSpectra.loadPklBin(ms2SpecPath.c_str()))
  {
    ERROR_MSG("Failed to re-load spectra from \'" << ms2SpecPath << "\'");
    return false;
  }

  return true;
}

bool performStatProtSeqs(ParameterList & ip)
{
  if (!ip.exists("FASTA_DATABASE") || !ip.exists("INPUT_SPEC_IDS"))
  {
    return false;
  }

  ParameterList statsParams;
  statsParams.addIfExists(ip, "TOLERANCE_PEAK");
  statsParams.setValue("INPUT_CONTIGS",
                       getProjPath(ip, "assembly/sps_seqs.pklbin"));
  statsParams.setValue("INPUT_STARS", getProjPath(ip, "spectra/stars.pklbin"));
  statsParams.setValue("INPUT_CONTIG_ABINFO",
                       getProjPath(ip, "assembly/component_info.bin"));
  statsParams.setValue("INPUT_MIDX",
                       getProjPath(ip, "homology/contigs_midx.pklbin"));
  statsParams.setValue("INPUT_MP", getProjPath(ip, "homology/contigs_mp.bin"));
  statsParams.setValue("INPUT_REF_INDICES",
                       getProjPath(ip, "spectra/contigs_indices.bin"));
  statsParams.setValue("INPUT_FASTA", ip.getValue("FASTA_DATABASE"));
  string aux = ip.getValue("EXE_DIR");
  aux += "/dancik_model.txt";
  statsParams.setValue("INPUT_ION_TYPES", aux);
  statsParams.addIfExists(ip, "INPUT_SPEC_IDS");
  statsParams.addIfExists(ip, "INPUT_CLUST_SPEC_IDS");
  statsParams.setValue("SPEC_ID_FORMAT",
                       ip.getValue("SPEC_ID_FORMAT", "MSGFDB"));

  statsParams.setValue("INPUT_FILE_INDEX",
                       getProjPath(ip, DEFAULT_USER_INPUT_FILES_LIST));
  statsParams.setValue("SPS_PROJECT_DIR", "");

  struct stat buf;
  if (stat(getProjPath(ip, "spectra/clusterData.bin").c_str(), &buf) != -1)
  {
    statsParams.setValue("INPUT_CLUSTERS_DIR", getProjPath(ip, ""));
  }
  statsParams.setValue("INPUT_FILE_MAPPING",
                       getProjPath(ip, DEFAULT_INPUT_MAPPING));

  statsParams.setValue("INPUT_SCAN_REF_FILES",
                       getProjPath(ip, DEFAULT_BIN_FILES_FILENAME));

  if (ip.exists("STATS_MIN_CONTIG_AA_TAG"))
  {
    statsParams.setValue("MIN_CONTIG_AA_TAG",
                         ip.getValue("STATS_MIN_CONTIG_AA_TAG"));
  }

  if (ip.exists("STATS_MIN_CONTIG_DB_MP"))
  {
    statsParams.setValue("MIN_CONTIG_DB_MP",
                         ip.getValue("STATS_MIN_CONTIG_DB_MP"));
  }

  if (ip.exists("STATS_ENDS_CHOP"))
  {
    statsParams.setValue("ENDS_CHOP", ip.getValue("STATS_ENDS_CHOP"));
  }

  if (ip.exists("STATS_TARGET_PROTEINS"))
  {
    statsParams.setValue("TARGET_PROTEINS",
                         ip.getValue("STATS_TARGET_PROTEINS"));
  }
  else
  {
    statsParams.setValue("TARGET_PROTEINS", "");
  }

  statsParams.setValue("TABLE_DELIM", "\t");
  statsParams.setValue("OUTPUT_SPS_STATS_FILE",
                       getProjPath(ip, "ReportData/tableStatsCummulative.tsv"));
  statsParams.setValue("OUTPUT_CONTIG_STATS_FILE",
                       getProjPath(ip, "ReportData/tableStatsContig.tsv"));

  statsParams.writeToFile(getProjPath(ip, "debug_statProtSeqs.params"));

  ExecReportSPSStats statsModule(statsParams);
  if (!statsModule.loadInputData())
  {
    return false;
  }
  if (!statsModule.invoke())
  {
    return false;
  }
  if (!statsModule.saveOutputData())
  {
    return false;
  }

  return true;
}
//-----------------------------------------------------------------------------
bool performGenoMS(ParameterList & ip)
{
  DEBUG_TRACE;

  ParameterList genoMSParams;

  genoMSParams.addIfExists(ip, "EXE_DIR");
  genoMSParams.addIfExists(ip, "OUTPUT_SPECTRA_PATH");

  genoMSParams.addIfExists(ip, "TOLERANCE_PEAK");
  genoMSParams.addIfExists(ip, "TOLERANCE_PM");
  genoMSParams.addIfExists(ip, "PEAK_PENALTY");

  genoMSParams.addIfExists(ip, "ALL2ALL_SIMILARITY");
  genoMSParams.addIfExists(ip, "HMM_LATE_ADD");
  genoMSParams.addIfExists(ip, "FDR_CUTOFF");
  genoMSParams.addIfExists(ip, "MUTATION_MODE");
  genoMSParams.addIfExists(ip, "DBCOMBINED");
  genoMSParams.addIfExists(ip, "DBROOTNAME");
  genoMSParams.addIfExists(ip, "GENOMESEQ");
  genoMSParams.addIfExists(ip, "DIGEST");
  genoMSParams.addIfExists(ip, "TEMPLATECONSTRAINTFILE");
  genoMSParams.addIfExists(ip, "FIXEDMOD");
  genoMSParams.addIfExists(ip, "PROJECT_DIR");
  genoMSParams.addIfExists(ip, "RUN_DBSEARCH");

  //new params for CID, ETD, HCD, or pairs/triplets
  genoMSParams.addIfExists(ip, "NUM_CONSECUTIVE");
  genoMSParams.addIfExists(ip, "INSTRUMENT_TYPE");
  genoMSParams.addIfExists(ip, "CLUSTER_MIN_SIZE");
  genoMSParams.addIfExists(ip, "CLUSTER_TOOL");
  genoMSParams.addIfExists(ip, "MERGE_SAME_PREC");

  //genoMSParams.setValue("PRMS",ip.getValue("PEPNOVO_OUTPUT_PRMS","spectra/specs_scored.prms"));
  //genoMSParams.setValue("PRMS", "spectra/pepnovo_out_CID.prms");
  genoMSParams.setValue("PRMS", SPECS_SCORED_MGF_PATH);
  //genoMSParams.setValue("SPECTRA",     "spectra/specs_ms_pepnovo.mgf");
  //genoMSParams.setValue("SPECTRA", "spectra/pepnovo_in_CID.mgf");
  genoMSParams.setValue("SPECTRA", SPECS_MS_MGF_PATH);
  genoMSParams.setValue("OUTPUT_FILE", "genoMS.out");
  genoMSParams.setValue("GENERATE_REPORTS", "1");
  genoMSParams.setValue("LOG_FILE", "genoMS.log");

  DEBUG_TRACE;
  stringstream aux;
  genoMSParams.print(aux);
  DEBUG_MSG(aux.str());
  DEBUG_TRACE;
  ExecGenoMS moduleGenoMS(genoMSParams);

  string errorString;
  bool isValid = moduleGenoMS.validateParams(errorString);
  TEST_VALID;

  DEBUG_TRACE;

  bool returnStatus = moduleGenoMS.invoke();
  // Test for return status
  TEST_RETURN_STATUS("moduleGenoMS");

  DEBUG_TRACE;

  return true;
}

bool performMergeOfCSPSAndGenoMS(ParameterList & ip)
{
  ///////////////////////////////////////////////
  //First we duplicate the input spectrum file
  ///////////////////////////////////////////////
  bool debug = false;
  DEBUG_TRACE;

  unsigned int cspsContigCount;
  unsigned int genomsContigCount;

  spsReports::ClusterData clusterData;

  SpecSet inputSpectra;
  SpecSet other;
  unsigned int specCount; // total (unique spectra)

  if (inputSpectra.loadPklBin("./spectra/specs_ms.pklbin") <= 0)
  {
    ERROR_MSG("Problem loading spectra from ./spectra/specs_ms.pklbin");
    return false;
  }

  DEBUG_MSG("First we have " << inputSpectra.size() << " inputSpectra!");
  if (other.loadPklBin("./spectra/specs_ms.pklbin") <= 0)
  {
    ERROR_MSG("Problem loading spectra from ./spectra/specs_ms.pklbin");
    return false;
  }

  inputSpectra.appendSpecSet(other);

  specCount = inputSpectra.size();

  DEBUG_MSG("Now we have " << specCount << " inputSpectra!");
  if (inputSpectra.savePklBin("./spectra/specs_ms.pklbin") <= 0)
  {
    ERROR_MSG("Problem saving doubled spectra to ./spectra/specs_ms.pklbin");
    return false;
  }
  if (debug)
  {
    inputSpectra.SaveSpecSet_mgf("./spectra/specs_ms.mgf");
  }

  ///////////////////////////////////////////////
  //Also duplicate the clusterData.bin
  ///////////////////////////////////////////////
  if ((ip.exists("CLUSTER_TOOL")
      && ip.getValue("CLUSTER_TOOL", "") == CLUSTER_PRMS)
      || (ip.exists("CLUSTER_MIN_SIZE")
          && ip.getValueInt("CLUSTER_MIN_SIZE") <= 0))
  {

    DEBUG_MSG("Doubling input_mapping.bin");
    vector<vector<int> > map;
    string fName = getProjPath(ip, DEFAULT_INPUT_MAPPING);
    Load_binArray(fName.c_str(), map);
    DEBUG_MSG("Loaded double");
    Save_binArrayDouble(fName.c_str(), map);
    DEBUG_MSG("Saved double");

  }
  else
  {
    DEBUG_MSG("Doubling clusterData.bin");
    clusterData.loadDataDouble(".");
    DEBUG_MSG("Loaded double");
    clusterData.saveData("..");
    DEBUG_MSG("Saved double");
  }
  ///////////////////////////////////////////////
  //Merge the star spectra
  ///////////////////////////////////////////////
  DEBUG_MSG("Merging star spectra, csps + genoms");
  if (inputSpectra.loadPklBin("./spectra/csps.stars.pklbin") <= 0)
  {
    ERROR_MSG("Problem loading spectra from ./spectra/csps.stars.pklbin");
    return false;
  }

  if (other.loadPklBin("./spectra/genoms.stars.pklbin") <= 0)
  {
    ERROR_MSG("Problem loading spectra from ./spectra/genoms.stars.pklbin");
    return false;
  }

  inputSpectra.appendSpecSet(other);
  if (inputSpectra.savePklBin("./spectra/stars.pklbin") <= 0)
  {
    ERROR_MSG("Problem saving doubled star spectra to ./spectra/stars.pklbin");
    return false;
  }

  DEBUG_MSG("Made " << inputSpectra.size() << " star spectra!");
  if (debug)
  {
    inputSpectra.SaveSpecSet_mgf("./spectra/stars.mgf");
  }
  ///////////////////////////////////////////////
  //Merge the contig spectra
  //////////////////////////////////////////////
  DEBUG_MSG("Merging contig spectra");
  if (inputSpectra.loadPklBin("./assembly/csps.sps_seqs.pklbin") <= 0)
  {
    ERROR_MSG("Problem loading spectra from ./assembly/csps.sps_seqs.pklbin");
    return false;
  }

  cspsContigCount = inputSpectra.size();

  if (other.loadPklBin("./assembly/genoms.sps_seqs.pklbin") <= 0)
  {
    ERROR_MSG("Problem loading spectra from ./assembly/genoms.sps_seqs.pklbin");
    return false;
  }
  genomsContigCount = other.size();
  inputSpectra.appendSpecSet(other);
  if (inputSpectra.savePklBin("./assembly/sps_seqs.pklbin") <= 0)
  {
    ERROR_MSG("Problem saving merged contig spectra to ./assembly/sps_seqs.pklbin");
    return false;
  }
  DEBUG_MSG("Made " << inputSpectra.size() << " contig spectra!");
  DEBUG_MSG("CSPS has contigs 0-" << (cspsContigCount - 1));
  DEBUG_MSG("GenoMS has contigs " << cspsContigCount << "-" << (genomsContigCount + cspsContigCount - 1));

  if (debug)
  {
    inputSpectra.SaveSpecSet_mgf("./assembly/sps_seqs.mgf");
  }
  ///////////////////////////////////////////////
  //Merge abruin infos
  ///////////////////////////////////////////////
  DEBUG_MSG("Merging abruijn info");
  abinfo_t info1;
  abinfo_t info2;
  string abinfoOutFile = "./assembly/component_info.bin";

  Load_abinfo("./assembly/csps.component_info.bin", info1);
  Load_abinfo("./assembly/genoms.component_info.bin", info2);

  if (!merge_abinfo(info1, info2, abinfoOutFile.c_str(), specCount / 2))
  {
    ERROR_MSG("Error merging abruijn info");
    return false;
  }
  if (debug)
  {
    abinfo_t info;
    Load_abinfo("./assembly/component_info.bin", info);
    dumpAbInfo("./assembly/component_info.txt", info);
  }
  ///////////////////////////////////////////////
  //merge contig_mp.bin
  ///////////////////////////////////////////////
  DEBUG_MSG("Merging contig_mp_all.bin");
  vector<vector<int> > *m_contigs_mp_sub1;
  vector<vector<int> > *m_contigs_mp_sub2;

  m_contigs_mp_sub1 = new vector<vector<int> >();
  m_contigs_mp_sub2 = new vector<vector<int> >();

  DEBUG_MSG("Initialized vectors..");
  // load the data
  if (Load_binArray<int, vector>("./homology/csps.contigs_mp.bin",
                                 *m_contigs_mp_sub1) < 0)
  {

    return false;
  }

  if (Load_binArray<int, vector>("./homology/genoms.contigs_mp.bin",
                                 *m_contigs_mp_sub2) < 0)
  {

    return false;
  }
  DEBUG_MSG("Loaded contig_mp_all vectors");
  if (Save_doubleBinArray<int>("./homology/contigs_mp.bin",
                               *m_contigs_mp_sub1,
                               *m_contigs_mp_sub2) < 0)
  {
    ERROR_MSG("Error merging arrays for ./homology/contigs_mp.bin");
    return false;

  }

  DEBUG_MSG("CSPS contig_mp.bin = " << (*m_contigs_mp_sub1).size() << "x" << (*m_contigs_mp_sub1)[0].size());
  DEBUG_MSG("GenoMS contig_mp.bin = " << (*m_contigs_mp_sub2).size() << "x" << (*m_contigs_mp_sub2)[0].size());
  if (debug)
  {

    cout << "DEBUG: homology/contigs_mp.bin" << endl;

    Load_binArray<int, vector>("./homology/contigs_mp.bin", *m_contigs_mp_sub1);
    for (int i = 0; i < (*m_contigs_mp_sub1).size(); ++i)
    {
      for (int j = 0; j < (*m_contigs_mp_sub1)[0].size(); ++j)
      {
        cout << " " << (*m_contigs_mp_sub1)[i][j];
      }
      cout << endl;
    }

  }
  ///////////////////////////////////////////////
  //Merge contigs_midx.pklbin
  ///////////////////////////////////////////////
  DEBUG_MSG("Merging contigs_midx.pklbin");
  if (inputSpectra.loadPklBin("./homology/csps.contigs_midx.pklbin") <= 0)
  {
    ERROR_MSG("Problem loading spectra from ./homology/csps.contigs_midx.pklbin");
    return false;
  }

  if (other.loadPklBin("./homology/genoms.contigs_midx.pklbin") <= 0)
  {
    ERROR_MSG("Problem loading spectra from ./homology/genoms.contigs_midx.pklbin");
    return false;
  }

  inputSpectra.appendSpecSet(other);
  if (inputSpectra.savePklBin("./homology/contigs_midx.pklbin") <= 0)
  {
    ERROR_MSG("Problem saving contig index spectra to ./homology/contigs_midx.pklbin");
    return false;
  }
  if (debug)
  {
    inputSpectra.SaveSpecSet_mgf("./homology/contigs_midx.mgf");
  }

  ///////////////////////////////////////////////
  //merge contig_mp_all.bin
  ///////////////////////////////////////////////
  DEBUG_MSG("Merging contig_mp_all.bin");
  vector<vector<int> > *m_contigs_mp1;
  vector<vector<int> > *m_contigs_mp2;

  m_contigs_mp1 = new vector<vector<int> >();
  m_contigs_mp2 = new vector<vector<int> >();

  DEBUG_MSG("Initialized vectors..");
  // load the data
  if (Load_binArray<int, vector>("./homology/csps.contigs_mp_all.bin",
                                 *m_contigs_mp1) < 0)
  {

    return false;
  }

  if (Load_binArray<int, vector>("./homology/genoms.contigs_mp_all.bin",
                                 *m_contigs_mp2) < 0)
  {

    return false;
  }
  DEBUG_MSG("Loaded contig_mp_all vectors");
  if (Save_doubleBinArray<int>("./homology/contigs_mp_all.bin",
                               *m_contigs_mp1,
                               *m_contigs_mp2) < 0)
  {
    ERROR_MSG("Error merging arrays for ./homology/contigs_mp_all.bin");
    return false;

  }

  DEBUG_MSG("CSPS contig_mp_all.bin = " << (*m_contigs_mp1).size() << "x" << (*m_contigs_mp1)[0].size());
  DEBUG_MSG("GenoMS contig_mp_all.bin = " << (*m_contigs_mp2).size() << "x" << (*m_contigs_mp2)[0].size());
  if (debug)
  {

    cout << "DEBUG: homology/contigs_mp_all.bin" << endl;

    Load_binArray<int, vector>("./homology/contigs_mp_all.bin", *m_contigs_mp1);
    for (int i = 0; i < (*m_contigs_mp1).size(); ++i)
    {
      for (int j = 0; j < (*m_contigs_mp1)[0].size(); ++j)
      {
        cout << " " << (*m_contigs_mp1)[i][j];
      }
      cout << endl;
    }

  }
  ///////////////////////////////////////////////
  //Merge contigs_midx_all.pklbin
  ///////////////////////////////////////////////
  /*DEBUG_MSG("Merging contigs_midx_all.pklbin");
   if (inputSpectra.loadPklBin("./homology/csps.contigs_midx_all.pklbin") <= 0)
   {
   ERROR_MSG("Problem loading spectra from ./homology/csps.contigs_midx_all.pklbin");
   return false;
   }

   if (other.loadPklBin("./homology/genoms.contigs_midx_all.pklbin") <= 0)
   {
   ERROR_MSG("Problem loading spectra from ./homology/genoms.contigs_midx_all.pklbin");
   return false;
   }

   inputSpectra.appendSpecSet(other);
   if (inputSpectra.savePklBin("./homology/contigs_midx_all.pklbin") <= 0)
   {
   ERROR_MSG("Problem saving contig index spectra to ./homology/contigs_midx_all.pklbin");
   return false;
   }
   if (debug)
   {
   inputSpectra.SaveSpecSet_mgf("./homology/contigs_midx_all.mgf");
   }
   */
  ///////////////////////////////////////////////
  //Merge contigs.pklbin
  ///////////////////////////////////////////////
  DEBUG_MSG("Merging contigs.pklbin");
  if (inputSpectra.loadPklBin("./spectra/csps.contigs.pklbin") <= 0)
  {
    ERROR_MSG("Problem loading spectra from ./spectra/csps.contigs.pklbin");
    return false;
  }
  unsigned int cspsMappedContigCount = inputSpectra.size();
  if (other.loadPklBin("./spectra/genoms.contigs.pklbin") <= 0)
  {
    ERROR_MSG("Problem loading spectra from ./spectra/genoms.contigs.pklbin");
    return false;
  }

  inputSpectra.appendSpecSet(other);
  if (inputSpectra.savePklBin("./spectra/contigs.pklbin") <= 0)
  {
    ERROR_MSG("Problem saving contig spectra to ./spectra/contigs.pklbin");
    return false;
  }
  DEBUG_MSG("Made " << inputSpectra.size() << " mapped contig spectra!");
  DEBUG_MSG("CSPS has contigs 0-" << (cspsMappedContigCount - 1));
  DEBUG_MSG("GenoMS has contigs " << cspsMappedContigCount << "-" << (inputSpectra.size() - 1));
  if (debug)
  {
    inputSpectra.SaveSpecSet_mgf("./spectra/contigs.mgf");
  }
  ///////////////////////////////////////////////
  //merge contig_indices.bin
  ///////////////////////////////////////////////
  DEBUG_MSG("Merging contig_indices.bin");

  vector<vector<int> > *m_contigIndex1;
  vector<vector<int> > *m_contigIndex2;

  m_contigIndex1 = new vector<vector<int> >();
  m_contigIndex2 = new vector<vector<int> >();

  // load the data
  if (Load_binArray<int, vector>("./spectra/csps.contigs_indices.bin",
                                 *m_contigIndex1) < 0)
  {

    ERROR_MSG("Error loading binary array ./spectra/csps.contigs_indices.bin");
    return false;
  }

  if (Load_binArray<int, vector>("./spectra/genoms.contigs_indices.bin",
                                 *m_contigIndex2) < 0)
  {

    ERROR_MSG("Error loading binary array ./spectra/genoms.contigs_indices.bin");
    return false;
  }
  DEBUG_MSG("CSPS contig_indices.bin = " << (*m_contigIndex1).size() << "x" << (*m_contigIndex1)[0].size());
  DEBUG_MSG("GenoMS contig_indices.bin = " << (*m_contigIndex2).size() << "x" << (*m_contigIndex2)[0].size());

  //Update component indices from genoms result
  for (int i = 0; i < m_contigIndex2->size(); ++i)
  {
    (*m_contigIndex2)[i][0] += cspsContigCount;

  }

  if (Save_doubleBinArray<int>("./spectra/contigs_indices.bin",
                               *m_contigIndex1,
                               *m_contigIndex2) < 0)
  {
    ERROR_MSG("Error merging arrays for ./spectra/contigs_indices.bin");
    return false;

  }

  if (debug)
  {

    cout << "DEBUG: spectra/contigs_indices.bin" << endl;

    Load_binArray<int, vector>("./spectra/contigs_indices.bin",
                               *m_contigIndex1);
    for (int i = 0; i < (*m_contigIndex1).size(); ++i)
    {
      for (int j = 0; j < (*m_contigIndex1)[0].size(); ++j)
      {
        cout << " " << (*m_contigIndex1)[i][j];
      }
      cout << endl;
    }

  }

  ///////////////////////////////////////////////
  //merge homglue_ref_mp.bin
  ///////////////////////////////////////////////
  DEBUG_MSG("Merging homglue_ref_mp.bin");
  vector<vector<int> > *m_homglueRefMp1;
  vector<vector<int> > *m_homglueRefMp2;

  m_homglueRefMp1 = new vector<vector<int> >();
  m_homglueRefMp2 = new vector<vector<int> >();

  // load the data
  if (Load_binArray<int, vector>("./homology/csps.homglue_ref_mp.bin",
                                 *m_homglueRefMp1) < 0)
  {

    return false;
  }

  if (Load_binArray<int, vector>("./homology/genoms.homglue_ref_mp.bin",
                                 *m_homglueRefMp2) < 0)
  {

    return false;
  }

  DEBUG_MSG("CSPS homglue_ref_mp.bin = " << (*m_homglueRefMp1).size() << "x" << (*m_homglueRefMp1)[0].size());
  DEBUG_MSG("GenoMS homglue_ref_mp.bin = " << (*m_homglueRefMp2).size() << "x" << (*m_homglueRefMp2)[0].size());

  if (Save_doubleBinArray("./homology/homglue_ref_mp.bin",
                          *m_homglueRefMp1,
                          *m_homglueRefMp2) < 0)
  {
    ERROR_MSG("Error merging arrays for ./homology/homglue_ref_mp.bin");
    return false;

  }

  if (debug)
  {

    cout << "DEBUG: homology/homglue_ref_mp.bin" << endl;

    Load_binArray<int, vector>("./homology/homglue_ref_mp.bin",
                               *m_homglueRefMp1);
    for (int i = 0; i < (*m_homglueRefMp1).size(); ++i)
    {
      for (int j = 0; j < (*m_homglueRefMp1)[0].size(); ++j)
      {
        cout << " " << (*m_homglueRefMp1)[i][j];
      }
      cout << endl;
    }

  }

  ///////////////////////////////////////////////
  //Merge homglue_ref_midx.pklbin
  ///////////////////////////////////////////////
  DEBUG_MSG("Merging homglue_ref_midx.pklbin");
  if (inputSpectra.loadPklBin("./homology/csps.homglue_ref_midx.pklbin") <= 0)
  {
    ERROR_MSG("Problem loading spectra from ./homology/csps.homglue_ref_midx.pklbin");
    return false;
  }

  if (other.loadPklBin("./homology/genoms.homglue_ref_midx.pklbin") <= 0)
  {
    ERROR_MSG("Problem loading spectra from ./homology/genoms.homglue_ref_midx.pklbin");
    return false;
  }

  inputSpectra.appendSpecSet(other);
  if (inputSpectra.savePklBin("./homology/homglue_ref_midx.pklbin") <= 0)
  {
    ERROR_MSG("Problem saving contig index spectra to ./homology/homglue_ref_midx.pklbin");
    return false;
  }

  if (debug)
  {
    DEBUG_MSG("homglue_ref_midx.pklbin has" << inputSpectra.size() << " spectra");
    inputSpectra.SaveSpecSet_mgf("./homology/homglue_ref_midx.mgf");
  }

  ///////////////////////////////////////////////
  //Merge homglue_matches.pklbin
  ///////////////////////////////////////////////
  DEBUG_MSG("Merging homglue_matches.pklbin");
  if (inputSpectra.loadPklBin("./homology/csps.homglue_matches.pklbin") <= 0)
  {
    ERROR_MSG("Problem loading spectra from ./homology/csps.homglue_matches.pklbin");
    return false;
  }

  if (other.loadPklBin("./homology/genoms.homglue_matches.pklbin") <= 0)
  {
    ERROR_MSG("Problem loading spectra from ./homology/genoms.homglue_matches.pklbin");
    return false;
  }

  inputSpectra.appendSpecSet(other);
  if (inputSpectra.savePklBin("./homology/homglue_matches.pklbin") <= 0)
  {
    ERROR_MSG("Problem saving contig index spectra to ./homology/homglue_matches.pklbin");
    return false;
  }
  if (debug)
  {
    DEBUG_MSG("homglue_matches.pklbin has" << inputSpectra.size() << " spectra");
    inputSpectra.SaveSpecSet_mgf("./homology/homglue_matches.mgf");
  }

  ///////////////////////////////////////////////
  //merge homglue_matches_mp.bin
  ///////////////////////////////////////////////
  DEBUG_MSG("Merging homglue_matches_mp.bin");
  vector<vector<int> > *m_homglueMatchMp1;
  vector<vector<int> > *m_homglueMatchMp2;

  m_homglueMatchMp1 = new vector<vector<int> >();
  m_homglueMatchMp2 = new vector<vector<int> >();

  // load the data
  if (Load_binArray<int, vector>("./homology/csps.homglue_matches_mp.bin",
                                 *m_homglueMatchMp1) < 0)
  {

    return false;
  }

  if (Load_binArray<int, vector>("./homology/genoms.homglue_matches_mp.bin",
                                 *m_homglueMatchMp2) < 0)
  {

    return false;
  }

  DEBUG_MSG("CSPS homglue_matches_mp.bin = " << (*m_homglueMatchMp1).size() << "x" << (*m_homglueMatchMp1)[0].size());
  DEBUG_MSG("GenoMS homglue_matches_mp.bin = " << (*m_homglueMatchMp2).size() << "x" << (*m_homglueMatchMp2)[0].size());

  if (Save_doubleBinArray("./homology/homglue_matches_mp.bin",
                          *m_homglueMatchMp1,
                          *m_homglueMatchMp2) < 0)
  {
    ERROR_MSG("Error merging arrays for ./homology/homglue_matches_mp.bin");
    return false;

  }

  ///////////////////////////////////////////////
  //Merge homglue_matches_midx.pklbin
  ///////////////////////////////////////////////
  DEBUG_MSG("Merging homglue_matches.pklbin");
  if (inputSpectra.loadPklBin("./homology/csps.homglue_matches_midx.pklbin")
      <= 0)
  {
    ERROR_MSG("Problem loading spectra from ./homology/csps.homglue_matches_midx.pklbin");
    return false;
  }

  if (other.loadPklBin("./homology/genoms.homglue_matches_midx.pklbin") <= 0)
  {
    ERROR_MSG("Problem loading spectra from ./homology/genoms.homglue_matches_midx.pklbin");
    return false;
  }

  inputSpectra.appendSpecSet(other);
  if (inputSpectra.savePklBin("./homology/homglue_matches_midx.pklbin") <= 0)
  {
    ERROR_MSG("Problem saving contig index spectra to ./homology/homglue_matches_midx.pklbin");
    return false;
  }

  if (debug)
  {
    DEBUG_MSG("homglue_matches_midx.pklbin has" << inputSpectra.size() << " spectra");
    inputSpectra.SaveSpecSet_mgf("./homology/homglue_matches_midx.mgf");
  }

  ///////////////////////////////////////////////
  //merge ref_sps_names.txt
  ///////////////////////////////////////////////
  DEBUG_MSG("Merging ref_sps_names.txt");
  vector<string> genoMScontigNames;
  int i, j;
  int contigNum = cspsContigCount + 1;
  //The GenoMS contig numbers need to be updated
  readFilesFromFile("./homology/genoms.ref_sps_names.txt", genoMScontigNames);

  FILE * fp = fopen("./homology/genoms.ref_sps_names.txt", "wb");
  //We assume the names are GenoMS:i, where i is the contig number
  for (i = 0; i < genoMScontigNames.size(); ++i)
  {
    DEBUG_MSG("Orig genoMS contig Name: " << genoMScontigNames[i]);
    string newName = "";
    for (j = 0; j < genoMScontigNames[i].size(); ++j)
    {
      if (genoMScontigNames[i][j] == ':')
        break;
      else
        newName += genoMScontigNames[i][j];
    }
    newName += ":";
    newName += intToString(contigNum);
    newName += "\n";
    genoMScontigNames[i] = newName;
    DEBUG_MSG("New genoMS contig Name: " << genoMScontigNames[i]);

    fwrite((void *)(genoMScontigNames[i].c_str()),
           genoMScontigNames[i].size(),
           1,
           fp);
    contigNum += 1;

  }

  fclose(fp);

  //Only concatenate if the csps version exists.
  ifstream tmpFile("./homology/csps.ref_sps_names.txt");
  if (tmpFile)
  {
    if (!concatenateFiles("./homology/csps.ref_sps_names.txt",
                          "./homology/genoms.ref_sps_names.txt",
                          "./homology/ref_sps_names.txt"))
      return false;
  }

  return true;

}

bool generateRelaunchScript(ParameterList & ip)
{
  string relauncherScriptName = "relauncher.sh";
  FILE* myfile;
  myfile = fopen(relauncherScriptName.c_str(), "rb");
  if (myfile)
  {
    fclose(myfile);
    return true;
  }

  string dbFileName = getProjPath(ip, "relaunch.protid.fasta");
  string paramsFile = "relaunch.params";
  string exeDir = ip.getValue("EXE_DIR", ".");

  string logLevel = ip.getValue("LOG_LEVEL", "0");
  if (strlen(logLevel.c_str()) == 0)
    logLevel = "0";
  //Create the command

  //Change to the project directory
  string relauncherScript = "cd " + getProjPath(ip, ".") + "\n";

  relauncherScript += "echo \"Running\" > status.txt\n";

  //If we were currently running in merge mode, then we need to move the old CSPS files back to their old places
  // if(ip.getValueInt("MERGE_FLAG",0) == 1)
  myfile = fopen("assembly/csps.sps_seqs.pklbin", "rb");
  if (myfile)
  {
    fclose(myfile);
    relauncherScript +=
        "mv assembly/csps.sps_seqs.pklbin assembly/sps_seqs.pklbin\n";
    relauncherScript +=
        "mv assembly/csps.component_info.bin assembly/component_info.bin\n";
    //relauncherScript +=
    //    "mv homology/csps.contigs_midx_all.pklbin homology/contigs_midx_all.pklbin\n";
    relauncherScript +=
        "mv homology/csps.contigs_midx.pklbin homology/contigs_midx.pklbin\n";
    relauncherScript +=
        "mv homology/csps.contigs_mp_all.bin homology/contigs_mp_all.bin\n";
    relauncherScript +=
        "mv homology/csps.contigs_mp.bin homology/contigs_mp.bin\n";
    relauncherScript +=
        "mv homology/csps.homglue_matches_midx.pklbin homology/homglue_matches_midx.pklbin\n";
    relauncherScript +=
        "mv homology/csps.homglue_matches_mp.bin homology/homglue_matches_mp.bin\n";
    relauncherScript +=
        "mv homology/csps.homglue_matches.pklbin homology/homglue_matches.pklbin\n";
    relauncherScript +=
        "mv homology/csps.homglue_ref_midx.pklbin homology/homglue_ref_midx.pklbin\n";
    relauncherScript +=
        "mv homology/csps.homglue_ref_mp.bin homology/homglue_ref_mp.bin\n";

    ifstream ifile("homology/csps.ref_sps_names.txt");
    if (ifile)
    {
      relauncherScript +=
          "mv homology/csps.ref_sps_names.txt homology/ref_sps_names.txt\n";
    }

    relauncherScript +=
        "mv spectra/csps.contigs_indices.bin spectra/contigs_indices.bin\n";
    relauncherScript +=
        "mv spectra/csps.contigs.pklbin spectra/contigs.pklbin\n";
    relauncherScript += "mv spectra/csps.stars.pklbin spectra/stars.pklbin\n";
  }

  //execute the command
  relauncherScript += exeDir + "/main_specnets " + paramsFile;
  relauncherScript += " -i tagsearch";
  relauncherScript += " -ll " + logLevel;
  relauncherScript += " -lf "
      + ip.getValue("LOG_FILE_NAME", "relauncher_log.txt");
  int val = ip.getValueInt("GRID_EXECUTION");
  if (val != 0)
    relauncherScript += " -g";
  //val = ip.getValueInt("GENOMS_FLAG");
  //if(val != 0)
  //  relauncherScript += " -q";
  //val = ip.getValueInt("MERGE_FLAG");
  //if(val != 0)
  //  relauncherScript += " -m";
  relauncherScript += "\n";
  //Write the relauncher.sh
  FILE * f = fopen(relauncherScriptName.c_str(), "wb");

  if (f == NULL)
  {
    ERROR_MSG("Unable to open relauncher file for writing");
    return false;
  }

  val = fwrite(relauncherScript.c_str(),
               sizeof(char),
               strlen(relauncherScript.c_str()),
               f);
  if (val != strlen(relauncherScript.c_str()))
  {
    ERROR_MSG("Problem encountered writing relauncher file");
    ERROR_MSG(val);
    ERROR_MSG(strlen(relauncherScript.c_str()));
    return false;
  }
  fclose(f);

  //Write the updated params file
  ip.setValue("FASTA_DATABASE", dbFileName);
  ip.setValue("GENOMS_FLAG", "0");
  ip.setValue("MERGE_FLAG", "0");
  ip.writeToFile(getProjPath(ip, paramsFile));

  DEBUG_MSG("Wrote relauncher.sh");
  return true;
}

//-----------------------------------------------------------------------------
bool performFilterPairs(ParameterList & ip,
                        SpecSet & inputSpectra,
                        SpecSet & inputSpectraMS2,
                        SpectrumPairSet & filteredPairs,
                        vector<TwoValues<float> > & ratios,
                        vector<TwoValues<float> > & means,
                        vector<float> & varTerms,
                        list<vector<float> > & alignStats,
                        vector<vector<float> > & specStats,
                        std::vector<unsigned int> & idxKept,
                        std::vector<TwoValues<float> > & pvalues,
                        bool gridExecutionFlag,
                        bool resume)
{
  ParameterList filterpairsParams;
  filterpairsParams.addIfExists(ip, "AMINO_ACID_MASSES");
  filterpairsParams.addIfExists(ip, "TOLERANCE_PEAK");
  filterpairsParams.addIfExists(ip, "TOLERANCE_PM");
  filterpairsParams.addIfExists(ip, "PARTIAL_OVERLAPS");
  filterpairsParams.addIfExists(ip, "AA_DIFF_COUNT");
  filterpairsParams.addIfExists(ip, "MIN_OVERLAP_AREA");
  filterpairsParams.addIfExists(ip, "MIN_RATIO");
  filterpairsParams.addIfExists(ip, "MIN_SHIFT");
  filterpairsParams.addIfExists(ip, "MAX_SHIFT");
  filterpairsParams.addIfExists(ip, "PAIRS_MATCH_MODE");
  filterpairsParams.addIfExists(ip, "PAIRS_MIN_COSINE");

  filterpairsParams.addIfExists(ip, "GRID_TYPE");
  filterpairsParams.addIfExists(ip, "GRID_NUMNODES");
  filterpairsParams.addIfExists(ip, "GRID_NUMCPUS");
  filterpairsParams.addIfExists(ip, "GRID_EXE_DIR");
  filterpairsParams.addIfExists(ip, "GRID_SGE_EXE_DIR");
  filterpairsParams.addIfExists(ip, "GRID_PARAMS");

  filterpairsParams.addIfExists(ip, "MAX_PVALUE");
  filterpairsParams.addIfExists(ip, "RESOLUTION");
  filterpairsParams.addIfExists(ip, "TAGS_MATCH_FLANK");
  filterpairsParams.addIfExists(ip, "TAGS_MATCH_COUNT");

  filterpairsParams.addIfExists(ip, "PROJECT_DIR");
  filterpairsParams.addIfExists(ip, "OUTPUT_SPECTRA_PATH");

  filterpairsParams.addIfDoesntExist("USE_MIN_DIST_57", "1");
  filterpairsParams.addIfDoesntExist("SPEC_TYPE_MSMS", "0");

  filterpairsParams.setValue("MIN_MATCHED_PEAKS",
                             ip.getValue("MIN_MATCHED_PEAKS"));
  filterpairsParams.setValue("MAX_SHIFT", ip.getValue("MAX_MOD_MASS"));

  filterpairsParams.setValue("INPUT_SPECS_PKLBIN",
                             getProjPath(ip, SPECS_SCORED_PATH));
  filterpairsParams.setValue("GRID_DATA_DIR", getProjPath(ip, ALIGNS_DIR));
  filterpairsParams.setValue("OUTPUT_ALIGNS",
                             getProjPath(ip, "./aligns/pairs_raw.bin"));
  filterpairsParams.setValue("OUTPUT_MEANS",
                             getProjPath(ip, "./aligns/means.bin"));
  filterpairsParams.setValue("OUTPUT_VARIANCE",
                             getProjPath(ip, "./aligns/vars.bin"));
  filterpairsParams.setValue("OUTPUT_RATIOS",
                             getProjPath(ip, "./aligns/ratios.bin"));

  filterpairsParams.writeToFile("debug_filterpairs.params");

  DEBUG_TRACE;

  stringstream aux;
  filterpairsParams.print(aux);
  DEBUG_MSG(aux.str());

  DEBUG_TRACE;

  if (filterpairsParams.getValue("PAIRS_MATCH_MODE", "") == "cosine")
  {
    filterpairsParams.setValue("OUTPUT_ALIGNS",
                               getProjPath(ip, "./aligns/pairs_cosine.bin"));
    for (unsigned int i = 0; i < inputSpectraMS2.size(); i++)
    {
      for (unsigned int j = 0; j < inputSpectraMS2[i].size(); j++)
        inputSpectraMS2[i][j][1] = sqrt(inputSpectraMS2[i][j][1]);
      inputSpectraMS2[i].normalize2();
    }
  }

  ExecFilterPairs moduleFilterPairs(filterpairsParams,
                                    &inputSpectra,
                                    &inputSpectraMS2,
                                    &inputSpectraMS2,
                                    &filteredPairs,
                                    &ratios,
                                    &means,
                                    &varTerms,
                                    &alignStats,
                                    &specStats,
                                    &idxKept,
                                    &pvalues);

  string errorString;
  bool isValid = moduleFilterPairs.validateParams(errorString);
  TEST_VALID;

  DEBUG_TRACE;

  bool returnStatus;
  if (!ip.exists("GRID_NUMNODES") or ip.getValueInt("GRID_NUMNODES") <= 0)
  {
    returnStatus = moduleFilterPairs.invoke();
  }
  else
  {
    DEBUG_TRACE;
    int numNodes = ip.getValueInt("GRID_NUMNODES");

    string gridType = ip.getValue("GRID_TYPE");
    if (gridType == "pbs")
    {
      ParallelPbsExecution exec(&moduleFilterPairs,
          gridExecutionFlag,
          !gridExecutionFlag,
          resume);
      returnStatus = exec.invoke(numNodes);
    }
    else if (gridType == "sge")
    {
      ParallelSgeExecution exec(&moduleFilterPairs,
          gridExecutionFlag,
          !gridExecutionFlag,
          resume);
      returnStatus = exec.invoke(numNodes);
    }
    else if (gridType == "threaded")
    {
      ParallelThreadedExecution exec(&moduleFilterPairs);
      returnStatus = exec.invoke(numNodes);
    }
  }

  // Test for return status
  TEST_RETURN_STATUS("moduleFilterPairs");

  DEBUG_TRACE;

  moduleFilterPairs.saveOutputData();

  return true;
}

//-----------------------------------------------------------------------------
bool performFilterAligns(ParameterList & ip,
                         SpectrumPairSet & filteredPairs,
                         vector<TwoValues<float> > & ratios,
                         vector<TwoValues<float> > & means,
                         vector<float> & varTerms,
                         std::vector<unsigned int> & idxKept,
                         std::vector<TwoValues<float> > & pvalues)
{
  ParameterList filteralignsParams;
  filteralignsParams.addIfExists(ip, "TOLERANCE_PM");
  filteralignsParams.addIfExists(ip, "MIN_RATIO");
  filteralignsParams.addIfExists(ip, "MAX_PVALUE");
  filteralignsParams.addIfExists(ip, "FILTER_TRIGS");

  filteralignsParams.setValue("INPUT_ALIGNS",
                              getProjPath(ip, "./aligns/pairs_raw.bin"));
  filteralignsParams.setValue("INPUT_MEANS",
                              getProjPath(ip, "./aligns/means.bin"));
  filteralignsParams.setValue("INPUT_VARIANCE",
                              getProjPath(ip, "./aligns/vars.bin"));
  filteralignsParams.setValue("INPUT_RATIOS",
                              getProjPath(ip, "./aligns/ratios.bin"));
  filteralignsParams.setValue("OUTPUT_ALIGNS",
                              getProjPath(ip, "./aligns/pairs.bin"));

  filteralignsParams.addIfExists(ip, "OUTPUT_SPECTRA_PATH");

  filteralignsParams.writeToFile("debug_filteraligns.params");

  DEBUG_TRACE;

  stringstream aux;
  filteralignsParams.print(aux);
  DEBUG_MSG(aux.str());

  DEBUG_TRACE;

  ExecFilterAligns moduleFilterAligns(filteralignsParams,
                                      &filteredPairs,
                                      &ratios,
                                      &means,
                                      &varTerms,
                                      &idxKept,
                                      &pvalues);

  string errorString;
  bool isValid = moduleFilterAligns.validateParams(errorString);
  TEST_VALID;

  DEBUG_TRACE;

  bool returnStatus = moduleFilterAligns.invoke();

  // Test for return status
  TEST_RETURN_STATUS("moduleFilterAligns");

  DEBUG_TRACE;

  // Save output data
  moduleFilterAligns.saveOutputData();

  return true;
}

//-----------------------------------------------------------------------------
bool performAlignment(ParameterList & ip,
                      AAJumps & jumps,
                      SpecSet &inputSpectra,
                      SpectrumPairSet &inputPairs,
                      SpecSet &pairAlignments,
                      SpecSet &starSpectraOnly,
                      SpecSet &starSpectra,
                      vector<unsigned int> &alignedSpectra)
{
  ParameterList alignParams;
  alignParams.addIfExists(ip, "TOLERANCE_PEAK");
  alignParams.addIfExists(ip, "TOLERANCE_PM");
  alignParams.addIfExists(ip, "RESOLUTION");
  alignParams.addIfExists(ip, "PARTIAL_OVERLAPS");
  alignParams.addIfExists(ip, "MAX_AA_JUMP");
  alignParams.addIfExists(ip, "PENALTY_PTM");
  alignParams.addIfExists(ip, "PENALTY_SAME_VERTEX");

  alignParams.addIfExists(ip, "PROJECT_DIR");
  alignParams.addIfExists(ip, "OUTPUT_SPECTRA_PATH");

  alignParams.setValue("PENALTY_PTM_PEAKS", "-1.0");
  alignParams.setValue("OUTPUT_SPECS",
                       getProjPath(ip, "./spectra/aligns_specs.pklbin"));
  alignParams.setValue("OUTPUT_STARS",
                       getProjPath(ip, "./spectra/stars_only.pklbin"));
  alignParams.setValue("OUTPUT_STARS_INDEX",
                       getProjPath(ip, "./spectra/stars_indices.bin"));
  alignParams.setValue("OUTPUT_STARS_ALL",
                       getProjPath(ip, "./spectra/stars.pklbin"));
  if (ip.exists("AMINO_ACID_MASSES"))
  {
    alignParams.setValue("AMINO_ACID_MASSES",
                         "../" + ip.getValue("AMINO_ACID_MASSES"));
  }

  alignParams.writeToFile("debug_align.params");

  DEBUG_TRACE;
  ExecAlignment moduleAlignment(alignParams,
                                &inputSpectra,
                                &inputPairs,
                                &pairAlignments,
                                &starSpectraOnly,
                                &starSpectra,
                                &alignedSpectra);

  string errorString;
  bool isValid = moduleAlignment.validateParams(errorString);
  TEST_VALID;

  DEBUG_TRACE;

#if 1
  bool returnStatus = moduleAlignment.invoke();
#else
  DEBUG_TRACE;
  ParallelThreadedExecution exec(&moduleAlignment);
  bool returnStatus = exec.invoke(1, 1);
#endif

  // Test for return status
  TEST_RETURN_STATUS("ExecAlignment");

  DEBUG_TRACE;
  // Save output data
  moduleAlignment.saveOutputData();

  return true;
}

//-----------------------------------------------------------------------------
bool performFilterStarPairs(ParameterList & ip,
                            SpectrumPairSet &inputPairs,
                            SpecSet &starSpectra,
                            vector<vector<float> > &ratios,
                            SpecSet &matchedPeaks)
{
  ParameterList filterstarpairsParams;
  filterstarpairsParams.addIfExists(ip, "TOLERANCE_PEAK");
  filterstarpairsParams.addIfExists(ip, "TOLERANCE_PM");
  filterstarpairsParams.addIfExists(ip, "MAX_MOD_MASS");
  filterstarpairsParams.addIfExists(ip, "PARTIAL_OVERLAPS");
  filterstarpairsParams.addIfExists(ip, "MIN_RATIO");
  filterstarpairsParams.addIfExists(ip, "MIN_MATCHED_PEAKS");
  filterstarpairsParams.addIfExists(ip, "MAX_AA_JUMP");
  filterstarpairsParams.addIfExists(ip, "PENALTY_PTM");
  filterstarpairsParams.addIfExists(ip, "PENALTY_SAME_VERTEX");

  filterstarpairsParams.addIfExists(ip, "PROJECT_DIR");
  filterstarpairsParams.addIfExists(ip, "OUTPUT_SPECTRA_PATH");

  filterstarpairsParams.setValue("SPEC_TYPE_MSMS", "0");
  filterstarpairsParams.setValue("OUTPUT_ALIGNS",
                                 getProjPath(ip, "aligns/pairs_stars.bin"));

  filterstarpairsParams.writeToFile(getProjPath(ip,
                                                "debug_filterstarpairs.params"));

  DEBUG_TRACE;
  ExecFilterStarPairs moduleFilterStarPairs(filterstarpairsParams,
                                            &inputPairs,
                                            &starSpectra,
                                            &ratios,
                                            &matchedPeaks);
  DEBUG_TRACE;

  string errorString;
  bool isValid = moduleFilterStarPairs.validateParams(errorString);
  TEST_VALID;

#if 1
  bool returnStatus = moduleFilterStarPairs.invoke();
#else
  ParallelSgeExecution exec(&moduleFilterStarPairs);
  //ParallelThreadedExecution exec(&moduleFilterStarPairs);
  bool returnStatus = exec.invoke(1, 1);
#endif

  // Test for return status
  TEST_RETURN_STATUS("ExecFilterStarPairs");

  // Save output data
  returnStatus = moduleFilterStarPairs.saveOutputData();
  //  TEST_SAVE_OUPUT_DATA("ExecFilterStarPairs");

  return true;
}

//-----------------------------------------------------------------------------

bool performAssembly(ParameterList & ip,
                     SpecSet & starSpectra,
                     SpectrumPairSet & starPairs,
                     Clusters & contigShifts,
                     abinfo_t & contigAbinfo)
{
  ParameterList assemblyParams;
  assemblyParams.addIfExists(ip, "TOLERANCE_PEAK");
  assemblyParams.addIfExists(ip, "TOLERANCE_PM");
  assemblyParams.addIfExists(ip, "MAX_MOD_MASS");
  assemblyParams.addIfExists(ip, "MAX_AA_JUMP");
  assemblyParams.addIfExists(ip, "PENALTY_PTM");
  assemblyParams.addIfExists(ip, "MIN_MATCHED_PEAKS");
  assemblyParams.addIfExists(ip, "PENALTY_SAME_VERTEX");
  assemblyParams.addIfExists(ip, "PENALTY_PTM");
  assemblyParams.addIfExists(ip, "PARALLEL_PATHS");
  assemblyParams.addIfExists(ip, "ADD_ENDPOINTS");

  assemblyParams.addIfExists(ip, "PROJECT_DIR");
  assemblyParams.addIfExists(ip, "OUTPUT_SPECTRA_PATH");

  assemblyParams.setValue("NO_SEQUENCING", "0");
  assemblyParams.setValue("EDGE_SCORE_TYPE", "1");
  assemblyParams.setValue("SPEC_TYPE_MSMS", "0");
  assemblyParams.setValue("GRAPH_TYPE", "2");
  assemblyParams.setValue("OUTPUT_COMPLETE_ABRUIJN", "");
  assemblyParams.setValue("PATH_MIN_PEAKS",
                          ip.getValue("SPSPATH_MIN_NUM_PEAKS"));
  assemblyParams.setValue("PATH_MIN_SPECS",
                          ip.getValue("SPSPATH_MIN_NUM_SPECS"));
  assemblyParams.setValue("MIN_EDGES_TO_COMPONENT",
                          ip.getValue("SPS_MIN_EDGES_TO_COMPONENT"));

  assemblyParams.setValue("OUTPUT_CLUSTERS",
                          getProjPath(ip,
                                      "assembly/path_spectra_as_cluster.txt"));
  assemblyParams.setValue("OUTPUT_SPECS",
                          getProjPath(ip, "assembly/sps_seqs.pklbin"));
  assemblyParams.setValue("OUTPUT_ABINFO",
                          getProjPath(ip, "assembly/component_info.bin"));

  assemblyParams.writeToFile(getProjPath(ip, "debug_assembly.params"));

  DEBUG_TRACE;
  ExecAssembly moduleAssembly(assemblyParams,
                              &starSpectra,
                              &starPairs,
                              &contigShifts,
                              &contigAbinfo);
  DEBUG_TRACE;

  string errorString;
  bool isValid = moduleAssembly.validateParams(errorString);
  TEST_VALID;

  bool returnStatus = moduleAssembly.invoke();
  // Test for return status
  TEST_RETURN_STATUS("ExecAssembly");

  returnStatus = moduleAssembly.saveOutputData();
  TEST_SAVE_OUPUT_DATA("ExecAssembly");

  return true;
}

//-----------------------------------------------------------------------------

bool performContigAlignment(ParameterList & ip,
                            Clusters & contigs,
                            SpectrumPairSet & contigPairs)
{
  ParameterList alignParams;
  alignParams.addIfExists(ip, "TOLERANCE_PEAK");
  alignParams.addIfExists(ip, "TOLERANCE_PM");
  alignParams.addIfExists(ip, "MIN_MATCHED_PEAKS");
  alignParams.addIfExists(ip, "MIN_METACONTIG_SCORE");

  alignParams.setValue("INPUT_SPECTRA",
                       getProjPath(ip, "assembly/sps_seqs.pklbin"));

  alignParams.setValue("OUTPUT_CONTIG_ALIGNS",
                       getProjPath(ip, "aligns/sps_seqs_aligns.bin"));

  alignParams.addIfExists(ip, "OUTPUT_SPECTRA_PATH");

  DEBUG_TRACE;
  ExecFilterContigPairs moduleAlign(alignParams,
                                    &contigs.consensus,
                                    &contigPairs);
  DEBUG_TRACE;

  string errorString;
  bool isValid = moduleAlign.validateParams(errorString);
  TEST_VALID;

  bool returnStatus = moduleAlign.invoke();
  // Test for return status
  TEST_RETURN_STATUS("ExecFilterContigPairs");

  returnStatus = moduleAlign.saveOutputData();
  TEST_SAVE_OUPUT_DATA("ExecFilterContigPairs");

  return true;
}

//-----------------------------------------------------------------------------

bool performMetaAssembly(ParameterList & ip,
                         Clusters & contigs,
                         SpectrumPairSet & contigAligns,
                         abinfo_t & contigAbruijn,
                         SpecSet & starSpectra,
                         Clusters & outputMetaContigs,
                         abinfo_t & metaContigAbruijn)
{
  ParameterList assemblyParams;
  assemblyParams.addIfExists(ip, "TOLERANCE_PEAK");
  assemblyParams.addIfExists(ip, "TOLERANCE_PM");
  assemblyParams.addIfExists(ip, "MIN_MATCHED_PEAKS");
  assemblyParams.addIfExists(ip, "MIN_METACONTIG_SCORE");
  assemblyParams.addIfExists(ip, "MIN_METACONTIG_SIZE");

  assemblyParams.addIfExists(ip, "OUTPUT_SPECTRA_PATH");

  assemblyParams.setValue("OUTPUT_CLUSTERS",
                          getProjPath(ip,
                                      "assembly/meta_path_spectra_as_cluster.txt"));
  assemblyParams.setValue("OUTPUT_SPECTRA",
                          getProjPath(ip, "assembly/sps_seqs.pklbin"));
  assemblyParams.setValue("OUTPUT_ABINFO",
                          getProjPath(ip, "assembly/component_info.bin"));

  assemblyParams.writeToFile(getProjPath(ip, "debug_metaAssembly.params"));

  DEBUG_TRACE;
  ExecMetaAssembly moduleAssembly(assemblyParams,
                                  &contigs.consensus,
                                  &contigAligns,
                                  &contigAbruijn,
                                  &starSpectra,
                                  &outputMetaContigs,
                                  &metaContigAbruijn);
  DEBUG_TRACE;

  string errorString;
  bool isValid = moduleAssembly.validateParams(errorString);
  TEST_VALID;

  bool returnStatus = moduleAssembly.invoke();
  // Test for return status
  TEST_RETURN_STATUS("ExecMetaAssembly");

  returnStatus = moduleAssembly.saveOutputData();
  TEST_SAVE_OUPUT_DATA("ExecMetaAssembly");

  return true;
}

//-----------------------------------------------------------------------------
bool performTagsearch(ParameterList & ip,
                      SpecSet & spectra,
                      DB_fasta &db,
                      vector<unsigned int> *specsToSearch,
                      PeptideSpectrumMatchSet & psmSet,
                      bool decoy)
{
  ParameterList tagsearchParams;
  tagsearchParams.addIfExists(ip, "RESOLUTION");
  tagsearchParams.addIfExists(ip, "MAX_MOD_MASS");
  tagsearchParams.addIfExists(ip, "TOLERANCE_PEAK");
  tagsearchParams.addIfExists(ip, "TOLERANCE_PM");

  tagsearchParams.addIfExists(ip, "PROJECT_DIR");
  tagsearchParams.addIfExists(ip, "OUTPUT_SPECTRA_PATH");

  tagsearchParams.addIfExists(ip, "MAX_PARSIMONY");
  tagsearchParams.addIfExists(ip, "TAG_LEN");
  tagsearchParams.addIfDoesntExist("TAG_LEN", "6");
  tagsearchParams.addIfExists(ip, "DOUBLE_AA_JUMPS");
  tagsearchParams.addIfDoesntExist("DOUBLE_AA_JUMPS", "1");
  tagsearchParams.addIfExists(ip, "MAX_NUM_TAGS");
  tagsearchParams.addIfDoesntExist("MAX_NUM_TAGS", "0");
  tagsearchParams.addIfExists(ip, "MATCH_TAG_FLANKING_MASSES");
  tagsearchParams.addIfDoesntExist("MATCH_TAG_FLANKING_MASSES", "0");
  tagsearchParams.addIfExists(ip, "TAG_MATCH_TOP_SCORING_ONLY");
  tagsearchParams.addIfExists(ip, "TAG_PEAK_SKIP_PENALTY");

  tagsearchParams.addIfExists(ip, "CORRECT_PEPTIDE_SPECS");
  tagsearchParams.addIfExists(ip, "CORRECT_PEPTIDE_PSMS");

  if (decoy)
  {
    DEBUG_TRACE;
    tagsearchParams.setValue("OUTPUT_PSM",
        getProjPath(ip, "assembly/tagsearchpsm_decoy.txt"));
    tagsearchParams.setValue("OUTPUT_PEPTIDES",
        getProjPath(ip, "spectra/tagsearch_decoy.txt"));
    //tagsearchParams.setValue("OUTPUT_SPECS",
    //                         getProjPath(ip, "assembly/meta_sps_seqs_decoy.pklbin"));
    tagsearchParams.writeToFile("debug_tagsearch_decoy.params");
  }
  else
  {
    tagsearchParams.setValue("OUTPUT_PSM",
        getProjPath(ip, "assembly/tagsearchpsm.txt"));
    tagsearchParams.setValue("OUTPUT_PEPTIDES",
        getProjPath(ip, "spectra/tagsearch.txt"));
    //tagsearchParams.setValue("OUTPUT_SPECS",
    //                         getProjPath(ip, "assembly/meta_sps_seqs.pklbin"));
    tagsearchParams.writeToFile("debug_tagsearch.params");
  }

  DEBUG_TRACE;
  ExecTagSearch moduleTagsearch(tagsearchParams,
                                &spectra,
                                &db,
                                specsToSearch,
                                &psmSet);
  DEBUG_TRACE;

  string errorString;
  bool isValid = moduleTagsearch.validateParams(errorString);
  TEST_VALID;

  bool returnStatus = moduleTagsearch.invoke();
  // Test for return status
  TEST_RETURN_STATUS("ExecTagSearch");

  // Saving output data
  returnStatus = moduleTagsearch.saveOutputData();
  TEST_SAVE_OUPUT_DATA("ExecTagSearch");

  return true;
}

//-----------------------------------------------------------------------------

bool performContigProtAlign(ParameterList & ip,
                            SpecSet & contigSpectra,
                            DB_fasta & db,
                            DB_fasta & dbDecoy,
                            PeptideSpectrumMatchSet & psmTag,
                            PeptideSpectrumMatchSet & psmTagDecoy,
                            PenaltyMatrix & penaltyMatrixBlosum,
                            PenaltyMatrix & penaltyMatrixMods,
                            SpecSet & matchedSpectraAll,
                            SpecSet & matchedSpectra,
                            PeptideSpectrumMatchSet & psmSet,
                            PeptideSpectrumMatchSet & psmSetDecoy,
                            PeptideSpectrumMatchSet & psmSetFdr,
                            bool gridExecutionFlag,
                            bool resume)
{
  ParameterList specProtAlignParams;
  specProtAlignParams.addIfExists(ip, "ALIGNMENT_RESOLUTION");
  specProtAlignParams.addIfExists(ip, "TOLERANCE_PEAK");
  specProtAlignParams.addIfExists(ip, "TOLERANCE_PM");

  specProtAlignParams.setValue("GRID_EXECUTION_FLAG",
                               gridExecutionFlag ? "1" : "0");
  specProtAlignParams.setValue("GRID_RESUME_FLAG", resume ? "1" : "0");
  specProtAlignParams.addIfExists(ip, "GRID_TYPE");

  //If we have an argument called "GRID_NUMNODES_SPECPROTALIGN", then replace GRID NUMNODES with that
  if (ip.exists("GRID_NUMNODES_SPECPROTALIGN"))
  {
    string val = ip.getValue("GRID_NUMNODES_SPECPROTALIGN");
    specProtAlignParams.setValue("GRID_NUMNODES", val);
  }
  else
    specProtAlignParams.addIfExists(ip, "GRID_NUMNODES");

  specProtAlignParams.addIfExists(ip, "GRID_NUMCPUS");
  specProtAlignParams.addIfExists(ip, "GRID_EXE_DIR");
  specProtAlignParams.addIfExists(ip, "GRID_SGE_EXE_DIR");
  specProtAlignParams.addIfExists(ip, "GRID_PARAMS");
  specProtAlignParams.setValue("GRID_DATA_DIR_TARGET",
                               getProjPath(ip, "./homology/grid_t"));
  specProtAlignParams.setValue("GRID_DATA_DIR_DECOY",
                               getProjPath(ip, "./homology/grid_d"));

  specProtAlignParams.addIfExists(ip, "PROJECT_DIR");
  specProtAlignParams.addIfExists(ip, "OUTPUT_SPECTRA_PATH");

  if (ip.exists("MIN_MATCHED_PEAKS_DB"))
    specProtAlignParams.setValue("MIN_MATCHED_PEAKS_DB",
                                 ip.getValue("MIN_MATCHED_PEAKS_DB"));
  else
    specProtAlignParams.addIfDoesntExist("MIN_MATCHED_PEAKS_DB", "7");
  specProtAlignParams.addIfExists(ip, "MAX_NUM_MODS");
  specProtAlignParams.addIfDoesntExist("MAX_NUM_MODS", "2");
  specProtAlignParams.addIfExists(ip, "MIN_MOD_MASS");
  specProtAlignParams.addIfDoesntExist("MIN_MOD_MASS", "-100");
  specProtAlignParams.addIfExists(ip, "MAX_MOD_MASS");
  specProtAlignParams.addIfDoesntExist("MAX_MOD_MASS", "100");
  specProtAlignParams.addIfExists(ip, "MAX_PARSIMONY");
  specProtAlignParams.addIfDoesntExist("MAX_PARSIMONY", "1");
  specProtAlignParams.addIfExists(ip, "ALIGNMENT_SCORE_THRESHOLD");
  specProtAlignParams.addIfDoesntExist("ALIGNMENT_SCORE_THRESHOLD", "0.90");
  specProtAlignParams.addIfExists(ip, "MAX_ALIGN_DB_GAP_AAS");
  specProtAlignParams.addIfDoesntExist("MAX_ALIGN_DB_GAP_AAS", "6");
  specProtAlignParams.addIfExists(ip, "MAX_ALIGN_SPECTRUM_GAP_DALTONS");
  specProtAlignParams.addIfDoesntExist("MAX_ALIGN_SPECTRUM_GAP_DALTONS",
                                       "1000");

  if (ip.exists("FDR_CONTIG_THRESHOLD"))
  {
    specProtAlignParams.setValue("FDR_THRESHOLD",
                                 ip.getValue("FDR_CONTIG_THRESHOLD"));
  }

  if (ip.exists("CONTIGPROTALIGN_START_IDX"))
  {
    specProtAlignParams.setValue("IDX_START",
                                 ip.getValue("CONTIGPROTALIGN_START_IDX"));
  }
  if (ip.exists("CONTIGPROTALIGN_END_IDX"))
  {
    specProtAlignParams.setValue("IDX_END",
                                 ip.getValue("CONTIGPROTALIGN_END_IDX"));
  }

  AAJumps jumps(1); // Amino acid masses
  jumps.saveJumps(getProjPath(ip, "homology/specprotalign_aa_masses.txt").c_str());
  specProtAlignParams.setValue("AMINO_ACID_MASSES",
                               getProjPath(ip,
                                           "homology/specprotalign_aa_masses.txt"));

  specProtAlignParams.setValue("BLOSUM_PENALTY_FILE",
                               getProjPath(ip,
                                           "homology/specprotalign_blosum_penal.txt"));
  specProtAlignParams.setValue("MODS_PENALTY_FILE",
                               getProjPath(ip,
                                           "homology/specprotalign_mod_penal.txt"));
  specProtAlignParams.setValue("KNOWN_MODS_FILE",
                               getProjPath(ip,
                                           "homology/specprotalign_mod_known.txt"));

  specProtAlignParams.setValue("OUTPUT_MATCHED_PEAKS_IDX_TGT",
                               getProjPath(ip,
                                           "homology/contigs_midx_tgt.pklbin"));
  specProtAlignParams.setValue("OUTPUT_MATCHED_SPECS_TGT",
                               getProjPath(ip, "homology/contigs_tgt.pklbin"));
  specProtAlignParams.setValue("OUTPUT_PSM_TGT",
                               getProjPath(ip, "homology/contigs_psm_tgt.txt"));
  specProtAlignParams.setValue("OUTPUT_MATCHED_PEAKS_IDX_DEC",
                               getProjPath(ip,
                                           "homology/contigs_midx_dec.pklbin"));
  specProtAlignParams.setValue("OUTPUT_MATCHED_SPECS_DEC",
                               getProjPath(ip, "homology/contigs_dec.pklbin"));
  specProtAlignParams.setValue("OUTPUT_PSM_DEC",
                               getProjPath(ip, "homology/contigs_psm_dec.txt"));

  specProtAlignParams.setValue("OUTPUT_MATCHED_PEAKS_IDX_ALL",
                               getProjPath(ip,
                                           "homology/contigs_midx_all.pklbin"));
  specProtAlignParams.setValue("OUTPUT_MATCHED_PEAKS_IDX",
                               getProjPath(ip, "homology/contigs_midx.pklbin"));
  specProtAlignParams.setValue("OUTPUT_MATCHED_SPECS",
                               getProjPath(ip, "spectra/contigs.pklbin"));
  specProtAlignParams.setValue("OUTPUT_SPECTRA_SPRINKLED",
                               getProjPath(ip,
                                           "homology/contigs_sprinkled.pklbin"));
  specProtAlignParams.setValue("OUTPUT_PSM",
                               getProjPath(ip, "homology/contigs_psm.txt"));

  specProtAlignParams.setValue("OUTPUT_MATCHED_PROTS",
                               getProjPath(ip, "homology/contigs_mp.bin"));
  specProtAlignParams.setValue("OUTPUT_MATCHED_PROTS_ALL",
                               getProjPath(ip, "homology/contigs_mp_all.bin"));
  specProtAlignParams.setValue("OUTPUT_MATCHED_SPECS_IDX",
                               getProjPath(ip, "spectra/contigs_indices.bin"));

  specProtAlignParams.setValue("OUTPUT_PSM_MOD_FREQS",
                               getProjPath(ip,
                                           "homology/contigs_psm_mod_freqs.txt"));
  specProtAlignParams.setValue("OUTPUT_PSM_FDR",
                               getProjPath(ip, "homology/contigs_psm_fdr.txt"));

  specProtAlignParams.setValue("ENFORCE_ENDPEAKS", "0"); // Do not enforce endpeaks for contigs

  specProtAlignParams.addIfExists(ip, "ALIGNMENT_SCORE_THRESHOLD");
  specProtAlignParams.addIfDoesntExist("ALIGNMENT_SCORE_THRESHOLD", "0.90");

  specProtAlignParams.addIfExists(ip, "PENALTY_ALIGNMENT");
  specProtAlignParams.addIfDoesntExist("PENALTY_ALIGNMENT", "0");
  specProtAlignParams.addIfExists(ip, "PENALTY_ALIGNMENT_ALPHA");
  specProtAlignParams.addIfDoesntExist("PENALTY_ALIGNMENT_ALPHA", "2.0");
  specProtAlignParams.addIfExists(ip, "PENALTY_ALIGNMENT_BETA");
  specProtAlignParams.addIfDoesntExist("PENALTY_ALIGNMENT_BETA", "2.0");
  specProtAlignParams.addIfExists(ip, "PENALTY_ALIGNMENT_UNKNOWN_PENALTY");
  specProtAlignParams.addIfDoesntExist("PENALTY_ALIGNMENT_UNKNOWN_PENALTY",
                                       "1.0");
  specProtAlignParams.addIfExists(ip, "PENALTY_ALIGNMENT_UNKNOWN_MULTIPLIER");
  specProtAlignParams.addIfDoesntExist("PENALTY_ALIGNMENT_UNKNOWN_MULTIPLIER",
                                       "2.0");

  specProtAlignParams.addIfExists(ip, "DEBUG_SPECPROTALIGN");
  specProtAlignParams.addIfExists(ip, "DEBUG_SPECPROTALIGN_SPECTRA");
  specProtAlignParams.addIfExists(ip, "DEBUG_SPECPROTALIGN_TIME");
  specProtAlignParams.addIfExists(ip, "DEBUG_SPECPROTALIGN_ANNO");
  specProtAlignParams.addIfExists(ip, "DEBUG_SPECPROTALIGN_SPRINKLE");

  specProtAlignParams.setValue("SPEC_TYPE_MSMS", "0");
  specProtAlignParams.setValue("MIN_RATIO", "0.4");

  specProtAlignParams.writeToFile(getProjPath(ip,
                                              "debug_contigprotalign.params"));

  DEBUG_VAR(contigSpectra.size());

  ExecSpecProtAlignTgtDecoy moduleSpecProtAlignTgtDecoy(specProtAlignParams,
                                                        &contigSpectra,
                                                        0x0, // No PRM spectra to "sprinkle" into contig spectra
                                                        &db,
                                                        &dbDecoy,
                                                        &psmTag,
                                                        &psmTagDecoy,
                                                        &penaltyMatrixBlosum,
                                                        &penaltyMatrixMods,
                                                        &matchedSpectraAll,
                                                        &matchedSpectra,
                                                        &psmSet,
                                                        &psmSetDecoy,
                                                        &psmSetFdr);
  DEBUG_TRACE;

  string errorString;
  bool isValid = moduleSpecProtAlignTgtDecoy.validateParams(errorString);
  TEST_VALID;

  bool returnStatus = moduleSpecProtAlignTgtDecoy.invoke();

  // Test for return status
  TEST_RETURN_STATUS("ExecSpecProtAlignTgtDecoy");

  DEBUG_VAR(matchedSpectra.size());

  returnStatus = moduleSpecProtAlignTgtDecoy.saveOutputData();
  DEBUG_VAR(returnStatus);
  TEST_SAVE_OUPUT_DATA("ExecSpecProtAlignTgtDecoy");

  return true;
}

//-----------------------------------------------------------------------------

bool performSpecProtAlign(ParameterList & ip,
                          SpecSet & contigSpectra,
                          SpecSet & prmSpectra,
                          DB_fasta & db,
                          DB_fasta & dbDecoy,
                          PeptideSpectrumMatchSet & psmTag,
                          PeptideSpectrumMatchSet & psmTagDecoy,
                          PenaltyMatrix & penaltyMatrixBlosum,
                          PenaltyMatrix & penaltyMatrixMods,
                          SpecSet & matchedSpectraAll,
                          SpecSet & matchedSpectra,
                          PeptideSpectrumMatchSet & psmSet,
                          PeptideSpectrumMatchSet & psmSetDecoy,
                          PeptideSpectrumMatchSet & psmSetFdr,
                          bool gridExecutionFlag,
                          bool resume)
{
  ParameterList specProtAlignParams;
  specProtAlignParams.addIfExists(ip, "ALIGNMENT_RESOLUTION");
  specProtAlignParams.addIfExists(ip, "MAX_MOD_MASS");
  specProtAlignParams.addIfExists(ip, "TOLERANCE_PEAK");
  specProtAlignParams.addIfExists(ip, "TOLERANCE_PM");

  // We need to reset endpoints for spectrum alignment
  specProtAlignParams.setValue("RESET_ENDPOINTS", "1");
  specProtAlignParams.addIfExists(ip, "SPEC_TYPE_MSMS");

  specProtAlignParams.setValue("GRID_EXECUTION_FLAG",
                               gridExecutionFlag ? "1" : "0");
  specProtAlignParams.setValue("GRID_RESUME_FLAG", resume ? "1" : "0");
  specProtAlignParams.addIfExists(ip, "GRID_TYPE");
  specProtAlignParams.addIfExists(ip, "GRID_NUMNODES");
  specProtAlignParams.addIfExists(ip, "GRID_NUMCPUS");
  specProtAlignParams.addIfExists(ip, "GRID_EXE_DIR");
  specProtAlignParams.addIfExists(ip, "GRID_SGE_EXE_DIR");
  specProtAlignParams.addIfExists(ip, "GRID_PARAMS");
  specProtAlignParams.setValue("GRID_DATA_DIR_TARGET",
                               getProjPath(ip, "./homology/grid2_t"));
  specProtAlignParams.setValue("GRID_DATA_DIR_DECOY",
                               getProjPath(ip, "./homology/grid2_d"));

  specProtAlignParams.addIfExists(ip, "PROJECT_DIR");
  specProtAlignParams.addIfExists(ip, "OUTPUT_SPECTRA_PATH");

  if (ip.exists("MIN_MATCHED_PEAKS_DB"))
    specProtAlignParams.setValue("MIN_MATCHED_PEAKS_DB",
                                 ip.getValue("MIN_MATCHED_PEAKS_DB"));
  else
    specProtAlignParams.addIfDoesntExist("MIN_MATCHED_PEAKS_DB", "7");
  specProtAlignParams.addIfExists(ip, "MAX_NUM_MODS");
  specProtAlignParams.addIfDoesntExist("MAX_NUM_MODS", "2");
  specProtAlignParams.addIfExists(ip, "MIN_MOD_MASS");
  specProtAlignParams.addIfDoesntExist("MIN_MOD_MASS", "-100");
  specProtAlignParams.addIfExists(ip, "MAX_MOD_MASS");
  specProtAlignParams.addIfDoesntExist("MAX_MOD_MASS", "100");
  specProtAlignParams.addIfExists(ip, "MAX_PARSIMONY");
  specProtAlignParams.addIfDoesntExist("MAX_PARSIMONY", "1");
  specProtAlignParams.addIfExists(ip, "MAX_ALIGN_DB_GAP_AAS");
  specProtAlignParams.addIfDoesntExist("MAX_ALIGN_DB_GAP_AAS", "6");
  specProtAlignParams.addIfExists(ip, "MAX_ALIGN_SPECTRUM_GAP_DALTONS");
  specProtAlignParams.addIfDoesntExist("MAX_ALIGN_SPECTRUM_GAP_DALTONS",
                                       "1000");

  specProtAlignParams.setValue("ENFORCE_ENDPEAKS", "1"); // Always enforce end peaks for spectrum matches

  specProtAlignParams.addIfExists(ip, "PENALTY_ALIGNMENT");
  specProtAlignParams.addIfDoesntExist("PENALTY_ALIGNMENT", "0");
  specProtAlignParams.addIfExists(ip, "PENALTY_ALIGNMENT_ALPHA");
  specProtAlignParams.addIfDoesntExist("PENALTY_ALIGNMENT_ALPHA", "2.0");
  specProtAlignParams.addIfExists(ip, "PENALTY_ALIGNMENT_BETA");
  specProtAlignParams.addIfDoesntExist("PENALTY_ALIGNMENT_BETA", "2.0");
  specProtAlignParams.addIfExists(ip, "PENALTY_ALIGNMENT_UNKNOWN_PENALTY");
  specProtAlignParams.addIfDoesntExist("PENALTY_ALIGNMENT_UNKNOWN_PENALTY",
                                       "1.0");
  specProtAlignParams.addIfExists(ip, "PENALTY_ALIGNMENT_UNKNOWN_MULTIPLIER");
  specProtAlignParams.addIfDoesntExist("PENALTY_ALIGNMENT_UNKNOWN_MULTIPLIER",
                                       "2.0");
  specProtAlignParams.addIfExists(ip,
                                  "SPECPROTALIGN_PRM_CONTRIBUTION_MULTIPLIER");
  specProtAlignParams.addIfExists(ip, "SPECPROTALIGN_ROUND_ANNOTATION_MAX");

  cout << "FDR_SPECTRUM_THRESHOLD = " << ip.getValue("FDR_SPECTRUM_THRESHOLD")
      << endl;
  if (ip.exists("FDR_SPECTRUM_THRESHOLD"))
  {
    specProtAlignParams.setValue("FDR_THRESHOLD",
                                 ip.getValue("FDR_SPECTRUM_THRESHOLD"));
  }

  if (ip.exists("SPECPROTALIGN_START_IDX"))
  {
    specProtAlignParams.setValue("IDX_START",
                                 ip.getValue("SPECPROTALIGN_START_IDX"));
  }
  if (ip.exists("SPECPROTALIGN_END_IDX"))
  {
    specProtAlignParams.setValue("IDX_END",
                                 ip.getValue("SPECPROTALIGN_END_IDX"));
  }

  AAJumps jumps(1); // Amino acid masses
  jumps.saveJumps(getProjPath(ip, "homology/specprotalign_aa_masses.txt").c_str());
  specProtAlignParams.setValue("AMINO_ACID_MASSES",
                               getProjPath(ip,
                                           "homology/specprotalign_aa_masses.txt"));

  specProtAlignParams.setValue("BLOSUM_PENALTY_FILE",
                               getProjPath(ip,
                                           "homology/specprotalign_blosum_penal.txt"));
  specProtAlignParams.setValue("MODS_PENALTY_FILE",
                               getProjPath(ip,
                                           "homology/specprotalign_mod_penal.txt"));
  specProtAlignParams.setValue("KNOWN_MODS_FILE",
                               getProjPath(ip,
                                           "homology/specprotalign_mod_known.txt"));

  specProtAlignParams.setValue("OUTPUT_MATCHED_PEAKS_IDX_TGT",
                               getProjPath(ip,
                                           "homology/spectrum_midx_tgt.pklbin"));
  specProtAlignParams.setValue("OUTPUT_MATCHED_SPECS_TGT",
                               getProjPath(ip, "homology/spectrum_tgt.pklbin"));
  specProtAlignParams.setValue("OUTPUT_PSM_TGT",
                               getProjPath(ip,
                                           "homology/spectrum_psm_tgt.txt"));
  specProtAlignParams.setValue("OUTPUT_MATCHED_PEAKS_IDX_DEC",
                               getProjPath(ip,
                                           "homology/spectrum_midx_dec.pklbin"));
  specProtAlignParams.setValue("OUTPUT_MATCHED_SPECS_DEC",
                               getProjPath(ip, "homology/spectrum_dec.pklbin"));
  specProtAlignParams.setValue("OUTPUT_PSM_DEC",
                               getProjPath(ip,
                                           "homology/spectrum_psm_dec.txt"));

  specProtAlignParams.setValue("OUTPUT_MATCHED_PEAKS_IDX",
                               getProjPath(ip,
                                           "homology/spectrum_midx.pklbin"));
  specProtAlignParams.setValue("OUTPUT_MATCHED_SPECS",
                               getProjPath(ip, "spectra/spectrum.pklbin"));
  specProtAlignParams.setValue("OUTPUT_SPECTRA_SPRINKLED",
                               getProjPath(ip,
                                           "homology/spectrum_sprinkled.pklbin"));
  specProtAlignParams.setValue("OUTPUT_PSM",
                               getProjPath(ip, "homology/spectrum_psm.txt"));

  specProtAlignParams.setValue("OUTPUT_MATCHED_PROTS",
                               getProjPath(ip, "homology/spectrum_mp.bin"));
  specProtAlignParams.setValue("OUTPUT_MATCHED_PROTS_ALL",
                               getProjPath(ip, "homology/spectrum_mp_all.bin"));
  specProtAlignParams.setValue("OUTPUT_MATCHED_SPECS_IDX",
                               getProjPath(ip, "spectra/spectrum_indices.bin"));

  specProtAlignParams.setValue("OUTPUT_PSM_MOD_FREQS",
                               getProjPath(ip,
                                           "homology/spectrum_psm_mod_freqs.txt"));
  specProtAlignParams.setValue("OUTPUT_PSM_FDR",
                               getProjPath(ip,
                                           "homology/spectrum_psm_fdr.txt"));

  specProtAlignParams.addIfExists(ip, "DEBUG_SPECPROTALIGN");
  specProtAlignParams.addIfExists(ip, "DEBUG_SPECPROTALIGN_SPECTRA");
  specProtAlignParams.addIfExists(ip, "DEBUG_SPECPROTALIGN_TIME");
  specProtAlignParams.addIfExists(ip, "DEBUG_SPECPROTALIGN_SPRINKLE");

  DEBUG_VAR(contigSpectra.size());

  DEBUG_TRACE;
  ExecSpecProtAlignTgtDecoy moduleSpecProtAlignTgtDecoy(specProtAlignParams,
                                                        &contigSpectra,
                                                        &prmSpectra,
                                                        &db,
                                                        &dbDecoy,
                                                        &psmTag,
                                                        &psmTagDecoy,
                                                        &penaltyMatrixBlosum,
                                                        &penaltyMatrixMods,
                                                        &matchedSpectraAll,
                                                        &matchedSpectra,
                                                        &psmSet,
                                                        &psmSetDecoy,
                                                        &psmSetFdr);
  DEBUG_TRACE;

  string errorString;
  bool isValid = moduleSpecProtAlignTgtDecoy.validateParams(errorString);
  TEST_VALID;

  bool returnStatus = moduleSpecProtAlignTgtDecoy.invoke();

  // Test for return status
  TEST_RETURN_STATUS("ExecSpecProtAlignTgtDecoy");

  DEBUG_VAR(matchedSpectra.size());

  returnStatus = moduleSpecProtAlignTgtDecoy.saveOutputData();
  DEBUG_VAR(returnStatus);
  TEST_SAVE_OUPUT_DATA("ExecSpecProtAlignTgtDecoy");

  return true;
}

//-----------------------------------------------------------------------------

bool performProtProtAlign(ParameterList & ip,
                          unsigned int refProtIdx,
                          set<unsigned int> dbIndexes,
                          DB_fasta & db)
{
  DEBUG_TRACE;
  ParameterList protProtAlignParams;

  protProtAlignParams.addIfExists(ip, "PROJECT_DIR");
  protProtAlignParams.addIfExists(ip, "OUTPUT_SPECTRA_PATH");

  //protProtAlignParams.addIfExists(ip,"CLUSTALW_EXE_DIR");
  if (ip.exists("EXE_DIR"))
  protProtAlignParams.setValue("CLUSTALW_EXE_DIR", ip.getValue("EXE_DIR"));

  ostringstream tmp;
  tmp << refProtIdx;
  protProtAlignParams.setValue("REFERENCE_PROTEIN_IDX", tmp.str());

  protProtAlignParams.addIfExists(ip, "CLUSTALW_MINSCORE");
  protProtAlignParams.addIfDoesntExist("CLUSTALW_MINSCORE", "250");

  protProtAlignParams.setValue("CLUSTALW_FILENAME_PREFIX",
      getProjPath(ip, "homology/cwseqs_"));
  protProtAlignParams.setValue("CLUSTALW_INDEX",
      getProjPath(ip, "homology/cwindex.txt"));
  protProtAlignParams.writeToFile("debug_protprotalign.params");

  DEBUG_TRACE;
  vector<string> alnFileNames;
  ExecProtProtAlign moduleProtProtAlign(protProtAlignParams,
      &db,
      &dbIndexes,
      &alnFileNames);
  DEBUG_TRACE;

  string errorString;
  bool isValid = moduleProtProtAlign.validateParams(errorString);
  TEST_VALID;

  bool returnStatus = moduleProtProtAlign.invoke();
  // Test for return status
  TEST_RETURN_STATUS("ExecProtProtAlign");
  // Save output data
  returnStatus = moduleProtProtAlign.saveOutputData();
  TEST_SAVE_OUPUT_DATA("ExecProtProtAlign");

  return true;
}

//-----------------------------------------------------------------------------

bool performHomologyAssembly(ParameterList & ip,
                             SpecSet & spectra,
                             DB_fasta & db,
                             SpecSet & contigShifts,
                             SpecSet & contigMatchedIndices,
                             vector<vector<int> > & contigMatchedProts)
{
  DEBUG_TRACE;
  ParameterList homologyAssemblyParams;
  homologyAssemblyParams.setValue("INPUT_HOMOLOGIES",
      getProjPath(ip, "homology/cwindex.txt"));
  //  homologyAssemblyParams.setValue("INPUT_SPECS_NAMES",      getProjPath(ip, "homology/ref_sps_names.txt") );
  homologyAssemblyParams.setValue("SPEC_TYPE_MSMS", "0");
  homologyAssemblyParams.setValue("GRAPH_TYPE", "2");
  homologyAssemblyParams.setValue("MIN_CONTIG_SET", "1");
  homologyAssemblyParams.setValue("EDGE_SCORE_TYPE", "1");
  homologyAssemblyParams.setValue("OUTPUT_SPECS",
      getProjPath(ip,
          "homology/homglue_matches.pklbin"));
  homologyAssemblyParams.setValue("OUTPUT_MATCHES_REF_MIDX",
      getProjPath(ip,
          "homology/homglue_ref_midx.pklbin"));
  homologyAssemblyParams.setValue("OUTPUT_MATCHES_REF_MP",
      getProjPath(ip, "homology/homglue_ref_mp.bin"));
  homologyAssemblyParams.setValue("OUTPUT_MATCHES_CSPS_MIDX",
      getProjPath(ip,
          "homology/homglue_matches_midx.pklbin"));
  homologyAssemblyParams.setValue("OUTPUT_MATCHES_CSPS_MP",
      getProjPath(ip,
          "homology/homglue_matches_mp.bin"));
  homologyAssemblyParams.setValue("OUTPUT_CSV",
      getProjPath(ip, "homglue_matches.txt"));
  homologyAssemblyParams.addIfExists(ip, "RESOLUTION");
  homologyAssemblyParams.addIfExists(ip, "MAX_MOD_MASS");
  homologyAssemblyParams.addIfExists(ip, "TOLERANCE_PEAK");
  homologyAssemblyParams.addIfExists(ip, "TOLERANCE_PM");
  homologyAssemblyParams.addIfExists(ip, "MAX_AA_JUMP");
  homologyAssemblyParams.addIfExists(ip, "PENALTY_PTM");
  homologyAssemblyParams.addIfExists(ip, "PENALTY_SAME_VERTEX");

  homologyAssemblyParams.addIfExists(ip, "PROJECT_DIR");
  homologyAssemblyParams.addIfExists(ip, "OUTPUT_SPECTRA_PATH");

  homologyAssemblyParams.writeToFile(getProjPath(ip, "debug_homology.params"));

  DEBUG_TRACE;
  ExecHomologyAssembly moduleHomologyAssembly(homologyAssemblyParams,
      &spectra,
      &db,
      &contigShifts,
      &contigMatchedIndices,
      &contigMatchedProts);
  DEBUG_TRACE;

  string errorString;
  bool isValid = moduleHomologyAssembly.validateParams(errorString);
  TEST_VALID;

  bool returnStatus = moduleHomologyAssembly.invoke();
  // Test for return status
  TEST_RETURN_STATUS("ExecHomologyAssembly");

  returnStatus = moduleHomologyAssembly.saveOutputData();
  TEST_SAVE_OUPUT_DATA("ExecHomologyAssembly");

  return true;
}

//-----------------------------------------------------------------------------

bool performReport(ParameterList & ip)
{
  DEBUG_MSG("Entering performReport()");

  string execMode = ip.getValue("EXECUTION_MODE", "sps");

  float resolution = ip.getValueDouble("RESOLUTION", 0.1);

  // Get current directory -- needed for absolute path composition
  string currentWorkindDir;
  currentWorkindDir = getCurrentDirectory();

  // auxiliary string used for file name composition using project name
  string aux;

  // Exec Spsplot section

  ParameterList reportSpsplotParams;

  reportSpsplotParams.addIfExists(ip, "PROJECT_DIR");
  reportSpsplotParams.addIfExists(ip, "EXE_DIR");
  reportSpsplotParams.addIfExists(ip, "OUTPUT_SPECTRA_PATH");

  reportSpsplotParams.addIfExists(ip, "REPORT_JOB");
  reportSpsplotParams.addIfExists(ip, "REPORT_USER");
  reportSpsplotParams.addIfExists(ip, "REPORT_PWD");
  // add dynamic report project directory (target upon relocation)
  reportSpsplotParams.addIfExists(ip, "REPORT_SERVER");
  reportSpsplotParams.addIfExists(ip, "REPORT_DIR_SERVER");
  reportSpsplotParams.addIfExists(ip, "REPORT_CELLS_PER_LINE");
  // specify if MS/MS images are shown
  reportSpsplotParams.addIfExists(ip, "REPORT_MSMS_IMAGES");
  // add PDF reports status
  reportSpsplotParams.addIfExists(ip, "REPORT_PDF");

  // add mass shift for report drawing
  reportSpsplotParams.addIfExists(ip, "REPORT_MASS_SHIFT");
  // add mass shift for report drawing
  reportSpsplotParams.addIfExists(ip, "MS2_SCORING_MODEL");
  // add amino acids file
  reportSpsplotParams.addIfExists(ip, "AMINO_ACID_MASSES");
  // add mass shift for report drawing
  //reportSpsplotParams.addIfExists(ip, "REPORT_MASS_SHIFT_PRM");
  // add mass shift for report drawing
  //reportSpsplotParams.addIfExists(ip, "REPORT_MS2_SCORING_MODEL_PRM");

  // add dynamic reports status
  //reportSpsplotParams.addIfExists(ip, "REPORT_DYNAMIC");
  if (ip.getValueInt("REPORT_DYNAMIC") == 1)
    reportSpsplotParams.setValue("REPORT_DYNAMIC", "1");

  if (ip.getValueInt("CLUSTER_MIN_SIZE") == 0
      && ip.getValueInt("MERGE_SAME_PREC", 0) > 1)
    reportSpsplotParams.setValue("NO_CLUSTERS", "1");

  // allow the realign in protein coverage pages
  if (ip.getValueInt("REPORT_REALIGN") == 1)
    reportSpsplotParams.setValue("ALLOW_REALIGN", "1");
  else
    reportSpsplotParams.setValue("ALLOW_REALIGN", "0");

  // allow interactivity with the server
  if (ip.getValueInt("REPORT_NOSERVER") == 1)
    reportSpsplotParams.setValue("DYNAMIC", "0");
  else
    reportSpsplotParams.setValue("DYNAMIC", "1");

  stringstream saux;
  saux << resolution;
  reportSpsplotParams.setValue("RESOLUTION", saux.str());

  int tool = 1;

  DEBUG_VAR(ip.getValueInt("GENOMS_FLAG"));
  DEBUG_VAR(ip.getValueInt("MERGE_FLAG"));

  if (ip.getValueInt("GENOMS_FLAG") == 1)
    tool = 2;

  if (ip.getValueInt("MERGE_FLAG") == 1)
    tool = 3;

  DEBUG_VAR(tool);

  reportSpsplotParams.setValue("TOOL", parseInt(tool));

  DEBUG_VAR(reportSpsplotParams.getValueInt("TOOL"));

  // add number of cpus
  if (ip.exists("CPUS"))
  {
    reportSpsplotParams.addIfExists(ip, "CPUS");
    //} else if(ip.exists("GRID_NUMCPUS")) {
    //  reportSpsplotParams.setValue("CPUS", ip.getValue("GRID_NUMCPUS") );
  }

  // Get report dir -- for the new reports
  string report_dir;
  string report_dir_aux = "report";
  if (ip.exists("REPORT_DIR"))
    report_dir_aux = ip.getValue("REPORT_DIR");

  // If report directory is a relative path, make it absolute
  if (report_dir_aux[0] != '/')
  {
    report_dir = currentWorkindDir;
    if (report_dir[report_dir.size() - 1] != '/')
      report_dir += '/';
  }
  report_dir += report_dir_aux;

  // set the report directory parameter
  reportSpsplotParams.setValue("OUTDIR", report_dir);

  aux = currentWorkindDir;
  reportSpsplotParams.setValue("FONT_PATH", aux);

  //aux = currentWorkindDir;
  //aux += "/spectra/stars_only.pklbin";
  aux = "spectra/stars_only.pklbin";
  reportSpsplotParams.setValue("OUTPUT_STARS", aux);
  //aux = currentWorkindDir;
  //aux += "/spectra/stars_indices.bin";
  aux = "spectra/stars_indices.bin";
  reportSpsplotParams.setValue("OUTPUT_STARS_INDEX", aux);
  //aux = currentWorkindDir ; aux += "/spectra/stars.pklbin";
  //reportSpsplotParams.setValue("OUTPUT_STARS_ALL",    aux);

  aux = ip.getValue("FASTA_DATABASE");
  string aux2;
  if (aux[0] != '/')
  {
    aux2 = currentWorkindDir;
    if (aux2[aux2.size() - 1] != '/')
      aux2 += '/';
  }
  aux2 += aux;
  reportSpsplotParams.setValue("FILE_FASTA", aux);

  if (execMode == string("sps"))
  {
    DEBUG_MSG("sps selected");

    //aux = currentWorkindDir;
    //aux += "/homology/ref_sps_names.txt";
    aux = "homology/ref_sps_names.txt";
    reportSpsplotParams.setValue("FILE_REFINDEX", aux);
    //aux = currentWorkindDir;
    //aux += "/assembly/component_info.bin";
    aux = "assembly/component_info.bin";
    reportSpsplotParams.setValue("FILE_COMP", aux);
    //aux = currentWorkindDir;
    //aux += "/assembly/sps_seqs.pklbin";
    aux = "assembly/sps_seqs.pklbin";
    reportSpsplotParams.setValue("FILE_SEQS", aux);

    // for the old reports
    //aux = currentWorkindDir;
    //aux += "/homology/contigs_mp_all.bin";
    aux = "homology/contigs_mp_all.bin";
    reportSpsplotParams.setValue("FILE_MP", aux);
    //aux = currentWorkindDir;
    //aux += "/homology/contigs_midx_all.pklbin";
    aux = "homology/contigs_midx_all.pklbin";
    reportSpsplotParams.setValue("FILE_MIDX", aux);

    // for the new reports
    //aux = currentWorkindDir;
    //aux += "/homology/contigs_mp.bin";
    aux = "homology/contigs_mp.bin";
    reportSpsplotParams.setValue("FILE_MP2", aux);
    //aux = currentWorkindDir;
    //aux += "/homology/contigs_midx.pklbin";
    aux = "homology/contigs_midx.pklbin";
    reportSpsplotParams.setValue("FILE_MIDX2", aux);

    //aux = currentWorkindDir;
    //aux += "/homology/homglue_ref_mp.bin";
    aux = "homology/homglue_ref_mp.bin";
    reportSpsplotParams.setValue("FILE_REFMP", aux);
    //aux = currentWorkindDir;
    //aux += "/homology/homglue_ref_midx.pklbin";
    aux = "homology/homglue_ref_midx.pklbin";
    reportSpsplotParams.setValue("FILE_REFMIDX", aux);

    //aux = currentWorkindDir;
    //aux += "/spectra/contigs.pklbin";
    aux = "spectra/contigs.pklbin";
    reportSpsplotParams.setValue("MATCHED_CONTIGS", aux);
    //aux = currentWorkindDir;
    //aux += "/spectra/contigs_indices.bin";
    aux = "spectra/contigs_indices.bin";
    reportSpsplotParams.setValue("MATCHED_CONTIGS_IDX", aux);

  }
  else
  {

    DEBUG_MSG("snets selected");

    //aux = currentWorkindDir;
    //aux += "/homology/ref_snets_names.txt";
    aux = "homology/ref_snets_names.txt";
    reportSpsplotParams.setValue("FILE_REFINDEX", aux);
    //aux = currentWorkindDir;
    //aux += "/homology/specnets_abruijn.bin";
    aux = "homology/specnets_abruijn.bin";
    reportSpsplotParams.setValue("FILE_COMP", aux);
    //    aux = currentWorkindDir ; aux += "/specnets/snets_specs.pklbin";
    //aux = currentWorkindDir;
    //aux += "/homology/homglue_matches.pklbin";
    aux = "homology/homglue_matches.pklbin";
    reportSpsplotParams.setValue("FILE_SEQS", aux);

    // for the old reports

    //    aux = currentWorkindDir ; aux += "/specnets/snets_mp.bin";
    //aux = currentWorkindDir;
    //aux += "/homology/homglue_matches_mp.bin";
    aux = "homology/homglue_matches_mp.bin";
    reportSpsplotParams.setValue("FILE_MP", aux);
    //    aux = currentWorkindDir ; aux += "/specnets/snets_midx.pklbin";
    //aux = currentWorkindDir;
    //aux += "/homology/homglue_matches_midx.pklbin";
    aux = "homology/homglue_matches_midx.pklbin";
    reportSpsplotParams.setValue("FILE_MIDX", aux);

    // for the new reports

    //    aux = currentWorkindDir ; aux += "/specnets/snets_mp.bin";
    //aux = currentWorkindDir;
    //aux += "/homology/homglue_matches_mp.bin";
    aux = "homology/homglue_matches_mp.bin";
    reportSpsplotParams.setValue("FILE_MP2", aux);
    //    aux = currentWorkindDir ; aux += "/specnets/snets_midx.pklbin";
    //aux = currentWorkindDir;
    //aux += "/homology/homglue_matches_midx.pklbin";
    aux = "homology/homglue_matches_midx.pklbin";
    reportSpsplotParams.setValue("FILE_MIDX2", aux);

    // Same as FILE_MP / FILE_MIDX - no homology mapping for spectral networks
    //aux = currentWorkindDir;
    //aux += "/homology/homglue_matches_mp.bin";
    aux = "homology/homglue_matches_mp.bin";
    reportSpsplotParams.setValue("FILE_REFMP", aux);
    //aux = currentWorkindDir;
    //aux += "/homology/homglue_matches_midx.pklbin";
    aux = "homology/homglue_matches_midx.pklbin";
    reportSpsplotParams.setValue("FILE_REFMIDX", aux);

    //aux = currentWorkindDir;
    //aux += "/homology/homglue_matches.pklbin";
    aux = "homology/homglue_matches.pklbin";
    reportSpsplotParams.setValue("MATCHED_CONTIGS", aux);
    //    aux = currentWorkindDir ; aux += "/spectra/contigs_indices.bin";
    //    reportSpsplotParams.setValue("MATCHED_CONTIGS_IDX",     aux);
  }

  //aux = currentWorkindDir;
  //aux += "/spectra/stars.pklbin";
  aux = "spectra/stars.pklbin";
  reportSpsplotParams.setValue("FILE_STARS", aux);
  //aux = currentWorkindDir;
  //aux += "/spectra/specs_ms.pklbin";
  if (ip.getValue("CLUSTER_TOOL", CLUSTER_MSCLUST) == CLUSTER_PRMS)
    aux = SPECS_SCORED_PATH;
  else
    aux = SPECS_MS_PATH;
  reportSpsplotParams.setValue("FILE_MS", aux);
  //aux = currentWorkindDir;
  //aux += "/spectra/input_index.txt";
  aux = "spectra/input_index.txt";
  reportSpsplotParams.setValue("FILE_INDEX", aux);

  //aux = currentWorkindDir;
  //aux += "/spectra/out/clust/clusters_0_";
  aux = "spectra/out/clust/clusters_0_";
  reportSpsplotParams.setValue("FILE_CLUSTER", aux);
  //aux = currentWorkindDir;
  //aux += '/';
  //aux += DEFAULT_INPUT_FILES_LIST;
  aux = DEFAULT_INPUT_FILES_LIST;
  reportSpsplotParams.setValue("FILE_CLUSTERMS", aux);
  //aux = currentWorkindDir;
  //aux += '/';
  //aux += DEFAULT_BIN_FILES_FILENAME;
  aux = DEFAULT_BIN_FILES_FILENAME;
  reportSpsplotParams.setValue("FILE_CLUSTERSCAN", aux);

  // write params to debug file
  reportSpsplotParams.writeToFile(getProjPath(ip, "debug_report.params"));

  // instatiate reports module
  DEBUG_TRACE;
  ExecReportSpsplot moduleReportSpsplot(reportSpsplotParams);
  DEBUG_TRACE;

  // test parameters
  string errorString;
  bool isValid = moduleReportSpsplot.validateParams(errorString);
  TEST_VALID;

  DEBUG_MSG("Invoking reports module");

  bool returnStatus = moduleReportSpsplot.invoke();
  // Test for return status
  TEST_RETURN_STATUS("ExecReportSpsplot");

  string source = ip.getValue("EXE_DIR");
  source += "/resources/css";
  string dest = ip.getValue("REPORT_DIR");
  CopyDirectory(source, dest);

  return true;
}

//-----------------------------------------------------------------------------
bool performExecMainSpecnets(ParameterList & ip,
                             SpecSet * msSpectra,
                             SpecSet * scoredSpectra,
                             SpecSet * starSpectra,
                             SpectrumPairSet * pairs,
                             DB_fasta * db,
                             PenaltyMatrix * penaltyMatrixBlosum,
                             PenaltyMatrix * penaltyMatrixMods,
                             PeptideSpectrumMatchSet * psms,
                             PeptideSpectrumMatchSet * origPsms,
                             SpecSet * psms_spectra,
                             SpecSet * psms_midx,
                             vector<vector<int> > * psms_mp,
                             SpecSet * snets_contigs,
                             SpecSet * snets_midx,
                             vector<vector<int> > * snets_mp)
{
  ParameterList specnetsParams;
  specnetsParams.addIfExists(ip, "OUTPUT_SPECTRA_PATH");
  specnetsParams.setValue("RESOLUTION", ip.getValue("RESOLUTION", "0.1"));
  specnetsParams.setValue("MIN_MOD_MASS", ip.getValue("MIN_MOD_MASS", "-100"));
  specnetsParams.setValue("MAX_MOD_MASS", ip.getValue("MAX_MOD_MASS", "100"));
  specnetsParams.setValue("TOLERANCE_PEAK",
                          ip.getValue("TOLERANCE_PEAK", "0.45"));
  specnetsParams.setValue("TOLERANCE_PM", ip.getValue("TOLERANCE_PM", "1.5"));

  //alignment params
  specnetsParams.addIfExists(ip, "MAX_PARSIMONY");
  specnetsParams.addIfExists(ip, "MAX_NUM_MODS");
  specnetsParams.addIfExists(ip, "MIN_MATCHED_PEAKS");
  specnetsParams.addIfExists(ip, "MATCHES_PER_MOD");
  specnetsParams.addIfExists(ip, "MIN_MATCHED_PEAKS_DB");
  specnetsParams.addIfExists(ip, "PENALTY_ALIGNMENT");
  specnetsParams.addIfExists(ip, "ENFORCE_ENDPEAKS");
  specnetsParams.addIfExists(ip, "PENALTY_ALIGNMENT_ALPHA");
  specnetsParams.addIfExists(ip, "PENALTY_ALIGNMENT_BETA");
  specnetsParams.addIfExists(ip, "MIN_PENALTY_FREQUENCY");
  specnetsParams.addIfExists(ip, "PEAK_EQUIVALENTS");
  specnetsParams.addIfExists(ip, "MAX_ALIGN_DB_GAP_AAS");
  specnetsParams.addIfExists(ip, "MAX_ALIGN_SPECTRUM_GAP_DALTONS");

  // Tagsearch params
  specnetsParams.addIfExists(ip, "INSPECT_PSMS");
  specnetsParams.addIfExists(ip, "INPUT_SPEC_IDS");
  specnetsParams.addIfExists(ip, "INPUT_CLUSTERS_DIR");
  specnetsParams.setValue("OUTPUT_PEPTIDES", "spectra/tagsearch_peptides.txt");
  specnetsParams.addIfExists(ip, "OUTPUT_RAW_RESULTS");
  specnetsParams.addIfExists(ip, "OUTPUT_SNETS_PTM_MATRIX");
  specnetsParams.addIfExists(ip, "OUTPUT_SPECS_PROJ");
  specnetsParams.addIfExists(ip, "OUTPUT_ALIGNS");
  specnetsParams.addIfExists(ip, "OUTPUT_ANNOTINFO");

  // SVM params
  specnetsParams.setValue("STARS_SCORED_PRM_OFFSET", "-1.0072763");
  specnetsParams.setValue("STARS_SCORED_SRM_OFFSET", "-1.0072763");
  specnetsParams.setValue("SCAN_ZERO_INDEX", "1");
  specnetsParams.addIfExists(ip, "SCAN_ZERO_INDEX");
  specnetsParams.addIfExists(ip, "MS2_SCORING_MODEL"); // To become "<trunk>/resources/dancik_model.txt
  specnetsParams.addIfExists(ip, "SPECS_MS_STATISTICS_CONFIG"); // To become "<trunk>/resources/spectra_stats_ms.txt
  specnetsParams.addIfExists(ip, "SPECS_SCORED_STATISTICS_CONFIG"); // To become "<trunk>/resources/spectra_stats_prm.txt
  specnetsParams.addIfExists(ip, "STARS_STATISTICS_CONFIG"); // To become "<trunk>/resources/spectra_stats_stars.txt
  specnetsParams.addIfExists(ip, "SVM_SCALE_PARAMS_CHARGE1"); // To become "<trunk>/resources/HEK293_charge1_model_svm_keys_range.txt
  specnetsParams.addIfExists(ip, "SVM_SCALE_PARAMS_CHARGE2"); // To become "<trunk>/resources/HEK293_charge2_model_svm_keys_range.txt
  specnetsParams.addIfExists(ip, "SVM_SCALE_PARAMS_CHARGE3"); // To become "<trunk>/resources/HEK293_charge3_model_svm_keys_range.txt
  specnetsParams.addIfExists(ip, "SVM_MODEL_CHARGE1"); // To become "<trunk>/resources/HEK293_charge1_model_SVM.model
  specnetsParams.addIfExists(ip, "SVM_MODEL_CHARGE2"); // To become "<trunk>/resources/HEK293_charge2_model_SVM.model
  specnetsParams.addIfExists(ip, "SVM_MODEL_CHARGE3"); // To become "<trunk>/resources/HEK293_charge3_model_SVM.model

  // Spectral networks params
  specnetsParams.setValue("SPECNETS_PROJ_TYPE",
                          ip.getValue("SPECNETS_PROJ_TYPE", "all")); // Possible values are "all" / "matched" to retain all/matched-only PRMs during propagation
  specnetsParams.setValue("SPECNETS_USE_SVM",
                          ip.getValue("SPECNETS_USE_SVM", "0"));

  // FDR params
  specnetsParams.setValue("INPUT_RESULTS_TYPE", "snets");
  specnetsParams.setValue("OUTPUT_FDR_RESULTS", "spectra/psms_fdr.txt");

  DEBUG_TRACE;
  SpectrumPairSet filteredPairs;

  specnetsParams.writeToFile("debug_specnets.params");

  MS2ScoringModel model;
  if (ip.exists("MS2_SCORING_MODEL"))
  {
    if (!model.LoadModel(ip.getValue("MS2_SCORING_MODEL").c_str()))
    {
      DEBUG_MSG("Could not load " << ip.getValue("MS2_SCORING_MODEL"));
      return false;
    }
  }
  DEBUG_MSG("SCORING MODEL LOADED");

  DEBUG_TRACE;
  ExecMainSpecnets mainSpecnets(specnetsParams,
                                msSpectra,
                                scoredSpectra,
                                starSpectra,
                                pairs,
                                &model,
                                db,
                                penaltyMatrixBlosum,
                                penaltyMatrixMods,
                                psms,
                                origPsms,
                                psms_spectra,
                                psms_midx,
                                psms_mp,
                                snets_contigs,
                                snets_midx,
                                snets_mp);
  DEBUG_TRACE;

  string errorString;
  bool isValid = mainSpecnets.validateParams(errorString);
  TEST_VALID;

  bool returnStatus = mainSpecnets.invoke();
  // Test for return status
  TEST_RETURN_STATUS("ExecMainSpecnets");

  // Saving output data
  returnStatus = mainSpecnets.saveOutputData();
  TEST_SAVE_OUPUT_DATA("ExecMainSpecnets");

  ofstream spsIndex(getProjPath(ip, "homology/ref_snets_names.txt").c_str(),
                    ios_base::out | ios_base::binary);
  for (unsigned int i = 0; i < snets_contigs->size(); i++)
    spsIndex << "snet:" << i + 1 << endl;
  spsIndex.close();

  //---------------------------------------------------------------------------
  // REPORT STAGE
  //---------------------------------------------------------------------------

  ofstream spsProj(getProjPath(ip, "sps_projects.txt").c_str(),
                   ios_base::out | ios_base::binary);
  spsProj << "sps;.;" << ip.getValue("TOLERANCE_PEAK") << ";"
      << ip.getValue("TOLERANCE_PM") << "\n";
  spsProj.close();

  if (!performReport(ip))
  {
    ERROR_MSG("Problem encountered during Report stage");
    exit(-100 - STAGE_REPORT);
  }

  return true;
}

#if 0
#include <sys/resource.h>
int setStackSize(int newSize)
{
  const rlim_t kStackSize = newSize * 1024 * 1024; // min stack size in MB
  struct rlimit rl;
  int result;

  result = getrlimit(RLIMIT_STACK, &rl);
  DEBUG_MSG("Stack soft limit: " << (int)(rl.rlim_cur) );
  DEBUG_MSG("Stack hard limit: " << (int)(rl.rlim_max) );

  if (result == 0)
  {
    if (rl.rlim_cur < kStackSize)
    {
      rl.rlim_cur = kStackSize;
      if(rl.rlim_cur > rl.rlim_max)
      rl.rlim_cur = rl.rlim_max;
      DEBUG_MSG("Changing stack soft limit to : " << rl.rlim_cur );
      result = setrlimit(RLIMIT_STACK, &rl);
      if (result != 0)
      {
        int aa = errno;
        WARN_MSG("setrlimit returned result = " << result);
        WARN_MSG("errno = " << aa);
        string aux = strerror(aa);
        WARN_MSG("error = " << aux);
        return 0;
      }
      DEBUG_MSG("Stack soft limit changed to : " << kStackSize );
    }
  }

  // ...

  return 1;
}
#endif

int showVersion(void)
{
  cout << PROGRAM_NAME << endl;
  cout << "main_specnets 3.0." << XSTR(SPS_VERSION) << endl;
  cout << "Build date: " << __DATE__ << " " << __TIME__ << endl;
  cout << "SH1: " << XSTR(GIT_SH1) << endl;

  cout << endl;

  cout << COPYRIGHT1 << endl;
  cout << COPYRIGHT2 << endl;
  cout << endl;

  return 0;
}

// MAIN
//-----------------------------------------------------------------------------
int main(int argc, char ** argv)
{
  int logLevel = 5;
  Logger::setDefaultLogger(Logger::getLogger(logLevel));

  string initialStageString("");
  string finalStageString("");
  string statusFileName("status.txt");

  bool resumeFlag = false;
  bool gridExecutionFlag = false;
  bool runGenoMSFlag = false;
  bool runMergedFlag = false;

  bool showHelp = false;
  for (size_t i = 0; i < argc; i++)
  {
    string arg(argv[i]);
    if (arg.compare("--help") == 0)
      showHelp = true;

    if (arg.compare("--version") == 0)
      return showVersion();
  }

  if (argc < 2 || showHelp)
  {
    if (showHelp)
    {
      cout << "Usage: main_specnets [PARAM FILE] [OPTION]..." << endl << endl;
      cout << "Optional arguments are listed below " << endl;
      cout << "  -i  <intialstage>   begin processing at specified stage:"
          << endl;
      cout
          << "                         ms2deconv, mscluster, scoring, filterpairs"
          << endl;
      cout
          << "                         filteraligns, alignment, filterstarpairs"
          << endl;
      cout
          << "                         assembly, metaassembly, tagsearch, specnets, contigprotalign"
          << endl;
      cout
          << "                         specprotalign, protprotalign, homologyassembly, genoms, statprotseqs, report"
          << endl;
      cout
          << "  -f  <finalstage>   end processing after completing specified stage:"
          << endl;
      cout << "                         mscluster, scoring, filterpairs"
          << endl;
      cout
          << "                         filteraligns, alignment, filterstarpairs"
          << endl;
      cout
          << "                         assembly, tagsearch, specnets, contigprotalign"
          << endl;
      cout
          << "                         specprotalign, protprotalign, homologyassembly, genoms, report"
          << endl;
      cout << "  -proteosafe         input params are passed by ProteoSAFe"
          << endl;
      cout << "  -g                  execution is on a grid" << endl;
      cout << "  -lf <filename>      name of log file for output" << endl;
      cout << "  -ll <loglevel>      log level for debug/warn/error output:"
          << endl;
      cout << "                         9 for errors only" << endl;
      cout << "                         5 for warnings and errors" << endl;
      cout << "                         0 for all debug output" << endl;
      cout << "  -s                   execute a single step then exit" << endl;
      cout << "  -z                   resume an earlier run and prepare"
          << endl;
      cout
          << "                       specprotalign for running on a remote grid"
          << endl;
      cout << "       NOTE: this option is only used when intermediate results"
          << endl;
      cout << "             were saved during a previous run where FilterPairs"
          << endl;
      cout << "             was executed on a remote grid" << endl;
      cout << "  -zno                 resume an earlier run and execute"
          << endl;
      cout << "                       remaining steps on a single thread"
          << endl;
      cout << "       NOTE: this option is only used when intermediate results"
          << endl;
      cout << "             were saved during a previous run where FilterPairs"
          << endl;
      cout << "             was executed on a remote grid" << endl;
      cout << "  -q                   Run in GenoMS mode" << endl;
      cout << "       NOTE: In this mode, only 4 stages are valid: mscluster, "
          << endl;
      cout << "             scoring, genoms, and report" << endl;
      cout << "  -m                   Run in integrated mode (GenoMS + CSPS)"
          << endl;
      cout
          << "  -x  <execmode>        execution mode (sps, specnets, or signatures):"
          << endl << endl;
    }
    else
    {
      cerr << "main_specnets: insufficient arguments" << endl;
      cerr << "Type \'main_specnets --help\' for more information." << endl
          << endl;
    }

    cout << PROGRAM_NAME << endl;
    cout << "main_specnets 3.0." << XSTR(SPS_VERSION) << endl;
    cout << endl;

    cout << COPYRIGHT1 << endl;
    cout << COPYRIGHT2 << endl;
    cout << endl;

    return -1;
  }

  // Parse the command line parameters
  vector<CommandLineParser::Option> listOptions;
  listOptions.push_back(CommandLineParser::Option("i", "INITIAL_STAGE", 1));
  listOptions.push_back(CommandLineParser::Option("f", "FINAL_STAGE", 1));
  listOptions.push_back(CommandLineParser::Option("g", "GRID_EXECUTION", 0));
  listOptions.push_back(CommandLineParser::Option("lf", "LOG_FILE_NAME", 1));
  listOptions.push_back(CommandLineParser::Option("ll", "LOG_LEVEL", 1));
  listOptions.push_back(CommandLineParser::Option("s", "SINGLE_STEP", 0));
  listOptions.push_back(CommandLineParser::Option("z", "RESUME_FLAG", 0));
  listOptions.push_back(CommandLineParser::Option("zno",
                                                  "RESUME_SINGLE_FLAG",
                                                  0));
  listOptions.push_back(CommandLineParser::Option("q", "GENOMS_FLAG", 0));
  listOptions.push_back(CommandLineParser::Option("m", "MERGE_FLAG", 0));
  listOptions.push_back(CommandLineParser::Option("x", "EXECUTION_MODE", 1));
  listOptions.push_back(CommandLineParser::Option("ccms_numnodes",
                                                  "GRID_NUMNODES",
                                                  1));

  CommandLineParser clp(argc, argv, 1, listOptions);
  string parserError = "";
  if (!clp.validate(parserError))
  {
    cerr << "main_specnets: " << parserError << endl;
    cerr << "Type \'main_specnets --help\' for more information." << endl
        << endl;
    cout << PROGRAM_NAME << endl;
    cout << "main_specnets 3.0." << XSTR(SPS_VERSION) << endl;
    cout << endl;
    cout << COPYRIGHT1 << endl;
    cout << COPYRIGHT2 << endl;
    cout << endl;
    return -1;
  }

  ParameterList commandLineParams;
  clp.getOptionsAsParameterList(commandLineParams);
  if (commandLineParams.exists("GENOMS_FLAG"))
    commandLineParams.setValue("GENOMS_FLAG", "1");

  if (commandLineParams.exists("MERGE_FLAG"))
    commandLineParams.setValue("MERGE_FLAG", "1");

  ParameterList ip;
  ip.readFromFile(argv[1]);
  ip.writeToFile("debug_sps.params");

  // Combine the command line parameters to the file ones
  //   Command line parameters take precedence (hence the overwrite flag set)
  ip.addList(commandLineParams, true);
  ip.writeToFile("debug_wcommand.params");

  string execMode = ip.getValue("EXECUTION_MODE", "sps");

  if (ip.getValueInt("GENOMS_FLAG", 0) == 1)
    runGenoMSFlag = true;

  if (ip.getValueInt("MERGE_FLAG", 0) == 1)
    runMergedFlag = true;

  logLevel = ip.getValueInt("LOG_LEVEL", 5);
  if (ip.exists("LOG_FILE_NAME"))
  {
    string logFileName = ip.getValue("LOG_FILE_NAME");
    Logger::setDefaultLogger(Logger::getLogger(logFileName, logLevel));
  }
  else
  {
    Logger::setDefaultLogger(Logger::getLogger(logLevel));
  }

  // add divert for segfault
  addSegFaultDivert();

  DEBUG_VAR(execMode);
  DEBUG_MSG("GenoMS_FLAG=" + ip.getValue("GENOMS_FLAG", "X"));
  DEBUG_MSG("Merge_FLAG=" + ip.getValue("MERGE_FLAG", "X"));

  DEBUG_TRACE;

  DEBUG_MSG("SH1:" << XSTR(GIT_SH1));

  if (!ip.exists("EXE_DIR"))
  {

    // extract EXE_DIR from command line
    string exeDir(argv[0]);

    // find last /, and remove from that point on
    size_t found = exeDir.find_last_of("/\\");
    string aux = exeDir.substr(0, found);

    //string mainSpecnetsStr = "/main_specnets";
    //exeDir.erase(exeDir.length() - mainSpecnetsStr.length(), mainSpecnetsStr.length());

    // remove /ExecFramework, if it exists
    string mainSpecnetsStr = "/ExecFramework";
    if (aux.rfind(mainSpecnetsStr) == aux.length() - mainSpecnetsStr.length())
      aux.erase(aux.length() - mainSpecnetsStr.length(),
                mainSpecnetsStr.length());
    else
    {
      mainSpecnetsStr = "\\ExecFramework";
      if (aux.rfind(mainSpecnetsStr) == aux.length() - mainSpecnetsStr.length())
        aux.erase(aux.length() - mainSpecnetsStr.length(),
                  mainSpecnetsStr.length());
    }

    ip.setValue("EXE_DIR", aux);
  }

  if (ip.exists("INITIAL_STAGE"))
  {
    initialStageString = commandLineParams.getValue("INITIAL_STAGE");
  }
  DEBUG_VAR(initialStageString);

  if (ip.exists("FINAL_STAGE"))
  {
    finalStageString = commandLineParams.getValue("FINAL_STAGE");
  }
  DEBUG_VAR(finalStageString);

  if (ip.exists("MIN_PVALUE"))
  {
    ERROR_MSG("MIN_PVALUE is a deprecated variable, must use MAX_PVALUE instead");
    return -1;
  }

  addDefaultParameterValues(ip);
  ip.writeToFile("debug_default.params");

  ip.addIfDoesntExist("OUTPUT_SPECTRA_PATH", getProjPath(ip, "./spectra"));
  if (ip.exists("EXE_DIR"))
  {
    string exeDir = ip.getValue("EXE_DIR");

    // if path begins with '~', exit program.
    if (exeDir[0] == '~')
    {
      cout
          << "EXE_DIR path begins with tilde (~). Paths beginning with tilde are not supported."
          << endl;
      exit(0);
    }

    // In case there is a "/" at the end of EXE_DIR.. remove it
    if (exeDir.length() > 2 && exeDir[exeDir.length() - 1] == '/')
    {
      exeDir = exeDir.substr(0, exeDir.length() - 1);
      ip.setValue("EXE_DIR", exeDir);
    }
  }

  if (ip.exists("GRID_EXE_DIR"))
  {
    // In case there is a "/" at the end of EXE_DIR.. remove it
    string gridExeDir = ip.getValue("GRID_EXE_DIR");
    if (gridExeDir.length() > 2 && gridExeDir[gridExeDir.length() - 1] == '/')
    {
      gridExeDir = gridExeDir.substr(0, gridExeDir.length() - 1);
      ip.setValue("GRID_EXE_DIR", gridExeDir);
    }
  }

  if (ip.exists("GRID_SGE_EXE_DIR"))
  {
    // In case there is a "/" at the end of GRID_SGE_EXE_DIR.. remove it
    string gridSgeExeDir = ip.getValue("GRID_SGE_EXE_DIR");
    if (gridSgeExeDir.length() > 2
        && gridSgeExeDir[gridSgeExeDir.length() - 1] == '/')
    {
      gridSgeExeDir = gridSgeExeDir.substr(0, gridSgeExeDir.length() - 1);
      ip.setValue("GRID_SGE_EXE_DIR", gridSgeExeDir);
    }
  }

  if (ip.exists("INITIAL_STAGE"))
  {
    initialStageString = ip.getValue("INITIAL_STAGE");
  }
  DEBUG_VAR(initialStageString);
  if (initialStageString.empty())
  {
    initialStageString = "begin";
  }
  DEBUG_VAR(initialStageString);

  if (ip.exists("FINAL_STAGE"))
  {
    finalStageString = ip.getValue("FINAL_STAGE");
  }
  DEBUG_VAR(finalStageString);
  if (finalStageString.empty())
  {
    finalStageString = "report";
  }
  DEBUG_VAR(finalStageString);

  map<string, Stage> map_stage;
  map_stage["begin"] = STAGE_BEGIN;
  map_stage["ms2deconv"] = STAGE_MS2DECONV;
  map_stage["mscluster"] = STAGE_MSCLUSTER;
  map_stage["scoring"] = STAGE_SCORING;
  map_stage["filterpairs"] = STAGE_FILTERPAIRS;
  map_stage["filteraligns"] = STAGE_FILTERALIGNS;
  map_stage["alignment"] = STAGE_ALIGNMENT;
  map_stage["filterstarpairs"] = STAGE_FILTERSTARPAIRS;
  map_stage["assembly"] = STAGE_ASSEMBLY;
  map_stage["metaassembly"] = STAGE_METAASSEMBLY;
  map_stage["tagsearch"] = STAGE_TAGSEARCH;
  map_stage["specnets"] = STAGE_SPECNETS;
  map_stage["contigprotalign"] = STAGE_CONTIGPROTALIGN;
  map_stage["specprotalign"] = STAGE_SPECPROTALIGN;
  map_stage["protprotalign"] = STAGE_PROTPROTALIGN;
  map_stage["homologyassembly"] = STAGE_HOMOLOGYASSEMBLY;
  map_stage["genoms"] = STAGE_GENOMS;
  map_stage["merge"] = STAGE_MERGE;
  map_stage["statprotseqs"] = STAGE_STATPROTSEQS;
  map_stage["report"] = STAGE_REPORT;

  if (map_stage.find(initialStageString) == map_stage.end())
  {
    ERROR_MSG("Unknown starting stage [" << initialStageString << "]");
    return -1;
  }

  if (map_stage.find(finalStageString) == map_stage.end())
  {
    ERROR_MSG("Unknown final stage [" << finalStageString << "]");
    return -1;
  }

  //Start the status as "running" and write itout
  writeStatusFile(statusFileName, "Running");

  int initialStage = map_stage[initialStageString];
  DEBUG_VAR(initialStage);
  int finalStage = map_stage[finalStageString];
  DEBUG_VAR(finalStage);

  if (commandLineParams.exists("RESUME_FLAG")
      || commandLineParams.exists("RESUME_SINGLE_FLAG"))
  {
    resumeFlag = true;
  }
  DEBUG_VAR(resumeFlag);

  if (commandLineParams.exists("GRID_EXECUTION"))
  {
    gridExecutionFlag = true;
  }
  DEBUG_VAR(gridExecutionFlag);

  //
  // Initialize environment
  //
  DEBUG_MSG("Spectral Networks 2.0.0: session started on: " << getCurrentTimeString());
  DEBUG_MSG("Starting stage: [" << initialStageString << "] on: " << getCurrentTimeString());

  bool res;
  res = mkdir_if_not_exist("spectra");
  if (res)
  {
    DEBUG_MSG("Made directory \'spectra\'");
  }
  res = mkdir_if_not_exist("aligns");
  if (res)
  {
    DEBUG_MSG("Made directory \'aligns\'");
  }
  if (commandLineParams.exists("PROTEOSAFE_FLAG"))
  {
    res = mkdir_if_not_exist("aligns/params");
    if (res)
    {
      DEBUG_MSG("Made directory \'aligns/params\'");
    }
  }
  res = mkdir_if_not_exist("specnets");
  if (res)
  {
    DEBUG_MSG("Made directory \'specnets\'");
  }
  res = mkdir_if_not_exist("homology");
  if (res)
  {
    DEBUG_MSG("Made directory \'homology\'");
  }
  res = mkdir_if_not_exist("homology/grid_t");
  if (res)
  {
    DEBUG_MSG("Made directory \'homology/grid_t\'");
  }
  res = mkdir_if_not_exist("homology/grid_d");
  if (res)
  {
    DEBUG_MSG("Made directory \'homology/grid_d\'");
  }
  res = mkdir_if_not_exist("homology/grid2_t");
  if (res)
  {
    DEBUG_MSG("Made directory \'homology/grid2_t\'");
  }
  res = mkdir_if_not_exist("homology/grid2_d");
  if (res)
  {
    DEBUG_MSG("Made directory \'homology/grid2_d\'");
  }
  res = mkdir_if_not_exist("report");
  if (res)
  {
    DEBUG_MSG("Made directory \'report\'");
  }
  res = mkdir_if_not_exist("ReportData");
  if (res)
  {
    DEBUG_MSG("Made directory \'ReportData\'");
  }
  res = mkdir_if_not_exist("assembly");
  if (res)
  {
    DEBUG_MSG("Made directory \'assembly\'");
  }

  // get LD_LIBRARY_PATH from system
  curLibPath = getenv("LD_LIBRARY_PATH");

  // Build the needed library path
  string libPath;
  libPath = ip.getValue("EXE_DIR");
  // set LD_LIBRARY_PATH to EXE_DIR + /libs.
  libPath += "/libs";

  string fullLibPath;
  // check if LD_LIBRARY_PATH is already defined.
  if (curLibPath)
  {
    // if it is, check if it contains the path we want.
    fullLibPath = curLibPath;
    // if the library path IS NOT contained in the path variable, add it, and set the environment variable.
    if (fullLibPath.find(libPath) == string::npos)
    {
      fullLibPath += ':';
      fullLibPath += libPath;
      mysetenv("LD_LIBRARY_PATH", fullLibPath.c_str());
    }
  }
  else
  {
    // if LD_LIBRARY_PATH is not defined,, define it with what we want.
    mysetenv("LD_LIBRARY_PATH", libPath.c_str());
  }

  //---------------------------------
  // Load amino acid masses
  //---------------------------------
  AAJumps jumps(1); // Amino acid masses
  if (ip.exists("AMINO_ACID_MASSES"))
  {
    DEBUG_MSG("Loading amino acid masses from [" << ip.getValue("AMINO_ACID_MASSES") << "]");
    if (!jumps.loadJumps(ip.getValue("AMINO_ACID_MASSES").c_str(), true))
    {
      ERROR_MSG("Unable to load amino acid jumps");
      ERROR_MSG("Aborting!");
      exit(-1);
    }
  }
  else
  {
    DEBUG_MSG("No amino acid masses loaded. Using defaults");
  }

  string exeDir = ip.getValue("EXE_DIR");
  string convertCmd = exeDir + "/convert ";

  vector<string> ms2Spectra_separate; // MS/MS spectra in separate files
  string inputPklbinFilesList = getProjPath(ip, DEFAULT_INPUT_FILES_LIST);

  if (initialStage == STAGE_BEGIN)
  {
    DEBUG_TRACE;

    // Build the needed library path
    string libPath = exeDir + "/libs";
    addEnvironmentVariable(convertCmd, "LD_LIBRARY_PATH", libPath);

    bool ret = loadInitialData(ip, ip.getValue("INPUT_SPECS_MS"), ms2Spectra_separate);

    if(!ret)
    {
      ERROR_MSG("Failed loading initial data.");
      writeStatusFile(statusFileName, "Error");
      exit(1);
    }

    if (commandLineParams.exists("SINGLE_STEP"))
    {
      DEBUG_MSG("Option -s given. Exiting after single step.");
      writeStatusFile(statusFileName, "Finished");
      exit(0);
    }
  }
  else
  {
    if (!readFilesFromFile(inputPklbinFilesList, ms2Spectra_separate))
    {
      ERROR_MSG("readFilesFromFile() failed for " << inputPklbinFilesList);
      exit(0);
    }
  }

  if (finalStage == STAGE_BEGIN)
  {
    DEBUG_MSG("Option -f given. Exiting after step: " << finalStageString);
    writeStatusFile(statusFileName, "Finished");
    exit(0);
  }

  //TEST_RETURN("moduleBegin", *ms2Spectra_separate);
  DEBUG_VAR(ms2Spectra_separate.size());

  if (ms2Spectra_separate.size() == 0)
  {
    ERROR_MSG("ms2Spectra_separate size is 0!");
    writeStatusFile(statusFileName, "Error");
    exit(-2);
  }

  vector<string> deconvMs2Spectra_separate; // MS/MS spectra in separate files
  string deconvPklbinFilesList = getProjPath(ip, DEFAULT_DECONV_FILES_LIST);

  if (initialStage <= STAGE_MS2DECONV && ip.getValueInt("DECONV_MS2", 0) > 0)
  {

    if (!performMS2Deconv(ip, ms2Spectra_separate, deconvMs2Spectra_separate))
    {
      exit(0);
    }

    DEBUG_VAR(ms2Spectra_separate.size());

    if (commandLineParams.exists("SINGLE_STEP"))
    {
      DEBUG_MSG("Option -s given. Exiting after single step.");
      writeStatusFile(statusFileName, "Finished");
      exit(0);
    }
  }
  else if (ip.getValueInt("DECONV_MS2", 0) > 0)
  {
    if (!readFilesFromFile(deconvPklbinFilesList, deconvMs2Spectra_separate))
    {
      ERROR_MSG("readFilesFromFile() failed for " << inputPklbinFilesList);
      exit(0);
    }
  }
  else
  {
    deconvMs2Spectra_separate = ms2Spectra_separate;
    deconvPklbinFilesList = inputPklbinFilesList;
  }

  if (finalStage == STAGE_MS2DECONV)
  {
    DEBUG_MSG("Option -f given. Exiting after step: " << finalStageString);
    writeStatusFile(statusFileName, "Finished");
    exit(0);
  }

  string clusterTool = ip.getValue("CLUSTER_TOOL", CLUSTER_MSCLUST);

  bool performClustering = ip.getValueInt("CLUSTER_MIN_SIZE", 0) > 0;
  SpecSet ms2Spectra;

  if (initialStage <= STAGE_MSCLUSTER)
  {

    DEBUG_MSG("Starting stage McCluster on: " << getCurrentTimeString());
    if (!performMsCluster(ip, deconvMs2Spectra_separate, ms2Spectra))
    {
      ERROR_MSG("Problem encountered during MsCluster stage");
      writeStatusFile(statusFileName, "Error");
      exit(-1);
    }

    if (commandLineParams.exists("SINGLE_STEP"))
    {
      DEBUG_MSG("Option -s given. Exiting after single step.");
      writeStatusFile(statusFileName, "Finished");
      exit(0);
    }
  }
  else
  {
    // load clustered spectra
    DEBUG_MSG("Bypassing MsCluster stage");
    if (!ms2Spectra.loadPklBin(getProjPath(ip, "./spectra/specs_ms.pklbin").c_str()))
    {
      ERROR_MSG("Could not find input MS/MS spectra!");
    }

    DEBUG_VAR(ms2Spectra.size());
  }

  // Test for empty return data structures
  TEST_RETURN("moduleMsCluster", ms2Spectra);

  DEBUG_VAR(ms2Spectra.size());

  if (ms2Spectra.size() == 0)
  {
    ERROR_MSG("ms2Spectra size is 0!");
    writeStatusFile(statusFileName, "Error");
    exit(-2);
  }

  if (finalStage == STAGE_MSCLUSTER)
  {
    DEBUG_MSG("Option -f given. Exiting after step: " << finalStageString);
    writeStatusFile(statusFileName, "Finished");
    exit(0);
  }

  //---------------------------------------------------------------------------
  // PEPNOVO/SCORING STAGE
  //---------------------------------------------------------------------------
  SpecSet prmSpectra;
  if (ip.getValue("PAIRS_MATCH_MODE", "") == "cosine")
  {
    prmSpectra = ms2Spectra;
  }
  else
  {
    if (initialStage <= STAGE_SCORING)
    {
      DEBUG_MSG("Starting stage Scoring on: " << getCurrentTimeString());
      if (!performScoring(ip, ms2Spectra, prmSpectra))
      {
        ERROR_MSG("Problem encountered during Scoring stage");
        writeStatusFile(statusFileName, "Error");
        exit(-3);
      }

      if (runGenoMSFlag || runMergedFlag)
      {
        // need to save a copy of the scored spectra as MGF so GenoMS works
        if (!ExecMergeConvert::saveSpecset(getProjPath(ip,
                                                       SPECS_SCORED_MGF_PATH),
                                           &prmSpectra))
        {
          ERROR_MSG("Failed to save to " << getProjPath(ip, SPECS_SCORED_MGF_PATH));
          exit(-3);
        }
      }

      if (commandLineParams.exists("SINGLE_STEP"))
      {
        DEBUG_MSG("Option -s given. Exiting after single step.");
        writeStatusFile(statusFileName, "Finished");
        exit(0);
      }
    }
    else
    {
      DEBUG_MSG("Bypassing Scoring stage");
      if (!prmSpectra.loadPklBin(getProjPath(ip, SPECS_SCORED_PATH).c_str()))
      {
        ERROR_MSG("Could not find input PRM spectra!");
      }
      DEBUG_VAR(prmSpectra.size());
    }
  }

  // Test for empty return data structures
  TEST_RETURN("moduleScoring", prmSpectra);

  if (finalStage == STAGE_SCORING)
  {
    DEBUG_MSG("Option -f given. Exiting after step: " << finalStageString);
    writeStatusFile(statusFileName, "Finished");
    exit(0);
  }

  DEBUG_TRACE;
  float peakTol = ip.getValueDouble("TOLERANCE_PEAK", 0.5);
  bool specTypeMSMS = ip.getValueBool("SPEC_TYPE_MSMS", false);
  float ionOffset = specTypeMSMS ? AAJumps::massHion : 0;
  for (unsigned int i = 0; i < prmSpectra.size(); i++)
  {
    prmSpectra[i].addZPMpeaks(peakTol, ionOffset, true);
  }

  //--------------------------------------------------------------------------
  // GenoMS STAGE
  //--------------------------------------------------------------------------

  if (runGenoMSFlag || runMergedFlag)
  {
    if (initialStage <= STAGE_GENOMS)
    {

#if INCLUDE_GENOMS == 0
      ERROR_MSG("GenoMS is not available. Returning in error.");
      exit(-STAGE_GENOMS);
      return false;
#endif

      DEBUG_MSG("Starting GenoMS");

      if (!performGenoMS(ip))
      {
        ERROR_MSG("Problem encountered during GenoMS");
        writeStatusFile(statusFileName, "Error");
        exit(-1);
      }
    }
    else
    {
      DEBUG_MSG("Bypassing GenoMS stage");
    }
    /*
     SpecSet testSpectra;

     testSpectra.loadPklBin( getProjPath(ip, "./spectra/stars.pklbin").c_str() );
     DEBUG_VAR(testSpectra.size());
     testSpectra.SaveSpecSet_mgf(getProjPath(ip, "./spectra/stars.mgf").c_str());
     DEBUG_TRACE;

     testSpectra.loadPklBin( getProjPath(ip, "./spectra/specs_ms.pklbin").c_str() );
     DEBUG_VAR(testSpectra.size());
     testSpectra.SaveSpecSet_mgf(getProjPath(ip, "./spectra/spec_ms2.mgf").c_str());
     DEBUG_TRACE;

     testSpectra.loadPklBin( getProjPath(ip, "./assembly/sps_seqs.pklbin").c_str() );
     DEBUG_VAR(testSpectra.size());
     testSpectra.SaveSpecSet_mgf(getProjPath(ip, "./assembly/sps_seqs.mgf").c_str());
     DEBUG_TRACE;

     testSpectra.loadPklBin( getProjPath(ip, "./homology/contigs_midx_all.pklbin").c_str() );
     DEBUG_VAR(testSpectra.size());
     testSpectra.SaveSpecSet_mgf(getProjPath(ip, "./homology/contigs_midx_all.mgf").c_str());
     DEBUG_TRACE;

     testSpectra.loadPklBin( getProjPath(ip, "./spectra/contigs.pklbin").c_str() );
     DEBUG_VAR(testSpectra.size());
     testSpectra.SaveSpecSet_mgf(getProjPath(ip, "./spectra/contigs.mgf").c_str());
     DEBUG_TRACE;

     testSpectra.loadPklBin( getProjPath(ip, "./homology/homglue_matches.pklbin").c_str() );
     DEBUG_VAR(testSpectra.size());
     testSpectra.SaveSpecSet_mgf(getProjPath(ip, "./homology/homglue_matches.mgf").c_str());
     DEBUG_TRACE;
     abinfo_t abinfo;
     Load_abinfo("./assembly/component_info.bin", abinfo);
     dumpAbInfo("./assembly/genoms_component_info.txt", abinfo);
     Save_abinfo_v1_0("./assembly/component_info.sps.bin",abinfo);
     */
  }

  if (runMergedFlag && initialStage <= STAGE_GENOMS)
  {

    DEBUG_MSG("Stage Genoms: " << STAGE_GENOMS);

    //We need to rename all of our results files
    if (rename("./spectra/stars.pklbin", "./spectra/genoms.stars.pklbin") < 0)
    {
      ERROR_MSG("Problem encountered renaming GenoMS files");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_GENOMS);
    }
    if (rename("./assembly/sps_seqs.pklbin",
               "./assembly/genoms.sps_seqs.pklbin") < 0)
    {
      ERROR_MSG("Problem encountered renaming GenoMS files");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_GENOMS);
    }
    if (rename("./assembly/component_info.bin",
               "./assembly/genoms.component_info.bin") < 0)
    {
      ERROR_MSG("Problem encountered renaming GenoMS files");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_GENOMS);
    }
    if (rename("./homology/contigs_mp.bin", "./homology/genoms.contigs_mp.bin")
        < 0)
    {
      ERROR_MSG("Problem encountered renaming GenoMS files");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_GENOMS);
    }
    if (rename("./homology/contigs_midx.pklbin",
               "./homology/genoms.contigs_midx.pklbin") < 0)
    {
      ERROR_MSG("Problem encountered renaming GenoMS files");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_GENOMS);
    }
    if (rename("./homology/contigs_mp_all.bin",
               "./homology/genoms.contigs_mp_all.bin") < 0)
    {
      ERROR_MSG("Problem encountered renaming GenoMS files");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_GENOMS);
    }
    if (rename("./homology/contigs_midx_all.pklbin",
               "./homology/genoms.contigs_midx_all.pklbin") < 0)
    {
      ERROR_MSG("Problem encountered renaming GenoMS files");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_GENOMS);
    }
    if (rename("./spectra/contigs.pklbin", "./spectra/genoms.contigs.pklbin")
        < 0)
    {
      ERROR_MSG("Problem encountered renaming GenoMS files");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_GENOMS);
    }
    if (rename("./spectra/contigs_indices.bin",
               "./spectra/genoms.contigs_indices.bin") < 0)
    {
      ERROR_MSG("Problem encountered renaming GenoMS files");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_GENOMS);
    }
    if (rename("./homology/homglue_ref_mp.bin",
               "./homology/genoms.homglue_ref_mp.bin") < 0)
    {
      ERROR_MSG("Problem encountered renaming GenoMS files");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_GENOMS);
    }
    if (rename("./homology/homglue_ref_midx.pklbin",
               "./homology/genoms.homglue_ref_midx.pklbin") < 0)
    {
      ERROR_MSG("Problem encountered renaming GenoMS files");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_GENOMS);
    }
    if (rename("./homology/homglue_matches.pklbin",
               "./homology/genoms.homglue_matches.pklbin") < 0)
    {
      ERROR_MSG("Problem encountered renaming GenoMS files");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_GENOMS);
    }
    if (rename("./homology/homglue_matches_mp.bin",
               "./homology/genoms.homglue_matches_mp.bin") < 0)
    {
      ERROR_MSG("Problem encountered renaming GenoMS files");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_GENOMS);
    }
    if (rename("./homology/homglue_matches_midx.pklbin",
               "./homology/genoms.homglue_matches_midx.pklbin") < 0)
    {
      ERROR_MSG("Problem encountered renaming GenoMS files");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_GENOMS);
    }
    if (rename("./homology/ref_sps_names.txt",
               "./homology/genoms.ref_sps_names.txt") < 0)
    {
      ERROR_MSG("Problem encountered renaming GenoMS files");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_GENOMS);
    }

    if (rename("./protid.fasta", "./genoms.protid.fasta") < 0)
    {
      ERROR_MSG("Problem encountered renaming GenoMS files");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_GENOMS);
    }

  }
  else if (runGenoMSFlag)
  {

    //If we are running GenoMS, then we need to update FASTA_DATABASE to point to our output

    ip.setValue("FASTA_DATABASE", "./protid.fasta");

    if (initialStage <= STAGE_REPORT)
    {
      DEBUG_MSG("Starting Report stage");
      ofstream spsProj(getProjPath(ip, "sps_projects.txt").c_str(),
                       ios_base::out | ios_base::binary);
      spsProj << "sps;.;" << ip.getValue("TOLERANCE_PEAK") << ";"
          << ip.getValue("TOLERANCE_PM") << "\n";
      spsProj.close();

      if (!performReport(ip))
      {
        ERROR_MSG("Problem encountered during Report stage");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_REPORT);
      }
      //--------------------------------------------------------------------------
      //Set up for relaunching
      //--------------------------------------------------------------------------
      if (!generateRelaunchScript(ip))
      {
        ERROR_MSG("Problem encountered during relaunch script creation");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_REPORT);
      }

    }
    else
    {
      DEBUG_MSG("Bypassing Report stage");
    }
    writeStatusFile(statusFileName, "Finished");
    return 0;
  }

  SpectrumPairSet filteredPairs;

  vector<TwoValues<float> > ratios;
  vector<TwoValues<float> > means;
  vector<float> varTerms;
  list<vector<float> > alignStats;
  vector<vector<float> > specStats;
  vector<unsigned int> idxKept;
  vector<TwoValues<float> > pvalues;

  //---------------------------------------------------------------------------
  // FLITERPAIRS STAGE
  //---------------------------------------------------------------------------
  if (initialStage <= STAGE_FILTERPAIRS)
  {
    DEBUG_MSG("Starting stage FilterPairs on: " << getCurrentTimeString());
    if (!performFilterPairs(ip,
                            prmSpectra,
                            ms2Spectra,
                            filteredPairs,
                            ratios,
                            means,
                            varTerms,
                            alignStats,
                            specStats,
                            idxKept,
                            pvalues,
                            gridExecutionFlag,
                            resumeFlag))
    {
      ERROR_MSG("Problem encountered during Filter Pairs stage");
      writeStatusFile(statusFileName, "Error");
      exit(-4);
    }

    if (commandLineParams.exists("PROTEOSAFE_FLAG"))
    {
      bool res = spsSystem("cp aligns/*.params aligns/params/");
      if (!res)
      {
        ERROR_MSG("Failed to copy params files to aligns/params !!!");
        exit(-4);
      }
    }

    if ((!ip.exists("GRID_NUMNODES") || ip.getValueInt("GRID_NUMNODES") > 0)
        && (!resumeFlag && !gridExecutionFlag))
    {
      // If we are doing a grid execution (and are not actually on the grid)
      //    and we aren't resuming... then exit (we'll resume execution later)
      DEBUG_MSG("Files for grid execution have been saved.");
      DEBUG_MSG("Restart with -z option when grid execution has been completed.");
      writeStatusFile(statusFileName, "Finished");
      exit(0);
    }

    if (ip.getValue("PAIRS_MATCH_MODE", "") == "cosine")
    {
      filteredPairs.saveToBinaryFile(getProjPath(ip,
                                                 "./aligns/pairs_cosine.bin"));
    }
    else
    {
      filteredPairs.saveToBinaryFile(getProjPath(ip, "./aligns/pairs_raw.bin"));
    }

    if (commandLineParams.exists("SINGLE_STEP"))
    {
      DEBUG_MSG("Option -s given. Exiting after single step.");
      writeStatusFile(statusFileName, "Finished");
      exit(0);
    }

  }
  else
  {
    DEBUG_MSG("Bypassing Filter Pairs stage");
    filteredPairs.loadFromBinaryFile(getProjPath(ip, "./aligns/pairs_raw.bin"));
    DEBUG_VAR(filteredPairs.size());
    Load_binArray("./aligns/ratios.bin", ratios);
    DEBUG_VAR(ratios.size());
    Load_binArray("./aligns/means.bin", means);
    DEBUG_VAR(means.size());
    Load_binArray("./aligns/vars.bin", varTerms);
    DEBUG_VAR(varTerms.size());
  }

  if (ip.getValue("PAIRS_MATCH_MODE", "") == "cosine")
  {
    DEBUG_MSG("Parameters include PAIRS_MATCH_MODE=cosine, exiting after ExecFilterPairs");
    return 0;
  }

  // Test for empty return data structures
  //  TEST_RETURN("moduleFilterPairs", ratios);
  //  TEST_RETURN("moduleFilterPairs", means);
  //  TEST_RETURN("moduleFilterPairs", varTerms);
  TEST_RETURN("moduleFilterPairs", filteredPairs);
  //  TEST_RETURN("moduleFilterPairs", pvalues);
  //  TEST_RETURN("moduleFilterPairs", idxKept);
  //  TEST_RETURN("moduleFilterPairs", alignStats);
  //  TEST_RETURN("moduleFilterPairs", specStats);

  if (finalStage == STAGE_FILTERPAIRS)
  {
    DEBUG_MSG("Option -f given. Exiting after step: " << finalStageString);
    writeStatusFile(statusFileName, "Finished");
    exit(0);
  }

  //---------------------------------------------------------------------------
  // FLITERALIGNS STAGE
  //---------------------------------------------------------------------------
  if (initialStage <= STAGE_FILTERALIGNS)
  {
    DEBUG_MSG("Starting stage FilterAligns on: " << getCurrentTimeString());
    if (!performFilterAligns(ip,
                             filteredPairs,
                             ratios,
                             means,
                             varTerms,
                             idxKept,
                             pvalues))
    {
      ERROR_MSG("Problem encountered during Filter Aligns stage");
      writeStatusFile(statusFileName, "Error");
      exit(-4);
    }

    if (commandLineParams.exists("SINGLE_STEP"))
    {
      DEBUG_MSG("Option -s given. Exiting after single step.");
      filteredPairs.saveToBinaryFile(getProjPath(ip, "./aligns/pairs.bin"));
      writeStatusFile(statusFileName, "Finished");
      exit(0);
    }
    else
      filteredPairs.saveToBinaryFile(getProjPath(ip, "./aligns/pairs.bin"));
  }
  else
  {
    DEBUG_MSG("Bypassing Filter Aligns stage");
    filteredPairs.loadFromBinaryFile(getProjPath(ip, "./aligns/pairs.bin"));
    DEBUG_VAR(filteredPairs.size());
  }

  // Test for empty return data structures
  TEST_RETURN("moduleFilterPairs", filteredPairs);

  if (finalStage == STAGE_FILTERALIGNS)
  {
    DEBUG_MSG("Option -f given. Exiting after step: " << finalStageString);
    writeStatusFile(statusFileName, "Finished");
    exit(0);
  }

  if (commandLineParams.exists("RESUME_SINGLE_FLAG"))
  {
    ip.setValue("GRID_NUMNODES", "-1");
  }

  SpecSet pairAlignments;
  SpecSet starSpectraOnly;
  SpecSet starSpectra;
  vector<unsigned int> alignedSpectra;

  //---------------------------------------------------------------------------
  // ALIGNMENT STAGE
  //---------------------------------------------------------------------------
  if (initialStage <= STAGE_ALIGNMENT)
  {
    DEBUG_MSG("Starting stage Alignment on: " << getCurrentTimeString());
    if (!performAlignment(ip,
                          jumps,
                          prmSpectra,
                          filteredPairs,
                          pairAlignments,
                          starSpectraOnly,
                          starSpectra,
                          alignedSpectra))
    {
      ERROR_MSG("Problem encountered during Alignment stage");
      writeStatusFile(statusFileName, "Error");
      exit(-5);
    }

    if (commandLineParams.exists("SINGLE_STEP"))
    {
      DEBUG_MSG("Option -s given. Exiting after single step.");
      writeStatusFile(statusFileName, "Finished");
      exit(0);
    }
  }
  else
  {
    DEBUG_MSG("Bypassing Alignment stage");
    starSpectra.loadPklBin(getProjPath(ip, "./spectra/stars.pklbin").c_str());
    DEBUG_VAR(starSpectra.size());
  }

  // Test for empty return data structures
  //TEST_RETURN("ExecAlignment", pairAlignments);
  //TEST_RETURN("ExecAlignment", starSpectraOnly);
  TEST_RETURN("ExecAlignment", starSpectra);
  //TEST_RETURN("ExecAlignment", alignedSpectra);

  if (finalStage == STAGE_ALIGNMENT)
  {
    DEBUG_MSG("Option -f given. Exiting after step: " << finalStageString);
    writeStatusFile(statusFileName, "Finished");
    exit(0);
  }

  vector<vector<float> > starRatios;
  SpecSet matchedPeaks;

  //---------------------------------------------------------------------------
  // FILTERSTARPAIRS STAGE
  //---------------------------------------------------------------------------
  if (initialStage <= STAGE_FILTERSTARPAIRS)
  {
    if (!performFilterStarPairs(ip,
                                filteredPairs,
                                starSpectra,
                                starRatios,
                                matchedPeaks))
    {
      ERROR_MSG("Problem encountered during Filter Star Pairs stage");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_FILTERSTARPAIRS);
    }

    if (commandLineParams.exists("SINGLE_STEP"))
    {
      DEBUG_MSG("Option -s given. Exiting after single step.");
      writeStatusFile(statusFileName, "Finished");
      exit(0);
    }
  }
  else
  {
    DEBUG_MSG("Bypassing Filter Star Pairs stage");
    //    DEBUG_MSG("Not implemented yet" );
    //    Load_binArray("./aligns/ratios.bin", starRatios);
    //    DEBUG_VAR(starRatios.size());
    //    matchedPeaks->loadPklBin("");
    //    DEBUG_VAR(matchedPeaks.size());

    if (!filteredPairs.loadFromBinaryFile("aligns/pairs_stars.bin"))
    {
      ERROR_MSG("Problem encountered during Filter Star Pairs stage (loading filteredPairs)");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_FILTERSTARPAIRS);
    }
  }

  if (finalStage == STAGE_FILTERSTARPAIRS)
  {
    DEBUG_MSG("Option -f given. Exiting after step: " << finalStageString);
    writeStatusFile(statusFileName, "Finished");
    exit(0);
  }

  DEBUG_MSG("Reloading Star specturm");
  // Need to reload the star spectra because they were altered by FilterStarPairs
  if (starSpectra.loadPklBin(getProjPath(ip, "./spectra/stars.pklbin").c_str())
      <= 0)
  {
    ERROR_MSG("Problem encountered while reloading star specturm");
    writeStatusFile(statusFileName, "Error");
    exit(-STAGE_ASSEMBLY);
  }
  DEBUG_VAR(starSpectra.size());

  if (execMode == string("signatures"))
  {
    return 0;
  }

  // Test for empty return data structures
  //TEST_RETURN("ExecFilterStarPairs", starRatios);
  //TEST_RETURN("ExecFilterStarPairs", matchedPeaks);

  //---------------------------------------------------------------------------
  // Spectral Networks STAGE
  //---------------------------------------------------------------------------
  DB_fasta dbAll;
  string dbFileName;
  DEBUG_MSG(ip.getValue("FASTA_DATABASE"));
  if (ip.exists("FASTA_DATABASE"))
  {

    dbFileName = ip.getValue("FASTA_DATABASE");

    //If we are integrating, then add GenoMS protein sequences to the DB
    if (runMergedFlag && initialStage <= STAGE_FILTERSTARPAIRS)
    {
      DEBUG_MSG("Concatenating " << dbFileName << " and genoms.protid.fasta");
      if (!concatenateFiles("./genoms.protid.fasta",
                            dbFileName,
                            "./protid.fasta"))
      {
        ERROR_MSG("Problem encountered concatenating fasta files");
        writeStatusFile(statusFileName, "Error");
        exit(-100);
      }
      dbFileName = "./protid.fasta";
      ip.setValue("FASTA_DATABASE", "./protid.fasta");
    }

    if (dbAll.Load(dbFileName.c_str()) <= 0)
    {
      ERROR_MSG("Problem encountered during Spectral Networks stage (loading " << dbFileName << ")");
      writeStatusFile(statusFileName, "Error");
      exit(-100);
    }
  }

  //---------------------------------------------------------------------------
  // Create penalty matrices for new alignment
  //---------------------------------------------------------------------------
  float resolution = ip.getValueDouble("ALIGNMENT_RESOLUTION", 0.1);
  float peakEquivalents = ip.getValueDouble("PEAK_EQUIVALENTS", 2.0);
  float minFrequency = ip.getValueDouble("MIN_PENALTY_FREQUENCY", 0.01);
  float unknownPenalty = ip.getValueDouble("PENALTY_ALIGNMENT_UNKNOWN_PENALTY",
                                           1.0);
  float unknownMultiplier =
      ip.getValueDouble("PENALTY_ALIGNMENT_UNKNOWN_MULTIPLIER", 2.0);
  DEBUG_VAR(resolution);
  DEBUG_VAR(peakEquivalents);
  PenaltyMatrix penaltyMatrixBlosum(jumps,
                                    resolution,
                                    unknownPenalty,
                                    unknownMultiplier);
  PenaltyMatrix penaltyMatrixMods(jumps,
                                  resolution,
                                  unknownPenalty,
                                  unknownMultiplier);

  bool penaltyAlign = (bool)ip.getValueInt("PENALTY_ALIGNMENT");
  if (penaltyAlign)
  {
    vector<pair<char, float> > vecRefAminoAcidsTest;
    if (jumps.getAllAArefs(vecRefAminoAcidsTest) == 0)
    {
      ERROR_MSG("Penalty alignment specified but no amino acid masses found.");
      exit(-1);
    }

    if (!ip.exists("BLOSUM_PENALTY_FILE"))
    {
      ERROR_MSG("Penalty alignment specified but BLOSUM_PENALTY_FILE not set");
      exit(-1);
    }

    string blosumFilename = ip.getValue("BLOSUM_PENALTY_FILE");
    DEBUG_VAR(blosumFilename);
    if (!penaltyMatrixBlosum.loadFromBlosum(blosumFilename, peakEquivalents))
    {
      ERROR_MSG("Unable to load BLOSUM_PENALTY_FILE");
      exit(-1);
    }

    if (!ip.exists("KNOWN_MODS_FILE"))
    {
      WARN_MSG("Penalty alignment specified but KNOWN_MODS_FILE not set");
    }
    string blosumMatFilename("homology/specprotalign_blosum_penal.txt");
    DEBUG_VAR(blosumMatFilename);
    penaltyMatrixBlosum.saveMatrix(blosumMatFilename);

    DEBUG_TRACE;
    map<float, float> modFreqs;
    filteredPairs.getModificationFrequencies(resolution, modFreqs);
    DEBUG_TRACE;
    float avgIntensity = prmSpectra.averageIntensity();
    DEBUG_VAR(avgIntensity);
    string knownmodFilename = ip.getValue("KNOWN_MODS_FILE");
    DEBUG_VAR(knownmodFilename);
    penaltyMatrixMods.loadKnownModifications(knownmodFilename);
    penaltyMatrixMods.createFromModificationFreqs(modFreqs,
                                                  peakEquivalents,
                                                  minFrequency,
                                                  avgIntensity);

    string matrixFilename("homology/specprotalign_mod_penal.txt");
    string modknownFilename("homology/specprotalign_mod_known.txt");
    penaltyMatrixMods.saveMatrix(matrixFilename);
    penaltyMatrixMods.saveKnownMods(modknownFilename);
  }

  if (execMode == string("specnets"))
  {
    // Entering Spectral Networks mode
    PeptideSpectrumMatchSet * input_psms = new PeptideSpectrumMatchSet;
    PeptideSpectrumMatchSet * output_psms = new PeptideSpectrumMatchSet;
    SpecSet * psms_spectra = new SpecSet;
    SpecSet * psms_midx = new SpecSet;
    vector<vector<int> > *psms_mp = new vector<vector<int> >;
    SpecSet * snets_contigs = new SpecSet;
    SpecSet * snets_midx = new SpecSet;
    float resolution = ip.getValueFloat("@", .1);

    vector<vector<int> > *snets_mp = new vector<vector<int> >;

    if (ip.exists("INSPECT_PSMS"))
    {
      //load in inspect seed PSMS
      if (!input_psms->loadInspectResultsFile(ip.getValue("INSPECT_PSMS").c_str(),
                                              ip.getValueBool("SCAN_ZERO_INDEX",
                                                              1)))
      {
        ERROR_MSG("Unable to load inspect file! " << ip.getValue("INSPECT_PSMS"));
        writeStatusFile(statusFileName, "Error");
        exit(-100);
      }
    }

    if (not performExecMainSpecnets(ip,
                                    &ms2Spectra,
                                    &prmSpectra,
                                    &starSpectra,
                                    &filteredPairs,
                                    &dbAll,
                                    &penaltyMatrixMods,
                                    &penaltyMatrixBlosum,
                                    output_psms,
                                    input_psms,
                                    psms_spectra,
                                    psms_midx,
                                    psms_mp,
                                    snets_contigs,
                                    snets_midx,
                                    snets_mp))
    {
      ERROR_MSG("Problem encountered during ExecMainSpecnets stage");
      writeStatusFile(statusFileName, "Error");
      exit(-100);
    }

    delete input_psms;
    delete output_psms;
    delete psms_spectra;
    delete psms_midx;
    delete psms_mp;
    delete snets_contigs;
    delete snets_midx;
    delete snets_mp;

    return 0;
  }

  Clusters contigShifts;
  abinfo_t contigAbinfo;
  //---------------------------------------------------------------------------
  // ASSEMBLY STAGE
  //---------------------------------------------------------------------------
  if (initialStage <= STAGE_ASSEMBLY and execMode == string("sps"))
  {
    if (!performAssembly(ip,
                         starSpectra,
                         filteredPairs,
                         contigShifts,
                         contigAbinfo))
    {
      ERROR_MSG("Problem encountered during Assembly stage");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_ASSEMBLY);
    }

    if (commandLineParams.exists("SINGLE_STEP"))
    {
      DEBUG_MSG("Option -s given. Exiting after single step.");
      writeStatusFile(statusFileName, "Finished");
      exit(0);
    }
    contigShifts.Save("assembly/path_spectra_as_cluster.txt");

    bool res = FileCopy("assembly/sps_seqs.pklbin",
                        "assembly/old_sps_seqs.pklbin");
    //spsSystem("cp assembly/sps_seqs.pklbin assembly/old_sps_seqs.pklbin");
    res = FileCopy("assembly/path_spectra_as_cluster.txt",
                   "assembly/old_path_spectra_as_cluster.txt");
    //= spsSystem("cp assembly/path_spectra_as_cluster.txt assembly/old_path_spectra_as_cluster.txt");
    res = FileCopy("assembly/component_info.bin",
                   "assembly/old_component_info.bin");
    //= spsSystem("cp assembly/component_info.bin assembly/old_component_info.bin");
  }
  else
  {
    DEBUG_MSG("Bypassing Assembly stage");
    if (contigShifts.Load("assembly/old_path_spectra_as_cluster.txt") <= 0
        and contigShifts.Load("assembly/path_spectra_as_cluster.txt") <= 0)
    {
      ERROR_MSG("Problem encountered while skipping Assembly stage");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_ASSEMBLY);
    }
    if (contigShifts.consensus.loadPklBin("assembly/old_sps_seqs.pklbin") <= 0
        and contigShifts.consensus.loadPklBin("assembly/sps_seqs.pklbin") <= 0)
    {
      ERROR_MSG("Problem encountered while skipping Assembly stage");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_ASSEMBLY);
    }
    DEBUG_MSG("Loading contigAbinfo...");
    if (!Load_abinfo("assembly/component_info.bin", contigAbinfo))
    {
      ERROR_MSG("Problem encountered while skipping Assembly stage");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_ASSEMBLY);
    }
  }

  // Test for empty return data structures
  TEST_RETURN("ExecAssembly", contigShifts);

  if (finalStage == STAGE_ASSEMBLY)
  {
    DEBUG_MSG("Option -f given. Exiting after step: " << finalStageString);
    writeStatusFile(statusFileName, "Finished");
    exit(0);
  }

  if (initialStage <= STAGE_METAASSEMBLY
      and ip.getValueInt("MIN_METACONTIG_SIZE", 0) > 0)
  {
    SpectrumPairSet contigPairs;
    Clusters metaContigs;
    abinfo_t metaContigAbinfo;

    if (!performContigAlignment(ip, contigShifts, contigPairs))
    {
      ERROR_MSG("Problem encountered during MetaAssembly alignment stage");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_METAASSEMBLY);
    }

    if (!performMetaAssembly(ip,
                             contigShifts,
                             contigPairs,
                             contigAbinfo,
                             starSpectra,
                             metaContigs,
                             metaContigAbinfo))
    {
      ERROR_MSG("Problem encountered during MetaAssembly stage");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_METAASSEMBLY);
    }

    if (commandLineParams.exists("SINGLE_STEP"))
    {
      DEBUG_MSG("Option -s given. Exiting after single step.");
      writeStatusFile(statusFileName, "Finished");
      exit(0);
    }

    contigShifts = metaContigs;
    contigAbinfo = metaContigAbinfo;
  }
  else if (ip.getValueInt("MIN_METACONTIG_SIZE", 0) > 0)
  {
    DEBUG_MSG("Bypassing MetaAssembly stage");
    if (contigShifts.Load("assembly/path_spectra_as_cluster.txt") <= 0)
    {
      ERROR_MSG("Problem encountered while skipping MetaAssembly stage");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_METAASSEMBLY);
    }
    if (!Load_abinfo("assembly/component_info.bin", contigAbinfo))
    {
      ERROR_MSG("Problem encountered while skipping MetaAssembly stage");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_METAASSEMBLY);
    }
  }

  if (finalStage == STAGE_METAASSEMBLY)
  {
    DEBUG_MSG("Option -f given. Exiting after step: " << finalStageString);
    writeStatusFile(statusFileName, "Finished");
    exit(0);
  }

  DEBUG_VAR(contigShifts.size());

  TEST_RETURN("ExecMetaAssembly", contigShifts);

  DB_fasta dbDecoys;
  if (ip.getValueFloat("FDR_CONTIG_THRESHOLD", -1.0) >= 0.0
      || ip.getValueFloat("FDR_SPECTRUM_THRESHOLD", -1.0) >= 0.0)
  {
    if (!ip.exists("DECOY_DATABASE"))
    {
      ERROR_MSG("FDR specfied for spectrum alignment but DECOY_DATABASE not set");
      exit(0);
    }
    string dbDecoyFileName = ip.getValue("DECOY_DATABASE");
    if (dbDecoys.Load(dbDecoyFileName.c_str()) <= 0)
    {
      ERROR_MSG("Unable to load DECOY_DATABASE [" << dbDecoyFileName << "]");
      exit(0);
    }
  }

  PeptideSpectrumMatchSet psmSetTag;
  PeptideSpectrumMatchSet psmSetTagDecoy;

  if (ip.exists("FASTA_DATABASE"))
  {

    //---------------------------------------------------------------------------
    // TAGSEARCH STAGE
    //---------------------------------------------------------------------------
    if (initialStage <= STAGE_TAGSEARCH)
    {
      bool ok;
      if (execMode == string("sps"))
      {
        ok = performTagsearch(ip,
                              contigShifts.consensus,
                              dbAll,
                              (vector<unsigned int> *)0, //specsToSearch,
                              psmSetTag,
                              false);
      }
      else
      {
        ok = performTagsearch(ip,
                              contigShifts.consensus,
                              dbAll,
                              &alignedSpectra,
                              psmSetTag,
                              false);
      }
      if (!ok)
      {
        ERROR_MSG("Problem encountered during TagSearch stage");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_TAGSEARCH);
      }

      if (dbDecoys.size() != 0)
      {
        if (!performTagsearch(ip,
                              contigShifts.consensus,
                              dbDecoys,
                              (vector<unsigned int> *)0, //specsToSearch,
                              psmSetTagDecoy,
                              true))
        {
          ERROR_MSG("Problem encountered during TagSearch stage");
          writeStatusFile(statusFileName, "Error");
          exit(-STAGE_TAGSEARCH);
        }
      }

      if (commandLineParams.exists("SINGLE_STEP"))
      {
        DEBUG_MSG("Option -s given. Exiting after single step.");
        writeStatusFile(statusFileName, "Finished");
        exit(0);
      }
    }
    else
    {
      // No need to load data here.. since we always have to reload after TagSearch
    }

    DEBUG_MSG("Reloading contig spectrum");
    // Clear the spectra so we don't double them (or the PSMs) on reload
    contigShifts.consensus.clear();
    // Need to reload the contigs because they were altered by TagSearch
    if (contigShifts.Load(getProjPath(ip,
                                      "assembly/path_spectra_as_cluster.txt").c_str(),
                          getProjPath(ip, "assembly/sps_seqs.pklbin").c_str(),
                          getProjPath(ip, "assembly/tagsearchpsm.txt").c_str())
        <= 0)
    {
      ERROR_MSG("Problem reloading contigs after TagSearch stage");
      writeStatusFile(statusFileName, "Error");
      exit(-1);
    }

    if (ip.getValueFloat("FDR_CONTIG_THRESHOLD", -1.0) >= 0.0
        || ip.getValueFloat("FDR_SPECTRUM_THRESHOLD", -1.0) >= 0.0)
    {
      psmSetTagDecoy.loadFromFile(getProjPath(ip,
                                              "assembly/tagsearchpsm_decoy.txt").c_str());
    }

    // Test for empty return data structures
    //TEST_RETURN("ExecTagSearch", alignedSpectra);
    DEBUG_VAR(contigShifts.size());

    if (finalStage == STAGE_TAGSEARCH)
    {
      DEBUG_MSG("Option -f given. Exiting after step: " << finalStageString);
      writeStatusFile(statusFileName, "Finished");
      exit(0);
    }

    DEBUG_VAR(contigShifts.size());

    //---------------------------------------------------------------------------
    // Contig/Protein Alignment STAGE
    //---------------------------------------------------------------------------
    SpecSet matchedContigsAll; // Matched versions of spectra significantly matched to a protein (All)
    SpecSet matchedContigs; // Matched versions of spectra significantly matched to a protein

    PeptideSpectrumMatchSet psmSetContig;
    PeptideSpectrumMatchSet psmSetContigDecoy;
    PeptideSpectrumMatchSet psmSetContigFdr;

    if (initialStage <= STAGE_CONTIGPROTALIGN)
    {
      DEBUG_VAR(contigShifts.size());
      if (!performContigProtAlign(ip,
                                  contigShifts.consensus,
                                  dbAll,
                                  dbDecoys,
                                  psmSetTag,
                                  psmSetTagDecoy,
                                  penaltyMatrixBlosum,
                                  penaltyMatrixMods,
                                  matchedContigsAll,
                                  matchedContigs,
                                  psmSetContig,
                                  psmSetContigDecoy,
                                  psmSetContigFdr,
                                  gridExecutionFlag,
                                  resumeFlag))
      {
        ERROR_MSG("Problem encountered during ContigProtAlign stage");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_CONTIGPROTALIGN);
      }

      if ((!ip.exists("GRID_NUMNODES") || ip.getValueInt("GRID_NUMNODES") > 0)
          && (!resumeFlag && !gridExecutionFlag))
      {
        // If we are doing a grid execution (and are not actually on the grid)
        //    and we aren't resuming... then exit (we'll resume execution later)
        DEBUG_MSG("Files for grid execution have been saved.");
        DEBUG_MSG("Restart with -z option when grid execution has been completed.");
        writeStatusFile(statusFileName, "Finished");
        exit(0);
      }

      if (commandLineParams.exists("SINGLE_STEP"))
      {
        DEBUG_MSG("Option -s given. Exiting after single step.");
        writeStatusFile(statusFileName, "Finished");
        exit(0);
      }
    }
    else
    {
      DEBUG_MSG("Bypassing Contig/Protein Alignment stage");
      if (matchedContigs.loadPklBin(getProjPath(ip, "spectra/contigs.pklbin").c_str(),
                                    getProjPath(ip, "homology/contigs_psm.txt").c_str(),
                                    getProjPath(ip,
                                                "homology/contigs_midx.pklbin").c_str())
          <= 0)
      {
        ERROR_MSG("Problem encountered while skipping ContigProtAlign stage");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_CONTIGPROTALIGN);
      }

      if (matchedContigsAll.loadPklBin(getProjPath(ip,
                                                   "homology/contigs_tgt.pklbin").c_str())
          <= 0)
      {
        ERROR_MSG("Problem encountered while skipping ContigProtAlign stage");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_CONTIGPROTALIGN);
      }

      if (!psmSetContig.loadFromFile(getProjPath(ip, "homology/contigs_psm.txt").c_str()))
      {
        ERROR_MSG("Problem encountered while skipping ContigProtAlign stage");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_CONTIGPROTALIGN);
      }

      if (ip.getValueFloat("FDR_CONTIG_THRESHOLD", -1.0) >= 0.0)
      {
        if (!psmSetContigDecoy.loadFromFile(getProjPath(ip,
                                                        "homology/contigs_psm_dec.txt").c_str()))
        {
          ERROR_MSG("Problem encountered while skipping ContigProtAlign stage");
          writeStatusFile(statusFileName, "Error");
          exit(-STAGE_CONTIGPROTALIGN);
        }
      }
    }

    //SpecSet matchedContigs; // Matched versions of spectra significantly matched to a protein
    //SpecSet matchedMassIndices; // Per-spectrum indices of matched spectrum/protein masses

    // Test for empty return data structures
    //TEST_RETURN("ExecContigProtAlign", matchedContigs);
    DEBUG_VAR(contigShifts.size());
    DEBUG_VAR(matchedContigs.size());
    DEBUG_VAR(matchedContigsAll.size());

    if (finalStage == STAGE_CONTIGPROTALIGN)
    {
      DEBUG_MSG("Option -f given. Exiting after step: " << finalStageString);
      writeStatusFile(statusFileName, "Finished");
      exit(0);
    }

    //---------------------------------------------------------------------------
    // Spectrum/Protein Alignment STAGE
    //---------------------------------------------------------------------------
    if (initialStage <= STAGE_SPECPROTALIGN
        && ip.getValueInt("ALIGN_STARS", 0) > 0)
    {
      SpecSet starAlignSpectra(starSpectra);
      SpecSet starAlignSpectraDecoy(starSpectra);
      for (int i = 0; i < starAlignSpectra.size(); i++)
      {
        starAlignSpectra[i] = starSpectra[i];
        starAlignSpectraDecoy[i] = starSpectra[i];
      }

      // Get the starting masses of the each spectrum in each contig
      vector<vector<sps::tuple<unsigned int, float, bool> > > outputAssembledShifts;
      DEBUG_VAR(contigAbinfo.size());
      getAssembledShifts(contigShifts.consensus,
                         contigAbinfo,
                         outputAssembledShifts);
      DEBUG_VAR(outputAssembledShifts.size());
      DEBUG_VAR(starSpectra.size());

      PeptideSpectrumMatchSet psmSetContigTag;
      PeptideSpectrumMatchSet psmSetContigTagDecoy;

      psmSetContig.addSpectra(&matchedContigsAll); // Associate target PSMs

      makeSpectrumTagsFromContig(matchedContigsAll,
                                 outputAssembledShifts,
                                 starSpectra,
                                 starAlignSpectra,
                                 psmSetContigTag);
      DEBUG_VAR(psmSetContigTag.size());

      matchedContigsAll.clearPsms(); // Clear all the target PSMs
      psmSetContigDecoy.addSpectra(&matchedContigsAll); // Associate decoy PSMs

      makeSpectrumTagsFromContig(matchedContigsAll,
                                 outputAssembledShifts,
                                 starSpectra,
                                 starAlignSpectraDecoy,
                                 psmSetContigTagDecoy);
      DEBUG_VAR(psmSetContigTagDecoy.size());

      matchedContigsAll.clearPsms(); // Clear all the decoy PSMs
      psmSetContig.addSpectra(&matchedContigsAll); // Associate target PSMs again (dumb I know.. fix later)

      SpecSet starPrmAll; // Matched versions of spectra significantly matched to a protein (All)
      SpecSet starPrm; // Matched versions of spectra significantly matched to a protein
      PeptideSpectrumMatchSet psmSetSpectra;
      PeptideSpectrumMatchSet psmSetSpectraDecoy;
      PeptideSpectrumMatchSet psmSetSpectraFdr;

      DEBUG_TRACE;
      if (!performSpecProtAlign(ip, starAlignSpectra, prmSpectra, // PRM spectra for "sprinkling" into star spectra
                                dbAll,
                                dbDecoys,
                                psmSetContigTag,
                                psmSetContigTagDecoy,
                                penaltyMatrixBlosum,
                                penaltyMatrixMods,
                                starPrmAll,
                                starPrm,
                                psmSetSpectra,
                                psmSetSpectraDecoy,
                                psmSetSpectraFdr,
                                gridExecutionFlag,
                                resumeFlag))
      {
        ERROR_MSG("Problem encountered during SpecProtAlign stage");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_SPECPROTALIGN);
      }

      if ((!ip.exists("GRID_NUMNODES") || ip.getValueInt("GRID_NUMNODES") > 0)
          && (!resumeFlag && !gridExecutionFlag))
      {
        // If we are doing a grid execution (and are not actually on the grid)
        //    and we aren't resuming... then exit (we'll resume execution later)
        DEBUG_MSG("Files for grid execution have been saved.");
        DEBUG_MSG("Restart with -z option when grid execution has been completed.");
        writeStatusFile(statusFileName, "Finished");
        exit(0);
      }

      if (commandLineParams.exists("SINGLE_STEP"))
      {
        DEBUG_MSG("Option -s given. Exiting after single step.");
        writeStatusFile(statusFileName, "Finished");
        exit(0);
      }
    }
    else
    {
      DEBUG_MSG("Bypassing Spectrum/Protein Alignment stage");
    }

    if (finalStage == STAGE_SPECPROTALIGN)
    {
      DEBUG_MSG("Option -f given. Exiting after step: " << finalStageString);
      writeStatusFile(statusFileName, "Finished");
      exit(0);
    }

    //---------------------------------------------------------------------------
    // Find reference protein for protein/protein sequence alignments
    //---------------------------------------------------------------------------
    unsigned int refProtIdx = 0, refProtMatchCount = 0;
    set<unsigned int> dbIndexes;

    //---------------------------------------------------------------------------
    // Assign the reference protein for protein/protein sequence alignments
    //---------------------------------------------------------------------------
    if (ip.exists("FORCE_REFERENCE"))
    {
      DEBUG_MSG("POINT A-prime");
      // Force use of the supplied reference
      refProtIdx = ip.getValueInt("FORCE_REFERENCE");
      DEBUG_MSG("Point B-prime");
      for (unsigned int i = 0; i < matchedContigs.size(); i++)
      {
        DEBUG_MSG("Point C-prime");
        dbIndexes.insert(matchedContigs[i].psmList.front()->m_dbIndex);
      }
      DEBUG_VAR(refProtIdx);

    }
    else
    {
      // Use the most common protein match as the reference
      DEBUG_MSG("Point A.a");
      vector<unsigned int> numMatches(dbAll.size());
      DEBUG_MSG("Point A.b");
      for (unsigned int i = 0; i < dbAll.size(); i++)
        numMatches[i] = 0;
      DEBUG_MSG("Point A");
      for (unsigned int i = 0; i < matchedContigs.size(); i++)
      {
        DEBUG_VAR(matchedContigs.size());
        DEBUG_VAR(i);
        DEBUG_VAR(matchedContigs[i].psmList.size());
        if (matchedContigs[i].psmList.size() > 0)
        {
          DEBUG_VAR(matchedContigs[i].psmList.front());
          int index = matchedContigs[i].psmList.front()->m_dbIndex;
          DEBUG_MSG("Point B");
          dbIndexes.insert(index);
          DEBUG_MSG("Point C");
          if (index >= 0 and index < dbAll.size())
          {
            DEBUG_MSG("Point D");
            //    numMatches[index]++;  // Select by highest number of matched contigs
            DEBUG_VAR(numMatches.size());
            DEBUG_VAR(index);

            DEBUG_VAR(matchedContigs[i].psmList.front()->m_matchedPeaks.size());
            numMatches[index] +=
                matchedContigs[i].psmList.front()->m_matchedPeaks.size(); // Select by highest number of matched masses
            if (numMatches[index] > refProtMatchCount)
            {
              DEBUG_MSG("Point E");
              refProtIdx = index;
              refProtMatchCount = numMatches[index];
            }
          }

          DEBUG_VAR(refProtIdx);
          DEBUG_VAR(refProtMatchCount);
        }
        //Since we didn't write out the PMS matches for these contigs (they are GenoMS contigs)
        else if (matchedContigs[i].psmList.size() == 0
            && (runGenoMSFlag || runMergedFlag))
        {

          int index = 0;
          DEBUG_MSG("Point B_NEC");
          dbIndexes.insert(index);
          DEBUG_MSG("Point C_NEC");
          if (index >= 0 and index < dbAll.size())
          {
            DEBUG_MSG("Point D_NEC");
            //    numMatches[index]++;  // Select by highest number of matched contigs
            DEBUG_VAR(numMatches.size());
            DEBUG_VAR(index);

            //DEBUG_VAR(matchedContigs[i].psmList.front()->m_matchedPeaks.size());
            numMatches[index] += matchedContigs[i].size(); // Select by highest number of matched masses
            if (numMatches[index] > refProtMatchCount)
            {
              DEBUG_MSG("Point E_NEC");
              refProtIdx = index;
              refProtMatchCount = numMatches[index];
            }
          }

          DEBUG_VAR(refProtIdx);
          DEBUG_VAR(refProtMatchCount);

        }
      }
    } // if (ip.exists("FORCE_REFERENCE")) {

    //---------------------------------------------------------------------------
    // Protein/Protein Alignment STAGE
    //---------------------------------------------------------------------------
    DEBUG_TRACE;
    string indexFileName; // File names of clustalw .aln output files
    if (initialStage <= STAGE_PROTPROTALIGN)
    {
      if (!performProtProtAlign(ip, refProtIdx, dbIndexes, dbAll))
      {
        ERROR_MSG("Problem encountered during ProtProtAlign stage");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_PROTPROTALIGN);
      }

      if (commandLineParams.exists("SINGLE_STEP"))
      {
        DEBUG_MSG("Option -s given. Exiting after single step.");
        writeStatusFile(statusFileName, "Finished");
        exit(0);
      }
    }
    else
    {
      DEBUG_MSG("Bypassing Protein/Protein Alignment stage");
    }

    DEBUG_VAR(dbAll.size());
    DEBUG_VAR(matchedContigs.size());

    //---------------------------------------------------------------------------
    // Homology-based Assembly STAGE
    //---------------------------------------------------------------------------
    SpecSet cspsContigs;
    SpecSet cspsMatchedIndices;
    vector<vector<int> > cspsMatchedProts;
    if (initialStage <= STAGE_HOMOLOGYASSEMBLY)
    {
      DEBUG_TRACE;
      if (!performHomologyAssembly(ip,
              matchedContigs,
              dbAll,
              cspsContigs,
              cspsMatchedIndices,
              cspsMatchedProts))
      {
        ERROR_MSG("Problem encountered during cSPS HomologyAssembly stage");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_HOMOLOGYASSEMBLY);
      }

      if (commandLineParams.exists("SINGLE_STEP"))
      {
        DEBUG_MSG("Option -s given. Exiting after single step.");
        writeStatusFile(statusFileName, "Finished");
        exit(0);
      }
    }
    else
    {
      DEBUG_MSG("Bypassing Homology Assembly stage");
      if (cspsContigs.loadPklBin(getProjPath(ip,
                  "homology/homglue_matches.pklbin").c_str())
          <= 0
          or cspsMatchedIndices.loadPklBin(getProjPath(ip,
                  "homology/homglue_matches_midx.pklbin").c_str())
          <= 0
          or Load_binArray<int> (getProjPath(ip,
                  "homology/homglue_matches_mp.bin").c_str(),
              cspsMatchedProts) <= 0)
      {
        ERROR_MSG("Problem encountered while skipping SpecProtAlign stage");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_HOMOLOGYASSEMBLY);
      }
    }

    // Test for empty return data structures
    TEST_RETURN("ExecHomologyAssembly", cspsContigs);
    TEST_RETURN("ExecHomologyAssembly", cspsMatchedIndices);
    TEST_RETURN("ExecHomologyAssembly", cspsMatchedProts);

    //---------------------------------------------------------------------------
    // MERGE RESULTS IF RUNNING IN INTEGRATIVE MODE
    //---------------------------------------------------------------------------
    if (runMergedFlag && initialStage <= STAGE_HOMOLOGYASSEMBLY)
    {
      //We need to rename all of our results files
      if (rename("./spectra/stars.pklbin", "./spectra/csps.stars.pklbin") < 0)
      {
        ERROR_MSG("Problem encountered renaming CSPS files");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_MERGE);
      }
      if (rename("./assembly/sps_seqs.pklbin",
                 "./assembly/csps.sps_seqs.pklbin") < 0)
      {
        ERROR_MSG("Problem encountered renaming CSPS files");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_MERGE);
      }
      if (rename("./assembly/component_info.bin",
                 "./assembly/csps.component_info.bin") < 0)
      {
        ERROR_MSG("Problem encountered renaming CSPS files");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_MERGE);
      }
      if (rename("./homology/contigs_mp.bin", "./homology/csps.contigs_mp.bin")
          < 0)
      {
        ERROR_MSG("Problem encountered renaming CSPS files");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_MERGE);
      }
      if (rename("./homology/contigs_midx.pklbin",
                 "./homology/csps.contigs_midx.pklbin") < 0)
      {
        ERROR_MSG("Problem encountered renaming CSPS files");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_MERGE);
      }
      if (rename("./homology/contigs_mp_all.bin",
                 "./homology/csps.contigs_mp_all.bin") < 0)
      {
        ERROR_MSG("Problem encountered renaming CSPS files");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_MERGE);
      }

      /*if (rename("./homology/contigs_midx_all.pklbin",
       "./homology/csps.contigs_midx_all.pklbin") < 0)
       {
       ERROR_MSG("Problem encountered renaming CSPS files");
       writeStatusFile(statusFileName, "Error");
       exit(-STAGE_MERGE);
       }*/

      if (rename("./spectra/contigs.pklbin", "./spectra/csps.contigs.pklbin")
          < 0)
      {
        ERROR_MSG("Problem encountered renaming CSPS files");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_MERGE);
      }

      if (rename("./spectra/contigs_indices.bin",
                 "./spectra/csps.contigs_indices.bin") < 0)
      {
        ERROR_MSG("Problem encountered renaming CSPS files");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_MERGE);
      }

      if (rename("./homology/homglue_ref_mp.bin",
                 "./homology/csps.homglue_ref_mp.bin") < 0)
      {
        ERROR_MSG("Problem encountered renaming CSPS files");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_MERGE);
      }

      if (rename("./homology/homglue_ref_midx.pklbin",
                 "./homology/csps.homglue_ref_midx.pklbin") < 0)
      {
        ERROR_MSG("Problem encountered renaming CSPS files");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_MERGE);
      }

      if (rename("./homology/homglue_matches.pklbin",
                 "./homology/csps.homglue_matches.pklbin") < 0)
      {
        ERROR_MSG("Problem encountered renaming CSPS files");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_MERGE);
      }

      if (rename("./homology/homglue_matches_mp.bin",
                 "./homology/csps.homglue_matches_mp.bin") < 0)
      {
        ERROR_MSG("Problem encountered renaming CSPS files");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_MERGE);
      }

      if (rename("./homology/homglue_matches_midx.pklbin",
                 "./homology/csps.homglue_matches_midx.pklbin") < 0)
      {
        ERROR_MSG("Problem encountered renaming CSPS files");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_MERGE);
      }

      ifstream ifile("homology/ref_sps_names.txt");
      if (ifile)
      {

        if (rename("./homology/ref_sps_names.txt",
                   "./homology/csps.ref_sps_names.txt") < 0)
        {
          ERROR_MSG("Problem encountered renaming CSPS files");
          writeStatusFile(statusFileName, "Error");
          exit(-STAGE_MERGE);
        }
      }

    } // if(runMergedFlag && initialStage <= STAGE_HOMOLOGYASSEMBLY)

    if (runMergedFlag && initialStage <= STAGE_MERGE)
    {
      if (!performMergeOfCSPSAndGenoMS(ip))
      {
        ERROR_MSG("Problem encounctered during merge of CSPS and GenoMS");
        writeStatusFile(statusFileName, "Error");
        exit(-STAGE_MERGE);
      }
    }
  }
  else // if(ip.exists("FASTA_DATABASE"))
  {
    //
    //  --- need to fix performReports to allow reporting de novo results only
    //
    //  WARN_MSG("No FASTA_DATABASE specified, reporting only de novo results.\n");
    WARN_MSG("No FASTA_DATABASE specified, returning without database matches.\n");
    return 0;
  }

  if (finalStage == STAGE_HOMOLOGYASSEMBLY)
  {
    DEBUG_MSG("Option -f given. Exiting after step: " << finalStageString);
    writeStatusFile(statusFileName, "Finished");
    exit(0);
  }

  bool generatedStats = false;
  if (initialStage <= STAGE_STATPROTSEQS
      && (ip.exists("INPUT_SPEC_IDS") || ip.exists("INPUT_CLUST_SPEC_IDS")))
  {
    if (!performStatProtSeqs(ip))
    {
      ERROR_MSG("Problem encountered during StatProtSeqs stage");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_STATPROTSEQS);
    }
    generatedStats = true;
    if (commandLineParams.exists("SINGLE_STEP"))
    {
      DEBUG_MSG("Option -s given. Exiting after single step.");
      writeStatusFile(statusFileName, "Finished");
      exit(0);
    }
  }

  if (generatedStats && finalStage == STAGE_STATPROTSEQS)
  {
    DEBUG_MSG("Option -f given. Exiting after step: " << finalStageString);
    writeStatusFile(statusFileName, "Finished");
    exit(0);
  }

  //---------------------------------------------------------------------------
  // REPORT STAGE
  //---------------------------------------------------------------------------
  if (initialStage <= STAGE_REPORT)
  {
    ofstream spsProj(getProjPath(ip, "sps_projects.txt").c_str(),
                     ios_base::out | ios_base::binary);
    spsProj << "sps;.;" << ip.getValue("TOLERANCE_PEAK") << ";"
        << ip.getValue("TOLERANCE_PM") << "\n";
    spsProj.close();

    if (!performReport(ip))
    {
      ERROR_MSG("Problem encountered during Report stage");
      writeStatusFile(statusFileName, "Error");
      exit(-STAGE_REPORT);
    }
    if (commandLineParams.exists("SINGLE_STEP"))
    {
      DEBUG_MSG("Option -s given. Exiting after single step.");
      writeStatusFile(statusFileName, "Finished");
      exit(0);
    }

  }
  else
  {
    DEBUG_MSG("Bypassing reports stage");
  }
  //--------------------------------------------------------------------------
  //Set up for relaunching
  //--------------------------------------------------------------------------
  if (!generateRelaunchScript(ip))
  {
    ERROR_MSG("Problem encountered during relaunch script creation");
    writeStatusFile(statusFileName, "Error");
    exit(-STAGE_REPORT);
  }

  //---------------------------------------------------------------------------
  // END
  //---------------------------------------------------------------------------
  writeStatusFile(statusFileName, "Finished");
  return 0;
} // END

