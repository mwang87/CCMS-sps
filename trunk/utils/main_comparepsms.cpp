//
//  main_specdump - stand alone executable for dumping spectrum data
//
#include "abruijn.h"
#include "ClusterData.h"
#include "CommandLineParser.h"
#include "db_fasta.h"
#include "Logger.h"
#include "ReportTableClusterConsensus.h"
#include "spectrum.h"
#include "SpecSet.h"
#include "SpectrumPairSet.h"
#include "tuple.h"

#include <stdlib.h>

#define DEBUG_ 0
#define DEBUG_PERCENT_MATCH 0

#define MATCH_THRESHOLD 0.8

using namespace std;
using namespace sps;
using namespace specnets;
using namespace spsReports;


struct ContigResult {
  int    contig;
  int    contigPeaks;
  bool   tagMatch;
  bool   tagMatchGood;
  bool   passFdr;
  int    contigMatch;
  int    maxContigGap;
  float  specProb;
  float  alignScore;
  float  fdr;
  bool   orient;

  ContigResult() {
    contig = -1;
    contigPeaks = -1;
    tagMatch = false;
    tagMatchGood = false;
    passFdr = false;
    contigMatch = -1;
    maxContigGap = -1;
    specProb = -1;
    alignScore = -1;
    fdr = -1;
    orient = 0;
  }
};


struct SpectrumResult {
  int contig;
  int spectrum;
  int specPeaks;
  bool contigPsmGood;
  bool passFdr;
  int specMatch;
  string specAnno;
  string specAnnoClean;
  int    specAnnoLength;
  int    maxSpecGap;
  int    maxDbGap;
  string msgfAnno;
  string msgfAnnoClean;
  string specProtein;
  string msgfProtein;
  float  specProb;
  float  alignScore;
  float  fdr;
  float  msgfProb;
  float  spChange;
  float  matchPercent;
  bool   orient;

  SpectrumResult() {
    contig = -1;
    spectrum = -1;
    specPeaks = -1;
    contigPsmGood = false;
    passFdr = false;
    specMatch = -1;
    specAnno = "---";
    specAnnoClean = "---";
    specAnnoLength = -1;
    maxSpecGap = -1;
    maxDbGap = -1;
    msgfAnno = "---";
    msgfAnnoClean = "---";
    specProtein = "---";
    msgfProtein = "---";
    specProb = -1;
    alignScore = -1;
    fdr = -1;
    msgfProb = -1;
    spChange = 0;
    matchPercent = 0;
    orient = 0;
  }
};

// -------------------------------------------------------------------------
string stripAnnotation(string & annotation)
{
  string cleanAnnotation; 
  static string aminoAcids("ABCDEFGHIJKLMNOPQRSTUVWXYZ");
  for (int iChar = 0; iChar < annotation.length(); iChar++)
  {
    if (aminoAcids.find_first_of(annotation[iChar]) != string::npos)
    {
      cleanAnnotation += annotation[iChar];
    }
  }
  return cleanAnnotation;
}

// -------------------------------------------------------------------------
float getPercentMatch(string & string1, string & string2)
{
  if (DEBUG_PERCENT_MATCH) DEBUG_VAR(string1.size());
  if (DEBUG_PERCENT_MATCH) DEBUG_VAR(string2.size());
  int dynaArray[string1.size()+1][string2.size()+1];
  for (int i = 0; i < string1.size() + 1; i++) {
    for (int j = 0; j < string2.size() + 1; j++) {
      dynaArray[i][j] = 0;
    }
  }
  for (int i = 1; i <= string1.size(); i++) {
    for (int j = 1; j <= string2.size(); j++) {
      int up = dynaArray[i-1][j];
      int left = dynaArray[i][j-1];
      int diag = dynaArray[i-1][j-1];
      if (DEBUG_PERCENT_MATCH) DEBUG_MSG(i << "  " << j);
      if (DEBUG_PERCENT_MATCH) DEBUG_MSG(up << "  " << left << "  " << diag);
      if (DEBUG_PERCENT_MATCH) DEBUG_MSG(string1[i-1] << "  " << string2[j-1]);
      if (string1[i-1] == string2[j-1]) {
        diag++;
      }
      if (DEBUG_PERCENT_MATCH) DEBUG_MSG(up << "  " << left << "  " << diag);
      int max1 = max(left, up);
      dynaArray[i][j] = max(max1, diag);
      if (DEBUG_PERCENT_MATCH) DEBUG_MSG(dynaArray[i][j]);
    }
  }

  if (DEBUG_PERCENT_MATCH) {
    cout << string1 << "  " << string2 << endl;
    for (int i = 0; i <= string1.size(); i++) {
      for (int j = 0; j <= string2.size(); j++) {
        cout << dynaArray[i][j] << "  ";
      }
      cout << endl;
    }
    cout << endl;
  }
  int match = dynaArray[string1.size()][string2.size()];
  if (DEBUG_PERCENT_MATCH) DEBUG_VAR(match);
  float percent = (float)match / (float)min(string1.size(), string2.size());
  if (DEBUG_PERCENT_MATCH) DEBUG_VAR(percent);
  percent = min(percent, (float)1.0);
  if (DEBUG_PERCENT_MATCH) DEBUG_VAR(percent);
  return percent;
}


// -------------------------------------------------------------------------
int main(int argc, char ** argv)
{
  Logger::setDefaultLogger(Logger::getLogger(0));

  if (argc != 6 && argc != 2) {
    cerr << "Usage: main_scorecontig msgf_psm_file [tag_psm_file contig_psm_fdr_file spec_psm_fdr_file abinfo_file]" << endl;
    return -1;
  }
  string tagfile("assembly/tagsearchpsm.txt");
  string contigsfdrcleanfile("contigs_psm_fdr_clean.txt");
  string specfdrfile("spectrum_psm_fdr_clean.txt");
  string compfile("assembly/component_info.bin");
  string specfulltgtfile("homology/spectrum_psm_tgt.txt");
  string specfulldecfile("homology/spectrum_psm_dec.txt");
  string contigspecfile("assembly/sps_seqs.pklbin");
  string starspecfile("spectra/stars.pklbin");
  string contigpeaksfile("homology/contigs_midx_tgt.pklbin");
  string starpeaksfile("homology/spectrum_midx_tgt.pklbin");
  string contigsallpsmfile("homology/contigs_psm_tgt.txt");

  if (argc != 2) {
    tagfile= argv[2];
    contigsfdrcleanfile = argv[3];
    specfdrfile = argv[4];
    compfile = argv[5];
  }
  DEBUG_VAR(tagfile);
  DEBUG_VAR(contigsfdrcleanfile);
  DEBUG_VAR(specfdrfile);
  DEBUG_VAR(compfile);

  DEBUG_MSG("Reading data...");

  PeptideSpectrumMatchSet	psmSetMsgf;
  if (!psmSetMsgf.loadFromFile(argv[1])) {
    ERROR_MSG("Loading MSGF PSM file [" << argv[1] << "]");
    return -2;
  }
  //psmSetMsgf.saveToFile("psm_msgf.psm");

  PeptideSpectrumMatchSet	psmSetContigTag;
  if (!psmSetContigTag.loadFromFile(tagfile.c_str())) {
    ERROR_MSG("Loading Tag PSM file [" << tagfile << "]");
    return -3;
  }
  //psmSetContigTag.saveToFile("psm_contig_tag.psm");

  map<int, set<string> > mapContigToProteinSet;
  for (int iPSM = 0; iPSM < psmSetContigTag.size(); iPSM++ ) {
    psmPtr psmTag = psmSetContigTag[iPSM];
    mapContigToProteinSet[psmTag->m_scanNum].insert(psmTag->m_protein);
  } // for (int iPSM = 0; iPSM < psmSetContigTag.size(); iPSM++ )

  PeptideSpectrumMatchSet	psmSetContigFdr;
  if (!psmSetContigFdr.loadFromFile(contigsfdrcleanfile.c_str())) {
    ERROR_MSG("Loading Contig PSM file [" << contigsfdrcleanfile << "]");
    return -4;
  }
  //psmSetContigFdr.saveToFile("psm_contig.psm");

  PeptideSpectrumMatchSet	psmSetSpecFdr;
  if (!psmSetSpecFdr.loadFromFile(specfdrfile.c_str())) {
    ERROR_MSG("Loading Spectra PSM file [" << specfdrfile << "]");
    return -5;
  }
  //psmSetSpecFdr.saveToFile("psm_spec.psm");

  PeptideSpectrumMatchSet	psmSetSpecFullTgt;
  if (!psmSetSpecFullTgt.loadFromFile(specfulltgtfile.c_str())) {
    ERROR_MSG("Loading Spectra Full PSM file [" << specfulltgtfile << "]");
    return -6;
  }
  //psmSetSpecFullTgt.saveToFile("psm_full_spec.psm");

  PeptideSpectrumMatchSet	psmSetSpecFullDec;
  if (!psmSetSpecFullDec.loadFromFile(specfulldecfile .c_str())) {
    ERROR_MSG("Loading Spectra Full PSM file [" << specfulldecfile << "]");
    return -6;
  }
  //psmSetSpecFullDec.saveToFile("psm_full_spec.psm");

  abinfo_t contigAbinfo;
  if (!Load_abinfo(compfile.c_str(), contigAbinfo)) {
    ERROR_MSG("Loading Abinfo file [" << compfile.c_str() << "]");
    return -7;
  }
  //DEBUG_VAR(contigAbinfo.size());

  SpecSet contigSpectra;
  if (contigSpectra.loadPklBin(contigspecfile.c_str()) <= 0) {
    ERROR_MSG("Loading Contig Spectra file [" << contigspecfile << "]");
    return -8;
  }
  //DEBUG_VAR(contigSpectra.size());

  SpecSet starSpectra;
  if (starSpectra.loadPklBin(starspecfile.c_str()) <= 0) {
    ERROR_MSG("Loading Star Spectra file [" << starspecfile << "]");
    return -9;
  }
  //DEBUG_VAR(starSpectra.size());

  SpecSet contigPeaks;
  if (contigPeaks.loadPklBin(contigpeaksfile.c_str()) <= 0) {
    ERROR_MSG("Loading contig peaks file [" << contigpeaksfile << "]");
    return -9;
  }
  //DEBUG_VAR(contigPeaks.size());

  SpecSet starPeaks;
  if (starPeaks.loadPklBin(starpeaksfile.c_str()) <= 0) {
    ERROR_MSG("Loading star peaks file [" << starspecfile << "]");
    return -9;
  }
  //DEBUG_VAR(starPeaks.size());

  PeptideSpectrumMatchSet	psmSetContigsAll;
  if (!psmSetContigsAll.loadFromFile(contigsallpsmfile.c_str())) {
    ERROR_MSG("Loading Contigs PSM file [" << contigsallpsmfile << "]");
    return -2;
  }
  //psmSetContigsAll.saveToFile("psm_msgf.psm");

  DEBUG_MSG("Iterating through all contigs/spectra...");

  //----------------------------------------------------------
  // Iterate through the abinfo and get all contigs and spectra
  //----------------------------------------------------------
  int totalSpecsInContigs = 0;
  map<int, int> specToContig;  
  map<int, set<int> > contigToSpecSet;  
  std::map<unsigned, // contig index
      std::pair<std::pair<vector<int> , vector<int> >, // spectrum index, flipped(1)/not-flipped(0)
          vector<std::pair< // ABruijn vertices
              vector<int> , vector<double> > // Spectrum index, peak mass
          > > >::iterator itr = contigAbinfo.begin();
  std::map<unsigned, // contig index
      std::pair<std::pair<vector<int> , vector<int> >, // spectrum index, flipped(1)/not-flipped(0)
          vector<std::pair< // ABruijn vertices
              vector<int> , vector<double> > // Spectrum index, peak mass
          > > >::iterator itr_end = contigAbinfo.end();
  for ( ; itr != itr_end; itr++) {
    int contigIndex = itr->first + 1;  // We prefer 1-based contigs
    //DEBUG_MSG("C = " << contigIndex);
    vector<int> specs = itr->second.first.first;
//    totalSpecsInContigs += specs.size();
    for (int i = 0; i < specs.size(); i++) {

      if (specs[i] >= starSpectra.size()) {
        WARN_MSG(specs[i] << " out of bounds!");
        continue;
      }
      int specIndex = specs[i] + 1;  // We prefer 1-based spectra
      //DEBUG_MSG("  S = " << specIndex);
      if (specToContig.find(specIndex) != specToContig.end()) {
        WARN_MSG(specIndex << " was found twice!");
        WARN_MSG("  " << contigIndex << " and " << specToContig[specIndex ]);
      } else {
        specToContig[specIndex] = contigIndex;
        contigToSpecSet[contigIndex].insert(specIndex);
        totalSpecsInContigs++;
      }
    }
  }
  DEBUG_VAR(specToContig.size());
  DEBUG_VAR(contigToSpecSet.size());


  DEBUG_MSG("Setting up results maps...");

  //----------------------------------------------------------
  // Setup the contig results vector
  //----------------------------------------------------------
  map<int, ContigResult> contigResults;
  map<int, set<int> >::iterator itrc = contigToSpecSet.begin();  
  map<int, set<int> >::iterator itrcEnd = contigToSpecSet.end();  
  for ( ; itrc != itrcEnd; itrc++) {
    //DEBUG_VAR(itrc->first);
    ContigResult newContigResult;
    contigResults[itrc->first] = newContigResult;
    contigResults[itrc->first].contig = itrc->first;
    contigResults[itrc->first].tagMatch = (mapContigToProteinSet.find(itrc->first) != mapContigToProteinSet.end());

    int contigIndex0 = itrc->first - 1;
    contigResults[itrc->first].contigPeaks = contigSpectra[contigIndex0].size();
    contigResults[itrc->first].maxContigGap = -2;
    if (contigPeaks[contigIndex0].size() > 0) {
      for (int i = 0; i < contigPeaks[contigIndex0].size() - 1; i++) {
        float peakGapContig = contigPeaks[contigIndex0][i+1][0] - contigPeaks[contigIndex0][i][0];
        if (peakGapContig > contigResults[itrc->first].maxContigGap) {
          contigResults[itrc->first].maxContigGap = peakGapContig;
        }
      }
    }
  }

  //----------------------------------------------------------
  // Setup the spectrum results vector
  //----------------------------------------------------------
  map<int, SpectrumResult> specResults;
  map<int, int>::iterator itrx = specToContig.begin();  
  map<int, int>::iterator itrxEnd = specToContig.end();  
  for ( ; itrx != itrxEnd; itrx++) {
    //DEBUG_VAR(itrx->first);
    SpectrumResult newSpecResult;
    specResults[itrx->first] = newSpecResult;
    specResults[itrx->first].contig = itrx->second;
    specResults[itrx->first].spectrum = itrx->first;

    int specIndex0 = itrx->first - 1;
    if (specIndex0 < 0 || specIndex0 >= starSpectra.size()) {
      continue;
    }

    specResults[itrx->first].specPeaks = starSpectra[specIndex0].size();

    specResults[itrx->first].maxSpecGap = -2;
    if (starPeaks[specIndex0].size() > 0) {
      for (int i = 0; i < starPeaks[specIndex0].size() - 1; i++) {
        if (itrx->first == 13057) DEBUG_MSG( i << "  " << starPeaks[specIndex0][i+1][0] << "  " <<  starPeaks[specIndex0][i][0]);
        float peakGapSpec = starPeaks[specIndex0][i+1][0] - starPeaks[specIndex0][i][0];
        if (peakGapSpec > specResults[itrx->first].maxSpecGap) {
          specResults[itrx->first].maxSpecGap = peakGapSpec;
        }
      }
    }
    specResults[itrx->first].maxDbGap  = -2;
    if (starPeaks[specIndex0].size() > 0) {
      for (int i = 0; i < starPeaks[specIndex0].size() - 1; i++) {
        if (itrx->first == 13057) DEBUG_MSG( i << "  " << starPeaks[specIndex0][i+1][1] << "  " <<  starPeaks[specIndex0][i][1]);
        float peakGapDb = starPeaks[specIndex0][i+1][1] - starPeaks[specIndex0][i][1];
        if (peakGapDb > specResults[itrx->first].maxDbGap) {
          specResults[itrx->first].maxDbGap = peakGapDb;
        }
      }
    }

  }

  DEBUG_MSG("Populating simple fields...");

  //----------------------------------------------------------
  // Create the map of contigs to their annotations
  //----------------------------------------------------------
  map<int, set<string> > mapContigToSetAnno;
  for (int iPSM = 0; iPSM < psmSetContigsAll.size(); iPSM++ ) {
    psmPtr psmSpec = psmSetContigsAll[iPSM];
    mapContigToSetAnno[psmSpec->m_scanNum].insert(psmSpec->m_annotation);
  }

  //----------------------------------------------------------
  // Populate some contig fields
  //----------------------------------------------------------
  for (int iPSM = 0; iPSM < psmSetContigFdr.size(); iPSM++ ) {
    psmPtr psmSpec = psmSetContigFdr[iPSM];
    contigResults[psmSpec->m_scanNum].passFdr = true;
    contigResults[psmSpec->m_scanNum].specProb = psmSpec->m_pValue;
    contigResults[psmSpec->m_scanNum].alignScore = psmSpec->m_score;
    contigResults[psmSpec->m_scanNum].fdr = psmSpec->m_fdr;
    contigResults[psmSpec->m_scanNum].orient = psmSpec->m_matchOrientation;
  }

  //----------------------------------------------------------
  // Set spectrum pass FDR flag
  //----------------------------------------------------------
  for (int iPSM = 0; iPSM < psmSetSpecFdr.size(); iPSM++ ) {
    psmPtr psmFdr = psmSetSpecFdr[iPSM];
    specResults[psmFdr->m_scanNum].passFdr = true;
  }

  //----------------------------------------------------------
  // Create the full PSM map
  //----------------------------------------------------------
  map<int, psmPtr> mapScanFullPsm;
  for (int iPSM = 0; iPSM < psmSetSpecFullTgt.size(); iPSM++ ) {
    psmPtr psmSpec = psmSetSpecFullTgt[iPSM];
    if (mapScanFullPsm.find(psmSpec->m_scanNum) == mapScanFullPsm.end()) {
      mapScanFullPsm[psmSpec->m_scanNum] = psmSpec;
    } else {
      if (psmSpec->m_score > mapScanFullPsm[psmSpec->m_scanNum]->m_score) {
        mapScanFullPsm[psmSpec->m_scanNum] = psmSpec;
      }
    }
  }
  for (int iPSM = 0; iPSM < psmSetSpecFullDec.size(); iPSM++ ) {
    psmPtr psmSpec = psmSetSpecFullDec[iPSM];
    if (mapScanFullPsm.find(psmSpec->m_scanNum) == mapScanFullPsm.end()) {
      mapScanFullPsm[psmSpec->m_scanNum] = psmSpec;
    } else {
      if (psmSpec->m_score > mapScanFullPsm[psmSpec->m_scanNum]->m_score) {
        mapScanFullPsm[psmSpec->m_scanNum] = psmSpec;
      }
    }
  }

  DEBUG_TRACE;
  //----------------------------------------------------------
  // Populate some spectrum fields
  //----------------------------------------------------------
  map<int, psmPtr>::iterator itrp = mapScanFullPsm.begin();
  map<int, psmPtr>::iterator itrpEnd = mapScanFullPsm.end();
  for ( ; itrp != itrpEnd; itrp++) {
    psmPtr psmSpec = itrp->second;
    specResults[psmSpec->m_scanNum].specAnno = psmSpec->m_annotation;
    specResults[psmSpec->m_scanNum].specAnnoClean = stripAnnotation(specResults[psmSpec->m_scanNum].specAnno);
    specResults[psmSpec->m_scanNum].specAnnoLength = specResults[psmSpec->m_scanNum].specAnnoClean.length();
    specResults[psmSpec->m_scanNum].specProtein = psmSpec->m_protein;
    specResults[psmSpec->m_scanNum].specAnno = psmSpec->m_annotation;
    specResults[psmSpec->m_scanNum].alignScore = psmSpec->m_score;
    specResults[psmSpec->m_scanNum].specProb = psmSpec->m_pValue;
    specResults[psmSpec->m_scanNum].fdr = psmSpec->m_fdr;
    specResults[psmSpec->m_scanNum].orient = psmSpec->m_matchOrientation;
  }

  DEBUG_TRACE;
  //----------------------------------------------------------
  // Create MSGF ID map
  //----------------------------------------------------------
  int stars_in_contigs_msgf_id = 0;
  map<int, psmPtr> mapScanMsgf;
  for (int iPSM = 0; iPSM < psmSetMsgf.size(); iPSM++ ) {

    psmPtr psmMsgf = psmSetMsgf[iPSM];
    mapScanMsgf[psmMsgf->m_scanNum] = psmMsgf;

    if (specToContig.find(psmMsgf->m_scanNum) != specToContig.end()) {
      stars_in_contigs_msgf_id++;
    }

    specResults[psmMsgf->m_scanNum].msgfAnno = psmMsgf->m_annotation;
    specResults[psmMsgf->m_scanNum].msgfAnnoClean = stripAnnotation(psmMsgf->m_annotation);
    specResults[psmMsgf->m_scanNum].msgfProtein = psmMsgf->m_protein;
    specResults[psmMsgf->m_scanNum].msgfProb = -log10(psmMsgf->m_score);

    if (specResults[psmMsgf->m_scanNum].alignScore == 0.0 || specResults[psmMsgf->m_scanNum].msgfProb == 0.0) {
      specResults[psmMsgf->m_scanNum].spChange = 0.0;
    } else {
      specResults[psmMsgf->m_scanNum].spChange = specResults[psmMsgf->m_scanNum].specProb / 
                                                 specResults[psmMsgf->m_scanNum].msgfProb;
    }

  } // for (int iPSM = 0; iPSM < psmSetMsgf.size(); iPSM++ )


  int nUnkonwnContigs = contigAbinfo.size();

  int totalContigsWithTags = 0;


  DEBUG_MSG("Finding contig/spectra tags...");
  //----------------------------------------------------------
  // Figure out if each contig has any tag and a good tag or not
  //----------------------------------------------------------
  for (int iTagPSM = 0; iTagPSM < psmSetContigTag.size(); iTagPSM++ ) {
    psmPtr psmTagContig = psmSetContigTag[iTagPSM];
    int contigIndex = psmTagContig->m_scanNum;

    if (contigResults[contigIndex].tagMatchGood == 1) {
      continue;
    }

    //DEBUG_MSG(psmTagContig->m_scanNum << "  " << psmTagContig->m_protein);
    set<int>::iterator itrs = contigToSpecSet[contigIndex].begin();  
    set<int>::iterator itrsEnd = contigToSpecSet[contigIndex].end();  
    for ( ; itrs != itrsEnd; itrs++) {

      int specIndex = *itrs;
      //DEBUG_MSG("  " << specIndex);
      if (mapScanMsgf.find(specIndex) == mapScanMsgf.end()) {
        continue;
      }
      string cleanContig = stripAnnotation(psmTagContig->m_annotation);
      //DEBUG_MSG("  " << cleanContig);

      if (specResults[specIndex].msgfAnnoClean.find(cleanContig) != string::npos) {
        contigResults[contigIndex].tagMatchGood = true;
        totalContigsWithTags++;
        //DEBUG_MSG("    YES");
        //break;  // only need one good tag
      }

      set<string>::iterator itrs2 = mapContigToSetAnno[contigIndex].begin();  
      set<string>::iterator itrs2End = mapContigToSetAnno[contigIndex].end();
      for ( ; itrs2 != itrs2End ; itrs2++) {
        string contigAnno = *itrs2;
        //DEBUG_MSG("  " << contigAnno);
        string cleanAnno = stripAnnotation(contigAnno);
        if (specResults[specIndex].msgfAnnoClean.find(cleanAnno) != string::npos) {
          specResults[specIndex].contigPsmGood = true;
          //DEBUG_MSG("    YES2");
          break; // Only need one
        }
      }

    } // for ( ; itrs != itrs_end; itrs++) {

  } // for (int iTagPSM = 0; iTagPSM < psmSetContigTag.size(); iTagPSM++ ) {

  DEBUG_VAR(totalContigsWithTags);

  int nNoMatch = 0;
  int nDiffProtein = 0;
  int nExactMatch = 0;
  int nMatchG80 = 0;
  int nDiffPeptide = 0;

  set<int> specPSMSet;  
  set<int> setNoMsgf;  
  set<int> setExactMatch;  
  set<int> setInexactMatch;  
  set<int> setAssociatedPeptide;
  set<int> setAssociatedBadPeptide;
  set<int> setDiffPeptide;
  map<int, int> mapSpecToMatchType;
  map<int, int> mapSpectrumToGoodMsgfID;

  DEBUG_MSG("Matching MSGF results...");
  //---------------------------------------------------------------------
  // Determine if the spectra from specnets match MSGFDB results
  //---------------------------------------------------------------------
  map<int, psmPtr>::iterator itrFull = mapScanFullPsm.begin();
  map<int, psmPtr>::iterator itrFullEnd = mapScanFullPsm.end();
  for ( ; itrFull != itrFullEnd; itrFull++ ) {

    const psmPtr & psmSpec = itrFull->second;
    int specIndex = psmSpec->m_scanNum;
    //DEBUG_VAR(specIndex);

    specPSMSet.insert(specIndex);
    int contigIndex = specToContig[specIndex];
    //DEBUG_VAR(contigIndex);

    string cleanAnnotationSpec = stripAnnotation(psmSpec ->m_annotation);
    //DEBUG_VAR(cleanAnnotationSpec);

    map<int, psmPtr>::iterator itrFind = mapScanMsgf.find(specIndex);
    if (itrFind == mapScanMsgf.end()) {
      if (specResults[specIndex].passFdr) setNoMsgf.insert(specIndex);
      continue;
    }

    psmPtr psmMsgf = itrFind->second;
    string cleanAnnotationMsgf = stripAnnotation(psmMsgf->m_annotation);
    //DEBUG_VAR(cleanAnnotationMsgf);

    mapSpectrumToGoodMsgfID[specIndex] = 0;

    map<int, set<string>  >::iterator itrTag = mapContigToProteinSet.find(contigIndex);
    if (itrTag != mapContigToProteinSet.end()) {
      set<string>::iterator itrProtein = itrTag->second.begin();
      set<string>::iterator itrProteinEnd = itrTag->second.end();
      for ( ; itrProtein != itrProteinEnd; itrProtein++) {
        if (itrProtein->find(psmMsgf->m_protein) != string::npos) {
          //DEBUG_MSG("  " << "found protein [" << contigIndex << "] protein [" << *itrProtein << "]");
          mapSpectrumToGoodMsgfID[specIndex] = 1;
          break;
        }
      }
    }

    float percentMatch = getPercentMatch(cleanAnnotationSpec, cleanAnnotationMsgf);
    specResults[specIndex].matchPercent = percentMatch;

    if (percentMatch >= MATCH_THRESHOLD) {
      //DEBUG_MSG("MATCH");
      nExactMatch++;
      if (specResults[specIndex].passFdr) {
        setExactMatch.insert(specIndex);
      }
      mapSpecToMatchType[specIndex] = 1;
      mapSpectrumToGoodMsgfID[specIndex] = 1;
      continue;
    }

    if (specResults[specIndex].passFdr) {
      //DEBUG_MSG("NON MATCH");
      setDiffPeptide.insert(specIndex);
    }
    mapSpecToMatchType[specIndex] = 3;
    //DEBUG_MSG(specIndex << "  " << mapSpecToMatchType[specIndex]);

  } // for (int iPSM = 0; iPSM < psmSetSpecFull.size(); iPSM++ )

  DEBUG_MSG("Finding similar spectra...");
  //---------------------------------------------------------------------
  // Find similar spectra in contig to figure out if unknowns are good or bad
  //    and if MSGFDB got the wrong answer
  //---------------------------------------------------------------------
  int nSimGood = 0;
  int nSimBad = 0;

  // Find spectra that can be matched to others in the contig
  map<int, psmPtr>::iterator itrm = mapScanFullPsm.begin();
  map<int, psmPtr>::iterator itrmEnd = mapScanFullPsm.end();
  for ( ; itrm != itrmEnd; itrm++) {
    const psmPtr & psmSpec = itrm->second;

//  for (int iPSM = 0; iPSM < psmSetSpecFull.size(); iPSM++ ) {
//    const psmPtr & psmSpec = psmSetSpecFull[iPSM];

    int specIndex = psmSpec->m_scanNum;
    int contigIndex = specToContig[specIndex];
    int matchType = mapSpecToMatchType[specIndex];
    string cleanAnnotationSpec1 = stripAnnotation(psmSpec ->m_annotation);
    // Only check previously unknown spectra in known contigs
    if (matchType != 0 && matchType != 3) {
      continue;
    }

    map<int, psmPtr>::iterator itrm2 = mapScanFullPsm.begin();
    map<int, psmPtr>::iterator itrm2End = mapScanFullPsm.end();
    for ( ; itrm2 != itrm2End; itrm2++) {
      const psmPtr & psmSpec2 = itrm2->second;

//    for (int iPSM2 = 0; iPSM2 < psmSetSpecFull.size(); iPSM2++ ) {
//      const psmPtr & psmSpec2 = psmSetSpecFull[iPSM2];

      int specIndex2 = psmSpec2->m_scanNum;
      if (specIndex == specIndex2) {
        continue;
      }
      int contigIndex2 = specToContig[specIndex2];
      if (contigIndex != contigIndex2) {
        continue;
      }
      int matchType2 = mapSpecToMatchType[specIndex2];
      if (matchType2 != 1 && matchType2 != 3) {
        continue;
      }

      string cleanAnnotationSpec2 = stripAnnotation(psmSpec2->m_annotation);
      float percentMatch = getPercentMatch(cleanAnnotationSpec1, cleanAnnotationSpec2);
      //DEBUG_MSG(contigIndex << "  " << contigIndex2 << "  " << specIndex << "  " << specIndex2);
      //DEBUG_MSG("  " << cleanAnnotationSpec1 << "  " << cleanAnnotationSpec2 << "  " << percentMatch << "  " << matchType2);
      specResults[specIndex].matchPercent = -percentMatch;
      if (percentMatch > MATCH_THRESHOLD) {
        if (matchType != 3) {
          if (matchType2 == 1 || matchType2 == 2) {
            mapSpecToMatchType[specIndex] = 4;
            mapSpectrumToGoodMsgfID[specIndex] = 1;
            if (specResults[specIndex].passFdr) nSimGood++;
          } else {
            mapSpecToMatchType[specIndex] = 5;
            if (specResults[specIndex].passFdr) nSimBad++;
          }
        } else if (matchType == 3 && matchType2 == 1) {
        }
        break;
      }
    }
  }

  DEBUG_MSG("Finding contig hybrids...");
  //---------------------------------------------------------------------
  // Figure out if the contigs are good, bad or hybrid
  //---------------------------------------------------------------------
  map<int, string> mapContigToMsgfProtein;
  map<int, int> mapContigMsgfHybrid;

  map<int, int>::iterator itrs = specToContig.begin();
  map<int, int>::iterator itrsEnd = specToContig.end();
  for ( ; itrs != itrsEnd; itrs++) {

    int specIndex = itrs->first;
    int contigIndex = itrs->second;

    map<int, psmPtr>::iterator itrFindMsgf =  mapScanMsgf.find(specIndex);
    if (itrFindMsgf == mapScanMsgf.end()) {
      map<int, string>::iterator itrFindProt =  mapContigToMsgfProtein.find(contigIndex);
      if (itrFindProt == mapContigToMsgfProtein.end()) {
      mapContigMsgfHybrid[contigIndex] = -1;
      }
      continue;
    }

    const psmPtr & psmSpec = itrFindMsgf->second;

    map<int, string>::iterator itrFindProt =  mapContigToMsgfProtein.find(contigIndex);
    if (itrFindProt == mapContigToMsgfProtein.end()) {
      mapContigToMsgfProtein[contigIndex] = psmSpec->m_protein;
      mapContigMsgfHybrid[contigIndex] = 0;
    } else {
      string prevProtein = itrFindProt->second;
      if (itrFindMsgf->second->m_protein != prevProtein) {
        mapContigMsgfHybrid[contigIndex] = 1;
      }
    }

  }

  //---------------------------------------------------------------------
  // Figure out if the contigs are good, bad or hybrid
  //---------------------------------------------------------------------
  map<int, vector<int> > mapContigVecType;
  map<int, int > mapContigHasTag;

  itrm = mapScanFullPsm.begin();
  for ( ; itrm != itrmEnd; itrm++) {
    const psmPtr & psmSpec = itrm->second;

//  for (int iPSM = 0; iPSM < psmSetSpecFull.size(); iPSM++ ) {
//    const psmPtr & psmSpec = psmSetSpecFull[iPSM];

    int specIndex = psmSpec->m_scanNum;
    int contigIndex = specToContig[specIndex];
    int matchType = mapSpecToMatchType[specIndex];
    int contigType = 0;
    if (matchType == 1 || matchType == 2 || matchType == 4 || matchType == 6) contigType = 1;
    if (matchType == 3 || matchType == 5) contigType = 3;
    mapContigVecType[contigIndex].push_back(contigType);
  }

  int nGoodContig = 0;
  int nBadContig = 0;
  int nHybridContig = 0;
  map<int, int> mapContigType;
  map<int, vector<int> >::iterator itrv = mapContigVecType.begin();
  map<int, vector<int> >::iterator itrvEnd = mapContigVecType.end();
  for (; itrv != itrvEnd; itrv++) {
    int contigIndex = itrv->first;
    vector<int> & vecType = itrv->second;
    for (int i = 0; i < vecType.size(); i++) {
      if (vecType[i] != 0) {
        if (mapContigType.find(contigIndex) == mapContigType.end()) {
          mapContigType[contigIndex] = vecType[i];
        } else if (mapContigType[contigIndex] != vecType[i]) {
          mapContigType[contigIndex] = 8;
        }
      }
    }
    if (mapContigType[contigIndex] == 1) nGoodContig++;
    if (mapContigType[contigIndex] == 3) nBadContig++;
    if (mapContigType[contigIndex] == 8) nHybridContig++;
  }

  DEBUG_MSG("Outputting results tables...");
  //---------------------------------------------------------------------
  // Output the spectrum results table
  //---------------------------------------------------------------------
  cout << "contig#\tscan#\tcontig_peaks\tspec_peaks\ttag_match\ttag_match_good\tgood_contig_psm\tpass_fdr\tcontig_match\tmsgf_hybrid\tspec_match\tspec_anno\tlength\tmax_contig_gap\tmax_spec_gap\tmax_db_gap\tmsgf_anno\tspec_prot\tmsgf_prot\tspec_prob\talign_score\tfdr\tmsgf_prob\tSP_change\tmatch%\torient" << endl;

  int totalSpectraWithGoodTags = 0;

  DEBUG_VAR(mapScanFullPsm.size());
  DEBUG_VAR(specToContig.size());

  int stars_id_msgf_noid_nocontigpsm = 0;
  int stars_id_msgf_noid_belowfdr = 0;
  int stars_id_msgf_noid_nogoodtag = 0;
  int stars_id_msgf_noid_hybridcontig = 0;
  int stars_id_msgf_noid_badcontig = 0;

  itrs = specToContig.begin();
  for ( ; itrs != itrsEnd; itrs++) {

    int specIndex = itrs->first;
    int contigIndex = specResults[specIndex].contig;

    cout << contigIndex << "\t";
    cout << specIndex << "\t";

    cout << contigResults[contigIndex].contigPeaks << "\t";
    cout << specResults[specIndex].specPeaks << "\t";

    cout << contigResults[contigIndex].tagMatch << "\t";
    cout << contigResults[contigIndex].tagMatchGood << "\t";
    cout << specResults[specIndex].contigPsmGood << "\t";
    cout << specResults[specIndex].passFdr << "\t";

    switch (mapContigType[contigIndex]) {
    case 1:
      cout << "good\t";
      break;
    case 3:
      cout << "bad\t";
      break;
    case 8:
      cout << "hybrid\t";
      break;
    default:
      cout << "unk\t";
      break;
    }

    switch (mapContigMsgfHybrid[contigIndex]) {
    case 0:
      cout << "homog\t";
      break;
    case 1:
      cout << "hybrid\t";
      break;
    case -1:
    default:
      cout << "unk\t";
      break;
    }

    //DEBUG_MSG(specIndex << "  " << mapSpecToMatchType[specIndex]);
    switch (mapSpecToMatchType[specIndex]) {
    case 1:
    case 2:
      cout << "match\t";
      break;
    case 3:
      cout << "nonmatch\t";
      break;
    case 4:
      cout << "simgood\t";
      break;
    case 5:
      cout << "simbad\t";
      break;
    case 6:
      cout << "msgfbad\t";
      break;
    default:
      cout << "unk\t";
      break;
    }

    cout << specResults[specIndex].specAnno << "\t"; 
    cout << specResults[specIndex].specAnnoLength << "\t"; 

    cout << contigResults[contigIndex].maxContigGap << "\t"; 
    cout << specResults[specIndex].maxSpecGap << "\t"; 
    cout << specResults[specIndex].maxDbGap << "\t"; 

    cout << specResults[specIndex].msgfAnno << "\t"; 
    cout << specResults[specIndex].specProtein << "\t";
    cout << specResults[specIndex].msgfProtein << "\t";
    cout << specResults[specIndex].specProb << "\t";
    cout << specResults[specIndex].alignScore << "\t";
    cout << specResults[specIndex].fdr << "\t";
    cout << specResults[specIndex].msgfProb << "\t";
    cout << specResults[specIndex].spChange << "\t";
    cout << specResults[specIndex].matchPercent << "\t";
    cout << specResults[specIndex].orient << endl;

    if (contigResults[contigIndex].tagMatchGood) {
      totalSpectraWithGoodTags++;
    }

    if (!specResults[specIndex].contigPsmGood && 
        !specResults[specIndex].passFdr &&
        (specResults[specIndex].msgfProb > 0)) {
      stars_id_msgf_noid_nocontigpsm++;
    }

    if (!specResults[specIndex].contigPsmGood && 
        !specResults[specIndex].passFdr &&
        (specResults[specIndex].msgfProb > 0) &&
        !contigResults[contigIndex].tagMatchGood) {
      stars_id_msgf_noid_nogoodtag++;
    }

    if (!specResults[specIndex].contigPsmGood && 
        !specResults[specIndex].passFdr &&
        (specResults[specIndex].msgfProb > 0) &&
        contigResults[contigIndex].tagMatchGood &&
        (mapContigType[contigIndex] == 8)) {
      stars_id_msgf_noid_hybridcontig++;
    }

    if (!specResults[specIndex].contigPsmGood && 
        !specResults[specIndex].passFdr &&
        (specResults[specIndex].msgfProb > 0) &&
        contigResults[contigIndex].tagMatchGood &&
        (mapContigType[contigIndex] == 3)) {
      stars_id_msgf_noid_badcontig++;
    }

    if (specResults[specIndex].contigPsmGood && 
        !specResults[specIndex].passFdr &&
        (specResults[specIndex].msgfProb > 0)) {
      stars_id_msgf_noid_belowfdr++;
    }

  }

  cout << endl;
  cout << "stars_in_contigs_all\t" << totalSpecsInContigs << endl;
  cout << "stars_in_contigs_good_tag\t" << totalSpectraWithGoodTags << endl;
  cout << "stars_in_contigs_fdr\t" << psmSetSpecFdr.size() << endl;
  cout << "stars_in_contigs_msgf_id\t" << stars_in_contigs_msgf_id << endl;
  cout << "stars_id_msgf_match\t" << setExactMatch.size() << endl;
  cout << "stars_id_msgf_mismatch\t" << setDiffPeptide.size() << endl;

  cout << "stars_id_msgf_noid\t" << stars_in_contigs_msgf_id - setExactMatch.size() - setDiffPeptide.size() << endl;
  cout << "stars_id_msgf_noid_nocontigpsm\t" << stars_id_msgf_noid_nocontigpsm << endl;
  cout << "stars_id_msgf_noid_nogoodtag\t" << stars_id_msgf_noid_nogoodtag << endl;
  cout << "stars_id_msgf_noid_hybridcontig\t" << stars_id_msgf_noid_hybridcontig << endl;
  cout << "stars_id_msgf_noid_badcontig\t" << stars_id_msgf_noid_badcontig << endl;
  cout << "stars_id_msgf_noid_belowfdr\t" << stars_id_msgf_noid_belowfdr << endl;

  cout << "stars_in_contigs_msgf_unknown\t" << setNoMsgf.size() << endl;
  cout << "stars_id_msgf_similar_good\t" << nSimGood << endl;
  cout << "stars_id_msgf_similar_bad\t" << nSimBad << endl;

  //---------------------------------------------------------------------
  // Output the contig results table
  //---------------------------------------------------------------------
  cout << endl;
  cout << "contig#\tcontig_peaks\thave_tag_match\ttag_match_good\tpass_fdr\tmsgf_contig\tcontig_type\tmax_contig_gap\talign_score\tfdr\torient" << endl;

  DEBUG_VAR(mapScanFullPsm.size());
  DEBUG_VAR(specToContig.size());

  itrc = contigToSpecSet.begin();  
  for ( ; itrc != itrcEnd; itrc++) {

    int contigIndex = itrc->first;

    cout << contigIndex << "\t";
    cout << contigResults[contigIndex].contigPeaks << "\t";
    cout << contigResults[contigIndex].tagMatch << "\t";
    cout << contigResults[contigIndex].tagMatchGood << "\t";
    cout << specResults[contigIndex].passFdr << "\t";

    switch (mapContigMsgfHybrid[contigIndex]) {
    case 0:
      cout << "homog\t";
      break;
    case 1:
      cout << "hybrid\t";
      break;
    case -1:
    default:
      cout << "unk\t";
      break;
    }


    switch (mapContigType[contigIndex]) {
    case 1:
      cout << "good\t";
      break;
    case 3:
      cout << "bad\t";
      break;
    case 8:
      cout << "hybrid\t";
      break;
    default:
      cout << "unk\t";
      break;
    }

    cout << contigResults[contigIndex].maxContigGap << "\t"; 
    cout << contigResults[contigIndex].alignScore << "\t";
    cout << contigResults[contigIndex].fdr << "\t";
    cout << contigResults[contigIndex].orient << endl;
  }


  cout << endl;
  cout << "contigs_total\t" << contigAbinfo.size() << endl;
  cout << "contigs_correct_tag\t" << totalContigsWithTags << endl;
  cout << "contigs_correct\t" << nGoodContig << endl;
  cout << "contigs_incorrect\t" << nBadContig << endl;
  cout << "contigs_hybrid\t" << nHybridContig << endl;
  cout << "contigs_unknown\t" << contigAbinfo.size() - nGoodContig - nBadContig - nHybridContig << endl;
  cout << endl;

  return 0;
}
