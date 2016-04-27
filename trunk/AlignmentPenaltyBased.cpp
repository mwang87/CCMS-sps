#include "AlignmentPenaltyBased.h"
#include "Logger.h"

#include <algorithm>
#include <iostream>
#include <string>
#include <limits.h>
#include <stdio.h>

#define DEBUG_RANGE 0
#define DEBUG_RANGE2 0
#define DEBUG_AAS 0
#define DEBUG_SPECS 0
#define DEBUG_ALIGN 0
#define DEBUG_ALIGN2 0
#define DEBUG_ALIGN3 0
#define DEBUG_ALIGN_NTERM 0
#define DEBUG_GAP 0
#define DEBUG_GAP_ANNO 0
#define DEBUG_ALIGN_GAP 0
#define DEBUG_CACHE 0


// Assuming vectors of size ~1000 then the cache size will end up in the MB region (not GB)
#define CACHE_LIMIT 10000

namespace specnets
{
  //--------------------------------------------------------------------------
  void debugAAs(PenaltyMatrix * modPenaltyMatrix,
		PenaltyMatrix * blossumPenaltyMatrix)
  {
    for (int i = 0; i < 26; i++) {
      string strAA;  // Create a dummy string
      strAA.append(1, char('A' + i)); // Append the DB AA char
      DEBUG_MSG(strAA << " = " << modPenaltyMatrix->getMass(strAA));
    }
    for (int i = 0; i < 26; i++) {
      string strAA;  // Create a dummy string
      strAA.append(1, char('A' + i)); // Append the DB AA char
      DEBUG_MSG(strAA << " = " << blossumPenaltyMatrix->getMass(strAA));
    }
    return;
  }

  //--------------------------------------------------------------------------
  void debugDBSpec(char * dbSeq, Spectrum & dbSpec)
  {
    DEBUG_MSG(0 << ", " << dbSpec[0][0] << ", " << dbSpec[0][0]);
    for (int dbSpecIdx = 1; dbSpecIdx < dbSpec.size(); dbSpecIdx++) {
      DEBUG_MSG(dbSeq[dbSpecIdx - 1] << ", " << dbSpecIdx << ", " << dbSpec[dbSpecIdx][0] << ", " << dbSpec[dbSpecIdx][0] - dbSpec[dbSpecIdx-1][0]);
    }
    return;
  }
  
  //--------------------------------------------------------------------------
  void debugSpec(Spectrum & spec)
  {
    DEBUG_MSG(0 << ", " << spec[0][0] << ", " << spec[0][1]);
    for (int specIdx = 1; specIdx < spec.size(); specIdx++) {
      DEBUG_MSG(specIdx << ", " << spec[specIdx][0] - spec[specIdx-1][0] << ", " << spec[specIdx][1]);
    }
    return;
  }

  //--------------------------------------------------------------------------
  void debugMatrices(Spectrum & 	spec,
		     Spectrum & 	dbSpec,
		     vector<vector< float > > & matchMatrix,
		     vector<vector< string > > & matchType)
  {
    for (int dbSpecIdx = 0; dbSpecIdx < dbSpec.size(); dbSpecIdx++) {
      bool outputRow = false;
      // Only output "interesting" rows
      for (int specIdx = 1; specIdx < spec.size(); specIdx++) {
        if (matchMatrix[dbSpecIdx][specIdx] != -(float)INT_MAX) {
          outputRow = true;
          break;
        }
      }
      if (outputRow) {
        cout << dbSpecIdx << "\t";
        for (int specIdx = 0; specIdx < spec.size(); specIdx++) {
          cout << matchMatrix[dbSpecIdx][specIdx] << "\t";
        }
        cout << endl;
      }
    }

    for (int dbSpecIdx = 0; dbSpecIdx < dbSpec.size(); dbSpecIdx++) {
      bool outputRow = false;
      // Only output "interesting" rows
      for (int specIdx = 1; specIdx < spec.size(); specIdx++) {
        if (matchMatrix[dbSpecIdx][specIdx] != -(float)INT_MAX) {
          outputRow = true;
          break;
        }
      }
      if (outputRow) {
        cout << dbSpecIdx << "\t";
        for (int specIdx = 0; specIdx < spec.size(); specIdx++) {
          cout << matchType[dbSpecIdx][specIdx] << "\t";
        }
        cout << endl;
      }
    }
    return;
  }

  //--------------------------------------------------------------------------
  AlignmentPenaltyBased::AlignmentPenaltyBased(int             maxSpecGap,
                                               PenaltyMatrix * modPenaltyMatrix,
                                               PenaltyMatrix * blossumPenaltyMatrix,
                                               float           penaltyAlpha,
                                               float           penaltyBeta,
                                               float           maxMod,
                                               float           minMod)
  : m_maxSpecGap(maxSpecGap + 1),
    m_modPenaltyMatrix(modPenaltyMatrix),
    m_blossumPenaltyMatrix(blossumPenaltyMatrix),
    m_penaltyAlpha(penaltyAlpha),
    m_penaltyBeta(penaltyBeta),
    m_maxMod(maxMod),
    m_minMod(minMod)
  {
    if (m_modPenaltyMatrix != 0x0 || m_modPenaltyMatrix != 0x0) {
      m_mapGapVectors.clear();
      createGapVectors();
    }
    m_cacheHitsTotal = 0;
    m_cacheMissesTotal = 0;
  }

  //--------------------------------------------------------------------------
  AlignmentPenaltyBased::~AlignmentPenaltyBased()
  {
    // Nothing here
  }

  //--------------------------------------------------------------------------
  void AlignmentPenaltyBased::clearCache(void)
  {
    m_mapCachedGapVectors.clear();
  }

  //--------------------------------------------------------------------------
  bool AlignmentPenaltyBased::isCached(string & aaString)
  {
    return m_mapCachedGapVectors.find(aaString) != m_mapCachedGapVectors.end();
  }
  
  //--------------------------------------------------------------------------
  void AlignmentPenaltyBased::sortStringLetters(string & toSort)
  {
    // Initialize the vector
    for (int i = 0; i < toSort.length() - 1; i++) {
      for (int j = i; j < toSort.length(); j++) {
        if (toSort[i] > toSort[j]) {
          char temp = toSort[i];
          toSort[i] = toSort[j];
          toSort[j] = temp;
        }
      }
    }
  }
  
  //---------------------------------------------------------------------------------------
  // This method adds  to the cache all the maxLength kMers that make up the protein string
  //---------------------------------------------------------------------------------------
  void AlignmentPenaltyBased::cacheProteinStrings(string & proteinString, int maxLength)
  {
    if (DEBUG_CACHE) DEBUG_VAR(proteinString);
    m_mapCachedGapVectors.clear();
    vector<float> scoresCombined;
    for (int iLength = maxLength - m_kMer + 1; iLength <= maxLength; iLength++) {
      if (DEBUG_CACHE) DEBUG_VAR(iLength);
      for (int iStart = 0; iStart < proteinString.length() - iLength; iStart++) {
	
	if (DEBUG_CACHE) DEBUG_VAR(iStart);
	string full = proteinString.substr(iStart, iLength);
	sortStringLetters(full);
	if (DEBUG_CACHE) DEBUG_VAR(full);
        if (isCached(full)) {
	  continue;
	}
	
	string left = full.substr(0, m_kMer);
	sortStringLetters(left);
	if (DEBUG_CACHE) DEBUG_VAR(left);
	vector<float> scoresLeft = m_mapGapVectors[left];

	string right = full.substr(m_kMer);
	sortStringLetters(right);
	if (DEBUG_CACHE) DEBUG_VAR(right);
	vector<float> scoresRight = m_mapGapVectors[left];
	
        initializeGapVector(scoresCombined);
        combineVectors(scoresCombined, scoresLeft, scoresRight);
	m_mapCachedGapVectors[full] = scoresCombined;
      }
    }
    if (DEBUG_CACHE) DEBUG_MSG("Done cacheing protein strings");
  }

  //--------------------------------------------------------------------------
  void AlignmentPenaltyBased::computeAlignment(Spectrum & 	spec,
                                               Spectrum & 	dbSpec,
                                               char * 		dbSeq,
                                               int 		dbIndex,
                                               int 			    matchOrientation,
                                               set<float> & startRange,
                                               int 		minMatchedPeaks,
                                               int 		maxGapSize,
                                               float 		pmTolerance,
                                               float 		tolerance,
                                               bool 		enforceEndpeaks)
  {
    DEBUG_MSG("START PENALTY ALIGNMENT");
    if (m_mapCachedGapVectors.size() > (float)CACHE_LIMIT * 0.8) {
      clearCache();
    }
    m_cacheHits = 0;
    m_cacheMisses = 0;
    
    if (spec.size() == 0) {
      WARN_MSG("Spectrum size is 0");
      return;
    }

    // Compute average peak intensity
    float avgPeakIntensity = 0;
    for (int i = 0; i < spec.size(); i++) {
      avgPeakIntensity += spec[i][1];
    }
    avgPeakIntensity /= spec.size();

#if DEBUG_ALIGN
    DEBUG_VAR(enforceEndpeaks);
    DEBUG_VAR(maxGapSize);
    DEBUG_VAR(startRange.size());
    DEBUG_VAR(pmTolerance);
    DEBUG_VAR(tolerance);
#endif
    if (DEBUG_ALIGN) DEBUG_VAR(avgPeakIntensity);
    if (DEBUG_AAS) debugAAs(m_modPenaltyMatrix, m_blossumPenaltyMatrix);
    if (DEBUG_SPECS) debugDBSpec(dbSeq, dbSpec);
    if (DEBUG_SPECS) debugSpec(spec);

    vector<vector< float > >         matchMatrix;
    vector<vector< pair<int,int> > > matchPtr;
    vector<vector< string > >        matchType;
    vector<vector<char> >            startFlags(spec.size());

    // Initialize the starting flags
    setStartingFlagArray(spec, dbSpec, startRange, startFlags, tolerance);

    matchMatrix.resize(dbSpec.size());
    matchPtr.resize(dbSpec.size());
    matchType.resize(dbSpec.size());

    char dbGapString[strlen(dbSeq) + 1];
    
    //---------------------------------
    //  Initialize the matrices
    //---------------------------------
    for (int dbSpecIdx = 0; dbSpecIdx < dbSpec.size(); dbSpecIdx++) {
      matchMatrix[dbSpecIdx].resize(spec.size());
      matchPtr[dbSpecIdx].resize(spec.size());
      matchType[dbSpecIdx].resize(spec.size());
      for (int specIdx = 0; specIdx < spec.size(); specIdx++) {
        matchMatrix[dbSpecIdx][specIdx] = -(float)INT_MAX;
        matchPtr[dbSpecIdx][specIdx] = make_pair<int,int>(-1,-1);
        matchType[dbSpecIdx][specIdx] = "?";
      }
    }

    //---------------------------------
    // Compute the match between the two spectra
    //---------------------------------
    for (int dbSpecIdx = 0; dbSpecIdx < dbSpec.size(); dbSpecIdx++) {

      for (int specIdx = 0; specIdx < spec.size(); specIdx++) {
    
        if (DEBUG_RANGE2) DEBUG_MSG(dbSpecIdx << ", " << specIdx << ", " << (int)startFlags[specIdx][dbSpecIdx]);

        if (!enforceEndpeaks or specIdx == 0) {
          // Matching can start on any mass, assume no predecessor 
          matchMatrix[dbSpecIdx][specIdx] = spec[specIdx][1];
        }

        if (DEBUG_RANGE2) DEBUG_MSG(matchMatrix[dbSpecIdx][specIdx]);

        for (int predDbSpecIdx = dbSpecIdx - 1;
             (predDbSpecIdx >= 0);
              predDbSpecIdx--) {

          for (int predSpecIdx = specIdx - 1;
        	   (predSpecIdx >= 0);
        	   predSpecIdx--) {

            if (!startFlags[predSpecIdx][predDbSpecIdx]) {
              // if location is not a valid start, skip scoring
              if (DEBUG_RANGE2) DEBUG_MSG("Skipping " << predDbSpecIdx << ", " << predSpecIdx);
              continue;
            }
            if (matchMatrix[predDbSpecIdx][predSpecIdx] == -(float)INT_MAX) {
              // Must build on a real score (can make anything out of -infinity)
              if (DEBUG_RANGE2) DEBUG_MSG("Skipping " << predDbSpecIdx << ", " << predSpecIdx);
              continue;
            }
        
            //makeDbString(dbGapString, dbSeq, predDbSpecIdx, dbSpecIdx);
            int dbGapLength = dbSpecIdx - predDbSpecIdx;//strlen(dbGapString);
            //int dbGapLength = dbSpecIdx - predDbSpecIdx;
            float dbGapMass = dbSpec[dbSpecIdx][0] - dbSpec[predDbSpecIdx][0];
            float dbAAMinMass = (dbSpecIdx - predDbSpecIdx) * 57.0;
            float specGapMass = spec[specIdx][0] - spec[predSpecIdx][0];

            if (DEBUG_ALIGN) DEBUG_MSG(predDbSpecIdx << "," << dbSpecIdx << "  " << predSpecIdx << "," << specIdx);
            if (DEBUG_ALIGN) DEBUG_MSG(dbGapLength << "  " << dbGapMass << "  " << specGapMass << "  " << dbAAMinMass);
            
            float deltaDbMass = dbSpec[dbSpecIdx][0] - dbSpec[predDbSpecIdx][0];
            float deltaSpecMass = spec[specIdx][0] - spec[predSpecIdx][0];
            float deltaMasses = deltaSpecMass - deltaDbMass;

            if (dbSpecIdx - predDbSpecIdx == 1) {

              if (deltaDbMass > 400.0) continue;
              if (deltaSpecMass > 400.0) continue;
            
              makeDbString(dbGapString, dbSeq, predDbSpecIdx, dbSpecIdx);
              string strAA(dbGapString);

              //if (DEBUG_ALIGN) DEBUG_MSG(strAA << "  " << deltaDbMass << "  " << deltaSpecMass << "  " << deltaMasses);

              if (deltaMasses < m_minMod || deltaMasses > m_maxMod) {
                continue;
              }

              if (abs(deltaMasses) < tolerance) {
              
                if (DEBUG_ALIGN) DEBUG_MSG("Exact Match");
                if (DEBUG_ALIGN) DEBUG_VAR(spec[specIdx][1]);
                float score = matchMatrix[predDbSpecIdx][predSpecIdx] +
                                      spec[specIdx][1];
                if (DEBUG_ALIGN) DEBUG_VAR(score);
                if (score > matchMatrix[dbSpecIdx][specIdx]) {
                  matchMatrix[dbSpecIdx][specIdx] = score;
                  if (DEBUG_ALIGN2) DEBUG_VAR(matchMatrix[predDbSpecIdx][predSpecIdx]);
                  if (DEBUG_ALIGN2) DEBUG_VAR(matchMatrix[dbSpecIdx][specIdx]);
                  matchPtr[dbSpecIdx][specIdx] = make_pair<int,int>(predDbSpecIdx,predSpecIdx);
                  matchType[dbSpecIdx][specIdx] = strAA;
                  if (DEBUG_ALIGN2) DEBUG_VAR(matchType[dbSpecIdx][specIdx]);
                }

              } else if (predSpecIdx == 0 && m_modPenaltyMatrix->isNterm(deltaMasses)) {
              
                if (DEBUG_ALIGN) DEBUG_MSG("Nterm Modification Match");
                // LARS: Figure out penalty for Nterm later
                float penaltyMod = -1.0;
                if (DEBUG_ALIGN) DEBUG_VAR(penaltyMod);
                float score = matchMatrix[predDbSpecIdx][predSpecIdx] +
                                        spec[specIdx][1] + penaltyMod;
                if (DEBUG_ALIGN) DEBUG_VAR(score);
                if (score > matchMatrix[dbSpecIdx][specIdx]) {
                  matchMatrix[dbSpecIdx][specIdx] = score;
                  if (DEBUG_ALIGN2) DEBUG_VAR(matchMatrix[predDbSpecIdx][predSpecIdx]);
                  if (DEBUG_ALIGN2) DEBUG_VAR(matchMatrix[dbSpecIdx][specIdx]);
                  matchPtr[dbSpecIdx][specIdx] = make_pair<int,int>(predDbSpecIdx,predSpecIdx);
                  matchType[dbSpecIdx][specIdx] = "[" + massToString(deltaMasses) + "]" + strAA;
                  if (DEBUG_ALIGN2) DEBUG_VAR(matchType[dbSpecIdx][specIdx]);
                }
              } else if (m_modPenaltyMatrix->isInMatrix(strAA, deltaMasses)) {
              
                if (DEBUG_ALIGN) DEBUG_MSG("Modification Match");
                float penaltyMod = (*m_modPenaltyMatrix)(strAA, deltaMasses, avgPeakIntensity);
                if (DEBUG_ALIGN) DEBUG_VAR(penaltyMod);
                float score = matchMatrix[predDbSpecIdx][predSpecIdx] +
                                        spec[specIdx][1] + penaltyMod;
                if (DEBUG_ALIGN) DEBUG_VAR(score);
                if (score > matchMatrix[dbSpecIdx][specIdx]) {
                  matchMatrix[dbSpecIdx][specIdx] = score;
                  if (DEBUG_ALIGN2) DEBUG_VAR(matchMatrix[predDbSpecIdx][predSpecIdx]);
                  if (DEBUG_ALIGN2) DEBUG_VAR(matchMatrix[dbSpecIdx][specIdx]);
                  matchPtr[dbSpecIdx][specIdx] = make_pair<int,int>(predDbSpecIdx,predSpecIdx);
                  matchType[dbSpecIdx][specIdx] = "(" + strAA + "," + massToString(deltaMasses) + ")";
                  if (DEBUG_ALIGN2) DEBUG_VAR(matchType[dbSpecIdx][specIdx]);
                }
#if 0
              } else if (m_blossumPenaltyMatrix->isInMatrix(strAA, deltaMasses)) {
              
                if (DEBUG_ALIGN) DEBUG_MSG("Blossum Match");
                float penaltyBlossum = (*m_blossumPenaltyMatrix)(strAA, deltaMasses, avgPeakIntensity);
                if (DEBUG_ALIGN) DEBUG_VAR(penaltyBlossum);
                float score = matchMatrix[predDbSpecIdx][predSpecIdx] +
                                        spec[specIdx][1] +
                                        m_penaltyBeta * penaltyBlossum;
                if (DEBUG_ALIGN) DEBUG_VAR(score);
                if (score > matchMatrix[dbSpecIdx][specIdx]) {
                  matchMatrix[dbSpecIdx][specIdx] = score;
                  if (DEBUG_ALIGN2) DEBUG_VAR(matchMatrix[predDbSpecIdx][predSpecIdx]);
                  if (DEBUG_ALIGN2) DEBUG_VAR(matchMatrix[dbSpecIdx][specIdx]);
                  matchPtr[dbSpecIdx][specIdx] = make_pair<int,int>(predDbSpecIdx,predSpecIdx);
                  matchType[dbSpecIdx][specIdx] = "(" + strAA + "," + massToString(deltaMasses) + ")";
                  if (DEBUG_ALIGN2) DEBUG_VAR(matchType[dbSpecIdx][specIdx]);
                }
#endif
              } else if (!m_modPenaltyMatrix->isInMatrix(strAA, deltaMasses) &&
                         m_modPenaltyMatrix->getMass(strAA) + deltaMasses >= 57.0) {

                if (DEBUG_ALIGN) DEBUG_MSG("Unknown Match");
                float penaltyMod = m_modPenaltyMatrix->getUnknownPenalty(avgPeakIntensity) * m_penaltyAlpha;
                if (DEBUG_ALIGN) DEBUG_VAR(penaltyMod);
                float score = matchMatrix[predDbSpecIdx][predSpecIdx] +
                                  spec[specIdx][1] +
                                  penaltyMod;
                if (DEBUG_ALIGN) DEBUG_VAR(score);
                if (score > matchMatrix[dbSpecIdx][specIdx]) {
                  matchMatrix[dbSpecIdx][specIdx] = score;
                  if (DEBUG_ALIGN2) DEBUG_VAR(matchMatrix[predDbSpecIdx][predSpecIdx]);
                  if (DEBUG_ALIGN2) DEBUG_VAR(matchMatrix[dbSpecIdx][specIdx]);
                  matchPtr[dbSpecIdx][specIdx] = make_pair<int,int>(predDbSpecIdx,predSpecIdx);
                  matchType[dbSpecIdx][specIdx] = "(" + strAA + "," + massToString(deltaMasses) + ")";
                  if (DEBUG_ALIGN2) DEBUG_VAR(matchType[dbSpecIdx][specIdx]);
                }

              }  // else if (m_modPenaltyMatrix->isInMatrix(strAA, deltaMasses)) && ...

            } else if (abs(deltaMasses) < tolerance) {
              
              if (DEBUG_ALIGN_GAP) DEBUG_MSG("Exact Gap Match");
              //if (DEBUG_ALIGN_GAP) DEBUG_VAR(dbGapString);
              if (DEBUG_ALIGN_GAP) DEBUG_VAR(spec[specIdx][1]);
              float score = matchMatrix[predDbSpecIdx][predSpecIdx] +
                                      spec[specIdx][1];
              if (DEBUG_ALIGN_GAP) DEBUG_VAR(score);
              if (score > matchMatrix[dbSpecIdx][specIdx]) {
                matchMatrix[dbSpecIdx][specIdx] = score;
                if (DEBUG_ALIGN_GAP) DEBUG_VAR(matchMatrix[predDbSpecIdx][predSpecIdx]);
                if (DEBUG_ALIGN_GAP) DEBUG_VAR(matchMatrix[dbSpecIdx][specIdx]);
                matchPtr[dbSpecIdx][specIdx] = make_pair<int,int>(predDbSpecIdx,predSpecIdx);
                makeDbString(dbGapString, dbSeq, predDbSpecIdx, dbSpecIdx);
                matchType[dbSpecIdx][specIdx] = dbGapString;
                if (DEBUG_ALIGN_GAP) DEBUG_VAR(matchType[dbSpecIdx][specIdx]);
              }


            } else if (dbAAMinMass < specGapMass && 
                       (dbSpecIdx - predDbSpecIdx <= maxGapSize) && 
                       (spec[specIdx][0] - spec[predSpecIdx][0] < m_maxSpecGap - 1) &&
            		  specGapMass > dbGapMass + m_minMod * dbGapLength &&
            		  specGapMass < dbGapMass + m_maxMod * dbGapLength &&
                       matchMatrix[predDbSpecIdx][predSpecIdx] + spec[specIdx][1] > matchMatrix[dbSpecIdx][specIdx]) {
 
              if (DEBUG_ALIGN_NTERM) DEBUG_MSG("Gap Match");
              string matchString;
              float penalty = -(float)INT_MAX;
              int specGapLength = (int)(specGapMass + 0.5);
              if (DEBUG_ALIGN_NTERM) DEBUG_VAR(specGapLength);
              makeDbString(dbGapString, dbSeq, predDbSpecIdx, dbSpecIdx);
              penalty = getGapPenalty(specGapLength, dbGapString, avgPeakIntensity);
              if (DEBUG_ALIGN_NTERM) DEBUG_VAR(penalty);
              if (penalty == 0.0) {
                matchString = dbGapString;
                if (DEBUG_ALIGN_NTERM) DEBUG_VAR(matchString);
              } else {
                char buffer[20];
                sprintf(buffer, "%0.2f", deltaMasses);
                matchString = "(" + string(dbGapString) + "," + string(buffer) + ")";  // Gap annotation (place holder because we don't compute exact)
                if (DEBUG_ALIGN_NTERM) DEBUG_VAR(matchString);
              }

              float bestGapPenalty = penalty;
              if (DEBUG_ALIGN_NTERM) DEBUG_VAR(bestGapPenalty);

              // If we are at nterm then loop over all the nterm mods and try them with gap
              if (predSpecIdx == 0) {
                const set<float> & setNtermMods = m_modPenaltyMatrix->getNtermMods();
                if (DEBUG_ALIGN_NTERM) DEBUG_VAR(setNtermMods.size());
                set<float>::const_iterator itrs = setNtermMods.begin();
                set<float>::const_iterator itrsEnd = setNtermMods.end();
                for ( ; itrs != itrsEnd; itrs++) {
                  float ntermMod = *itrs;
                  if (DEBUG_ALIGN_NTERM) DEBUG_VAR(ntermMod);
                  float newSpecGapMass = specGapMass - ntermMod;
                  if (DEBUG_ALIGN_NTERM) DEBUG_VAR(newSpecGapMass);
                  if (DEBUG_ALIGN_NTERM) DEBUG_VAR(dbAAMinMass);
                  if (DEBUG_ALIGN_NTERM) DEBUG_VAR(dbGapMass);
                  if (DEBUG_ALIGN_NTERM) DEBUG_VAR(m_minMod * dbGapLength);
                  if (DEBUG_ALIGN_NTERM) DEBUG_VAR(m_maxMod * dbGapLength);
                  if (dbAAMinMass < newSpecGapMass && 
            	        newSpecGapMass > dbGapMass + m_minMod * dbGapLength &&
            	        newSpecGapMass < dbGapMass + m_maxMod * dbGapLength) {

                    if (DEBUG_ALIGN_NTERM) DEBUG_TRACE;
                    int specGapLength = (int)(newSpecGapMass  + 0.5);
                    if (DEBUG_ALIGN_NTERM) DEBUG_VAR(specGapLength);
                    float ntermPenalty = getGapPenalty(specGapLength, dbGapString, avgPeakIntensity);
                    if (ntermPenalty > bestGapPenalty) {
                      if (DEBUG_ALIGN_NTERM) DEBUG_VAR(ntermPenalty);
                      if (ntermPenalty == 0.0) {
                        matchString = "[" + massToString(ntermMod)+ "]" + dbGapString;
                        if (DEBUG_ALIGN_NTERM) DEBUG_VAR(matchString);
                      } else {
                        char buffer[20];
                        sprintf(buffer, "%0.2f", deltaMasses - ntermMod);
                        matchString = "[" + massToString(ntermMod)+ "]" + "(" + dbGapString + "," + string(buffer) + ")";
                        if (DEBUG_ALIGN_NTERM) DEBUG_VAR(matchString);
                      }
                      bestGapPenalty = ntermPenalty;
                      if (DEBUG_ALIGN_NTERM) DEBUG_VAR(bestGapPenalty);
                    }
                  }
                }
              }

              if (DEBUG_ALIGN_NTERM) DEBUG_VAR(bestGapPenalty);
              float score = matchMatrix[predDbSpecIdx][predSpecIdx] +
                                spec[specIdx][1] + bestGapPenalty;
              if (DEBUG_ALIGN_NTERM) DEBUG_VAR(score);
              if (score > matchMatrix[dbSpecIdx][specIdx]) {
                matchMatrix[dbSpecIdx][specIdx] = score;
                if (DEBUG_ALIGN_NTERM) DEBUG_VAR(matchMatrix[predDbSpecIdx][predSpecIdx]);
                if (DEBUG_ALIGN_NTERM) DEBUG_VAR(matchMatrix[dbSpecIdx][specIdx]);
                matchPtr[dbSpecIdx][specIdx] = make_pair<int,int>(predDbSpecIdx,predSpecIdx);
                matchType[dbSpecIdx][specIdx] = matchString;
                if (DEBUG_ALIGN_NTERM) DEBUG_VAR(matchType[dbSpecIdx][specIdx]);
              }

            }

          } // predSpecIdx

        } // predDbSpecIdx

      } // specIdxS

    } // dbSpecIdx

    if (DEBUG_ALIGN3) debugMatrices(spec, dbSpec, matchMatrix, matchType);
       
    //---------------------------------
    // Find best scoring match
    //---------------------------------
    float bestScore = -(float)INT_MAX;
    pair<int,int> bestMatchPtr;

    if (enforceEndpeaks) {
      for (int dbSpecIdx = 0; dbSpecIdx < dbSpec.size(); dbSpecIdx++) {
        if (matchMatrix[dbSpecIdx][spec.size()-1] > bestScore) {
          bestScore = matchMatrix[dbSpecIdx][spec.size()-1];
          bestMatchPtr = make_pair<int,int>(dbSpecIdx,spec.size()-1);
          if (DEBUG_ALIGN2) DEBUG_MSG(bestMatchPtr.first << ", " << bestMatchPtr.second);
        }
      }
    } else {
      for (int dbSpecIdx = 0; dbSpecIdx < dbSpec.size(); dbSpecIdx++) {
        for (int specIdx = 0; specIdx < spec.size(); specIdx++) {
          if (matchMatrix[dbSpecIdx][specIdx] > bestScore) {
            bestScore = matchMatrix[dbSpecIdx][specIdx];
            bestMatchPtr = make_pair<int,int>(dbSpecIdx,specIdx);
            if (DEBUG_ALIGN2) DEBUG_MSG(bestMatchPtr.first << ", " << bestMatchPtr.second);
          }
        }
      }
    }
    
    if (DEBUG_ALIGN2) DEBUG_MSG(bestMatchPtr.first << ", " << bestMatchPtr.second);

    //---------------------------------
    // Create PSM for the match
    //---------------------------------
    if (bestScore != -(float)INT_MAX) {

      // Create a PSM for this match      
      psmPtr p(new PeptideSpectrumMatch);
      p->m_score = bestScore;
      p->m_dbIndex = dbIndex;
      p->m_matchOrientation = matchOrientation;
      p->m_spectrum = &spec;

      int lastDbIndex = bestMatchPtr.first;
      
      if (DEBUG_ALIGN2) DEBUG_VAR(p->m_score);
      if (DEBUG_ALIGN2) DEBUG_VAR(p->m_matchOrientation);
      if (DEBUG_ALIGN2) DEBUG_VAR(bestMatchPtr.first);
      if (DEBUG_ALIGN2) DEBUG_VAR(bestMatchPtr.second);
      if (DEBUG_ALIGN2) DEBUG_VAR(matchType[bestMatchPtr.first][bestMatchPtr.second]);

      p->m_annotation.clear();
      // Check for gap at end
      if (bestMatchPtr.second != spec.size() - 1) {
        float endMass = spec[spec.size() - 1][0] - spec[bestMatchPtr.second][0];
        p->m_annotation = "[" + massToString(endMass) + "]";
      }

      pair<int,int> prevMatchPtr = bestMatchPtr;
      TwoValues<int> prevIndices;
      prevIndices[0] = -1;
      prevIndices[1] = -1;
      bool wasGap = false;
      while (bestMatchPtr.second != -1) {
        TwoValues<int> tmpIndices;
        tmpIndices[0] = bestMatchPtr.second;
        tmpIndices[1] = bestMatchPtr.first;
        p->m_matchedPeaks.push_back(tmpIndices);

        string matchAnnotation = matchType[bestMatchPtr.first][bestMatchPtr.second];

        if (matchAnnotation != "?") {
          p->m_annotation = matchAnnotation + p->m_annotation;
          prevIndices = tmpIndices;
	 }
        prevMatchPtr = bestMatchPtr;
        bestMatchPtr = matchPtr[bestMatchPtr.first][bestMatchPtr.second];

#if DEBUG_ALIGN2
        DEBUG_VAR(bestMatchPtr.first);
        DEBUG_VAR(bestMatchPtr.second);
        if (bestMatchPtr.first != -1 || bestMatchPtr.second != -1) {
          DEBUG_VAR(matchType[bestMatchPtr.first][bestMatchPtr.second]);
        }
#endif
      }
      int firstDbIndex = prevMatchPtr.first;
      p->m_startMass = dbSpec[firstDbIndex][0];

      // Check for gap at beginning
      if (prevMatchPtr.second != 0) {
        p->m_annotation = "[" + massToString(spec[prevMatchPtr.second][0]) + "]" + p->m_annotation;
        p->m_startMass -= spec[prevMatchPtr.second][0]; // Adjust start mass backward to account for gap
      }

      if (DEBUG_ALIGN2) DEBUG_VAR(firstDbIndex);
      if (DEBUG_ALIGN2) DEBUG_VAR(lastDbIndex);

      if (lastDbIndex - firstDbIndex > 2) {
        if (DEBUG_ALIGN2) DEBUG_VAR(p->m_annotation);
        if (DEBUG_ALIGN2) DEBUG_VAR(p->m_matchedPeaks.size());
        reverse(p->m_matchedPeaks.begin(),p->m_matchedPeaks.end());
        spec.psmList.push_back(p);
      } 

    }  // if (bestScore != -(float)INT_MAX)

    if (DEBUG_CACHE) DEBUG_VAR(m_cacheHits);
    if (DEBUG_CACHE) DEBUG_VAR(m_cacheMisses);
    m_cacheHitsTotal += m_cacheHits;
    m_cacheMissesTotal += m_cacheMisses;

    if (DEBUG_CACHE) DEBUG_VAR(m_cacheHitsTotal);
    if (DEBUG_CACHE) DEBUG_VAR(m_cacheMissesTotal);

    return;
  }

  //--------------------------------------------------------------------------
  float AlignmentPenaltyBased::getPeptideMass(string & stringPeptide)
  {
    float totalMass = 0.0;
    for (int i = 0; i < stringPeptide.length(); i++) {
      string strAA;
      strAA.append(1, stringPeptide[i]);
      float mass = m_modPenaltyMatrix->getMass(strAA) * 0.9995;
      totalMass += mass;
    }
    return totalMass;
  }

  //--------------------------------------------------------------------------
  void AlignmentPenaltyBased::computeAllGapAnnotations(string & stringAnnotationIn,
                                                       string & stringAnnotationOut)
  {
    DEBUG_VAR(stringAnnotationIn);
    stringAnnotationOut.clear();
    if (stringAnnotationIn.empty()) return;

    bool leftParen = false;
    bool comma = false;
    string gapString;
    string gapMass;
    for (int i = 0; i < stringAnnotationIn.length(); i++) {
      if (stringAnnotationIn[i] == ')') {
        DEBUG_VAR(gapString);
        DEBUG_VAR(gapMass);
        string annotation;
        float gapMassFloat;
        char buffer[20];
        sscanf(gapMass.c_str(), "%f", &gapMassFloat);
        DEBUG_VAR(gapMassFloat);
        getGapAnnotation(getPeptideMass(gapString) + gapMassFloat, gapString, annotation);
        DEBUG_VAR(annotation);
        stringAnnotationOut += annotation;
        leftParen = false;
        comma = false;
      } else if (stringAnnotationIn[i] == '(') {
        leftParen = true;
        gapString.clear();
        gapMass.clear();
      } else if (leftParen) {
        if (stringAnnotationIn[i] == ',') {
          comma = true;
        } else if (comma) {
          gapMass += stringAnnotationIn[i];
        } else {
          gapString += stringAnnotationIn[i];
        }
      } else {
        stringAnnotationOut += stringAnnotationIn[i];
      }
    }
    DEBUG_VAR(stringAnnotationOut );
    return;
  }

  //--------------------------------------------------------------------------
  void AlignmentPenaltyBased::getGapAnnotation(float    specGapLengthFloat,
                                               string & aaString,
                                               string & stringAnnotation)
  {
    if (DEBUG_GAP_ANNO) DEBUG_VAR(specGapLengthFloat);
    int specGapLength = (int)(specGapLengthFloat + 0.5);
    if (DEBUG_GAP_ANNO) DEBUG_VAR(specGapLength);
    if (aaString.empty()) {
      return;
    }
    string strAA;
    strAA.append(1, aaString[0]);
    if (DEBUG_GAP_ANNO) DEBUG_VAR(strAA);

    vector<float> scoresSoFar = m_mapGapVectors[strAA];

    int aaMass = (int)m_modPenaltyMatrix->getMass(strAA);
    if (DEBUG_GAP_ANNO) DEBUG_VAR(aaMass);

    vector<float> scoresCombined(m_maxSpecGap);
    initializeGapVector(scoresCombined);
    vector<string> vecAnnotations1(m_maxSpecGap);
    vecAnnotations1[aaMass] = strAA;
    vector<string> vecAnnotations2(m_maxSpecGap);

    // Initialize first annotation vector
    for (int i = 57; i < aaMass; i++) {
      vecAnnotations1[i] = makeAnnotation(strAA, i - aaMass);
    }
    for (int i = aaMass + 1; i < m_maxSpecGap; i++) {
      vecAnnotations1[i] = makeAnnotation(strAA, i - aaMass);
    }

#if 0
    for (int i = 0; i < m_maxSpecGap; i++) {
      DEBUG_MSG(i << "  " << vecAnnotations1[i]);
    }
#endif

    vector<string> & prevAnno = vecAnnotations1;
    vector<string> & newAnno = vecAnnotations2;
    // Add annotations for each additional character
    for (int iChar = 1; iChar < aaString.length(); iChar++) {
      string strNewAA;
      strNewAA.append(1, aaString[iChar]);
      //if (DEBUG_GAP_ANNO) DEBUG_VAR(strNewAA);
      int newAAMass = (int)m_modPenaltyMatrix->getMass(strNewAA);
      //if (DEBUG_GAP_ANNO) DEBUG_VAR(newAAMass);
      combineVectorsWithAnnotation(scoresSoFar, m_mapGapVectors[strNewAA], scoresCombined,
                                   vecAnnotations1, vecAnnotations2, strNewAA, newAAMass);
      scoresSoFar = scoresCombined;

      // swap the annotation vectors to avoid copies
      vector<string> & tempAnno = prevAnno;
      prevAnno = newAnno;
      newAnno = tempAnno;
    }

#if 0
    for (int i = 0; i < prevAnno.size(); i++) {
      DEBUG_MSG(i << "  " << prevAnno[i]);
    }
#endif

    if (specGapLength > prevAnno.size()) {
      WARN_MSG("Gap size too large!");
      // This should never happen.. but just in case
      stringAnnotation = "(" + aaString +",x.xxx)";
      return;
    }
    stringAnnotation = prevAnno[specGapLength];
    return;
  }

  //--------------------------------------------------------------------------
  float AlignmentPenaltyBased::getGapPenalty(int      specGapLength,
                                             string   aaString,
                                             float    avgPeakIntensity)
  {
    if (DEBUG_GAP) DEBUG_VAR(specGapLength);
    if (DEBUG_GAP) DEBUG_VAR(aaString);
    if (aaString.empty()) {
      WARN_MSG("Empty string passed to getGapPenalty!")
      return -(float)INT_MAX;
    }

    // If the length of the string is small enough then we have the answer in map
    string sortedString = aaString;
    sortStringLetters(sortedString);
    if (aaString.length() <= m_kMer) {

      // See if this string is in the map
      if (m_mapGapVectors.find(sortedString) == m_mapGapVectors.end()) {
        return -(float)INT_MAX;
      }

      return m_mapGapVectors[sortedString][specGapLength] * avgPeakIntensity;
    }

    if (isCached(sortedString)) {
      //if (DEBUG_CACHE) DEBUG_MSG("Retrieved answer for [" << sortedString << "] from the cache");
      m_cacheHits++;
      return m_mapCachedGapVectors[sortedString][specGapLength] * avgPeakIntensity;
    }
    //if (DEBUG_CACHE) DEBUG_MSG("Answer for [" << sortedString << "] not in cache");
    m_cacheMisses++;

    vector<float> scoresCombined(m_maxSpecGap);
    initializeGapVector(scoresCombined);

    string left = sortedString.substr(0, m_kMer);
    sortStringLetters(left);
    if (DEBUG_GAP) DEBUG_VAR(left);
    vector<float> scoresSoFar = m_mapGapVectors[left];

    string remaining = sortedString.substr(m_kMer);
    if (DEBUG_GAP) DEBUG_VAR(remaining);
    while (!remaining.empty()) {
      string left = remaining.substr(0, m_kMer);
      //DEBUG_VAR(left);
      sortStringLetters(left);
      //DEBUG_VAR(left);

      combineVectors(scoresSoFar, m_mapGapVectors[left], scoresCombined);
      scoresSoFar = scoresCombined;

      if (remaining.size() > m_kMer) {
        remaining = remaining.substr(m_kMer);
      } else {
        remaining = "";
      }
      //DEBUG_VAR(remaining);
    }

    // Cache this string in case user needs it again
    m_mapCachedGapVectors[sortedString] = scoresSoFar;

    // Return the value
    return scoresSoFar[specGapLength] * avgPeakIntensity;
  }

  //--------------------------------------------------------------------------
  void AlignmentPenaltyBased::combineVectors(vector<float> & scores1,
                                             vector<float> & scores2,
                                             vector<float> & outputScores)
  {
    vector<float>::iterator itrOut = outputScores.begin();
    vector<float>::iterator itr1 = scores1.begin();
    for (int i = 0; i < scores1.size(); i++, itrOut++, itr1++) {
      vector<float>::iterator itrOut2 = itrOut;
      vector<float>::iterator itr2 = scores2.begin();
      for (int j = 0; j < scores2.size() - i; j++, itrOut2++, itr2++) {
        float sum = *itr1 + *itr2;
        if (*itrOut2 < sum) {
          *itrOut2 = sum;
        }
      } // for (int j = 0; j < scores2.size(); j++) {
    } // for (int i = 0; i < scores1.size(); i++) {

    return;
  }

  //--------------------------------------------------------------------------
  void AlignmentPenaltyBased::combineVectorsWithAnnotation(vector<float>  & scores1,
                                                           vector<float>  & scores2,
                                                           vector<float>  & outputScores,
                                                           vector<string> & annotationsOld,
                                                           vector<string> & annotationsNew,
                                                           string         & newAA,
                                                           int            & newAAMass)
  {
    vector<char> vecFlags(outputScores.size());
    vector<char>::iterator itrFlags = vecFlags.begin();
    vector<string>::iterator itrAnnoOld = annotationsOld.begin();
    vector<string>::iterator itrAnnoNew = annotationsNew.begin();
    vector<float>::iterator itrOut = outputScores.begin();
    vector<float>::iterator itr1 = scores1.begin();
    for (int i = 0; i < scores1.size(); i++, itrOut++, itr1++, itrFlags++, itrAnnoOld++, itrAnnoNew++) {
      vector<float>::iterator itrOut2 = itrOut;
      vector<float>::iterator itr2 = scores2.begin();
      vector<char>::iterator itrFlags2 = itrFlags;
      vector<string>::iterator itrAnnoNew2 = itrAnnoNew;
      for (int j = 0; j < scores2.size() - i; j++, itrOut2++, itr2++, itrFlags2++, itrAnnoNew2++) {
        float sum = *itr1 + *itr2;
        if ((*itrFlags2 == 0) || (*itrOut2 < sum)) {
          *itrFlags2 = 1;   //vecFlags[i + j] = 1;
          *itrOut2 = sum;   //outputScores[i + j] = scores1[i] + scores2[j];
          if (*itr2 == 0) {  // if scores2[j] == 0
            *itrAnnoNew2 = *itrAnnoOld + newAA; //annotationsNew[i + j] = annotationsOld[i] + newAA;
          } else {
            *itrAnnoNew2 = *itrAnnoOld + makeAnnotation(newAA, j - newAAMass);  //annotationsNew[i + j] = annotationsOld[i] + makeAnnotation(newAA, j - newAAMass);
          }
        }
      } // for (int j = 0; j < scores2.size(); j++) {

    } // for (int i = 0; i < scores1.size(); i++) {

    return;
  }

  //--------------------------------------------------------------------------
  void AlignmentPenaltyBased::createGapVectors(void)
  {
    // Set this so we know which vectors are pre-created
    //   kMer = 4 means all combinations of 4 AA's (sorted)
    m_kMer = 4;

    vector<string> aaVec;
    m_modPenaltyMatrix->getAminoAcids(aaVec);
    if (DEBUG_GAP) DEBUG_VAR(aaVec.size());
    sort(aaVec.begin(), aaVec.end());

    // Create all the single AA gap vectors
    for (int a1 = 0; a1 < aaVec.size(); a1++) {
      string strAA1 = aaVec[a1];
      //DEBUG_VAR(strAA1);
      vector<float> scores1(m_maxSpecGap);
      initializeGapVector(scores1);
      fillVectorSingle(scores1, strAA1);
      m_mapGapVectors[strAA1] = scores1;
    }


    // Create all the AA gap vectors (up to 3mers for right now)
    for (int a1 = 0; a1 < aaVec.size(); a1++) {

      string strAA1 = aaVec[a1];
      //DEBUG_VAR(strAA1);

      for (int a2 = 0; a2 < aaVec.size(); a2++) {

        string strAA2 = aaVec[a2];
        string strAACombined2 = strAA1 + strAA2;
        sortStringLetters(strAACombined2);

        // Check to see if we've already done this one
        if (m_mapGapVectors.find(strAACombined2) != m_mapGapVectors.end()) {
          //DEBUG_MSG("Already Done");
          continue;
        }
        //DEBUG_VAR(strAACombined2);

        vector<float> scoresCombined2(m_maxSpecGap);
        initializeGapVector(scoresCombined2);

        combineVectors(m_mapGapVectors[strAA1], m_mapGapVectors[strAA2], scoresCombined2);
        m_mapGapVectors[strAACombined2] = scoresCombined2;

        for (int a3 = 0; a3 < aaVec.size(); a3++) {

          string strAA3 = aaVec[a3];
          string strAACombined3 = strAACombined2 + strAA3;
          sortStringLetters(strAACombined3);

          // Check to see if we've already done this one
          if (m_mapGapVectors.find(strAACombined3) != m_mapGapVectors.end()) {
            //DEBUG_MSG("Already Done");
            continue;
          }
          //DEBUG_VAR(strAACombined3);

          vector<float> scoresCombined3(m_maxSpecGap);
          initializeGapVector(scoresCombined3);
          combineVectors(m_mapGapVectors[strAACombined2], m_mapGapVectors[strAA3], scoresCombined3);
          m_mapGapVectors[strAACombined3] = scoresCombined3;
#if 1
          for (int a4 = 0; a4 < aaVec.size(); a4++) {

            string strAA4 = aaVec[a4];
            string strAACombined4 = strAACombined3 + strAA4;
            sortStringLetters(strAACombined4);

            // Check to see if we've already done this one
            if (m_mapGapVectors.find(strAACombined4) != m_mapGapVectors.end()) {
              //DEBUG_MSG("Already Done");
              continue;
            }
            //DEBUG_VAR(strAACombined4);

            vector<float> scoresCombined4(m_maxSpecGap);
            initializeGapVector(scoresCombined4);
            combineVectors(m_mapGapVectors[strAACombined3], m_mapGapVectors[strAA4], scoresCombined4);
            m_mapGapVectors[strAACombined4] = scoresCombined4;
          }
#endif
        }  // for (int a3 = 0; a3 < aaVec.size(); a3++) {

      }  // for (int a2 = 0; a2 < aaVec.size(); a2++) {

    } // for (int a1 = 0; a1 < aaVec.size(); a1++) {

  }

  //--------------------------------------------------------------------------
  void AlignmentPenaltyBased::fillVectorSingle(vector<float> & scores,
                                              string           strAA)
  {
    //DEBUG_VAR(strAA);
    // No penalty for hitting DB mass exactly
    float aaMass = m_modPenaltyMatrix->getMass(strAA);
    //DEBUG_VAR(aaMass);
    size_t index = (size_t)aaMass;
    scores[index] = 0.0;

    std::map<float, float> penalties;
    m_modPenaltyMatrix->getPenalties(strAA, penalties, 1.0);
    std::map<float, float>::iterator itrm = penalties.begin();
    std::map<float, float>::iterator itrm_end = penalties.end();
    // Loop over all the possible mods (and their penalties)
    for (;itrm != itrm_end; itrm++) {
      float massDiff = itrm->first;
      float penalty = itrm->second;
      if (penalty == 0.0) {
        continue;
      }
      size_t indexPenalty = (size_t)(index + massDiff);
      // Penalties are negative so a "higher" penalty is better
      float newPenalty = penalty;
      if (!m_modPenaltyMatrix->isKnown(strAA, massDiff)) { // Known mods don't get multiplied by alpha
        newPenalty *= m_penaltyAlpha;
      }
      if (indexPenalty >= 0 &&
          indexPenalty < scores.size()) {
        scores[indexPenalty] = newPenalty;
      }
    }

    // Fill in all possible unknown modifications
    int minRealDelta = (int)(m_modPenaltyMatrix->getMass(strAA) - 57.0);
    int startUnk = max((int)(index - minRealDelta), (int)(index + m_minMod)); // min mod is negative
    startUnk = max(0, startUnk);  // Don't go beyond beginning of array
    int endUnk = min((int)(index + m_maxMod), (int)scores.size());
    float unkPenalty = m_modPenaltyMatrix->getUnknownPenalty(m_penaltyAlpha); // multiply by alpha
    for (int iUnk = startUnk; iUnk < endUnk; iUnk++) {
      float newPenalty = unkPenalty;
      if (newPenalty > scores[iUnk]) {
        scores[iUnk] = newPenalty;
      }
    }

    return;
  }

  //--------------------------------------------------------------------------
  void AlignmentPenaltyBased::initializeGapVector(vector<float> & scores)
  {
    // Initialize the vector
    for (int i = 0; i < scores.size(); i++) {
      scores[i] = -(float)INT_MAX;
    }
  }

  //--------------------------------------------------------------------------
  string AlignmentPenaltyBased::makeAnnotation(string & strAA, int mass)
  {
    char strMod[100];
    sprintf(strMod, "%.3f", (float)mass);
    return "(" + strAA + "," + strMod + ")";
  }

  //--------------------------------------------------------------------------
  void AlignmentPenaltyBased::makeDbString(char * dbGapStringOut, char * dbSeq, int start, int end)
  {
    int length = end - start;
    strncpy(dbGapStringOut, dbSeq + start, length);
    dbGapStringOut[length] = '\0';
    return;
  }

  //--------------------------------------------------------------------------
  string AlignmentPenaltyBased::massToString(float mass)
  {
    char strMod[100];
    sprintf(strMod, "%.3f", mass);
    return string(strMod);
  }

  //--------------------------------------------------------------------------
  void AlignmentPenaltyBased::setStartingFlagArray(Spectrum &              spec,
                                                   Spectrum &              dbSpec,
                                                   set<float> &            startRange,
                                                   vector<vector<char> > & startFlags,
                                                   float                   tolerance)
  {
    //---------------------------------
    //  INITIALIZE THE STARTING FLAGS
    //---------------------------------
    //setStartingFlagArray(spec, dbSpec, startRange, startFlags, tolerance);

    for (int i = 0; i < spec.size(); i++) {
      vector<char> newArray(dbSpec.size());
      startFlags[i] = newArray;
      for (int j = 0; j < dbSpec.size(); j++) {
	    // If the startRange set is empty then allow any start position
	    if (startRange.size() == 0) {
	      startFlags[i][j] = 1;
	    } else {
	      startFlags[i][j] = 0;
	    }
      }
    }

    if (DEBUG_RANGE) DEBUG_VAR(m_minMod);
    if (DEBUG_RANGE) DEBUG_VAR(m_maxMod);
    // Find all the valid starting points in the matrix
    set<float>::iterator itr = startRange.begin();
    set<float>::iterator itrEnd = startRange.end();
    for (; itr != itrEnd; itr++) {
      int minIdx1, maxIdx1, minIdx2, maxIdx2;
      if (DEBUG_RANGE) DEBUG_VAR(*itr);
      float minStartMass = *itr + m_minMod;
      float maxStartMass = *itr + m_maxMod;
      if (DEBUG_RANGE) DEBUG_VAR(minStartMass);
      if (DEBUG_RANGE) DEBUG_VAR(maxStartMass);
      float lastMass = 0;
      for (int specIdx = 0; specIdx < spec.size(); specIdx++) {
	    float massDiff = spec[specIdx][0] - lastMass;
	    if (DEBUG_RANGE) DEBUG_VAR(massDiff);
	    minStartMass += massDiff;
	    maxStartMass += massDiff;
	    if (DEBUG_RANGE) DEBUG_VAR(minStartMass);
	    if (DEBUG_RANGE) DEBUG_VAR(maxStartMass);
	    lastMass = spec[specIdx][0];
	    // Find peaks close to desired bounds
	    float T = 200.0;

	    list<int> matches1;
	    list<int> matches2;
	    dbSpec.setPeakTolerance(0);
	    dbSpec.findPeaks(minStartMass, T, &matches1);
	    dbSpec.findPeaks(maxStartMass, T, &matches2);
	    if (DEBUG_RANGE) DEBUG_VAR(matches1.size());
	    if (DEBUG_RANGE) DEBUG_VAR(matches1.front());
	    if (DEBUG_RANGE) DEBUG_VAR(matches1.back());
	    if (DEBUG_RANGE) DEBUG_VAR(matches2.size());
	    if (DEBUG_RANGE) DEBUG_VAR(matches2.front());
	    if (DEBUG_RANGE) DEBUG_VAR(matches2.back());
	    minIdx1 = (matches1.size() == 0) ? -1 : matches1.front();
	    maxIdx1 = (matches1.size() == 0) ? -1 : matches1.back();
	    minIdx2 = (matches2.size() == 0) ? -1 : matches2.front();
	    maxIdx2 = (matches2.size() == 0) ? -1 : matches2.back();

	    // If peak is not found, index could be -1.. so make sure at least 0
	    minIdx1 = max<int>(minIdx1,0);
	    minIdx2 = max<int>(minIdx2,0);
	    // Make sure min peak is INSIDE the desired range
	    while (minIdx1 < dbSpec.size() &&
			dbSpec[minIdx1][0] < minStartMass - tolerance - AAJumps::massH2O) minIdx1++;

	    // Make sure max peak is OUTSIDE (above) the desired range
	    while (minIdx2 < dbSpec.size() &&
			dbSpec[minIdx2][0] < maxStartMass + tolerance + AAJumps::massH2O) minIdx2++;

	    // Mark all peaks in range as valid starting points in matrix
	    if (DEBUG_RANGE) DEBUG_VAR(minIdx1);
	    if (DEBUG_RANGE) DEBUG_VAR(minIdx2);
	    for (int i = minIdx1; i < minIdx2 && i < dbSpec.size(); i++) {
	      startFlags[specIdx][i] = 1;
	    }
      } // for (int specIdx = 0; specIdx < spec.size(); specIdx++)
    } // for (; itr != itrEnd; itr++) {

#if DEBUG_RANGE
    for (int j = 0; j < dbSpec.size(); j++) {
      cout << j << "\t";
      for (int i = 0; i < spec.size(); i++) {
	cout << (int)startFlags[i][j] << "\t";
      }
      cout << endl;
    }
#endif
    return;
  }


} // namespace specnets
