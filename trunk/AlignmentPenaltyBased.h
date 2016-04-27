#ifndef ALIGNMENT_PENALTY_BASED2
#define ALIGNMENT_PENALTY_BASED2

#include "spectrum.h"
#include "PenaltyMatrix.h"
#include <set>

namespace specnets
{
  using namespace std;

  class AlignmentPenaltyBased
  {
  public:

    AlignmentPenaltyBased(int             maxSpecGap,
                          PenaltyMatrix * modPenaltyMatrix,
                          PenaltyMatrix * blossumPenaltyMatrix,
                          float           penaltyAlpha,
                          float           penaltyBeta,
                          float           maxMod,
                          float           minMod);
    ~AlignmentPenaltyBased();

    void clearCache(void);
    bool isCached(string & aaString);
    void cacheProteinStrings(string & proteinString, int maxLength);
    void sortStringLetters(string & toSort);

    void computeAlignment(Spectrum &	spec,
                          Spectrum &	dbSpec,
                          char * 		dbSeq,
                          int 		dbIndex,
                          int 		matchOrientation,
                          set<float> & 	startRange,
                          int 		minMatchedPeaks,
                          int 		maxGapSize,
                          float 		pmTolerance,
                          float 		tolerance,
                          bool 		enforceEndpeaks = false);

    void computeAllGapAnnotations(string & stringAnnotationIn,
                                  string & stringAnnotationOut);

    void getGapAnnotation(float    specGapLengthFloat,
                          string & aaString,
                          string & stringAnnotation);

    float getGapPenalty(int      specGapLength,
                        string   aaString,
                        float    avgPeakIntensity);

    int 			      m_cacheHits;
    int 			      m_cacheMisses;
    int 			      m_cacheHitsTotal;
    int 			      m_cacheMissesTotal;
  private:
    int 			      m_kMer;
    int 			      m_maxSpecGap;
    PenaltyMatrix * m_modPenaltyMatrix;
    PenaltyMatrix * m_blossumPenaltyMatrix;
    float           m_penaltyAlpha;
    float           m_penaltyBeta;
    float           m_maxMod;
    float           m_minMod;
    map<string, vector<float> > m_mapGapVectors;
    map<string, vector<float> > m_mapCachedGapVectors;

    void createGapVectors();

    void combineVectors(vector<float> & scores1,
                        vector<float> & scores2,
                        vector<float> & outputScores);

    void combineVectorsWithAnnotation(vector<float>  & scores1,
                                      vector<float>  & scores2,
                                      vector<float>  & outputScores,
                                      vector<string> & annotationsOld,
                                      vector<string> & annotationsNew,
                                      string         & newAA,
                                      int            & newAAMass);

    void fillVectorSingle(vector<float> & scores,
                          string          strAA);

    float getPeptideMass(string & stringPeptide);

    void initializeGapVector(vector<float> & scores);

    string makeAnnotation(string & strAA, int mass);

    void makeDbString(char * dbGapStringOut, char * dbSeq, int start, int end);

    string massToString(float mass);

    void setStartingFlagArray(Spectrum &   			      spec,
                              Spectrum &  			      dbSpec,
                              set<float> & 			      startRange,
                              vector<vector<char> > & startFlags,
                              float                   tolerance);
  };

}

#endif

