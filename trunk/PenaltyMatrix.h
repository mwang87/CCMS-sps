#ifndef __PenaltyMatrix_H__
#define __PenaltyMatrix_H__

#include "aminoacid.h"

// System Includes
#include <map>
#include <set>
#include <string>
#include <vector>


namespace specnets
{
  /*! \brief Holds a penalty matrix.

   A penalty matrix is a two-dimensional matrix indexed by Amino Acid Sequence
   and mass <string, float>.
  */
  class PenaltyMatrix
  {
  public:
    //! \name CONSTRUCTORS
    //@{
    /*! \brief The default constructor.
    */
    PenaltyMatrix(AAJumps & jumps, 
                  float resolution = 1.0, 
                  float unknownPenalty = 1.0, 
                  float unknownMultiplier = 2.0);

    //! \name DESTRUCTOR
    //@{
    virtual ~PenaltyMatrix(void);
    //@}

    //! \name ACCESSORS
    //@{
    bool isInMatrix(string & seq, float mass);
    bool isKnown(string & seq, float mass);
    bool isNterm(float mass);
    float operator()(string & seq, float mass, float averagePeakIntensity = 0.0);

    void  getPenalties(string & seq, std::map<float, float> & penaltyMap, float averagePeakIntensity = 0.0);
    float getMass(string aa);
    float getUnknownPenalty(float averagePeakIntensity = 0.0);
    const map<string, set<float> > & getKnownMods(void);
    const set<float> & getNtermMods(void);

    bool saveMatrix(std::string & filename);
    bool saveKnownMods(string & filename);
    bool saveAminoAcids(string & filename);

    void getAminoAcids(vector<string> & aaVec);
    //@}

    //! \name MODIFIERS
    //@{
    bool load(std::string & filename);
    bool loadAminoAcids(string & filename);

    bool loadFromBlosum(std::string & filename, float peakEquivalents);

    bool createFromModificationFreqs(map<float, float> & modFreq,
                                     float peakEquivalents,
                                     float minFrequency,
                                     float averagePeakIntensity = 1.0);

    bool loadKnownModifications(std::string & filename);
    //@}

  protected:

    std::map<std::string, std::map<float, float> > penalties;
    std::map<std::string, float> mapCharMods;
    std::map<float, std::string> mapModChars;
    float m_resolution;
    float m_unknownPenalty;
    float m_unknownMultiplier;
    float m_allSpectraAveragePeakIntensity;
    
  private:
    //! \name NOT IMPLEMENTED
    //@{
    PenaltyMatrix(const PenaltyMatrix & that);
    //@}

    float roundMass(float mass);

    map<string, set<float> > m_knownMods;
    set<float> m_knownNtermMods;
  };

} // namespace specnets

#endif // __PenaltyMatrix_H__

