/*
 * GFTable.h
 *
 *  Created on: Jul 19, 2013
 *      Author: aguthals
 */

#ifndef GFTABLE_H_
#define GFTABLE_H_

#include "Cluster.h"
#include "spectrum.h"

#include <vector>
#include <list>


using namespace std;

namespace specnets
{

  struct GFCell
  {
    // Cummulative probability of all peptides ending at this node
    double probability;

    // Does this node patch a the prefix-residue mass of a valid peptide? (0/1)
    char isValid;
  };

  class GFTable : public SpectrumItem
  {

  protected:

    // ith mass -> jth score -> GFCell
    vector<vector<GFCell> > m_table;

    // ith mass -> score of node @ mass i
    vector<int> m_nodeScores;

    // Lower mass range of parent mass
    unsigned int m_firstPM;

    // Integer-valued amino acid masses loaded from AAJumps
    vector<int> m_aaMasses;

  public:

    // Default constructor
    GFTable();

    /** Constructs graph from a spectrum
     * @param spectrum
     * @param jumps set of AA masses to impose edges
     **/
    GFTable(const Spectrum &spectrum, const AAJumps &jumps);

    /** Constructs graph from a spectrum
     * @param spectrum
     * @param jumps set of AA masses to impose edges
     **/
    void initialize(const Spectrum &spectrum, const AAJumps &jumps);

    // Access the vector of GFNodes at a given mass
    inline const vector<GFCell> &operator[](unsigned int mass) const
    {
      return m_table[mass];
    }

    // Access the vector of GFNodes at a given mass
    inline vector<GFCell> &operator[](unsigned int mass)
    {
      return m_table[mass];
    }

    inline const vector<int> &getAAMasses() const
    {
      return m_aaMasses;
    }

    // Access the score at a given mass
    inline unsigned int getScore(unsigned int mass) const
    {
      return m_nodeScores[mass];
    }

    // Get the number of mass bins encoded in the graph
    inline unsigned int size() const
    {
      return m_table.size();
    }

    inline const unsigned int & getPmLowerBound() const
    {
      return m_firstPM;
    }

    /** Removes all nodes/edges that do not encode peptides that have probability within the given threshold.
     * @param probThreshold
     **/
    void encodeDictionary(double probThreshold);

  private:

    /**
     * Backtracks m_table from a given cell and marks any node found as valid in the input table. If a
     *   valid node is encountered, the exploration ends
     * @param mass
     * @param score
     * @param validNodes
     */
    void exploreBackwards(const unsigned int &mass,
                          const unsigned int &score,
                          vector<vector<char> > &validNodes);
  };
}

#endif /* GFTABLE_H_ */
