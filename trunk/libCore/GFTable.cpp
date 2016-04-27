/*
 * GFTable.cpp
 *
 *  Created on: Jul 19, 2013
 *      Author: aguthals
 */

#include "GFTable.h"
#include "aminoacid.h"
#include "Logger.h"

using namespace std;

namespace specnets
{

  GFTable::GFTable() :
    SpectrumItem(), m_table(0), m_nodeScores(0), m_firstPM(0), m_aaMasses(0)
  {
  }

  GFTable::GFTable(const Spectrum &spectrum, const AAJumps &jumps) :
    SpectrumItem(), m_table(spectrum.size()), m_nodeScores(spectrum.size()),
        m_firstPM(0), m_aaMasses(0)
  {
    initialize(spectrum, jumps);
  }

  void GFTable::initialize(const Spectrum &spectrum, const AAJumps &jumps)
  {
    if (spectrum.parentMass < AAJumps::massMH)
    {
      WARN_MSG("bad spectrum with parent mass " << spectrum.parentMass
          << ", skipping GF calculation");
      return;
    }
    if (spectrum.size() == 0)
    {
      m_table.resize(0);
      m_nodeScores.resize(0);
      m_firstPM = 0;
      m_aaMasses.resize(0);
      return;
    }

    m_scan = spectrum.scan;
    m_filename = spectrum.fileName;
    m_index = -1;

    // compute size of node and edge vectors so we can lookup the parent mass at the end each
    unsigned int massBinsSize = round(spectrum.parentMass
        + spectrum.parentMassTol - AAJumps::massMH) + 1;

    m_firstPM = round(spectrum.parentMass - spectrum.parentMassTol
        - AAJumps::massMH);

    DEBUG_VAR(massBinsSize);
    DEBUG_VAR(m_firstPM);

    // initialize vector sizes
    m_nodeScores.assign(massBinsSize, 0);
    m_table.resize(massBinsSize);
    unsigned int totalScore = 0;

    for (unsigned int i = 0; i < spectrum.size(); i++)
    {
      unsigned int roundedMass = round(spectrum[i][0]);
      unsigned int nodeScore = round(spectrum[i][1]);
      totalScore += nodeScore;

      // If we have higher resolution we need to spread the score by peak tolerance
      unsigned int lowerMass = round(spectrum[i][0] - spectrum.getTolerance(i));
      unsigned int upperMass = round(spectrum[i][0] + spectrum.getTolerance(i));
      for (unsigned int j = lowerMass; j <= upperMass; j++)
      {
        m_nodeScores[j] = max(nodeScore, (unsigned int&)m_nodeScores[j]);
      }
    }

    DEBUG_TRACE;

    GFCell defaultNode;
    defaultNode.probability = 0;
    defaultNode.isValid = (char)0;

    if (m_nodeScores[0] == 0)
    {
      m_nodeScores[0] = 1;
    }

    DEBUG_TRACE;

    m_aaMasses.resize(jumps.size());
    for (unsigned int aa = 0; aa < jumps.size(); aa++)
    {
      // Load all rounded AA masses
      m_aaMasses[aa] = round(jumps[aa]);
    }

    unsigned int runningTotalScore = 0;
    int prevScore, currentNodeScore, prevMass;
    double totalProb;

    double aaProb = 1.0 / (double)jumps.size();

    for (int m = 0; m < massBinsSize; m++)
    {
      currentNodeScore = m_nodeScores[m];

      // Each 2nd degree vector only needs to be as large as the total score in the spectrum to this peak mass
      runningTotalScore += currentNodeScore;

      if (m < m_firstPM)
      {
        m_table[m].assign(runningTotalScore + 1, defaultNode);
      }
      else
      {
        m_table[m].assign(totalScore + 1, defaultNode);
      }

      if (m == 0)
      {
        // Initialize base case of the recursion
        m_table[m][0].probability = 1.0;
        m_table[m][0].isValid = (char)1;
      }
      for (int s = 0; s < m_table[m].size(); s++)
      {
        prevScore = s - currentNodeScore;
        if (prevScore < 0)
        {
          continue;
        }
        totalProb = 0;
        for (int aa = 0; aa < m_aaMasses.size(); aa++)
        {
          prevMass = m - m_aaMasses[aa];
          if (prevMass < 0 || (m_table[prevMass][prevScore].isValid == (char)0))
          {
            continue;
          }
          // add up the probability
          totalProb += m_table[prevMass][prevScore].probability * aaProb;
        }
        m_table[m][s].probability = totalProb;
        m_table[m][s].isValid = (char)1;
      }
    }
  }

  void GFTable::encodeDictionary(double probThreshold)
  {
    double totalSpecProb = 0;
    unsigned int maxScore = m_table[m_firstPM].size() - 1;
    unsigned int cutoffScore = 0;

    // Find the maximum score that has a total probability >= to the threshold
    for (unsigned int s = maxScore; s >= 0; s--)
    {
      cutoffScore = s;
      for (unsigned int m = m_firstPM; m < m_table.size(); m++)
      {
        totalSpecProb += m_table[m][s].probability;
      }
      if (totalSpecProb >= probThreshold)
      {
        break;
      }
    }

    DEBUG_VAR(totalSpecProb);

    // Mass m -> Score s -> Are peptides ending at (m,s) prefixes of at least one peptide in the dictionary? (0/1)
    vector<vector<char> > validNodes(m_table.size());
    for (unsigned int m = m_firstPM; m < m_table.size(); m++)
    {
      // Initialize all nodes as invalid
      validNodes[m].assign(m_table[m].size(), (char)0);
      for (unsigned int s = maxScore; s >= cutoffScore; s--)
      {
        // Set all nodes in paths ending at this node as valid
        exploreBackwards(m, s, validNodes);
      }
    }

    if (totalSpecProb > 0 && validNodes[0].size() == 0)
    {
      ERROR_MSG("Failed to backtrack to origin!!! (should not happen)");
      return;
    }

    GFCell invalidNode;
    invalidNode.probability = 0;
    invalidNode.isValid = (char)0;

    // For any node that we marked invalid, reflect that invalidity in our class data structure
    for (unsigned int m = 0; m < m_table.size(); m++)
    {
      if (validNodes[m].size() == 0)
      {
        m_table[m].assign(m_table[m].size(), invalidNode);
        continue;
      }

      for (unsigned int s = 0; s < m_table[m].size(); s++)
      {
        if (validNodes[m][s] == (char)0)
        {
          m_table[m][s] = invalidNode;
        }
      }
    }

  }

  void GFTable::exploreBackwards(const unsigned int &mass,
                                 const unsigned int &score,
                                 vector<vector<char> > &validNodes)
  {
    list<pair<unsigned int, unsigned int> > nextNodes;

    // Initialize the queue if we are starting from a reachable node
    if (m_table[mass][score].isValid == (char)1)
    {
      validNodes[mass][score] = (char)1;
      nextNodes.push_back(pair<unsigned int, unsigned int> (mass, score));
    }

    unsigned int nextMass;
    unsigned int nextScore;
    int m, s;
    while (nextNodes.size() > 0)
    {
      // explore by breadth-first search (so we pop from the front)
      nextMass = nextNodes.front().first;
      nextScore = nextNodes.front().second;
      nextNodes.pop_front();

      for (unsigned int aa = 0; aa < m_aaMasses.size(); aa++)
      {
        m = nextMass - m_aaMasses[aa];
        s = nextScore - m_nodeScores[nextMass];
        if (m < 0 || s < 0 || m_table[m][s].isValid == (char)0)
        {// skip nodes that don't exist in the GF tables
          continue;
        }

        if (validNodes[m].size() > 0 && validNodes[m][s] == (char)1)
        {// skip nodes we already marked as valid
          continue;
        }
        else if (validNodes[m].size() == 0)
        {// sometimes we reach an un-explored mass and we need to initialize the column
          validNodes[m].assign(m_table[m].size(), (char)0);
        }
        // Upon reaching an invalid node, mark it as valid and keep exploring from it.
        validNodes[m][s] = (char)1;
        nextNodes.push_back(pair<unsigned int, unsigned int> (m, s));
      }
    }
  }

}
