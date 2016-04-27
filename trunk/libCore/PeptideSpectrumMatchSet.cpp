#include "PeptideSpectrumMatchSet.h"

namespace specnets
{
  // -------------------------------------------------------------------------
  psmPtr & PeptideSpectrumMatchSet::operator[](unsigned int i)
  {
    return m_psmSet[i];
  }

  // -------------------------------------------------------------------------
  const psmPtr & PeptideSpectrumMatchSet::operator[](unsigned int i) const
  {
    return m_psmSet[i];
  }

  // -------------------------------------------------------------------------
  void PeptideSpectrumMatchSet::removePsmSetItem(psmPtr item)
  {
    for (int i = 0; i < m_psmSet.size(); i++)
    {

      std::tr1::shared_ptr<PeptideSpectrumMatch> p =
          std::tr1::dynamic_pointer_cast<PeptideSpectrumMatch>(m_psmSet[i]); //I AM A PSMNETWORK POINTER

      // 	  		PeptideSpectrumMatch * p = (PeptideSpectrumMatch *)m_psmSet[i].get();

      if (p == item)
      {
        m_psmSet.erase(m_psmSet.begin() + i, m_psmSet.begin() + i + 1);
        return;
      }
    }
    DEBUG_MSG("PeptideSpectrumMatchSet::removePsmSetItem:: pointer not found!");
  }
  // -------------------------------------------------------------------------
  PeptideSpectrumMatchSet & PeptideSpectrumMatchSet::operator=(PeptideSpectrumMatchSet &other)
  {
    if (this == &other)
    {
      return *this;
    }
    m_psmSet.clear();
    for (int i = 0; i < other.m_psmSet.size(); i++)
    {
      psmPtr currMatch(new PeptideSpectrumMatch);
      *currMatch = *(other.m_psmSet[i]);
      m_psmSet.push_back(currMatch);
    }
    return *this;
  }
  // -------------------------------------------------------------------------
  unsigned int PeptideSpectrumMatchSet::push_back(const psmPtr &other)
  {
    psmPtr currMatch(new PeptideSpectrumMatch);
    *currMatch = *other;
    m_psmSet.push_back(currMatch);
    return m_psmSet.size();
  }

  // -------------------------------------------------------------------------
  unsigned int PeptideSpectrumMatchSet::resize(unsigned int newSize)
  {
    m_psmSet.resize(newSize);
    return m_psmSet.size();
  }
  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatchSet::parseLine(psmPtr currMatch,
                                          vector<string> &line,
                                          int zeroIndexed,
                                          int scanIndex,
                                          int specIdxIndex,
                                          int spectrumFileIndex,
                                          int annotationIndex,
                                          int origAnnotationIndex,
                                          int proteinNameIndex,
                                          int dbIndexIndex,
                                          int numModsIndex,
                                          int matchOrientationIndex,
                                          int startMassIndex,
                                          int chargeIndex,
                                          int scoreIndex,
                                          int pvalueIndex,
                                          int strictenvIndex,
                                          int unstrictenvIndex,
                                          int compoundNameIndex,
                                          int fileScanUniqueIDIndex,
                                          int fdrIndex,
                                          bool isInspect,
                                          int decoyIndex,
                                          int organismIndex,
                                          int libraryNameIndex,
                                          int smilesIndex,
                                          int inchiIndex,
                                          int inchiAuxIndex,
                                          int libmetadataIndex,
                                          int sharedpeakIndex,
                                          int ppmerrIndex,
                                          int abundanceIndex,
                                          int FDRIndex,
                                          int massdiffIndex,
                                          int spectrumidIndex)
  {

    bool haveScanNum = false;

    if (scanIndex != -1)
    {
      if (!sscanf(line[scanIndex].c_str(), "%d", &(currMatch->m_scanNum)))
      {
        ERROR_MSG("Unable to get scan number!");
        return false;
      }
      if (currMatch->m_scanNum >= 0)
      {
        haveScanNum = true;
      }
    }

    if (!haveScanNum && specIdxIndex != -1)
    {
      if (!sscanf(line[specIdxIndex].c_str(), "%d", &(currMatch->m_scanNum)))
      {
        ERROR_MSG("Unable to get scan number from spectrum index!");
        return false;
      }
    }

    if (currMatch->m_scanNum >= 0 && zeroIndexed)
    {
      currMatch->m_scanNum += 1;
    }

    if (spectrumFileIndex != -1)
    {
      currMatch->m_spectrumFile = line[spectrumFileIndex];
    }

    if (annotationIndex != -1)
    {
      if (isInspect)
      {
        currMatch->inspectToSpecNets(line[annotationIndex],
                                     currMatch->m_annotation);
      }
      else
      {
        currMatch->m_annotation = line[annotationIndex];
      }
    }

    if (origAnnotationIndex != -1)
    {
      if (isInspect)
      {
        currMatch->inspectToSpecNets(line[origAnnotationIndex],
                                     currMatch->m_origAnnotation);
      }
      else
      {
        currMatch->m_origAnnotation = line[origAnnotationIndex];
      }
    }

    if (proteinNameIndex != -1)
    {
      currMatch->m_protein = line[proteinNameIndex];
      if (currMatch->m_protein.substr(0, 3).compare("XXX") == 0)
      {
        currMatch->m_isDecoy = true;
      }
      else
      {
        currMatch->m_isDecoy = false;
      }
    }

    if (dbIndexIndex != -1)
    {
      if (!sscanf(line[dbIndexIndex].c_str(), "%d", &(currMatch->m_dbIndex)))
      {
        ERROR_MSG("Unable to get DB index!");
        return false;
      }
    }
    if (numModsIndex != -1)
    {
      if (!sscanf(line[numModsIndex].c_str(), "%d", &(currMatch->m_numMods)))
      {
        ERROR_MSG("Unable to get num mods!");
        return false;
      }
    }

    if (matchOrientationIndex != -1)
    {
      if (!sscanf(line[matchOrientationIndex].c_str(),
                  "%d",
                  &(currMatch->m_matchOrientation)))
      {
        ERROR_MSG("Unable to get match orientation!");
        return false;
      }
    }

    if (startMassIndex != -1)
    {
      if (!sscanf(line[startMassIndex].c_str(), "%f", &(currMatch->m_startMass)))
      {
        ERROR_MSG("Unable to get start mass!");
        return false;
      }
    }
    if (chargeIndex != -1)
    {
      if (!sscanf(line[chargeIndex].c_str(), "%d", &(currMatch->m_charge)))
      {
        ERROR_MSG("Unable to get charge!");
        return false;
      }
    }
    if (scoreIndex != -1)
    {
      if (!sscanf(line[scoreIndex].c_str(), "%f", &(currMatch->m_score)))
      {
        ERROR_MSG("Unable to get score!");
        return false;
      }
    }
    if (pvalueIndex != -1)
    {
      float temp;
      if (!sscanf(line[pvalueIndex].c_str(), "%f", &(temp)))
      {
        ERROR_MSG("Unable to get p-value");
        return false;
      }
      currMatch->m_pValue = temp;
    }
    if (fdrIndex != -1)
    {
      if (!sscanf(line[fdrIndex].c_str(), "%f", &(currMatch->m_fdr)))
      {
        ERROR_MSG("Unable to get FDR index!");
        return false;
      }
    }

    if (strictenvIndex != -1)
    {
      if (!sscanf(line[strictenvIndex].c_str(),
                  "%f",
                  &(currMatch->m_strict_envelope_score)))
      {
        ERROR_MSG("Unable to get strict envelope score");
        return false;
      }
    }
    if (unstrictenvIndex != -1)
    {
      if (!sscanf(line[unstrictenvIndex].c_str(),
                  "%f",
                  &(currMatch->m_unstrict_envelope_score)))
      {
        ERROR_MSG("Unable to get unstrict envelope score");
        return false;
      }
    }

    if (decoyIndex != -1)
    {
      if (!sscanf(line[decoyIndex].c_str(), "%i", &(currMatch->m_isDecoy)))
      {
        ERROR_MSG("Unable to get unstrict envelope score");
        return false;
      }
    }

    if (organismIndex != -1)
    {
      currMatch->m_organism = line[organismIndex];
    }

    if (libraryNameIndex != -1)
    {
      currMatch->m_library_name = line[libraryNameIndex];
    }

    if (compoundNameIndex != -1)
    {
      currMatch->m_compound_name = (line[compoundNameIndex]);
    }

    if (smilesIndex != -1)
    {
      currMatch->m_smiles = (line[smilesIndex]);
    }

    if (inchiIndex != -1)
    {
      currMatch->m_InChI = (line[inchiIndex]);
    }

    if (inchiAuxIndex != -1)
    {
      currMatch->m_InChI_Aux = (line[inchiAuxIndex]);
    }

    if (libmetadataIndex != -1)
    {
      currMatch->m_notes = (line[libmetadataIndex]);
    }

    if (sharedpeakIndex != -1)
    {
      if (!sscanf(line[sharedpeakIndex].c_str(),
                  "%d",
                  &(currMatch->m_shared_peaks)))
      {
        ERROR_MSG("Unable to get shared peak count!");
        return false;
      }
    }

    if (sharedpeakIndex != -1)
    {
      if (!sscanf(line[ppmerrIndex].c_str(), "%f", &(currMatch->m_mz_error_ppm)))
      {
        ERROR_MSG("Unable to get ppm error!");
        return false;
      }
    }

    if (abundanceIndex != -1)
    {
      if (!sscanf(line[abundanceIndex].c_str(), "%d", &(currMatch->m_abundance)))
      {
        ERROR_MSG("Unable to get ppm error!");
        return false;
      }
    }

    if (FDRIndex != -1)
    {
      if (!sscanf(line[FDRIndex].c_str(), "%f", &(currMatch->m_fdr)))
      {
        ERROR_MSG("Unable to get FDR!");
        return false;
      }
    }

    if (massdiffIndex != -1)
    {
      if (!sscanf(line[massdiffIndex].c_str(),
                  "%f",
                  &(currMatch->m_parentmass_difference)))
      {
        ERROR_MSG("Unable to get Mass Difference!");
        return false;
      }
    }
    
    if (spectrumidIndex != -1)
    {
      currMatch->m_spectrumID = (line[spectrumidIndex]);
    }
    
    

    return true;
  }

  // -------------------------------------------------------------------------
  void PeptideSpectrumMatchSet::getPSMSet(SpecSet * spectra)
  {
    m_psmSet.clear(); //empty any current psm annotations

    for (int i = 0; i < spectra->size(); i++)
    {
      Spectrum spectrum = (*spectra)[i];
      if (spectrum.psmList.size() > 0)
      {
        list<psmPtr>::iterator it;

        for (it = spectrum.psmList.begin(); it != spectrum.psmList.end(); it++)
        {
          psmPtr currMatch(new PeptideSpectrumMatch);
          (*currMatch) = **it;
          currMatch->m_scanNum = spectrum.scan;
          m_psmSet.push_back(currMatch);
        }
      }
    }
  }
  // -------------------------------------------------------------------------
  void PeptideSpectrumMatchSet::addSpectraByFilename(SpecSet * spectra,
                                                     bool addToSpectra)
  {
    vector<psmPtr>::iterator it;
    spectra->index();
    for (it = m_psmSet.begin(); it != m_psmSet.end(); it++)
    {
      Spectrum temp;
      temp.fileName = (**it).m_spectrumFile;
      temp.scan = (**it).m_scanNum;
      string key = temp.getUniqueID();
      Spectrum * spectrum = spectra->getIndex(key);
      if (spectrum == NULL)
      {
        WARN_MSG("Unable to find spectrum for key " << key);
      }
      else
      {
        if (addToSpectra)
        {
          spectrum->psmList.push_back(*it);
        }
        (*it)->m_spectrum = spectrum;
      }
    }
  }
  // -------------------------------------------------------------------------
  void PeptideSpectrumMatchSet::addSpectra(SpecSet * spectra, bool addToSpectra)
  {
    vector<psmPtr>::iterator it;
    DEBUG_TRACE;
    for (it = m_psmSet.begin(); it != m_psmSet.end(); it++)
    {
      int scanNum = (*it)->m_scanNum;
      Spectrum * spectrum = spectra->getScan(scanNum); //getScan returns null if spectrum is not found

      if (spectrum == NULL)
      {
        WARN_MSG("Unable to find spectrum for scan " << scanNum);
      }
      else
      {
        if (addToSpectra)
        {
          spectrum->psmList.push_back(*it);
        }
      }
      (*it)->m_spectrum = spectrum;
    }
  }
  // -------------------------------------------------------------------------
  void PeptideSpectrumMatchSet::addSpectra(SpecSet * spectra,
                                           string filename,
                                           bool addToSpectra)
  {
    vector<psmPtr>::iterator it;
    for (it = m_psmSet.begin(); it != m_psmSet.end(); it++)
    {
      string psmFile;
      extractFileName((*it)->m_spectrumFile, psmFile, false);

      if (filename.compare(psmFile) == 0)
      {
        int scanNum = (*it)->m_scanNum;

        Spectrum * spectrum = spectra->getScan(scanNum);
        if (spectrum == NULL)
        {
          WARN_MSG("Unable to find spectrum for scan " << scanNum);
        }
        else
        {
          if (addToSpectra)
          {
            spectrum->psmList.push_back(*it);
          }
        }
        (*it)->m_spectrum = spectrum;
      }
    }
  }

  int PeptideSpectrumMatchSet::cluster(map<int, list<int> >& clusterInfo,
                                       short mergeType)
  {
    WARN_MSG("Calling deprecated version of cluster(), assuming one spectrum file in PSM set...");
    tr1::unordered_set < string > filenames;
    for (unsigned int i = 0; i < m_psmSet.size(); i++)
    {
      FilenameManager mngr(m_psmSet[i]->m_spectrumFile);
      filenames.insert(mngr.filename);
    }
    if (filenames.size() > 1)
    {
      ERROR_MSG("Found multiple spectrum files in PSM set, must use non-deprecated cluster method");
      abort();
    }
    map<int, list<pair<int, string> > > newClustInfo;
    for (map<int, list<int> >::iterator clustIt = clusterInfo.begin(); clustIt
        != clusterInfo.end(); clustIt++)
    {
      list < pair<int, string> > newL;
      for (list<int>::const_iterator childIt = clustIt->second.begin(); childIt
          != clustIt->second.end(); childIt++)
      {
        newL.push_back(pair<int, string> (*childIt, *filenames.begin()));
      }
      newClustInfo[clustIt->first] = newL;
    }
    return cluster(newClustInfo, mergeType);
  }

  // -------------------------------------------------------------------------
  int PeptideSpectrumMatchSet::cluster(map<int, list<pair<int, string> > >& clusterInfo,
                                       short mergeType)
  {
    tr1::unordered_map < string, map<int, unsigned int> > scanToPSMIdx;
    for (unsigned int i = 0; i < m_psmSet.size(); i++)
    {
      FilenameManager mngr(m_psmSet[i]->m_spectrumFile);
      //DEBUG_MSG("Have PSM for (" << m_psmSet[i]->m_scanNum << ", \'" << mngr.filename << "\')");
      if (scanToPSMIdx.count(mngr.filename) > 0)
      {
        scanToPSMIdx[mngr.filename][m_psmSet[i]->m_scanNum] = i;
      }
      else
      {
        map<int, unsigned int> tempMap;
        tempMap[m_psmSet[i]->m_scanNum] = i;
        scanToPSMIdx[mngr.filename] = tempMap;
      }

    }

    //DEBUG_VAR(scanToPSMIdx.size());

    if (m_psmSet.size() < 2)
    {
      return 0;
    }

    vector<psmPtr> newPSMSet(0);
    set<int> clusteredScans;

    float numClusters = 0;
    float numAgreeClusters = 0;

    for (map<int, list<pair<int, string> > >::iterator clustIt =
        clusterInfo.begin(); clustIt != clusterInfo.end(); clustIt++)
    {
      int newScan = clustIt->first;

      //DEBUG_VAR(newScan);

      psmPtr clusteredPSM;
      bool foundPSM = false;
      map < string, pair<int, psmPtr> > peptideCounts;
      set < string > unModPeptides;
      if (mergeType == 0)
      {
        for (list<pair<int, string> >::const_iterator childIt =
            clustIt->second.begin(); childIt != clustIt->second.end(); childIt++)
        {
          FilenameManager mngr(childIt->second);
          if (scanToPSMIdx.count(mngr.filename) == 0
              || scanToPSMIdx[mngr.filename].count(childIt->first) == 0)
          {
            //DEBUG_MSG("Failed to locate (" << childIt->first << ", \'" << mngr.filename << "\')");
            continue;
          }
          foundPSM = true;
          clusteredScans.insert(childIt->first);
          unsigned int psmIdx = scanToPSMIdx[mngr.filename][childIt->first];
          psmPtr nextPSM(m_psmSet[psmIdx]);
          unModPeptides.insert(AAJumps::stripMods(nextPSM->m_annotation));

          /*
           DEBUG_VAR(childIt->first);
           DEBUG_VAR(mngr.filename);
           DEBUG_VAR(nextPSM->m_annotation);
           */

          if (peptideCounts.count(nextPSM->m_annotation) == 0)
          {
            pair<int, psmPtr> newPair(1, nextPSM);
            peptideCounts[nextPSM->m_annotation] = newPair;
          }
          else
          {
            psmPtr prevPSM(peptideCounts[nextPSM->m_annotation].second);
            if (prevPSM->m_pValue > nextPSM->m_pValue)
            {
              peptideCounts[nextPSM->m_annotation].second = nextPSM;
            }
            peptideCounts[nextPSM->m_annotation].first++;
          }
        }

        bool initialized = false;
        int countSoFar = 0;

        for (map<string, pair<int, psmPtr> >::const_iterator countIt =
            peptideCounts.begin(); countIt != peptideCounts.end(); countIt++)
        {

          if (!initialized)
          {
            clusteredPSM = countIt->second.second;
            countSoFar = countIt->second.first;
            initialized = true;
            continue;
          }

          if (countSoFar < countIt->second.first)
          {
            clusteredPSM = countIt->second.second;
            countSoFar = countIt->second.first;
            continue;
          }

          if (clusteredPSM->m_pValue > countIt->second.second->m_pValue)
          {
            clusteredPSM = countIt->second.second;
            countSoFar = countIt->second.first;
            continue;
          }
        }
      }
      if (!foundPSM)
      {
        continue;
      }
      /*
       if (clustIt->second.size() > 3)
       {
       numClusters += 1.0;
       if (unModPeptides.size() == 1)
       {
       numAgreeClusters += 1.0;
       }
       else
       {
       DEBUG_VAR(peptideCounts.size());
       for (map<string, pair<int, psmPtr> >::const_iterator countIt =
       peptideCounts.begin(); countIt != peptideCounts.end(); countIt++)
       {
       DEBUG_VAR(countIt->first);
       }
       }
       }
       */
      //DEBUG_VAR(clusteredPSM->m_annotation);
      clusteredPSM->m_scanNum = newScan;
      newPSMSet.push_back(clusteredPSM);
    }
    /*
     for (unsigned int i = 0; i < m_psmSet.size(); i++)
     {
     if (clusteredScans.count(m_psmSet[i]->m_scanNum) > 0)
     {
     continue;
     }
     newPSMSet.push_back(m_psmSet[i]);
     }
     */

    m_psmSet = newPSMSet;

    //DEBUG_MSG((100.0 * numAgreeClusters)/numClusters << " percent of all clusters have agreeing PSMs (over " << numClusters << " clusters)");

    return clusteredScans.size();
  }

  // -------------------------------------------------------------------------
  void PeptideSpectrumMatchSet::getModCounts(map<string, map<float, float> > & mapModCount,
                                             bool useOrig)
  {
    vector<psmPtr>::iterator it;
    for (it = m_psmSet.begin(); it != m_psmSet.end(); it++)
    {
      vector<float> modifications;
      vector<unsigned int> positions;
      vector<unsigned int> lengths;
      (*it)->getModificationsAndPositions(modifications, positions, lengths);
      
      string annotation;
      if (useOrig) {
        annotation = (*it)->m_origAnnotation;
      } else {
        annotation = (*it)->m_annotation;
      }
      
      string cleanAnnotation; 
      static string aminoAcids("ABCDEFGHIJKLMNOPQRSTUVWXYZ");
      for (int iChar = 0; iChar < annotation.length(); iChar++) {
	if (aminoAcids.find_first_of(annotation[iChar]) != string::npos) {
	  cleanAnnotation += annotation[iChar];
	}
      }
      
      for (int i = 0; i < positions.size(); i++) {
	float length = (float)lengths[i];
	float modValue = (float)modifications[i];
	float intModValue = 0.0;
	if (modValue < 0) {
	  intModValue = (float)int(modValue - 0.5);
	} else {
	  intModValue = (float)int(modValue + 0.5);
	}
	if (positions[i] == 0) {
	  mapModCount["nTerm"][intModValue] = 
	      mapModCount["nTerm"][intModValue] + 1.0;
	}
	
        for (int k = 0; k < lengths[i]; k++) {
	  int pos = positions[i] + k;
#if 0 // Leave these out for now	  
	  if (pos == 1) {
            mapModCount["nTerm"][intModValue] = 
	        mapModCount["nTerm"][intModValue] + 1.0 / length;
	  }
#endif	  
	  string modChar("X");
	  modChar[0] = cleanAnnotation[pos-1];
	  mapModCount[modChar][intModValue] = 
	      mapModCount[modChar][intModValue] + 1.0 / length;
        }
	
      } // for (int iChar = 0; iChar < annotation.length(); iChar++) {

    } // for (it = m_psmSet.begin(); it != m_psmSet.end(); it++)

    return;
  }

  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatchSet::saveModMatrix(const char * filename,
                                              bool useOrig)
  {
    ofstream ofs(filename, ios::out);
    if (!ofs)
    {
      return false;
    }

    set<float> setAllMods;

    map < string, map<float, float> > mapModCount;
    getModCounts(mapModCount, useOrig);

    map<string, map<float, float> >::iterator itrc1 = mapModCount.begin();
    map<string, map<float, float> >::iterator itrc_end1 = mapModCount.end();
    ofs << "\t";
    for (; itrc1 != itrc_end1; itrc1++)
    {
      ofs << itrc1->first << "\t";

      map<float, float>::iterator itrc2 = itrc1->second.begin();
      map<float, float>::iterator itrc_end2 = itrc1->second.end();
      for (; itrc2 != itrc_end2; itrc2++)
      {
        setAllMods.insert(itrc2->first);
      }
    }
    ofs << "Total" << endl;

    set<float>::iterator itrs = setAllMods.begin();
    set<float>::iterator itrs_end = setAllMods.end();
    for (; itrs != itrs_end; itrs++)
    {
      ofs << *itrs;
      float rowTotal = 0;
      map<string, map<float, float> >::iterator itrc1 = mapModCount.begin();
      map<string, map<float, float> >::iterator itrc_end1 = mapModCount.end();
      for (; itrc1 != itrc_end1; itrc1++)
      {
        ofs << "\t" << mapModCount[itrc1->first][*itrs];
        rowTotal += mapModCount[itrc1->first][*itrs];
      }
      ofs << "\t" << rowTotal << endl;
    }
    return true;
  }

  // -------------------------------------------------------------------------
  void PeptideSpectrumMatchSet::maximumParsimony(void)
  {
    DEBUG_VAR(m_psmSet.size());

    // Initialize assignment vector
    int maxProtIdx = 0;
    int maxSpecIdx = 0;
    // Find maximum protein and maximum spectrum index of a matched protein
    for (size_t i = 0; i < m_psmSet.size(); i++)
    {
      maxProtIdx = max(maxProtIdx, m_psmSet[i]->m_dbIndex);
      maxSpecIdx = max(maxSpecIdx, m_psmSet[i]->m_scanNum);
    }

    //DEBUG_VAR(maxProtIdx);
    //DEBUG_VAR(maxSpecIdx);

    vector<psmPtr> assignments(maxSpecIdx);
    for (int pepIdx = 0; pepIdx < assignments.size(); pepIdx++)
    {
      psmPtr p(new PeptideSpectrumMatch);
      assignments[pepIdx] = p;
      assignments[pepIdx]->m_dbIndex = -1;
    }

    // Sort all the PSMs into scan based lists
    vector < list<psmPtr> > scanPsmLists(maxSpecIdx);
    for (size_t i = 0; i < m_psmSet.size(); i++)
    {
      scanPsmLists[m_psmSet[i]->m_scanNum - 1].push_back(m_psmSet[i]);
    }

    vector < list<psmPtr> > pepsPerProt(maxProtIdx + 1);
    //DEBUG_VAR(pepsPerProt.size());

    unsigned int maxHits, maxHitsIdx;
    float maxHitsMass;
    do
    {
      maxHits = 0;
      maxHitsIdx = 0;

      // Clear the peptides per protein vector
      for (int protIdx = 0; protIdx < pepsPerProt.size(); protIdx++)
      {
        pepsPerProt[protIdx].clear();
      }

      // Create the lists of peptides per protein and compute the max
      for (int scan = 0; scan < scanPsmLists.size(); scan++)
      {
        for (list<psmPtr>::iterator iter = scanPsmLists[scan].begin(); iter
            != scanPsmLists[scan].end(); iter++)
        {
          int dbIndex = (*iter)->m_dbIndex;
          pepsPerProt[dbIndex].push_back(*iter);
          if (pepsPerProt[dbIndex].size() > maxHits)
          {
            maxHits = pepsPerProt[dbIndex].size();
            maxHitsIdx = dbIndex;
          }
        }
      }
      //DEBUG_VAR(maxHits);
      //DEBUG_VAR(maxHitsIdx);

      if (maxHits > 0)
      {
        for (list<psmPtr>::iterator iter = pepsPerProt[maxHitsIdx].begin(); iter
            != pepsPerProt[maxHitsIdx].end(); iter++)
        {
          scanPsmLists[(*iter)->m_scanNum - 1].clear();
          assignments[(*iter)->m_scanNum - 1] = (*iter);
        }
      }

    } while (maxHits > 0);

    m_psmSet.clear();

    // Copy the assignments back to the vector
    for (int pepIdx = 0; pepIdx < assignments.size(); pepIdx++)
    {
      if (assignments[pepIdx]->m_dbIndex != -1)
      {
        m_psmSet.push_back(assignments[pepIdx]);
      }
    }

    DEBUG_VAR(m_psmSet.size());

    return;
  }

  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatchSet::loadFromFile(const char * filename)
  {
    // Clear any previously existing PSMs
    m_psmSet.clear();

    vector < string > header;
    vector < vector<string> > lines;
    vector < string > requiredHeader;
    vector<int> requiredHeaderIndex;

    //set up ALL headers
    requiredHeader.push_back("#Scan#");
    requiredHeader.push_back("SpectrumFile");
    requiredHeader.push_back("Annotation");
    requiredHeader.push_back("OrigAnnotation");
    requiredHeader.push_back("Protein");
    requiredHeader.push_back("dbIndex");
    requiredHeader.push_back("numMods");
    requiredHeader.push_back("matchOrientation");
    requiredHeader.push_back("startMass");
    requiredHeader.push_back("Charge");
    requiredHeader.push_back("MQScore");
    requiredHeader.push_back("p-value");
    /*
     requiredHeader.push_back("isDecoy");
     requiredHeader.push_back("StrictEnvelopeScore");
     requiredHeader.push_back("UnstrictEvelopeScore");
     requiredHeader.push_back("CompoundName");
     requiredHeader.push_back("FileScanUniqueID");
     */

    if (!DelimitedTextReader::loadDelimitedFile(filename,
                                                "\t",
                                                "",
                                                header,
                                                lines,
                                                requiredHeader,
                                                requiredHeaderIndex))
    {
      ERROR_MSG("Unable to open results! " << filename);
      return false;
    }

    int scanIndex = -1;
    int spectrumFileIndex = -1;
    int annotationIndex = -1;
    int origAnnotationIndex = -1;
    int proteinIndex = -1;
    int dbIndex = -1;
    int numModsIndex = -1;
    int matchOrientationIndex = -1;
    int startMassIndex = -1;
    int chargeIndex = -1;
    int scoreIndex = -1;
    int pvalueIndex = -1;
    int strictenvIndex = -1;
    int unstrictenvIndex = -1;
    int compoundNameIndex = -1;
    int fileScanUniqueIDIndex = -1;
    int fdrIndex = -1;
    int decoyIndex = -1;
    int organismIndex = -1;
    int libraryNameIndex = -1;
    int inchiIndex = -1;
    ;
    int inchiAuxIndex = -1;
    int smilesIndex = -1;
    int metadataIndex = -1;
    int sharedpeakIndex = -1;
    int ppmerrIndex = -1;
    int abundanceIndex = -1;
    int FDRIndex = -1;
    int massdiffIndex = -1;
    int spectrumidIndex = -1;

    //map headers
    for (int i = 0; i < header.size(); i++)
    {
      if (header[i].compare("#Scan#") == 0)
      {
        scanIndex = i;
      }
      else if (header[i].compare("SpectrumFile") == 0)
      {
        spectrumFileIndex = i;
      }
      else if (header[i].compare("Annotation") == 0)
      {
        annotationIndex = i;
      }
      else if (header[i].compare("OrigAnnotation") == 0)
      {
        origAnnotationIndex = i;
      }
      else if (header[i].compare("Protein") == 0)
      {
        proteinIndex = i;
      }
      else if (header[i].compare("dbIndex") == 0)
      {
        dbIndex = i;
      }
      else if (header[i].compare("numMods") == 0)
      {
        numModsIndex = i;
      }
      else if (header[i].compare("matchOrientation") == 0)
      {
        matchOrientationIndex = i;
      }
      else if (header[i].compare("startMass") == 0)
      {
        startMassIndex = i;
      }
      else if (header[i].compare("Charge") == 0)
      {
        chargeIndex = i;
      }
      else if (header[i].compare("MQScore") == 0)
      {
        scoreIndex = i;
      }
      else if (header[i].compare("p-value") == 0)
      {
        pvalueIndex = i;
      }
      else if (header[i].compare("isDecoy") == 0)
      {
        decoyIndex = i;
      }
      else if (header[i].compare("ScrictEnvelope") == 0)
      {
        strictenvIndex = i;
      }
      else if (header[i].compare("UnstrictEnvelope") == 0)
      {
        unstrictenvIndex = i;
      }
      else if (header[i].compare("CompoundName") == 0)
      {
        compoundNameIndex = i;
      }
      else if (header[i].compare("Organism") == 0)
      {
        organismIndex = i;
      }
      else if (header[i].compare("FileScanUniqueID") == 0)
      {
        fileScanUniqueIDIndex = i;
      }
      else if (header[i].compare("FDR") == 0)
      {
        fdrIndex = i;
      }
      else if (header[i].compare("LibraryName") == 0)
      {
        libraryNameIndex = i;
      }
      else if (header[i].compare("LibMetaData") == 0)
      {
        metadataIndex = i;
      }
      else if (header[i].compare("Smiles") == 0)
      {
        smilesIndex = i;
      }
      else if (header[i].compare("Inchi") == 0)
      {
        inchiIndex = i;
      }
      else if (header[i].compare("LibSearchSharedPeaks") == 0)
      {
        sharedpeakIndex = i;
      }
      else if (header[i].compare("mzErrorPPM") == 0)
      {
        ppmerrIndex = i;
      }
      else if (header[i].compare("Abundance") == 0)
      {
        abundanceIndex = i;
      }
      else if (header[i].compare("StrictEnvelopeScore") == 0)
      {
        strictenvIndex = i;
      }
      else if (header[i].compare("UnstrictEvelopeScore") == 0)
      {
        unstrictenvIndex = i;
      }
      else if (header[i].compare("FDR") == 0)
      {
        FDRIndex = i;
      }
      else if (header[i].compare("ParentMassDiff") == 0)
      {
        massdiffIndex = i;
      }
      else if (header[i].compare("LibrarySpectrumID") == 0)
      {
        spectrumidIndex = i;
      }
    }

    //parse results into m_psmSet
    for (int i = 0; i < lines.size(); i++)
    {
      psmPtr currMatch(new PeptideSpectrumMatch);

      if (!parseLine(currMatch,
                     lines[i],
                     false,
                     scanIndex,
                     -1,
                     spectrumFileIndex,
                     annotationIndex,
                     origAnnotationIndex,
                     proteinIndex,
                     dbIndex,
                     numModsIndex,
                     matchOrientationIndex,
                     startMassIndex,
                     chargeIndex,
                     scoreIndex,
                     pvalueIndex,
                     strictenvIndex,
                     unstrictenvIndex,
                     compoundNameIndex,
                     fileScanUniqueIDIndex,
                     fdrIndex,
                     false,
                     decoyIndex,
                     organismIndex,
                     libraryNameIndex,
                     smilesIndex,
                     inchiIndex,
                     inchiAuxIndex,
                     metadataIndex,
                     sharedpeakIndex,
                     ppmerrIndex,
                     abundanceIndex,
                     FDRIndex,
                     massdiffIndex,
                     spectrumidIndex))
      {
        WARN_MSG("Unable to parse line" << i);
      }
      else
      {
        /*cout << currMatch->m_spectrumFile << "\t" << currMatch->m_scanNum << "\t"
         << currMatch->m_annotation << endl; */
        m_psmSet.push_back(currMatch);
      }
    }

    DEBUG_VAR(m_psmSet.size());
    return true;
  }

  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatchSet::loadFromFiles(const char * resultsFileList)
  {
    vector < vector<string> > fileList;
    map<string, unsigned int> fileListHeader;
    vector < string > requiredHeader;
    requiredHeader.push_back("Path");
    vector<int> requiredHeaderIndex;

    if (!DelimitedTextReader::loadDelimitedFile(resultsFileList,
                                                "\t",
                                                "#",
                                                fileListHeader,
                                                fileList,
                                                requiredHeader,
                                                requiredHeaderIndex))
    {
      ERROR_MSG("Unable to load files!" << resultsFileList);
      return false;
    }

    for (int i = 0; i < fileList.size(); i++)
    {
      //get path
      const char * path = fileList[i][requiredHeaderIndex[0]].c_str();
      DEBUG_VAR(path);

      if (!this->loadFromFile(path))
      {
        ERROR_MSG("Unable to load file! "
            << fileList[i][requiredHeaderIndex[0]]);
        return false;
      }
    }
    return true;
  }
  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatchSet::saveToFile(const char * filename,
                                           bool includeHeader)
  {
    ofstream ofs(filename);
    if (!ofs || !ofs.good())
    {
      ERROR_MSG("Unable to open file [" << filename << "]");
      return false;
    }

    if (includeHeader)
    {
      ofs
          << "#Scan#\tSpectrumFile\tAnnotation\tOrigAnnotation\tProtein\tdbIndex\tnumMods\tmatchOrientation\tstartMass\tCharge\tMQScore\tp-value\tisDecoy\tStrictEnvelopeScore\tUnstrictEvelopeScore\tCompoundName\tOrganism\tFileScanUniqueID\tFDR\tLibraryName\tmzErrorPPM\tLibMetaData\tSmiles\tInchi\tLibSearchSharedPeaks\tAbundance\tParentMassDiff\tSpecMZ\tExactMass\tLibrarySpectrumID"
          << endl;
    }

    for (int i = 0; i < size(); i++)
    {
      m_psmSet[i]->saveToFile(ofs);
    }
    return true;
  }

  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatchSet::loadSpecnetsResultsFile(const char * resultsFile,
                                                        bool zeroIndexed)
  {
    vector < string > header;
    vector < vector<string> > lines;
    vector < string > requiredHeader;
    vector<int> requiredHeaderIndex;

    //set up inspect headers
    requiredHeader.push_back("Scan#");
    requiredHeader.push_back("Annotation");
    requiredHeader.push_back("Charge");
    requiredHeader.push_back("Protein");

    if (!DelimitedTextReader::loadDelimitedFile(resultsFile,
                                                "\t",
                                                "",
                                                header,
                                                lines,
                                                requiredHeader,
                                                requiredHeaderIndex))
    {
      ERROR_MSG("Unable to open results! " << resultsFile);
      return false;
    }

    int scanIndex = -1;
    int spectrumFileIndex = -1;
    int annotationIndex = -1;
    int chargeIndex = -1;
    int scoreIndex = -1;
    int pvalueIndex = -1;
    int proteinIndex = -1;

    //map headers
    for (int i = 0; i < header.size(); i++)
    {
      if (header[i].compare("Scan#") == 0)
      {
        scanIndex = i;
      }
      else if (header[i].compare("#SpectrumFile") == 0)
      {
        spectrumFileIndex = i;
      }
      else if (header[i].compare("Annotation") == 0)
      {
        annotationIndex = i;
      }
      else if (header[i].compare("Charge") == 0)
      {
        chargeIndex = i;
      }
      else if (header[i].compare("MQScore") == 0)
      {
        scoreIndex = i;
      }
      else if (header[i].compare("p-value") == 0)
      {
        pvalueIndex = i;
      }
      else if (header[i].compare("Protein") == 0)
      {
        proteinIndex = i;
      }
    }

    //parse results into m_psmSet
    for (int i = 0; i < lines.size(); i++)
    {
      psmPtr currMatch(new PeptideSpectrumMatch);

      if (!parseLine(currMatch,
                     lines[i],
                     zeroIndexed,
                     scanIndex,
                     -1,
                     spectrumFileIndex,
                     annotationIndex,
                     -1,
                     proteinIndex,
                     -1,
                     -1,
                     -1,
                     -1,
                     chargeIndex,
                     scoreIndex,
                     pvalueIndex,
                     -1,
                     -1,
                     -1,
                     -1,
                     -1,
                     false,
                     -1))
      {
        WARN_MSG("Unable to parse line" << i);
      }
      else
      {
        /*cout << currMatch->m_spectrumFile << "\t" << currMatch->m_scanNum << "\t"
         << currMatch->m_annotation << endl; */
        m_psmSet.push_back(currMatch);
      }
    }

    return true;
  }
  // -------------------------------------------------------------------------

  bool PeptideSpectrumMatchSet::loadInspectResultsFiles(const char * resultsFileList)
  {
    vector < vector<string> > inspectFileList;
    map<string, unsigned int> inspectFileListHeader;
    vector < string > requiredHeader;
    requiredHeader.push_back("Path");
    vector<int> requiredHeaderIndex;

    if (!DelimitedTextReader::loadDelimitedFile(resultsFileList,
                                                "\t",
                                                "#",
                                                inspectFileListHeader,
                                                inspectFileList,
                                                requiredHeader,
                                                requiredHeaderIndex))
    {
      ERROR_MSG("Unable to load inspect files!" << resultsFileList);
      return false;
    }

    //set up column numbers for optional parameters.
    int isZeroIndexedColumn = -1;

    map<string, unsigned int>::iterator it;

    it = inspectFileListHeader.find("isZeroIndexed");

    if (it != inspectFileListHeader.end())
    {
      isZeroIndexedColumn = it->second;
    }

    for (int i = 0; i < inspectFileList.size(); i++)
    {
      //get path
      const char * path = inspectFileList[i][requiredHeaderIndex[0]].c_str();
      DEBUG_VAR(path);

      //make sure results are loaded with correct indices.
      if (isZeroIndexedColumn != -1)
      {
        //make sure that sscanf doesn't overrun. Set to int then cast to bool.
        int isZeroIndexed = 0;
        if (!sscanf(inspectFileList[i][isZeroIndexedColumn].c_str(),
                    "%d",
                    &isZeroIndexed))
        {
          WARN_MSG("Unable to convert isZeroIndexed to bool!");
        }

        isZeroIndexed = (bool)isZeroIndexed;

        if (!this->loadInspectResultsFile(path, isZeroIndexed))
        {
          ERROR_MSG("Unable to load inspect file! "
              << inspectFileList[i][requiredHeaderIndex[0]]);
          return false;
        }
      }
      else
      {
        if (!this->loadInspectResultsFile(path))
        {
          ERROR_MSG("Unable to load inspect file! "
              << inspectFileList[i][requiredHeaderIndex[0]]);
          return false;
        }
      }
    }
    return true;
  }

  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatchSet::loadInspectResultsFile(const char * resultsFile,
                                                       bool zeroIndexed)
  {
    vector < string > header;
    vector < vector<string> > lines;
    vector < string > requiredHeader;
    vector<int> requiredHeaderIndex;

    //set up inspect headers
    requiredHeader.push_back("Scan#");
    requiredHeader.push_back("Annotation");
    requiredHeader.push_back("Charge");
    requiredHeader.push_back("Protein");

    if (!DelimitedTextReader::loadDelimitedFile(resultsFile,
                                                "\t",
                                                "",
                                                header,
                                                lines,
                                                requiredHeader,
                                                requiredHeaderIndex))
    {
      ERROR_MSG("Unable to open results! " << resultsFile);
      return false;
    }

    int scanIndex = -1;
    int spectrumFileIndex = -1;
    int annotationIndex = -1;
    int chargeIndex = -1;
    int scoreIndex = -1;
    int pvalueIndex = -1;
    int proteinIndex = -1;

    //map headers
    for (int i = 0; i < header.size(); i++)
    {
      if (header[i].compare("Scan#") == 0)
      {
        scanIndex = i;
      }
      else if (header[i].compare("#SpectrumFile") == 0)
      {
        spectrumFileIndex = i;
      }
      else if (header[i].compare("Annotation") == 0)
      {
        annotationIndex = i;
      }
      else if (header[i].compare("Charge") == 0)
      {
        chargeIndex = i;
      }
      else if (header[i].compare("MQScore") == 0)
      {
        scoreIndex = i;
      }
      else if (header[i].compare("p-value") == 0)
      {
        pvalueIndex = i;
      }
      else if (header[i].compare("Protein") == 0)
      {
        proteinIndex = i;
      }
    }

    //parse results into m_psmSet
    for (int i = 0; i < lines.size(); i++)
    {
      psmPtr currMatch(new PeptideSpectrumMatch);

      if (!parseLine(currMatch,
                     lines[i],
                     zeroIndexed,
                     scanIndex,
                     -1,
                     spectrumFileIndex,
                     annotationIndex,
                     -1,
                     proteinIndex,
                     -1,
                     -1,
                     -1,
                     -1,
                     chargeIndex,
                     scoreIndex,
                     pvalueIndex,
                     -1,
                     -1,
                     -1,
                     -1,
                     -1,
                     true,
                     -1))
      {
        WARN_MSG("Unable to parse line" << i);
      }
      else
      {
        /*cout << currMatch->m_spectrumFile << "\t" << currMatch->m_scanNum << "\t"
         << currMatch->m_annotation << "\t" << currMatch->m_isDecoy << endl;*/
        m_psmSet.push_back(currMatch);
      }
    }

    return true;
  }

  bool PeptideSpectrumMatchSet::loadMSGFPlusResultsFile(const char * resultsFile,
                                                        bool zeroIndexed)
  {
    vector < string > header;
    vector < vector<string> > lines;
    vector < string > requiredHeader;
    vector<int> requiredHeaderIndex;

    //set up MSGFDB headers
    requiredHeader.push_back("#SpecFile");
    requiredHeader.push_back("SpecID");
    requiredHeader.push_back("ScanNum");
    requiredHeader.push_back("Charge");
    requiredHeader.push_back("Protein");
    requiredHeader.push_back("Peptide");
    requiredHeader.push_back("QValue");

    if (!DelimitedTextReader::loadDelimitedFile(resultsFile,
                                                "\t",
                                                "",
                                                header,
                                                lines,
                                                requiredHeader,
                                                requiredHeaderIndex))
    {
      ERROR_MSG("Unable to open results! " << resultsFile);
      return false;
    }

    int scanIndex = -1;
    int specIdxIndex = -1;
    int spectrumFileIndex = -1;
    int annotationIndex = -1;
    int chargeIndex = -1;
    int scoreIndex = -1;
    int pvalueIndex = -1;
    int proteinIndex = -1;
    int precusorIndex = -1;
    int fragIndex = -1;

    //map headers
    for (int i = 0; i < header.size(); i++)
    {
      if (header[i].compare("SpecID") == 0)
      {
        specIdxIndex = i;
      }
      else if (header[i].compare("ScanNum") == 0)
      {
        scanIndex = i;
      }
      else if (header[i].compare("#SpecFile") == 0)
      {
        spectrumFileIndex = i;
      }
      else if (header[i].compare("Peptide") == 0)
      {
        annotationIndex = i;
      }
      else if (header[i].compare("Charge") == 0)
      {
        chargeIndex = i;
      }
      else if (header[i].compare("SpecEValue") == 0)
      {
        scoreIndex = i;
      }
      else if (header[i].compare("QValue") == 0)
      {
        pvalueIndex = i;
      }
      else if (header[i].compare("Protein") == 0)
      {
        proteinIndex = i;
      }
      else if (header[i].compare("Precursor") == 0)
      {
        precusorIndex = i;
      }
      else if (header[i].compare("FragMethod") == 0)
      {
        fragIndex = i;
      }
    }

    vector < string > specIDStrs;
    const char* specIDSep = "=";

    //parse results into m_psmSet
    for (int i = 0; i < lines.size(); i++)
    {
      psmPtr currMatch(new PeptideSpectrumMatch);

      if (!parseLine(currMatch,
                     lines[i],
                     zeroIndexed,
                     scanIndex,
                     -1,
                     spectrumFileIndex,
                     annotationIndex,
                     -1,
                     proteinIndex,
                     -1,
                     -1,
                     -1,
                     -1,
                     chargeIndex,
                     scoreIndex,
                     pvalueIndex,
                     -1,
                     -1,
                     -1,
                     -1,
                     -1,
                     true,
                     -1))
      {
        WARN_MSG("Unable to parse line" << i);
      }
      else
      {

        if (currMatch->m_scanNum <= 0)
        {
          splitText(lines[i][specIdxIndex].c_str(), specIDStrs, specIDSep);

          if (specIDStrs.size() > 1)
          {
            if (!sscanf(specIDStrs.back().c_str(),
                        "%d",
                        &(currMatch->m_scanNum)))
            {
              ERROR_MSG("Unable to get scan number from spectrum index!");
              return false;
            }
          }
        }
        /* cout << currMatch->m_spectrumFile << "\t" << currMatch->m_scanNum << "\t"
         << currMatch->m_annotation << endl; */
        m_psmSet.push_back(currMatch);
      }
    }

    return true;
  }

  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatchSet::loadMSGFDBResultsFile(const char * resultsFile,
                                                      bool zeroIndexed)
  {
    vector < string > header;
    vector < vector<string> > lines;
    vector < string > requiredHeader;
    vector<int> requiredHeaderIndex;

    //set up MSGFDB headers
    requiredHeader.push_back("SpecIndex");
    requiredHeader.push_back("Scan#");
    requiredHeader.push_back("Peptide");
    requiredHeader.push_back("Charge");
    requiredHeader.push_back("Protein");
    requiredHeader.push_back("P-value");
    // requiredHeader.push_back("FDR");

    if (!DelimitedTextReader::loadDelimitedFile(resultsFile,
                                                "\t",
                                                "",
                                                header,
                                                lines,
                                                requiredHeader,
                                                requiredHeaderIndex))
    {
      ERROR_MSG("Unable to open results! " << resultsFile);
      return false;
    }

    int scanIndex = -1;
    int specIdxIndex = -1;
    int spectrumFileIndex = -1;
    int annotationIndex = -1;
    int chargeIndex = -1;
    int scoreIndex = -1;
    int pvalueIndex = -1;
    int proteinIndex = -1;
    int precusorIndex = -1;
    int fragIndex = -1;
    int fdrIndex = -1;

    //map headers
    for (int i = 0; i < header.size(); i++)
    {
      if (header[i].compare("SpecIndex") == 0)
      {
        specIdxIndex = i;
      }
      else if (header[i].compare("Scan#") == 0)
      {
        scanIndex = i;
      }
      else if (header[i].compare("#SpecFile") == 0)
      {
        spectrumFileIndex = i;
      }
      else if (header[i].compare("Peptide") == 0)
      {
        annotationIndex = i;
      }
      else if (header[i].compare("Charge") == 0)
      {
        chargeIndex = i;
      }
      else if (header[i].compare("SpecProb") == 0)
      {
        scoreIndex = i;
      }
      else if (header[i].compare("P-value") == 0)
      {
        pvalueIndex = i;
      }
      else if (header[i].compare("Protein") == 0)
      {
        proteinIndex = i;
      }
      else if (header[i].compare("Precursor") == 0)
      {
        precusorIndex = i;
      }
      else if (header[i].compare("FragMethod") == 0)
      {
        fragIndex = i;
      }
      else if (header[i].compare("FDR") == 0)
      {
        fdrIndex = i;
      }
    }

    //parse results into m_psmSet
    string scanStr("");
    vector < string > clusteredScans(0);
    string specIdxStr("");
    vector < string > clusteredIdxs(0);
    // string fragStr("");
    //vector<string> clusteredFrags(0);

    const char* clustSep = "/";

    for (int i = 0; i < lines.size(); i++)
    {
      scanStr = lines[i][scanIndex];
      specIdxStr = lines[i][specIdxIndex];
      //fragStr = lines[i][fragIndex];

      splitText(scanStr.c_str(), clusteredScans, clustSep);
      splitText(specIdxStr.c_str(), clusteredIdxs, clustSep);

      if (clusteredScans.size() != clusteredIdxs.size())
      {
        ERROR_MSG("Number of clustered scans does not match the number of clustered indices for line "
            << i);
        DEBUG_MSG( lines[i][scanIndex]);
        DEBUG_MSG(lines[i][specIdxIndex]);
        return false;
      }

      vector < string > lineCopy(lines[i]);
      for (int j = 0; j < clusteredScans.size(); j++)
      {
        lineCopy[scanIndex] = clusteredScans[j];
        lineCopy[specIdxIndex] = clusteredIdxs[j];

        psmPtr currMatch(new PeptideSpectrumMatch);

        if (!parseLine(currMatch,
                       lineCopy,
                       zeroIndexed,
                       scanIndex,
                       specIdxIndex,
                       spectrumFileIndex,
                       annotationIndex,
                       -1,
                       proteinIndex,
                       -1,
                       -1,
                       -1,
                       -1,
                       chargeIndex,
                       scoreIndex,
                       pvalueIndex,
                       -1,
                       -1,
                       -1,
                       -1,
                       -1,
                       true,
                       -1,
                       -1,
                       -1,
                       -1,
                       -1,
                       -1,
                       -1,
                       -1,
                       -1,
                       -1,
                       fdrIndex))
        {
          WARN_MSG("Unable to parse line" << i);
        }
        else
        {
          /* cout << currMatch->m_spectrumFile << "\t" << currMatch->m_scanNum << "\t"
           << currMatch->m_annotation << endl; */
          m_psmSet.push_back(currMatch);
        }
      }
    }

    return true;
  }

  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatchSet::loadMSGFDBResultsFiles(const char * resultsFileList)
  {
    vector < vector<string> > msgfFileList;
    map<string, unsigned int> msgfFileListHeader;
    vector < string > requiredHeader;
    requiredHeader.push_back("Path");
    vector<int> requiredHeaderIndex;

    if (!DelimitedTextReader::loadDelimitedFile(resultsFileList,
                                                "\t",
                                                "#",
                                                msgfFileListHeader,
                                                msgfFileList,
                                                requiredHeader,
                                                requiredHeaderIndex))
    {
      ERROR_MSG("Unable to load msgf files!" << resultsFileList);
      return false;
    }

    //set up column numbers for optional parameters.
    int isZeroIndexedColumn = -1;

    map<string, unsigned int>::iterator it;

    it = msgfFileListHeader.find("isZeroIndexed");

    if (it != msgfFileListHeader.end())
    {
      isZeroIndexedColumn = it->second;
    }

    for (int i = 0; i < msgfFileList.size(); i++)
    {
      //get path
      const char * path = msgfFileList[i][requiredHeaderIndex[0]].c_str();
      DEBUG_VAR(path);

      //make sure results are loaded with correct indices.
      if (isZeroIndexedColumn != -1)
      {
        //make sure that sscanf doesn't overrun. Set to int then cast to bool.
        int isZeroIndexed = 0;
        if (!sscanf(msgfFileList[i][isZeroIndexedColumn].c_str(),
                    "%d",
                    &isZeroIndexed))
        {
          WARN_MSG("Unable to convert isZeroIndexed to bool!");
        }

        isZeroIndexed = (bool)isZeroIndexed;

        if (!this->loadMSGFDBResultsFile(path, isZeroIndexed))
        {
          ERROR_MSG("Unable to load msgf file! "
              << msgfFileList[i][requiredHeaderIndex[0]]);
          return false;
        }
      }
      else
      {
        if (!this->loadMSGFDBResultsFile(path))
        {
          ERROR_MSG("Unable to load msgf file! "
              << msgfFileList[i][requiredHeaderIndex[0]]);
          return false;
        }
      }
    }
    return true;
  }

  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatchSet::loadModaResultsFile(const char * resultsFile,
                                                    bool zeroIndexed)
  {
    vector < string > header;
    vector < vector<string> > lines;
    vector < string > requiredHeader;
    vector<int> requiredHeaderIndex;

    //set up MSGFDB headers
    //requiredHeader.push_back("Index");
    requiredHeader.push_back("ScanNo");
    requiredHeader.push_back("Peptide");
    requiredHeader.push_back("Charge");
    requiredHeader.push_back("Protein");
    requiredHeader.push_back("Probability");

    if (!DelimitedTextReader::loadDelimitedFile(resultsFile,
                                                "\t",
                                                "",
                                                header,
                                                lines,
                                                requiredHeader,
                                                requiredHeaderIndex))
    {
      ERROR_MSG("Unable to open results! " << resultsFile);
      return false;
    }

    int scanIndex = -1;
    int specIdxIndex = -1;
    int spectrumFileIndex = -1;
    int annotationIndex = -1;
    int chargeIndex = -1;
    int scoreIndex = -1;
    int pvalueIndex = -1;
    int proteinIndex = -1;
    int precusorIndex = -1;
    int fragIndex = -1;

    //map headers
    for (int i = 0; i < header.size(); i++)
    {
      if (header[i].compare("Index") == 0)
      {
        scanIndex = i;
        //    } else if (header[i].compare("ScanNo") == 0) {
        //      scanIndex = i;
      }
      else if (header[i].compare("SpectrumFile") == 0)
      {
        spectrumFileIndex = i;
      }
      else if (header[i].compare("Peptide") == 0)
      {
        annotationIndex = i;
      }
      else if (header[i].compare("Charge") == 0)
      {
        chargeIndex = i;
      }
      else if (header[i].compare("SpecProb") == 0)
      {
        scoreIndex = i;
      }
      else if (header[i].compare("Probability") == 0)
      {
        pvalueIndex = i;
      }
      else if (header[i].compare("Protein") == 0)
      {
        proteinIndex = i;
      }
      else if (header[i].compare("DeltaMass") == 0)
      { // LARS - Not sure about this
        precusorIndex = i;
      }
    }

    //parse results into m_psmSet
    for (int i = 0; i < lines.size(); i++)
    {
      psmPtr currMatch(new PeptideSpectrumMatch);

      if (!parseLine(currMatch,
                     lines[i],
                     zeroIndexed,
                     scanIndex,
                     -1,
                     spectrumFileIndex,
                     annotationIndex,
                     -1,
                     proteinIndex,
                     -1,
                     -1,
                     -1,
                     -1,
                     chargeIndex,
                     scoreIndex,
                     pvalueIndex,
                     -1,
                     -1,
                     -1,
                     -1,
                     -1,
                     true,
                     -1))
      {
        WARN_MSG("Unable to parse line" << i);
      }
      else
      {
        /*cout << currMatch->m_spectrumFile << "\t" << currMatch->m_scanNum << "\t"
         << currMatch->m_annotation << "\t" << currMatch->m_isDecoy << endl;*/
        m_psmSet.push_back(currMatch);
      }
    }

    return true;
  }

  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatchSet::loadSpecnetsReportFile(const char * resultsFile)
  {
    vector < vector<string> > lines;

    if (!DelimitedTextReader::loadDelimitedFileNoHeader(resultsFile,
                                                        ";",
                                                        "",
                                                        lines))
    {
      ERROR_MSG("Unable to open results! " << resultsFile);
      return false;
    }

    int scanIndex = 0;
    int annotationIndex = 4;
    int chargeIndex = 8;
    int proteinIndex = -1;

    //parse results into m_psmSet
    for (int i = 0; i < lines.size(); i++)
    {
      psmPtr currMatch(new PeptideSpectrumMatch);

      if (!parseLine(currMatch,
                     lines[i],
                     false,
                     scanIndex,
                     -1,
                     -1,
                     annotationIndex,
                     -1,
                     proteinIndex,
                     -1,
                     -1,
                     -1,
                     -1,
                     chargeIndex,
                     -1,
                     -1,
                     -1,
                     -1,
                     -1,
                     -1,
                     -1,
                     false,
                     -1))
      {
        WARN_MSG("Unable to parse line" << i);
      }
      else
      {
        /*cout << currMatch->m_spectrumFile << "\t" << currMatch->m_scanNum << "\t"
         << currMatch->m_annotation << "\t" << currMatch->m_isDecoy << endl;*/
        m_psmSet.push_back(currMatch);
      }
    }

    return true;
  }

  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatchSetSpectralLibraryLoader::loadSpecnetsResultsFile(const char * resultsFile,
                                                                             bool zeroIndexed)
  {
    vector < string > header;
    vector < vector<string> > lines;
    vector < string > requiredHeader;
    vector<int> requiredHeaderIndex;

    //set up inspect headers
    requiredHeader.push_back("Scan#");
    requiredHeader.push_back("Charge");

    if (!DelimitedTextReader::loadDelimitedFile(resultsFile,
                                                "\t",
                                                "",
                                                header,
                                                lines,
                                                requiredHeader,
                                                requiredHeaderIndex))
    {
      ERROR_MSG("Unable to open results! " << resultsFile);
      return false;
    }

    int scanIndex = -1;
    int spectrumFileIndex = -1;
    int annotationIndex = -1;
    int chargeIndex = -1;
    int scoreIndex = -1;
    int pvalueIndex = -1;
    int proteinIndex = -1;
    int strictenvIndex = -1;
    int unstrictenvIndex = -1;

    //map headers
    for (int i = 0; i < header.size(); i++)
    {
      if (header[i].compare("Scan#") == 0)
      {
        scanIndex = i;
      }
      else if ((header[i].compare("#SpectrumFile") == 0)
          or (header[i].compare("#SpecFile") == 0))
      {
        spectrumFileIndex = i;
      }
      else if ((header[i].compare("Annotation") == 0)
          or (header[i].compare("Peptide") == 0))
      {
        annotationIndex = i;
      }
      else if (header[i].compare("Charge") == 0)
      {
        chargeIndex = i;
      }
      else if (header[i].compare("MSGFScore") == 0)
      {
        scoreIndex = i;
      }
      //else if ((header[i].compare("p-value") == 0) or (header[i].compare("P-value") == 0))
      else if ((header[i].compare("SpecProb") == 0))
      {
        pvalueIndex = i;
      }
      else if (header[i].compare("Protein") == 0)
      {
        proteinIndex = i;
      }
      else if (header[i].compare("ScrictEnvelope") == 0)
      {
        strictenvIndex = i;
      }
      else if (header[i].compare("UnstrictEnvelope") == 0)
      {
        unstrictenvIndex = i;
      }
    }

    //parse results into m_psmSet
    for (int i = 0; i < lines.size(); i++)
    {
      psmPtr currMatch(new PeptideSpectrumMatch);

      if (!parseLine(currMatch,
                     lines[i],
                     zeroIndexed,
                     scanIndex,
                     spectrumFileIndex,
                     annotationIndex,
                     -1,
                     proteinIndex,
                     -1,
                     -1,
                     -1,
                     -1,
                     chargeIndex,
                     scoreIndex,
                     pvalueIndex,
                     strictenvIndex,
                     unstrictenvIndex))
      {
        WARN_MSG("Unable to parse line" << i);
      }
      else
      {
        /*cout << currMatch->m_spectrumFile << "\t" << currMatch->m_scanNum << "\t"
         << currMatch->m_annotation << endl; */
        m_psmSet.push_back(currMatch);
      }
    }
    return true;
  }

  bool PeptideSpectrumMatchSetSpectralLibraryLoader::parseLine(psmPtr currMatch,
                                                               vector<string> &line,
                                                               int zeroIndexed,
                                                               int scanIndex,
                                                               int spectrumFileIndex,
                                                               int annotationIndex,
                                                               int origAnnotationIndex,
                                                               int proteinIndex,
                                                               int dbIndexIndex,
                                                               int numModsIndex,
                                                               int matchOrientationIndex,
                                                               int startMassIndex,
                                                               int chargeIndex,
                                                               int scoreIndex,
                                                               int pvalueIndex,
                                                               int strictenvIndex,
                                                               int unstrictenvIndex)
  {

    if (scanIndex != -1)
    {
      if (!sscanf(line[scanIndex].c_str(), "%d", &(currMatch->m_scanNum)))
      {
        ERROR_MSG("Unable to get scan number!");
        return false;
      }
    }

    if (spectrumFileIndex != -1)
    {
      currMatch->m_spectrumFile = line[spectrumFileIndex];
    }
    if (annotationIndex != -1)
    {

      currMatch->m_annotation = line[annotationIndex];
    }
    if (origAnnotationIndex != -1)
    {
      currMatch->m_origAnnotation = line[origAnnotationIndex];
    }
    if (proteinIndex != -1)
    {
      currMatch->m_protein = line[proteinIndex];
      if (currMatch->m_protein.substr(0, 3).compare("XXX") == 0)
      {
        currMatch->m_isDecoy = true;
      }
      else
      {
        currMatch->m_isDecoy = false;
      }
    }
    if (dbIndexIndex != -1)
    {
      if (!sscanf(line[dbIndexIndex].c_str(), "%d", &(currMatch->m_dbIndex)))
      {
        ERROR_MSG("Unable to get DB index!");
        return false;
      }
    }
    if (numModsIndex != -1)
    {
      if (!sscanf(line[numModsIndex].c_str(), "%d", &(currMatch->m_numMods)))
      {
        ERROR_MSG("Unable to get num mods!");
        return false;
      }
    }
    if (matchOrientationIndex != -1)
    {
      if (!sscanf(line[matchOrientationIndex].c_str(),
                  "%d",
                  &(currMatch->m_matchOrientation)))
      {
        ERROR_MSG("Unable to get match orientation!");
        return false;
      }
    }
    if (startMassIndex != -1)
    {
      if (!sscanf(line[startMassIndex].c_str(), "%f", &(currMatch->m_startMass)))
      {
        ERROR_MSG("Unable to get start mass!");
        return false;
      }
    }
    if (chargeIndex != -1)
    {
      if (!sscanf(line[chargeIndex].c_str(), "%d", &(currMatch->m_charge)))
      {
        ERROR_MSG("Unable to get charge!");
        return false;
      }
    }
    if (scoreIndex != -1)
    {
      if (!sscanf(line[scoreIndex].c_str(), "%f", &(currMatch->m_score)))
      {
        ERROR_MSG("Unable to get score!");
        return false;
      }
    }
    if (pvalueIndex != -1)
    {
      float temp;
      if (!sscanf(line[pvalueIndex].c_str(), "%f", &(temp)))
      {
        ERROR_MSG("Unable to get p-value");
        return false;
      }
      currMatch->m_pValue = temp;
    }
    if (strictenvIndex != -1)
    {
      if (!sscanf(line[strictenvIndex].c_str(),
                  "%f",
                  &(currMatch->m_strict_envelope_score)))
      {
        ERROR_MSG("Unable to get strict envelope score");
        return false;
      }
    }
    if (unstrictenvIndex != -1)
    {
      if (!sscanf(line[unstrictenvIndex].c_str(),
                  "%f",
                  &(currMatch->m_unstrict_envelope_score)))
      {
        ERROR_MSG("Unable to get unstrict envelope score");
        return false;
      }
    }

    return true;
  }

} // namespace specnets
