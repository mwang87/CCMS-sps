///////////////////////////////////////////////////////////////////////////////
#include "SpecplotInterface2.h"
#include "mzxml.h"
#include "PWizInterface.h"

///////////////////////////////////////////////////////////////////////////////
namespace specnets {
///////////////////////////////////////////////////////////////////////////////
SpecplotInterface2::SpecplotInterface2()
{
}

SpecplotInterface2::~SpecplotInterface2()
{
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
int SpecplotInterface2::load(string &fn, string &ext, bool spectrumScanPresent, bool spectrumIdPresent, bool usePwizFirst)
{
  bool fileLoaded = false;
  // load te file
  if(!specSet) {
    // create object
    specSet = new SpecSet();
    
    // load using pwiz
    if(usePwizFirst && !fileLoaded && (ext.compare("mzxml") || ext.compare("mzml"))) {
      PWizInterface pwiz;
      fileLoaded = pwiz.loadDataUsingPWiz(fn, *specSet, 2);
    }

    // load MGF, given and id or scan   
    if(!fileLoaded && false) {
      if(ext.compare("mgf") == 0) {
        int idx = (m_spectrumIndex > 0 ? m_spectrumIndex : 0);
        int scan = getInt(m_spectrumScan.c_str());
        scan = (scan > 0 ? scan : 0);
        fileLoaded = specSet->LoadSpecSet_mgf(fn.c_str(), scan, idx);
      }
    }

    // load mzxml file. Load using first the SPS mzxml loader. Use the appropriate method in case scan # is specified
    if(!fileLoaded) {
      if(ext.compare("mzxml") == 0) {
        // auxilizary vector needed
        vector<short> msLevel;
        // hold the return value
        int ret;
        // load method depends on the presence of spectrumscan specification
        try {
          if(m_spectrumScan.empty()) {
            fileLoaded = LoadMzxml( (char * const)(fn.c_str()), fn, *specSet, & msLevel, 2);
          } else {
            fileLoaded = LoadMzxml(fn.c_str(), fn, *specSet, m_spectrumScan.c_str(), & msLevel, 2);
          }
        } catch (...) {
          fileLoaded = false;
        }
      }
    }
    
    // If not file was loaded at this point, load it using the specsset load methods.
    if(!fileLoaded) {
      // load other formats
      fileLoaded = specSet->Load(fn.c_str(), NULL);
    }
    
    // load using pwiz - final attempt
    if(!fileLoaded) {
      PWizInterface pwiz;
      fileLoaded = pwiz.loadDataUsingPWiz(fn, *specSet, 2);
    }
    
  }

  // if no file was loaded at this point, fail
  if(!fileLoaded) {
    stringstream err;
    err << "Error loading file: " << fn;
    return error(err.str());
  }  

  // Check if spectrum ID is present. If it is, find the corresponding spectrum index 
  if(spectrumIdPresent) {
    bool found = false;
    for(int i = 0 ; i < specSet->size() ; i++){
      if(  (*specSet)[i].psmList.size() == 1 ){
        if ( (*specSet)[i].psmList.front()->m_spectrumID == m_spectrumID) {
          m_spectrumIndex = i+1;
          plotSpectrum.setSpectrumIndex(m_spectrumIndex);
          found = true;
        }
      }
    }
    // in case we didn't find it, exit in error
    if(!found) {
        stringstream err; err << "ERROR: Spectrum ID " << m_spectrumID << " not found.";
        return error(err.str());
    }
  }
    
  // If spectrumscan was specified,   
  if(spectrumScanPresent) {
    // get spectrumscan from index
    bool found = false;
    for(int i = 0 ; i < specSet->size() ; i++){
      if((*specSet)[i].scan == getInt(m_spectrumScan.c_str())) {
        m_spectrumIndex = i+1;
        plotSpectrum.setSpectrumIndex(m_spectrumIndex);
        found = true;
      }
    }
    // in case we didn't find it, exit in error
    if(!found) {
      stringstream err; err << "ERROR: Spectrumscan not found.";
      return error(err.str());
    }
  }
  
  if(fileLoaded)
    m_inputSpectraFilename = fn;

  return fileLoaded;
}
///////////////////////////////////////////////////////////////////////////////
}; // namespace specnets
///////////////////////////////////////////////////////////////////////////////
