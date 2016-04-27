//
//  main_specdump - stand alone executable for dumping spectrum data
//
#include "CommandLineParser.h"
#include "Logger.h"
#include "spectrum.h"
#include "SpecSet.h"
#include "SpectrumPairSet.h"
#include <stdlib.h>

using namespace specnets;
using namespace std;

// -------------------------------------------------------------------------
int main(int argc, char ** argv)
{
  Logger::setDefaultLogger(Logger::getLogger(0));

  DEBUG_TRACE;
  if (argc != 3)
  {
    cerr << "Usage: main_specdump type specfile" << endl;
    cerr << "       Valid types are: mgf, prms, specset, specpairset, pklbin" << endl;
    return -1;
  }
  
  string type = argv[1];

  if (type == "mgf" || type == "prms" || type == "specset" || type == "pklbin")
  {
    SpecSet spectra1;
    size_t size1 = 0;


    if (type == "mgf")
    {
      DEBUG_TRACE;
      spectra1.LoadSpecSet_mgf(argv[2]);
      size1 = spectra1.size();
      DEBUG_VAR(size1);
    }
    else if (type == "prms")
    {
      DEBUG_TRACE;
      spectra1.LoadSpecSet_prmsv3(argv[2]);
      size1 = spectra1.size();
      DEBUG_VAR(size1);
    }
    else if (type == "specset" || type == "pklbin")
    {
      DEBUG_TRACE;
      spectra1.loadPklBin(argv[2]);
      size1 = spectra1.size();
      DEBUG_VAR(size1);
    }

    for (size_t i = 0; i < size1; i++)
    {
      cout << "i = " << i << endl;
      
      //DEBUG_VAR(i);
      size_t peakSize1 = spectra1[i].size();
      //DEBUG_VAR(peakSize1);

      cout << "Parent Mass = " << spectra1[i].parentMass << endl;
      cout << "Scan Number = " << spectra1[i].scan << endl;
      cout << "MS Level = " << spectra1[i].msLevel << endl;
      
      for (size_t j = 0; j < peakSize1; j++)
      {
        TwoValues<float> peak1 = spectra1[i][j];
        //DEBUG_VAR(peak1.values[0]);
        //DEBUG_VAR(peak1.values[1]);
        cout << peak1.values[0] << ", " << peak1.values[1] << endl;
      }  
    }
  }
  else if (type == "specpairset")
  {
    SpectrumPairSet specpairset1;
    specpairset1.loadFromBinaryFile(argv[2]);
    size_t size1 = specpairset1.size();
    DEBUG_VAR(size1);

    for (size_t i = 0; i < size1; i++)
    {
      cout << "i =" << i << endl;
      SpectrumPair pair1 = specpairset1[i];

      cout << "pair1.spec1 = " << pair1.spec1 << endl;
      cout << "pair1.spec2 = " << pair1.spec2 << endl;
      cout << "pair1.score1 = " << pair1.score1 << endl;
      cout << "pair1.score2 = " << pair1.score2 << endl;
      cout << "pair1.shift1 = " << pair1.shift1 << endl;
      cout << "pair1.shift2 = " << pair1.shift2 << endl;
      cout << "pair1.specC = " << pair1.specC << endl;
      cout << "pair1.spec2rev = " << pair1.spec2rev << endl;
    }
  }
  else
  {
    cerr << "Usage: main_specdump type specfile" << endl;
    cerr << "       Valid types are: mgf, prms, specset, specpairset" << endl;
    return -1;
  }  
  
  return 0;
}

