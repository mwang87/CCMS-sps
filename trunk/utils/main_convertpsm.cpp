//
// util_parsimony - Create parsimonious PSM set
//
#include "Logger.h"
#include "PeptideSpectrumMatchSet.h"
#include <stdlib.h>

using namespace specnets;
using namespace std;

// -------------------------------------------------------------------------
int main(int argc, char ** argv)
{
  Logger::setDefaultLogger(Logger::getLogger(0));

  if (argc != 4)
  {
    cerr << "Usage: convertpsm input_type PSM_file_in PSM_file_out" << endl;
    cerr << "       Valid types are: inspect msgfdb moda" << endl << endl;
    cerr << "       Converts various PSM file types to the Specnets standard format" << endl;
    return -1;
  }
  
  PeptideSpectrumMatchSet psmSet;

  string type(argv[1]);
  DEBUG_VAR(type);

  DEBUG_MSG("Loading file [" << argv[2] << "] as type [" << type << "]");
  if (type == "inspect") {
    if (!psmSet.loadInspectResultsFile(argv[2])) {
      ERROR_MSG("Unable to load [" << argv[2] << "] as type [inspect]");
      return -1;
    }
  } else if (type == "msgfdb") {
    if (!psmSet.loadMSGFDBResultsFile(argv[2])) {
      ERROR_MSG("Unable to load [" << argv[2] << "] as type [msgfdb]");
      return -1;
    }
  } else if (type == "moda") {
    if (!psmSet.loadModaResultsFile(argv[2])) {
      ERROR_MSG("Unable to load [" << argv[2] << "] as type [moda]");
      return -1;
    }
  } else {
      ERROR_MSG("Unknown type [" << argv[2] << "]");
      return -1;
  }
  //DEBUG_VAR(psmSet.size());

  DEBUG_MSG("Saving file [" << argv[3] << "]");
  psmSet.saveToFile(argv[3]);

  return 0;
}
