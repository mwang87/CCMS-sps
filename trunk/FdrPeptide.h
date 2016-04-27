#ifndef FDR_PEPTIDE_H_
#define FDR_PEPTIDE_H_

// External Includes
#include "PeptideSpectrumMatch.h"
#include "PeptideSpectrumMatchSet.h"
#include "utils.h"

// Module Includes
#include "Logger.h"

// System Includes
#include <stdio.h>
#include <string>
#include <vector>

//TR1 includes. GCC 4.0 and above only!
#ifdef __GLIBCXX__
#  include <tr1/unordered_map>
#else
#  ifdef __IBMCPP__
#    define __IBMCPP_TR1__
#  endif
#  include <unordered_map>
#endif

namespace specnets
{
  /*! \brief false discovery rates class for peptides. Input can be single PSM set or separate decoy
   * and target peptides.
   *
   */
  class FdrPeptide
  {
  public:
    /** Calculates p-values for input psms when decoy and target search is done together
     *
     * @param inputPeptides - input psm set. Needs to have m_isDecoy flag set for decoy hits!
     * @param outputPeptides - output psm set. Calculates m_pValue.
     * @param scaling factor - True database fraction / false database fraction
     * @param comparison function - Function used to compare PSMs
     * @return - true if function completes correctly.
     */
   static bool concatenatedTargetDecoy(PeptideSpectrumMatchSet &inputPeptides,
                                       PeptideSpectrumMatchSet &outputPeptides,
                                       double scalingFactor = 1,
                                       bool (*fdrCompareFunction)(psmPtr, psmPtr) = 0);
    /** Calculates p-values for input psms when decoy and target search is done separately
        *
        * @param inputPeptides - input psm set. Needs to have m_isDecoy flag set for decoy hits!
        * @param outputPeptides - output psm set. Calculates m_pValue.
        * @param scaling factor - True database fraction / false database fraction
        * @param comparison function - Function used to compare PSMs
        * @return - true if function completes correctly.
        */
    static bool separateTargetDecoy(PeptideSpectrumMatchSet &inputPeptides,
                                    PeptideSpectrumMatchSet &outputPeptides,
                                    double scalingFactor = 1,
                                    bool (*fdrCompareFunction)(psmPtr, psmPtr) = 0);
    /** Calculates p-values for input psms when decoy and target search is done together
     *
     * @param inputPeptides - input psm set. Needs to have m_isDecoy flag set for decoy hits!
     * @param outputPeptides - output psm set. Calculates m_pValue.
     * @param scaling factor - True database fraction / false database fraction
     * @param comparison function - Function used to compare PSMs
     * @return - true if function completes correctly.
     */
    static bool calculatePValues(PeptideSpectrumMatchSet &inputPeptides,
                                 PeptideSpectrumMatchSet &outputPeptides,
                                 double scalingFactor = 1,
                                 bool (*fdrCompareFunction)(psmPtr, psmPtr) = 0);

    /** Merges psm sets when decoy and target search is done separately
     *
     * @param targetPeptides - input psm set for target matches.
     * @param decoyPeptides - input psm set for decoy matches.
     * @param outputPeptides - output psm set. Includes both decoy and target hits.
     * @return - true if function completes correctly.
     */
    static bool mergeTargetDecoy(PeptideSpectrumMatchSet &targetPeptides,
                                 PeptideSpectrumMatchSet &decoyPeptides,
                                 PeptideSpectrumMatchSet &outputPeptides);

    /** filters psm set to top scoring hit
     *
     * @param inputPeptides - input psm set.
     * @param outputPeptides - output psm set.
     * @return - true if function completes correctly.
     */
    static bool filterToTopHit(PeptideSpectrumMatchSet &inputPeptides,
                               PeptideSpectrumMatchSet &outputPeptides);

    /** concatenates target and decoy psm set.
     *
     * @param targetPeptides - input psm set for target matches.
     * @param decoyPeptides - input psm set for decoy matches.
     * @param outputPeptides - output psm set. Includes both decoy and target hits.
     * @return - true if function completes correctly.
     */
    static void concatenateTargetDecoy(PeptideSpectrumMatchSet &targetPeptides,
                                       PeptideSpectrumMatchSet &decoyPeptides,
                                       PeptideSpectrumMatchSet &outputPeptides);

    /** Takes in input psm set with pValues already set and filters based on cutoff
     *
     * @param inputPeptides - input psm set. m_pValue must be set!
     * @param cutoff - p-value cutoff.
     * @return - true if function completes correctly. False if m_pValue not set on input.
     */
    static bool filterByPValue(PeptideSpectrumMatchSet &inputPeptides,
                               double cutoff);
  };
}

#endif
