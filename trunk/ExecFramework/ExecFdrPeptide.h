#ifndef __EXECFDRPEPTIDE_H__
#define __EXECFDRPEPTIDE_H__

// External Includes
#include "ExecBase.h"
#include "ParameterList.h"
#include "DelimitedTextReader.h"
#include "PeptideSpectrumMatch.h"
#include "PeptideSpectrumMatchSet.h"
#include "utils.h"
#include "FdrPeptide.h"

// Module Includes
#include "Logger.h"
#include "FileUtils.h"

// System Includes
#include <stdio.h>
#include <string>
#include <vector>

namespace specnets
{
  /*! \brief Class to generate statistics on a per spectrum basis
   *
   This class is used to generate the statistical input to SVM models
   */
  class ExecFdrPeptide : public ExecBase
  {
  public:
    //! \name CONSTRUCTORS
    //@{
    /*! \brief The exemplar constructor.

     Generally this constructor should not be used. It is only used by the
     module execution factory in order to create an exemplar object (without
     valid parameters) which is then used to create a real (valid) object
     using the clone() method.
     @sa clone()
     */
    ExecFdrPeptide(void);

    /*! \brief The default constructor

     This is the default constructor. A valid set of parameters must be
     supplied in the input parameter. The parameters can then be verified
     using the validateParams() method.
     @sa validateParams()
     @param inputParams structure containing all input parameters necessary for execution
     */
    ExecFdrPeptide(const ParameterList & inputParams);

    /*! \brief Constructor which takes both the parameters and the data structures

     Constructor which takes both the parameters and the data structures necessary for
     generating FDR output. When using this constructor no external data need be read in.
     The parameters can be verified using the validateParams() method.
     @sa validateParams()
     @param inputParams structure containing all input parameters necessary for execution
     @param peptideResults Input PSMSet to be considered
     @param fdrPeptides Output PSMSet for FDR values
     */
    ExecFdrPeptide(const ParameterList & inputParams,
                   PeptideSpectrumMatchSet * peptideResults,
                   PeptideSpectrumMatchSet * fdrPeptides);

    /*! \brief Constructor which takes both the parameters and the data structures

     Constructor which takes both the parameters and the data structures necessary for
     generating FDR output. When using this constructor no external data need be read in.
     The parameters can be verified using the validateParams() method. This constructor
     takes in separate target and decoy results. Note: target and decoy results MUST have
     m_spectrum parameter set for all spectra!
     @sa validateParams()
     @param inputParams structure containing all input parameters necessary for execution
     @param targetResults Input target PSMSet to be considered
     @param decoyResults Input decoy PSMSet to be considered
     @param fdrPeptides Output PSMSet for FDR values
     */
    ExecFdrPeptide(const ParameterList & inputParams,
                   PeptideSpectrumMatchSet * targetResults,
                   PeptideSpectrumMatchSet * decoyResults,
                   PeptideSpectrumMatchSet * fdrPeptides);


    //@}

    //! \name DESTRUCTOR
    //@{
    virtual ~ExecFdrPeptide(void);
    //@}

    //! \name ACCESSORS
    //@{
    /*! \brief Creates a new module of the virtual class with the given params.

     The parameters should be sufficient for the derived class to be invoked properly.

     @return A pointer to the newly created object
     */
    virtual ExecBase * clone(const ParameterList & inputParams) const;
    //! \name MODIFIERS
    //@{
    /*! \brief Executes the module.

     In order to call this method succesfully all the necessary data for
     execution must already be loaded into memory (data members). This can
     be accomplished using the loadInputData() method.

     @return True if execution finished successfully, false otherwise.
     @sa loadInputData()
     */
    virtual bool invoke(void);

    /*! \brief Loads the input data from the files specified in the params.

     Loads all the data from files into the data members. This method is
     primarily used by the execution module to load necessary data when
     executing in a separate process..

     @return True if data was loaded successfully, false otherwise.
     @sa ExecFdrPeptide(const ParameterList & input_params), saveOutputData()
     */
    virtual bool loadInputData(void);

    /*! \brief Saves all the result data to files specified in the params.

     Saves all the result data generated by the module to files specified
     by values in the params. This is used to either save the data permanantly,
     or to be loaded back in after remote execution is finished and results
     need to be merged by the merge() method.

     @param filenames A list of file names that contain the data necessary to run the module
     @return True if data was saved successfully, false otherwise.
     @sa ExecFdrPeptide(const ParameterList & input_params), loadInputData(), merge()
     */
    virtual bool saveOutputData(void);

    /*! \brief Saves all the internal data into files specified in the params.

     Saves all the data required for an external process to execute the
     module. The external process would call loadInputData() to reload the
     data into the members before calling invoke(). The user passes a vector
     that will contain the names of all files necessary to run the module as
     a separate process. The first file in this list will always be the main
     parameter file.

     @param filenames A list of file names that contain the data necessary to run the module
     @return True if data was saved successfully, false otherwise.
     @sa ExecFdrPeptide(const ParameterList & input_params), loadInputData()
     */
    virtual bool saveInputData(std::vector<std::string> & filenames);

    /*! \brief Loads the output data from the files specified in the params.

     Loads all the data from output files into the data members. The purpose
     of this is to ready "child" modules to be merged back together after
     being executed separately.

     @return True if data was loaded successfully, false otherwise.
     @sa ExecFdrPeptide(const ParameterList & input_params), saveOutputData()
     */
    virtual bool loadOutputData(void);

    /*! \brief Splits the module into multiple "children" for parallel execution.

     Divides the work required by the module into a vector of sub-modules
     that can be executed in parallel (by the the ParallelExecution() class.
     The split method should divide the work into "nodes" and sub-divide into
     "cpus" for each node. Because the precise implementation of parallel execution
     is not known apriori, nodes may be thought of as large batches, while the
     number of CPUs can be thought of as modules which will be run simultaneously
     within those batches. For example, split may be called with numNodes = 1
     and numCpus = 4 on a quad-processor machine to prep for a threaded implementation.
     Or it may be called with numNodes = 4 and numCpus = 1 to prepare for
     distributed execution on 4 separate machines on a network.

     @param numNodes Number of separate nodes
     @param numCpus Number of CPUs per node
     @return The set of sub-modules (children) that the original module has been split into.
     An empty vector implies an error.
     @sa merge()
     */
    virtual std::vector<ExecBase *> const & split(int numSplit);

    /*! \brief Merges the child modules back into a complete result

     This method is only called when split is used for parallel execution
     and after the execution of all children has been completed. This method
     will merge the results generated by each of the child modules into
     one cohesive result as if the module had been run as a single entity.

     @return True if merge could be performed successfully, false otherwise.
     @sa split(int numNodes, int numCpus)
     */
    virtual bool merge(void);

    /*! \brief Performs validation of the input parameters.

     Checks the parameters structure provided in the constructor to see if they
     are sufficient and correct to invoke the module. Also sets the internal
     validity flag so that isValid() will return the correct result.

     @param error A description of the error (if any occurs)
     @return True if the parameters for the module are valid, false otherwise.
     @sa isValid()
     */
    virtual bool validateParams(std::string & error);
    //@}

  private:
    bool ownInput; //! Does this object "own" the input data pointers (and hence have to free them)

    PeptideSpectrumMatchSet * m_peptideResults; //! The peptide results we're considering.
    PeptideSpectrumMatchSet * m_targetResults; //! Target peptide results
    PeptideSpectrumMatchSet * m_decoyResults; //! Decoy peptide results
    PeptideSpectrumMatchSet * m_fdrPeptides; //! FDR output

    bool ownOutput; //! Does this object "own" the output data pointers (and hence have to free them)

  };

} // namespace specnets

#endif // __ExecFdrPeptide_H__
